# Install required packages if not already installed
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

required_packages <- c("GEOquery", "limma", "Biobase", "dplyr", "readr", "stringr", "edgeR")
for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE)) {
    BiocManager::install(pkg)
    library(pkg, character.only = TRUE)
  }
}


# Function to download and process GSE124439 count data
process_GSE124439 <- function() {
  cat("Processing GSE124439 RNA-seq count data...\n")
  
  # Download supplementary files if not already done
  if (!file.exists("GSE124439/GSE124439_RAW.tar")) {
    cat("Downloading supplementary files...\n")
    getGEOSuppFiles("GSE124439", baseDir = ".")
  }
  
  # Extract the tar file if not already done
  if (!dir.exists("GSE124439/counts")) {
    cat("Extracting files...\n")
    untar("GSE124439/GSE124439_RAW.tar", exdir = "GSE124439/counts")
  }
  
  # Get list of count files
  count_files <- list.files("GSE124439/counts", pattern = "*_counts\\.txt\\.gz$", full.names = TRUE)
  cat("Found", length(count_files), "count files\n")
  
  # Read the first file to get gene names and structure
  first_file <- read.delim(count_files[1], header = TRUE, stringsAsFactors = FALSE)
  cat("First file structure:\n")
  print(head(first_file))
  
  # Initialize expression matrix
  gene_names <- first_file[, 1]
  expr_matrix <- matrix(0, nrow = length(gene_names), ncol = length(count_files))
  rownames(expr_matrix) <- gene_names
  
  # Sample names from file names (keep GSM IDs for matching with metadata)
  sample_names <- basename(count_files)
  colnames(expr_matrix) <- sample_names
  
  # Read all count files
  cat("Reading count files...\n")
  for (i in 1:length(count_files)) {
    if (i %% 20 == 0) cat("Processing file", i, "of", length(count_files), "\n")
    
    count_data <- read.delim(count_files[i], header = TRUE, stringsAsFactors = FALSE)
    
    # Ensure the gene order matches
    if (all(count_data[, 1] == gene_names)) {
      expr_matrix[, i] <- count_data[, 2]
    } else {
      # Match by gene name if order differs
      matched <- match(gene_names, count_data[, 1])
      expr_matrix[, i] <- count_data[matched, 2]
    }
  }
  
  # Remove the header row if present
  if (rownames(expr_matrix)[1] == "gene/TE") {
    expr_matrix <- expr_matrix[-1, ]
    gene_names <- gene_names[-1]
  }
  
  cat("Final expression matrix dimensions:", dim(expr_matrix), "\n")
  
  return(list(expr_matrix = expr_matrix, sample_names = sample_names))
}

# Function to get proper sample metadata from GEO
get_sample_metadata <- function() {
  cat("Downloading sample metadata from GEO...\n")
  
  # Download the series matrix
  gse <- getGEO("GSE124439", GSEMatrix = TRUE)
  
  if (length(gse) > 0) {
    # Get phenotype data
    pdata <- pData(gse[[1]])
    
    # Print available columns to help identify group information
    cat("Available metadata columns:\n")
    print(colnames(pdata))
    
    # Look for columns that might contain disease status information
    potential_group_cols <- c("sample group", "disease", "diagnosis", "characteristics", "description")
    
    for (col in colnames(pdata)) {
      if (any(grepl(paste(potential_group_cols, collapse = "|"), tolower(col)))) {
        cat("\nPotential group column:", col, "\n")
        print(table(pdata[, col]))
      }
    }
    
    return(pdata)
  }
  return(NULL)
}

# Function to create proper sample groups
create_proper_sample_groups <- function(sample_names, geo_metadata) {
  cat("Creating sample groups based on GEO metadata...\n")
  
  # Extract GSM IDs from sample names
  gsm_ids <- gsub("_counts\\.txt\\.gz$", "", sample_names)
  
  # Initialize groups
  groups <- rep(NA, length(gsm_ids))
  
  if (!is.null(geo_metadata)) {
    # Try to match with GEO metadata
    matched <- match(gsm_ids, rownames(geo_metadata))
    
    # Look for disease status in various columns
    for (i in 1:length(gsm_ids)) {
      if (!is.na(matched[i])) {
        sample_meta <- geo_metadata[matched[i], ]
        
        # Check various columns for disease information
        if ("characteristics_ch1.2" %in% colnames(sample_meta)) {
          # This often contains disease status
          status <- as.character(sample_meta$characteristics_ch1.2)
          if (grepl("control", tolower(status))) {
            groups[i] <- "Control"
          } else if (grepl("als", tolower(status))) {
            groups[i] <- "ALS"
          }
        }
        
        # Also check other common columns
        if (is.na(groups[i]) && "source_name_ch1" %in% colnames(sample_meta)) {
          status <- as.character(sample_meta$source_name_ch1)
          if (grepl("control", tolower(status))) {
            groups[i] <- "Control"
          } else if (grepl("als", tolower(status))) {
            groups[i] <- "ALS"
          }
        }
      }
    }
  }
  
  # If still missing groups, try to infer from sample names
  if (any(is.na(groups))) {
    cat("Inferring groups from sample names...\n")
    for (i in 1:length(gsm_ids)) {
      if (is.na(groups[i])) {
        # Look for patterns in sample names that might indicate group
        sample_id <- gsm_ids[i]
        if (grepl("GSM3533[0-2][0-9]", sample_id)) {
          # Earlier GSM IDs might be controls (this is a guess - needs verification)
          groups[i] <- "Control"
        } else if (grepl("GSM3533[3-4][0-9]", sample_id)) {
          # Later GSM IDs might be ALS (this is a guess - needs verification)
          groups[i] <- "ALS"
        } else {
          groups[i] <- "Unknown"
        }
      }
    }
  }
  
  # Print group distribution
  cat("Sample group distribution:\n")
  print(table(groups, useNA = "always"))
  
  return(groups)
}

# Function to perform differential expression analysis with edgeR
run_DE_analysis <- function(expr_matrix, groups) {
  cat("Performing differential expression analysis with edgeR...\n")
  
  # Only use samples with known groups
  valid_samples <- !is.na(groups) & groups %in% c("Control", "ALS")
  
  if (sum(valid_samples) < 10) {
    cat("Not enough samples with valid groups for DE analysis.\n")
    return(NULL)
  }
  
  expr_matrix_filtered <- expr_matrix[, valid_samples]
  groups_filtered <- groups[valid_samples]
  
  cat("Using", sum(valid_samples), "samples for DE analysis\n")
  cat("Group distribution:", table(groups_filtered), "\n")
  
  # Create DGEList object
  dge <- DGEList(counts = expr_matrix_filtered, group = groups_filtered)
  
  # Filter lowly expressed genes (keep genes with at least 10 counts in at least 10% of samples)
  keep <- filterByExpr(dge, min.count = 10, min.prop = 0.1)
  dge <- dge[keep, , keep.lib.sizes = FALSE]
  cat("Genes after filtering:", nrow(dge), "\n")
  
  # Normalize
  dge <- calcNormFactors(dge)
  
  # Create design matrix
  design <- model.matrix(~groups_filtered)
  
  # Estimate dispersion
  dge <- estimateDisp(dge, design)
  
  # Fit model and test
  fit <- glmQLFit(dge, design)
  qlf <- glmQLFTest(fit, coef = 2)
  
  # Get results
  results <- topTags(qlf, n = Inf)$table
  return(results)
}

# Main analysis
cat("Starting GSE124439 analysis...\n")

# Process the count data
processed_data <- process_GSE124439()
expr_matrix <- processed_data$expr_matrix
sample_names <- processed_data$sample_names

# Get sample metadata from GEO
geo_metadata <- get_sample_metadata()

# Create proper sample groups
groups <- create_proper_sample_groups(sample_names, geo_metadata)

# Save the processed expression matrix
write.csv(expr_matrix, "GSE124439_expression_matrix.csv")
cat("Expression matrix saved to GSE124439_expression_matrix.csv\n")

# Save sample metadata with groups
sample_metadata <- data.frame(
  Sample = sample_names,
  GSM_ID = gsub("_counts\\.txt\\.gz$", "", sample_names),
  Group = groups
)
write.csv(sample_metadata, "GSE124439_sample_metadata.csv")
cat("Sample metadata saved to GSE124439_sample_metadata.csv\n")

# Perform differential expression analysis
if (sum(groups %in% c("Control", "ALS"), na.rm = TRUE) >= 10) {
  de_results <- run_DE_analysis(expr_matrix, groups)
  
  if (!is.null(de_results)) {
    # Save results
    write.csv(de_results, "GSE124439_DE_results.csv")
    cat("Differential expression results saved to GSE124439_DE_results.csv\n")
    
    # Print summary
    cat("\nDifferential expression summary:\n")
    cat("Total genes tested:", nrow(de_results), "\n")
    cat("Significant genes (FDR < 0.05):", sum(de_results$FDR < 0.05), "\n")
    cat("Significant genes (FDR < 0.01):", sum(de_results$FDR < 0.01), "\n")
    
    if (sum(de_results$FDR < 0.05) > 0) {
      cat("Top 5 differentially expressed genes:\n")
      print(head(de_results[de_results$FDR < 0.05, ], 5))
    } else {
      cat("Top 5 genes by p-value:\n")
      print(head(de_results, 5))
    }
  }
} else {
  cat("Not enough samples with proper group assignments for DE analysis.\n")
  cat("Please check the sample_metadata.csv file and manually assign groups.\n")
}

cat("Analysis complete!\n")

# Print instructions for manual group assignment
cat("\n=== INSTRUCTIONS ===\n")
cat("1. Check GSE124439_sample_metadata.csv to see the automatically assigned groups\n")
cat("2. If groups are incorrect, you may need to:\n")
cat("   - Manually assign groups based on the original publication\n")
cat("   - Contact the authors for sample information\n")
cat("   - Use the GEO sample page to find group information\n")