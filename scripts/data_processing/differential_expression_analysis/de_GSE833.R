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


# Function to download and process GSE833 data
process_GSE833 <- function() {
  cat("Processing GSE833 data...\n")
  
  # Download supplementary files if not already done
  if (!file.exists("GSE833/GSE833_RAW.tar")) {
    cat("Downloading supplementary files...\n")
    getGEOSuppFiles("GSE833", baseDir = ".")
  }
  
  # Check what files were downloaded
  downloaded_files <- list.files("GSE833", pattern = "\\.tar$|\\.zip$|\\.txt\\.gz$", full.names = TRUE)
  cat("Downloaded files:", downloaded_files, "\n")
  
  # Extract files if they exist
  if (length(downloaded_files) > 0) {
    for (file in downloaded_files) {
      if (grepl("\\.tar$", file)) {
        cat("Extracting tar file:", file, "\n")
        untar(file, exdir = "GSE833/raw_data")
      } else if (grepl("\\.zip$", file)) {
        cat("Extracting zip file:", file, "\n")
        unzip(file, exdir = "GSE833/raw_data")
      }
    }
  }
  
  # Look for data files in the extracted directory
  data_files <- list.files("GSE833/raw_data", pattern = "\\.txt$|\\.csv$|\\.xls$|\\.cel$", 
                          full.names = TRUE, recursive = TRUE)
  
  # If no extracted files, check for direct data files
  if (length(data_files) == 0) {
    data_files <- list.files("GSE833", pattern = "\\.txt\\.gz$|\\.csv\\.gz$", 
                            full.names = TRUE)
  }
  
  cat("Found", length(data_files), "data files\n")
  
  if (length(data_files) == 0) {
    cat("No raw data files found. Trying to use series matrix data...\n")
    return(NULL)
  }
  
  # Check the first file to understand the format
  if (length(data_files) > 0) {
    first_file <- data_files[1]
    cat("First file:", first_file, "\n")
    
    # Try to read the file
    tryCatch({
      if (grepl("\\.cel$", first_file, ignore.case = TRUE)) {
        cat("CEL file detected - using affy package\n")
        if (!require("affy", quietly = TRUE)) {
          BiocManager::install("affy")
          library(affy)
        }
        # Read CEL files
        cel_data <- ReadAffy(celfile.path = "GSE833/raw_data")
        expr_matrix <- exprs(cel_data)
        sample_names <- sampleNames(cel_data)
        
      } else if (grepl("\\.txt$|\\.csv$", first_file)) {
        # Try to read as text file
        first_data <- read.delim(first_file, header = TRUE, stringsAsFactors = FALSE, nrows = 10)
        cat("First file structure:\n")
        print(head(first_data))
        
        # This will need to be customized based on the actual file format
        # For now, return NULL and we'll use the series matrix
        return(NULL)
      }
    }, error = function(e) {
      cat("Error reading data files:", e$message, "\n")
      return(NULL)
    })
  }
  
  return(NULL)
}

# Function to get expression data from series matrix
get_series_matrix_data <- function() {
  cat("Getting expression data from series matrix...\n")
  
  # Download the series matrix
  gse <- getGEO("GSE833", GSEMatrix = TRUE)
  
  if (length(gse) > 0) {
    gse_data <- gse[[1]]
    
    # Get expression matrix
    expr_matrix <- exprs(gse_data)
    cat("Expression matrix dimensions:", dim(expr_matrix), "\n")
    
    # Get sample names
    sample_names <- colnames(expr_matrix)
    
    return(list(expr_matrix = expr_matrix, sample_names = sample_names, gse_data = gse_data))
  }
  
  return(NULL)
}

# Function to get proper sample metadata from GEO
get_sample_metadata <- function(gse_data = NULL) {
  cat("Downloading sample metadata from GEO...\n")
  
  if (is.null(gse_data)) {
    # Download the series matrix if not provided
    gse <- getGEO("GSE833", GSEMatrix = TRUE)
    if (length(gse) > 0) {
      gse_data <- gse[[1]]
    }
  }
  
  if (!is.null(gse_data)) {
    # Get phenotype data
    pdata <- pData(gse_data)
    
    # Print available columns to help identify group information
    cat("Available metadata columns:\n")
    print(colnames(pdata))
    
    # Look for columns that might contain group information
    potential_group_cols <- c("disease", "diagnosis", "characteristics", "description", 
                             "tissue", "cell type", "treatment", "condition")
    
    for (col in colnames(pdata)) {
      if (any(grepl(paste(potential_group_cols, collapse = "|"), tolower(col)))) {
        cat("\nPotential group column:", col, "\n")
        unique_vals <- unique(pdata[, col])
        if (length(unique_vals) <= 20) {  # Only print if not too many values
          print(table(pdata[, col]))
        } else {
          cat("Too many unique values (", length(unique_vals), "), showing first 10:\n")
          print(head(unique_vals, 10))
        }
      }
    }
    
    return(pdata)
  }
  return(NULL)
}

# Function to create proper sample groups
create_proper_sample_groups <- function(sample_names, geo_metadata) {
  cat("Creating sample groups based on GEO metadata...\n")
  
  # Initialize groups
  groups <- rep(NA, length(sample_names))
  
  if (!is.null(geo_metadata)) {
    # Try to match with GEO metadata
    matched <- match(sample_names, rownames(geo_metadata))
    
    # Look for group information in various columns
    for (i in 1:length(sample_names)) {
      if (!is.na(matched[i])) {
        sample_meta <- geo_metadata[matched[i], ]
        
        # Check various columns for group information
        for (col in colnames(sample_meta)) {
          if (!is.na(groups[i])) break  # Stop if group already assigned
          
          value <- as.character(sample_meta[, col])
          if (!is.na(value) && value != "") {
            # Look for common group patterns
            if (grepl("control|normal|healthy", tolower(value))) {
              groups[i] <- "Control"
            } else if (grepl("als|motor neuron|amyotrophic", tolower(value))) {
              groups[i] <- "ALS"
            } else if (grepl("alzheimer|ad", tolower(value))) {
              groups[i] <- "Alzheimer"
            } else if (grepl("parkinson|pd", tolower(value))) {
              groups[i] <- "Parkinson"
            } else if (grepl("huntington|hd", tolower(value))) {
              groups[i] <- "Huntington"
            }
          }
        }
      }
    }
  }
  
  # If still missing groups, try to infer from sample names
  if (any(is.na(groups))) {
    cat("Inferring groups from sample names...\n")
    for (i in 1:length(sample_names)) {
      if (is.na(groups[i])) {
        sample_id <- sample_names[i]
        # Add any known patterns for GSE833 here
        if (grepl("control|ctrl|normal", tolower(sample_id))) {
          groups[i] <- "Control"
        } else if (grepl("als|motor", tolower(sample_id))) {
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

# Function to perform differential expression analysis
run_DE_analysis <- function(expr_matrix, groups, data_type = "microarray") {
  cat("Performing differential expression analysis...\n")
  
  # Only use samples with known groups
  valid_groups <- unique(groups[!is.na(groups) & groups != "Unknown"])
  if (length(valid_groups) < 2) {
    cat("Need at least 2 groups for DE analysis. Found:", valid_groups, "\n")
    return(NULL)
  }
  
  # Use the two most common groups for comparison
  group_counts <- table(groups[groups %in% valid_groups])
  top_groups <- names(sort(group_counts, decreasing = TRUE))[1:2]
  
  valid_samples <- groups %in% top_groups
  expr_matrix_filtered <- expr_matrix[, valid_samples]
  groups_filtered <- groups[valid_samples]
  
  cat("Using", sum(valid_samples), "samples for DE analysis\n")
  cat("Group distribution:", table(groups_filtered), "\n")
  cat("Comparison:", paste(top_groups, collapse = " vs "), "\n")
  
  if (data_type == "rnaseq") {
    # Use edgeR for RNA-seq data
    if (!require("edgeR", quietly = TRUE)) {
      BiocManager::install("edgeR")
      library(edgeR)
    }
    
    # Create DGEList object
    dge <- DGEList(counts = expr_matrix_filtered, group = groups_filtered)
    
    # Filter lowly expressed genes
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
    
  } else {
    # Use limma for microarray data
    # Create design matrix
    design <- model.matrix(~0 + groups_filtered)
    colnames(design) <- gsub("groups_filtered", "", colnames(design))
    
    # Make contrasts
    contrast_matrix <- makeContrasts(
      contrasts = paste(top_groups[1], top_groups[2], sep = "-"),
      levels = design
    )
    
    # Fit linear model
    fit <- lmFit(expr_matrix_filtered, design)
    fit2 <- contrasts.fit(fit, contrast_matrix)
    fit2 <- eBayes(fit2)
    
    # Get results
    results <- topTable(fit2, number = Inf, adjust.method = "BH")
  }
  
  return(results)
}

# Main analysis
cat("Starting GSE833 analysis...\n")

# Try to get data from series matrix first (most reliable for older GEO datasets)
series_data <- get_series_matrix_data()

if (!is.null(series_data)) {
  expr_matrix <- series_data$expr_matrix
  sample_names <- series_data$sample_names
  gse_data <- series_data$gse_data
  data_type <- "microarray"  # GSE833 is likely microarray data
  
  cat("Successfully obtained expression data from series matrix\n")
} else {
  # Try to process raw data files
  processed_data <- process_GSE833()
  if (!is.null(processed_data)) {
    expr_matrix <- processed_data$expr_matrix
    sample_names <- processed_data$sample_names
    data_type <- "rnaseq"  # Assuming RNA-seq if raw files are available
  } else {
    stop("Could not obtain expression data for GSE833")
  }
}

# Get sample metadata
geo_metadata <- get_sample_metadata(if(exists("gse_data")) gse_data else NULL)

# Create proper sample groups
groups <- create_proper_sample_groups(sample_names, geo_metadata)

# Save the processed expression matrix
write.csv(expr_matrix, "GSE833_expression_matrix.csv")
cat("Expression matrix saved to GSE833_expression_matrix.csv\n")

# Save sample metadata with groups
sample_metadata <- data.frame(
  Sample = sample_names,
  Group = groups
)
write.csv(sample_metadata, "GSE833_sample_metadata.csv")
cat("Sample metadata saved to GSE833_sample_metadata.csv\n")

# Perform differential expression analysis
if (length(unique(groups[!is.na(groups) & groups != "Unknown"])) >= 2) {
  de_results <- run_DE_analysis(expr_matrix, groups, data_type)
  
  if (!is.null(de_results)) {
    # Save results
    write.csv(de_results, "GSE833_DE_results.csv")
    cat("Differential expression results saved to GSE833_DE_results.csv\n")
    
    # Print summary
    cat("\nDifferential expression summary:\n")
    cat("Total genes tested:", nrow(de_results), "\n")
    
    if ("adj.P.Val" %in% colnames(de_results)) {
      # limma results
      sig_genes <- sum(de_results$adj.P.Val < 0.05, na.rm = TRUE)
      cat("Significant genes (FDR < 0.05):", sig_genes, "\n")
      cat("Significant genes (FDR < 0.01):", sum(de_results$adj.P.Val < 0.01, na.rm = TRUE), "\n")
      
      if (sig_genes > 0) {
        cat("Top 5 differentially expressed genes:\n")
        print(head(de_results[de_results$adj.P.Val < 0.05, ], 5))
      } else {
        cat("Top 5 genes by p-value:\n")
        print(head(de_results, 5))
      }
    } else if ("FDR" %in% colnames(de_results)) {
      # edgeR results
      sig_genes <- sum(de_results$FDR < 0.05, na.rm = TRUE)
      cat("Significant genes (FDR < 0.05):", sig_genes, "\n")
      cat("Significant genes (FDR < 0.01):", sum(de_results$FDR < 0.01, na.rm = TRUE), "\n")
      
      if (sig_genes > 0) {
        cat("Top 5 differentially expressed genes:\n")
        print(head(de_results[de_results$FDR < 0.05, ], 5))
      } else {
        cat("Top 5 genes by p-value:\n")
        print(head(de_results, 5))
      }
    }
  }
} else {
  cat("Not enough samples with proper group assignments for DE analysis.\n")
  cat("Please check the sample_metadata.csv file and manually assign groups.\n")
}

cat("Analysis complete!\n")

# Print instructions for manual group assignment
cat("\n=== INSTRUCTIONS ===\n")
cat("1. Check GSE833_sample_metadata.csv to see the automatically assigned groups\n")
cat("2. If groups are incorrect, you may need to:\n")
cat("   - Manually assign groups based on the original publication\n")
cat("   - Check the GEO sample page for detailed sample information\n")
cat("   - Modify the create_proper_sample_groups() function with known group patterns\n")