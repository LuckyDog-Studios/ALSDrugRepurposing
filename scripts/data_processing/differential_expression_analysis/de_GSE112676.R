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

# Function to get expression data from series matrix
get_series_matrix_data <- function() {
  cat("Getting expression data from series matrix...\n")
  
  # Download the series matrix
  gse <- getGEO("GSE112676", GSEMatrix = TRUE)
  
  if (length(gse) > 0) {
    gse_data <- gse[[1]]
    
    # Get expression matrix
    expr_matrix <- exprs(gse_data)
    cat("Expression matrix dimensions:", dim(expr_matrix), "\n")
    
    # Get sample names
    sample_names <- colnames(expr_matrix)
    
    # Determine data type based on expression values
    if (all(expr_matrix == floor(expr_matrix))) {
      data_type <- "rnaseq"
    } else {
      data_type <- "microarray"
    }
    
    return(list(expr_matrix = expr_matrix, sample_names = sample_names, 
                gse_data = gse_data, data_type = data_type))
  }
  
  return(NULL)
}

# Function to get proper sample metadata from GEO
get_sample_metadata <- function(gse_data = NULL) {
  cat("Downloading sample metadata from GEO...\n")
  
  if (is.null(gse_data)) {
    # Download the series matrix if not provided
    gse <- getGEO("GSE112676", GSEMatrix = TRUE)
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
    potential_group_cols <- c("diagnosis", "disease", "characteristics", "group")
    
    for (col in colnames(pdata)) {
      if (any(grepl(paste(potential_group_cols, collapse = "|"), tolower(col)))) {
        cat("\nPotential group column:", col, "\n")
        unique_vals <- unique(pdata[, col])
        if (length(unique_vals) <= 20) {
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

# Function to create proper sample groups - FIXED VERSION
create_proper_sample_groups <- function(sample_names, geo_metadata) {
  cat("Creating sample groups based on GEO metadata...\n")
  
  # Initialize groups
  groups <- rep(NA, length(sample_names))
  
  if (!is.null(geo_metadata)) {
    # Try to match with GEO metadata
    matched <- match(sample_names, rownames(geo_metadata))
    
    # Look for group information in various columns - PRIORITIZE DIAGNOSIS
    for (i in 1:length(sample_names)) {
      if (!is.na(matched[i])) {
        sample_meta <- geo_metadata[matched[i], ]
        
        # FIRST PRIORITY: Look for diagnosis information
        diagnosis_cols <- grep("diagnosis", colnames(sample_meta), ignore.case = TRUE, value = TRUE)
        for (col in diagnosis_cols) {
          if (!is.na(groups[i])) break
          
          value <- as.character(sample_meta[, col])
          if (!is.na(value) && value != "") {
            if (grepl("control|normal|healthy|non.neuronal|con", tolower(value))) {
              groups[i] <- "Control"
            } else if (grepl("als|motor.neuron|amyotrophic", tolower(value))) {
              groups[i] <- "ALS"
            }
          }
        }
        
        # SECOND PRIORITY: Look in characteristics columns
        if (is.na(groups[i])) {
          chars_cols <- grep("characteristics", colnames(sample_meta), ignore.case = TRUE, value = TRUE)
          for (col in chars_cols) {
            if (!is.na(groups[i])) break
            
            value <- as.character(sample_meta[, col])
            if (!is.na(value) && value != "") {
              if (grepl("control|normal|healthy|non.neuronal|con", tolower(value))) {
                groups[i] <- "Control"
              } else if (grepl("als|motor.neuron|amyotrophic", tolower(value))) {
                groups[i] <- "ALS"
              } else if (grepl("c9orf72|c9", tolower(value))) {
                groups[i] <- "C9orf72_ALS"
              } else if (grepl("sod1", tolower(value))) {
                groups[i] <- "SOD1_ALS"
              }
            }
          }
        }
        
        # THIRD PRIORITY: Look in source name
        if (is.na(groups[i]) && "source_name_ch1" %in% colnames(sample_meta)) {
          value <- as.character(sample_meta$source_name_ch1)
          if (grepl("control|normal|healthy|non.neuronal", tolower(value))) {
            groups[i] <- "Control"
          } else if (grepl("als|motor.neuron|amyotrophic", tolower(value))) {
            groups[i] <- "ALS"
          }
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
  valid_samples <- !is.na(groups) & groups %in% c("Control", "ALS")
  
  if (sum(valid_samples) < 10) {
    cat("Not enough samples with valid groups for DE analysis.\n")
    cat("Control samples:", sum(groups == "Control", na.rm = TRUE), "\n")
    cat("ALS samples:", sum(groups == "ALS", na.rm = TRUE), "\n")
    return(NULL)
  }
  
  expr_matrix_filtered <- expr_matrix[, valid_samples]
  groups_filtered <- groups[valid_samples]
  
  cat("Using", sum(valid_samples), "samples for DE analysis\n")
  cat("Group distribution - Control:", sum(groups_filtered == "Control"), 
      "ALS:", sum(groups_filtered == "ALS"), "\n")
  
  if (data_type == "rnaseq") {
    # Use edgeR for RNA-seq data
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
    
    # Make contrasts (ALS vs Control)
    contrast_matrix <- makeContrasts(
      contrasts = "ALS-Control",
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
cat("Starting GSE112676 analysis...\n")

# Get expression data from series matrix
series_data <- get_series_matrix_data()

if (!is.null(series_data)) {
  expr_matrix <- series_data$expr_matrix
  sample_names <- series_data$sample_names
  gse_data <- series_data$gse_data
  data_type <- series_data$data_type
  cat("Successfully obtained expression data from series matrix\n")
} else {
  stop("Could not obtain expression data for GSE112676")
}

# Get sample metadata
geo_metadata <- get_sample_metadata(gse_data)

# Create proper sample groups - FIXED VERSION
groups <- create_proper_sample_groups(sample_names, geo_metadata)

# Save the processed expression matrix
write.csv(expr_matrix, "GSE112676_expression_matrix.csv")
cat("Expression matrix saved to GSE112676_expression_matrix.csv\n")

# Save sample metadata with groups
sample_metadata <- data.frame(
  Sample = sample_names,
  Group = groups
)
write.csv(sample_metadata, "GSE112676_sample_metadata.csv")
cat("Sample metadata saved to GSE112676_sample_metadata.csv\n")

# Perform differential expression analysis
if (sum(groups %in% c("Control", "ALS"), na.rm = TRUE) >= 10) {
  de_results <- run_DE_analysis(expr_matrix, groups, data_type)
  
  if (!is.null(de_results)) {
    # Save results
    write.csv(de_results, "GSE112676_DE_results.csv")
    cat("Differential expression results saved to GSE112676_DE_results.csv\n")
    
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
  cat("Control samples:", sum(groups == "Control", na.rm = TRUE), "\n")
  cat("ALS samples:", sum(groups == "ALS", na.rm = TRUE), "\n")
  cat("Please check the sample_metadata.csv file and manually assign groups.\n")
}

cat("Analysis complete!\n")