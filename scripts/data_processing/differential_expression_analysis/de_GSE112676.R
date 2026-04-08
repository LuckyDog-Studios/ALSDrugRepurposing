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

# Function to create proper sample groups 
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
    return(NULL)
  }
  
  expr_matrix_filtered <- expr_matrix[, valid_samples]
  
  # Explicitly set factor levels so Control is ALWAYS the reference
  groups_filtered <- factor(groups[valid_samples], levels = c("Control", "ALS"))
  
  cat("Using", sum(valid_samples), "samples for DE analysis\n")
  
  if (data_type == "rnaseq") {
    dge <- DGEList(counts = expr_matrix_filtered, group = groups_filtered)
    keep <- filterByExpr(dge, min.count = 10, min.prop = 0.1)
    dge <- dge[keep, , keep.lib.sizes = FALSE]
    dge <- calcNormFactors(dge)
    design <- model.matrix(~groups_filtered)
    dge <- estimateDisp(dge, design)
    fit <- glmQLFit(dge, design)
    qlf <- glmQLFTest(fit, coef = 2)
    results <- topTags(qlf, n = Inf)$table
    
  } else {
    # Check for and apply log2 transformation if needed (Standard GEO2R heuristic)
    qx <- as.numeric(quantile(expr_matrix_filtered, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
    LogC <- (qx[5] > 100) || (qx[6]-qx[1] > 50 && qx[2] > 0)
    
    if (LogC) {
      cat("Data appears to be raw intensities. Applying log2 transformation...\n")
      expr_matrix_filtered[expr_matrix_filtered <= 0] <- NaN 
      expr_matrix_filtered <- log2(expr_matrix_filtered)
    } else {
      cat("Data appears to already be log2-transformed.\n")
    }

    # Create design matrix
    design <- model.matrix(~0 + groups_filtered)
    colnames(design) <- levels(groups_filtered)
    
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

# Create proper sample groups
groups <- create_proper_sample_groups(sample_names, geo_metadata)

# Save the processed expression matrix
write.csv(expr_matrix, "GSE112676_expression_matrix.csv")
cat("Expression matrix saved to GSE112676_expression_matrix.csv\n")

# Save sample metadata with groups
sample_metadata <- data.frame(
  Sample = sample_names,
  Group = groups
)
write.csv(sample_metadata, "GSE112676_sample_metadata.csv", row.names = FALSE)
cat("Sample metadata saved to GSE112676_sample_metadata.csv\n")

# Perform differential expression analysis
if (sum(groups %in% c("Control", "ALS"), na.rm = TRUE) >= 10) {
  de_results <- run_DE_analysis(expr_matrix, groups, data_type)
  
  if (!is.null(de_results)) {
    
    # ==========================================
    # NEW: Annotate Probe IDs to Gene Symbols
    # ==========================================
    cat("Annotating probe IDs to Gene Symbols...\n")
    
    # Extract feature data directly from the downloaded GEO object
    feature_data <- fData(gse_data)
    
    # Use regex to find the column containing gene symbols (usually "Symbol" or "ILMN_Gene")
    symbol_cols <- grep("symbol|ilmn_gene", colnames(feature_data), ignore.case = TRUE, value = TRUE)
    
    # Move rownames (Probe IDs) into their own column
    de_results$Probe_ID <- rownames(de_results)
    
    if (length(symbol_cols) > 0) {
      symbol_col <- symbol_cols[1] # Take the first matching column
      cat("Found gene symbols in column:", symbol_col, "\n")
      # Map the symbols using the rownames
      de_results$Gene_Symbol <- feature_data[rownames(de_results), symbol_col]
    } else {
      cat("Warning: Could not automatically identify a Gene Symbol column in the metadata.\n")
      de_results$Gene_Symbol <- NA
    }
    
    # Reorder columns so Probe_ID and Gene_Symbol are at the front for readability
    cols_to_move <- c("Probe_ID", "Gene_Symbol")
    other_cols <- setdiff(colnames(de_results), cols_to_move)
    de_results <- de_results[, c(cols_to_move, other_cols)]
    # ==========================================
    
    # Create directory safely if it doesn't exist
    if (!dir.exists("differential_expression")) {
      dir.create("differential_expression", recursive = TRUE)
    }
    
    # Save results (using row.names = FALSE since we made Probe_ID a proper column)
    write.csv(de_results, "differential_expression/GSE112676_DE_results.csv", row.names = FALSE)
    cat("Differential expression results saved to differential_expression/GSE112676_DE_results.csv\n")
    
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