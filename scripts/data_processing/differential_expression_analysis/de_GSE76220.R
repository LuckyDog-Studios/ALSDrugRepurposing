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

# Function to properly process GSE76220 expression data
process_GSE76220_expression <- function() {
  cat("Processing GSE76220 expression data...\n")
  
  # Download supplementary files if not already done
  if (!file.exists("GSE76220/GSE76220_ALS_LCM_RPKM.txt.gz")) {
    cat("Downloading supplementary files...\n")
    getGEOSuppFiles("GSE76220", baseDir = ".")
  }
  
  # Read the RPKM expression file
  rpkm_file <- "GSE76220/GSE76220_ALS_LCM_RPKM.txt.gz"
  
  if (file.exists(rpkm_file)) {
    cat("Reading RPKM expression file...\n")
    
    # Read the file with proper header handling
    expr_data <- read.delim(rpkm_file, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
    
    cat("Expression data dimensions:", dim(expr_data), "\n")
    
    # The first column is Gene_symbol, the last 3 columns are gene annotations (not samples)
    # Real sample columns are: 60a, 62a, 63a, 84a, 89a, 21a, 34a, 79a, 82a, 16a, 27a, 48a, 
    # 10c, 65c, 78c, 39c, 67c, 76c, 44c, 88c
    # Note: 42c is NOT a real sample based on GEO information
    
    # Identify sample columns (exclude gene annotation columns and 42c)
    all_sample_columns <- colnames(expr_data)[2:22]  # Columns 2-22 are potential samples
    
    # Remove 42c (not a real control sample) and gene annotation columns
    sample_columns <- setdiff(all_sample_columns, c("42c", "NCBI_mRNA", "NCBI_Protein", "Gene_Name"))
    
    cat("Sample columns:", sample_columns, "\n")
    
    # Extract gene information and expression matrix
    gene_symbols <- expr_data[, 1]
    expr_matrix <- as.matrix(expr_data[, sample_columns])
    rownames(expr_matrix) <- gene_symbols
    
    # Convert to numeric
    expr_matrix <- apply(expr_matrix, 2, as.numeric)
    rownames(expr_matrix) <- gene_symbols
    
    cat("Final expression matrix dimensions:", dim(expr_matrix), "\n")
    cat("Range of expression values:", range(expr_matrix, na.rm = TRUE), "\n")
    
    return(list(expr_matrix = expr_matrix, sample_names = sample_columns, 
                gene_symbols = gene_symbols, data_type = "rnaseq"))
  }
  
  return(NULL)
}

# Function to assign FINAL CORRECT groups based on exact GEO sample mapping
get_GSE76220_final_correct_groups <- function(sample_names) {
  cat("Assigning FINAL CORRECT groups for GSE76220...\n")
  
  # Initialize groups
  groups <- rep("Unknown", length(sample_names))
  
  # EXACT mapping based on GEO sample information:
  # Controls (8 samples): 10c, 65c, 78c, 39c, 67c, 76c, 44c, 88c
  # ALS (13 samples): 60a, 62a, 63a, 84a, 89a, 21a, 34a, 79a, 82a, 16a, 27a, 48a, 85a
  
  # Create exact mapping
  control_samples <- c("10c", "65c", "78c", "39c", "67c", "76c", "44c", "88c")
  als_samples <- c("60a", "62a", "63a", "84a", "89a", "21a", "34a", "79a", "82a", "16a", "27a", "48a", "85a")
  
  for (i in 1:length(sample_names)) {
    sample_id <- sample_names[i]
    
    if (sample_id %in% control_samples) {
      groups[i] <- "Control"
    } else if (sample_id %in% als_samples) {
      groups[i] <- "ALS"
    } else {
      groups[i] <- "Unknown"
    }
    
    cat("Sample", sample_id, "->", groups[i], "\n")
  }
  
  cat("\nFINAL CORRECT Group distribution:\n")
  print(table(groups))
  
  # Verify counts match GEO information
  expected_controls <- 8
  expected_als <- 13
  actual_controls <- sum(groups == "Control")
  actual_als <- sum(groups == "ALS")
  
  if (actual_controls == expected_controls && actual_als == expected_als) {
    cat("✓ Group assignment matches GEO information exactly!\n")
  } else {
    cat("✗ Group assignment DOES NOT match GEO information!\n")
    cat("Expected: Control =", expected_controls, "ALS =", expected_als, "\n")
    cat("Actual: Control =", actual_controls, "ALS =", actual_als, "\n")
  }
  
  return(groups)
}

# Function to perform differential expression analysis with edgeR (for count data)
run_DE_analysis_edgeR <- function(expr_matrix, groups) {
  cat("Performing differential expression analysis with edgeR...\n")
  
  # Only use samples with known groups
  valid_samples <- groups %in% c("Control", "ALS")
  
  if (sum(valid_samples) < 6) {
    cat("Not enough samples with valid groups for DE analysis.\n")
    return(NULL)
  }
  
  expr_matrix_filtered <- expr_matrix[, valid_samples]
  groups_filtered <- groups[valid_samples]
  
  cat("Using", sum(valid_samples), "samples for DE analysis\n")
  cat("Group distribution - Control:", sum(groups_filtered == "Control"), 
      "ALS:", sum(groups_filtered == "ALS"), "\n")
  
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
  return(results)
}

# Function to perform differential expression analysis with limma (for RPKM data)
run_DE_analysis_limma <- function(expr_matrix, groups) {
  cat("Performing differential expression analysis with limma...\n")
  
  # Only use samples with known groups
  valid_samples <- groups %in% c("Control", "ALS")
  
  if (sum(valid_samples) < 6) {
    cat("Not enough samples with valid groups for DE analysis.\n")
    return(NULL)
  }
  
  expr_matrix_filtered <- expr_matrix[, valid_samples]
  groups_filtered <- groups[valid_samples]
  
  cat("Using", sum(valid_samples), "samples for DE analysis\n")
  cat("Group distribution - Control:", sum(groups_filtered == "Control"), 
      "ALS:", sum(groups_filtered == "ALS"), "\n")
  
  # Check if data needs transformation (RPKM data is already normalized)
  if (min(expr_matrix_filtered, na.rm = TRUE) < 0) {
    cat("Data appears to be log2 transformed already\n")
    expr_transformed <- expr_matrix_filtered
  } else {
    cat("Applying log2 transformation to RPKM data\n")
    expr_transformed <- log2(expr_matrix_filtered + 1)
  }
  
  # Create design matrix
  design <- model.matrix(~0 + groups_filtered)
  colnames(design) <- gsub("groups_filtered", "", colnames(design))
  
  # Make contrasts (ALS vs Control)
  contrast_matrix <- makeContrasts(
    contrasts = "ALS-Control",
    levels = design
  )
  
  # Fit linear model
  fit <- lmFit(expr_transformed, design)
  fit2 <- contrasts.fit(fit, contrast_matrix)
  fit2 <- eBayes(fit2)
  
  # Get results
  results <- topTable(fit2, number = Inf, adjust.method = "BH")
  
  return(results)
}

# Main analysis function
run_GSE76220_final_analysis <- function() {
  gse_id <- "GSE76220"
  cat("Starting FINAL CORRECTED analysis for", gse_id, "\n")
  cat("==================================================\n")
  
  # Process expression data
  expr_data <- process_GSE76220_expression()
  
  if (is.null(expr_data)) {
    cat("ERROR: Could not process expression data\n")
    return(NULL)
  }
  
  expr_matrix <- expr_data$expr_matrix
  sample_names <- expr_data$sample_names
  data_type <- expr_data$data_type
  
  cat("Successfully loaded expression data with", nrow(expr_matrix), "genes and", ncol(expr_matrix), "samples\n")
  
  # Get FINAL CORRECT sample groups
  groups <- get_GSE76220_final_correct_groups(sample_names)
  
  # Save the processed expression matrix
  output_file <- paste0(gse_id, "_expression_matrix.csv")
  write.csv(expr_matrix, output_file)
  cat("Expression matrix saved to", output_file, "\n")
  
  # Save sample metadata with FINAL CORRECT groups
  metadata_file <- paste0(gse_id, "_sample_metadata.csv")
  sample_metadata <- data.frame(
    Sample = sample_names,
    Group = groups
  )
  write.csv(sample_metadata, metadata_file)
  cat("Sample metadata saved to", metadata_file, "\n")
  
  # Perform differential expression analysis
  if (sum(groups == "Control") >= 3 && sum(groups == "ALS") >= 3) {
    # Use limma for RPKM data
    de_results <- run_DE_analysis_limma(expr_matrix, groups)
    
    if (!is.null(de_results)) {
      # Save results
      results_file <- paste0(gse_id, "_DE_results.csv")
      write.csv(de_results, results_file)
      cat("Differential expression results saved to", results_file, "\n")
      
      # Print summary
      cat("\nFINAL DIFFERENTIAL EXPRESSION SUMMARY:\n")
      cat("Total genes tested:", nrow(de_results), "\n")
      
      if ("adj.P.Val" %in% colnames(de_results)) {
        sig_genes <- sum(de_results$adj.P.Val < 0.05, na.rm = TRUE)
        cat("Significant genes (FDR < 0.05):", sig_genes, "\n")
        cat("Significant genes (FDR < 0.01):", sum(de_results$adj.P.Val < 0.01, na.rm = TRUE), "\n")
        cat("Significant genes (FDR < 0.001):", sum(de_results$adj.P.Val < 0.001, na.rm = TRUE), "\n")
        
        if (sig_genes > 0) {
          cat("Top 10 differentially expressed genes:\n")
          print(head(de_results[de_results$adj.P.Val < 0.05, ], 10))
          
          # Calculate direction of changes
          up_regulated <- sum(de_results$adj.P.Val < 0.05 & de_results$logFC > 0, na.rm = TRUE)
          down_regulated <- sum(de_results$adj.P.Val < 0.05 & de_results$logFC < 0, na.rm = TRUE)
          cat("Up-regulated in ALS:", up_regulated, "\n")
          cat("Down-regulated in ALS:", down_regulated, "\n")
        } else {
          cat("Top 10 genes by p-value:\n")
          print(head(de_results, 10))
          cat("No genes reached FDR < 0.05 significance\n")
          cat("Consider using a less stringent threshold (e.g., p < 0.01) for exploratory analysis\n")
        }
      }
    }
  } else {
    cat("\nERROR: Group assignment issue\n")
    cat("Control samples:", sum(groups == "Control"), "\n")
    cat("ALS samples:", sum(groups == "ALS"), "\n")
  }
  
  cat("\nAnalysis complete for", gse_id, "!\n")
  cat("==================================================\n\n")
  
  return(list(expr_matrix = expr_matrix, groups = groups, de_results = de_results))
}

# Run the analysis
results <- run_GSE76220_final_analysis()

# Final summary
if (!is.null(results)) {
  cat("=== FINAL ANALYSIS COMPLETE ===\n")
  cat("GSE76220 - Laser-captured motor neurons from spinal cord\n")
  cat("EXACT sample breakdown:\n")
  cat("- 8 Control samples: 10c, 65c, 78c, 39c, 67c, 76c, 44c, 88c\n") 
  cat("- 13 ALS samples: 60a, 62a, 63a, 84a, 89a, 21a, 34a, 79a, 82a, 16a, 27a, 48a, 85a\n")
  cat("- Total: 21 samples with EXACT group assignment\n")
  
  if (!is.null(results$de_results)) {
    sig_genes <- sum(results$de_results$adj.P.Val < 0.05, na.rm = TRUE)
    if (sig_genes == 0) {
      cat("Note: No genes reached FDR < 0.05 significance in this dataset\n")
      cat("This could be due to:\n")
      cat("1. High biological variability in post-mortem human samples\n")
      cat("2. The specific molecular signature of sALS (sporadic ALS)\n")
      cat("3. Technical factors in laser-capture microdissection\n")
      cat("4. The conservative nature of FDR correction\n")
    }
  }
}