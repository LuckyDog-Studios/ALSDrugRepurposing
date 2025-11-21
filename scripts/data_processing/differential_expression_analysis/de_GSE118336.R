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
get_series_matrix_data <- function(gse_id) {
  cat("Getting expression data from series matrix for", gse_id, "...\n")
  
  # Download the series matrix
  gse <- getGEO(gse_id, GSEMatrix = TRUE)
  
  if (length(gse) > 0) {
    gse_data <- gse[[1]]
    
    # Get expression matrix
    expr_matrix <- exprs(gse_data)
    cat("Expression matrix dimensions:", dim(expr_matrix), "\n")
    
    # Get sample names
    sample_names <- colnames(expr_matrix)
    
    # Determine data type based on expression values
    if (all(expr_matrix == floor(expr_matrix)) && max(expr_matrix) > 1000) {
      data_type <- "rnaseq"
    } else {
      data_type <- "microarray"
    }
    
    cat("Detected data type:", data_type, "\n")
    
    return(list(expr_matrix = expr_matrix, sample_names = sample_names, 
                gse_data = gse_data, data_type = data_type))
  }
  
  return(NULL)
}

# Function to get proper sample metadata from GEO
get_sample_metadata <- function(gse_data = NULL, gse_id = NULL) {
  cat("Downloading sample metadata from GEO...\n")
  
  if (is.null(gse_data)) {
    # Download the series matrix if not provided
    gse <- getGEO(gse_id, GSEMatrix = TRUE)
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
    potential_group_cols <- c("diagnosis", "disease", "characteristics", "group", "condition", "treatment", "genotype")
    
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

# Function to create proper sample groups for GSE118336 - FIXED VERSION
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
        
        # FIRST PRIORITY: Look for GENOTYPE information (key for GSE118336)
        genotype_cols <- grep("genotype", colnames(sample_meta), ignore.case = TRUE, value = TRUE)
        for (col in genotype_cols) {
          if (!is.na(groups[i])) break
          
          value <- as.character(sample_meta[, col])
          if (!is.na(value) && value != "") {
            # GSE118336 specific genotype patterns
            if (grepl("fuswt/wt|wild.type|control", tolower(value))) {
              groups[i] <- "Control"
            } else if (grepl("fush517d/h517d|h517d/h517d|mutant", tolower(value))) {
              groups[i] <- "FUS_H517D_Mutant"
            } else if (grepl("fuswt/h517d|heterozygous", tolower(value))) {
              groups[i] <- "FUS_Heterozygous"
            }
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

# Function to perform differential expression analysis for a specific comparison
run_DE_comparison <- function(expr_matrix, groups, group1, group2, data_type = "microarray", comparison_name = NULL) {
  if (is.null(comparison_name)) {
    comparison_name <- paste0(group1, "_vs_", group2)
  }
  
  cat("\n=== Performing DE analysis:", comparison_name, "===\n")
  
  # Select samples for this comparison
  valid_samples <- groups %in% c(group1, group2)
  expr_matrix_filtered <- expr_matrix[, valid_samples]
  groups_filtered <- groups[valid_samples]
  
  cat("Using", sum(valid_samples), "samples for DE analysis\n")
  cat("Group distribution:", table(groups_filtered), "\n")
  
  if (sum(valid_samples) < 6) {
    cat("Not enough samples for this comparison (need at least 3 per group)\n")
    return(NULL)
  }
  
  if (data_type == "rnaseq") {
    # Use edgeR for RNA-seq data
    dge <- DGEList(counts = expr_matrix_filtered, group = groups_filtered)
    keep <- filterByExpr(dge, min.count = 10, min.prop = 0.1)
    dge <- dge[keep, , keep.lib.sizes = FALSE]
    cat("Genes after filtering:", nrow(dge), "\n")
    
    dge <- calcNormFactors(dge)
    design <- model.matrix(~groups_filtered)
    dge <- estimateDisp(dge, design)
    fit <- glmQLFit(dge, design)
    qlf <- glmQLFTest(fit, coef = 2)
    results <- topTags(qlf, n = Inf)$table
    
  } else {
    # Use limma for microarray data
    design <- model.matrix(~0 + groups_filtered)
    colnames(design) <- gsub("groups_filtered", "", colnames(design))
    
    # Make contrasts (group1 vs group2)
    contrast_formula <- paste(group1, group2, sep = "-")
    cat("Contrast:", contrast_formula, "\n")
    
    contrast_matrix <- makeContrasts(
      contrasts = contrast_formula,
      levels = design
    )
    
    fit <- lmFit(expr_matrix_filtered, design)
    fit2 <- contrasts.fit(fit, contrast_matrix)
    fit2 <- eBayes(fit2)
    results <- topTable(fit2, number = Inf, adjust.method = "BH")
  }
  
  return(results)
}

# Function to perform comprehensive differential expression analysis
run_comprehensive_DE_analysis <- function(expr_matrix, groups, data_type = "microarray") {
  cat("Performing comprehensive differential expression analysis...\n")
  
  # Find all valid groups
  valid_groups <- unique(groups[!is.na(groups) & groups != "Unknown"])
  cat("Available groups:", paste(valid_groups, collapse = ", "), "\n")
  
  if (length(valid_groups) < 2) {
    cat("Need at least 2 groups for DE analysis\n")
    return(NULL)
  }
  
  # Generate all possible pairwise comparisons
  comparisons <- list()
  
  if ("Control" %in% valid_groups) {
    # Compare each mutant group against control
    mutant_groups <- setdiff(valid_groups, "Control")
    for (mutant in mutant_groups) {
      comparison_name <- paste0(mutant, "_vs_Control")
      comparisons[[comparison_name]] <- list(group1 = mutant, group2 = "Control")
    }
    
    # If there are multiple mutant groups, compare them against each other
    if (length(mutant_groups) >= 2) {
      for (i in 1:(length(mutant_groups)-1)) {
        for (j in (i+1):length(mutant_groups)) {
          comparison_name <- paste0(mutant_groups[i], "_vs_", mutant_groups[j])
          comparisons[[comparison_name]] <- list(group1 = mutant_groups[i], group2 = mutant_groups[j])
        }
      }
    }
  } else {
    # If no control, compare all groups pairwise
    for (i in 1:(length(valid_groups)-1)) {
      for (j in (i+1):length(valid_groups)) {
        comparison_name <- paste0(valid_groups[i], "_vs_", valid_groups[j])
        comparisons[[comparison_name]] <- list(group1 = valid_groups[i], group2 = valid_groups[j])
      }
    }
  }
  
  # Run all comparisons
  all_results <- list()
  
  for (comp_name in names(comparisons)) {
    group1 <- comparisons[[comp_name]]$group1
    group2 <- comparisons[[comp_name]]$group2
    
    cat("\n", strrep("-", 50), "\n")
    results <- run_DE_comparison(expr_matrix, groups, group1, group2, data_type, comp_name)
    
    if (!is.null(results)) {
      all_results[[comp_name]] <- results
      
      # Print summary for this comparison
      if ("adj.P.Val" %in% colnames(results)) {
        sig_genes <- sum(results$adj.P.Val < 0.05, na.rm = TRUE)
        cat("Significant genes (FDR < 0.05):", sig_genes, "\n")
        if (sig_genes > 0) {
          cat("Top gene:", rownames(results)[1], "logFC:", round(results$logFC[1], 3), 
              "FDR:", formatC(results$adj.P.Val[1], format = "e", digits = 2), "\n")
        }
      }
    }
  }
  
  return(all_results)
}

# Main analysis function for any GSE
run_GSE_analysis <- function(gse_id) {
  cat("Starting", gse_id, "analysis...\n")
  cat("==================================================\n")
  
  # Get expression data from series matrix
  series_data <- get_series_matrix_data(gse_id)
  
  if (!is.null(series_data)) {
    expr_matrix <- series_data$expr_matrix
    sample_names <- series_data$sample_names
    gse_data <- series_data$gse_data
    data_type <- series_data$data_type
    cat("Successfully obtained expression data from series matrix\n")
  } else {
    cat("Could not obtain expression data from series matrix\n")
    return(NULL)
  }
  
  # Get sample metadata
  geo_metadata <- get_sample_metadata(gse_data, gse_id)
  
  # Create proper sample groups
  groups <- create_proper_sample_groups(sample_names, geo_metadata)
  
  # Save the processed expression matrix
  output_file <- paste0(gse_id, "_expression_matrix.csv")
  write.csv(expr_matrix, output_file)
  cat("Expression matrix saved to", output_file, "\n")
  
  # Save sample metadata with groups
  metadata_file <- paste0(gse_id, "_sample_metadata.csv")
  sample_metadata <- data.frame(
    Sample = sample_names,
    Group = groups
  )
  write.csv(sample_metadata, metadata_file)
  cat("Sample metadata saved to", metadata_file, "\n")
  
  # Perform comprehensive differential expression analysis
  valid_groups <- unique(groups[!is.na(groups) & groups != "Unknown"])
  if (length(valid_groups) >= 2) {
    all_de_results <- run_comprehensive_DE_analysis(expr_matrix, groups, data_type)
    
    if (!is.null(all_de_results) && length(all_de_results) > 0) {
      # Save all results
      for (comp_name in names(all_de_results)) {
        results_file <- paste0(gse_id, "_DE_", comp_name, ".csv")
        write.csv(all_de_results[[comp_name]], results_file)
        cat("Results saved to", results_file, "\n")
      }
      
      # Print overall summary
      cat("\n")
      cat(strrep("=", 60), "\n")
      cat("OVERALL DIFFERENTIAL EXPRESSION SUMMARY\n")
      cat(strrep("=", 60), "\n")
      
      for (comp_name in names(all_de_results)) {
        results <- all_de_results[[comp_name]]
        if ("adj.P.Val" %in% colnames(results)) {
          sig_genes <- sum(results$adj.P.Val < 0.05, na.rm = TRUE)
          cat(comp_name, ":", sig_genes, "significant genes (FDR < 0.05)\n")
        }
      }
    }
  } else {
    cat("Not enough samples with proper group assignments for DE analysis.\n")
    cat("Available groups:", paste(valid_groups, collapse = ", "), "\n")
  }
  
  cat("\nAnalysis complete for", gse_id, "!\n")
  cat("==================================================\n\n")
  
  return(list(expr_matrix = expr_matrix, groups = groups, de_results = all_de_results))
}

# Run analysis for GSE118336
results <- run_GSE_analysis("GSE118336")

# Print final instructions
cat("=== INSTRUCTIONS ===\n")
cat("1. Multiple comparisons have been analyzed:\n")
cat("   - FUS_Heterozygous vs Control\n")
cat("   - FUS_H517D_Mutant vs Control\n")
cat("   - FUS_Heterozygous vs FUS_H517D_Mutant\n")
cat("2. Each comparison has been saved to a separate CSV file\n")
cat("3. Positive logFC means higher in the FIRST group of the comparison name\n")
cat("4. Negative logFC means higher in the SECOND group\n")