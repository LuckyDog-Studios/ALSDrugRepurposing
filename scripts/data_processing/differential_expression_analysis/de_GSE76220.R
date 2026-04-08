# ==============================================================================
# GSE76220 Differential Expression Analysis (Optimized for RPKM)
# ==============================================================================

# 1. Install/Load Required Packages
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
required_packages <- c("GEOquery", "limma", "Biobase", "dplyr", "readr", "stringr")
for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE)) {
    BiocManager::install(pkg)
    library(pkg, character.only = TRUE)
  }
}

# ------------------------------------------------------------------------------
# Function: Process Expression Data
# ------------------------------------------------------------------------------
process_GSE76220_expression <- function() {
  cat("--- Step 1: Processing GSE76220 RPKM Data ---\n")
  
  if (!file.exists("GSE76220/GSE76220_ALS_LCM_RPKM.txt.gz")) {
    cat("Downloading supplementary files from GEO...\n")
    getGEOSuppFiles("GSE76220", baseDir = ".")
  }
  
  rpkm_file <- "GSE76220/GSE76220_ALS_LCM_RPKM.txt.gz"
  expr_data <- read.delim(rpkm_file, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
  
  # Identify sample columns
  all_cols <- colnames(expr_data)
  annotation_cols <- c("Gene_symbol", "NCBI_mRNA", "NCBI_Protein", "Gene_Name", "42c")
  sample_columns <- setdiff(all_cols[1:22], annotation_cols) 
  
  gene_symbols <- expr_data$Gene_symbol
  expr_matrix <- as.matrix(expr_data[, sample_columns])
  expr_matrix <- apply(expr_matrix, 2, as.numeric)
  rownames(expr_matrix) <- gene_symbols
  
  # Handle Duplicate Gene Symbols
  if (any(duplicated(rownames(expr_matrix)))) {
    cat("Handling duplicate gene symbols by averaging...\n")
    expr_matrix <- aggregate(expr_matrix, list(Symbol = rownames(expr_matrix)), FUN = mean)
    rownames(expr_matrix) <- expr_matrix$Symbol
    expr_matrix$Symbol <- NULL
    expr_matrix <- as.matrix(expr_matrix)
  }
  
  cat("Cleaned matrix dimensions:", dim(expr_matrix), "\n")
  return(list(expr_matrix = expr_matrix, sample_names = sample_columns))
}

# ------------------------------------------------------------------------------
# Function: Define Sample Groups
# ------------------------------------------------------------------------------
get_groups <- function(sample_names) {
  control_samples <- c("10c", "65c", "78c", "39c", "67c", "76c", "44c", "88c")
  als_samples <- c("60a", "62a", "63a", "84a", "89a", "21a", "34a", "79a", "82a", "16a", "27a", "48a", "85a")
  
  groups <- ifelse(sample_names %in% control_samples, "Control", 
                   ifelse(sample_names %in% als_samples, "ALS", "Unknown"))
  
  # Validation check
  if(any(groups == "Unknown")) stop("Error: Some samples were not correctly assigned to groups.")
  
  cat("\nGroup distribution:\n")
  print(table(groups))
  return(factor(groups, levels = c("Control", "ALS")))
}

# ------------------------------------------------------------------------------
# Function: Limma Differential Expression
# ------------------------------------------------------------------------------
run_limma_analysis <- function(expr_matrix, groups) {
  cat("\n--- Step 2: Running Differential Expression (limma-trend) ---\n")
  
  # 1. Filter Low Expression (RPKM > 0.5 in at least 5 samples)
  keep <- rowSums(expr_matrix > 0.5) >= 5
  expr_filtered <- expr_matrix[keep, ]
  cat("Genes remaining after filtering:", nrow(expr_filtered), "\n")
  
  # 2. Log2 Transformation 
  # Using a +1 offset is standard for RPKM to handle zeros
  expr_transformed <- log2(expr_filtered + 1)
  
  # 3. Setup Design and Contrasts
  # We use the group names directly as column headers for clarity
  design <- model.matrix(~0 + groups)
  colnames(design) <- levels(groups)
  
  # Define the comparison: ALS minus Control
  # Positive logFC = Up in ALS; Negative logFC = Down in ALS
  contrast_matrix <- makeContrasts(ALS_vs_Control = ALS - Control, levels = design)
  
  # 4. Linear Modeling with Empirical Bayes Trend
  # 'trend = TRUE' is vital when using RPKM/CPM instead of voom(counts)
  fit <- lmFit(expr_transformed, design)
  fit2 <- contrasts.fit(fit, contrast_matrix)
  fit2 <- eBayes(fit2, trend = TRUE) 
  
  # 5. Extract Results
  results <- topTable(fit2, coef = "ALS_vs_Control", number = Inf, adjust.method = "BH")
  
  # Add Mean Expression (AveExpr) as a log2 value for plotting
  return(results)
}

# ------------------------------------------------------------------------------
# Main Execution
# ------------------------------------------------------------------------------
run_full_analysis <- function() {
  if(!dir.exists("differential_expression")) dir.create("differential_expression")
  
  data_load <- process_GSE76220_expression()
  groups <- get_groups(data_load$sample_names)
  
  de_results <- run_limma_analysis(data_load$expr_matrix, groups)
  
  # Ensure numeric columns are stored as double (Inf values from limma are
  # written as the string "Inf" by write.csv, making the column varchar)
  de_results$logFC     <- as.double(de_results$logFC)
  de_results$AveExpr   <- as.double(de_results$AveExpr)
  de_results$t         <- as.double(de_results$t)
  de_results$P.Value   <- as.double(de_results$P.Value)
  de_results$adj.P.Val <- as.double(de_results$adj.P.Val)
  de_results$B         <- as.double(de_results$B)

  # Replace any ±Inf (e.g. from zero-variance genes) with NA
  de_results[sapply(de_results, is.numeric)] <- lapply(
    de_results[sapply(de_results, is.numeric)],
    function(x) { x[is.infinite(x)] <- NA; x }
  )

  write.csv(de_results, "differential_expression/GSE76220_DE_Results.csv", na = "")
  
  cat("\n--- Final Summary ---\n")
  cat("Total genes tested:", nrow(de_results), "\n")
  cat("Significant genes (FDR < 0.05):", sum(de_results$adj.P.Val < 0.05, na.rm = TRUE), "\n")
  
  return(de_results)
}

final_results <- run_full_analysis()