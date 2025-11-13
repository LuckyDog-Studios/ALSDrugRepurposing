# Add diagnostic function to check each dataset
diagnose_datasets <- function(file_paths, dataset_names) {
  cat("\n=== DATASET DIAGNOSIS ===\n")
  
  for (i in 1:length(file_paths)) {
    cat("\n--- Diagnosing:", dataset_names[i], "---\n")
    
    # Read raw data
    raw_data <- read.csv(file_paths[i], stringsAsFactors = FALSE, check.names = FALSE)
    cat("  Raw data dimensions:", dim(raw_data), "\n")
    cat("  Column names:", paste(colnames(raw_data), collapse = ", "), "\n")
    
    # Check for NA values
    cat("  NA counts in key columns:\n")
    if ("logFC" %in% colnames(raw_data)) cat("    logFC:", sum(is.na(raw_data$logFC)), "\n")
    if ("P.Value" %in% colnames(raw_data)) cat("    P.Value:", sum(is.na(raw_data$P.Value)), "\n")
    if ("adj.P.Val" %in% colnames(raw_data)) cat("    adj.P.Val:", sum(is.na(raw_data$adj.P.Val)), "\n")
    
    # Check value ranges
    if ("logFC" %in% colnames(raw_data)) {
      cat("  logFC range:", range(raw_data$logFC, na.rm = TRUE), "\n")
    }
    if ("P.Value" %in% colnames(raw_data)) {
      cat("  P.Value range:", range(raw_data$P.Value, na.rm = TRUE), "\n")
    }
    if ("adj.P.Val" %in% colnames(raw_data)) {
      cat("  adj.P.Val range:", range(raw_data$adj.P.Val, na.rm = TRUE), "\n")
    }
    
    # Test filtering with current thresholds
    if (all(c("logFC", "adj.P.Val") %in% colnames(raw_data))) {
      sig_genes <- raw_data[abs(raw_data$logFC) >= 0.6 & raw_data$adj.P.Val <= 0.05, ]
      cat("  Genes passing current filters (logFC>=0.6, adj.P.Val<=0.05):", nrow(sig_genes), "\n")
      
      # Test with relaxed filters
      sig_relaxed <- raw_data[abs(raw_data$logFC) >= 0.5 & raw_data$adj.P.Val <= 0.1, ]
      cat("  Genes with relaxed filters (logFC>=0.5, adj.P.Val<=0.1):", nrow(sig_relaxed), "\n")
    }
  }
}


file_paths <- c(
  "datasets/processed/GSE112680/GSE112680_annotated.csv",
  "datasets/processed/GSE112676/GSE112676_annotated.csv", 
  "datasets/processed/GSE76220/GSE76220_properly_annotated.csv",
  "datasets/processed/GSE56500/results_sALS_genes.csv",
  "datasets/processed/GSE56500/results_C9_genes.csv"
)
dataset_names <- c("GSE112680", "GSE112676", "GSE76220", "GSE56500_sALS", "GSE56500_C9")

# Run diagnosis
diagnose_datasets(file_paths, dataset_names)