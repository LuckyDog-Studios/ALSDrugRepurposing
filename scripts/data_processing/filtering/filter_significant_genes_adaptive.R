# Load required libraries
library(dplyr)
library(tidyr)
library(ggplot2)
library(VennDiagram)
library(UpSetR)
library(here)
source(here("scripts","data_processing", "visualization", "comprehensive_visualizations.R"))

# Enhanced format detection function
detect_and_standardize_format <- function(data, dataset_name) {
  cat("Detecting format for:", dataset_name, "\n")
  
  # Create a standardized data frame
  standardized_data <- data
  
  # Convert column names to lowercase for easier matching
  colnames(standardized_data) <- tolower(colnames(standardized_data))
  
  # Map different column name variations to standard names
  # Check for logFC columns
  if ("logfc" %in% colnames(standardized_data)) {
    standardized_data$logFC <- standardized_data$logfc
  } else {
    fc_cols <- grep("logfc|log2fc|fc|foldchange", colnames(standardized_data), ignore.case = TRUE, value = TRUE)
    if (length(fc_cols) > 0) {
      standardized_data$logFC <- standardized_data[[fc_cols[1]]]
      cat("  Using column for logFC:", fc_cols[1], "\n")
    }
  }
  
  # Check for P.Value columns
  if ("p.value" %in% colnames(standardized_data)) {
    standardized_data$P.Value <- standardized_data$p.value
  } else {
    pval_cols <- grep("p.value|pval|p_value", colnames(standardized_data), ignore.case = TRUE, value = TRUE)
    if (length(pval_cols) > 0) {
      standardized_data$P.Value <- standardized_data[[pval_cols[1]]]
      cat("  Using column for P.Value:", pval_cols[1], "\n")
    }
  }
  
  # Check for adj.P.Val columns
  if ("adj.p.val" %in% colnames(standardized_data)) {
    standardized_data$adj.P.Val <- standardized_data$adj.p.val
  } else {
    adjp_cols <- grep("adj.p.val|fdr|qval|padj", colnames(standardized_data), ignore.case = TRUE, value = TRUE)
    if (length(adjp_cols) > 0) {
      standardized_data$adj.P.Val <- standardized_data[[adjp_cols[1]]]
      cat("  Using column for adj.P.Val:", adjp_cols[1], "\n")
    }
  }
  
  # Gene symbol columns
  if ("hgnc_symbol" %in% colnames(standardized_data)) {
    standardized_data$GeneSymbol <- standardized_data$hgnc_symbol
  } else if ("genesymbol" %in% colnames(standardized_data)) {
    standardized_data$GeneSymbol <- standardized_data$genesymbol
  } else if ("gene_symbol" %in% colnames(standardized_data)) {
    standardized_data$GeneSymbol <- standardized_data$gene_symbol
  } else {
    symbol_cols <- grep("genesymbol|hgnc_symbol|symbol|gene_name|genename", colnames(standardized_data), ignore.case = TRUE, value = TRUE)
    if (length(symbol_cols) > 0) {
      standardized_data$GeneSymbol <- standardized_data[[symbol_cols[1]]]
      cat("  Using column for GeneSymbol:", symbol_cols[1], "\n")
    }
  }
  
  # Validate that we have the required columns
  required_cols <- c("logFC", "P.Value", "adj.P.Val", "GeneSymbol")
  missing_cols <- setdiff(required_cols, colnames(standardized_data))
  
  if (length(missing_cols) > 0) {
    cat("  WARNING: Missing columns:", paste(missing_cols, collapse = ", "), "\n")
    cat("  Available columns:", paste(colnames(standardized_data), collapse = ", "), "\n")
  }
  
  return(standardized_data)
}

# Function to read and preprocess dataset
read_expression_data <- function(file_path, dataset_name) {
  cat("Reading dataset:", dataset_name, "\n")
  
  # Read the CSV file
  data <- read.csv(file_path, stringsAsFactors = FALSE, check.names = FALSE)
  
  # Detect format and standardize column names
  data <- detect_and_standardize_format(data, dataset_name)
  
  # Required columns check
  required_cols <- c("logFC", "P.Value", "adj.P.Val", "GeneSymbol")
  if (!all(required_cols %in% colnames(data))) {
    stop(paste("Missing required columns in", dataset_name, 
               ". Required:", paste(required_cols, collapse = ", "),
               "\nAvailable:", paste(colnames(data), collapse = ", ")))
  }
  
  # Clean GeneSymbol column (remove ///, NA, and empty values)
  data$GeneSymbol_clean <- sapply(data$GeneSymbol, function(x) {
    if (is.na(x) || x == "" || x == "NA") {
      return(NA)
    }
    # Extract first symbol if multiple separated by " /// "
    clean_symbol <- strsplit(as.character(x), " /// ")[[1]][1]
    return(clean_symbol)
  })
  
  # Remove rows with NA or empty GeneSymbols
  data <- data[!is.na(data$GeneSymbol_clean) & data$GeneSymbol_clean != "", ]
  
  # Add dataset identifier
  data$Dataset <- dataset_name
  
  cat("  Total rows:", nrow(data), "\n")
  cat("  Unique genes:", length(unique(data$GeneSymbol_clean)), "\n")
  cat("  logFC range:", round(range(data$logFC, na.rm = TRUE), 3), "\n")
  cat("  adj.P.Val range:", round(range(data$adj.P.Val, na.rm = TRUE), 3), "\n")
  
  return(data)
}

# Enhanced filtering function with dataset-specific thresholds
filter_significant_genes <- function(data, dataset_name,
                                    pval_threshold = 0.05, 
                                    adj_pval_threshold = 0.05,
                                    logFC_threshold = 0.5,
                                    method = "adj.P.Val") {
  
  cat("Filtering", dataset_name, "with thresholds:\n")
  cat("  P-value:", pval_threshold, "\n")
  cat("  adj.P.Val:", adj_pval_threshold, "\n") 
  cat("  logFC:", logFC_threshold, "\n")
  cat("  Method:", method, "\n")
  
  
  if (method == "adj.P.Val") {
    sig_adjp <- sum(data$adj.P.Val <= adj_pval_threshold, na.rm = TRUE)
    sig_fc <- sum(abs(data$logFC) >= logFC_threshold, na.rm = TRUE)
    sig_both <- sum(data$adj.P.Val <= adj_pval_threshold & 
                     abs(data$logFC) >= logFC_threshold, na.rm = TRUE)
    
    cat("  Genes passing adj.P.Val threshold:", sig_adjp, "\n")
    cat("  Genes passing logFC threshold:", sig_fc, "\n")
    cat("  Genes passing both:", sig_both, "\n")
    
    significant_data <- data %>%
      filter(adj.P.Val <= adj_pval_threshold & abs(logFC) >= logFC_threshold)
    
  } else if (method == "P.Value") {
    sig_p <- sum(data$P.Value <= pval_threshold, na.rm = TRUE)
    sig_fc <- sum(abs(data$logFC) >= logFC_threshold, na.rm = TRUE)
    sig_both <- sum(data$P.Value <= pval_threshold & 
                     abs(data$logFC) >= logFC_threshold, na.rm = TRUE)
    
    cat("  Genes passing P.Value threshold:", sig_p, "\n")
    cat("  Genes passing logFC threshold:", sig_fc, "\n")
    cat("  Genes passing both:", sig_both, "\n")
    
    significant_data <- data %>%
      filter(P.Value <= pval_threshold & abs(logFC) >= logFC_threshold)
  }
  
  # For genes with multiple probes, keep the most significant one
  if (nrow(significant_data) > 0) {
    significant_genes <- significant_data %>%
      group_by(GeneSymbol_clean) %>%
      slice_min(order_by = P.Value, n = 1) %>%
      ungroup()
  } else {
    significant_genes <- significant_data
  }
  
  cat("  Final significant genes:", nrow(significant_genes), "\n")
  
  return(significant_genes)
}

# Function to create consensus gene list
find_consensus_genes <- function(significant_lists, min_datasets = 2) {
  cat("\nFinding consensus genes...\n")
  cat("  Minimum datasets required:", min_datasets, "\n")
  
  # Remove empty datasets
  non_empty_lists <- significant_lists[sapply(significant_lists, nrow) > 0]
  
  if (length(non_empty_lists) == 0) {
    cat("  WARNING: No datasets contain significant genes!\n")
    return(list(
      consensus_genes = character(0),
      presence_matrix = matrix(0, nrow = 0, ncol = 0),
      gene_counts = integer(0)
    ))
  }
  
  # Get all unique genes across all datasets
  all_genes <- unique(unlist(lapply(non_empty_lists, function(x) unique(x$GeneSymbol_clean))))
  
  # Create presence/absence matrix
  presence_matrix <- matrix(0, nrow = length(all_genes), ncol = length(significant_lists))
  rownames(presence_matrix) <- all_genes
  colnames(presence_matrix) <- names(significant_lists)
  
  # Fill the matrix
  for (i in 1:length(significant_lists)) {
    dataset_genes <- unique(significant_lists[[i]]$GeneSymbol_clean)
    presence_matrix[rownames(presence_matrix) %in% dataset_genes, i] <- 1
  }
  
  # Calculate how many datasets each gene appears in
  gene_counts <- rowSums(presence_matrix)
  
  # Get consensus genes
  consensus_genes <- names(gene_counts[gene_counts >= min_datasets])
  
  cat("  Total unique genes across all datasets:", length(all_genes), "\n")
  cat("  Consensus genes (in >=", min_datasets, "datasets):", length(consensus_genes), "\n")
  
  # Print distribution
  cat("  Distribution of genes by dataset count:\n")
  print(table(gene_counts))
  
  return(list(
    consensus_genes = consensus_genes,
    presence_matrix = presence_matrix,
    gene_counts = gene_counts
  ))
}

# Function to save results
save_results <- function(significant_lists, consensus_results, output_dir = "results") {
  
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Save individual significant gene lists
  for (ds_name in names(significant_lists)) {
    if (nrow(significant_lists[[ds_name]]) > 0) {
      write.csv(significant_lists[[ds_name]],
                file.path(output_dir, paste0("significant_genes_", ds_name, ".csv")),
                row.names = FALSE)
    }
  }
  
  # Save consensus gene list
  if (length(consensus_results$consensus_genes) > 0) {
    consensus_df <- data.frame(
      GeneSymbol = consensus_results$consensus_genes,
      DatasetCount = consensus_results$gene_counts[consensus_results$consensus_genes]
    )
    write.csv(consensus_df,
              file.path(output_dir, "consensus_genes.csv"),
              row.names = FALSE)
  }
  
  # Save summary statistics
  summary_stats <- data.frame(
    Dataset = names(significant_lists),
    Significant_Genes = sapply(significant_lists, nrow),
    Unique_Genes = sapply(significant_lists, function(x) length(unique(x$GeneSymbol_clean)))
  )
  write.csv(summary_stats,
            file.path(output_dir, "summary_statistics.csv"),
            row.names = FALSE)
  
  cat("Results saved to:", output_dir, "\n")
}

# MAIN ANALYSIS PIPELINE with adaptive thresholds
run_analysis_adaptive <- function(file_paths, dataset_names, 
                                 output_dir = "results", min_datasets=2) {
  
  cat("Starting adaptive significant gene analysis pipeline...\n")
  cat("=====================================================\n")
  
  # 1. Read all datasets
  all_data <- list()
  for (i in 1:length(file_paths)) {
    all_data[[dataset_names[i]]] <- read_expression_data(file_paths[i], dataset_names[i])
  }
  
  # 2. Filter significant genes with adaptive thresholds
  significant_lists <- list()
  
  for (ds_name in dataset_names) {
    cat("\n--- Processing", ds_name, "---\n")
    
    # Set adaptive thresholds based on dataset characteristics
    if (ds_name == "GSE112680") {
      # Small effect sizes, use lenient thresholds
      sig_genes <- filter_significant_genes(
        all_data[[ds_name]], ds_name,
        pval_threshold = 0.05,
        adj_pval_threshold = 0.01,
        logFC_threshold = 0.3,
        method = "adj.P.Val"
      )
    } else if (ds_name == "GSE76220") {
      sig_genes <- filter_significant_genes(
        all_data[[ds_name]], ds_name,
        pval_threshold = 0.01,     # Stricter p-value
        adj_pval_threshold = 0.05,
        logFC_threshold = 0.6,
        method = "P.Value"        
      )
    } else {
      # Default thresholds for other datasets
      sig_genes <- filter_significant_genes(
        all_data[[ds_name]], ds_name,
        pval_threshold = 0.05,
        adj_pval_threshold = 0.05,
        logFC_threshold = 0.5,     # Moderate FC threshold
        method = "adj.P.Val"
      )
    }
    
    significant_lists[[ds_name]] <- sig_genes
  }
  
  # 3. Find consensus genes with adaptive minimum
  datasets_with_genes <- sum(sapply(significant_lists, nrow) > 0)
  
  cat("\nAdaptive consensus settings:\n")
  cat("  Datasets with significant genes:", datasets_with_genes, "\n")
  cat("  Minimum datasets for consensus:", min_datasets, "\n")
  
  consensus_results <- find_consensus_genes(significant_lists, min_datasets = min_datasets)
  
  # 4. Create visualizations
  # Define your dataset-specific thresholds
  dataset_thresholds <- list(
  GSE112680 = list(
    adj_pval_threshold = 0.01,
    logFC_threshold = 0.3
  ),
  GSE112676 = list(
    adj_pval_threshold = 0.05,
    logFC_threshold = 0.5
  ),
  GSE76220 = list(
    adj_pval_threshold = 0.05,
    logFC_threshold = 0.6
  ),
  GSE56500_sALS = list(
    adj_pval_threshold = 0.05,
    logFC_threshold = 0.5
  ),
  GSE56500_C9 = list(
    adj_pval_threshold = 0.05,
    logFC_threshold = 0.5
  )
)

# Use in your function call
create_comprehensive_visualizations(
  significant_lists = results$significant_lists,
  consensus_results = results$consensus_results,
  all_data = all_data,
  output_dir = "visualization",
  dataset_thresholds = dataset_thresholds
)
  
  # 5. Save results
  save_results(significant_lists, consensus_results, output_dir)
  
  # 6. Print final summary
  cat("\n=== ANALYSIS COMPLETE ===\n")
  cat("Final consensus gene count:", length(consensus_results$consensus_genes), "\n")
  
  if (length(consensus_results$consensus_genes) > 0) {
    cat("Consensus genes for PPI network:\n")
    print(consensus_results$consensus_genes)
  } else {
    cat("\nNo consensus genes found. Individual dataset results saved.\n")
    cat("Consider using individual dataset results or further relaxing thresholds.\n")
  }
  
  return(list(
    significant_lists = significant_lists,
    consensus_genes = consensus_results$consensus_genes,
    consensus_results = consensus_results
  ))
}

# Define file paths and dataset names
file_paths <- c(
  "datasets/processed/GSE112680/GSE112680_annotated.csv",
  "datasets/processed/GSE112676/GSE112676_annotated.csv", 
  "datasets/processed/GSE76220/GSE76220_properly_annotated.csv",
  "datasets/processed/GSE56500/results_sALS_genes.csv",
  "datasets/processed/GSE56500/results_C9_genes.csv"
)
dataset_names <- c("GSE112680", "GSE112676", "GSE76220", "GSE56500_sALS", "GSE56500_C9")

# Run the adaptive analysis
results <- run_analysis_adaptive(
  file_paths = file_paths,
  dataset_names = dataset_names,
  output_dir = "datasets/filtered",
  min_datasets = 1
)

# Print summary
cat("\n=== FINAL SUMMARY ===\n")
for (ds_name in dataset_names) {
  cat(sprintf("%-15s: %4d significant genes\n", ds_name, nrow(results$significant_lists[[ds_name]])))
}
cat(sprintf("%-15s: %4d consensus genes\n", "CONSENSUS", length(results$consensus_genes)))