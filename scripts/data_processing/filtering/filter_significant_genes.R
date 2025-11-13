# Load required libraries
library(dplyr)
library(tidyr)
library(ggplot2)
library(VennDiagram)
library(UpSetR)

#Detects the proper column names/format
detect_and_standardize_format <- function(data, dataset_name) {
  cat("Detecting format for:", dataset_name, "\n")
  
  # Create a standardized data frame
  standardized_data <- data
  
  # Check column names and standardize
  colnames(standardized_data) <- tolower(colnames(standardized_data))
  
  # Map different column name variations to standard names
  if ("logfc" %in% colnames(standardized_data)) {
    # Already has standard names (like your second example)
    cat("  Detected: Standard format\n")
    standardized_data$logFC <- standardized_data$logfc
    standardized_data$P.Value <- standardized_data$p.value
    standardized_data$adj.P.Val <- standardized_data$adj.p.val
    
    # Handle gene symbol columns
    if ("hgnc_symbol" %in% colnames(standardized_data)) {
      standardized_data$GeneSymbol <- standardized_data$hgnc_symbol
    } else if ("genesymbol" %in% colnames(standardized_data)) {
      standardized_data$GeneSymbol <- standardized_data$genesymbol
    }
    
  } else if ("logfc" %in% colnames(standardized_data)) {
    # First format (your original example)
    cat("  Detected: Original format\n")
    standardized_data$logFC <- standardized_data$logfc
    standardized_data$P.Value <- standardized_data$p.value
    standardized_data$adj.P.Val <- standardized_data$adj.p.val
    standardized_data$GeneSymbol <- standardized_data$genesymbol
    
  } else {
    cat("  Warning: Unknown format, attempting to use as-is\n")
    # Try to use whatever is there
    if (!"logfc" %in% colnames(standardized_data)) {
      # Find potential logFC column
      fc_cols <- grep("logfc|log2fc|fc|foldchange", colnames(standardized_data), ignore.case = TRUE, value = TRUE)
      if (length(fc_cols) > 0) {
        standardized_data$logFC <- standardized_data[[fc_cols[1]]]
        cat("  Using column for logFC:", fc_cols[1], "\n")
      }
    }
    
    if (!"p.value" %in% colnames(standardized_data)) {
      pval_cols <- grep("p.value|pval|p_value", colnames(standardized_data), ignore.case = TRUE, value = TRUE)
      if (length(pval_cols) > 0) {
        standardized_data$P.Value <- standardized_data[[pval_cols[1]]]
        cat("  Using column for P.Value:", pval_cols[1], "\n")
      }
    }
    
    if (!"adj.p.val" %in% colnames(standardized_data)) {
      adjp_cols <- grep("adj.p.val|fdr|qval|padj", colnames(standardized_data), ignore.case = TRUE, value = TRUE)
      if (length(adjp_cols) > 0) {
        standardized_data$adj.P.Val <- standardized_data[[adjp_cols[1]]]
        cat("  Using column for adj.P.Val:", adjp_cols[1], "\n")
      }
    }
    
    # Gene symbol columns
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
  cat("  Column names:", paste(colnames(data), collapse = ", "), "\n")
  
  return(data)
}

# Function to filter significant genes
filter_significant_genes <- function(data, 
                                    pval_threshold = 0.05, 
                                    adj_pval_threshold = 0.05,
                                    logFC_threshold = 1.0,
                                    method = "adj.P.Val") {
  
  cat("Filtering significant genes...\n")
  cat("  P-value threshold:", pval_threshold, "\n")
  cat("  adj.P.Val threshold:", adj_pval_threshold, "\n") 
  cat("  logFC threshold:", logFC_threshold, "\n")
  cat("  Method:", method, "\n")
  
  if (method == "adj.P.Val") {
    significant_data <- data %>%
      filter(adj.P.Val <= adj_pval_threshold & abs(logFC) >= logFC_threshold)
  } else if (method == "P.Value") {
    significant_data <- data %>%
      filter(P.Value <= pval_threshold & abs(logFC) >= logFC_threshold)
  } else {
    stop("Method must be either 'adj.P.Val' or 'P.Value'")
  }
  
  # For genes with multiple probes, keep the most significant one
  significant_genes <- significant_data %>%
    group_by(GeneSymbol_clean) %>%
    slice_min(order_by = P.Value, n = 1) %>%
    ungroup()
  
  cat("  Significant genes found:", nrow(significant_genes), "\n")
  cat("  Unique significant genes:", length(unique(significant_genes$GeneSymbol_clean)), "\n")
  
  return(significant_genes)
}

# Function to create consensus gene list
find_consensus_genes <- function(significant_lists, min_datasets = 2) {
  cat("\nFinding consensus genes...\n")
  cat("  Minimum datasets required:", min_datasets, "\n")
  
  # Get all unique genes across all datasets
  all_genes <- unique(unlist(lapply(significant_lists, function(x) unique(x$GeneSymbol_clean))))
  
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

# Function to create visualization
create_visualizations <- function(significant_lists, consensus_results, output_dir = "results") {
  
  # Create output directory
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # 1. Upset plot
  pdf(file.path(output_dir, "upset_plot.pdf"), width = 10, height = 6)
  upset_data <- as.data.frame(consensus_results$presence_matrix)
  upset(upset_data, nsets = ncol(upset_data), main.bar.color = "steelblue",
        sets.bar.color = "darkred", mb.ratio = c(0.6, 0.4),
        order.by = "freq")
  dev.off()
  
  # 2. Dataset overlap bar plot
  overlap_counts <- table(consensus_results$gene_counts)
  pdf(file.path(output_dir, "dataset_overlap.pdf"), width = 8, height = 6)
  barplot(overlap_counts, col = "lightblue", 
          main = "Number of Genes by Dataset Count",
          xlab = "Number of Datasets", ylab = "Number of Genes")
  dev.off()
  
  cat("Visualizations saved to:", output_dir, "\n")
}

# Function to save results
save_results <- function(significant_lists, consensus_results, output_dir = "results") {
  
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Save individual significant gene lists
  for (ds_name in names(significant_lists)) {
    write.csv(significant_lists[[ds_name]],
              file.path(output_dir, paste0("significant_genes_", ds_name, ".csv")),
              row.names = FALSE)
  }
  
  # Save consensus gene list
  consensus_df <- data.frame(
    GeneSymbol = consensus_results$consensus_genes,
    DatasetCount = consensus_results$gene_counts[consensus_results$consensus_genes]
  )
  write.csv(consensus_df,
            file.path(output_dir, "consensus_genes.csv"),
            row.names = FALSE)
  
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

# MAIN ANALYSIS PIPELINE
run_analysis <- function(file_paths, dataset_names, 
                        pval_threshold = 0.05,
                        adj_pval_threshold = 0.05, 
                        logFC_threshold = 1.0,
                        min_datasets = 2,
                        output_dir = "results") {
  
  cat("Starting significant gene analysis pipeline...\n")
  cat("=============================================\n")
  
  # 1. Read all datasets
  all_data <- list()
  for (i in 1:length(file_paths)) {
    all_data[[dataset_names[i]]] <- read_expression_data(file_paths[i], dataset_names[i])
  }
  
  # 2. Filter significant genes from each dataset
  significant_lists <- list()
  for (ds_name in dataset_names) {
    cat("\n--- Processing", ds_name, "---\n")
    significant_lists[[ds_name]] <- filter_significant_genes(
      all_data[[ds_name]],
      pval_threshold = pval_threshold,
      adj_pval_threshold = adj_pval_threshold,
      logFC_threshold = logFC_threshold,
      method = "adj.P.Val"
    )
  }
  
  # 3. Find consensus genes
  consensus_results <- find_consensus_genes(significant_lists, min_datasets = min_datasets)
  
  # 4. Create visualizations
  create_visualizations(significant_lists, consensus_results, output_dir)
  
  # 5. Save results
  save_results(significant_lists, consensus_results, output_dir)
  
  # 6. Print final summary
  cat("\n=== ANALYSIS COMPLETE ===\n")
  cat("Final consensus gene count:", length(consensus_results$consensus_genes), "\n")
  cat("These genes can be used for PPI network analysis.\n")
  
  return(list(
    significant_lists = significant_lists,
    consensus_genes = consensus_results$consensus_genes,
    consensus_results = consensus_results
  ))
}

# ====================
file_paths <- c(
  "datasets/processed/GSE112680/GSE112680_annotated.csv",
  "datasets/processed/GSE112676/GSE112676_annotated.csv", 
  "datasets/processed/GSE76220/GSE76220_properly_annotated.csv",
  "datasets/processed/GSE56500/results_sALS_genes.csv",
  "datasets/processed/GSE56500/results_C9_genes.csv"
)
dataset_names <- c("GSE112680", "GSE112676", "GSE76220", "GSE56500_sALS", "GSE56500_C9")

# Run the analysis
results <- run_analysis(
  file_paths = file_paths,
  dataset_names = dataset_names,
  pval_threshold = 0.05,        
  adj_pval_threshold = 0.05, 
  logFC_threshold = 0.6,        
  min_datasets = 1,             
  output_dir = "datasets/filtered"
)

# results
print("Consensus genes for PPI network:")
print(results$consensus_genes)
