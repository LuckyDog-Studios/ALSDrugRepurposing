#!/usr/bin/env Rscript

# Script to perform comprehensive diagnostics on DE results
# Helps identify issues with statistical significance, effect sizes, and data quality

library(data.table)
library(dplyr)
library(ggplot2)
library(patchwork)

# Function to perform comprehensive diagnostics
diagnose_de_results <- function(de_file, output_dir = "./de_diagnostics") {
  
  # Create output directory
  if(!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  cat("========================================\n")
  cat("COMPREHENSIVE DE RESULTS DIAGNOSTICS\n")
  cat("========================================\n\n")
  
  # Load data
  cat("1. Loading data...\n")
  de_data <- fread(de_file)
  cat(paste0("   Total rows: ", nrow(de_data), "\n"))
  cat(paste0("   Total columns: ", ncol(de_data), "\n"))
  cat("   Columns:\n")
  print(colnames(de_data))
  
  # Check for required columns
  required_cols <- c("SYMBOL", "logFC", "PValue", "adj.P.Val")
  missing_cols <- setdiff(required_cols, colnames(de_data))
  if(length(missing_cols) > 0) {
    cat(paste0("\n⚠ WARNING: Missing required columns: ", paste(missing_cols, collapse = ", "), "\n"))
  } else {
    cat("\n✅ All required columns present\n")
  }
  
  # ========== BASIC STATISTICS ==========
  cat("\n2. BASIC STATISTICS:\n")
  cat("   ---------------------------------\n")
  
  # Unique genes
  unique_genes <- unique(de_data$SYMBOL)
  unique_genes <- unique_genes[!is.na(unique_genes) & unique_genes != ""]
  cat(paste0("   Unique genes: ", length(unique_genes), "\n"))
  
  # Dataset information
  if("n_datasets" %in% colnames(de_data)) {
    cat("\n   Dataset distribution:\n")
    dataset_dist <- de_data[, .(count = .N, unique_genes = length(unique(SYMBOL))), by = n_datasets][order(-n_datasets)]
    print(dataset_dist)
    
    # Save dataset distribution
    fwrite(dataset_dist, file.path(output_dir, "dataset_distribution.tsv"), sep = "\t")
  }
  
  if("Dataset" %in% colnames(de_data)) {
    dataset_counts <- de_data[, .N, by = Dataset][order(-N)]
    cat("\n   Rows per dataset:\n")
    print(dataset_counts)
  }
  
  # ========== STATISTICAL QUALITY CHECKS ==========
  cat("\n3. STATISTICAL QUALITY CHECKS:\n")
  cat("   ---------------------------------\n")
  
  # Check for NA values
  na_counts <- sapply(de_data[, ..required_cols], function(x) sum(is.na(x)))
  cat("   NA counts in key columns:\n")
  for(col in names(na_counts)) {
    cat(paste0("     ", col, ": ", na_counts[col], " (", 
               round(100 * na_counts[col]/nrow(de_data), 1), "%)\n"))
  }
  
  # Check P-value distribution
  cat("\n   P-value summary:\n")
  pval_summary <- summary(de_data$PValue)
  print(pval_summary)
  
  # Check adjusted P-value distribution
  cat("\n   Adjusted P-value summary:\n")
  adj_pval_summary <- summary(de_data$adj.P.Val)
  print(adj_pval_summary)
  
  # Check logFC distribution
  cat("\n   logFC summary:\n")
  logfc_summary <- summary(de_data$logFC)
  print(logfc_summary)
  
  # Count by significance thresholds
  cat("\n   Significance by different thresholds:\n")
  thresholds <- list(
    "adj.P.Val < 0.001" = sum(de_data$adj.P.Val < 0.001, na.rm = TRUE),
    "adj.P.Val < 0.01" = sum(de_data$adj.P.Val < 0.01, na.rm = TRUE),
    "adj.P.Val < 0.05" = sum(de_data$adj.P.Val < 0.05, na.rm = TRUE),
    "adj.P.Val < 0.1" = sum(de_data$adj.P.Val < 0.1, na.rm = TRUE),
    "adj.P.Val < 0.2" = sum(de_data$adj.P.Val < 0.2, na.rm = TRUE)
  )
  
  for(thresh in names(thresholds)) {
    cat(paste0("     ", thresh, ": ", thresholds[[thresh]], 
               " (", round(100 * thresholds[[thresh]]/nrow(de_data), 1), "%)\n"))
  }
  
  # Count by effect size thresholds
  cat("\n   Effect size distribution:\n")
  effect_sizes <- list(
    "|logFC| > 2" = sum(abs(de_data$logFC) > 2, na.rm = TRUE),
    "|logFC| > 1" = sum(abs(de_data$logFC) > 1, na.rm = TRUE),
    "|logFC| > 0.5" = sum(abs(de_data$logFC) > 0.5, na.rm = TRUE),
    "|logFC| > 0.2" = sum(abs(de_data$logFC) > 0.2, na.rm = TRUE)
  )
  
  for(effect in names(effect_sizes)) {
    cat(paste0("     ", effect, ": ", effect_sizes[[effect]], 
               " (", round(100 * effect_sizes[[effect]]/nrow(de_data), 1), "%)\n"))
  }
  
  # ========== GENE-LEVEL ANALYSIS ==========
  cat("\n4. GENE-LEVEL ANALYSIS:\n")
  cat("   ---------------------------------\n")
  
  # Create gene-level summary
  if(all(c("SYMBOL", "adj.P.Val", "logFC", "n_datasets") %in% colnames(de_data))) {
    gene_summary <- de_data[, .(
      n_observations = .N,
      min_adj_p = min(adj.P.Val, na.rm = TRUE),
      max_logfc = max(abs(logFC), na.rm = TRUE),
      mean_logfc = mean(logFC, na.rm = TRUE),
      sd_logfc = sd(logFC, na.rm = TRUE),
      n_datasets = first(n_datasets)
    ), by = SYMBOL][order(min_adj_p)]
    
    cat("   Top 20 most significant genes:\n")
    print(head(gene_summary, 20))
    
    # Save gene summary
    fwrite(gene_summary, file.path(output_dir, "gene_summary.tsv"), sep = "\t")
    
    # Count genes by number of datasets
    if("n_datasets" %in% colnames(gene_summary)) {
      cat("\n   Genes by dataset appearance:\n")
      gene_dataset_counts <- gene_summary[, .(count = .N), by = n_datasets][order(-n_datasets)]
      print(gene_dataset_counts)
    }
    
    # Most differentially expressed genes
    cat("\n   Top 20 genes by absolute logFC:\n")
    top_logfc <- gene_summary[order(-max_logfc), .(SYMBOL, max_logfc, min_adj_p, n_datasets)][1:20]
    print(top_logfc)
    
  } else {
    cat("   ⚠ Cannot perform gene-level analysis - missing required columns\n")
  }
  
  # ========== PLOTS ==========
  cat("\n5. GENERATING DIAGNOSTIC PLOTS...\n")
  
  tryCatch({
    # Create plots directory
    plots_dir <- file.path(output_dir, "plots")
    if(!dir.exists(plots_dir)) dir.create(plots_dir)
    
    # Plot 1: P-value histogram
    p1 <- ggplot(de_data, aes(x = PValue)) +
      geom_histogram(bins = 50, fill = "steelblue", alpha = 0.7) +
      labs(title = "Distribution of P-values",
           x = "P-value", y = "Count") +
      theme_minimal() +
      geom_vline(xintercept = 0.05, color = "red", linetype = "dashed")
    
    # Plot 2: Adjusted P-value histogram
    p2 <- ggplot(de_data, aes(x = adj.P.Val)) +
      geom_histogram(bins = 50, fill = "darkorange", alpha = 0.7) +
      labs(title = "Distribution of Adjusted P-values",
           x = "Adjusted P-value", y = "Count") +
      theme_minimal() +
      geom_vline(xintercept = 0.05, color = "red", linetype = "dashed")
    
    # Plot 3: logFC distribution
    p3 <- ggplot(de_data, aes(x = logFC)) +
      geom_histogram(bins = 100, fill = "darkgreen", alpha = 0.7) +
      labs(title = "Distribution of logFC values",
           x = "logFC", y = "Count") +
      theme_minimal() +
      geom_vline(xintercept = c(-1, 1), color = "red", linetype = "dashed")
    
    # Plot 4: Volcano plot (if we have enough data)
    if(nrow(de_data) > 0 && sum(!is.na(de_data$logFC) & !is.na(de_data$adj.P.Val)) > 0) {
      volcano_data <- de_data[, .(logFC, adj.P.Val, SYMBOL)]
      volcano_data <- volcano_data[!is.na(logFC) & !is.na(adj.P.Val)]
      
      # Add significance labels
      volcano_data$significance <- "Not significant"
      volcano_data$significance[volcano_data$adj.P.Val < 0.05 & abs(volcano_data$logFC) > 1] <- "Significant"
      
      p4 <- ggplot(volcano_data, aes(x = logFC, y = -log10(adj.P.Val), color = significance)) +
        geom_point(alpha = 0.5, size = 1) +
        scale_color_manual(values = c("Not significant" = "gray", "Significant" = "red")) +
        labs(title = "Volcano Plot",
             x = "logFC", y = "-log10(adj.P.Val)") +
        theme_minimal() +
        geom_vline(xintercept = c(-1, 1), linetype = "dashed", alpha = 0.5) +
        geom_hline(yintercept = -log10(0.05), linetype = "dashed", alpha = 0.5)
      
      # Save volcano plot
      ggsave(file.path(plots_dir, "volcano_plot.png"), p4, width = 8, height = 6, dpi = 300)
    }
    
    # Plot 5: Dataset distribution (if available)
    if("n_datasets" %in% colnames(de_data)) {
      dataset_plot_data <- de_data[, .N, by = n_datasets][order(n_datasets)]
      p5 <- ggplot(dataset_plot_data, aes(x = factor(n_datasets), y = N)) +
        geom_bar(stat = "identity", fill = "purple", alpha = 0.7) +
        labs(title = "Number of Observations by Dataset Count",
             x = "Number of Datasets", y = "Count") +
        theme_minimal()
      
      ggsave(file.path(plots_dir, "dataset_distribution.png"), p5, width = 8, height = 6, dpi = 300)
    }
    
    # Save basic plots
    ggsave(file.path(plots_dir, "pvalue_distribution.png"), p1, width = 8, height = 6, dpi = 300)
    ggsave(file.path(plots_dir, "adj_pvalue_distribution.png"), p2, width = 8, height = 6, dpi = 300)
    ggsave(file.path(plots_dir, "logfc_distribution.png"), p3, width = 8, height = 6, dpi = 300)
    
    # Combined plot
    combined_plot <- (p1 | p2) / p3
    ggsave(file.path(plots_dir, "combined_diagnostics.png"), combined_plot, 
           width = 12, height = 9, dpi = 300)
    
    cat("   ✅ Diagnostic plots saved to: ", plots_dir, "\n")
    
  }, error = function(e) {
    cat("   ⚠ Could not create plots: ", e$message, "\n")
  })
  
  # ========== QUALITY METRICS ==========
  cat("\n6. QUALITY METRICS:\n")
  cat("   ---------------------------------\n")
  
  # Calculate quality scores
  quality_metrics <- list()
  
  # 1. Proportion of significant genes
  sig_genes <- unique(de_data[adj.P.Val < 0.05, SYMBOL])
  quality_metrics$prop_significant <- length(sig_genes) / length(unique_genes)
  
  # 2. Proportion with large effect sizes
  large_effect <- unique(de_data[abs(logFC) > 1, SYMBOL])
  quality_metrics$prop_large_effect <- length(large_effect) / length(unique_genes)
  
  # 3. Consistency across datasets (if available)
  if("n_datasets" %in% colnames(de_data)) {
    multi_dataset_genes <- unique(de_data[n_datasets >= 2, SYMBOL])
    quality_metrics$prop_multi_dataset <- length(multi_dataset_genes) / length(unique_genes)
  }
  
  # 4. Signal-to-noise ratio (approximate)
  median_logfc <- median(abs(de_data$logFC), na.rm = TRUE)
  quality_metrics$median_abs_logfc <- median_logfc
  
  cat("   Quality metrics:\n")
  for(metric in names(quality_metrics)) {
    if(is.numeric(quality_metrics[[metric]])) {
      cat(paste0("     ", metric, ": ", round(quality_metrics[[metric]], 4), "\n"))
    }
  }
  
  # ========== PROBLEM IDENTIFICATION ==========
  cat("\n7. PROBLEM IDENTIFICATION:\n")
  cat("   ---------------------------------\n")
  
  problems <- list()
  
  # Check for suspicious P-value distribution
  prop_pval_1 <- sum(de_data$PValue == 1, na.rm = TRUE) / nrow(de_data)
  if(prop_pval_1 > 0.1) {
    problems$high_pval_1 <- paste0("High proportion of P-value = 1: ", 
                                   round(100 * prop_pval_1, 1), "%")
  }
  
  # Check for suspicious adjusted P-value distribution
  prop_adj_pval_1 <- sum(de_data$adj.P.Val == 1, na.rm = TRUE) / nrow(de_data)
  if(prop_adj_pval_1 > 0.1) {
    problems$high_adj_pval_1 <- paste0("High proportion of adj.P.Val = 1: ", 
                                       round(100 * prop_adj_pval_1, 1), "%")
  }
  
  # Check for too few significant genes
  prop_sig <- sum(de_data$adj.P.Val < 0.05, na.rm = TRUE) / nrow(de_data)
  if(prop_sig < 0.01) {
    problems$few_significant <- paste0("Very few significant genes: ", 
                                       round(100 * prop_sig, 2), "%")
  }
  
  # Check for small effect sizes
  median_abs_logfc <- median(abs(de_data$logFC), na.rm = TRUE)
  if(median_abs_logfc < 0.1) {
    problems$small_effects <- paste0("Small median effect size: ", 
                                     round(median_abs_logfc, 3))
  }
  
  if(length(problems) > 0) {
    cat("   ⚠ POTENTIAL PROBLEMS IDENTIFIED:\n")
    for(prob in names(problems)) {
      cat(paste0("     • ", problems[[prob]], "\n"))
    }
  } else {
    cat("   ✅ No major problems identified\n")
  }
  
  # ========== RECOMMENDATIONS ==========
  cat("\n8. RECOMMENDATIONS:\n")
  cat("   ---------------------------------\n")
  
  recommendations <- list()
  
  # Based on significance distribution
  prop_sig_05 <- sum(de_data$adj.P.Val < 0.05, na.rm = TRUE) / nrow(de_data)
  if(prop_sig_05 < 0.05) {
    recommendations$threshold <- "Use lenient thresholds: adj.P.Val < 0.1 or logFC > 0.5"
  } else if(prop_sig_05 > 0.2) {
    recommendations$threshold <- "Use strict thresholds: adj.P.Val < 0.01 and logFC > 1"
  } else {
    recommendations$threshold <- "Use moderate thresholds: adj.P.Val < 0.05 and logFC > 0.5"
  }
  
  # Based on effect sizes
  prop_large_effect <- sum(abs(de_data$logFC) > 1, na.rm = TRUE) / nrow(de_data)
  if(prop_large_effect < 0.1) {
    recommendations$effect_size <- "Focus on genes with |logFC| > 0.5 (moderate effects)"
  }
  
  # Based on dataset consistency
  if("n_datasets" %in% colnames(de_data)) {
    prop_multi_dataset <- sum(de_data$n_datasets >= 2, na.rm = TRUE) / nrow(de_data)
    if(prop_multi_dataset > 0.3) {
      recommendations$consistency <- "Require genes in ≥2 datasets for reproducibility"
    }
  }
  
  cat("   Recommendations for filtering:\n")
  for(rec in names(recommendations)) {
    cat(paste0("     • ", recommendations[[rec]], "\n"))
  }
  
  # ========== SUMMARY REPORT ==========
  cat("\n9. SUMMARY REPORT:\n")
  cat("   ---------------------------------\n")
  
  summary_report <- list(
    "File" = de_file,
    "Total rows" = nrow(de_data),
    "Unique genes" = length(unique_genes),
    "Significant genes (adj.P.Val < 0.05)" = length(sig_genes),
    "Large effect genes (|logFC| > 1)" = length(large_effect),
    "Median |logFC|" = round(median_abs_logfc, 3),
    "Proportion significant" = round(prop_sig, 3)
  )
  
  if("n_datasets" %in% colnames(de_data)) {
    summary_report[["Genes in ≥2 datasets"]] <- length(multi_dataset_genes)
    summary_report[["Max datasets per gene"]] <- max(de_data$n_datasets, na.rm = TRUE)
  }
  
  for(item in names(summary_report)) {
    cat(paste0("     ", item, ": ", summary_report[[item]], "\n"))
  }
  
  # Save summary report
  summary_df <- data.frame(
    Metric = names(summary_report),
    Value = unlist(summary_report),
    stringsAsFactors = FALSE
  )
  
  fwrite(summary_df, file.path(output_dir, "summary_report.tsv"), sep = "\t")
  
  # ========== FINAL OUTPUT ==========
  cat("\n" , strrep("=", 40), "\n", sep = "")
  cat("DIAGNOSTICS COMPLETE\n")
  cat(strrep("=", 40), "\n")
  cat(paste0("Output directory: ", output_dir, "\n"))
  cat("Files created:\n")
  cat("  - dataset_distribution.tsv\n")
  cat("  - gene_summary.tsv\n")
  cat("  - summary_report.tsv\n")
  cat("  - plots/ (directory with diagnostic plots)\n")
  
  return(list(
    de_data = de_data,
    gene_summary = if(exists("gene_summary")) gene_summary else NULL,
    quality_metrics = quality_metrics,
    problems = problems,
    recommendations = recommendations
  ))
}

# ========== RUN DIAGNOSTICS ==========
cat("Starting DE Results Diagnostics...\n\n")

# Define your file
de_file <- "C:/Users/noahm/PycharmProjects/MarbleProject/de_results/combined/combined_DE_results_annotated_cleaned.csv"
output_dir <- "./de_diagnostics"

# Run diagnostics
if(file.exists(de_file)) {
  results <- diagnose_de_results(de_file, output_dir)
  
  # Additional specific checks
  cat("\n" , strrep("=", 40), "\n", sep = "")
  cat("ADDITIONAL SPECIFIC CHECKS\n")
  cat(strrep("=", 40), "\n")
  
  # Load the data for additional checks
  de_data <- fread(de_file)
  
  # Check for ALS-relevant genes
  als_genes <- c("SOD1", "TARDBP", "C9ORF72", "FUS", "VCP", "OPTN", "UBQLN2")
  present_als <- intersect(als_genes, unique(de_data$SYMBOL))
  
  cat(paste0("\nALS-relevant genes in your data: ", length(present_als), "/", length(als_genes), "\n"))
  if(length(present_als) > 0) {
    cat("Present: ", paste(present_als, collapse = ", "), "\n")
    
    # Check their significance
    als_significance <- de_data[SYMBOL %in% present_als, .(
      SYMBOL, 
      min_adj_p = min(adj.P.Val),
      max_logfc = max(abs(logFC)),
      n_datasets = if("n_datasets" %in% colnames(de_data)) first(n_datasets) else NA
    ), by = SYMBOL][order(min_adj_p)]

    cat("\nSignificance of ALS genes:\n")
    print(als_significance)
  } else {
    cat("⚠ WARNING: No key ALS genes found in your DE results!\n")
  }
  
  # Check top genes for biological relevance
  cat("\nTop 20 genes by significance and their potential relevance:\n")
  top_genes <- de_data[, .(
    min_adj_p = min(adj.P.Val),
    max_logfc = max(abs(logFC)),
    n_datasets = ifelse("n_datasets" %in% colnames(.), first(n_datasets), NA)
  ), by = SYMBOL][order(min_adj_p)][1:20]
  
  # Categorize by potential relevance
  neurodegeneration_keywords <- c("APOE", "MAPT", "SNCA", "APP", "PSEN", "PARK", "LRRK", "HTT")
  inflammation_keywords <- c("IL", "TNF", "NFKB", "STAT", "CXCL", "CCL")
  stress_keywords <- c("HSP", "CHOP", "ATF", "XBP", "IRE", "PERK")
  
  top_genes$potential_relevance <- "Other"
  top_genes$potential_relevance[grepl(paste(neurodegeneration_keywords, collapse = "|"), top_genes$SYMBOL)] <- "Neurodegeneration"
  top_genes$potential_relevance[grepl(paste(inflammation_keywords, collapse = "|"), top_genes$SYMBOL)] <- "Inflammation"
  top_genes$potential_relevance[grepl(paste(stress_keywords, collapse = "|"), top_genes$SYMBOL)] <- "Stress Response"
  
  print(top_genes)
  
} else {
  cat(paste0("ERROR: File not found: ", de_file, "\n"))
}