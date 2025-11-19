
library(limma)
library(dplyr)
library(ggplot2)

# Your current results
results <- read.csv("../../datasets/processed/GSE76220/GSE76220_limma_results_CORRECTED.csv")

cat("=== CURRENT RESULTS DIAGNOSIS ===\n")
cat("Total genes:", nrow(results), "\n")
cat("Genes with raw p < 0.05:", sum(results$P.Value < 0.05, na.rm = TRUE), "\n")
cat("Genes with FDR < 0.05:", sum(results$adj.P.Val < 0.05, na.rm = TRUE), "\n")
cat("Genes with FDR < 0.1:", sum(results$adj.P.Val < 0.1, na.rm = TRUE), "\n")

# Check the distribution
cat("\nP-value distribution:\n")
print(quantile(results$P.Value, probs = seq(0, 1, 0.1), na.rm = TRUE))

cat("\nFold change distribution:\n")
print(quantile(results$logFC, probs = seq(0, 1, 0.1), na.rm = TRUE))

# Since no genes pass FDR, let's use a less stringent approach
# Option 1: Use raw p-value cutoff
sig_raw <- results[results$P.Value < 0.001 & !is.na(results$P.Value), ]
sig_raw <- sig_raw[order(sig_raw$P.Value), ]

cat("\n=== GENES WITH P < 0.001 ===\n")
cat("Number of genes:", nrow(sig_raw), "\n")

if(nrow(sig_raw) > 0) {
  print(head(sig_raw[, c("Gene_symbol", "logFC", "P.Value", "adj.P.Val")], 20))
}

# Option 2: Use fold change filtering
sig_fc <- results[abs(results$logFC) > 0.5 & results$P.Value < 0.05 & !is.na(results$P.Value), ]
sig_fc <- sig_fc[order(sig_fc$P.Value), ]

cat("\n=== GENES WITH |FC| > 1.4 AND P < 0.05 ===\n")
cat("Number of genes:", nrow(sig_fc), "\n")

if(nrow(sig_fc) > 0) {
  print(head(sig_fc[, c("Gene_symbol", "logFC", "P.Value", "adj.P.Val")], 20))
}

# Create better visualization
create_enhanced_plots <- function(results, output_dir = "./datasets/processed/") {

  # Remove NA p-values
  plot_data <- results[!is.na(results$P.Value), ]

  p1 <- ggplot(plot_data, aes(x = logFC, y = -log10(P.Value))) +
    geom_point(aes(color = adj.P.Val < 0.1), alpha = 0.6) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
    geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", color = "blue") +
    scale_color_manual(values = c("gray", "red")) +
    labs(title = "Volcano Plot - ALS vs Control",
         subtitle = paste("Total genes:", nrow(plot_data),
                         "| P < 0.05:", sum(plot_data$P.Value < 0.05),
                         "| FDR < 0.1:", sum(plot_data$adj.P.Val < 0.1))) +
    theme_minimal()

  p2 <- ggplot(plot_data, aes(x = P.Value)) +
    geom_histogram(bins = 50, fill = "lightblue", color = "black") +
    labs(title = "P-value Distribution",
         x = "P-value", y = "Frequency") +
    theme_minimal()

  # MA plot
  p3 <- ggplot(plot_data, aes(x = AveExpr, y = logFC)) +
    geom_point(aes(color = P.Value < 0.05), alpha = 0.6) +
    scale_color_manual(values = c("gray", "red")) +
    geom_smooth(method = "loess", color = "blue") +
    labs(title = "MA Plot",
         x = "Average Expression", y = "log2 Fold Change") +
    theme_minimal()

  # Combine plots
  library(patchwork)
  combined_plot <- (p1 | p2) / p3

  ggsave(file.path(output_dir, "GSE76220_enhanced_plots.png"),
         combined_plot, width = 12, height = 10, dpi = 300)
}

create_enhanced_plots(results)

# Try alternative analysis approach
cat("\n=== ALTERNATIVE ANALYSIS SUGGESTIONS ===\n")
cat("1. The weak signal could be due to:\n")
cat("   - Small sample size (8 vs 13)\n")
cat("   - High biological variability\n")
cat("   - Subtle disease effects\n")
cat("   - Need for covariate adjustment\n\n")

cat("2. Consider:\n")
cat("   - Using a less stringent FDR cutoff (e.g., 0.1 or 0.2)\n")
cat("   - Focusing on genes with both p < 0.01 and |FC| > 0.5\n")
cat("   - Pathway analysis with relaxed thresholds\n")
cat("   - Checking for batch effects\n")

# Save relaxed threshold results
relaxed_sig <- results[results$adj.P.Val < 0.2 & !is.na(results$adj.P.Val), ]
relaxed_sig <- relaxed_sig[order(relaxed_sig$P.Value), ]

write.csv(relaxed_sig,
          "../../datasets/processed/GSE76220/GSE76220_relaxed_FDR_0.2.csv",
          row.names = FALSE)

cat("\nGenes with FDR < 0.2:", nrow(relaxed_sig), "\n")

# For gene symbol conversion (if needed)
if(all(grepl("^[0-9]+$", results$Gene_symbol))) {
  cat("\n=== GENE ID CONVERSION ===\n")
  cat("Your gene symbols appear to be database IDs. Consider converting to gene names.\n")
  cat("You can use biomaRt or similar tools to convert these to gene symbols.\n")
}

# Final summary
cat("\n=== FINAL RECOMMENDATIONS ===\n")
cat("1. Your analysis is technically correct but shows weak differential expression\n")
cat("2. Consider reporting genes with:\n")
cat("   - FDR < 0.2 (less stringent)\n")
cat("   - OR p < 0.001 + |logFC| > 0.5\n")
cat("3. Focus on biological interpretation of top hits rather than strict FDR cutoff\n")
cat("4. Consider if additional covariates (age, sex, batch) could improve power\n")