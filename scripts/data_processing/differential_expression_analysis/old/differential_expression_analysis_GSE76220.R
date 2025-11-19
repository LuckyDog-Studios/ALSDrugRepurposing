library(limma)
library(dplyr)

# Read and inspect the data more carefully
expr_data <- read.table("../../datasets/raw/GSE76220/GSE76220_ALS_LCM_RPKM.txt",
                        header = TRUE, sep = "\t",
                        check.names = FALSE,
                        stringsAsFactors = FALSE,
                        fill = TRUE,
                        quote = "")

print("=== DATA INSPECTION ===")
print(dim(expr_data))
print("First few columns:")
print(head(colnames(expr_data)))

# Let's examine the structure more carefully
print("Data types by column:")
print(sapply(expr_data[1:10, 1:10], class))

# Identify sample columns (those with numeric data)
is_numeric_col <- sapply(expr_data, function(x) {
  num_vals <- suppressWarnings(as.numeric(as.character(x)))
  sum(!is.na(num_vals)) > nrow(expr_data) * 0.5  # More than 50% numeric
})

sample_cols <- which(is_numeric_col)
print("Sample columns identified:")
print(colnames(expr_data)[sample_cols])

# Extract just the expression data
expr_matrix <- as.matrix(expr_data[, sample_cols])
mode(expr_matrix) <- "numeric"

# Remove rows with all zeros or NAs
row_sums <- rowSums(expr_matrix, na.rm = TRUE)
expr_matrix <- expr_matrix[row_sums > 0, ]

print("Filtered expression matrix dimensions:")
print(dim(expr_matrix))

# Check data distribution
print("Expression value summary:")
print(summary(as.vector(expr_matrix)))

# Create groups
sample_names <- colnames(expr_matrix)
groups <- ifelse(grepl("a$", sample_names), "ALS", "Control")
groups <- factor(groups, levels = c("Control", "ALS"))

print("Group distribution:")
print(table(groups))

# Apply robust log2 transformation
# First, check if we need to handle zeros
zero_count <- sum(expr_matrix == 0, na.rm = TRUE)
total_count <- length(expr_matrix)
print(paste("Zero values:", zero_count, "/", total_count,
            "(", round(zero_count/total_count*100, 1), "%)"))

# Use a smaller pseudocount if many zeros
if (zero_count/total_count > 0.1) {
  pseudocount <- 0.1
} else {
  pseudocount <- 1
}

expr_matrix_log2 <- log2(expr_matrix + pseudocount)

print("After log2 transformation:")
print(summary(as.vector(expr_matrix_log2)))

# Create design matrix
design <- model.matrix(~ groups)

# Filter lowly expressed genes
# Keep genes with expression above 1 RPKM in at least 20% of samples
min_samples <- ceiling(ncol(expr_matrix_log2) * 0.2)
keep_genes <- rowSums(expr_matrix_log2 > 1) >= min_samples
expr_matrix_filtered <- expr_matrix_log2[keep_genes, ]

print(paste("Genes kept after filtering:", sum(keep_genes), "/", nrow(expr_matrix_log2)))

# Perform quantile normalization (common for RNA-seq)
library(preprocessCore)
expr_matrix_norm <- normalizeQuantiles(expr_matrix_filtered)

print("After quantile normalization:")
print(summary(as.vector(expr_matrix_norm)))

# Fit the model
fit <- lmFit(expr_matrix_norm, design)
fit <- eBayes(fit)

# Get results
results <- topTable(fit, coef = "groupsALS", number = Inf, adjust.method = "BH")

print("=== RESULTS SUMMARY ===")
print(paste("Total genes tested:", nrow(results)))
print(paste("Genes with adj.P.Val < 0.05:", sum(results$adj.P.Val < 0.05, na.rm = TRUE)))
print(paste("Genes with adj.P.Val < 0.01:", sum(results$adj.P.Val < 0.01, na.rm = TRUE)))

# Add gene symbols
results$Gene_symbol <- rownames(results)

# Extract significant results
sig_results <- results[results$adj.P.Val < 0.05 & !is.na(results$adj.P.Val), ]
sig_results <- sig_results[order(sig_results$P.Value), ]

print("Top 10 most significant genes:")
print(head(sig_results, 10))

# Check effect sizes
print("Fold change distribution in significant genes:")
if(nrow(sig_results) > 0) {
  print(summary(sig_results$logFC))

  cat("\n=== STRONGLY REGULATED GENES ===\n")
  strong_up <- sig_results[sig_results$logFC > 1 & sig_results$adj.P.Val < 0.05, ]
  strong_down <- sig_results[sig_results$logFC < -1 & sig_results$adj.P.Val < 0.05, ]

  cat("Strongly upregulated (logFC > 1):", nrow(strong_up), "\n")
  cat("Strongly downregulated (logFC < -1):", nrow(strong_down), "\n")

  if(nrow(strong_up) > 0) {
    cat("\nTop upregulated genes:\n")
    print(head(strong_up[order(-strong_up$logFC), c("Gene_symbol", "logFC", "P.Value", "adj.P.Val")], 10))
  }

  if(nrow(strong_down) > 0) {
    cat("\nTop downregulated genes:\n")
    print(head(strong_down[order(strong_down$logFC), c("Gene_symbol", "logFC", "P.Value", "adj.P.Val")], 10))
  }
} else {
  cat("No significant genes found at FDR < 0.05\n")
}

# Save results
write.csv(results, "../../datasets/processed/GSE76220/GSE76220_limma_results_CORRECTED.csv", row.names = FALSE)

# Create diagnostic plots
png("./datasets/processed/GSE76220_QC_plots.png", width = 1000, height = 800)
par(mfrow = c(2, 2))

# Distribution of expression values
hist(expr_matrix_norm, breaks = 50, main = "Expression Distribution",
     xlab = "log2(RPKM + pseudocount)")

# Mean-variance relationship
plot(rowMeans(expr_matrix_norm), apply(expr_matrix_norm, 1, sd),
     xlab = "Mean expression", ylab = "Standard deviation",
     main = "Mean-Variance Relationship")

# Volcano plot
volcano_data <- results[!is.na(results$P.Value), ]
plot(volcano_data$logFC, -log10(volcano_data$P.Value),
     xlab = "log2 Fold Change", ylab = "-log10(p-value)",
     main = "Volcano Plot", pch = 20, cex = 0.6,
     col = ifelse(volcano_data$adj.P.Val < 0.05, "red", "gray60"))
abline(h = -log10(0.05), lty = 2, col = "blue")
abline(v = c(-1, 1), lty = 2, col = "blue")

# P-value distribution
hist(volcano_data$P.Value, breaks = 50, main = "P-value Distribution",
     xlab = "P-value")

dev.off()

cat("\nAnalysis complete! Check the QC plots to verify data quality.\n")