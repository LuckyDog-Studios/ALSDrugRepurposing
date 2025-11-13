# Load required libraries
library(GEOquery)
library(limma)
library(ggplot2)

# Download and load GEO data
gse <- getGEO("GSE112680", GSEMatrix = TRUE)

# Extract the first element (usually the main dataset)
if (class(gse) == "list") {
  gse <- gse[[1]]
}

# Now extract expression data and phenotype data
expr_data <- exprs(gse)
pheno_data <- pData(gse)

# Explore the phenotype data to find the diagnosis column
print(colnames(pheno_data))

# Look for columns containing diagnosis information
diagnosis_cols <- grep("diagnosis", colnames(pheno_data), value = TRUE, ignore.case = TRUE)
print(diagnosis_cols)

# Use the correct diagnosis column - it's "diagnosis:ch1"
print(head(pheno_data$`diagnosis:ch1`))

# Extract diagnosis information from the correct column
groups <- pheno_data$`diagnosis:ch1`
groups <- gsub("diagnosis: ", "", groups)
groups <- gsub(" ", "_", groups)  # Replace spaces with underscores
groups <- gsub(":", "_", groups)  # Replace colons with underscores
groups <- make.names(groups)      # Ensure syntactically valid names

groups <- factor(groups)

# Check the group distribution
print(table(groups))

# Remove any samples with missing diagnosis
keep <- !is.na(groups) & groups != ""
expr_data <- expr_data[, keep]
groups <- groups[keep]
groups <- factor(groups)

print(table(groups))

# Create design matrix
design <- model.matrix(~0 + groups)
colnames(design) <- levels(groups)

print("Design matrix columns:")
print(colnames(design))

# Fit linear model
fit <- lmFit(expr_data, design)

# Create contrast matrix - compare ALS vs CON (control)
# Use the actual cleaned group names from your design matrix
cont.matrix <- makeContrasts(
  ALS_vs_CON = ALS - CON,
  MIM_vs_CON = MIM - CON, 
  ALS_vs_MIM = ALS - MIM,
  levels = design
)

# Fit contrasts
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)

# Get results for ALS vs CON (main comparison)
results_ALS_vs_CON <- topTable(fit2, coef = "ALS_vs_CON", number = Inf, adjust.method = "BH")
head(results_ALS_vs_CON)

# Save results
write.csv(results_ALS_vs_CON, "GSE112680_ALS_vs_CON_limma_results.csv", row.names = TRUE)

# Summary statistics
cat("Summary of differential expression (ALS vs CON):\n")
cat("Total genes:", nrow(results_ALS_vs_CON), "\n")
cat("Significant genes (adj.P.Val < 0.05):", sum(results_ALS_vs_CON$adj.P.Val < 0.05), "\n")
cat("Up-regulated in ALS (adj.P.Val < 0.05 & logFC > 1):", 
    sum(results_ALS_vs_CON$adj.P.Val < 0.05 & results_ALS_vs_CON$logFC > 1), "\n")
cat("Down-regulated in ALS (adj.P.Val < 0.05 & logFC < -1):", 
    sum(results_ALS_vs_CON$adj.P.Val < 0.05 & results_ALS_vs_CON$logFC < -1), "\n")