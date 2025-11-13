#gse56500_c9 <- read.csv("datasets\\processed\\GSE56500\\results_C9_genes.csv")
#gse56500_sALS <- read.csv("datasets\\processed\\GSE56500\\results_sALS_genes.csv")
#gse112676 <- read.csv("datasets\\processed\\GSE112676\\results_GSE112676.csv")
#gse112680 <- read.csv("datasets\\processed\\GSE112680\\GSE112680_ALS_vs_CON_limma_results.csv")

# Clear environment and start fresh
rm(list = ls())

# Load libraries
library(AnnotationDbi)
library(org.Hs.eg.db)
library(dplyr)

# Read the original CSV file
gse76220 <- read.csv("datasets/processed/GSE76220/GSE76220_limma_results_CORRECTED.csv")

# Inspect the original structure
cat("=== ORIGINAL DATA STRUCTURE ===\n")
dim(gse76220)
head(gse76220)
str(gse76220$Gene_symbol)

# Check for basic data issues
cat("=== DATA QUALITY CHECK ===\n")
cat("Total rows:", nrow(gse76220), "\n")
cat("Missing Gene_symbol:", sum(is.na(gse76220$Gene_symbol)), "\n")
cat("Duplicate Gene_symbol:", sum(duplicated(gse76220$Gene_symbol)), "\n")
cat("Unique Gene_symbol:", length(unique(gse76220$Gene_symbol)), "\n")

# Check statistical columns
stat_cols <- c("logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B")
cat("Rows with all NA statistics:", 
    sum(apply(gse76220[, stat_cols], 1, function(x) all(is.na(x)))), "\n")

# Remove rows with all NA statistics
gse76220_clean <- gse76220[!apply(gse76220[, stat_cols], 1, function(x) all(is.na(x))), ]

# Remove duplicate Gene_symbols (keep first occurrence)
gse76220_clean <- gse76220_clean[!duplicated(gse76220_clean$Gene_symbol), ]

cat("=== AFTER CLEANING ===\n")
cat("Clean rows:", nrow(gse76220_clean), "\n")

# Convert to character (required for mapIds)
gse76220_clean$Gene_symbol <- as.character(gse76220_clean$Gene_symbol)

# Verify conversion
cat("Gene_symbol class:", class(gse76220_clean$Gene_symbol), "\n")

# Map Entrez IDs to HGNC symbols
cat("=== MAPPING WITH org.Hs.eg.db ===\n")
gse76220_clean$hgnc_symbol <- mapIds(org.Hs.eg.db,
                                   keys = gse76220_clean$Gene_symbol,
                                   column = "SYMBOL",
                                   keytype = "ENTREZID",
                                   multiVals = "first")

# Check initial mapping results
cat("Genes with HGNC symbols:", sum(!is.na(gse76220_clean$hgnc_symbol)), "\n")
cat("Genes without HGNC symbols:", sum(is.na(gse76220_clean$hgnc_symbol)), "\n")

# Create clean display names
gse76220_clean$Gene_display <- ifelse(
  !is.na(gse76220_clean$hgnc_symbol),
  gse76220_clean$hgnc_symbol,                    # Use HGNC symbol if available
  paste0("ENTREZ:", gse76220_clean$Gene_symbol)  # Otherwise use Entrez ID
)

# Ensure no empty strings
gse76220_clean$Gene_display[gse76220_clean$Gene_display == ""] <- 
  paste0("ENTREZ:", gse76220_clean$Gene_symbol[gse76220_clean$Gene_display == ""])

cat("=== FINAL VERIFICATION ===\n")
cat("Total genes in clean dataset:", nrow(gse76220_clean), "\n")
cat("Genes with HGNC symbols:", sum(!is.na(gse76220_clean$hgnc_symbol)), "\n")
cat("Genes using Entrez fallback:", sum(is.na(gse76220_clean$hgnc_symbol)), "\n")
cat("Mapping success rate:", 
    round(sum(!is.na(gse76220_clean$hgnc_symbol)) / nrow(gse76220_clean) * 100, 1), "%\n")

# Check samples of both mapped and unmapped
cat("\n=== SAMPLE OF MAPPED GENES ===\n")
mapped_sample <- head(gse76220_clean[!is.na(gse76220_clean$hgnc_symbol), 
                                    c("Gene_symbol", "hgnc_symbol", "Gene_display")])
print(mapped_sample)

cat("\n=== SAMPLE OF UNMAPPED GENES ===\n")
unmapped_sample <- head(gse76220_clean[is.na(gse76220_clean$hgnc_symbol), 
                                      c("Gene_symbol", "hgnc_symbol", "Gene_display")])
print(unmapped_sample)

# Save the final cleaned and annotated dataset
write.csv(gse76220_clean, 
          "datasets/processed/GSE76220/GSE76220_properly_annotated.csv", 
          row.names = FALSE)

cat("Properly annotated dataset saved!\n")