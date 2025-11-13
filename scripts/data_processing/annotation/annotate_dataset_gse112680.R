# Clear environment and start fresh
rm(list = ls())

# Load libraries
library(AnnotationDbi)
library(org.Hs.eg.db)
library(dplyr)

# Read the GSE112676 CSV file
gse112680 <- read.csv("datasets/processed/GSE112680/GSE112680_ALS_vs_CON_limma_results.csv")

# Inspect the original structure
cat("=== GSE112680 DATA STRUCTURE ===\n")
dim(gse112680)
head(gse112680)
colnames(gse112680)

# Check what the first column contains (it's the probe IDs)
cat("First column name:", colnames(gse112680)[1], "\n")
cat("Sample probe IDs:\n")
print(head(gse112680[,1]))

# Extract probe IDs and rename the column
probe_ids <- gse112680[,1]
colnames(gse112680)[1] <- "ProbeID"

cat("Number of probe IDs:", length(probe_ids), "\n")

# Install and load the appropriate Illumina annotation package
# Based on the "ILMN_" prefix, this is likely from Illumina HumanHT-12 arrays

cat("=== USING ILLUMINA ANNOTATION PACKAGE ===\n")

# Try different Illumina annotation packages
illumina_packages <- c("illuminaHumanv4.db", "illuminaHumanv3.db", "illuminaHumanv2.db")

selected_package <- NULL
for(pkg in illumina_packages) {
  cat("Trying", pkg, "...\n")
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    cat("Installing", pkg, "...\n")
    tryCatch({
      BiocManager::install(pkg)
      library(pkg, character.only = TRUE)
      selected_package <- pkg
      cat("Successfully loaded", pkg, "\n")
      break
    }, error = function(e) {
      cat("Failed to install/load", pkg, "\n")
    })
  } else {
    selected_package <- pkg
    cat("Successfully loaded", pkg, "\n")
    break
  }
}

if(is.null(selected_package)) {
  stop("Could not load any Illumina annotation package. Please install manually.")
}

cat("Using annotation package:", selected_package, "\n")

# Map probes to Entrez IDs
cat("Mapping probes to Entrez IDs...\n")
probe_entrez <- mapIds(get(selected_package),
                      keys = probe_ids,
                      column = "ENTREZID",
                      keytype = "PROBEID",
                      multiVals = "first")

# Map probes to gene symbols directly
cat("Mapping probes to gene symbols...\n")
probe_symbols <- mapIds(get(selected_package),
                       keys = probe_ids,
                       column = "SYMBOL",
                       keytype = "PROBEID",
                       multiVals = "first")

# Create mapping dataframe
probe_mapping <- data.frame(
  ProbeID = names(probe_entrez),
  entrezgene_id = as.character(probe_entrez),
  hgnc_symbol = probe_symbols,
  stringsAsFactors = FALSE
)

# Remove rows where ProbeID is NA
probe_mapping <- probe_mapping[!is.na(probe_mapping$ProbeID), ]

cat("Probe mapping results:", nrow(probe_mapping), "probes mapped\n")

# Merge with the original data
gse112680_annotated <- gse112680 %>%
  left_join(probe_mapping, by = "ProbeID")

# Create display names
gse112680_annotated$Gene_display <- ifelse(
  !is.na(gse112680_annotated$hgnc_symbol),
  gse112680_annotated$hgnc_symbol,
  ifelse(!is.na(gse112680_annotated$entrezgene_id),
         paste0("ENTREZ:", gse112680_annotated$entrezgene_id),
         gse112680_annotated$ProbeID)  # Fallback to probe ID
)

# Final verification
cat("=== FINAL GSE112676 ANNOTATION RESULTS ===\n")
cat("Total probes:", nrow(gse112680_annotated), "\n")
cat("Probes with HGNC symbols:", sum(!is.na(gse112680_annotated$hgnc_symbol)), "\n")
cat("Probes with Entrez IDs:", sum(!is.na(gse112680_annotated$entrezgene_id)), "\n")
cat("Probes using probe ID fallback:", sum(is.na(gse112680_annotated$hgnc_symbol) & 
                                           is.na(gse112680_annotated$entrezgene_id)), "\n")

# Check samples
cat("\n=== SAMPLE OF ANNOTATED PROBES ===\n")
sample_data <- head(gse112680_annotated[, c("ProbeID", "entrezgene_id", "hgnc_symbol", "Gene_display")])
print(sample_data)

# Check if there are any significant probes and their annotation
significant_probes <- gse112680_annotated[gse112680_annotated$significant == TRUE, ]
cat("\n=== SIGNIFICANT PROBES ANNOTATION ===\n")
cat("Significant probes:", nrow(significant_probes), "\n")
cat("Significant probes with HGNC symbols:", sum(!is.na(significant_probes$hgnc_symbol)), "\n")

if(nrow(significant_probes) > 0) {
  cat("Sample of significant annotated probes:\n")
  print(head(significant_probes[, c("ProbeID", "logFC", "adj.P.Val", "hgnc_symbol", "Gene_display")]))
}

# Save the annotated dataset
write.csv(gse112680_annotated, 
          "datasets/processed/GSE112680/GSE112680_annotated.csv", 
          row.names = FALSE)

cat("GSE112680 annotated dataset saved!\n")