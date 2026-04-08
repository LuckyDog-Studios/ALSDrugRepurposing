#!/usr/bin/env Rscript
# GSE833 Differential Expression & Annotation Validator
# Optimized for HG-U95A (Affymetrix) datasets

# 1. SETUP & DEPENDENCIES
required_packages <- c("dplyr", "stringr", "data.table", "GEOquery", "limma", "Biobase")
new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) {
  if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
  BiocManager::install(new_packages)
}

library(dplyr)
library(stringr)
library(data.table)
library(GEOquery)
library(limma)

# 2. ANNOTATION CHECKER (Your original logic, streamlined)
check_annotation_format <- function(annotations) {
  cat("\n=== Running Annotation Format Analysis ===\n")
  annotations <- annotations[!is.na(annotations) & annotations != ""]
  
  format_tests <- list(
    "AFFYMETRIX_PROBE_ID" = "^\\d+_([asxibf])_at$|^AFFX-",
    "ENTREZ_GENE_ID"      = "^\\d+$",
    "ENSEMBL_GENE_ID"     = "^ENS[A-Z]{0,3}G\\d{11}$",
    "GENE_SYMBOL"         = "^[A-Za-z][A-Za-z0-9-]*[A-Za-z0-9]$"
  )
  
  results <- sapply(format_tests, function(p) {
    sum(grepl(p, annotations)) / length(annotations) * 100
  })
  
  best_match <- names(results)[which.max(results)]
  cat("Primary ID Type Detected:", best_match, "(", round(max(results), 1), "%)\n")
  return(best_match)
}

# 3. MAIN ANALYSIS PIPELINE
run_gse833_analysis <- function() {
  cat("Step 1: Downloading GSE833 from GEO...\n")
  gse <- getGEO("GSE833", destdir = ".", getGPL = TRUE)
  eset <- gse[[1]]
  
  # Step 2: Clean Metadata & Grouping
  groups <- pData(eset)$title
  group_assignment <- ifelse(grepl("Control", groups, ignore.case = TRUE), "Control", "ALS")
  group_factor <- factor(group_assignment, levels = c("Control", "ALS"))
  
  cat("\nGroup counts:\n")
  print(table(Group = group_factor))
  
  # Step 3: Expression Pre-processing
  ex <- exprs(eset)
  qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
  LogC <- (qx[5] > 100) || (qx[6]-qx[1] > 50 && qx[2] > 0)
  if (LogC) { 
    ex[ex <= 0] <- NaN
    ex <- log2(ex) 
    cat("Applied log2 transformation.\n")
  }
  
  # Step 4: Differential Expression with Limma
  design <- model.matrix(~0 + group_factor)
  colnames(design) <- levels(group_factor)
  fit <- lmFit(ex, design)
  cont.matrix <- makeContrasts(ALS - Control, levels=design)
  fit2 <- contrasts.fit(fit, cont.matrix)
  fit2 <- eBayes(fit2)
  
  # Step 5: Extract Results
  stats_results <- topTable(fit2, adjust="fdr", number=Inf) %>%
    as.data.frame() %>%
    tibble::rownames_to_column("Probe_ID")
  
  # Step 6: Smart Mapping of Platform Data
  cat("Mapping Probe IDs to Gene Symbols...\n")
  annot_data <- fData(eset)
  
  # Find the symbol column using a case-insensitive search
  sym_col <- grep("Symbol", colnames(annot_data), ignore.case = TRUE, value = TRUE)[1]
  
  if(is.na(sym_col)) {
    cat("⚠️ Warning: Could not find a 'Symbol' column. Using Probe IDs only.\n")
    final_df <- stats_results %>% mutate(Gene_Symbol = Probe_ID)
  } else {
    annot_clean <- annot_data %>% 
      select(Probe_ID = ID, Gene_Symbol = !!sym(sym_col))
    final_df <- left_join(stats_results, annot_clean, by = "Probe_ID")
  }
  
  return(final_df)
}

main <- function() {
  cat("GSE833 DE Analysis & ID Checker\n")
  cat("===============================\n")
  
  results <- tryCatch({
    run_gse833_analysis()
  }, error = function(e) {
    message("\n❌ Analysis failed: ", e$message)
    return(NULL)
  })
  
  if (!is.null(results)) {
    cat("\nTop 5 Differentially Expressed Genes (by P-Value):\n")
    # Using the new standardized column name: Gene_Symbol
    print(head(results %>% select(Probe_ID, Gene_Symbol, logFC, adj.P.Val), 5))
    
    write.csv(results, "differential_expression/GSE833_DE_Results.csv", row.names = FALSE)
    cat("\nResults saved to: ", getwd(), "/GSE833_DE_Results.csv\n", sep="")
  }
}

main()