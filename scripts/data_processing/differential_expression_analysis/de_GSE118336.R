# Load libraries
library(GEOquery)
library(oligo)     # Essential for CEL files
library(limma)     # Standard for Microarray DE
library(dplyr)
library(data.table)
library(pd.hta.2.0) # The specific platform package your log identified

# -----------------------------
# 1. Robust Data Loader (CEL -> RMA)
# -----------------------------
load_gse118336_rma <- function(gse_id) {
  # 1. Download if missing
  if (!dir.exists(gse_id)) {
    cat("Downloading RAW.tar for", gse_id, "...\n")
    getGEOSuppFiles(gse_id, makeDirectory = TRUE, baseDir = getwd())
  }
  
  # 2. Extract TAR
  raw_tar <- list.files(gse_id, pattern = "\\.tar$", full.names = TRUE)
  untar_dir <- file.path(gse_id, "CEL_files")
  
  if (!dir.exists(untar_dir)) {
    cat("Extracting CEL files...\n")
    dir.create(untar_dir)
    untar(raw_tar, exdir = untar_dir)
  }
  
  # 3. Read and Normalize
  cel_files <- list.files(untar_dir, pattern = "\\.CEL\\.gz$", full.names = TRUE)
  if (length(cel_files) == 0) stop("No CEL files found!")

  cat("Normalizing", length(cel_files), "samples via RMA (this may take a minute)...\n")
  raw_data <- read.celfiles(cel_files)
  eset <- rma(raw_data) # Background correction, normalization, summarization
  
  expr_matrix <- exprs(eset)
  
  # Clean column names to just GSM IDs (e.g., GSM3325490)
  colnames(expr_matrix) <- gsub("^([^_]+)_.*$", "\\1", basename(colnames(expr_matrix)))
  
  return(expr_matrix)
}

# -----------------------------
# 2. Metadata & Mapping (Fixed Matching)
# -----------------------------
get_metadata_gse118336 <- function(gse_id, count_gsm_ids) {
  cat("Downloading metadata for", gse_id, "...\n")
  gse <- getGEO(gse_id, GSEMatrix = TRUE)[[1]]
  pdata <- pData(gse)
  
  # Identify genotype column
  char_col <- grep("characteristics_ch1", colnames(pdata), value = TRUE)[1]
  
  # Assign Groups
  pdata$Genotype <- "Other"
  pdata$Genotype[grepl("wt/wt", pdata[[char_col]], ignore.case = TRUE)] <- "Control"
  pdata$Genotype[grepl("wt/H517D", pdata[[char_col]], ignore.case = TRUE)] <- "FUS_Het"
  pdata$Genotype[grepl("H517D/H517D", pdata[[char_col]], ignore.case = TRUE)] <- "FUS_Mutant"

  mapping_df <- data.frame(
    GSM_ID = as.character(rownames(pdata)),
    Condition = pdata$Genotype,
    stringsAsFactors = FALSE
  )
  
  # Reorder metadata to match matrix columns exactly
  final_metadata <- mapping_df[match(count_gsm_ids, mapping_df$GSM_ID), ]
  rownames(final_metadata) <- count_gsm_ids
  
  cat("Mapping Summary:\n")
  print(table(final_metadata$Condition, useNA = "always"))
  
  return(final_metadata)
}

# -----------------------------
# 3. Optimized DE Analysis (limma + Annotation)
# -----------------------------
run_DE_analysis_limma <- function(expr_matrix, metadata, group1, group2) {
  cat("\nRunning DE:", group1, "vs", group2, "\n")
  
  # Subset samples
  keep_samples <- metadata$Condition %in% c(group1, group2)
  expr_sub <- expr_matrix[, keep_samples]
  meta_sub <- metadata[keep_samples, ]
  
  # Check if we have enough samples
  if (ncol(expr_sub) < 2) {
    warning("Insufficient samples for comparison: ", group1, " vs ", group2)
    return(NULL)
  }

  # Setup groups (Reference group is group2/Control)
  group <- factor(meta_sub$Condition, levels = c(group2, group1))
  design <- model.matrix(~group)
  colnames(design) <- c("Intercept", "Contrast")
  
  # Limma Pipeline
  fit <- lmFit(expr_sub, design)
  fit <- eBayes(fit)
  res <- topTable(fit, coef = "Contrast", number = Inf, adjust.method = "BH")
  
  # Optional: Map Probes to Gene Symbols (Uses the platform package)
  # This makes the data much easier to read for ALS research
  if (require("hta20sttranscriptcluster.db", quietly = TRUE)) {
      anno <- select(hta20sttranscriptcluster.db, keys = rownames(res), 
                     columns = c("SYMBOL", "GENENAME"), keytype = "PROBEID")
      anno <- anno[!duplicated(anno$PROBEID), ]
      res <- cbind(ProbeID = rownames(res), res)
      res <- left_join(res, anno, by = c("ProbeID" = "PROBEID"))
  }
  
  return(res)
}

# -----------------------------
# Execution
# -----------------------------
options(timeout = 3600)

# 1. Load Data (CEL -> RMA)
# This will use the files already extracted in your folder
expr_matrix <- load_gse118336_rma("GSE118336")

# 2. Map Metadata
meta <- get_metadata_gse118336("GSE118336", colnames(expr_matrix))

# 3. Run and Save Comparisons
if(!dir.exists("differential_expression")) dir.create("differential_expression")

comparisons <- list(
  Mutant_vs_Control = c("FUS_Mutant", "Control"),
  Het_vs_Control    = c("FUS_Het", "Control")
)

for (name in names(comparisons)) {
  results <- run_DE_analysis_limma(expr_matrix, meta, comparisons[[name]][1], comparisons[[name]][2])
  
  if (!is.null(results)) {
    write.csv(results, paste0("differential_expression/GSE118336_", name, "_RMA.csv"), row.names = FALSE)
    cat("Saved results for:", name, "\n")
  }
}

cat("\nAnalysis Complete.\n")