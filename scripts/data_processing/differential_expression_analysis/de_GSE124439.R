# Load libraries
library(GEOquery)
library(edgeR)
library(dplyr)
library(data.table) # Faster for large count files

# -----------------------------
# 1. Robust Count Loader
# -----------------------------
# Using data.table for speed and a more efficient merge approach
read_count_files <- function(count_dir) {
  count_files <- list.files(count_dir, pattern = "_counts\\.txt\\.gz$", full.names = TRUE)
  if (length(count_files) == 0) stop("No count files found!")

  cat("Loading", length(count_files), "files...\n")

  # 1. Read the first file to initialize the matrix and get Gene IDs
  first_df <- fread(count_files[1], data.table = FALSE)
  if (first_df[1,1] == "gene/TE") first_df <- first_df[-1, ]
  
  gene_names <- first_df[[1]]
  sample_name <- gsub("_counts\\.txt\\.gz$", "", basename(count_files[1]))
  
  # 2. Pre-allocate a matrix for speed (Rows = Genes, Cols = Samples)
  expr_matrix <- matrix(NA, nrow = length(gene_names), ncol = length(count_files))
  colnames(expr_matrix) <- gsub("_counts\\.txt\\.gz$", "", basename(count_files))
  rownames(expr_matrix) <- gene_names
  
  # 3. Fill the matrix column by column
  for (i in seq_along(count_files)) {
    # Only read the second column (the counts) to save RAM
    temp_counts <- fread(count_files[i], select = 2, data.table = FALSE)
    
    # Handle that pesky header row if it exists
    if (nrow(temp_counts) > length(gene_names)) {
        temp_counts <- temp_counts[-1, , drop = FALSE]
    }
    
    expr_matrix[, i] <- as.numeric(temp_counts[[1]])
    
    if (i %% 20 == 0) cat("Processed", i, "files...\n")
  }

  return(expr_matrix)
}

# -----------------------------
# 2. Metadata & Mapping
# -----------------------------
get_metadata_and_groups <- function(gse_id, file_names) {
  gse <- getGEO(gse_id, GSEMatrix = TRUE)[[1]]
  pdata <- pData(gse)
  
  # 1. Collapse all characteristic columns into one string per sample to search everything
  all_chars <- apply(pdata[, grep("characteristics", colnames(pdata))], 1, paste, collapse = " ")

  # 2. Create the mapping table
  mapping_df <- data.frame(
    GSM_ID = rownames(pdata),
    Title_ID = as.character(pdata$title), 
    Condition = ifelse(grepl("ALS", all_chars, ignore.case = TRUE), "ALS", 
                ifelse(grepl("Control", all_chars, ignore.case = TRUE), "Control", "Other")),
    stringsAsFactors = FALSE
  )

  # 3. Fuzzy match filenames (same as before)
  matched_indices <- sapply(file_names, function(x) {
    idx <- which(sapply(mapping_df$Title_ID, function(id) grepl(id, x)))
    if(length(idx) > 0) return(idx[1]) else return(NA)
  })

  final_metadata <- mapping_df[matched_indices, ]
  rownames(final_metadata) <- file_names
  
  cat("New Mapping Summary:\n")
  print(table(final_metadata$Condition, useNA = "always"))
  
  return(final_metadata)
}

# -----------------------------
# 3. Optimized DE Analysis
# -----------------------------
run_DE_analysis <- function(counts, metadata) {
  # 1. Subset to only relevant samples immediately
  keep_samples <- !is.na(metadata$Condition) & metadata$Condition %in% c("ALS", "Control")
  
  if (sum(keep_samples) < 2) {
    stop("Error: Insufficient samples (less than 2) matched 'ALS' or 'Control'.")
  }

  counts_sub <- counts[, keep_samples]
  meta_sub <- metadata[keep_samples, ]
  
  # 2. Ensure Factor Levels
  group <- factor(meta_sub$Condition, levels = c("Control", "ALS"))
  
  # 3. edgeR Workflow
  dge <- DGEList(counts = counts_sub, group = group)
  
  # Use filterByExpr which is smarter about group sizes
  keep_genes <- filterByExpr(dge, design = model.matrix(~group))
  dge <- dge[keep_genes, , keep.lib.sizes = FALSE]
  
  dge <- calcNormFactors(dge)
  design <- model.matrix(~group)
  dge <- estimateDisp(dge, design)
  
  fit <- glmQLFit(dge, design)
  qlf <- glmQLFTest(fit, coef = 2) 
  
  # 4. Explicitly sort by p-value
  res <- topTags(qlf, n = Inf, adjust.method = "BH")$table
  return(res)
}

# -----------------------------
# Execution
# -----------------------------
# 1. Load Data
counts <- read_count_files("GSE124439/counts")

# 2. Map Metadata
# Ensure your working directory contains the files or update path
meta <- get_metadata_and_groups("GSE124439", colnames(counts))

# 3. Run DE
results <- run_DE_analysis(counts, meta)

# 4. Save
write.csv(results, "differential_expression/GSE124439_ALS_vs_Control_Results_2.csv")
cat("Analysis Complete. Results saved.\n")