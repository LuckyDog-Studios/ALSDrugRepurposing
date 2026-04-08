# Differential expression analysis for PRIDE dataset PXD065424
# G3BP1 co-IP LFQ proteomics in C9ORF72-ALS vs Control iPSC-derived motor neurons
# Conditions: Basal and Arsenite-induced stress granule assembly
# Comparison: C9ALS vs Control (within each condition)
#
# Dataset downloads automatically from PRIDE on first run.
# Output saved to differential_expression/PXD065424_DE_*.csv

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

required_packages <- c("rpx", "DEP", "SummarizedExperiment", "limma",
                       "dplyr", "tidyr", "readr", "vsn")
for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    BiocManager::install(pkg)
    library(pkg, character.only = TRUE)
  }
}

# ── Configuration ─────────────────────────────────────────────────────────────
PRIDE_ID    <- "PXD065424"
DATA_DIR    <- file.path("pride_data", PRIDE_ID)
OUTPUT_DIR  <- "differential_expression"
MIN_VALID_FRAC <- 0.70   # min fraction of valid values per condition for filtering

# ── Download from PRIDE ───────────────────────────────────────────────────────
if (!dir.exists(DATA_DIR)) dir.create(DATA_DIR, recursive = TRUE)

pg_file <- file.path(DATA_DIR, "proteinGroups.txt")

if (!file.exists(pg_file)) {
  px         <- PXDataset(PRIDE_ID)
  base_url   <- pxurl(px)   # e.g. ftp://ftp.pride.ebi.ac.uk/.../PXD065424
  all_files  <- pxfiles(px)

  pride_download <- function(filename, dest_path) {
    url <- paste0(base_url, "/", filename)
    cat("Downloading", filename, "\n  from:", url, "\n")
    download.file(url, dest_path, mode = "wb")
  }

  # Try a direct proteinGroups.txt first; fall back to the zip
  pg_remote <- grep("proteinGroups", all_files, value = TRUE, ignore.case = TRUE)

  if (length(pg_remote) > 0) {
    pride_download(pg_remote[1], pg_file)

  } else {
    zip_remote <- grep("\\.zip$", all_files, value = TRUE, ignore.case = TRUE)
    if (length(zip_remote) == 0)
      stop("Cannot find proteinGroups.txt or a zip in PRIDE. Files: ",
           paste(all_files, collapse = ", "))

    zip_path <- file.path(DATA_DIR, zip_remote[1])
    pride_download(zip_remote[1], zip_path)

    cat("Extracting", zip_path, "...\n")
    unzip(zip_path, exdir = DATA_DIR)

    found <- list.files(DATA_DIR, pattern = "proteinGroups\\.txt$",
                        recursive = TRUE, full.names = TRUE, ignore.case = TRUE)
    if (length(found) == 0)
      stop("proteinGroups.txt not found inside the zip. Extracted to: ", DATA_DIR)

    cat("Found:", found[1], "\n")
    if (normalizePath(found[1]) != normalizePath(pg_file))
      file.copy(found[1], pg_file)
  }
}

# ── Load proteinGroups ────────────────────────────────────────────────────────
cat("Loading proteinGroups.txt ...\n")
pg <- read.delim(pg_file, stringsAsFactors = FALSE, check.names = FALSE)
cat("Loaded", nrow(pg), "protein groups,", ncol(pg), "columns\n")

# Print all LFQ column names so the user can verify group assignments
lfq_all <- grep("^LFQ intensity ", colnames(pg), value = TRUE)
cat("\nAll LFQ intensity columns found:\n")
print(sub("^LFQ intensity ", "", lfq_all))

# ── Assign experimental groups from sample names ──────────────────────────────
# Actual naming: {genotype}_{condition}_{replicate}
#   Genotype:  WT (control), CS29 / CS30 (C9ORF72 patient lines → C9ALS)
#   Condition: CTRL (basal), ARS (arsenite stress),
#              RNase_CTRL (RNase-treated basal), RNase_ARS (RNase-treated stress)
#   FLAG samples = IP negative controls → excluded (condition left NA)
infer_groups <- function(sample_names) {
  genotype <- dplyr::case_when(
    grepl("^WT_",    sample_names) ~ "Control",
    grepl("^CS[0-9]+_", sample_names) ~ "C9ALS",
    TRUE ~ NA_character_
  )
  condition <- dplyr::case_when(
    grepl("_RNase_ARS_",  sample_names) ~ "RNase_ARS",
    grepl("_RNase_CTRL_", sample_names) ~ "RNase_CTRL",
    grepl("_ARS_",        sample_names) ~ "ARS",
    grepl("_CTRL_",       sample_names) ~ "CTRL",
    TRUE ~ NA_character_   # FLAG and anything else excluded
  )
  data.frame(sample = sample_names, genotype = genotype, condition = condition,
             stringsAsFactors = FALSE)
}

sample_names_all <- sub("^LFQ intensity ", "", lfq_all)
group_df <- infer_groups(sample_names_all)

cat("\nInferred group assignments:\n")
print(group_df)

cat("\nGroup assignments (FLAG samples excluded as IP negative controls):\n")
print(group_df)

excluded <- is.na(group_df$condition)
if (any(excluded))
  cat("\nExcluded samples:", paste(group_df$sample[excluded], collapse = ", "), "\n")

# ── Helper: DEP workflow for one condition ─────────────────────────────────────
run_dep_analysis <- function(pg, group_df, condition_label, min_valid_frac) {

  cat("\n=== Condition:", condition_label,
      "| Comparison: C9ALS vs Control ===\n")

  # Select samples for this condition with known genotype (exclude NA conditions)
  sel <- !is.na(group_df$condition) & group_df$condition == condition_label & !is.na(group_df$genotype)
  sub_df <- group_df[sel, ]

  if (nrow(sub_df) == 0) {
    cat("  No samples found for this condition – skipping.\n")
    return(NULL)
  }

  n_c9  <- sum(sub_df$genotype == "C9ALS")
  n_ctl <- sum(sub_df$genotype == "Control")
  cat("  C9ALS:", n_c9, "| Control:", n_ctl, "\n")

  if (n_c9 < 2 || n_ctl < 2) {
    cat("  Need at least 2 replicates per group – skipping.\n")
    return(NULL)
  }

  # Build the data.frame for make_se
  sample_cols <- sub_df$sample
  data_wide <- pg %>%
    select(`Protein IDs`, `Gene names`,
           `Potential contaminant`, `Reverse`, `Only identified by site`,
           all_of(paste0("LFQ intensity ", sample_cols))) %>%
    filter(is.na(`Potential contaminant`) | `Potential contaminant` != "+",
           is.na(`Reverse`)               | `Reverse`               != "+",
           is.na(`Only identified by site`) | `Only identified by site` != "+") %>%
    rename_with(~ sub("^LFQ intensity ", "", .x), starts_with("LFQ intensity")) %>%
    mutate(
      name = sub(";.*", "", `Gene names`),
      name = ifelse(name == "" | is.na(name), sub(";.*", "", `Protein IDs`), name),
      name = make.unique(name)
    ) %>%
    mutate(across(all_of(sample_cols), ~ ifelse(.x == 0, NA, .x)))

  # Build input for make_se (avoid S4Vectors masking dplyr::rename)
  se_input <- as.data.frame(data_wide[, c("Protein IDs", "name", sample_cols)],
                             stringsAsFactors = FALSE)
  colnames(se_input)[colnames(se_input) == "Protein IDs"] <- "ID"
  col_indices <- which(colnames(se_input) %in% sample_cols)

  experimental_design <- data.frame(
    label     = sample_cols,
    condition = sub_df$genotype,
    replicate = seq_len(nrow(sub_df)),
    stringsAsFactors = FALSE
  )

  se <- make_se(se_input, col_indices, experimental_design)
  cat("  SummarizedExperiment:", nrow(se), "proteins,", ncol(se), "samples\n")

  # Filter missing values
  min_valid <- min(table(experimental_design$condition))
  thr <- floor((1 - min_valid_frac) * min_valid)
  se_filt <- filter_missval(se, thr = thr)
  cat("  After missing-value filter:", nrow(se_filt), "proteins\n")

  if (nrow(se_filt) < 10) {
    cat("  Too few proteins after filtering – skipping.\n")
    return(NULL)
  }

  # Normalise and impute
  se_norm <- normalize_vsn(se_filt)
  se_imp  <- impute(se_norm, fun = "MinProb", q = 0.01)

  # DE test: C9ALS vs Control
  contrast_name <- "C9ALS_vs_Control"
  se_diff <- test_diff(se_imp, type = "manual", test = contrast_name)
  dep     <- add_rejections(se_diff, alpha = 0.05, lfc = log2(1.5))
  results <- get_results(dep)

  cat("  Result columns:", paste(colnames(results), collapse = ", "), "\n")

  # Dynamically detect column names
  logfc_col <- grep("_ratio$",       colnames(results), value = TRUE)[1]
  pval_col  <- grep("_p\\.val$",     colnames(results), value = TRUE)[1]
  padj_col  <- grep("_p\\.adj$",     colnames(results), value = TRUE)[1]
  sig_col   <- grep("_significant$", colnames(results), value = TRUE)[1]

  out <- results %>%
    transmute(
      Protein_ID  = ID,
      Gene_Symbol = name,
      logFC       = as.double(.data[[logfc_col]]),
      P.Value     = as.double(.data[[pval_col]]),
      adj.P.Val   = as.double(.data[[padj_col]]),
      significant = .data[[sig_col]],
      condition   = condition_label,
      comparison  = contrast_name
    ) %>%
    mutate(across(where(is.double), ~ { .x[is.infinite(.x)] <- NA; .x })) %>%
    arrange(adj.P.Val)

  return(out)
}

# ── Run for each condition ─────────────────────────────────────────────────────
if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)

conditions <- unique(na.omit(group_df$condition))
all_results <- list()

for (cond in conditions) {
  res <- run_dep_analysis(pg, group_df, cond, MIN_VALID_FRAC)
  if (!is.null(res)) {
    all_results[[cond]] <- res
    out_path <- file.path(OUTPUT_DIR,
                          paste0("PXD065424_DE_", cond, "_C9ALS_vs_Control.csv"))
    write.csv(res, out_path, row.names = FALSE, na = "")
    cat("\n", cond, "results saved to", out_path, "\n")
    cat("  Proteins tested:", nrow(res), "\n")
    cat("  Significant (FDR < 0.05, |logFC| > log2(1.5)):",
        sum(res$significant, na.rm = TRUE), "\n")
  }
}

# Combined output
if (length(all_results) > 0) {
  combined <- bind_rows(all_results)
  out_path <- file.path(OUTPUT_DIR, "PXD065424_DE_results.csv")
  write.csv(combined, out_path, row.names = FALSE, na = "")
  cat("\nCombined results saved to", out_path, "\n")
}

cat("\nAnalysis complete!\n")
