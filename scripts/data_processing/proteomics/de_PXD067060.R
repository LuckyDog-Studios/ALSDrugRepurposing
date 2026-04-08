# Differential expression analysis for PRIDE dataset PXD067060
# ALS proteomics - MaxQuant LFQ data from cortex (CTX) and spinal cord (SPC)
# Comparison: Late-stage (timepoint 7) vs Early-stage (timepoint 1)
#
# Data location: sark_data/MQ_txt_Sark/proteinGroups.txt

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

required_packages <- c("DEP", "SummarizedExperiment", "limma", "dplyr", "tidyr", "readr", "vsn")
for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    BiocManager::install(pkg)
    library(pkg, character.only = TRUE)
  }
}

# ── Configuration ─────────────────────────────────────────────────────────────
MQ_TXT_DIR     <- "sark_data/MQ_txt_Sark"
OUTPUT_DIR     <- "differential_expression"
EARLY_TP       <- "1"   # timepoint label treated as "early / reference"
LATE_TP        <- "7"   # timepoint label treated as "late / ALS-advanced"
MIN_VALID_FRAC <- 0.70  # minimum fraction of valid (non-zero) values per condition

# ── Load proteinGroups ────────────────────────────────────────────────────────
cat("Loading proteinGroups.txt ...\n")
pg <- read.delim(file.path(MQ_TXT_DIR, "proteinGroups.txt"),
                 stringsAsFactors = FALSE, check.names = FALSE)
cat("Loaded", nrow(pg), "protein groups,", ncol(pg), "columns\n")

# ── Helper: run DEP workflow for one tissue ───────────────────────────────────
run_dep_analysis <- function(pg, tissue, early_tp, late_tp,
                             min_valid_frac = MIN_VALID_FRAC) {

  cat("\n=== Tissue:", tissue, "| Comparison: TP", late_tp, "vs TP", early_tp, "===\n")

  # --- Select LFQ columns for this tissue and the two timepoints of interest --
  lfq_cols <- grep(paste0("^LFQ intensity ", tissue, "_P"),
                   colnames(pg), value = TRUE)

  # Keep only samples whose name ends in _<early_tp> or _<late_tp>
  tp_pattern <- paste0("_(", early_tp, "|", late_tp, ")$")
  lfq_cols   <- grep(tp_pattern, lfq_cols, value = TRUE)

  if (length(lfq_cols) == 0) {
    cat("  No LFQ columns found – skipping.\n")
    return(NULL)
  }
  cat("  Using", length(lfq_cols), "samples\n")

  # --- Build the wide data frame expected by DEP --------------------------------
  # DEP needs: "name" (unique protein ID), "ID" (gene name), and one column per sample
  data_wide <- pg %>%
    select(`Protein IDs`, `Gene names`,
           `Potential contaminant`, `Reverse`, `Only identified by site`,
           all_of(lfq_cols)) %>%
    # Remove contaminants / reverse / only-by-site hits
    filter(is.na(`Potential contaminant`) | `Potential contaminant` != "+",
           is.na(`Reverse`)               | `Reverse`               != "+",
           is.na(`Only identified by site`) | `Only identified by site` != "+") %>%
    # Shorten column names: remove "LFQ intensity " prefix
    rename_with(~ sub("^LFQ intensity ", "", .x), starts_with("LFQ intensity")) %>%
    # Gene name column
    mutate(Gene_names = `Gene names`,
           # Use first gene name for proteins with multiple; fall back to Protein IDs
           name = sub(";.*", "", `Gene names`),
           name = ifelse(name == "" | is.na(name),
                         sub(";.*", "", `Protein IDs`),
                         name)) %>%
    # Make protein names unique (DEP requirement)
    mutate(name = make.unique(name)) %>%
    # Replace 0 with NA (MaxQuant encodes missing as 0)
    mutate(across(all_of(sub("^LFQ intensity ", "", lfq_cols)),
                  ~ ifelse(.x == 0, NA, .x)))

  sample_names <- sub("^LFQ intensity ", "", lfq_cols)

  # --- Experimental design table -----------------------------------------------
  # Sample name format:  <TISSUE>_P<patient>_<timepoint>
  experimental_design <- data.frame(
    label     = sample_names,
    condition = ifelse(grepl(paste0("_", early_tp, "$"), sample_names),
                       paste0("TP", early_tp),
                       paste0("TP", late_tp)),
    replicate = sub(paste0("^", tissue, "_P(\\d+)_.*"), "\\1", sample_names),
    stringsAsFactors = FALSE
  )
  cat("  Experimental design:\n")
  print(table(experimental_design$condition))

  # --- Create SummarizedExperiment via DEP -------------------------------------
  # Build input data.frame with base R to avoid S4Vectors masking dplyr::rename
  se_input <- as.data.frame(data_wide[, c("Protein IDs", "name", sample_names)],
                             stringsAsFactors = FALSE)
  colnames(se_input)[colnames(se_input) == "Protein IDs"] <- "ID"
  col_indices <- which(colnames(se_input) %in% sample_names)

  se <- make_se(se_input, col_indices, experimental_design)

  cat("  SummarizedExperiment:", nrow(se), "proteins,", ncol(se), "samples\n")

  # --- Filter: keep proteins with >= min_valid_frac valid values per condition --
  se_filt <- filter_missval(se, thr = floor((1 - min_valid_frac) *
                                              min(table(experimental_design$condition))))
  cat("  After missing-value filter:", nrow(se_filt), "proteins\n")

  if (nrow(se_filt) < 10) {
    cat("  Too few proteins after filtering – skipping.\n")
    return(NULL)
  }

  # --- Normalise (variance stabilising normalisation) --------------------------
  se_norm <- normalize_vsn(se_filt)

  # --- Impute missing values (MinProb: draws from left tail of distribution) ---
  se_imp <- impute(se_norm, fun = "MinProb", q = 0.01)

  # --- Differential expression test (limma) ------------------------------------
  contrast_name <- paste0("TP", late_tp, "_vs_TP", early_tp)
  se_diff <- test_diff(se_imp, type = "manual",
                       test = contrast_name)

  # --- Add significance calls ---------------------------------------------------
  dep <- add_rejections(se_diff, alpha = 0.05, lfc = log2(1.5))

  # --- Extract results table ---------------------------------------------------
  results <- get_results(dep)

  cat("  Result columns:", paste(colnames(results), collapse = ", "), "\n")

  # --- Dynamically detect column names (DEP naming can vary by version) --------
  logfc_col <- grep("_ratio$",       colnames(results), value = TRUE)[1]
  pval_col  <- grep("_p\\.val$",     colnames(results), value = TRUE)[1]
  padj_col  <- grep("_p\\.adj$",     colnames(results), value = TRUE)[1]
  sig_col   <- grep("_significant$", colnames(results), value = TRUE)[1]

  cat("  Detected columns: logFC =", logfc_col, "| p.val =", pval_col,
      "| p.adj =", padj_col, "| sig =", sig_col, "\n")

  out <- results %>%
    transmute(
      Protein_ID  = ID,
      Gene_Symbol = name,
      logFC       = .data[[logfc_col]],
      P.Value     = .data[[pval_col]],
      adj.P.Val   = .data[[padj_col]],
      significant = .data[[sig_col]],
      tissue      = tissue,
      comparison  = contrast_name
    ) %>%
    arrange(adj.P.Val)

  return(out)
}

# ── Run analysis for each tissue ──────────────────────────────────────────────
if (!dir.exists(OUTPUT_DIR)) {
  dir.create(OUTPUT_DIR, recursive = TRUE)
}

results_ctx <- run_dep_analysis(pg, tissue = "CTX",
                                 early_tp = EARLY_TP, late_tp = LATE_TP)
results_spc <- run_dep_analysis(pg, tissue = "SPC",
                                 early_tp = EARLY_TP, late_tp = LATE_TP)

# ── Save per-tissue CSVs ──────────────────────────────────────────────────────
if (!is.null(results_ctx)) {
  out_path <- file.path(OUTPUT_DIR, "PXD067060_DE_CTX_TP7_vs_TP1.csv")
  write.csv(results_ctx, out_path, row.names = FALSE)
  cat("\nCTX results saved to", out_path, "\n")
  cat("  Proteins tested:", nrow(results_ctx), "\n")
  cat("  Significant (FDR < 0.05, |logFC| > log2(1.5)):",
      sum(results_ctx$significant, na.rm = TRUE), "\n")
}

if (!is.null(results_spc)) {
  out_path <- file.path(OUTPUT_DIR, "PXD067060_DE_SPC_TP7_vs_TP1.csv")
  write.csv(results_spc, out_path, row.names = FALSE)
  cat("\nSPC results saved to", out_path, "\n")
  cat("  Proteins tested:", nrow(results_spc), "\n")
  cat("  Significant (FDR < 0.05, |logFC| > log2(1.5)):",
      sum(results_spc$significant, na.rm = TRUE), "\n")
}

# ── Combined output (both tissues) ───────────────────────────────────────────
all_results <- bind_rows(results_ctx, results_spc)
if (nrow(all_results) > 0) {
  out_path <- file.path(OUTPUT_DIR, "PXD067060_DE_results.csv")
  write.csv(all_results, out_path, row.names = FALSE)
  cat("\nCombined results saved to", out_path, "\n")
}

cat("\nAnalysis complete!\n")
