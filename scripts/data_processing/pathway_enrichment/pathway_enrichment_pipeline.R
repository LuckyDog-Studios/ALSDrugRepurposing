# =============================================================================
# Multi-Omic Pathway Enrichment Pipeline — ALS Drug Repurposing
# =============================================================================
# Identifies pathways consistently dysregulated at both mRNA and protein level
# in ALS spinal cord and motor cortex, then builds a CMap input signature.
#
# Inputs  : spinal_cord_mRNA_merged.csv, GSE124439_clean.csv (mRNA)
#           41598_2025_11466_MOESM1_ESM.xlsx  (protein, PXD062542)
# Outputs : scripts/data_processing/pathway_enrichment/
# =============================================================================

set.seed(42)

# ── resolve script directory so relative paths work from any working dir ──────
SCRIPT_DIR <- tryCatch(
  dirname(normalizePath(sys.frame(1)$ofile, winslash = "/")),
  error = function(e) getwd()
)
# Fall back to the canonical path when sourced interactively
if (!grepl("pathway_enrichment", SCRIPT_DIR)) {
  SCRIPT_DIR <- "C:/Users/noahm/PycharmProjects/MarbleProject/scripts/data_processing/pathway_enrichment"
}
OUT_DIR    <- SCRIPT_DIR
PROJECT    <- normalizePath(file.path(SCRIPT_DIR, "../../.."), winslash = "/")

# Input paths (actual locations confirmed on disk)
MRNA_SC_PATH  <- file.path(PROJECT, "mRNA_protein_corr_results/spinal_cord_mRNA_merged.csv")
MRNA_CTX_PATH <- file.path(PROJECT, "differential_expression/gene expression/GSE124439_clean.csv")
EXCEL_PATH    <- file.path(PROJECT, "mRNA_protein_corr_results/41598_2025_11466_MOESM1_ESM.xlsx")

# ── print manifest and halt if any file is missing ────────────────────────────
cat(strrep("=", 70), "\n")
cat("Expected input files:\n")
cat(strrep("=", 70), "\n")
inputs <- list(
  "spinal cord mRNA"  = MRNA_SC_PATH,
  "cortex mRNA"       = MRNA_CTX_PATH,
  "proteomics Excel"  = EXCEL_PATH
)
all_ok <- TRUE
for (nm in names(inputs)) {
  exists <- file.exists(inputs[[nm]])
  cat(sprintf("  [%s]  %s\n       %s\n", if (exists) "OK" else "MISSING", nm, inputs[[nm]]))
  if (!exists) all_ok <- FALSE
}
if (!all_ok) stop("One or more input files not found. Aborting.")
cat("\nAll inputs confirmed. Proceeding.\n\n")


# =============================================================================
# STEP 0  Install and load packages
# =============================================================================
cat(strrep("=", 70), "\n")
cat("STEP 0: Installing / loading packages\n")
cat(strrep("=", 70), "\n")

# CRAN packages
cran_pkgs <- c("tidyverse", "readxl", "writexl", "ggplot2", "pheatmap", "ggrepel")
for (pkg in cran_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message(sprintf("  Installing %s ...", pkg))
    install.packages(pkg, repos = "https://cloud.r-project.org", quiet = TRUE)
  }
}

# Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager", repos = "https://cloud.r-project.org", quiet = TRUE)

bioc_pkgs <- c("clusterProfiler", "org.Hs.eg.db", "enrichplot", "ReactomePA", "fgsea", "DOSE")
for (pkg in bioc_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message(sprintf("  Installing %s ...", pkg))
    tryCatch(
      BiocManager::install(pkg, ask = FALSE, update = FALSE, quiet = TRUE),
      error = function(e) stop(sprintf("Failed to install %s: %s", pkg, conditionMessage(e)))
    )
  }
}

suppressPackageStartupMessages({
  library(tidyverse)
  library(readxl)
  library(ggplot2)
  library(pheatmap)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(enrichplot)
  library(ReactomePA)
  library(fgsea)
  library(DOSE)
})
cat("  All packages loaded.\n\n")

# Helper: save a data frame as CSV in OUT_DIR
save_csv <- function(df, fname) {
  write.csv(df, file.path(OUT_DIR, fname), row.names = FALSE)
  invisible(df)
}

# Helper: gracefully extract @result from enrichResult objects
get_result <- function(res) {
  if (is.null(res)) return(NULL)
  tryCatch(as.data.frame(res@result), error = function(e) NULL)
}


# =============================================================================
# STEP 1  Load and prepare all four gene lists
# =============================================================================
cat(strrep("=", 70), "\n")
cat("STEP 1: Loading and preparing gene lists\n")
cat(strrep("=", 70), "\n")

# ── generic cleaner ───────────────────────────────────────────────────────────
clean_gene_list <- function(df, sym_col, lfc_col, pval_col, label) {
  # Use base R renaming to avoid tidyselect/all_of() version conflicts
  colnames(df)[colnames(df) == sym_col]  <- "Gene_Symbol"
  colnames(df)[colnames(df) == lfc_col]  <- "logFC"
  colnames(df)[colnames(df) == pval_col] <- "adjPval"

  df <- df[, c("Gene_Symbol", "logFC", "adjPval"), drop = FALSE]
  df$Gene_Symbol <- as.character(df$Gene_Symbol)
  df$logFC       <- as.numeric(df$logFC)
  df$adjPval     <- as.numeric(df$adjPval)

  df <- df[!is.na(df$Gene_Symbol) & df$Gene_Symbol != "", ]
  df <- df[!is.na(df$logFC) & !is.na(df$adjPval), ]
  df <- df[order(df$adjPval), ]
  df <- df[!duplicated(df$Gene_Symbol), ]          # keep lowest adjPval per symbol
  df <- df[order(-df$logFC), ]

  n_sig <- sum(df$adjPval < 0.05, na.rm = TRUE)
  cat(sprintf("  %s: %d genes | logFC range [%.2f, %.2f] | %d with adjPval < 0.05\n",
              label, nrow(df), min(df$logFC), max(df$logFC), n_sig))
  df
}

# 1a — Spinal cord mRNA
sc_mrna_raw <- read.csv(MRNA_SC_PATH)
sc_mrna <- clean_gene_list(sc_mrna_raw, "Gene_Symbol", "mRNA_logFC_mean",
                           "mRNA_adjPval_max", "SC mRNA")
save_csv(sc_mrna, "SC_mRNA_input.csv")

# 1b — Motor cortex mRNA
ctx_mrna_raw <- read.csv(MRNA_CTX_PATH)
ctx_mrna <- clean_gene_list(ctx_mrna_raw, "Gene_Symbol", "mRNA_logFC",
                            "mRNA_adjPval", "CTX mRNA")
save_csv(ctx_mrna, "CTX_mRNA_input.csv")

# 1c — Spinal cord protein (lumbar — matches mRNA source tissue)
prot_raw <- read_excel(EXCEL_PATH, sheet = "Human limma results")
cat(sprintf("  Excel loaded: %d rows x %d cols\n", nrow(prot_raw), ncol(prot_raw)))

sc_prot <- dplyr::select(prot_raw,
         Gene_Symbol = Gene.Symbol,
         logFC       = M.lumbar_vs_C.lumbar_diff,
         adjPval     = M.lumbar_vs_C.lumbar_p.adj)
sc_prot <- clean_gene_list(sc_prot, "Gene_Symbol", "logFC", "adjPval", "SC protein (lumbar)")
save_csv(sc_prot, "SC_protein_input.csv")

# 1d — Motor cortex protein
ctx_prot <- dplyr::select(prot_raw,
         Gene_Symbol = Gene.Symbol,
         logFC       = M.Mcx_vs_C.Mcx_diff,
         adjPval     = M.Mcx_vs_C.Mcx_p.adj)
ctx_prot <- clean_gene_list(ctx_prot, "Gene_Symbol", "logFC", "adjPval", "CTX protein")
save_csv(ctx_prot, "CTX_protein_input.csv")

cat("\n")


# =============================================================================
# STEP 2  Gene symbol -> Entrez ID conversion
# =============================================================================
cat(strrep("=", 70), "\n")
cat("STEP 2: Gene symbol -> Entrez ID conversion\n")
cat(strrep("=", 70), "\n")

convert_ids <- function(df, label) {
  mapped <- bitr(df$Gene_Symbol, fromType = "SYMBOL", toType = "ENTREZID",
                 OrgDb = org.Hs.eg.db)
  colnames(mapped)[colnames(mapped) == "SYMBOL"]   <- "Gene_Symbol"
  colnames(mapped)[colnames(mapped) == "ENTREZID"] <- "Entrez"
  result <- df %>% inner_join(mapped, by = "Gene_Symbol")
  cat(sprintf("  %s: %d/%d symbols mapped to Entrez IDs\n",
              label, nrow(result), nrow(df)))
  result
}

sc_mrna_ids  <- convert_ids(sc_mrna,  "SC mRNA")
ctx_mrna_ids <- convert_ids(ctx_mrna, "CTX mRNA")
sc_prot_ids  <- convert_ids(sc_prot,  "SC protein")
ctx_prot_ids <- convert_ids(ctx_prot, "CTX protein")

# Build a unified symbol<->Entrez lookup table from all four lists
lookup <- bind_rows(
  dplyr::select(sc_mrna_ids,  Gene_Symbol, Entrez),
  dplyr::select(ctx_mrna_ids, Gene_Symbol, Entrez),
  dplyr::select(sc_prot_ids,  Gene_Symbol, Entrez),
  dplyr::select(ctx_prot_ids, Gene_Symbol, Entrez)
) %>% distinct()
save_csv(lookup, "entrez_symbol_lookup.csv")
cat(sprintf("  Lookup table: %d unique gene-Entrez pairs\n\n", nrow(lookup)))

# Helper: create a named, sorted ranked vector from a mapped data frame
make_ranked <- function(df) {
  ranked <- setNames(df$logFC, df$Entrez)
  sort(ranked, decreasing = TRUE)
}

ranked_sc_mrna  <- make_ranked(sc_mrna_ids)
ranked_ctx_mrna <- make_ranked(ctx_mrna_ids)
ranked_sc_prot  <- make_ranked(sc_prot_ids)
ranked_ctx_prot <- make_ranked(ctx_prot_ids)


# =============================================================================
# STEP 3  GSEA enrichment — KEGG, Reactome, GO-BP for all 4 lists
# =============================================================================
cat(strrep("=", 70), "\n")
cat("STEP 3: GSEA enrichment (KEGG / Reactome / GO-BP)\n")
cat(strrep("=", 70), "\n")

GSEA_PARAMS <- list(
  pvalueCutoff  = 0.05,
  pAdjustMethod = "BH",
  minGSSize     = 10,
  maxGSSize     = 500,
  eps           = 0,
  seed          = TRUE,
  nPermSimple   = 1000
)

# Columns to retain in saved CSV
KEEP_COLS <- c("ID", "Description", "NES", "pvalue", "p.adjust",
               "setSize", "core_enrichment")

run_gsea_safe <- function(label, ranked_vec, type) {
  cat(sprintf("  Running %s for %s ...\n", type, label))
  res <- tryCatch({
    if (type == "KEGG") {
      gseKEGG(geneList      = ranked_vec,
              organism      = "hsa",
              pvalueCutoff  = GSEA_PARAMS$pvalueCutoff,
              pAdjustMethod = GSEA_PARAMS$pAdjustMethod,
              minGSSize     = GSEA_PARAMS$minGSSize,
              maxGSSize     = GSEA_PARAMS$maxGSSize,
              eps           = GSEA_PARAMS$eps,
              seed          = GSEA_PARAMS$seed,
              nPermSimple   = GSEA_PARAMS$nPermSimple)
    } else if (type == "Reactome") {
      gsePathway(geneList      = ranked_vec,
                 organism      = "human",
                 pvalueCutoff  = GSEA_PARAMS$pvalueCutoff,
                 pAdjustMethod = GSEA_PARAMS$pAdjustMethod,
                 minGSSize     = GSEA_PARAMS$minGSSize,
                 maxGSSize     = GSEA_PARAMS$maxGSSize,
                 eps           = GSEA_PARAMS$eps,
                 seed          = GSEA_PARAMS$seed,
                 nPermSimple   = GSEA_PARAMS$nPermSimple)
    } else {  # GO-BP
      gseGO(geneList      = ranked_vec,
            OrgDb         = org.Hs.eg.db,
            ont           = "BP",
            pvalueCutoff  = GSEA_PARAMS$pvalueCutoff,
            pAdjustMethod = GSEA_PARAMS$pAdjustMethod,
            minGSSize     = GSEA_PARAMS$minGSSize,
            maxGSSize     = GSEA_PARAMS$maxGSSize,
            eps           = GSEA_PARAMS$eps,
            seed          = GSEA_PARAMS$seed,
            nPermSimple   = GSEA_PARAMS$nPermSimple)
    }
  }, error = function(e) {
    message(sprintf("  WARNING: %s %s failed: %s", label, type, conditionMessage(e)))
    NULL
  })

  df <- get_result(res)
  if (is.null(df) || nrow(df) == 0) {
    warning(sprintf("  WARNING: %s %s returned 0 significant pathways.", label, type))
    df <- data.frame(ID = character(), Description = character(), NES = numeric(),
                     pvalue = numeric(), p.adjust = numeric(),
                     setSize = integer(), core_enrichment = character())
  } else {
    # Retain only columns that exist
    existing <- intersect(KEEP_COLS, colnames(df))
    df <- df[, existing, drop = FALSE]
    cat(sprintf("    -> %d significant pathways\n", nrow(df)))
  }
  list(result = df, object = res)
}

# Run all 12 combinations
datasets <- list(
  SC_mRNA   = ranked_sc_mrna,
  CTX_mRNA  = ranked_ctx_mrna,
  SC_protein  = ranked_sc_prot,
  CTX_protein = ranked_ctx_prot
)

gsea_results <- list()
for (ds_name in names(datasets)) {
  gsea_results[[ds_name]] <- list()
  for (db in c("KEGG", "Reactome", "GO")) {
    key <- paste0(ds_name, "_", db)
    res_list <- run_gsea_safe(ds_name, datasets[[ds_name]], db)
    gsea_results[[ds_name]][[db]] <- res_list

    # Save CSV
    fname <- paste0(ds_name, "_", db, ".csv")
    save_csv(res_list$result, fname)
  }
}
cat("\n")


# =============================================================================
# STEP 4  Find multi-omic confirmed pathways per tissue
# =============================================================================
cat(strrep("=", 70), "\n")
cat("STEP 4: Finding multi-omic confirmed pathways per tissue\n")
cat(strrep("=", 70), "\n")

# Find pathways significant in both mRNA and protein, with matching NES direction
confirm_pathways <- function(mrna_list, prot_list, tissue_label) {
  confirmed_all <- list()

  for (db in c("KEGG", "Reactome", "GO")) {
    mrna_df <- mrna_list[[db]]$result
    prot_df <- prot_list[[db]]$result

    if (nrow(mrna_df) == 0 || nrow(prot_df) == 0) {
      message(sprintf("  Skipping %s %s — no significant pathways in one or both inputs", tissue_label, db))
      next
    }

    # Inner join on Description (more stable than ID across databases)
    both <- inner_join(
      dplyr::select(mrna_df, ID, Description, NES, p.adjust) %>%
        dplyr::rename(mRNA_NES = NES, mRNA_padj = p.adjust),
      dplyr::select(prot_df, Description, NES, p.adjust) %>%
        dplyr::rename(protein_NES = NES, protein_padj = p.adjust),
      by = "Description"
    ) %>%
      # Keep only pathways where NES directions match
      filter(sign(mRNA_NES) == sign(protein_NES)) %>%
      mutate(
        direction = if_else(mRNA_NES > 0, "UP", "DOWN"),
        Database  = db
      ) %>%
      dplyr::select(Database, ID, Description, direction,
                    mRNA_NES, mRNA_padj, protein_NES, protein_padj)

    if (nrow(both) > 0) {
      cat(sprintf("  %s %s: %d confirmed pathways (%d UP, %d DOWN)\n",
                  tissue_label, db, nrow(both),
                  sum(both$direction == "UP"), sum(both$direction == "DOWN")))
      confirmed_all[[db]] <- both
    } else {
      message(sprintf("  %s %s: 0 pathways with matching direction", tissue_label, db))
    }
  }

  if (length(confirmed_all) == 0) return(NULL)

  combined <- bind_rows(confirmed_all) %>%
    arrange(pmin(mRNA_padj, protein_padj)) %>%
    # Remove cross-database duplicates by Description, keep lowest combined padj
    distinct(Description, .keep_all = TRUE)

  combined
}

sc_confirmed  <- confirm_pathways(gsea_results[["SC_mRNA"]],  gsea_results[["SC_protein"]],  "SC")
ctx_confirmed <- confirm_pathways(gsea_results[["CTX_mRNA"]], gsea_results[["CTX_protein"]], "CTX")

if (!is.null(sc_confirmed)) {
  save_csv(sc_confirmed, "SC_confirmed_pathways.csv")
  cat(sprintf("  SC total: %d confirmed (%d UP, %d DOWN)\n",
              nrow(sc_confirmed), sum(sc_confirmed$direction == "UP"),
              sum(sc_confirmed$direction == "DOWN")))
} else {
  cat("  SC: no confirmed multi-omic pathways\n")
  sc_confirmed <- data.frame()
}

if (!is.null(ctx_confirmed)) {
  save_csv(ctx_confirmed, "CTX_confirmed_pathways.csv")
  cat(sprintf("  CTX total: %d confirmed (%d UP, %d DOWN)\n",
              nrow(ctx_confirmed), sum(ctx_confirmed$direction == "UP"),
              sum(ctx_confirmed$direction == "DOWN")))
} else {
  cat("  CTX: no confirmed multi-omic pathways\n")
  ctx_confirmed <- data.frame()
}
cat("\n")


# =============================================================================
# STEP 5  Cross-tissue intersection
# =============================================================================
cat(strrep("=", 70), "\n")
cat("STEP 5: Cross-tissue intersection\n")
cat(strrep("=", 70), "\n")

if (nrow(sc_confirmed) == 0 || nrow(ctx_confirmed) == 0) {
  warning("  WARNING: one or both tissue confirmed tables are empty — cross-tissue intersection will be empty.")
  cross_tissue <- data.frame(
    Database = character(), ID = character(), Description = character(),
    SC_direction = character(), CTX_direction = character(), consistent = logical(),
    SC_mRNA_NES = numeric(), SC_protein_NES = numeric(),
    CTX_mRNA_NES = numeric(), CTX_protein_NES = numeric(),
    SC_mRNA_padj = numeric(), SC_protein_padj = numeric(),
    CTX_mRNA_padj = numeric(), CTX_protein_padj = numeric()
  )
} else {
  cross_tissue <- inner_join(
    sc_confirmed %>%
      dplyr::rename(SC_direction  = direction,
                    SC_mRNA_NES   = mRNA_NES,   SC_protein_NES  = protein_NES,
                    SC_mRNA_padj  = mRNA_padj,  SC_protein_padj = protein_padj),
    dplyr::select(ctx_confirmed,
                  Description, direction, mRNA_NES, protein_NES, mRNA_padj, protein_padj) %>%
      dplyr::rename(CTX_direction  = direction,
                    CTX_mRNA_NES   = mRNA_NES,  CTX_protein_NES  = protein_NES,
                    CTX_mRNA_padj  = mRNA_padj, CTX_protein_padj = protein_padj),
    by = "Description"
  ) %>%
    mutate(consistent = (SC_direction == CTX_direction))
}

save_csv(cross_tissue, "cross_tissue_confirmed_pathways.csv")
cat(sprintf("  Cross-tissue pathways: %d total\n", nrow(cross_tissue)))
if (nrow(cross_tissue) > 0) {
  cat(sprintf("    Consistent direction (both tissues): %d\n",  sum(cross_tissue$consistent)))
  cat(sprintf("    Tissue-discordant:                   %d\n",  sum(!cross_tissue$consistent)))
}
cat("\n")


# =============================================================================
# STEP 6  Flag UPS / proteostasis-relevant pathways
# =============================================================================
cat(strrep("=", 70), "\n")
cat("STEP 6: Flagging UPS/proteostasis-relevant pathways\n")
cat(strrep("=", 70), "\n")

UPS_KEYWORDS <- paste(
  c("proteasom", "ubiquitin", "autophagy", "protein folding",
    "unfolded protein", "stress granule", "protein degradation",
    "protein homeostasis", "proteostasis", "chaperone",
    "heat shock", "lysosom", "aggrephagy"),
  collapse = "|"
)

if (nrow(cross_tissue) > 0) {
  cross_tissue <- cross_tissue %>%
    mutate(UPS_relevant = grepl(UPS_KEYWORDS, Description, ignore.case = TRUE))
  n_ups <- sum(cross_tissue$UPS_relevant)
  cat(sprintf("  UPS-relevant cross-tissue pathways: %d\n", n_ups))
  if (n_ups > 0) {
    cat("  Pathways:\n")
    cross_tissue %>% filter(UPS_relevant) %>% pull(Description) %>%
      sort() %>% paste("   ", .) %>% cat(sep = "\n")
  }
} else {
  cross_tissue$UPS_relevant <- logical(0)
  n_ups <- 0
  cat("  No cross-tissue pathways to flag.\n")
}

save_csv(cross_tissue, "cross_tissue_confirmed_pathways.csv")   # overwrite with UPS column
cat("\n")


# =============================================================================
# STEP 7  Build CMap input gene lists
# =============================================================================
cat(strrep("=", 70), "\n")
cat("STEP 7: Building CMap input gene lists\n")
cat(strrep("=", 70), "\n")

# Helper: split core_enrichment (Entrez IDs separated by "/"), map to symbols
parse_core_genes <- function(core_str) {
  ids <- unlist(strsplit(as.character(core_str), "/"))
  ids[ids != "" & !is.na(ids)]
}

entrez_to_symbol <- function(entrez_ids) {
  # Use lookup table first; fall back to bitr for any unmapped
  mapped <- lookup %>% filter(Entrez %in% entrez_ids)
  still_missing <- setdiff(entrez_ids, mapped$Entrez)
  if (length(still_missing) > 0) {
    extra <- tryCatch(
      {
        tmp <- bitr(still_missing, fromType = "ENTREZID", toType = "SYMBOL", OrgDb = org.Hs.eg.db)
        colnames(tmp)[colnames(tmp) == "ENTREZID"] <- "Entrez"
        colnames(tmp)[colnames(tmp) == "SYMBOL"]   <- "Gene_Symbol"
        tmp
      },
      error = function(e) data.frame(Entrez = character(), Gene_Symbol = character())
    )
    mapped <- bind_rows(mapped, extra)
  }
  mapped
}

build_cmap_list <- function(pathway_df, direction_filter) {
  if (nrow(pathway_df) == 0) return(tibble(Gene_Symbol = character(), n_pathways_supporting = integer()))

  target_rows <- pathway_df %>%
    filter(consistent == TRUE, SC_direction == direction_filter)

  if (nrow(target_rows) == 0) return(tibble(Gene_Symbol = character(), n_pathways_supporting = integer()))

  # Collect core enrichment genes from BOTH SC and CTX GSEA results for this direction
  # Pull from the corresponding enrichment CSVs using Description matching
  all_db_results <- list()
  for (db in c("KEGG", "Reactome", "GO")) {
    for (tissue in c("SC_mRNA", "SC_protein", "CTX_mRNA", "CTX_protein")) {
      db_res <- gsea_results[[tissue]][[db]]$result
      if (!is.null(db_res) && nrow(db_res) > 0 && "core_enrichment" %in% colnames(db_res)) {
        matching <- db_res %>%
          filter(Description %in% target_rows$Description)
        all_db_results[[paste(tissue, db)]] <- matching
      }
    }
  }

  if (length(all_db_results) == 0) return(tibble(Gene_Symbol = character(), n_pathways_supporting = integer()))

  gene_pathway_counts <- bind_rows(all_db_results) %>%
    rowwise() %>%
    mutate(entrez_list = list(parse_core_genes(core_enrichment))) %>%
    unnest(entrez_list) %>%
    dplyr::rename(Entrez = entrez_list) %>%
    group_by(Entrez) %>%
    summarise(n_pathways_supporting = n_distinct(Description), .groups = "drop")

  # Map Entrez -> Symbol
  sym_map <- entrez_to_symbol(gene_pathway_counts$Entrez)
  result <- gene_pathway_counts %>%
    inner_join(sym_map, by = "Entrez") %>%
    dplyr::select(Gene_Symbol, n_pathways_supporting) %>%
    distinct(Gene_Symbol, .keep_all = TRUE) %>%
    arrange(desc(n_pathways_supporting), Gene_Symbol)

  result
}

cmap_up   <- build_cmap_list(cross_tissue, "UP")    # ALS upregulated -> CMap knock-down
cmap_down <- build_cmap_list(cross_tissue, "DOWN")  # ALS downregulated -> CMap activate

cat(sprintf("  CMap upregulated gene list (to reverse DOWN): %d genes\n",   nrow(cmap_up)))
cat(sprintf("  CMap downregulated gene list (to reverse UP): %d genes\n\n", nrow(cmap_down)))

save_csv(cmap_up   %>% dplyr::select(Gene_Symbol), "CMap_input_upregulated_genes.csv")
save_csv(cmap_down %>% dplyr::select(Gene_Symbol), "CMap_input_downregulated_genes.csv")

# Combined signature with n_pathways_supporting and direction
cmap_sig <- bind_rows(
  cmap_up   %>% mutate(direction = "UP"),
  cmap_down %>% mutate(direction = "DOWN")
) %>%
  dplyr::select(Gene_Symbol, direction, n_pathways_supporting) %>%
  arrange(desc(n_pathways_supporting), Gene_Symbol)

save_csv(cmap_sig, "CMap_input_signature.csv")

if (nrow(cmap_sig) > 0) {
  cat("  Top 10 genes by pathway support:\n")
  print(head(cmap_sig, 10))
}
cat("\n")


# =============================================================================
# STEP 8  Visualizations
# =============================================================================
cat(strrep("=", 70), "\n")
cat("STEP 8: Generating visualizations\n")
cat(strrep("=", 70), "\n")

# ── 8a: Dot plots for each of 12 enrichment results ──────────────────────────
dot_plot_gsea <- function(res_df, title, fname) {
  if (is.null(res_df) || nrow(res_df) < 5) {
    message(sprintf("  Skipping dot plot for '%s' (< 5 significant pathways)", title))
    return(invisible(NULL))
  }

  plot_df <- res_df %>%
    arrange(desc(abs(NES))) %>%
    distinct(Description, .keep_all = TRUE) %>%   # guard against duplicate descriptions
    slice_head(n = 15) %>%
    mutate(Description = str_trunc(Description, 55)) %>%
    distinct(Description, .keep_all = TRUE) %>%   # guard after truncation too
    mutate(Description = factor(Description, levels = rev(unique(Description))))

  p <- ggplot(plot_df, aes(x = NES, y = Description, size = setSize,
                            color = p.adjust)) +
    geom_point() +
    scale_color_gradient(low = "red", high = "blue", name = "adj. p-value") +
    scale_size_continuous(name = "Gene set size", range = c(2, 8)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    labs(title = title, x = "NES", y = NULL) +
    theme_bw(base_size = 10) +
    theme(axis.text.y = element_text(size = 8))

  ggsave(file.path(OUT_DIR, fname), p, width = 9, height = 6, dpi = 150)
  cat(sprintf("  Saved: %s\n", fname))
}

for (ds_name in names(datasets)) {
  for (db in c("KEGG", "Reactome", "GO")) {
    res_df <- gsea_results[[ds_name]][[db]]$result
    title  <- paste(ds_name, db, "GSEA")
    fname  <- paste0(ds_name, "_", db, "_dotplot.png")
    dot_plot_gsea(res_df, title, fname)
  }
}

# ── 8b: Cross-tissue NES heatmap ─────────────────────────────────────────────
if (nrow(cross_tissue) >= 2) {
  nes_mat <- cross_tissue %>%
    dplyr::select(Description, SC_mRNA_NES, SC_protein_NES, CTX_mRNA_NES, CTX_protein_NES) %>%
    column_to_rownames("Description") %>%
    as.matrix()

  # Shorten long pathway names for readability
  rownames(nes_mat) <- str_trunc(rownames(nes_mat), 60)

  # Row annotation: UPS_relevant flag
  row_ann <- data.frame(
    UPS = if_else(cross_tissue$UPS_relevant, "Yes", "No"),
    Direction = cross_tissue$SC_direction,
    Consistent = if_else(cross_tissue$consistent, "Yes", "No"),
    row.names = str_trunc(cross_tissue$Description, 60)
  )
  ann_colors <- list(
    UPS       = c(Yes = "#E53935", No = "#B0BEC5"),
    Direction = c(UP = "#1565C0", DOWN = "#6A1B9A"),
    Consistent = c(Yes = "#2E7D32", No = "#F57F17")
  )

  png(file.path(OUT_DIR, "cross_tissue_heatmap.png"),
      width = 1400, height = max(600, nrow(nes_mat) * 28), res = 150)
  pheatmap(
    nes_mat,
    annotation_row  = row_ann,
    annotation_colors = ann_colors,
    color           = colorRampPalette(c("#1565C0", "white", "#C62828"))(100),
    breaks          = seq(-3, 3, length.out = 101),
    cluster_rows    = TRUE,
    cluster_cols    = FALSE,
    fontsize_row    = 8,
    fontsize_col    = 10,
    main            = "Cross-tissue Confirmed Pathways — NES Heatmap",
    border_color    = NA
  )
  dev.off()
  cat("  Saved: cross_tissue_heatmap.png\n")
} else {
  message("  Skipping heatmap — fewer than 2 cross-tissue pathways")
}

# ── 8c: UPS pathway highlight bar chart ──────────────────────────────────────
ups_df <- cross_tissue %>% filter(UPS_relevant == TRUE)
if (nrow(ups_df) > 0) {
  ups_long <- ups_df %>%
    dplyr::select(Description,
                  SC_mRNA     = SC_mRNA_NES,
                  SC_protein  = SC_protein_NES,
                  CTX_mRNA    = CTX_mRNA_NES,
                  CTX_protein = CTX_protein_NES) %>%
    pivot_longer(-Description, names_to = "Condition", values_to = "NES") %>%
    mutate(Description = str_trunc(Description, 50),
           Condition   = factor(Condition,
                                levels = c("SC_mRNA", "SC_protein",
                                           "CTX_mRNA", "CTX_protein")))

  p_ups <- ggplot(ups_long, aes(x = Condition, y = NES, fill = NES > 0)) +
    geom_col(position = "dodge", show.legend = FALSE) +
    facet_wrap(~ Description, ncol = 2, labeller = label_wrap_gen(30)) +
    scale_fill_manual(values = c("TRUE" = "#C62828", "FALSE" = "#1565C0")) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    labs(title = "UPS/Proteostasis-Relevant Pathways — NES across Conditions",
         x = NULL, y = "NES") +
    theme_bw(base_size = 9) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))

  ggsave(file.path(OUT_DIR, "UPS_pathways_plot.png"), p_ups,
         width = max(8, ceiling(nrow(ups_df) / 2) * 3.5),
         height = ceiling(nrow(ups_df) / 2) * 3.5,
         dpi = 150, limitsize = FALSE)
  cat("  Saved: UPS_pathways_plot.png\n")
} else {
  message("  No UPS-relevant cross-tissue pathways — skipping UPS plot")
}
cat("\n")


# =============================================================================
# STEP 9  Summary report
# =============================================================================
cat(strrep("=", 70), "\n")
cat("STEP 9: Writing summary report\n")
cat(strrep("=", 70), "\n")

# Count significant pathways per enrichment
sig_counts <- sapply(names(datasets), function(ds) {
  sapply(c("KEGG", "Reactome", "GO"), function(db) {
    nrow(gsea_results[[ds]][[db]]$result)
  }, USE.NAMES = TRUE)
}, USE.NAMES = TRUE)

ups_pathway_names <- if (nrow(ups_df) > 0) paste(ups_df$Description, collapse = "\n    ") else "(none)"

top10_genes <- if (nrow(cmap_sig) > 0) {
  head(cmap_sig, 10) %>%
    mutate(line = sprintf("    %s (%s, %d pathways)", Gene_Symbol, direction, n_pathways_supporting)) %>%
    pull(line) %>% paste(collapse = "\n")
} else "(none)"

report <- sprintf(
"Multi-Omic Pathway Enrichment Pipeline — Summary Report
%s

STEP 1 — Gene list sizes after cleaning
  SC mRNA:      %d genes  (logFC range: [%.2f, %.2f])  %d with adjPval < 0.05
  CTX mRNA:     %d genes  (logFC range: [%.2f, %.2f])  %d with adjPval < 0.05
  SC protein:   %d genes  (logFC range: [%.2f, %.2f])  %d with adjPval < 0.05
  CTX protein:  %d genes  (logFC range: [%.2f, %.2f])  %d with adjPval < 0.05

STEP 2 — Entrez ID mapping
  SC mRNA:      %d mapped
  CTX mRNA:     %d mapped
  SC protein:   %d mapped
  CTX protein:  %d mapped

STEP 3 — Significant GSEA pathways (p.adjust < 0.05)
  SC_mRNA    : KEGG=%d  Reactome=%d  GO=%d
  CTX_mRNA   : KEGG=%d  Reactome=%d  GO=%d
  SC_protein : KEGG=%d  Reactome=%d  GO=%d
  CTX_protein: KEGG=%d  Reactome=%d  GO=%d

STEP 4 — Multi-omic confirmed pathways (same direction, both omics)
  Spinal cord:   %d  (UP: %d, DOWN: %d)
  Motor cortex:  %d  (UP: %d, DOWN: %d)

STEP 5 — Cross-tissue confirmed pathways
  Total cross-tissue:         %d
  Consistent direction:       %d
  Tissue-discordant:          %d

STEP 6 — UPS/proteostasis-relevant cross-tissue pathways: %d
  %s

STEP 7 — CMap input signature
  Genes in upregulated list:   %d
  Genes in downregulated list: %d
  Total unique genes:          %d

  Top 10 genes by pathway support:
%s

%s
",
  strrep("=", 70),
  nrow(sc_mrna),  min(sc_mrna$logFC),  max(sc_mrna$logFC),  sum(sc_mrna$adjPval < 0.05),
  nrow(ctx_mrna), min(ctx_mrna$logFC), max(ctx_mrna$logFC), sum(ctx_mrna$adjPval < 0.05),
  nrow(sc_prot),  min(sc_prot$logFC),  max(sc_prot$logFC),  sum(sc_prot$adjPval < 0.05),
  nrow(ctx_prot), min(ctx_prot$logFC), max(ctx_prot$logFC), sum(ctx_prot$adjPval < 0.05),
  nrow(sc_mrna_ids), nrow(ctx_mrna_ids), nrow(sc_prot_ids), nrow(ctx_prot_ids),
  sig_counts["KEGG",   "SC_mRNA"],   sig_counts["Reactome", "SC_mRNA"],   sig_counts["GO",  "SC_mRNA"],
  sig_counts["KEGG",   "CTX_mRNA"],  sig_counts["Reactome", "CTX_mRNA"],  sig_counts["GO",  "CTX_mRNA"],
  sig_counts["KEGG",   "SC_protein"],sig_counts["Reactome", "SC_protein"],sig_counts["GO",  "SC_protein"],
  sig_counts["KEGG",   "CTX_protein"],sig_counts["Reactome","CTX_protein"],sig_counts["GO", "CTX_protein"],
  if (nrow(sc_confirmed)  > 0) nrow(sc_confirmed)  else 0,
  if (nrow(sc_confirmed)  > 0) sum(sc_confirmed$direction  == "UP") else 0,
  if (nrow(sc_confirmed)  > 0) sum(sc_confirmed$direction  == "DOWN") else 0,
  if (nrow(ctx_confirmed) > 0) nrow(ctx_confirmed) else 0,
  if (nrow(ctx_confirmed) > 0) sum(ctx_confirmed$direction == "UP") else 0,
  if (nrow(ctx_confirmed) > 0) sum(ctx_confirmed$direction == "DOWN") else 0,
  nrow(cross_tissue),
  sum(cross_tissue$consistent),
  sum(!cross_tissue$consistent),
  n_ups,
  ups_pathway_names,
  nrow(cmap_up), nrow(cmap_down),
  n_distinct(c(cmap_up$Gene_Symbol, cmap_down$Gene_Symbol)),
  top10_genes,
  strrep("=", 70)
)

cat(report)
writeLines(report, file.path(OUT_DIR, "pathway_enrichment_summary.txt"))
cat("Saved: pathway_enrichment_summary.txt\n")
cat("\nPipeline complete.\n")
