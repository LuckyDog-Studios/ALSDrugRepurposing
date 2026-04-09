"""
SetA_pipeline.py
================
Multi-omic Set A construction for ALS drug repurposing STRING network.

Set A = genes with concordant differential expression in BOTH mRNA AND protein
        data from the SAME tissue (spinal cord or motor cortex).

Inputs:
  mRNA_protein_corr_results/spinal_cord_mRNA_merged.csv           (SC mRNA)
  differential_expression/gene expression/GSE124439_clean.csv     (CTX mRNA)
  differential_expression/protein expression/PXD062542_MOESM1_ESM.xlsx (both protein)

Outputs saved to: scripts/data_processing/multiomics_filtering/
"""

# ── pip installs (non-standard) ──────────────────────────────────────────────
# pip install pandas numpy openpyxl

import sys
import os
import re
import textwrap

import pandas as pd
import numpy as np

# ── Paths ────────────────────────────────────────────────────────────────────
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.abspath(os.path.join(SCRIPT_DIR, "..", "..", ".."))
OUT_DIR = SCRIPT_DIR  # all outputs go here

INPUT_PATHS = {
    "SC_mRNA":   os.path.join(PROJECT_ROOT, "mRNA_protein_corr_results",
                              "spinal_cord_mRNA_merged.csv"),
    "CTX_mRNA":  os.path.join(PROJECT_ROOT, "differential_expression",
                              "gene expression", "GSE124439_clean.csv"),
    "PROTEIN":   os.path.join(PROJECT_ROOT, "differential_expression",
                              "protein expression",
                              "PXD062542_MOESM1_ESM.xlsx"),
}

# ── Thresholds ────────────────────────────────────────────────────────────────
MRNA_FC_THRESH   = 0.5
MRNA_PADJ_THRESH = 0.05
PROT_FC_THRESH   = 0.3
PROT_PADJ_THRESH = 0.05

# ── Reference gene lists (Step 7) ─────────────────────────────────────────────
KNOWN_ALS = {
    "SOD1","TARDBP","FUS","C9orf72","UBQLN2","OPTN","TBK1","VCP","SQSTM1",
    "HNRNPA1","HNRNPA2B1","MATR3","TIA1","ATXN2","NEK1","CCNF","KIF5A",
    "SETX","ANG","VAPB","SIGMAR1",
}
UPS_RELATED = {
    "PSMD1","PSMD2","PSMD3","PSMD4","PSMD7","PSMD11","PSMD14",
    "PSMB1","PSMB2","PSMB3","PSMB4","PSMB5","PSMB6","PSMB7",
    "PSMA1","PSMA2","PSMA3","PSMA4","PSMA5","PSMA6","PSMA7",
    "UBA1","UBA2","UBA3","UBB","UBC","UCHL1","UCHL5",
    "SQSTM1","UBQLN1","UBQLN2","HDAC6","VCP","p62",
}
DISCORDANT_LOW = {
    "DBT","GCN1L1","ILK","IMPA1","MBP","MPZ","MX1","MYL6","NCAM1",
    "PEA15","PPP6C","PRPH","SFPQ","ANPEP","CD34","CFB","COL5A3",
    "CRABP2","FCHSD1","GGT5","GPX3","GSTM1","HLA-C","ICAM1",
    "IGF2BP2","ITGB3","ITIH4","MDK","PODN","SAA2","SVIL",
    "THNSL2","TNXB","TYMP","VTN",
}
STRESS_GRANULE = {
    "SPON1","ITIH4","CFL2","DBN1","VGF","CKAP5","NTN1",
    "STIP1","GNAS","IGHG1","RBM6",
}

# ── Helpers ───────────────────────────────────────────────────────────────────
HGNC_RE = re.compile(r'^[A-Za-z][A-Za-z0-9\-]*$')

def is_valid_symbol(s):
    """Return True if s looks like a valid HGNC gene symbol."""
    if not isinstance(s, str):
        return False
    s = s.strip()
    if not s:
        return False
    return bool(HGNC_RE.match(s))

def clean_symbols(df, col="Gene_Symbol"):
    """Strip whitespace, drop NA / empty / invalid HGNC symbols."""
    df = df.copy()
    df[col] = df[col].astype(str).str.strip()
    mask = df[col].apply(is_valid_symbol)
    n_dropped = (~mask).sum()
    df = df[mask].copy()
    return df, n_dropped

def direction(logfc_series):
    return np.where(logfc_series > 0, "UP", "DOWN")

def out(fname):
    return os.path.join(OUT_DIR, fname)

def sep(title="", width=60):
    print("\n" + "=" * width)
    if title:
        print(f"  {title}")
        print("=" * width)

# ── Preflight: confirm inputs exist ───────────────────────────────────────────
print("\nExpected input file paths:")
all_ok = True
for label, path in INPUT_PATHS.items():
    exists = os.path.isfile(path)
    status = "OK" if exists else "MISSING"
    print(f"  [{status}] {label}: {path}")
    if not exists:
        all_ok = False

if not all_ok:
    print("\nERROR: One or more input files are missing. "
          "Please check paths and re-run.")
    sys.exit(1)

print("\nAll input files found. Starting pipeline...\n")

# ─────────────────────────────────────────────────────────────────────────────
# STEP 1 — Load and inspect all inputs
# ─────────────────────────────────────────────────────────────────────────────
sep("STEP 1 — Load & inspect inputs")

report_lines = []

def inspect(label, df, sym_col="Gene_Symbol"):
    lines = [
        f"\n{'-'*60}",
        f"Source : {label}",
        f"Shape  : {df.shape}",
        f"Columns: {list(df.columns)}",
        "First 3 rows:",
        df.head(3).to_string(index=False),
        f"N rows missing {sym_col}: "
        f"{df[sym_col].isna().sum() if sym_col in df.columns else 'N/A (col absent)'}",
    ]
    for l in lines:
        print(l)
    report_lines.extend(lines)

# Load SC mRNA
sc_mrna_raw = pd.read_csv(INPUT_PATHS["SC_mRNA"])
inspect("spinal_cord_mRNA_merged.csv", sc_mrna_raw)

# Load CTX mRNA
ctx_mrna_raw = pd.read_csv(INPUT_PATHS["CTX_mRNA"])
inspect("GSE124439_clean.csv", ctx_mrna_raw)

# Load protein Excel (Human limma results sheet)
try:
    prot_raw = pd.read_excel(
        INPUT_PATHS["PROTEIN"],
        sheet_name="Human limma results",
        engine="openpyxl",
    )
except Exception as e:
    print(f"\nERROR loading protein Excel: {e}")
    sys.exit(1)

inspect("PXD062542 — spinal cord protein", prot_raw, sym_col="Gene.Symbol")
inspect("PXD062542 — motor cortex protein", prot_raw, sym_col="Gene.Symbol")

# Save inspection report
with open(out("00_input_inspection.txt"), "w") as fh:
    fh.write("SetA Pipeline — Input Inspection Report\n")
    fh.write("=" * 60 + "\n")
    fh.write("\n".join(report_lines))

print("\n[Saved] 00_input_inspection.txt")

# ─────────────────────────────────────────────────────────────────────────────
# STEP 2 — Clean and filter mRNA data
# ─────────────────────────────────────────────────────────────────────────────
sep("STEP 2 — mRNA filtering")

# ── Spinal cord mRNA ─────────────────────────────────────────────────────────
# adjPval strategy: use MIN of the two dataset-level adjusted p-values
# (at least one of GSE833 / GSE76220 must reach FDR < 0.05).
# Using the MAX would require both datasets to be independently significant,
# which is too conservative given differing sample sizes and statistical power.
print("\n--- Spinal cord mRNA ---")
sc_mrna = sc_mrna_raw[["Gene_Symbol", "mRNA_logFC_mean",
                        "mRNA_adjPval_833", "mRNA_adjPval_76220"]].copy()
sc_mrna["adjPval"] = sc_mrna[["mRNA_adjPval_833", "mRNA_adjPval_76220"]].min(axis=1)
sc_mrna = sc_mrna[["Gene_Symbol", "mRNA_logFC_mean", "adjPval"]].copy()
sc_mrna.columns = ["Gene_Symbol", "logFC", "adjPval"]
n_before = len(sc_mrna)

sc_mrna, n_sym_drop = clean_symbols(sc_mrna)
sc_mrna["logFC"]   = pd.to_numeric(sc_mrna["logFC"],   errors="coerce")
sc_mrna["adjPval"] = pd.to_numeric(sc_mrna["adjPval"], errors="coerce")
sc_mrna.dropna(subset=["logFC", "adjPval"], inplace=True)

sc_mrna = sc_mrna[
    (sc_mrna["logFC"].abs() > MRNA_FC_THRESH) &
    (sc_mrna["adjPval"] < MRNA_PADJ_THRESH)
].copy()

sc_mrna["tissue"]    = "spinal_cord"
sc_mrna["omic"]      = "mRNA"
sc_mrna["direction"] = direction(sc_mrna["logFC"])
n_after = len(sc_mrna)
n_up    = (sc_mrna["direction"] == "UP").sum()
n_down  = (sc_mrna["direction"] == "DOWN").sum()

sc_mrna.to_csv(out("SC_mRNA_filtered.csv"), index=False)
print(f"  Before filtering : {n_before} genes")
print(f"  Invalid symbols  : {n_sym_drop} dropped")
print(f"  After filtering  : {n_after} genes  (UP={n_up}, DOWN={n_down})")
print("[Saved] SC_mRNA_filtered.csv")

# ── Motor cortex mRNA ─────────────────────────────────────────────────────────
print("\n--- Motor cortex mRNA ---")
ctx_mrna = ctx_mrna_raw[["Gene_Symbol", "mRNA_logFC", "mRNA_adjPval"]].copy()
ctx_mrna.columns = ["Gene_Symbol", "logFC", "adjPval"]
n_before = len(ctx_mrna)

ctx_mrna, n_sym_drop = clean_symbols(ctx_mrna)
ctx_mrna["logFC"]   = pd.to_numeric(ctx_mrna["logFC"],   errors="coerce")
ctx_mrna["adjPval"] = pd.to_numeric(ctx_mrna["adjPval"], errors="coerce")
ctx_mrna.dropna(subset=["logFC", "adjPval"], inplace=True)

ctx_mrna = ctx_mrna[
    (ctx_mrna["logFC"].abs() > MRNA_FC_THRESH) &
    (ctx_mrna["adjPval"] < MRNA_PADJ_THRESH)
].copy()

ctx_mrna["tissue"]    = "cortex"
ctx_mrna["omic"]      = "mRNA"
ctx_mrna["direction"] = direction(ctx_mrna["logFC"])
n_after = len(ctx_mrna)
n_up    = (ctx_mrna["direction"] == "UP").sum()
n_down  = (ctx_mrna["direction"] == "DOWN").sum()

ctx_mrna.to_csv(out("CTX_mRNA_filtered.csv"), index=False)
print(f"  Before filtering : {n_before} genes")
print(f"  Invalid symbols  : {n_sym_drop} dropped")
print(f"  After filtering  : {n_after} genes  (UP={n_up}, DOWN={n_down})")
print("[Saved] CTX_mRNA_filtered.csv")

# ─────────────────────────────────────────────────────────────────────────────
# STEP 3 — Clean and filter protein data
# ─────────────────────────────────────────────────────────────────────────────
sep("STEP 3 — Protein filtering")

def filter_protein(df_raw, sym_col, logfc_col, padj_col, tissue_label, omic_label):
    """Generic protein filter for one tissue."""
    df = df_raw[[sym_col, logfc_col, padj_col]].copy()
    df.columns = ["Gene_Symbol", "logFC", "adjPval"]
    n_before = len(df)

    df, n_sym_drop = clean_symbols(df)
    df["logFC"]   = pd.to_numeric(df["logFC"],   errors="coerce")
    df["adjPval"] = pd.to_numeric(df["adjPval"], errors="coerce")
    df.dropna(subset=["logFC", "adjPval"], inplace=True)

    # Deduplicate: keep lowest adjPval per symbol
    df = df.sort_values("adjPval").drop_duplicates(subset="Gene_Symbol", keep="first")

    df = df[
        (df["logFC"].abs() > PROT_FC_THRESH) &
        (df["adjPval"] < PROT_PADJ_THRESH)
    ].copy()

    df["tissue"]    = tissue_label
    df["omic"]      = omic_label
    df["direction"] = direction(df["logFC"])
    n_after = len(df)
    n_up    = (df["direction"] == "UP").sum()
    n_down  = (df["direction"] == "DOWN").sum()

    print(f"  Before filtering : {n_before} genes")
    print(f"  Invalid symbols  : {n_sym_drop} dropped")
    print(f"  After filtering  : {n_after} genes  (UP={n_up}, DOWN={n_down})")
    return df

print("\n--- Spinal cord protein (lumbar) ---")
sc_prot = filter_protein(
    prot_raw,
    sym_col   = "Gene.Symbol",
    logfc_col = "M.lumbar_vs_C.lumbar_diff",
    padj_col  = "M.lumbar_vs_C.lumbar_p.adj",
    tissue_label = "spinal_cord",
    omic_label   = "protein",
)
sc_prot.to_csv(out("SC_protein_filtered.csv"), index=False)
print("[Saved] SC_protein_filtered.csv")

print("\n--- Motor cortex protein ---")
ctx_prot = filter_protein(
    prot_raw,
    sym_col   = "Gene.Symbol",
    logfc_col = "M.Mcx_vs_C.Mcx_diff",
    padj_col  = "M.Mcx_vs_C.Mcx_p.adj",
    tissue_label = "cortex",
    omic_label   = "protein",
)
ctx_prot.to_csv(out("CTX_protein_filtered.csv"), index=False)
print("[Saved] CTX_protein_filtered.csv")

# ─────────────────────────────────────────────────────────────────────────────
# STEP 4 — Build Set A: same-tissue multi-omic intersection
# ─────────────────────────────────────────────────────────────────────────────
sep("STEP 4 — Set A intersections")

def build_set_a(mrna_df, prot_df, tissue_prefix):
    """
    Inner join mRNA and protein on Gene_Symbol.
    Keep only concordant direction rows.
    Returns annotated DataFrame.
    """
    m = mrna_df[["Gene_Symbol", "logFC", "adjPval", "direction"]].copy()
    m.columns = [
        "Gene_Symbol",
        f"{tissue_prefix}_mRNA_logFC",
        f"{tissue_prefix}_mRNA_adjPval",
        "dir_mrna",
    ]
    p = prot_df[["Gene_Symbol", "logFC", "adjPval", "direction"]].copy()
    p.columns = [
        "Gene_Symbol",
        f"{tissue_prefix}_protein_logFC",
        f"{tissue_prefix}_protein_adjPval",
        "dir_prot",
    ]
    merged = pd.merge(m, p, on="Gene_Symbol", how="inner")

    # Direction concordance
    concordant = merged["dir_mrna"] == merged["dir_prot"]
    n_discordant = (~concordant).sum()
    if n_discordant:
        print(f"  Discarding {n_discordant} gene(s) with conflicting mRNA/protein direction")
    merged = merged[concordant].copy()

    merged["direction"]      = merged["dir_mrna"]
    merged["evidence_count"] = 2
    merged["tissue"]         = tissue_prefix.lower()
    merged.drop(columns=["dir_mrna", "dir_prot"], inplace=True)

    n_up   = (merged["direction"] == "UP").sum()
    n_down = (merged["direction"] == "DOWN").sum()
    print(f"  Set A genes: {len(merged)}  (UP={n_up}, DOWN={n_down})")
    return merged

print("\n--- Spinal cord Set A ---")
setA_sc  = build_set_a(sc_mrna,  sc_prot,  tissue_prefix="SC")
setA_sc.to_csv(out("SetA_spinal_cord.csv"), index=False)
print("[Saved] SetA_spinal_cord.csv")

print("\n--- Motor cortex Set A ---")
setA_ctx = build_set_a(ctx_mrna, ctx_prot, tissue_prefix="CTX")
setA_ctx.to_csv(out("SetA_cortex.csv"), index=False)
print("[Saved] SetA_cortex.csv")

# ─────────────────────────────────────────────────────────────────────────────
# STEP 5 — Combine and annotate full Set A
# ─────────────────────────────────────────────────────────────────────────────
sep("STEP 5 — Combine into SetA_combined")

# Rename columns for merge compatibility
sc_cols  = setA_sc.rename(columns={"direction": "direction_SC",
                                    "tissue": "tissue_SC",
                                    "evidence_count": "ec_SC"})
ctx_cols = setA_ctx.rename(columns={"direction": "direction_CTX",
                                     "tissue": "tissue_CTX",
                                     "evidence_count": "ec_CTX"})

# Outer join both tissues
combined = pd.merge(
    sc_cols.drop(columns=["tissue_SC", "ec_SC"]),
    ctx_cols.drop(columns=["tissue_CTX", "ec_CTX"]),
    on="Gene_Symbol",
    how="outer",
)

# Presence flags
combined["present_in_SC"]  = combined["Gene_Symbol"].isin(setA_sc["Gene_Symbol"])
combined["present_in_CTX"] = combined["Gene_Symbol"].isin(setA_ctx["Gene_Symbol"])

# Cross-tissue logic
both = combined["present_in_SC"] & combined["present_in_CTX"]
same_dir = combined["direction_SC"] == combined["direction_CTX"]

combined["cross_tissue"]            = False
combined["cross_tissue_discordant"] = False
combined.loc[both &  same_dir, "cross_tissue"]            = True
combined.loc[both & ~same_dir, "cross_tissue_discordant"] = True

# evidence_count: 4 if both tissues same direction, 2 otherwise
combined["evidence_count"] = 2
combined.loc[combined["cross_tissue"], "evidence_count"] = 4

# Final column order
col_order = [
    "Gene_Symbol",
    "direction_SC", "direction_CTX",
    "SC_mRNA_logFC", "SC_mRNA_adjPval",
    "SC_protein_logFC", "SC_protein_adjPval",
    "CTX_mRNA_logFC", "CTX_mRNA_adjPval",
    "CTX_protein_logFC", "CTX_protein_adjPval",
    "evidence_count", "cross_tissue", "cross_tissue_discordant",
    "present_in_SC", "present_in_CTX",
]
combined = combined[col_order]
combined.sort_values(["evidence_count", "Gene_Symbol"],
                     ascending=[False, True], inplace=True)
combined.reset_index(drop=True, inplace=True)

combined.to_csv(out("SetA_combined.csv"), index=False)

n_total    = len(combined)
n_sc_only  = (combined["present_in_SC"]  & ~combined["present_in_CTX"]).sum()
n_ctx_only = (combined["present_in_CTX"] & ~combined["present_in_SC"]).sum()
n_both_same= combined["cross_tissue"].sum()
n_both_disc= combined["cross_tissue_discordant"].sum()

print(f"  Total unique genes in Set A : {n_total}")
print(f"  Spinal cord only            : {n_sc_only}")
print(f"  Motor cortex only           : {n_ctx_only}")
print(f"  Both tissues, same direction: {n_both_same}  (evidence_count=4)")
print(f"  Both tissues, discordant    : {n_both_disc}")
print("[Saved] SetA_combined.csv")

# ─────────────────────────────────────────────────────────────────────────────
# STEP 6 — Prepare STRING input files
# ─────────────────────────────────────────────────────────────────────────────
sep("STEP 6 — STRING input files")

# Plain gene list
string_genes = combined["Gene_Symbol"].unique()
with open(out("SetA_STRING_input.txt"), "w") as fh:
    fh.write("\n".join(sorted(string_genes)) + "\n")
print(f"[Saved] SetA_STRING_input.txt  ({len(string_genes)} genes)")

# Ranked version
ranked = combined[["Gene_Symbol", "evidence_count",
                    "present_in_SC", "present_in_CTX"]].copy()
ranked.sort_values(["evidence_count", "Gene_Symbol"],
                   ascending=[False, True], inplace=True)
ranked.reset_index(drop=True, inplace=True)
ranked.to_csv(out("SetA_STRING_input_ranked.csv"), index=False)
print("[Saved] SetA_STRING_input_ranked.csv (pre-annotation)")

# ─────────────────────────────────────────────────────────────────────────────
# STEP 7 — Biological annotation
# ─────────────────────────────────────────────────────────────────────────────
sep("STEP 7 — Biological annotation")

gene_set = set(combined["Gene_Symbol"])

def overlap_report(label, ref_set):
    hits = sorted(gene_set & ref_set)
    print(f"\n  {label}:")
    print(f"    Overlap: {len(hits)} gene(s)")
    if hits:
        print(f"    Genes  : {', '.join(hits)}")
    return hits

hits_als  = overlap_report("Known ALS genes",           KNOWN_ALS)
hits_ups  = overlap_report("UPS-related genes",         UPS_RELATED)
hits_disc = overlap_report("Discordant_LOW (prev.)",    DISCORDANT_LOW)
hits_sg   = overlap_report("Stress granule interactome",STRESS_GRANULE)

# Add annotation columns to ranked file
ranked["known_ALS_gene"]             = ranked["Gene_Symbol"].isin(KNOWN_ALS)
ranked["UPS_related"]                = ranked["Gene_Symbol"].isin(UPS_RELATED)
ranked["discordant_LOW"]             = ranked["Gene_Symbol"].isin(DISCORDANT_LOW)
ranked["stress_granule_interactome"] = ranked["Gene_Symbol"].isin(STRESS_GRANULE)
ranked["n_annotation_categories"]   = (
    ranked["known_ALS_gene"].astype(int) +
    ranked["UPS_related"].astype(int) +
    ranked["discordant_LOW"].astype(int) +
    ranked["stress_granule_interactome"].astype(int)
)

ranked.to_csv(out("SetA_STRING_input_ranked.csv"), index=False)
print("\n[Saved] SetA_STRING_input_ranked.csv (with annotations)")

# ─────────────────────────────────────────────────────────────────────────────
# STEP 8 — Summary report
# ─────────────────────────────────────────────────────────────────────────────
sep("STEP 8 — Summary report")

cross_tissue_genes = sorted(
    combined.loc[combined["cross_tissue"], "Gene_Symbol"].tolist()
)
high_priority_genes = sorted(
    ranked.loc[ranked["n_annotation_categories"] >= 2, "Gene_Symbol"].tolist()
)

summary_lines = [
    "SetA Pipeline — Summary Report",
    "=" * 60,
    "",
    "── Filtered inputs ──────────────────────────────────────────",
    f"  SC mRNA  genes passing filters : {len(sc_mrna)}",
    f"  CTX mRNA genes passing filters : {len(ctx_mrna)}",
    f"  SC protein genes passing filters: {len(sc_prot)}",
    f"  CTX protein genes passing filters: {len(ctx_prot)}",
    "",
    "── Set A per-tissue ─────────────────────────────────────────",
    f"  SetA spinal cord : {len(setA_sc)} genes  "
    f"(UP={( setA_sc['direction']=='UP').sum()}, DOWN={( setA_sc['direction']=='DOWN').sum()})",
    f"  SetA motor cortex: {len(setA_ctx)} genes  "
    f"(UP={(setA_ctx['direction']=='UP').sum()}, DOWN={(setA_ctx['direction']=='DOWN').sum()})",
    "",
    "── SetA combined ────────────────────────────────────────────",
    f"  Total unique genes              : {n_total}",
    f"  Spinal cord only                : {n_sc_only}",
    f"  Motor cortex only               : {n_ctx_only}",
    f"  Both tissues, same direction    : {n_both_same}  (evidence_count=4)",
    f"  Both tissues, discordant        : {n_both_disc}",
    "",
    "── Biological annotation overlaps ───────────────────────────",
    f"  Known ALS genes                 : {len(hits_als)} — {', '.join(hits_als) if hits_als else 'none'}",
    f"  UPS-related genes               : {len(hits_ups)} — {', '.join(hits_ups) if hits_ups else 'none'}",
    f"  Discordant_LOW (prev. analysis) : {len(hits_disc)} — {', '.join(hits_disc) if hits_disc else 'none'}",
    f"  Stress granule interactome      : {len(hits_sg)} — {', '.join(hits_sg) if hits_sg else 'none'}",
    "",
    "── Cross-tissue confirmed genes (evidence_count = 4) ────────",
    f"  N = {len(cross_tissue_genes)}",
]
if cross_tissue_genes:
    summary_lines.append("  " + ", ".join(cross_tissue_genes))
else:
    summary_lines.append("  (none)")

summary_lines += [
    "",
    "── Highest-priority STRING seeds (n_annotation_categories ≥ 2) ─",
    f"  N = {len(high_priority_genes)}",
]
if high_priority_genes:
    # Wrap long list for readability
    wrapped = textwrap.fill(", ".join(high_priority_genes), width=70,
                            initial_indent="  ", subsequent_indent="  ")
    summary_lines.append(wrapped)
else:
    summary_lines.append("  (none)")

summary_text = "\n".join(summary_lines) + "\n"
print("\n" + summary_text)

with open(out("SetA_summary.txt"), "w", encoding="utf-8") as fh:
    fh.write(summary_text)
print("[Saved] SetA_summary.txt")

sep("Pipeline complete")
print(f"\nAll outputs in: {OUT_DIR}\n")
