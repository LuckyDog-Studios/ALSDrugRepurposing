"""
mRNA-Protein Correlation Pipeline for ALS Drug Repurposing
===========================================================
Identifies discordant genes (high mRNA / low protein or vice versa)
as candidates for UPS dysfunction in ALS.

Differential expression has already been computed — this pipeline
performs cleaning, correlation, and prioritisation only.

All outputs are saved to the same directory as this script.
"""

import subprocess, sys

# ── install dependencies ──────────────────────────────────────────────────────
for pkg in ["pandas", "numpy", "scipy", "matplotlib", "seaborn", "mygene"]:
    subprocess.check_call([sys.executable, "-m", "pip", "install", pkg, "--quiet"])

import os, re, textwrap, io
from pathlib import Path

# Force UTF-8 output on Windows to avoid cp1252 encoding errors
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
sys.stderr = io.TextIOWrapper(sys.stderr.buffer, encoding='utf-8', errors='replace')

import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

# ── path setup ────────────────────────────────────────────────────────────────
SCRIPT_DIR  = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parents[2]          # …/MarbleProject
MRNA_DIR    = PROJECT_ROOT / "differential_expression" / "gene expression"
PROT_DIR    = PROJECT_ROOT / "differential_expression" / "protein expression"
OUT_DIR     = SCRIPT_DIR                      # all outputs land here

def out(fname):
    return OUT_DIR / fname

# ── file manifest (actual filenames on disk) ──────────────────────────────────
FILES = {
    # mRNA – included
    "GSE833":   MRNA_DIR / "GSE833_DE_Results.csv",
    "GSE76220": MRNA_DIR / "GSE76220_DE_Results.csv",
    "GSE124439": MRNA_DIR / "GSE124439_ALS_vs_Control_Results.csv",
    # mRNA – skipped
    "GSE118336_Homo":   MRNA_DIR / "GSE118336_annotated_de_FUS_H517D_Mutant_vs_Control.csv",
    "GSE118336_Hetero": MRNA_DIR / "GSE118336_annotated_de_FUS_Heterozygous_vs_Control.csv",
    "GSE112676":        MRNA_DIR / "GSE112676_DE_results.csv",
    # protein – main correlation
    "PXD067060_CTX": PROT_DIR / "PXD067060_DE_CTX_TP7_vs_TP1.csv",
    "PXD067060_SPC": PROT_DIR / "PXD067060_DE_SPC_TP7_vs_TP1.csv",
    "PXD067060_ALL": PROT_DIR / "PXD067060_DE_results.csv",
    # protein – network expansion
    "PXD065424_ARS":       PROT_DIR / "PXD065424_DE_ARS_C9ALS_vs_Control.csv",
    "PXD065424_CTRL":      PROT_DIR / "PXD065424_DE_CTRL_C9ALS_vs_Control.csv",
    "PXD065424_RNase_ARS":  PROT_DIR / "PXD065424_DE_RNase_ARS_C9ALS_vs_Control.csv",
    "PXD065424_RNase_CTRL": PROT_DIR / "PXD065424_DE_RNase_CTRL_C9ALS_vs_Control.csv",
    "PXD065424_ALL": PROT_DIR / "PXD065424_DE_results.csv",
}

print("\n" + "="*70)
print("Expected input files:")
print("="*70)
missing = []
for name, path in FILES.items():
    status = "OK" if path.exists() else "MISSING"
    print(f"  [{status}]  {name}: {path.name}")
    if not path.exists():
        missing.append(name)

if missing:
    print(f"\nERROR: {len(missing)} file(s) not found: {missing}")
    print("Please verify paths and rerun.")
    sys.exit(1)
else:
    print("\nAll files found. Proceeding.\n")


# ─────────────────────────────────────────────────────────────────────────────
# STEP 1  Inspect all input files
# ─────────────────────────────────────────────────────────────────────────────
print("="*70)
print("STEP 1: Inspecting all input files")
print("="*70)

report_lines = ["mRNA-Protein Correlation Pipeline — Inspection Report", "="*70, ""]

def inspect(label, path):
    """Load a CSV and return inspection summary + the dataframe."""
    df = pd.read_csv(path)
    lines = [
        f"File: {label} ({path.name})",
        f"  Shape:   {df.shape[0]} rows × {df.shape[1]} columns",
        f"  Columns: {df.columns.tolist()}",
        "  First 3 rows:",
    ]
    for _, row in df.head(3).iterrows():
        lines.append(f"    {row.to_dict()}")
    lines.append("  Missing values per column:")
    for col, n in df.isnull().sum().items():
        lines.append(f"    {col}: {n}")
    for col in ["condition", "comparison", "tissue"]:
        if col in df.columns:
            lines.append(f"  Unique '{col}' values: {df[col].unique().tolist()}")
    lines.append("")
    print("\n".join(lines))
    return df, lines

raw_data = {}
for label, path in FILES.items():
    df, lines = inspect(label, path)
    raw_data[label] = df
    report_lines.extend(lines)

with open(out("00_inspection_report.txt"), "w", encoding="utf-8") as fh:
    fh.write("\n".join(report_lines))
print("Saved: 00_inspection_report.txt\n")


# ─────────────────────────────────────────────────────────────────────────────
# STEP 2  Clean and standardise mRNA DEG files
# ─────────────────────────────────────────────────────────────────────────────
print("="*70)
print("STEP 2: Cleaning and standardising mRNA DEG files")
print("="*70)

# Valid HGNC symbol: starts with a letter, contains only letters/digits/hyphens
VALID_SYMBOL = re.compile(r'^[A-Za-z][A-Za-z0-9\-]*$')

def is_valid_symbol(s):
    if pd.isna(s) or str(s).strip() == "":
        return False
    return bool(VALID_SYMBOL.match(str(s).strip()))

def clean_mrna(df, label, gene_col, logfc_col, adjpval_col):
    """
    Standardise a mRNA DEG dataframe to Gene_Symbol / mRNA_logFC / mRNA_adjPval.
    Applies HGNC symbol validation, removes duplicates (keep min adj.P.Val).
    """
    n_raw = len(df)

    # Rename gene column to Gene_Symbol
    df = df.rename(columns={gene_col: "Gene_Symbol",
                             logfc_col: "mRNA_logFC",
                             adjpval_col: "mRNA_adjPval"})

    # Keep only the three columns we need
    df = df[["Gene_Symbol", "mRNA_logFC", "mRNA_adjPval"]].copy()

    # Strip whitespace from gene symbols
    df["Gene_Symbol"] = df["Gene_Symbol"].astype(str).str.strip()

    # Drop rows where logFC or adj.P.Val is missing
    n_before_na = len(df)
    df = df.dropna(subset=["mRNA_logFC", "mRNA_adjPval"])
    n_dropped_na = n_before_na - len(df)

    # Filter to valid HGNC symbols
    valid_mask = df["Gene_Symbol"].apply(is_valid_symbol)
    n_invalid = (~valid_mask).sum()
    df = df[valid_mask].copy()

    # Deduplicate: keep row with lowest adj.P.Val per gene
    df = df.sort_values("mRNA_adjPval").drop_duplicates(subset="Gene_Symbol", keep="first")
    df = df.reset_index(drop=True)

    print(f"  {label}: {n_raw} raw rows → "
          f"{n_dropped_na} dropped (missing logFC/adjPval) → "
          f"{n_invalid} dropped (invalid symbol) → "
          f"{len(df)} genes remain")
    return df


# GSE833: Gene_Symbol is the last column
df833 = raw_data["GSE833"].copy()
df833_clean = clean_mrna(df833, "GSE833",
                          gene_col="Gene_Symbol",
                          logfc_col="logFC",
                          adjpval_col="adj.P.Val")
df833_clean.to_csv(out("GSE833_clean.csv"), index=False)
print("  Saved: GSE833_clean.csv")

# GSE76220: gene symbols in unnamed first column (row-index from R)
df76220 = raw_data["GSE76220"].copy()
# The unnamed column holds gene symbols
unnamed_col = df76220.columns[0]   # 'Unnamed: 0'
df76220_clean = clean_mrna(df76220, "GSE76220",
                            gene_col=unnamed_col,
                            logfc_col="logFC",
                            adjpval_col="adj.P.Val")
df76220_clean.to_csv(out("GSE76220_clean.csv"), index=False)
print("  Saved: GSE76220_clean.csv")

# GSE124439: standard layout; raw p-value column is 'PValue' (not used here)
df124439 = raw_data["GSE124439"].copy()
df124439_clean = clean_mrna(df124439, "GSE124439",
                             gene_col="Gene_Symbol",
                             logfc_col="logFC",
                             adjpval_col="adj.P.Val")
df124439_clean.to_csv(out("GSE124439_clean.csv"), index=False)
print("  Saved: GSE124439_clean.csv\n")


# ─────────────────────────────────────────────────────────────────────────────
# STEP 3  Merge spinal cord mRNA datasets (GSE833 + GSE76220)
# ─────────────────────────────────────────────────────────────────────────────
print("="*70)
print("STEP 3: Merging spinal cord mRNA datasets (GSE833 + GSE76220)")
print("="*70)

# Inner join on Gene_Symbol → only genes present in BOTH datasets
sc_mrna = pd.merge(
    df833_clean.rename(columns={"mRNA_logFC": "mRNA_logFC_833",
                                "mRNA_adjPval": "mRNA_adjPval_833"}),
    df76220_clean.rename(columns={"mRNA_logFC": "mRNA_logFC_76220",
                                  "mRNA_adjPval": "mRNA_adjPval_76220"}),
    on="Gene_Symbol", how="inner"
)

# Mean logFC across both datasets
sc_mrna["mRNA_logFC_mean"] = (sc_mrna["mRNA_logFC_833"] + sc_mrna["mRNA_logFC_76220"]) / 2

# Most conservative (highest) adj.P.Val
sc_mrna["mRNA_adjPval_max"] = sc_mrna[["mRNA_adjPval_833", "mRNA_adjPval_76220"]].max(axis=1)

sc_mrna = sc_mrna[["Gene_Symbol", "mRNA_logFC_mean", "mRNA_adjPval_max",
                    "mRNA_logFC_833", "mRNA_logFC_76220",
                    "mRNA_adjPval_833", "mRNA_adjPval_76220"]]

sc_mrna.to_csv(out("spinal_cord_mRNA_merged.csv"), index=False)
print(f"  Spinal cord merged mRNA: {len(sc_mrna)} genes (inner join GSE833 ∩ GSE76220)")
print("  Saved: spinal_cord_mRNA_merged.csv\n")


# ─────────────────────────────────────────────────────────────────────────────
# STEP 4  Prepare protein DEG files (PXD067060)
# ─────────────────────────────────────────────────────────────────────────────
print("="*70)
print("STEP 4: Preparing PXD067060 protein DEG files")
print("="*70)

def clean_protein_pxd067060(df, label):
    """
    PXD067060 already has Gene_Symbol. Standardise columns and deduplicate.
    """
    print(f"  {label} columns: {df.columns.tolist()}")

    # Both CTX and SPC files have Gene_Symbol, logFC, adj.P.Val
    df = df[["Gene_Symbol", "logFC", "adj.P.Val"]].copy()
    df = df.rename(columns={"logFC": "protein_logFC", "adj.P.Val": "protein_adjPval"})
    df["Gene_Symbol"] = df["Gene_Symbol"].astype(str).str.strip()

    n_raw = len(df)
    df = df.dropna(subset=["protein_logFC", "protein_adjPval"])
    df = df[df["Gene_Symbol"].apply(is_valid_symbol)]
    df = df.sort_values("protein_adjPval").drop_duplicates(subset="Gene_Symbol", keep="first")
    df = df.reset_index(drop=True)

    print(f"  {label}: {n_raw} raw rows → {len(df)} genes after cleaning")
    return df

df067_ctx = raw_data["PXD067060_CTX"].copy()
df067_spc = raw_data["PXD067060_SPC"].copy()

prot_ctx = clean_protein_pxd067060(df067_ctx, "PXD067060_CTX")
prot_spc = clean_protein_pxd067060(df067_spc, "PXD067060_SPC")

prot_ctx.to_csv(out("PXD067060_CTX_clean.csv"), index=False)
prot_spc.to_csv(out("PXD067060_SPC_clean.csv"), index=False)
print("  Saved: PXD067060_CTX_clean.csv, PXD067060_SPC_clean.csv\n")


# ─────────────────────────────────────────────────────────────────────────────
# Helpers for Steps 5 & 6
# ─────────────────────────────────────────────────────────────────────────────

QUADRANT_COLORS = {
    "concordant_up":   "#2196F3",   # blue
    "concordant_down": "#9C27B0",   # purple
    "discordant_HIGH": "#FF9800",   # orange  (high mRNA, low protein)
    "discordant_LOW":  "#F44336",   # red     (low mRNA, high protein — UPS candidates)
    "ambiguous":       "#BDBDBD",   # grey
}

def classify_quadrant(mrna_lfc, prot_lfc, thresh=0.5):
    if mrna_lfc > thresh and prot_lfc > thresh:
        return "concordant_up"
    elif mrna_lfc < -thresh and prot_lfc < -thresh:
        return "concordant_down"
    elif mrna_lfc > thresh and prot_lfc < -thresh:
        return "discordant_HIGH"
    elif mrna_lfc < -thresh and prot_lfc > thresh:
        return "discordant_LOW"
    else:
        return "ambiguous"


def run_correlation(mrna_df, prot_df, mrna_lfc_col,
                    tissue_label, plot_fname, results_fname):
    """
    Inner-join mRNA and protein data, compute correlation and quadrants,
    generate scatter plot, save results CSV.
    Returns the merged dataframe and correlation stats.
    """
    merged = pd.merge(mrna_df, prot_df, on="Gene_Symbol", how="inner")
    print(f"  {tissue_label}: {len(merged)} genes in overlap")

    # Drop any rows where either logFC is NaN (defensive)
    merged = merged.dropna(subset=[mrna_lfc_col, "protein_logFC"])
    n = len(merged)

    # Discordance score = protein_logFC − mRNA_logFC
    merged["discordance_score"] = merged["protein_logFC"] - merged[mrna_lfc_col]

    # Quadrant classification
    merged["quadrant"] = merged.apply(
        lambda r: classify_quadrant(r[mrna_lfc_col], r["protein_logFC"]), axis=1
    )

    # Pearson and Spearman correlations
    pearson_r, pearson_p = pearsonr(merged[mrna_lfc_col], merged["protein_logFC"])
    spearman_r, spearman_p = spearmanr(merged[mrna_lfc_col], merged["protein_logFC"])
    print(f"  {tissue_label}: Pearson r={pearson_r:.4f} (p={pearson_p:.3e}), "
          f"Spearman r={spearman_r:.4f} (p={spearman_p:.3e})")

    # Quadrant counts
    for q, cnt in merged["quadrant"].value_counts().items():
        print(f"    {q}: {cnt}")

    # ── scatter plot ──────────────────────────────────────────────────────────
    fig, ax = plt.subplots(figsize=(9, 8))

    for quad, color in QUADRANT_COLORS.items():
        mask = merged["quadrant"] == quad
        ax.scatter(
            merged.loc[mask, mrna_lfc_col],
            merged.loc[mask, "protein_logFC"],
            c=color, label=f"{quad} (n={mask.sum()})",
            alpha=0.65, s=25, linewidths=0
        )

    # Diagonal reference line y = x
    lim = max(
        abs(merged[mrna_lfc_col]).max(),
        abs(merged["protein_logFC"]).max()
    ) * 1.1
    ax.plot([-lim, lim], [-lim, lim], "k--", linewidth=0.8, alpha=0.5, label="y = x")

    # Label 10 most discordant genes (largest |discordance_score|)
    top10 = merged.nlargest(10, "discordance_score", keep="first")
    bot10 = merged.nsmallest(10, "discordance_score", keep="first")
    label_genes = pd.concat([top10, bot10]).drop_duplicates(subset="Gene_Symbol")
    for _, row in label_genes.iterrows():
        ax.annotate(
            row["Gene_Symbol"],
            xy=(row[mrna_lfc_col], row["protein_logFC"]),
            xytext=(5, 5), textcoords="offset points",
            fontsize=7, alpha=0.85,
        )

    ax.set_xlabel("mRNA log₂FC", fontsize=12)
    ax.set_ylabel("Protein log₂FC", fontsize=12)
    ax.set_title(
        f"{tissue_label} mRNA–Protein Correlation\n"
        f"Pearson r = {pearson_r:.3f}  |  Spearman r = {spearman_r:.3f}  |  n = {n}",
        fontsize=13
    )
    ax.legend(fontsize=8, loc="upper left")
    ax.axhline(0, color="grey", linewidth=0.5, linestyle=":")
    ax.axvline(0, color="grey", linewidth=0.5, linestyle=":")
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    plt.tight_layout()
    plt.savefig(out(plot_fname), dpi=150)
    plt.close()
    print(f"  Saved plot: {plot_fname}")

    # ── results CSV ──────────────────────────────────────────────────────────
    # Build output columns based on whether this is spinal cord (has adjPval_max) or cortex
    if "mRNA_adjPval_max" in merged.columns:
        adjpval_col = "mRNA_adjPval_max"
    else:
        adjpval_col = "mRNA_adjPval"

    out_df = merged[["Gene_Symbol", mrna_lfc_col, adjpval_col,
                      "protein_logFC", "protein_adjPval",
                      "discordance_score", "quadrant"]].copy()
    # Standardise column names
    out_df = out_df.rename(columns={mrna_lfc_col: "mRNA_logFC_mean",
                                    adjpval_col: "mRNA_adjPval_max"})
    out_df.to_csv(out(results_fname), index=False)
    print(f"  Saved results: {results_fname}\n")

    stats = {
        "tissue": tissue_label,
        "n_genes": n,
        "pearson_r": pearson_r,
        "pearson_p": pearson_p,
        "spearman_r": spearman_r,
        "spearman_p": spearman_p,
        "quadrant_counts": merged["quadrant"].value_counts().to_dict(),
    }
    return out_df, stats


# ─────────────────────────────────────────────────────────────────────────────
# STEP 5  mRNA–protein correlation — spinal cord
# ─────────────────────────────────────────────────────────────────────────────
print("="*70)
print("STEP 5: mRNA–protein correlation — spinal cord")
print("="*70)

sc_results, sc_stats = run_correlation(
    mrna_df=sc_mrna,
    prot_df=prot_spc,
    mrna_lfc_col="mRNA_logFC_mean",
    tissue_label="Spinal Cord",
    plot_fname="spinal_cord_correlation_plot.png",
    results_fname="spinal_cord_correlation_results.csv",
)


# ─────────────────────────────────────────────────────────────────────────────
# STEP 6  mRNA–protein correlation — cortex
# ─────────────────────────────────────────────────────────────────────────────
print("="*70)
print("STEP 6: mRNA–protein correlation — cortex")
print("="*70)

# For cortex, single dataset: rename for consistent internal use
ctx_mrna = df124439_clean.rename(columns={"mRNA_logFC": "mRNA_logFC_mean",
                                          "mRNA_adjPval": "mRNA_adjPval"})

ctx_results, ctx_stats = run_correlation(
    mrna_df=ctx_mrna,
    prot_df=prot_ctx,
    mrna_lfc_col="mRNA_logFC_mean",
    tissue_label="Cortex",
    plot_fname="cortex_correlation_plot.png",
    results_fname="cortex_correlation_results.csv",
)


# ─────────────────────────────────────────────────────────────────────────────
# STEP 7  Cross-tissue prioritisation for STRING
# ─────────────────────────────────────────────────────────────────────────────
print("="*70)
print("STEP 7: Cross-tissue prioritisation for STRING")
print("="*70)

# Extract discordant_LOW genes per tissue
sc_low = sc_results[sc_results["quadrant"] == "discordant_LOW"][
    ["Gene_Symbol", "mRNA_logFC_mean", "protein_logFC", "discordance_score"]
].rename(columns={
    "mRNA_logFC_mean": "spinal_mRNA_logFC",
    "protein_logFC": "spinal_protein_logFC",
    "discordance_score": "spinal_discordance",
})

ctx_low = ctx_results[ctx_results["quadrant"] == "discordant_LOW"][
    ["Gene_Symbol", "mRNA_logFC_mean", "protein_logFC", "discordance_score"]
].rename(columns={
    "mRNA_logFC_mean": "cortex_mRNA_logFC",
    "protein_logFC": "cortex_protein_logFC",
    "discordance_score": "cortex_discordance",
})

print(f"  discordant_LOW in spinal cord: {len(sc_low)}")
print(f"  discordant_LOW in cortex:      {len(ctx_low)}")

# ── outer join to combine all discordant_LOW genes across tissues ─────────────
all_low = pd.merge(sc_low, ctx_low, on="Gene_Symbol", how="outer")
all_low["present_in_spinal"] = all_low["spinal_mRNA_logFC"].notna()
all_low["present_in_cortex"]  = all_low["cortex_mRNA_logFC"].notna()

# High-confidence: discordant_LOW in BOTH tissues
high_conf = all_low[all_low["present_in_spinal"] & all_low["present_in_cortex"]].copy()

# Extended: discordant_LOW in at least ONE tissue (i.e., the full outer join)
extended = all_low.copy()

print(f"  High-confidence (both tissues): {len(high_conf)} genes")
print(f"  Extended (either tissue):       {len(extended)} genes")

high_conf.to_csv(out("ALS_Ubiquitome_STRING_input_highconfidence.csv"), index=False)
extended.to_csv(out("ALS_Ubiquitome_STRING_input_extended.csv"), index=False)
print("  Saved: ALS_Ubiquitome_STRING_input_highconfidence.csv")
print("  Saved: ALS_Ubiquitome_STRING_input_extended.csv\n")


# ─────────────────────────────────────────────────────────────────────────────
# STEP 8  PXD065424 network expansion prep
# ─────────────────────────────────────────────────────────────────────────────
print("="*70)
print("STEP 8: PXD065424 network expansion (G3BP1 interactome)")
print("="*70)

df_ars      = raw_data["PXD065424_ARS"].copy()
df_rnase_ars = raw_data["PXD065424_RNase_ARS"].copy()
df_ctrl     = raw_data["PXD065424_CTRL"].copy()
df_rnase_ctrl = raw_data["PXD065424_RNase_CTRL"].copy()

def filter_significant(df, label):
    """Keep rows where significant == TRUE (handles bool or string variants)."""
    sig_col = "significant"
    # Normalise: True / 'TRUE' / 'True' → True
    df[sig_col] = df[sig_col].astype(str).str.upper().map({"TRUE": True, "FALSE": False})
    sig_df = df[df[sig_col] == True][["Gene_Symbol", "logFC", "adj.P.Val"]].copy()
    sig_df["Gene_Symbol"] = sig_df["Gene_Symbol"].astype(str).str.strip()
    sig_df = sig_df[sig_df["Gene_Symbol"].apply(is_valid_symbol)]
    sig_df = sig_df.dropna(subset=["logFC", "adj.P.Val"])
    sig_df = sig_df.sort_values("adj.P.Val").drop_duplicates(subset="Gene_Symbol", keep="first")
    print(f"  {label}: {len(sig_df)} significant genes")
    return sig_df

sig_ars      = filter_significant(df_ars,      "ARS_C9ALS_vs_Control")
sig_rnase_ars = filter_significant(df_rnase_ars, "RNase_ARS_C9ALS_vs_Control")

# Intersection: significant in BOTH ARS conditions
# = RNA-independent stress granule interactors (highest confidence)
common_ars_genes = set(sig_ars["Gene_Symbol"]) & set(sig_rnase_ars["Gene_Symbol"])
print(f"\n  Genes significant in BOTH ARS conditions: {len(common_ars_genes)}")

hc_interactome = pd.merge(
    sig_ars[["Gene_Symbol", "logFC", "adj.P.Val"]].rename(
        columns={"logFC": "logFC_ARS", "adj.P.Val": "adj.P.Val_ARS"}),
    sig_rnase_ars[["Gene_Symbol", "logFC", "adj.P.Val"]].rename(
        columns={"logFC": "logFC_RNase_ARS", "adj.P.Val": "adj.P.Val_RNase_ARS"}),
    on="Gene_Symbol", how="inner"
)
hc_interactome.to_csv(out("PXD065424_highconfidence_interactome.csv"), index=False)
print("  Saved: PXD065424_highconfidence_interactome.csv")

# Overlap with STRING high-confidence list
if len(high_conf) > 0:
    overlap = set(hc_interactome["Gene_Symbol"]) & set(high_conf["Gene_Symbol"])
    print(f"\n  Overlap (STRING high-conf ∩ PXD065424 interactome): {len(overlap)} genes")
    if overlap:
        print(f"  Priority nodes: {sorted(overlap)}")

    priority_genes = set(hc_interactome["Gene_Symbol"]) & set(extended["Gene_Symbol"])
    priority_df = pd.merge(
        hc_interactome,
        extended[["Gene_Symbol", "spinal_mRNA_logFC", "spinal_protein_logFC",
                  "spinal_discordance", "cortex_mRNA_logFC", "cortex_protein_logFC",
                  "cortex_discordance", "present_in_spinal", "present_in_cortex"]],
        on="Gene_Symbol", how="inner"
    )
    priority_df.to_csv(out("ALS_Ubiquitome_priority_nodes.csv"), index=False)
    print(f"  Saved: ALS_Ubiquitome_priority_nodes.csv ({len(priority_df)} genes)")
else:
    # No STRING high-conf genes → priority nodes = interactome ∩ extended list
    priority_df = pd.merge(hc_interactome,
                           extended[["Gene_Symbol", "spinal_mRNA_logFC",
                                     "spinal_protein_logFC", "spinal_discordance",
                                     "cortex_mRNA_logFC", "cortex_protein_logFC",
                                     "cortex_discordance", "present_in_spinal",
                                     "present_in_cortex"]],
                           on="Gene_Symbol", how="inner")
    priority_df.to_csv(out("ALS_Ubiquitome_priority_nodes.csv"), index=False)
    print(f"  Saved: ALS_Ubiquitome_priority_nodes.csv ({len(priority_df)} genes)")

print()


# ─────────────────────────────────────────────────────────────────────────────
# STEP 9  Summary report
# ─────────────────────────────────────────────────────────────────────────────
print("="*70)
print("STEP 9: Summary report")
print("="*70)

def fmt_qcounts(qdict):
    return "  |  ".join(f"{k}: {v}" for k, v in sorted(qdict.items()))

summary = textwrap.dedent(f"""
mRNA–Protein Correlation Pipeline — Summary Report
{"="*70}

▌ SPINAL CORD (GSE833 + GSE76220 mRNA  ×  PXD067060_SPC protein)
  Genes analysed:      {sc_stats['n_genes']}
  Pearson r:           {sc_stats['pearson_r']:.4f}  (p = {sc_stats['pearson_p']:.3e})
  Spearman r:          {sc_stats['spearman_r']:.4f}  (p = {sc_stats['spearman_p']:.3e})
  Quadrant breakdown:  {fmt_qcounts(sc_stats['quadrant_counts'])}

▌ CORTEX (GSE124439 mRNA  ×  PXD067060_CTX protein)
  Genes analysed:      {ctx_stats['n_genes']}
  Pearson r:           {ctx_stats['pearson_r']:.4f}  (p = {ctx_stats['pearson_p']:.3e})
  Spearman r:          {ctx_stats['spearman_r']:.4f}  (p = {ctx_stats['spearman_p']:.3e})
  Quadrant breakdown:  {fmt_qcounts(ctx_stats['quadrant_counts'])}

▌ STRING INPUT LISTS
  High-confidence (discordant_LOW in BOTH tissues): {len(high_conf)} genes
  Extended (discordant_LOW in EITHER tissue):       {len(extended)} genes

▌ PXD065424 NETWORK EXPANSION
  Significant in ARS condition:                     {len(sig_ars)} genes
  Significant in RNase+ARS condition:               {len(sig_rnase_ars)} genes
  High-confidence interactome (both ARS):           {len(hc_interactome)} genes

▌ PRIORITY NODES (interactome ∩ discordant_LOW):    {len(priority_df)} genes
{"="*70}

mRNA files cleaned:
  GSE833   → {len(df833_clean)} genes
  GSE76220 → {len(df76220_clean)} genes
  GSE124439 → {len(df124439_clean)} genes
  Spinal cord merged (inner join): {len(sc_mrna)} genes

Protein files cleaned:
  PXD067060 CTX → {len(prot_ctx)} genes
  PXD067060 SPC → {len(prot_spc)} genes

Skipped datasets:
  GSE118336 (iPSC FUS lines)  — skipped per pipeline design
  GSE112676 (blood)           — wrong tissue, excluded

Output files in {OUT_DIR}:
  00_inspection_report.txt
  GSE833_clean.csv
  GSE76220_clean.csv
  GSE124439_clean.csv
  spinal_cord_mRNA_merged.csv
  PXD067060_CTX_clean.csv
  PXD067060_SPC_clean.csv
  spinal_cord_correlation_results.csv
  spinal_cord_correlation_plot.png
  cortex_correlation_results.csv
  cortex_correlation_plot.png
  ALS_Ubiquitome_STRING_input_highconfidence.csv
  ALS_Ubiquitome_STRING_input_extended.csv
  PXD065424_highconfidence_interactome.csv
  ALS_Ubiquitome_priority_nodes.csv
  09_summary_report.txt
""")

print(summary)
with open(out("09_summary_report.txt"), "w", encoding="utf-8") as fh:
    fh.write(summary)
print("Saved: 09_summary_report.txt")
print("\nPipeline complete.")
