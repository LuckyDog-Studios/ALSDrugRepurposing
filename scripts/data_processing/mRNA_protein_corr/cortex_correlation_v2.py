"""
Cortex mRNA-Protein Correlation — v2 (PXD062542 replacement)
=============================================================
Replaces the invalid PXD067060 cortex correlation (longitudinal TP7 vs TP1)
with PXD062542, a true ALS vs control TMT proteomics dataset from post-mortem
motor cortex and spinal cord (Pan et al., Scientific Reports 2025).

Runs Steps 1-4 only. Does NOT touch any spinal cord output files.
"""

import subprocess, sys, io

for pkg in ["pandas", "numpy", "scipy", "matplotlib", "seaborn", "openpyxl"]:
    subprocess.check_call([sys.executable, "-m", "pip", "install", pkg, "--quiet"],
                          stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

import re
import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

# Force UTF-8 output on Windows
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")
sys.stderr = io.TextIOWrapper(sys.stderr.buffer, encoding="utf-8", errors="replace")

# ── paths (all relative to this script's directory) ───────────────────────────
WORK_DIR = __file__  # resolve below
from pathlib import Path
WORK_DIR = Path(__file__).resolve().parent

EXCEL_FILE  = WORK_DIR / "41598_2025_11466_MOESM1_ESM.xlsx"
MRNA_FILE   = WORK_DIR / "GSE124439_clean.csv"
SC_RESULTS  = WORK_DIR / "spinal_cord_correlation_results.csv"
INTERACTOME = WORK_DIR / "archive_v1" / "PXD065424_highconfidence_interactome.csv"

def out(fname):
    return WORK_DIR / fname

# ── validate inputs ───────────────────────────────────────────────────────────
print("="*70)
print("Cortex Correlation v2 — Input file check")
print("="*70)
inputs = {
    "Excel proteomics (PXD062542)": EXCEL_FILE,
    "mRNA DEG (GSE124439_clean)":   MRNA_FILE,
    "Spinal cord results (read-only)": SC_RESULTS,
    "PXD065424 interactome (read-only)": INTERACTOME,
}
missing = []
for label, path in inputs.items():
    status = "OK" if path.exists() else "MISSING"
    print(f"  [{status}]  {label}: {path.name}")
    if not path.exists():
        missing.append(str(path))

if missing:
    print(f"\nERROR: {len(missing)} required file(s) not found:")
    for m in missing:
        print(f"  {m}")
    sys.exit(1)
print("\nAll inputs confirmed. Proceeding.\n")

# ── gene symbol validation ────────────────────────────────────────────────────
VALID_SYMBOL = re.compile(r'^[A-Za-z][A-Za-z0-9\-]*$')

def is_valid_symbol(s):
    if pd.isna(s) or str(s).strip() == "":
        return False
    return bool(VALID_SYMBOL.match(str(s).strip()))


# =============================================================================
# STEP 1  Extract and clean motor cortex protein data from PXD062542
# =============================================================================
print("="*70)
print("STEP 1: Extract motor cortex protein data from PXD062542")
print("="*70)

# Load the "Human limma results" sheet — ignore "Mouse limma results"
print('  Loading sheet "Human limma results" from Excel...')
prot_raw = pd.read_excel(EXCEL_FILE, sheet_name="Human limma results")
print(f"  Raw shape: {prot_raw.shape[0]} rows x {prot_raw.shape[1]} columns")

# Confirm required motor cortex columns are present
required_cols = {
    "Gene.Symbol":              "Gene_Symbol",
    "M.Mcx_vs_C.Mcx_diff":     "protein_logFC",
    "M.Mcx_vs_C.Mcx_p.adj":    "protein_adjPval",
    "M.Mcx_vs_C.Mcx_significant": "protein_significant",
}
for col in required_cols:
    if col not in prot_raw.columns:
        print(f"\nERROR: Expected column '{col}' not found in sheet.")
        print(f"  Available columns: {prot_raw.columns.tolist()}")
        sys.exit(1)

# Extract and rename the four columns we need
prot = prot_raw[list(required_cols.keys())].copy()
prot = prot.rename(columns=required_cols)

n_raw = len(prot)

# Strip whitespace from gene symbols
prot["Gene_Symbol"] = prot["Gene_Symbol"].astype(str).str.strip()

# Drop rows where logFC or adj.P.Val is missing
prot = prot.dropna(subset=["protein_logFC", "protein_adjPval"])
n_after_na = len(prot)

# Filter to valid HGNC gene symbols
valid_mask = prot["Gene_Symbol"].apply(is_valid_symbol)
n_invalid = (~valid_mask).sum()
prot = prot[valid_mask].copy()

# Deduplicate: keep row with lowest adj.P.Val per gene
prot = (prot
        .sort_values("protein_adjPval")
        .drop_duplicates(subset="Gene_Symbol", keep="first")
        .reset_index(drop=True))

print(f"  {n_raw} raw rows")
print(f"  -{n_raw - n_after_na} dropped (missing logFC or adj.P.Val)")
print(f"  -{n_invalid} dropped (invalid HGNC symbol)")
print(f"  {len(prot)} proteins remain after cleaning")

prot.to_csv(out("PXD062542_CTX_clean.csv"), index=False)
print("  Saved: PXD062542_CTX_clean.csv\n")


# =============================================================================
# STEP 2  Cortex mRNA-protein correlation
# =============================================================================
print("="*70)
print("STEP 2: Cortex mRNA-protein correlation (GSE124439 x PXD062542)")
print("="*70)

mrna = pd.read_csv(MRNA_FILE)

# Inner join on Gene_Symbol
merged = pd.merge(mrna, prot, on="Gene_Symbol", how="inner")
print(f"  mRNA genes: {len(mrna)}")
print(f"  Protein genes: {len(prot)}")
print(f"  Overlap (inner join): {len(merged)} genes")

# Drop any rows where either logFC is NaN (defensive)
merged = merged.dropna(subset=["mRNA_logFC", "protein_logFC"])
n = len(merged)

# Discordance score = protein_logFC - mRNA_logFC
merged["discordance_score"] = merged["protein_logFC"] - merged["mRNA_logFC"]

# Quadrant classification (threshold = 0.5)
THRESH = 0.5

def classify_quadrant(mrna_lfc, prot_lfc, thresh=THRESH):
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

merged["quadrant"] = merged.apply(
    lambda r: classify_quadrant(r["mRNA_logFC"], r["protein_logFC"]), axis=1
)

# Pearson and Spearman correlations
pearson_r, pearson_p = pearsonr(merged["mRNA_logFC"], merged["protein_logFC"])
spearman_r, spearman_p = spearmanr(merged["mRNA_logFC"], merged["protein_logFC"])

print(f"\n  Pearson r  = {pearson_r:.4f}  (p = {pearson_p:.3e})")
print(f"  Spearman r = {spearman_r:.4f}  (p = {spearman_p:.3e})")
print("\n  Quadrant counts:")
for q, cnt in merged["quadrant"].value_counts().items():
    marker = "  <-- UPS candidates" if q == "discordant_LOW" else ""
    print(f"    {q}: {cnt}{marker}")

# ── scatter plot ──────────────────────────────────────────────────────────────
QUADRANT_COLORS = {
    "concordant_up":   "#2196F3",
    "concordant_down": "#9C27B0",
    "discordant_HIGH": "#FF9800",
    "discordant_LOW":  "#F44336",
    "ambiguous":       "#BDBDBD",
}

fig, ax = plt.subplots(figsize=(9, 8))

for quad, color in QUADRANT_COLORS.items():
    mask = merged["quadrant"] == quad
    ax.scatter(
        merged.loc[mask, "mRNA_logFC"],
        merged.loc[mask, "protein_logFC"],
        c=color, label=f"{quad} (n={mask.sum()})",
        alpha=0.65, s=25, linewidths=0,
    )

# Diagonal reference line y = x
lim = max(
    merged["mRNA_logFC"].abs().max(),
    merged["protein_logFC"].abs().max()
) * 1.1
ax.plot([-lim, lim], [-lim, lim], "k--", linewidth=0.8, alpha=0.5, label="y = x")

# Label 10 most discordant genes by absolute discordance_score
top_disc = merged.reindex(
    merged["discordance_score"].abs().nlargest(10).index
)
for _, row in top_disc.iterrows():
    ax.annotate(
        row["Gene_Symbol"],
        xy=(row["mRNA_logFC"], row["protein_logFC"]),
        xytext=(5, 5), textcoords="offset points",
        fontsize=7, alpha=0.85,
    )

ax.set_xlabel("mRNA log2FC (GSE124439, motor cortex ALS vs control)", fontsize=11)
ax.set_ylabel("Protein log2FC (PXD062542, motor cortex ALS vs control)", fontsize=11)
ax.set_title(
    f"Motor Cortex mRNA-Protein Correlation (PXD062542 v2)\n"
    f"Pearson r = {pearson_r:.3f}  |  Spearman r = {spearman_r:.3f}  |  n = {n}",
    fontsize=12,
)
ax.legend(fontsize=8, loc="upper left")
ax.axhline(0, color="grey", linewidth=0.5, linestyle=":")
ax.axvline(0, color="grey", linewidth=0.5, linestyle=":")
ax.set_xlim(-lim, lim)
ax.set_ylim(-lim, lim)
plt.tight_layout()
plt.savefig(out("cortex_correlation_plot_v2.png"), dpi=150)
plt.close()
print("\n  Saved plot: cortex_correlation_plot_v2.png")

# ── results CSV ───────────────────────────────────────────────────────────────
ctx_results = merged[[
    "Gene_Symbol", "mRNA_logFC", "mRNA_adjPval",
    "protein_logFC", "protein_adjPval", "protein_significant",
    "discordance_score", "quadrant"
]].copy()
ctx_results.to_csv(out("cortex_correlation_results_v2.csv"), index=False)
print("  Saved results: cortex_correlation_results_v2.csv\n")


# =============================================================================
# STEP 3  Update cross-tissue STRING inputs
# =============================================================================
print("="*70)
print("STEP 3: Updating cross-tissue STRING input lists")
print("="*70)

# Load spinal cord results (read-only — do not modify)
sc = pd.read_csv(SC_RESULTS)

# Extract discordant_LOW genes per tissue
sc_low = sc[sc["quadrant"] == "discordant_LOW"][[
    "Gene_Symbol", "mRNA_logFC_mean", "protein_logFC", "discordance_score"
]].rename(columns={
    "mRNA_logFC_mean":  "spinal_mRNA_logFC",
    "protein_logFC":    "spinal_protein_logFC",
    "discordance_score": "spinal_discordance",
})

ctx_low = ctx_results[ctx_results["quadrant"] == "discordant_LOW"][[
    "Gene_Symbol", "mRNA_logFC", "protein_logFC", "discordance_score"
]].rename(columns={
    "mRNA_logFC":       "cortex_mRNA_logFC",
    "protein_logFC":    "cortex_protein_logFC",
    "discordance_score": "cortex_discordance",
})

print(f"  discordant_LOW spinal cord: {len(sc_low)} genes")
print(f"  discordant_LOW cortex (v2): {len(ctx_low)} genes")

# Overlap between the two discordant_LOW sets
sc_low_genes  = set(sc_low["Gene_Symbol"])
ctx_low_genes = set(ctx_low["Gene_Symbol"])
both_tissues  = sc_low_genes & ctx_low_genes
print(f"\n  Overlap (discordant_LOW in both tissues): {len(both_tissues)} genes")
if both_tissues:
    print(f"    {sorted(both_tissues)}")

# Outer join to combine all discordant_LOW genes across tissues
all_low = pd.merge(sc_low, ctx_low, on="Gene_Symbol", how="outer")
all_low["present_in_spinal"] = all_low["spinal_mRNA_logFC"].notna()
all_low["present_in_cortex"]  = all_low["cortex_mRNA_logFC"].notna()

# High-confidence: discordant_LOW in BOTH tissues
high_conf = all_low[all_low["present_in_spinal"] & all_low["present_in_cortex"]].copy()

# Extended: discordant_LOW in EITHER tissue
extended = all_low.copy()

print(f"\n  High-confidence STRING input (both tissues): {len(high_conf)} genes")
print(f"  Extended STRING input (either tissue):       {len(extended)} genes")

high_conf.to_csv(out("ALS_Ubiquitome_STRING_input_highconfidence.csv"), index=False)
extended.to_csv(out("ALS_Ubiquitome_STRING_input_extended.csv"), index=False)
print("  Saved (overwritten): ALS_Ubiquitome_STRING_input_highconfidence.csv")
print("  Saved (overwritten): ALS_Ubiquitome_STRING_input_extended.csv")

# ── overlap with PXD065424 interactome ───────────────────────────────────────
print("\n  Checking overlap with PXD065424 high-confidence interactome...")
interactome = pd.read_csv(INTERACTOME)
print(f"  Interactome genes: {len(interactome)}")

# Overlap: interactome genes that appear in the extended (either-tissue) discordant_LOW list
priority_df = pd.merge(
    interactome,
    extended[[
        "Gene_Symbol", "spinal_mRNA_logFC", "spinal_protein_logFC", "spinal_discordance",
        "cortex_mRNA_logFC", "cortex_protein_logFC", "cortex_discordance",
        "present_in_spinal", "present_in_cortex",
    ]],
    on="Gene_Symbol", how="inner"
)
print(f"  Priority nodes (interactome ∩ discordant_LOW): {len(priority_df)} genes")
if len(priority_df) > 0:
    print(f"    {priority_df['Gene_Symbol'].tolist()}")

priority_df.to_csv(out("ALS_Ubiquitome_priority_nodes.csv"), index=False)
print("  Saved (overwritten): ALS_Ubiquitome_priority_nodes.csv\n")


# =============================================================================
# STEP 4  Summary
# =============================================================================
print("="*70)
print("STEP 4: Summary")
print("="*70)

# Cross-tissue discordant_LOW comparison
ctx_low_vs_sc = ctx_low_genes & sc_low_genes
sc_only = sc_low_genes - ctx_low_genes
ctx_only = ctx_low_genes - sc_low_genes

summary = f"""
Cortex Correlation v2 — Summary Report
{"="*70}

Source datasets
  mRNA:    GSE124439 (motor cortex, ALS vs control, RNA-seq)
  Protein: PXD062542 (motor cortex, ALS vs control, TMT proteomics)
           Pan et al., Scientific Reports 2025
           Sheet: "Human limma results", columns: M.Mcx_vs_C.Mcx_*

STEP 1 — PXD062542 motor cortex protein cleaning
  Proteins after cleaning:  {len(prot)}

STEP 2 — Cortex mRNA-protein correlation
  Genes in overlap (inner join):  {n}
  Pearson r:   {pearson_r:.4f}  (p = {pearson_p:.3e})
  Spearman r:  {spearman_r:.4f}  (p = {spearman_p:.3e})

  Quadrant counts:
"""
for q, cnt in merged["quadrant"].value_counts().items():
    summary += f"    {q}: {cnt}\n"

summary += f"""
  discordant_LOW genes (cortex v2):
{chr(10).join("    " + g for g in sorted(ctx_low_genes)) if ctx_low_genes else "    (none)"}

STEP 3 — Updated cross-tissue STRING inputs
  High-confidence (discordant_LOW in BOTH tissues):  {len(high_conf)} genes
{chr(10).join("    " + g for g in sorted(both_tissues)) if both_tissues else "    (none)"}

  Extended (discordant_LOW in EITHER tissue):         {len(extended)} genes

  PXD065424 interactome x discordant_LOW (priority nodes): {len(priority_df)} genes
{chr(10).join("    " + g for g in priority_df['Gene_Symbol'].tolist()) if len(priority_df) > 0 else "    (none)"}

STEP 4 — Cross-tissue discordant_LOW comparison
  discordant_LOW spinal cord (v1, unchanged): {len(sc_low)} genes
    {sorted(sc_low_genes)}
  discordant_LOW cortex (v2):                 {len(ctx_low)} genes
    {sorted(ctx_low_genes) if ctx_low_genes else "(none)"}

  Overlap (in both tissues):     {len(ctx_low_vs_sc)} genes  -> {sorted(ctx_low_vs_sc) if ctx_low_vs_sc else "(none)"}
  Spinal cord only:              {len(sc_only)} genes
  Cortex only:                   {len(ctx_only)} genes

{"="*70}
Output files (scripts/data_processing/mRNA_protein_corr/):
  PXD062542_CTX_clean.csv
  cortex_correlation_results_v2.csv
  cortex_correlation_plot_v2.png
  ALS_Ubiquitome_STRING_input_highconfidence.csv  [overwritten]
  ALS_Ubiquitome_STRING_input_extended.csv        [overwritten]
  ALS_Ubiquitome_priority_nodes.csv               [overwritten]
  cortex_v2_summary.txt
"""

print(summary)
with open(out("cortex_v2_summary.txt"), "w", encoding="utf-8") as fh:
    fh.write(summary)
print("Saved: cortex_v2_summary.txt")
print("\nDone.")
