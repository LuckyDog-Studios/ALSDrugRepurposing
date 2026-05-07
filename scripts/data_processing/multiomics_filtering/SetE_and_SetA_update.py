# pip install pandas numpy openpyxl
import os
import re
import sys
import numpy as np
import pandas as pd

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
BASE_DIR = os.path.join(os.path.dirname(__file__), "..", "..", "..")

INPUT_LCM = {
    "increased_FDR05": os.path.join(BASE_DIR, "differential_expression", "gene expression", "LCM_meta_SC_increased_FDR05.XLSX"),
    "decreased_FDR05": os.path.join(BASE_DIR, "differential_expression", "gene expression", "LCM_meta_SC_decreased_FDR05.XLSX"),
    "increased_FDR10": os.path.join(BASE_DIR, "differential_expression", "gene expression", "LCM_meta_SC_increased_FDR10.XLSX"),
    "decreased_FDR10": os.path.join(BASE_DIR, "differential_expression", "gene expression", "LCM_meta_SC_decreased_FDR10.XLSX"),
}
INPUT_SC_MRNA    = os.path.join(BASE_DIR, "mRNA_protein_corr_results", "spinal_cord_mRNA_merged.csv")
INPUT_SC_PROTEIN = os.path.join(BASE_DIR, "multiomics_filtering_results", "SC_protein_filtered.csv")
INPUT_CTX_MRNA   = os.path.join(BASE_DIR, "differential_expression", "gene expression", "GSE124439_clean.csv")
INPUT_CTX_PROTEIN= os.path.join(BASE_DIR, "differential_expression", "protein expression", "PXD062542_MOESM1_ESM.xlsx")
INPUT_SETA_V1    = os.path.join(BASE_DIR, "multiomics_filtering_results", "SetA_combined.csv")

OUT_DIR = os.path.join(BASE_DIR, "multiomics_filtering_results")
os.makedirs(OUT_DIR, exist_ok=True)

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
HGNC_RE = re.compile(r"^[A-Za-z][A-Za-z0-9\-]*$")

def is_hgnc_column(series: pd.Series) -> bool:
    """Return True if ≥80 % of non-null values look like HGNC symbols."""
    valid = series.dropna().astype(str)
    if len(valid) == 0:
        return False
    return (valid.str.match(r"^[A-Za-z][A-Za-z0-9\-]*$").sum() / len(valid)) >= 0.80

def find_column(df: pd.DataFrame, candidates: list[str], label: str) -> str:
    """Case-insensitive column finder; raises with helpful message on failure."""
    lower_map = {c.lower(): c for c in df.columns}
    for cand in candidates:
        if cand.lower() in lower_map:
            return lower_map[cand.lower()]
    raise ValueError(f"Cannot find {label} column. Available columns: {list(df.columns)}")

def auto_gene_col(df: pd.DataFrame) -> str:
    """Find the gene-symbol column automatically."""
    # Try known names first
    known = ["Gene.Symbol", "Gene_Symbol", "Symbol", "gene", "GeneName", "gene_name", "SYMBOL"]
    lower_map = {c.lower(): c for c in df.columns}
    for k in known:
        if k.lower() in lower_map:
            col = lower_map[k.lower()]
            if is_hgnc_column(df[col]):
                return col
    # Fallback: scan all string columns
    for col in df.columns:
        if df[col].dtype == object and is_hgnc_column(df[col]):
            return col
    raise ValueError(f"Cannot identify gene symbol column. Columns are: {list(df.columns)}\n"
                     "Please specify which column contains HGNC gene symbols.")

def clean_gene_col(series: pd.Series) -> pd.Series:
    s = series.astype(str).str.strip()
    s = s.replace({"nan": np.nan, "": np.nan, "NA": np.nan, "None": np.nan})
    valid_mask = s.str.match(r"^[A-Za-z][A-Za-z0-9\-]*$", na=False)
    s[~valid_mask] = np.nan
    return s

# ---------------------------------------------------------------------------
# Pre-flight: confirm all inputs exist
# ---------------------------------------------------------------------------
print("\n=== Pre-flight: checking input files ===")
all_inputs = {
    **{k: v for k, v in INPUT_LCM.items()},
    "SC mRNA":       INPUT_SC_MRNA,
    "SC protein":    INPUT_SC_PROTEIN,
    "CTX mRNA":      INPUT_CTX_MRNA,
    "CTX protein":   INPUT_CTX_PROTEIN,
}
missing = []
for label, path in all_inputs.items():
    exists = os.path.isfile(path)
    status = "OK" if exists else "MISSING"
    print(f"  [{status}] {label}: {path}")
    if not exists:
        missing.append(path)

if missing:
    print("\nERROR: Missing input files — aborting.")
    for p in missing:
        print(f"  {p}")
    sys.exit(1)
print("All inputs found.\n")

# ---------------------------------------------------------------------------
# PART 1 — Extract Set E (LCM meta-analysis gene lists)
# ---------------------------------------------------------------------------
print("=" * 60)
print("PART 1: Extracting Set E (LCM meta-analysis DEGs)")
print("=" * 60)

# ---------------------------------------------------------------------------
# Step 1: Inspect all 4 LCM Excel files
# ---------------------------------------------------------------------------
print("\n--- Step 1: Inspecting LCM Excel files ---")
inspection_lines = []

raw_lcm = {}
for key, path in INPUT_LCM.items():
    print(f"  Loading {key} …")
    xl = pd.ExcelFile(path)
    # Sheet may be named "A" or "(A) DEG List" etc.
    sheet = None
    for s in xl.sheet_names:
        if s == "A" or s.startswith("(A)"):
            sheet = s
            break
    if sheet is None:
        print(f"  WARNING: no sheet 'A' or '(A)*' in {key}. Available: {xl.sheet_names}")
        sheet = xl.sheet_names[0]
    df = xl.parse(sheet)
    raw_lcm[key] = df

    header = f"\n{'='*50}\nFile: {key}\nPath: {path}\nSheet: {sheet}\n"
    cols   = f"Columns: {list(df.columns)}\n"
    nrows  = f"N rows: {len(df)}\n"
    sample = f"First 3 rows:\n{df.head(3).to_string()}\n"
    block  = header + cols + nrows + sample
    inspection_lines.append(block)
    print(f"    {len(df)} rows, columns: {list(df.columns)}")

inspection_path = os.path.join(OUT_DIR, "SetE_inspection.txt")
with open(inspection_path, "w") as f:
    f.write("LCM Excel File Inspection\n")
    f.writelines(inspection_lines)
print(f"\nInspection saved -> {inspection_path}")

# ---------------------------------------------------------------------------
# Step 2: Extract and clean each LCM file
# ---------------------------------------------------------------------------
print("\n--- Step 2: Cleaning LCM gene lists ---")

SMD_CANDIDATES = ["SMD", "smd", "Hedges_g", "effect_size", "EffectSize"]
FDR_CANDIDATES = ["FDR", "fdr", "adj.P.Val", "adjPval", "padj", "p.adj", "BH", "adjusted_pvalue"]

def clean_lcm_df(df: pd.DataFrame, key: str) -> pd.DataFrame:
    gene_col = auto_gene_col(df)
    smd_col  = find_column(df, SMD_CANDIDATES, "SMD")
    fdr_col  = find_column(df, FDR_CANDIDATES, "FDR")

    print(f"    {key}: gene='{gene_col}', SMD='{smd_col}', FDR='{fdr_col}'")

    out = df[[gene_col, smd_col, fdr_col]].copy()
    out.columns = ["Gene_Symbol", "LCM_SMD", "LCM_FDR"]
    out["Gene_Symbol"] = clean_gene_col(out["Gene_Symbol"])
    out = out.dropna(subset=["Gene_Symbol"])
    out["LCM_SMD"] = pd.to_numeric(out["LCM_SMD"], errors="coerce")
    out["LCM_FDR"] = pd.to_numeric(out["LCM_FDR"], errors="coerce")
    out = out.dropna(subset=["LCM_SMD"])

    # De-duplicate: keep highest |SMD|
    out["abs_SMD"] = out["LCM_SMD"].abs()
    out = out.sort_values("abs_SMD", ascending=False).drop_duplicates(subset="Gene_Symbol").drop(columns="abs_SMD")

    out["direction"] = np.where(out["LCM_SMD"] > 0, "UP", "DOWN")
    out = out[["Gene_Symbol", "direction", "LCM_SMD", "LCM_FDR"]]
    return out.reset_index(drop=True)

cleaned_lcm = {k: clean_lcm_df(df, k) for k, df in raw_lcm.items()}
for k, df in cleaned_lcm.items():
    print(f"    {k}: {len(df)} genes after cleaning")

# Primary Set E (FDR05)
sete_fdr05 = pd.concat(
    [cleaned_lcm["increased_FDR05"], cleaned_lcm["decreased_FDR05"]],
    ignore_index=True
)
# De-duplicate across files (keep highest |SMD|)
sete_fdr05["abs_SMD"] = sete_fdr05["LCM_SMD"].abs()
sete_fdr05 = sete_fdr05.sort_values("abs_SMD", ascending=False)\
                        .drop_duplicates(subset="Gene_Symbol")\
                        .drop(columns="abs_SMD")\
                        .reset_index(drop=True)

sete_fdr05_path = os.path.join(OUT_DIR, "SetE_LCM_FDR05_combined.csv")
sete_fdr05.to_csv(sete_fdr05_path, index=False)
n_up_05   = (sete_fdr05["direction"] == "UP").sum()
n_down_05 = (sete_fdr05["direction"] == "DOWN").sum()
print(f"\nSet E FDR05: {len(sete_fdr05)} total | UP={n_up_05} | DOWN={n_down_05}")
print(f"  Saved -> {sete_fdr05_path}")

# Secondary Set E (FDR10)
sete_fdr10 = pd.concat(
    [cleaned_lcm["increased_FDR10"], cleaned_lcm["decreased_FDR10"]],
    ignore_index=True
)
sete_fdr10["abs_SMD"] = sete_fdr10["LCM_SMD"].abs()
sete_fdr10 = sete_fdr10.sort_values("abs_SMD", ascending=False)\
                        .drop_duplicates(subset="Gene_Symbol")\
                        .drop(columns="abs_SMD")\
                        .reset_index(drop=True)

sete_fdr10_path = os.path.join(OUT_DIR, "SetE_LCM_FDR10_combined.csv")
sete_fdr10.to_csv(sete_fdr10_path, index=False)

fdr05_genes = set(sete_fdr05["Gene_Symbol"])
sete_fdr10_additional = sete_fdr10[~sete_fdr10["Gene_Symbol"].isin(fdr05_genes)].reset_index(drop=True)
sete_fdr10_additional_path = os.path.join(OUT_DIR, "SetE_LCM_FDR10_additional.csv")
sete_fdr10_additional.to_csv(sete_fdr10_additional_path, index=False)

print(f"Set E FDR10: {len(sete_fdr10)} total | {len(sete_fdr10_additional)} additional beyond FDR05")
print(f"  Saved -> {sete_fdr10_path}")
print(f"  Additional -> {sete_fdr10_additional_path}")

# ---------------------------------------------------------------------------
# PART 2 — Rebuild Set A
# ---------------------------------------------------------------------------
print("\n" + "=" * 60)
print("PART 2: Rebuilding Set A with expanded SC mRNA")
print("=" * 60)

# ---------------------------------------------------------------------------
# Step 3: Build expanded SC mRNA list
# ---------------------------------------------------------------------------
print("\n--- Step 3: Building expanded SC mRNA list ---")

geo = pd.read_csv(INPUT_SC_MRNA)
# The file is already filtered; add required columns
geo["mRNA_source"] = "GEO_merged"
geo["direction"] = np.where(geo["mRNA_logFC_mean"] > 0, "UP", "DOWN")
# Standardise to working columns
geo_work = geo[["Gene_Symbol", "direction", "mRNA_source", "mRNA_logFC_mean"]].copy()
geo_work["LCM_SMD"] = np.nan
print(f"  GEO SC mRNA: {len(geo_work)} genes")

lcm_work = sete_fdr05[["Gene_Symbol", "direction", "LCM_SMD"]].copy()
lcm_work["mRNA_source"] = "LCM_meta"
lcm_work["mRNA_logFC_mean"] = np.nan
print(f"  LCM meta SC mRNA: {len(lcm_work)} genes")

# Merge
geo_genes = set(geo_work["Gene_Symbol"])
lcm_genes = set(lcm_work["Gene_Symbol"])
both_genes = geo_genes & lcm_genes
geo_only   = geo_genes - lcm_genes
lcm_only   = lcm_genes - geo_genes

rows = []
# GEO-only
rows.append(geo_work[geo_work["Gene_Symbol"].isin(geo_only)].copy())

# LCM-only
rows.append(lcm_work[lcm_work["Gene_Symbol"].isin(lcm_only)].copy())

# Both sources
concordant_count  = 0
discordant_count  = 0
discordant_genes  = set()
for gene in both_genes:
    g_dir = geo_work.loc[geo_work["Gene_Symbol"] == gene, "direction"].iloc[0]
    l_dir = lcm_work.loc[lcm_work["Gene_Symbol"] == gene, "direction"].iloc[0]
    geo_row = geo_work[geo_work["Gene_Symbol"] == gene].iloc[0]
    lcm_row = lcm_work[lcm_work["Gene_Symbol"] == gene].iloc[0]

    if g_dir == l_dir:
        concordant_count += 1
        row = {
            "Gene_Symbol":    gene,
            "direction":      g_dir,
            "mRNA_source":    "GEO_and_LCM",
            "mRNA_logFC_mean": geo_row["mRNA_logFC_mean"],
            "LCM_SMD":        lcm_row["LCM_SMD"],
        }
    else:
        discordant_count += 1
        discordant_genes.add(gene)
        row = {
            "Gene_Symbol":    gene,
            "direction":      "DISCORDANT",
            "mRNA_source":    "GEO_and_LCM",
            "mRNA_logFC_mean": geo_row["mRNA_logFC_mean"],
            "LCM_SMD":        lcm_row["LCM_SMD"],
        }
    rows.append(pd.DataFrame([row]))

sc_mrna_expanded = pd.concat(rows, ignore_index=True)
sc_mrna_expanded = sc_mrna_expanded[["Gene_Symbol", "direction", "mRNA_source",
                                      "mRNA_logFC_mean", "LCM_SMD"]]

sc_mrna_exp_path = os.path.join(OUT_DIR, "SC_mRNA_expanded.csv")
sc_mrna_expanded.to_csv(sc_mrna_exp_path, index=False)

print(f"  GEO only:              {len(geo_only)}")
print(f"  LCM meta only:         {len(lcm_only)}")
print(f"  Both — concordant:     {concordant_count}")
print(f"  Both — discordant:     {discordant_count} (excluded from Set A)")
print(f"  Total in expanded list:{len(sc_mrna_expanded)}")
print(f"  Saved -> {sc_mrna_exp_path}")

# ---------------------------------------------------------------------------
# Step 4: Rebuild Set A — spinal cord
# ---------------------------------------------------------------------------
print("\n--- Step 4: Rebuilding Set A — spinal cord ---")

sc_mrna_filt = sc_mrna_expanded[sc_mrna_expanded["direction"] != "DISCORDANT"].copy()
sc_prot = pd.read_csv(INPUT_SC_PROTEIN)

sc_joined = pd.merge(sc_mrna_filt, sc_prot, on="Gene_Symbol", how="inner")
sc_joined = sc_joined[sc_joined["direction_x"] == sc_joined["direction_y"]].copy()

seta_sc = pd.DataFrame({
    "Gene_Symbol":      sc_joined["Gene_Symbol"],
    "SC_mRNA_logFC":    sc_joined["mRNA_logFC_mean"],
    "SC_LCM_SMD":       sc_joined["LCM_SMD"],
    "SC_protein_logFC": sc_joined["logFC"],
    "SC_protein_adjPval": sc_joined["adjPval"],
    "SC_mRNA_source":   sc_joined["mRNA_source"],
    "direction":        sc_joined["direction_x"],
    "evidence_count":   2,
    "tissue":           "spinal_cord",
})

seta_sc_path = os.path.join(OUT_DIR, "SetA_spinal_cord_v2.csv")
seta_sc.to_csv(seta_sc_path, index=False)

src_counts = seta_sc["SC_mRNA_source"].value_counts()
print(f"  SetA spinal cord v2: {len(seta_sc)} genes")
print(f"    UP:            {(seta_sc['direction']=='UP').sum()}")
print(f"    DOWN:          {(seta_sc['direction']=='DOWN').sum()}")
print(f"    GEO_merged:    {src_counts.get('GEO_merged', 0)}")
print(f"    LCM_meta:      {src_counts.get('LCM_meta', 0)}")
print(f"    GEO_and_LCM:   {src_counts.get('GEO_and_LCM', 0)}")
print(f"  Saved -> {seta_sc_path}")

# ---------------------------------------------------------------------------
# Step 5: Rebuild Set A — motor cortex
# ---------------------------------------------------------------------------
print("\n--- Step 5: Rebuilding Set A — motor cortex ---")

ctx_mrna_raw = pd.read_csv(INPUT_CTX_MRNA)
ctx_mrna = ctx_mrna_raw.copy()
# Find required columns
mrna_fc_col  = find_column(ctx_mrna, ["mRNA_logFC", "logFC", "log2FoldChange", "log2FC"], "mRNA_logFC")
mrna_pv_col  = find_column(ctx_mrna, ["mRNA_adjPval", "adjPval", "padj", "adj.P.Val", "FDR"], "mRNA_adjPval")
gene_col_ctx = find_column(ctx_mrna, ["Gene_Symbol", "Gene.Symbol", "gene", "Symbol"], "Gene_Symbol")

ctx_mrna = ctx_mrna.rename(columns={gene_col_ctx: "Gene_Symbol",
                                     mrna_fc_col:  "mRNA_logFC",
                                     mrna_pv_col:  "mRNA_adjPval"})
ctx_mrna["Gene_Symbol"] = clean_gene_col(ctx_mrna["Gene_Symbol"])
ctx_mrna = ctx_mrna.dropna(subset=["Gene_Symbol", "mRNA_logFC", "mRNA_adjPval"])
ctx_mrna["mRNA_logFC"]   = pd.to_numeric(ctx_mrna["mRNA_logFC"],   errors="coerce")
ctx_mrna["mRNA_adjPval"] = pd.to_numeric(ctx_mrna["mRNA_adjPval"], errors="coerce")
ctx_mrna = ctx_mrna[ctx_mrna["mRNA_logFC"].abs() > 0.5]
ctx_mrna = ctx_mrna[ctx_mrna["mRNA_adjPval"] < 0.05]
ctx_mrna["direction"] = np.where(ctx_mrna["mRNA_logFC"] > 0, "UP", "DOWN")
ctx_mrna = ctx_mrna[["Gene_Symbol", "mRNA_logFC", "mRNA_adjPval", "direction"]]\
           .drop_duplicates(subset="Gene_Symbol").reset_index(drop=True)
print(f"  CTX mRNA passing filter: {len(ctx_mrna)} genes")

ctx_prot_xl = pd.ExcelFile(INPUT_CTX_PROTEIN)
print(f"  CTX protein sheets: {ctx_prot_xl.sheet_names}")
ctx_prot_raw = ctx_prot_xl.parse("Human limma results")
ctx_prot = ctx_prot_raw[["Gene.Symbol", "M.Mcx_vs_C.Mcx_diff", "M.Mcx_vs_C.Mcx_p.adj"]].copy()
ctx_prot.columns = ["Gene_Symbol", "logFC", "adjPval"]
ctx_prot["Gene_Symbol"] = clean_gene_col(ctx_prot["Gene_Symbol"])
ctx_prot = ctx_prot.dropna(subset=["Gene_Symbol"])
ctx_prot["logFC"]   = pd.to_numeric(ctx_prot["logFC"],   errors="coerce")
ctx_prot["adjPval"] = pd.to_numeric(ctx_prot["adjPval"], errors="coerce")
ctx_prot = ctx_prot.dropna(subset=["logFC", "adjPval"])
# De-duplicate: keep lowest adjPval
ctx_prot = ctx_prot.sort_values("adjPval").drop_duplicates(subset="Gene_Symbol").reset_index(drop=True)
ctx_prot = ctx_prot[ctx_prot["logFC"].abs() > 0.3]
ctx_prot = ctx_prot[ctx_prot["adjPval"] < 0.05]
ctx_prot["direction"] = np.where(ctx_prot["logFC"] > 0, "UP", "DOWN")
print(f"  CTX protein passing filter: {len(ctx_prot)} genes")

ctx_joined = pd.merge(ctx_mrna, ctx_prot, on="Gene_Symbol", how="inner")
ctx_joined = ctx_joined[ctx_joined["direction_x"] == ctx_joined["direction_y"]].copy()

seta_ctx = pd.DataFrame({
    "Gene_Symbol":        ctx_joined["Gene_Symbol"],
    "CTX_mRNA_logFC":     ctx_joined["mRNA_logFC"],
    "CTX_mRNA_adjPval":   ctx_joined["mRNA_adjPval"],
    "CTX_protein_logFC":  ctx_joined["logFC"],
    "CTX_protein_adjPval":ctx_joined["adjPval"],
    "direction":          ctx_joined["direction_x"],
    "evidence_count":     2,
    "tissue":             "cortex",
    "SC_mRNA_source":     np.nan,
})

seta_ctx_path = os.path.join(OUT_DIR, "SetA_cortex_v2.csv")
seta_ctx.to_csv(seta_ctx_path, index=False)
print(f"  SetA cortex v2: {len(seta_ctx)} genes")
print(f"    UP:   {(seta_ctx['direction']=='UP').sum()}")
print(f"    DOWN: {(seta_ctx['direction']=='DOWN').sum()}")
print(f"  Saved -> {seta_ctx_path}")

# ---------------------------------------------------------------------------
# Step 6: Combine into SetA_combined_v2
# ---------------------------------------------------------------------------
print("\n--- Step 6: Building SetA_combined_v2 ---")

sc_genes  = set(seta_sc["Gene_Symbol"])
ctx_genes = set(seta_ctx["Gene_Symbol"])

all_genes = sc_genes | ctx_genes
records = []

sc_idx  = seta_sc.set_index("Gene_Symbol")
ctx_idx = seta_ctx.set_index("Gene_Symbol")

for gene in sorted(all_genes):
    in_sc  = gene in sc_genes
    in_ctx = gene in ctx_genes

    sc_row  = sc_idx.loc[gene]  if in_sc  else None
    ctx_row = ctx_idx.loc[gene] if in_ctx else None

    dir_sc  = sc_row["direction"]  if in_sc  else np.nan
    dir_ctx = ctx_row["direction"] if in_ctx else np.nan

    cross       = False
    cross_disc  = False
    ev_count    = 0
    if in_sc:
        ev_count += 2
    if in_ctx:
        ev_count += 2
    if in_sc and in_ctx:
        if dir_sc == dir_ctx:
            cross = True
        else:
            cross_disc = True

    rec = {
        "Gene_Symbol":          gene,
        "present_in_SC":        in_sc,
        "present_in_CTX":       in_ctx,
        "direction_SC":         dir_sc,
        "direction_CTX":        dir_ctx,
        "SC_mRNA_logFC":        sc_row["SC_mRNA_logFC"]      if in_sc  else np.nan,
        "SC_LCM_SMD":           sc_row["SC_LCM_SMD"]         if in_sc  else np.nan,
        "SC_mRNA_source":       sc_row["SC_mRNA_source"]     if in_sc  else np.nan,
        "SC_protein_logFC":     sc_row["SC_protein_logFC"]   if in_sc  else np.nan,
        "SC_protein_adjPval":   sc_row["SC_protein_adjPval"] if in_sc  else np.nan,
        "CTX_mRNA_logFC":       ctx_row["CTX_mRNA_logFC"]      if in_ctx else np.nan,
        "CTX_mRNA_adjPval":     ctx_row["CTX_mRNA_adjPval"]    if in_ctx else np.nan,
        "CTX_protein_logFC":    ctx_row["CTX_protein_logFC"]   if in_ctx else np.nan,
        "CTX_protein_adjPval":  ctx_row["CTX_protein_adjPval"] if in_ctx else np.nan,
        "evidence_count":       ev_count,
        "cross_tissue":         cross,
        "cross_tissue_discordant": cross_disc,
    }
    records.append(rec)

seta_v2 = pd.DataFrame(records)
seta_v2_path = os.path.join(OUT_DIR, "SetA_combined_v2.csv")
seta_v2.to_csv(seta_v2_path, index=False)

n_sc_only    = seta_v2[seta_v2["present_in_SC"] & ~seta_v2["present_in_CTX"]].shape[0]
n_ctx_only   = seta_v2[~seta_v2["present_in_SC"] & seta_v2["present_in_CTX"]].shape[0]
n_cross      = seta_v2["cross_tissue"].sum()
n_cross_disc = seta_v2["cross_tissue_discordant"].sum()

print(f"  SetA combined v2: {len(seta_v2)} total unique genes")
print(f"    SC only:              {n_sc_only}")
print(f"    CTX only:             {n_ctx_only}")
print(f"    Cross-tissue same dir:{n_cross} (evidence_count=4)")
print(f"    Cross-tissue discord: {n_cross_disc}")
print(f"  Saved -> {seta_v2_path}")

# Compare to v1
n_v1 = None
new_in_v2 = None
if os.path.isfile(INPUT_SETA_V1):
    seta_v1 = pd.read_csv(INPUT_SETA_V1)
    n_v1 = len(seta_v1)
    v1_genes = set(seta_v1["Gene_Symbol"])
    v2_genes = set(seta_v2["Gene_Symbol"])
    new_in_v2_set = v2_genes - v1_genes
    new_in_v2 = len(new_in_v2_set)
    print(f"  v1 had {n_v1} genes; v2 adds {new_in_v2} new genes")
else:
    print("  (SetA_combined.csv not found — no v1 comparison)")

# ---------------------------------------------------------------------------
# Step 7: Summary report
# ---------------------------------------------------------------------------
print("\n--- Step 7: Writing summary report ---")

lines = [
    "SetA v2 Summary Report",
    "=" * 50,
    "",
    "--- Set E: LCM meta-analysis ---",
]
for key, df in cleaned_lcm.items():
    lines.append(f"  {key}: {len(df)} genes after cleaning")

lines += [
    "",
    f"  SetE FDR05 (primary): {len(sete_fdr05)} genes  (UP={n_up_05}, DOWN={n_down_05})",
    f"  SetE FDR10 total:     {len(sete_fdr10)} genes",
    f"  SetE FDR10 additional (not in FDR05): {len(sete_fdr10_additional)} genes",
    "",
    "--- SC mRNA expanded ---",
    f"  GEO_merged only:      {len(geo_only)}",
    f"  LCM_meta only:        {len(lcm_only)}",
    f"  Both concordant:      {concordant_count}",
    f"  Both discordant:      {discordant_count} (excluded from Set A)",
    f"  Total in expanded:    {len(sc_mrna_expanded)}",
    "",
    "--- Set A v2 gene counts ---",
    f"  SetA spinal cord v2:  {len(seta_sc)} genes  (UP={(seta_sc['direction']=='UP').sum()}, DOWN={(seta_sc['direction']=='DOWN').sum()})",
    f"  SetA cortex v2:       {len(seta_ctx)} genes  (UP={(seta_ctx['direction']=='UP').sum()}, DOWN={(seta_ctx['direction']=='DOWN').sum()})",
    f"  SetA combined v2:     {len(seta_v2)} total unique genes",
    f"    SC only:              {n_sc_only}",
    f"    CTX only:             {n_ctx_only}",
    f"    Cross-tissue same dir:{n_cross} (evidence_count=4)",
    f"    Cross-tissue discord: {n_cross_disc}",
    "",
    "--- Comparison to v1 ---",
]
if n_v1 is not None:
    lines.append(f"  SetA combined v1: {n_v1} genes")
    lines.append(f"  New genes in v2:  {new_in_v2}")
else:
    lines.append("  SetA v1 not found — no comparison available")

summary_path = os.path.join(OUT_DIR, "SetA_v2_summary.txt")
with open(summary_path, "w") as f:
    f.write("\n".join(lines) + "\n")
print(f"  Saved -> {summary_path}")

print("\n=== All done ===")
