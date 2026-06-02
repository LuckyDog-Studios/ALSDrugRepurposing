"""
CMap Query Input Pipeline for ALS Drug Repurposing
Builds two CMap query files from multi-omic evidence layers:
  Query 1: Inflammation-focused (pathway enrichment genes only)
  Query 2: Multi-omic signature (all evidence layers combined)
"""

import os
import sys
import re
import pandas as pd
import numpy as np

# ---------------------------------------------------------------------------
# Output directory
# ---------------------------------------------------------------------------
OUT_DIR = "cmap_prep_results"
os.makedirs(OUT_DIR, exist_ok=True)

# ---------------------------------------------------------------------------
# Input paths
# NOTE: pathway enrichment folder is pathway_enrichment_results/ (not pathway_enrichment/)
# ---------------------------------------------------------------------------
PATHS = {
    "pathway_up":    "pathway_enrichment_results/CMap_input_upregulated_genes.csv",
    "setA":          "multiomics_filtering_results/SetA_combined_v2.csv",
    "setE":          "multiomics_filtering_results/SetE_LCM_FDR05_combined.csv",
    "disc_sc":       "mRNA_protein_corr_results/spinal_cord_correlation_results.csv",
    "disc_ctx":      "mRNA_protein_corr_results/cortex_correlation_results_v2.csv",
    "cytoscape":     "string_results/primary/cytoscape_notes.txt",
}

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def valid_gene(sym):
    """Return True if sym looks like a real HGNC gene symbol."""
    if pd.isna(sym) or str(sym).strip() == "":
        return False
    s = str(sym).strip()
    # Must start with letter, contain only alphanumeric / hyphen / dot
    return bool(re.match(r'^[A-Za-z][A-Za-z0-9\-\.]*$', s))


def clean_genes(series):
    return series.dropna().astype(str).str.strip().loc[lambda s: s.apply(valid_gene)]


def merge_sources(df, gene_col, source_col, n_col):
    """
    Given a DataFrame with Gene_Symbol and source columns, collapse duplicate
    genes by concatenating sources and counting them.
    """
    grouped = (
        df.groupby(gene_col)[source_col]
        .apply(lambda s: "+".join(sorted(set(s))))
        .reset_index()
    )
    grouped[n_col] = grouped[source_col].apply(lambda s: len(s.split("+")))
    return grouped


# ---------------------------------------------------------------------------
# Step 0: Pre-flight check - confirm all input files exist
# ---------------------------------------------------------------------------
print("=" * 70)
print("CMap Input Pipeline - ALS Drug Repurposing")
print("=" * 70)
print("\nExpected input paths:")
all_found = True
for key, path in PATHS.items():
    exists = os.path.exists(path)
    status = "OK" if exists else "MISSING"
    print(f"  [{status}] {path}")
    if not exists:
        all_found = False

if not all_found:
    print("\nERROR: One or more input files are missing. Stopping.")
    sys.exit(1)

print("\nAll input files found. Proceeding.\n")

# ---------------------------------------------------------------------------
# Step 1: Load and inspect all input files
# ---------------------------------------------------------------------------
print("=" * 70)
print("STEP 1: Load and inspect input files")
print("=" * 70)

inspection_lines = []

def inspect(label, df, gene_col):
    lines = [
        f"\n--- {label} ---",
        f"  Shape:   {df.shape}",
        f"  Columns: {list(df.columns)}",
        f"  Head (3 rows):\n{df.head(3).to_string(index=False)}",
        f"  N unique gene symbols: {df[gene_col].dropna().nunique()}",
    ]
    for ln in lines:
        print(ln)
    inspection_lines.extend(lines)


# Pathway enrichment
df_pathway = pd.read_csv(PATHS["pathway_up"])
inspect("pathway_enrichment_results/CMap_input_upregulated_genes.csv", df_pathway, "Gene_Symbol")

# Set A
df_setA = pd.read_csv(PATHS["setA"])
inspect("multiomics_filtering_results/SetA_combined_v2.csv", df_setA, "Gene_Symbol")

# Set E
df_setE = pd.read_csv(PATHS["setE"])
inspect("multiomics_filtering_results/SetE_LCM_FDR05_combined.csv", df_setE, "Gene_Symbol")

# Discordant LOW - spinal cord
df_sc = pd.read_csv(PATHS["disc_sc"])
inspect("mRNA_protein_corr_results/spinal_cord_correlation_results.csv", df_sc, "Gene_Symbol")

# Discordant LOW - cortex
df_ctx = pd.read_csv(PATHS["disc_ctx"])
inspect("mRNA_protein_corr_results/cortex_correlation_results_v2.csv", df_ctx, "Gene_Symbol")

# Cytoscape notes - parse separately
with open(PATHS["cytoscape"], "r") as fh:
    cyto_text = fh.read()

def parse_gene_list(text, header_pattern):
    """
    Extract gene symbols following a header line matching header_pattern.
    Stops at the next blank line after at least one gene is collected,
    or at the next section header.
    """
    lines = text.splitlines()
    in_section = False
    genes = []
    for line in lines:
        stripped = line.strip()
        if re.search(header_pattern, stripped, re.IGNORECASE):
            in_section = True
            continue
        if in_section:
            if stripped == "":
                if genes:          # blank line after collecting genes = end of section
                    break
                continue           # blank line before first gene - keep going
            # Stop if we hit another section header
            if re.search(r'^top\s+\d+', stripped, re.IGNORECASE):
                break
            if valid_gene(stripped):
                genes.append(stripped)
    return genes

betweenness_genes = parse_gene_list(cyto_text, r'betweeness|betweenness')
degree_genes      = parse_gene_list(cyto_text, r'top\s+10\s+degree')
hub_genes         = sorted(set(betweenness_genes) | set(degree_genes))

cyto_lines = [
    f"\n--- {PATHS['cytoscape']} ---",
    f"  Betweenness centrality genes ({len(betweenness_genes)}): {betweenness_genes}",
    f"  Degree genes ({len(degree_genes)}): {degree_genes}",
    f"  Hub gene union ({len(hub_genes)}): {hub_genes}",
]
for ln in cyto_lines:
    print(ln)
inspection_lines.extend(cyto_lines)

# Save inspection report
with open(os.path.join(OUT_DIR, "00_inspection.txt"), "w") as fh:
    fh.write("\n".join(inspection_lines))
print(f"\nInspection saved to {OUT_DIR}/00_inspection.txt")

# ---------------------------------------------------------------------------
# Step 2: Build Query 1 - Inflammation-focused
# ---------------------------------------------------------------------------
print("\n" + "=" * 70)
print("STEP 2: Build Query 1 - Inflammation-focused")
print("=" * 70)

# Upregulated - pathway enrichment genes
q1_up = df_pathway[["Gene_Symbol"]].copy()
q1_up["Gene_Symbol"] = q1_up["Gene_Symbol"].astype(str).str.strip()
q1_up = q1_up[q1_up["Gene_Symbol"].apply(valid_gene)].drop_duplicates("Gene_Symbol")
q1_up["direction"] = "UP"
q1_up["source"]    = "pathway_enrichment"
q1_up = q1_up[["Gene_Symbol", "direction", "source"]].reset_index(drop=True)

q1_up.to_csv(os.path.join(OUT_DIR, "Query1_upregulated.csv"), index=False)
with open(os.path.join(OUT_DIR, "Query1_upregulated.txt"), "w") as fh:
    fh.write("\n".join(q1_up["Gene_Symbol"].tolist()))
print(f"Query 1 UP:   {len(q1_up)} genes  ->  Query1_upregulated.csv / .txt")

# Downregulated - none for Query 1
with open(os.path.join(OUT_DIR, "Query1_downregulated.txt"), "w") as fh:
    fh.write("# No downregulated genes in Query 1")
print("Query 1 DOWN: 0 genes  ->  Query1_downregulated.txt (placeholder)")

# ---------------------------------------------------------------------------
# Step 3: Build Query 2 - Multi-omic signature
# ---------------------------------------------------------------------------
print("\n" + "=" * 70)
print("STEP 3: Build Query 2 - Multi-omic signature")
print("=" * 70)

# ---- Direction helper for Set A ----
def get_direction(row):
    sc = row.get("direction_SC", np.nan)
    ctx = row.get("direction_CTX", np.nan)
    if pd.notna(sc) and str(sc).strip() not in ("", "nan"):
        return str(sc).strip()
    if pd.notna(ctx) and str(ctx).strip() not in ("", "nan"):
        return str(ctx).strip()
    return np.nan

df_setA["direction"] = df_setA.apply(get_direction, axis=1)

# ---- Query 2 UPREGULATED ----

# Source 1: pathway enrichment (same as Q1)
src1_up = q1_up[["Gene_Symbol"]].copy()
src1_up["source"] = "pathway_enrichment"

# Source 2: Set A UP genes
setA_up = df_setA[df_setA["direction"] == "UP"][["Gene_Symbol"]].copy()
setA_up["Gene_Symbol"] = setA_up["Gene_Symbol"].astype(str).str.strip()
setA_up = setA_up[setA_up["Gene_Symbol"].apply(valid_gene)].drop_duplicates("Gene_Symbol")
setA_up["source"] = "SetA_multiomics_UP"

# Source 3: Network hub genes
src3_up = pd.DataFrame({"Gene_Symbol": hub_genes, "source": "network_hub"})

# Combine
all_up = pd.concat([src1_up, setA_up, src3_up], ignore_index=True)
all_up["Gene_Symbol"] = all_up["Gene_Symbol"].astype(str).str.strip()
all_up = all_up[all_up["Gene_Symbol"].apply(valid_gene)]

q2_up = merge_sources(all_up, "Gene_Symbol", "source", "n_up_sources")
q2_up["direction"] = "UP"
q2_up = q2_up[["Gene_Symbol", "direction", "source", "n_up_sources"]]
q2_up = q2_up.sort_values("n_up_sources", ascending=False).reset_index(drop=True)

n_pe_only   = q2_up[q2_up["source"] == "pathway_enrichment"].shape[0]
n_setA_only = q2_up[q2_up["source"] == "SetA_multiomics_UP"].shape[0]
n_hub_only  = q2_up[q2_up["source"] == "network_hub"].shape[0]
n_multi_up  = q2_up[q2_up["n_up_sources"] >= 2].shape[0]

print(f"\nQuery 2 UP - {len(q2_up)} total genes")
print(f"  Pathway enrichment only: {n_pe_only}")
print(f"  Set A only:              {n_setA_only}")
print(f"  Network hub only:        {n_hub_only}")
print(f"  Supported by 2+ sources: {n_multi_up}")

q2_up.to_csv(os.path.join(OUT_DIR, "Query2_upregulated.csv"), index=False)
with open(os.path.join(OUT_DIR, "Query2_upregulated.txt"), "w") as fh:
    fh.write("\n".join(q2_up["Gene_Symbol"].tolist()))

# ---- Query 2 DOWNREGULATED ----

# Source 1: Set A DOWN genes
setA_down = df_setA[df_setA["direction"] == "DOWN"][["Gene_Symbol"]].copy()
setA_down["Gene_Symbol"] = setA_down["Gene_Symbol"].astype(str).str.strip()
setA_down = setA_down[setA_down["Gene_Symbol"].apply(valid_gene)].drop_duplicates("Gene_Symbol")
setA_down["source"] = "SetA_multiomics_DOWN"

# Source 2: Set E LCM DOWN genes
setE_down = df_setE[df_setE["direction"] == "DOWN"][["Gene_Symbol"]].copy()
setE_down["Gene_Symbol"] = setE_down["Gene_Symbol"].astype(str).str.strip()
setE_down = setE_down[setE_down["Gene_Symbol"].apply(valid_gene)].drop_duplicates("Gene_Symbol")
setE_down["source"] = "LCM_meta_DOWN"

# Source 3: Discordant LOW - spinal cord
disc_sc_down = df_sc[df_sc["quadrant"] == "discordant_LOW"][["Gene_Symbol"]].copy()
disc_sc_down["Gene_Symbol"] = disc_sc_down["Gene_Symbol"].astype(str).str.strip()
disc_sc_down = disc_sc_down[disc_sc_down["Gene_Symbol"].apply(valid_gene)].drop_duplicates("Gene_Symbol")
disc_sc_down["source"] = "discordant_LOW_SC"

# Source 4: Discordant LOW - cortex
disc_ctx_down = df_ctx[df_ctx["quadrant"] == "discordant_LOW"][["Gene_Symbol"]].copy()
disc_ctx_down["Gene_Symbol"] = disc_ctx_down["Gene_Symbol"].astype(str).str.strip()
disc_ctx_down = disc_ctx_down[disc_ctx_down["Gene_Symbol"].apply(valid_gene)].drop_duplicates("Gene_Symbol")
disc_ctx_down["source"] = "discordant_LOW_CTX"

# Combine
all_down = pd.concat([setA_down, setE_down, disc_sc_down, disc_ctx_down], ignore_index=True)
all_down["Gene_Symbol"] = all_down["Gene_Symbol"].astype(str).str.strip()
all_down = all_down[all_down["Gene_Symbol"].apply(valid_gene)]

q2_down = merge_sources(all_down, "Gene_Symbol", "source", "n_down_sources")
q2_down["direction"] = "DOWN"
q2_down = q2_down[["Gene_Symbol", "direction", "source", "n_down_sources"]]
q2_down = q2_down.sort_values("n_down_sources", ascending=False).reset_index(drop=True)

n_setA_d   = q2_down[q2_down["source"] == "SetA_multiomics_DOWN"].shape[0]
n_lcm_d    = q2_down[q2_down["source"] == "LCM_meta_DOWN"].shape[0]
n_disc_sc  = q2_down[q2_down["source"] == "discordant_LOW_SC"].shape[0]
n_disc_ctx = q2_down[q2_down["source"] == "discordant_LOW_CTX"].shape[0]
n_multi_down = q2_down[q2_down["n_down_sources"] >= 2].shape[0]

print(f"\nQuery 2 DOWN - {len(q2_down)} total genes")
print(f"  Set A only:             {n_setA_d}")
print(f"  LCM meta only:          {n_lcm_d}")
print(f"  Discordant LOW SC only: {n_disc_sc}")
print(f"  Discordant LOW CTX only:{n_disc_ctx}")
print(f"  Supported by 2+ sources:{n_multi_down}")

q2_down.to_csv(os.path.join(OUT_DIR, "Query2_downregulated.csv"), index=False)
with open(os.path.join(OUT_DIR, "Query2_downregulated.txt"), "w") as fh:
    fh.write("\n".join(q2_down["Gene_Symbol"].tolist()))

# ---------------------------------------------------------------------------
# Step 4: Conflict resolution
# ---------------------------------------------------------------------------
print("\n" + "=" * 70)
print("STEP 4: Conflict resolution")
print("=" * 70)

up_genes   = set(q2_up["Gene_Symbol"])
down_genes = set(q2_down["Gene_Symbol"])
conflicts  = up_genes & down_genes

if conflicts:
    print(f"\n{len(conflicts)} conflicting gene(s) found (present in both UP and DOWN):")
    for gene in sorted(conflicts):
        up_src   = q2_up.loc[q2_up["Gene_Symbol"] == gene, "source"].values[0]
        down_src = q2_down.loc[q2_down["Gene_Symbol"] == gene, "source"].values[0]
        print(f"  {gene}: UP sources=[{up_src}]  DOWN sources=[{down_src}]  -> REMOVED")
else:
    print("\nNo conflicts found.")

# Remove conflicts from both lists
q2_up_clean   = q2_up[~q2_up["Gene_Symbol"].isin(conflicts)].reset_index(drop=True)
q2_down_clean = q2_down[~q2_down["Gene_Symbol"].isin(conflicts)].reset_index(drop=True)

# Overwrite files with conflict-free versions
q2_up_clean.to_csv(os.path.join(OUT_DIR, "Query2_upregulated.csv"), index=False)
with open(os.path.join(OUT_DIR, "Query2_upregulated.txt"), "w") as fh:
    fh.write("\n".join(q2_up_clean["Gene_Symbol"].tolist()))

q2_down_clean.to_csv(os.path.join(OUT_DIR, "Query2_downregulated.csv"), index=False)
with open(os.path.join(OUT_DIR, "Query2_downregulated.txt"), "w") as fh:
    fh.write("\n".join(q2_down_clean["Gene_Symbol"].tolist()))

print(f"\nAfter conflict removal:")
print(f"  Query 2 UP:   {len(q2_up_clean)} genes")
print(f"  Query 2 DOWN: {len(q2_down_clean)} genes")

# ---------------------------------------------------------------------------
# Step 5: Overlap analysis - Query 1 vs Query 2 upregulated
# ---------------------------------------------------------------------------
print("\n" + "=" * 70)
print("STEP 5: Overlap analysis - Query 1 vs Query 2 upregulated")
print("=" * 70)

q1_genes = set(q1_up["Gene_Symbol"])
q2_genes = set(q2_up_clean["Gene_Symbol"])

overlap   = q1_genes & q2_genes
unique_q2 = q2_genes - q1_genes

print(f"  Genes in both Q1 and Q2 UP:  {len(overlap)}")
print(f"  Unique to Q2 UP:             {len(unique_q2)}")

# Build overlap table
overlap_rows = []
all_genes_union = q1_genes | q2_genes
for gene in sorted(all_genes_union):
    in_q1 = gene in q1_genes
    in_q2 = gene in q2_genes
    q2_src = q2_up_clean.loc[q2_up_clean["Gene_Symbol"] == gene, "source"].values
    q2_src_str = q2_src[0] if len(q2_src) > 0 else ""
    overlap_rows.append({
        "Gene_Symbol":     gene,
        "direction":       "UP",
        "present_in_Q1":   in_q1,
        "present_in_Q2":   in_q2,
        "Q2_sources":      q2_src_str,
    })

df_overlap = pd.DataFrame(overlap_rows)
df_overlap.to_csv(os.path.join(OUT_DIR, "Query1_Query2_overlap.csv"), index=False)
print(f"  Overlap table saved -> {OUT_DIR}/Query1_Query2_overlap.csv")

# ---------------------------------------------------------------------------
# Step 6: CMap-ready summary files
# ---------------------------------------------------------------------------
print("\n" + "=" * 70)
print("STEP 6: Generate CMap-ready input files")
print("=" * 70)

# Query 1 input
q1_cmap_path = os.path.join(OUT_DIR, "CMAP_Query1_input.txt")
with open(q1_cmap_path, "w") as fh:
    fh.write("# Query 1: ALS Inflammation Signature (pathway enrichment)\n")
    fh.write(f"# Upregulated genes: {len(q1_up)}\n")
    fh.write("# Downregulated genes: 0\n")
    fh.write("\n")
    fh.write("\n".join(q1_up["Gene_Symbol"].tolist()))
print(f"  CMAP_Query1_input.txt              ({len(q1_up)} UP genes)")

# Query 2 upregulated
q2_up_cmap_path = os.path.join(OUT_DIR, "CMAP_Query2_upregulated.txt")
with open(q2_up_cmap_path, "w") as fh:
    fh.write("\n".join(q2_up_clean["Gene_Symbol"].tolist()))
print(f"  CMAP_Query2_upregulated.txt        ({len(q2_up_clean)} genes)")

# Query 2 downregulated
q2_down_cmap_path = os.path.join(OUT_DIR, "CMAP_Query2_downregulated.txt")
with open(q2_down_cmap_path, "w") as fh:
    fh.write("\n".join(q2_down_clean["Gene_Symbol"].tolist()))
print(f"  CMAP_Query2_downregulated.txt      ({len(q2_down_clean)} genes)")

# ---------------------------------------------------------------------------
# Step 7: Summary report
# ---------------------------------------------------------------------------
print("\n" + "=" * 70)
print("STEP 7: Summary report")
print("=" * 70)

top10_up   = q2_up_clean.nlargest(10, "n_up_sources")[["Gene_Symbol", "source", "n_up_sources"]]
top10_down = q2_down_clean.nlargest(10, "n_down_sources")[["Gene_Symbol", "source", "n_down_sources"]]

# Recompute per-source counts on clean lists (for accurate per-source singles)
def count_source_only(df, src_col, src_name):
    return df[df[src_col] == src_name].shape[0]

n_pe_only_c   = count_source_only(q2_up_clean,   "source", "pathway_enrichment")
n_setA_uc     = count_source_only(q2_up_clean,   "source", "SetA_multiomics_UP")
n_hub_uc      = count_source_only(q2_up_clean,   "source", "network_hub")
n_multi_up_c  = q2_up_clean[q2_up_clean["n_up_sources"] >= 2].shape[0]

n_setA_dc     = count_source_only(q2_down_clean,  "source", "SetA_multiomics_DOWN")
n_lcm_dc      = count_source_only(q2_down_clean,  "source", "LCM_meta_DOWN")
n_disc_sc_c   = count_source_only(q2_down_clean,  "source", "discordant_LOW_SC")
n_disc_ctx_c  = count_source_only(q2_down_clean,  "source", "discordant_LOW_CTX")
n_multi_down_c = q2_down_clean[q2_down_clean["n_down_sources"] >= 2].shape[0]

summary = f"""
CMap Input Pipeline - Summary Report
=====================================
Date: {pd.Timestamp.today().strftime('%Y-%m-%d')}

QUERY 1 - ALS Inflammation Signature (pathway enrichment)
  Upregulated genes:   {len(q1_up)}
  Downregulated genes: 0
  Sources: pathway_enrichment (neutrophil degranulation pathway)

QUERY 2 - Multi-omic ALS Signature
  Upregulated genes:   {len(q2_up_clean)}
    Pathway enrichment only:   {n_pe_only_c}
    Set A multiomics only:     {n_setA_uc}
    Network hub only:          {n_hub_uc}
    Supported by 2+ sources:   {n_multi_up_c}  (highest confidence UP)

  Downregulated genes: {len(q2_down_clean)}
    Set A multiomics only:     {n_setA_dc}
    LCM meta-analysis only:    {n_lcm_dc}
    Discordant LOW SC only:    {n_disc_sc_c}
    Discordant LOW CTX only:   {n_disc_ctx_c}
    Supported by 2+ sources:   {n_multi_down_c}  (highest confidence DOWN)

CONFLICT RESOLUTION
  Conflicts removed from Query 2: {len(conflicts)}
  {"Removed genes: " + ", ".join(sorted(conflicts)) if conflicts else "No conflicts detected."}

QUERY OVERLAP
  Genes shared between Q1 UP and Q2 UP: {len(overlap)}
  Genes unique to Q2 UP (multi-omic addition): {len(unique_q2)}

TOP 10 MOST SUPPORTED UPREGULATED GENES (Query 2):
{top10_up.to_string(index=False)}

TOP 10 MOST SUPPORTED DOWNREGULATED GENES (Query 2):
{top10_down.to_string(index=False)}

SUBMISSION INSTRUCTIONS
  Submit CMAP_Query1_input.txt and
  CMAP_Query2_upregulated.txt + CMAP_Query2_downregulated.txt
  separately to clue.io L1000 query tool. Compare top negative
  connectivity score drugs between both queries. Drugs appearing
  in top hits of both queries are highest priority candidates.
"""

print(summary)
with open(os.path.join(OUT_DIR, "cmap_prep_summary.txt"), "w") as fh:
    fh.write(summary)

print(f"Summary saved -> {OUT_DIR}/cmap_prep_summary.txt")
print("\nPipeline complete. All outputs written to:", OUT_DIR)
