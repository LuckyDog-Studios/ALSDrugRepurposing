"""
network_basis_pipeline.py

Compiles the final STRING network seed list from four evidence layers for ALS
drug-repurposing analysis.

Evidence sets:
  Set A  -- multi-omic concordant hits (mRNA + protein, same tissue)
  Set B  -- proteostatic failure / discordant_LOW UPS signature
  Set C  -- RNA-independent G3BP1 stress-granule interactors (C9orf72 iPSC MNs)
  Set E  -- LCM motor-neuron meta-analysis (Ziff et al. 2024)
           *** Set E ONLY boosts evidence for genes already in A, B, or C ***
           Set E genes not found in A/B/C are excluded from the network.
"""

# pip install pandas numpy  (both are standard in most scientific Python environments)

import os
import sys
import textwrap
import pandas as pd
import numpy as np

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(
    os.path.abspath(__file__)))))   # repo root

INPUT_PATHS = {
    "SetA": os.path.join(ROOT, "multiomics_filtering_results", "SetA_combined_v2.csv"),
    "SetB": os.path.join(ROOT, "mRNA_protein_corr_results", "ALS_Ubiquitome_STRING_input_extended.csv"),
    "SetC": os.path.join(ROOT, "mRNA_protein_corr_results", "archive_v1", "PXD065424_highconfidence_interactome.csv"),
    "SetE": os.path.join(ROOT, "multiomics_filtering_results", "SetE_LCM_FDR05_combined.csv"),
}

OUT_DIR = os.path.join(ROOT, "network_basis_results")
os.makedirs(OUT_DIR, exist_ok=True)

# ---------------------------------------------------------------------------
# Biological reference lists (Step 6)
# ---------------------------------------------------------------------------
UPS_GENES = {
    "PSMD1","PSMD2","PSMD3","PSMD4","PSMD7","PSMD11","PSMD14",
    "PSMB1","PSMB2","PSMB3","PSMB4","PSMB5","PSMB6","PSMB7",
    "PSMA1","PSMA2","PSMA3","PSMA4","PSMA5","PSMA6","PSMA7",
    "UBA1","UBA2","UBA3","UBB","UBC","UCHL1","UCHL5",
    "SQSTM1","UBQLN1","UBQLN2","HDAC6","VCP","OPTN","TBK1",
}
ALS_GENES = {
    "SOD1","TARDBP","FUS","C9orf72","UBQLN2","OPTN","TBK1","VCP",
    "SQSTM1","HNRNPA1","HNRNPA2B1","MATR3","TIA1","ATXN2","NEK1",
    "CCNF","KIF5A","SETX","ANG","VAPB","SIGMAR1",
}
SG_RBP_GENES = {
    "G3BP1","G3BP2","TIA1","TIAL1","FUS","TARDBP","HNRNPA1",
    "HNRNPA2B1","ATXN2","STMN2","SFPQ","MATR3",
}
NEUROINFLAM_GENES = {
    "CFI","CP","GBP2","ITGB2","SERPINA3",
    "ANPEP","CD34","CFB","ICAM1","HLA-C","GPX3","GSTM1",
}


# ===========================================================================
# Helpers
# ===========================================================================

def load_csv(label: str, path: str) -> pd.DataFrame:
    """Load a CSV; exit with a clear message if the file is missing."""
    if not os.path.exists(path):
        print(f"\n[ERROR] Cannot find {label} at:\n  {path}\nExiting.")
        sys.exit(1)
    df = pd.read_csv(path)
    return df


def write_lines(path: str, items) -> None:
    """Write one item per line, no header."""
    with open(path, "w") as fh:
        for item in items:
            fh.write(str(item) + "\n")


def bool_col(series: pd.Series) -> pd.Series:
    """Normalise a column to True/False (handles bool and string representations)."""
    if series.dtype == bool:
        return series
    return series.map(lambda x: str(x).strip().upper() in {"TRUE", "1", "YES"})


# ===========================================================================
# Step 0: Confirm input files exist
# ===========================================================================
print("=" * 70)
print("NETWORK BASIS PIPELINE -- ALS STRING seed list")
print("=" * 70)
print("\nExpected input files:")
all_found = True
for label, path in INPUT_PATHS.items():
    rel = os.path.relpath(path, ROOT)
    exists = os.path.exists(path)
    status = "OK" if exists else "MISSING"
    print(f"  [{status}]  {label}: {rel}")
    if not exists:
        all_found = False

if not all_found:
    print("\n[ERROR] One or more input files are missing. "
          "Please check paths and re-run.")
    sys.exit(1)

print("\nAll input files found. Proceeding.\n")

# ===========================================================================
# Step 1: Load and inspect all input files
# ===========================================================================
print("=" * 70)
print("STEP 1 -- Load and inspect input files")
print("=" * 70)

inspection_lines = []

def inspect(label, df):
    lines = [
        f"=== {label} ===",
        f"Shape: {df.shape}",
        f"Columns: {list(df.columns)}",
        "First 3 rows:",
        df.head(3).to_string(index=False),
        f"N unique Gene_Symbol: {df['Gene_Symbol'].nunique()}",
        "",
    ]
    for line in lines:
        print(line)
    inspection_lines.extend(lines)

raw = {label: load_csv(label, path) for label, path in INPUT_PATHS.items()}
for label, df in raw.items():
    inspect(label, df)

inspection_path = os.path.join(OUT_DIR, "00_input_inspection.txt")
with open(inspection_path, "w") as fh:
    fh.write("\n".join(inspection_lines))
print(f"Saved inspection -> {os.path.relpath(inspection_path, ROOT)}\n")


# ===========================================================================
# Step 2: Extract gene lists with metadata
# ===========================================================================
print("=" * 70)
print("STEP 2 -- Extract gene lists")
print("=" * 70)

# -- Set A -----------------------------------------------------------------
dfA_raw = raw["SetA"].copy()

setA = pd.DataFrame({
    "Gene_Symbol": dfA_raw["Gene_Symbol"],
    "SetA_present_in_SC":  bool_col(dfA_raw["present_in_SC"]),
    "SetA_present_in_CTX": bool_col(dfA_raw["present_in_CTX"]),
    # cross-tissue defined as evidence_count == 4
    "SetA_cross_tissue":   (dfA_raw["evidence_count"] == 4),
    "SetA_direction_SC":   dfA_raw["direction_SC"].where(dfA_raw["present_in_SC"], other=pd.NA),
    "SetA_direction_CTX":  dfA_raw["direction_CTX"].where(dfA_raw["present_in_CTX"], other=pd.NA),
    "SetA_mRNA_source":    dfA_raw["SC_mRNA_source"],
})
setA = setA.drop_duplicates("Gene_Symbol").reset_index(drop=True)
print(f"Set A: {len(setA)} genes loaded")

# -- Set B -----------------------------------------------------------------
dfB_raw = raw["SetB"].copy()

setB = pd.DataFrame({"Gene_Symbol": dfB_raw["Gene_Symbol"]})
if "present_in_spinal" in dfB_raw.columns and "present_in_cortex" in dfB_raw.columns:
    setB["SetB_in_SC"]  = bool_col(dfB_raw["present_in_spinal"])
    setB["SetB_in_CTX"] = bool_col(dfB_raw["present_in_cortex"])
    print("Set B: columns 'present_in_spinal' and 'present_in_cortex' found -- "
          "using for SetB_in_SC / SetB_in_CTX")
else:
    print(f"Set B: tissue columns not found. Available: {list(dfB_raw.columns)}")
    print("  -> Setting SetB_in_SC = NA and SetB_in_CTX = NA")
    setB["SetB_in_SC"]  = pd.NA
    setB["SetB_in_CTX"] = pd.NA

setB = setB.drop_duplicates("Gene_Symbol").reset_index(drop=True)
print(f"Set B: {len(setB)} genes loaded")

# -- Set C -----------------------------------------------------------------
dfC_raw = raw["SetC"].copy()
setC = dfC_raw[["Gene_Symbol"]].drop_duplicates().reset_index(drop=True)
print(f"Set C: {len(setC)} genes loaded")

# -- Set E -----------------------------------------------------------------
dfE_raw = raw["SetE"].copy()
setE = pd.DataFrame({
    "Gene_Symbol":    dfE_raw["Gene_Symbol"],
    "SetE_direction": dfE_raw["direction"],
    "SetE_SMD":       dfE_raw["LCM_SMD"],
    "SetE_FDR":       dfE_raw["LCM_FDR"],
})
setE = setE.drop_duplicates("Gene_Symbol").reset_index(drop=True)
print(f"Set E: {len(setE)} genes loaded (total before overlap filter)\n")


# ===========================================================================
# Step 3: Determine Set E contribution
# ===========================================================================
print("=" * 70)
print("STEP 3 -- Set E overlap with Sets A, B, C")
print("=" * 70)

genesA = set(setA["Gene_Symbol"])
genesB = set(setB["Gene_Symbol"])
genesC = set(setC["Gene_Symbol"])
genesE = set(setE["Gene_Symbol"])

E_in_A = genesE & genesA
E_in_B = genesE & genesB
E_in_C = genesE & genesC
SetE_in_network = genesE & (genesA | genesB | genesC)
SetE_excluded   = genesE - (genesA | genesB | genesC)

print(f"Set E genes overlapping Set A: {len(E_in_A)}")
print(f"Set E genes overlapping Set B: {len(E_in_B)}")
print(f"Set E genes overlapping Set C: {len(E_in_C)}")
print(f"Set E genes in network (>=1 of A/B/C): {len(SetE_in_network)}")
print(f"Set E genes EXCLUDED (not in A/B/C): {len(SetE_excluded)}\n")


# ===========================================================================
# Step 4: Build master network seed list
# ===========================================================================
print("=" * 70)
print("STEP 4 -- Build master network seed list")
print("=" * 70)

# Union of A, B, C as base
all_genes = sorted(genesA | genesB | genesC)
master = pd.DataFrame({"Gene_Symbol": all_genes})

# Set membership flags
master["in_SetA"] = master["Gene_Symbol"].isin(genesA)
master["in_SetB"] = master["Gene_Symbol"].isin(genesB)
master["in_SetC"] = master["Gene_Symbol"].isin(genesC)
master["in_SetE"] = master["Gene_Symbol"].isin(SetE_in_network)

# Merge Set A metadata
master = master.merge(setA, on="Gene_Symbol", how="left")

# Merge Set B tissue columns
master = master.merge(
    setB[["Gene_Symbol", "SetB_in_SC", "SetB_in_CTX"]],
    on="Gene_Symbol", how="left"
)

# Merge Set E metadata
master = master.merge(
    setE[["Gene_Symbol", "SetE_direction", "SetE_SMD", "SetE_FDR"]],
    on="Gene_Symbol", how="left"
)

# n_sets and priority tier
master["n_sets"] = (
    master["in_SetA"].astype(int) +
    master["in_SetB"].astype(int) +
    master["in_SetC"].astype(int) +
    master["in_SetE"].astype(int)
)

def assign_tier(n):
    if n >= 3:
        return "Tier1"
    elif n == 2:
        return "Tier2"
    else:
        return "Tier3"

master["priority_tier"] = master["n_sets"].map(assign_tier)

# Sort: n_sets descending, then alphabetical
master = master.sort_values(["n_sets", "Gene_Symbol"], ascending=[False, True])
master = master.reset_index(drop=True)

# Reorder columns to spec
col_order = [
    "Gene_Symbol",
    "in_SetA", "in_SetB", "in_SetC", "in_SetE",
    "SetA_present_in_SC", "SetA_present_in_CTX", "SetA_cross_tissue",
    "SetA_direction_SC", "SetA_direction_CTX", "SetA_mRNA_source",
    "SetB_in_SC", "SetB_in_CTX",
    "SetE_direction", "SetE_SMD", "SetE_FDR",
    "n_sets", "priority_tier",
]
master = master[col_order]

# Save
master_path = os.path.join(OUT_DIR, "network_seed_master.csv")
master.to_csv(master_path, index=False)

tier1 = master[master["priority_tier"] == "Tier1"]
tier2 = master[master["priority_tier"] == "Tier2"]
tier3 = master[master["priority_tier"] == "Tier3"]

print(f"Total unique genes in master list: {len(master)}")
print(f"  Genes from Set A contributing:  {master['in_SetA'].sum()}")
print(f"  Genes from Set B contributing:  {master['in_SetB'].sum()}")
print(f"  Genes from Set C contributing:  {master['in_SetC'].sum()}")
print(f"  Genes from Set E contributing:  {master['in_SetE'].sum()}")
print(f"  Tier 1 (n_sets >= 3): {len(tier1)}")
print(f"  Tier 2 (n_sets == 2): {len(tier2)}")
print(f"  Tier 3 (n_sets == 1): {len(tier3)}")

print("\n--- Tier 1 genes ---")
for _, row in tier1.iterrows():
    sets = "+".join([s for s, f in [("A", row["in_SetA"]), ("B", row["in_SetB"]),
                                     ("C", row["in_SetC"]), ("E", row["in_SetE"])] if f])
    print(f"  {row['Gene_Symbol']:15s}  sets={sets}  n={row['n_sets']}")

print("\n--- Tier 2 genes ---")
for _, row in tier2.iterrows():
    sets = "+".join([s for s, f in [("A", row["in_SetA"]), ("B", row["in_SetB"]),
                                     ("C", row["in_SetC"]), ("E", row["in_SetE"])] if f])
    print(f"  {row['Gene_Symbol']:15s}  sets={sets}  n={row['n_sets']}")

print(f"\nSaved master list -> {os.path.relpath(master_path, ROOT)}\n")


# ===========================================================================
# Step 5: Generate STRING input files
# ===========================================================================
print("=" * 70)
print("STEP 5 -- Generate STRING input files")
print("=" * 70)

# Primary: Tier 1 + Tier 2 (n_sets >= 2)
primary = master[master["n_sets"] >= 2].copy()
primary_genes = primary["Gene_Symbol"].tolist()

primary_txt = os.path.join(OUT_DIR, "STRING_input_primary.txt")
write_lines(primary_txt, primary_genes)
print(f"Primary STRING input ({len(primary_genes)} genes) -> "
      f"{os.path.relpath(primary_txt, ROOT)}")

# Annotated primary (biological annotations added in Step 6, saved there)
# Extended: all genes
extended_genes = master["Gene_Symbol"].tolist()
extended_txt = os.path.join(OUT_DIR, "STRING_input_extended.txt")
write_lines(extended_txt, extended_genes)
print(f"Extended STRING input ({len(extended_genes)} genes) -> "
      f"{os.path.relpath(extended_txt, ROOT)}\n")


# ===========================================================================
# Step 6: Biological annotation check
# ===========================================================================
print("=" * 70)
print("STEP 6 -- Biological annotation check")
print("=" * 70)

primary_set = set(primary_genes)

def report_overlap(category: str, ref_genes: set, gene_set: set) -> list:
    hits = sorted(gene_set & ref_genes)
    pct  = 100 * len(hits) / len(gene_set) if gene_set else 0
    print(f"\n{category}:")
    print(f"  N overlapping: {len(hits)} / {len(gene_set)} "
          f"({pct:.1f}% of primary input)")
    if hits:
        print("  Genes: " + ", ".join(hits))
    return hits

ups_hits    = report_overlap("UPS / proteostasis",    UPS_GENES,      primary_set)
als_hits    = report_overlap("Known ALS genes",       ALS_GENES,      primary_set)
sg_hits     = report_overlap("Stress granule / RBP",  SG_RBP_GENES,   primary_set)
neuro_hits  = report_overlap("Neuroinflammation",     NEUROINFLAM_GENES, primary_set)

# Add annotation columns to primary annotated table
primary = primary.copy()
primary["UPS_proteostasis"]      = primary["Gene_Symbol"].isin(UPS_GENES)
primary["known_ALS_gene"]        = primary["Gene_Symbol"].isin(ALS_GENES)
primary["stress_granule_RBP"]    = primary["Gene_Symbol"].isin(SG_RBP_GENES)
primary["neuroinflammation"]     = primary["Gene_Symbol"].isin(NEUROINFLAM_GENES)
primary["n_annotation_categories"] = (
    primary["UPS_proteostasis"].astype(int) +
    primary["known_ALS_gene"].astype(int) +
    primary["stress_granule_RBP"].astype(int) +
    primary["neuroinflammation"].astype(int)
)

primary_ann_path = os.path.join(OUT_DIR, "STRING_input_primary_annotated.csv")
primary.to_csv(primary_ann_path, index=False)
print(f"\nAnnotated primary table -> {os.path.relpath(primary_ann_path, ROOT)}\n")


# ===========================================================================
# Step 7: Summary report
# ===========================================================================
print("=" * 70)
print("STEP 7 -- Writing summary report")
print("=" * 70)

def sets_str(row):
    return "+".join([s for s, f in [("A", row["in_SetA"]), ("B", row["in_SetB"]),
                                     ("C", row["in_SetC"]), ("E", row["in_SetE"])] if f])

summary_lines = [
    "=" * 70,
    "NETWORK BASIS PIPELINE -- Summary Report",
    "=" * 70,
    "",
    "INPUT SETS",
    "-" * 40,
    f"Set A (multi-omic concordant):            {len(setA):>4} genes",
    f"Set B (proteostatic failure / UPS):       {len(setB):>4} genes",
    f"Set C (stress-granule interactors):       {len(setC):>4} genes",
    f"Set E (LCM meta-analysis, total):         {len(setE):>4} genes",
    "",
    "SET E CONTRIBUTION",
    "-" * 40,
    f"Set E overlapping Set A:                  {len(E_in_A):>4}",
    f"Set E overlapping Set B:                  {len(E_in_B):>4}",
    f"Set E overlapping Set C:                  {len(E_in_C):>4}",
    f"Set E contributing to network (>=1 of A/B/C): {len(SetE_in_network):>3}",
    f"Set E excluded (not in A/B/C):            {len(SetE_excluded):>4}",
    "",
    "MASTER SEED LIST",
    "-" * 40,
    f"Total unique genes:                       {len(master):>4}",
    f"  from Set A:                             {master['in_SetA'].sum():>4}",
    f"  from Set B:                             {master['in_SetB'].sum():>4}",
    f"  from Set C:                             {master['in_SetC'].sum():>4}",
    f"  from Set E (boosting):                  {master['in_SetE'].sum():>4}",
    "",
    f"Tier 1 (n_sets >= 3, highest confidence): {len(tier1):>4}",
    f"Tier 2 (n_sets == 2):                     {len(tier2):>4}",
    f"Tier 3 (n_sets == 1):                     {len(tier3):>4}",
    "",
    "STRING INPUT FILES",
    "-" * 40,
    f"Primary input (Tier 1+2, n_sets >= 2):    {len(primary_genes):>4} genes",
    f"Extended input (all tiers):               {len(extended_genes):>4} genes",
    "",
    "TIER 1 GENES (n_sets >= 3)",
    "-" * 40,
]
for _, row in tier1.iterrows():
    summary_lines.append(
        f"  {row['Gene_Symbol']:15s}  sets={sets_str(row):10s}  n={row['n_sets']}"
    )

summary_lines += ["", "TIER 2 GENES (n_sets == 2)", "-" * 40]
for _, row in tier2.iterrows():
    summary_lines.append(
        f"  {row['Gene_Symbol']:15s}  sets={sets_str(row):10s}  n={row['n_sets']}"
    )

summary_lines += [
    "",
    "BIOLOGICAL ANNOTATION OVERLAPS (primary STRING input)",
    "-" * 40,
    f"UPS / proteostasis ({len(ups_hits)} genes):   " +
    (", ".join(ups_hits) if ups_hits else "none"),
    f"Known ALS genes ({len(als_hits)} genes):      " +
    (", ".join(als_hits) if als_hits else "none"),
    f"Stress granule / RBP ({len(sg_hits)} genes):  " +
    (", ".join(sg_hits) if sg_hits else "none"),
    f"Neuroinflammation ({len(neuro_hits)} genes):  " +
    (", ".join(neuro_hits) if neuro_hits else "none"),
    "",
    "USAGE NOTE",
    "-" * 40,
    textwrap.fill(
        "Submit STRING_input_primary.txt to STRING database at string-db.org. "
        "Use confidence threshold 0.7, organism Human (9606), network type full "
        "STRING network. After export, import STRING_input_primary_annotated.csv "
        "into Cytoscape as node annotation table for coloring nodes by evidence "
        "layer (in_SetA, in_SetB, in_SetC, in_SetE) and biological category.",
        width=70,
    ),
    "",
]

summary_path = os.path.join(OUT_DIR, "network_basis_summary.txt")
with open(summary_path, "w") as fh:
    fh.write("\n".join(summary_lines))

print("\n".join(summary_lines))
print(f"\nSaved summary -> {os.path.relpath(summary_path, ROOT)}")

print("\n" + "=" * 70)
print("Pipeline complete.")
print(f"All outputs in: {os.path.relpath(OUT_DIR, ROOT)}/")
print("=" * 70)
