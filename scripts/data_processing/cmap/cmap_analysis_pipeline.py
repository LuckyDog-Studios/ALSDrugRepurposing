"""
CMap L1000 Analysis Pipeline for ALS Drug Repurposing
Parses two query_result.gct files and performs cross-query analysis
to identify drug repurposing candidates for ALS.
"""

# pip install pandas numpy matplotlib seaborn scipy

import os
import sys
import re
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
from collections import Counter

# ---------------------------------------------------------------------------
# Config
# ---------------------------------------------------------------------------
CMAP_DIR = "cmap_results"
OUT_DIR  = "cmap_results/parsed"
os.makedirs(OUT_DIR, exist_ok=True)

HUB_GENES = [
    "C3", "CYBB", "ITGB2", "PTPRC",
    "HLA-DRA", "CP", "HMGCS1", "NDUFA12",
    "ABAT", "FECH",
    "ICAM1", "ITGB3", "VTN", "ANPEP",
    "GGT5", "CD34",
]

OUT_COLS = ["pert_iname", "pert_id", "cell_iname", "raw_cs", "norm_cs",
            "fdr_q_nlog10", "tas", "moa", "target_name"]

# ---------------------------------------------------------------------------
# Directory scan
# ---------------------------------------------------------------------------
print("=" * 70)
print("CMap Analysis Pipeline -- ALS Drug Repurposing")
print("=" * 70)
print(f"\nScanning: {CMAP_DIR}/\n")

scan_lines = [f"Scan of {CMAP_DIR}/\n{'='*60}"]

all_gct_files = []
query_result_files = []

for root, dirs, files in os.walk(CMAP_DIR):
    # Skip our own output directory
    dirs[:] = [d for d in dirs if d != "parsed"]
    rel_root = os.path.relpath(root, CMAP_DIR)
    for fname in sorted(files):
        if fname.endswith(".gct") or fname.endswith(".gctx"):
            full = os.path.join(root, fname)
            rel  = os.path.relpath(full, CMAP_DIR)
            all_gct_files.append(full)
            print(f"  [GCT] {rel}")
            scan_lines.append(f"  {rel}")
            # query_result.gct inside arfs/TAG/
            if fname == "query_result.gct" and "arfs" in root and "TAG" in root:
                query_result_files.append(full)

print(f"\nTotal GCT/GCTX files found: {len(all_gct_files)}")
print(f"\nIdentified query_result.gct files ({len(query_result_files)} found):")
for i, p in enumerate(query_result_files, 1):
    print(f"  Query {i}: {os.path.relpath(p, CMAP_DIR)}")

scan_lines += [
    f"\nTotal GCT/GCTX files: {len(all_gct_files)}",
    f"\nIdentified query_result.gct files:",
]
for i, p in enumerate(query_result_files, 1):
    scan_lines.append(f"  Query {i}: {os.path.relpath(p, CMAP_DIR)}")

with open(os.path.join(OUT_DIR, "00_directory_scan.txt"), "w") as fh:
    fh.write("\n".join(scan_lines))

if len(query_result_files) != 2:
    print(f"\nERROR: Expected 2 query_result.gct files, found {len(query_result_files)}.")
    print("Check that both query result folders are present inside cmap_results/")
    sys.exit(1)

# Sort so that query_1 comes before query_2 by folder name
query_result_files = sorted(query_result_files)

print("\n" + "-" * 70)
print("Please confirm the above file paths are correct.")
confirm = input("Proceed? (yes/no): ").strip().lower()
if confirm not in ("yes", "y"):
    print("Aborted by user.")
    sys.exit(0)

Q1_GCT = query_result_files[0]
Q2_GCT = query_result_files[1]
print(f"\nQuery 1 file: {Q1_GCT}")
print(f"Query 2 file: {Q2_GCT}")

# ---------------------------------------------------------------------------
# Helper: GCT parser
# ---------------------------------------------------------------------------
def parse_gct(path):
    """
    Parse a CMap GCT 1.3 file.
    - Skips the version tag and dimensions line (first 2 lines)
    - Treats line 3 as the header
    - Drops the 'desc' metadata row (line 4 in the file, row 0 after read)
    - Replaces -666 sentinel values with NaN
    - Converts numeric columns
    """
    df = pd.read_csv(path, sep="\t", skiprows=2, low_memory=False)

    raw_shape = df.shape
    print(f"  Raw shape (incl. desc row): {raw_shape}")
    print(f"  Columns: {list(df.columns)}")

    # Confirm and drop the desc/metadata row
    first_id = str(df.iloc[0]["id"]).strip().lower()
    assert first_id == "desc", (
        f"Expected 'desc' as first row id, got '{first_id}'. "
        "Check GCT format -- desc row may be missing or offset."
    )
    df = df.iloc[1:].reset_index(drop=True)
    print(f"  Shape after dropping desc row: {df.shape}")

    # Replace CMap null sentinel -666 (string and numeric)
    df = df.replace({"-666": np.nan, -666: np.nan, -666.0: np.nan})
    # Also handle quoted empty strings
    df = df.replace({"": np.nan, '""': np.nan})

    # Convert numeric columns
    numeric_cols = [
        "raw_cs", "norm_cs", "fdr_q_nlog10", "tas", "cc_q75",
        "qc_pass", "is_hiq", "is_ncs_sig", "is_exemplar_sig",
        "is_null_sig", "nsample", "ss_ngene",
    ]
    for col in numeric_cols:
        if col in df.columns:
            before_na = df[col].isna().sum()
            df[col] = pd.to_numeric(df[col], errors="coerce")
            after_na = df[col].isna().sum()
            if after_na > before_na:
                print(f"  WARNING: {after_na - before_na} new NaN in '{col}' "
                      f"after numeric conversion (non-numeric values present)")

    # Verify raw_cs converted cleanly (no unexpected NaN surge)
    n_raw_cs_na = df["raw_cs"].isna().sum()
    print(f"  raw_cs NaN count after conversion: {n_raw_cs_na} "
          f"({'OK' if n_raw_cs_na == 0 else 'CHECK'})")

    return df, raw_shape[0]


# ---------------------------------------------------------------------------
# Helper: quality filters
# ---------------------------------------------------------------------------
def apply_filters(df, label):
    n0 = len(df)

    df_cp = df[df["pert_type"] == "trt_cp"].copy()
    n_cp = len(df_cp)

    df_ex = df_cp[df_cp["is_exemplar_sig"] == 1].copy()
    n_ex = len(df_ex)

    df_qc = df_ex[df_ex["qc_pass"] == 1].copy()
    n_qc = len(df_qc)

    df_tas = df_qc[df_qc["tas"] >= 0.2].copy()
    n_tas = len(df_tas)

    print(f"\n{label} -- quality filtering:")
    print(f"  Total rows in raw GCT:      {n0:>8,}")
    print(f"  After trt_cp:               {n_cp:>8,}")
    print(f"  After is_exemplar_sig==1:   {n_ex:>8,}")
    print(f"  After qc_pass==1:           {n_qc:>8,}")
    print(f"  After tas >= 0.2:           {n_tas:>8,}")

    df_named = df_tas[df_tas["moa"].notna()].copy()
    n_named = len(df_named)
    print(f"  After moa known (named):    {n_named:>8,}")

    print(f"\n  raw_cs range:  {df_tas['raw_cs'].min():.4f} to {df_tas['raw_cs'].max():.4f}")
    print(f"  norm_cs range: {df_tas['norm_cs'].min():.4f} to {df_tas['norm_cs'].max():.4f}")

    return (
        df_tas[OUT_COLS].sort_values("raw_cs").reset_index(drop=True),
        df_named[OUT_COLS].sort_values("raw_cs").reset_index(drop=True),
        {"n_raw": n0, "n_cp": n_cp, "n_ex": n_ex, "n_qc": n_qc,
         "n_tas": n_tas, "n_named": n_named},
    )


# ---------------------------------------------------------------------------
# Helper: hub gene annotation
# ---------------------------------------------------------------------------
def annotate_hubs(df):
    hub_set = set(HUB_GENES)

    def check(tgt):
        if pd.isna(tgt):
            return [], False, 0
        parts = [t.strip() for t in str(tgt).split("|")]
        matched = [t for t in parts if t in hub_set]
        return matched, bool(matched), len(matched)

    results = df["target_name"].apply(check)
    df = df.copy()
    df["matched_hubs"]    = results.apply(lambda x: "|".join(x[0]) if x[0] else "")
    df["hub_target"]      = results.apply(lambda x: x[1])
    df["n_hubs_targeted"] = results.apply(lambda x: x[2])
    return df


# ---------------------------------------------------------------------------
# Step 1: Parse both GCT files
# ---------------------------------------------------------------------------
print("\n" + "=" * 70)
print("STEP 1: Parse GCT files")
print("=" * 70)

print(f"\nParsing Query 1: {Q1_GCT}")
df_q1_raw, n_q1_raw = parse_gct(Q1_GCT)

print(f"\nParsing Query 2: {Q2_GCT}")
df_q2_raw, n_q2_raw = parse_gct(Q2_GCT)

q1_all, q1_named, q1_stats = apply_filters(df_q1_raw, "Query 1")
q2_all, q2_named, q2_stats = apply_filters(df_q2_raw, "Query 2")

q1_all.to_csv(os.path.join(OUT_DIR, "Query1_all_compounds.csv"), index=False)
q2_all.to_csv(os.path.join(OUT_DIR, "Query2_all_compounds.csv"), index=False)
q1_named.to_csv(os.path.join(OUT_DIR, "Query1_named_drugs.csv"), index=False)
q2_named.to_csv(os.path.join(OUT_DIR, "Query2_named_drugs.csv"), index=False)

print("\nStep 1 outputs saved.")

# ---------------------------------------------------------------------------
# Step 2: Score distribution analysis
# ---------------------------------------------------------------------------
print("\n" + "=" * 70)
print("STEP 2: Score distribution analysis")
print("=" * 70)


def score_stats(df, label):
    cs = df["raw_cs"]
    n_strong_neg  = (cs < -0.3).sum()
    n_mod_neg     = (cs < -0.2).sum()
    n_strong_pos  = (cs >  0.3).sum()
    print(f"\n{label}:")
    print(f"  N with raw_cs < -0.3 (strong negative): {n_strong_neg}")
    print(f"  N with raw_cs < -0.2 (moderate negative): {n_mod_neg}")
    print(f"  N with raw_cs > +0.3 (strong positive):  {n_strong_pos}")
    print(f"  Median raw_cs: {cs.median():.4f}")
    print(f"  Mean raw_cs:   {cs.mean():.4f}")
    return n_strong_neg, n_mod_neg, n_strong_pos


q1_sn, q1_mn, q1_sp = score_stats(q1_all, "Query 1")
q2_sn, q2_mn, q2_sp = score_stats(q2_all, "Query 2")


def plot_distribution(df, label, out_path):
    fig, ax = plt.subplots(figsize=(9, 5))
    ax.hist(df["raw_cs"], bins=80, color="steelblue", edgecolor="white", linewidth=0.3)
    ax.axvline(-0.3, color="red",    linestyle="--", linewidth=1.5, label="raw_cs = -0.3")
    ax.axvline(-0.2, color="orange", linestyle="--", linewidth=1.5, label="raw_cs = -0.2")
    ax.set_xlabel("raw_cs (Connectivity Score)", fontsize=12)
    ax.set_ylabel("Count", fontsize=12)
    ax.set_title(f"CMap raw_cs Distribution -- {label}", fontsize=13)
    ax.legend(fontsize=10)
    sns.despine(ax=ax)
    plt.tight_layout()
    plt.savefig(out_path, dpi=150)
    plt.close()
    print(f"  Saved: {out_path}")


plot_distribution(q1_all, "Query 1: ALS Inflammation Signature",
                  os.path.join(OUT_DIR, "Query1_score_distribution.png"))
plot_distribution(q2_all, "Query 2: ALS Multi-omic Signature",
                  os.path.join(OUT_DIR, "Query2_score_distribution.png"))

# ---------------------------------------------------------------------------
# Step 3: Top drug hits per query
# ---------------------------------------------------------------------------
print("\n" + "=" * 70)
print("STEP 3: Top drug hits per query")
print("=" * 70)


def top_hits(df, n, direction, label, query_tag):
    if direction == "neg":
        sub = df.nsmallest(n, "raw_cs").copy()
    else:
        sub = df.nlargest(n, "raw_cs").copy()
    sub.insert(0, "rank", range(1, len(sub) + 1))
    sub["query"] = query_tag
    return sub


q1_top50_neg = top_hits(q1_named, 50, "neg", "Query1", "Q1")
q2_top50_neg = top_hits(q2_named, 50, "neg", "Query2", "Q2")
q1_top50_pos = top_hits(q1_named, 50, "pos", "Query1", "Q1")
q2_top50_pos = top_hits(q2_named, 50, "pos", "Query2", "Q2")

q1_top50_neg.to_csv(os.path.join(OUT_DIR, "Query1_top50_negative.csv"), index=False)
q2_top50_neg.to_csv(os.path.join(OUT_DIR, "Query2_top50_negative.csv"), index=False)
q1_top50_pos.to_csv(os.path.join(OUT_DIR, "Query1_top50_positive.csv"), index=False)
q2_top50_pos.to_csv(os.path.join(OUT_DIR, "Query2_top50_positive.csv"), index=False)

print(f"  Q1 top 50 negative: raw_cs range "
      f"{q1_top50_neg['raw_cs'].min():.4f} to {q1_top50_neg['raw_cs'].max():.4f}")
print(f"  Q2 top 50 negative: raw_cs range "
      f"{q2_top50_neg['raw_cs'].min():.4f} to {q2_top50_neg['raw_cs'].max():.4f}")

# ---------------------------------------------------------------------------
# Step 4: Cross-query overlap analysis
# ---------------------------------------------------------------------------
print("\n" + "=" * 70)
print("STEP 4: Cross-query overlap analysis")
print("=" * 70)


def build_overlap(q1_df, q2_df, n, q1_named, q2_named):
    """
    Find pert_inames in top-N of both queries.
    For each match return the best (most negative) row from each query.
    """
    q1_drugs = set(q1_df["pert_iname"].str.upper())
    q2_drugs = set(q2_df["pert_iname"].str.upper())
    common   = q1_drugs & q2_drugs

    # Best row per drug = most negative raw_cs
    q1_best = (q1_named.copy()
               .assign(_upper=q1_named["pert_iname"].str.upper())
               .sort_values("raw_cs")
               .drop_duplicates("_upper")
               .set_index("_upper"))
    q2_best = (q2_named.copy()
               .assign(_upper=q2_named["pert_iname"].str.upper())
               .sort_values("raw_cs")
               .drop_duplicates("_upper")
               .set_index("_upper"))

    rows = []
    for drug_upper in sorted(common):
        r1 = q1_best.loc[drug_upper]
        r2 = q2_best.loc[drug_upper]
        # Rank within the original top-N list
        q1_rank_series = q1_df["pert_iname"].str.upper().reset_index(drop=True)
        q2_rank_series = q2_df["pert_iname"].str.upper().reset_index(drop=True)
        q1_rank = q1_rank_series[q1_rank_series == drug_upper].index[0] + 1
        q2_rank = q2_rank_series[q2_rank_series == drug_upper].index[0] + 1
        rows.append({
            "pert_iname":   r1["pert_iname"],
            "Q1_raw_cs":    r1["raw_cs"],
            "Q2_raw_cs":    r2["raw_cs"],
            "mean_cs":      (r1["raw_cs"] + r2["raw_cs"]) / 2,
            "moa":          r1["moa"],
            "target_name":  r1["target_name"],
            "Q1_rank":      int(q1_rank),
            "Q2_rank":      int(q2_rank),
        })

    df_overlap = pd.DataFrame(rows).sort_values("mean_cs").reset_index(drop=True)
    return df_overlap


# Need top-N drug lists for ranking
def top_n_drugs(df_named, n):
    return df_named.nsmallest(n, "raw_cs")

q1_top25  = top_n_drugs(q1_named, 25)
q2_top25  = top_n_drugs(q2_named, 25)
q1_top50  = top_n_drugs(q1_named, 50)
q2_top50  = top_n_drugs(q2_named, 50)
q1_top100 = top_n_drugs(q1_named, 100)
q2_top100 = top_n_drugs(q2_named, 100)

overlap25  = build_overlap(q1_top25,  q2_top25,  25,  q1_named, q2_named)
overlap50  = build_overlap(q1_top50,  q2_top50,  50,  q1_named, q2_named)
overlap100 = build_overlap(q1_top100, q2_top100, 100, q1_named, q2_named)

overlap25.to_csv(os.path.join(OUT_DIR, "overlap_top25_both_queries.csv"),  index=False)
overlap50.to_csv(os.path.join(OUT_DIR, "overlap_top50_both_queries.csv"),  index=False)
overlap100.to_csv(os.path.join(OUT_DIR, "overlap_top100_both_queries.csv"), index=False)

print(f"  Drugs in top 25 overlap:  {len(overlap25)}")
print(f"  Drugs in top 50 overlap:  {len(overlap50)}")
print(f"  Drugs in top 100 overlap: {len(overlap100)}")

if len(overlap25):
    print(f"\n  Top 25 overlap drugs: {', '.join(overlap25['pert_iname'].tolist())}")

# ---------------------------------------------------------------------------
# Step 5: Hub gene annotation
# ---------------------------------------------------------------------------
print("\n" + "=" * 70)
print("STEP 5: Network hub gene cross-reference")
print("=" * 70)

overlap100_hub = annotate_hubs(overlap100)
q1_top50_hub   = annotate_hubs(q1_top50_neg)
q2_top50_hub   = annotate_hubs(q2_top50_neg)

overlap100_hub.to_csv(os.path.join(OUT_DIR, "overlap_top100_hub_annotated.csv"), index=False)
q1_top50_hub.to_csv(os.path.join(OUT_DIR, "Query1_top50_hub_annotated.csv"),    index=False)
q2_top50_hub.to_csv(os.path.join(OUT_DIR, "Query2_top50_hub_annotated.csv"),    index=False)


def report_hub_hits(df, label):
    hits = df[df["hub_target"] == True]
    print(f"\n  {label}: {len(hits)} drug(s) targeting hub genes")
    for _, r in hits.iterrows():
        name = r.get("pert_iname", r.get("Gene_Symbol", "?"))
        print(f"    {name:30s}  hubs={r['matched_hubs']}  n={r['n_hubs_targeted']}")


report_hub_hits(overlap100_hub, "Overlap top 100")
report_hub_hits(q1_top50_hub,   "Query 1 top 50")
report_hub_hits(q2_top50_hub,   "Query 2 top 50")

# ---------------------------------------------------------------------------
# Step 6: MOA enrichment analysis
# ---------------------------------------------------------------------------
print("\n" + "=" * 70)
print("STEP 6: MOA enrichment analysis")
print("=" * 70)


def moa_enrichment(df_all, df_top50):
    def split_moa(series):
        counts = Counter()
        for val in series.dropna():
            for term in str(val).split("|"):
                term = term.strip()
                if term:
                    counts[term] += 1
        return counts

    all_counts  = split_moa(df_all["moa"])
    top_counts  = split_moa(df_top50["moa"])
    n_all   = len(df_all)
    n_top   = len(df_top50)

    rows = []
    for term, cnt_all in all_counts.items():
        cnt_top = top_counts.get(term, 0)
        # enrichment = fraction in top50 / fraction in all
        enrich = (cnt_top / n_top) / (cnt_all / n_all) if cnt_all > 0 else 0.0
        rows.append({
            "moa_term":        term,
            "count_top50":     cnt_top,
            "count_all":       cnt_all,
            "enrichment_ratio": round(enrich, 4),
        })

    return (pd.DataFrame(rows)
            .sort_values("enrichment_ratio", ascending=False)
            .reset_index(drop=True))


q1_moa = moa_enrichment(q1_named, q1_top50_neg)
q2_moa = moa_enrichment(q2_named, q2_top50_neg)

q1_moa.to_csv(os.path.join(OUT_DIR, "Query1_MOA_enrichment.csv"), index=False)
q2_moa.to_csv(os.path.join(OUT_DIR, "Query2_MOA_enrichment.csv"), index=False)

print("\n  Top 10 enriched MOA classes -- Query 1:")
for _, r in q1_moa[q1_moa["count_top50"] > 0].head(10).iterrows():
    print(f"    {r['moa_term']:40s}  top50={r['count_top50']:2d}  "
          f"all={r['count_all']:4d}  ratio={r['enrichment_ratio']:.3f}")

print("\n  Top 10 enriched MOA classes -- Query 2:")
for _, r in q2_moa[q2_moa["count_top50"] > 0].head(10).iterrows():
    print(f"    {r['moa_term']:40s}  top50={r['count_top50']:2d}  "
          f"all={r['count_all']:4d}  ratio={r['enrichment_ratio']:.3f}")

# ---------------------------------------------------------------------------
# Step 7: Final drug shortlist
# ---------------------------------------------------------------------------
print("\n" + "=" * 70)
print("STEP 7: Final drug shortlist")
print("=" * 70)

# Gather all drugs in overlap50 and overlap100 with hub annotation
overlap50_hub  = annotate_hubs(overlap50)
overlap100_full = overlap100_hub.copy()  # already annotated

overlap50_names  = set(overlap50_hub["pert_iname"].str.upper())
overlap100_names = set(overlap100_full["pert_iname"].str.upper())

def assign_tier(name_upper, hub_target):
    in50  = name_upper in overlap50_names
    in100 = name_upper in overlap100_names
    if in50 and hub_target:
        return 1
    if in50:
        return 2
    if in100 and hub_target:
        return 3
    return None  # not eligible

# Build shortlist from overlap50 union (overlap100 & hub)
shortlist_names = overlap50_names | {
    r["pert_iname"].upper()
    for _, r in overlap100_full.iterrows()
    if r["hub_target"]
}

rows = []
for _, r in overlap100_full.iterrows():
    uname = r["pert_iname"].upper()
    if uname not in shortlist_names:
        continue
    tier = assign_tier(uname, r["hub_target"])
    if tier is None:
        continue

    # Pull Q1/Q2 scores (use best row from full named list)
    q1_row = q1_named[q1_named["pert_iname"].str.upper() == uname]
    q2_row = q2_named[q2_named["pert_iname"].str.upper() == uname]
    q1_cs  = q1_row["raw_cs"].min() if len(q1_row) else np.nan
    q2_cs  = q2_row["raw_cs"].min() if len(q2_row) else np.nan

    # Rank within top 100
    q1_rank_series = q1_top100["pert_iname"].str.upper().reset_index(drop=True)
    q2_rank_series = q2_top100["pert_iname"].str.upper().reset_index(drop=True)
    q1_rank_idx = q1_rank_series[q1_rank_series == uname].index
    q2_rank_idx = q2_rank_series[q2_rank_series == uname].index
    q1_rank = int(q1_rank_idx[0]) + 1 if len(q1_rank_idx) else np.nan
    q2_rank = int(q2_rank_idx[0]) + 1 if len(q2_rank_idx) else np.nan

    rows.append({
        "pert_iname":      r["pert_iname"],
        "Q1_raw_cs":       round(q1_cs, 4) if not np.isnan(q1_cs) else np.nan,
        "Q2_raw_cs":       round(q2_cs, 4) if not np.isnan(q2_cs) else np.nan,
        "mean_cs":         round(np.nanmean([q1_cs, q2_cs]), 4),
        "Q1_rank":         q1_rank,
        "Q2_rank":         q2_rank,
        "moa":             r["moa"],
        "target_name":     r["target_name"],
        "hub_target":      r["hub_target"],
        "matched_hubs":    r["matched_hubs"],
        "n_hubs_targeted": r["n_hubs_targeted"],
        "evidence_tier":   tier,
        "als_literature_flag": "CHECK",
    })

shortlist = (pd.DataFrame(rows)
             .drop_duplicates("pert_iname")
             .sort_values(["evidence_tier", "mean_cs"])
             .reset_index(drop=True))

shortlist.to_csv(os.path.join(OUT_DIR, "final_drug_shortlist.csv"), index=False)

tier_counts = shortlist["evidence_tier"].value_counts().sort_index()
print(f"\n  Final shortlist: {len(shortlist)} drugs")
for tier, cnt in tier_counts.items():
    label = {1: "Tier 1 (top50 both + hub target)",
             2: "Tier 2 (top50 both)",
             3: "Tier 3 (top100 both + hub target)"}.get(tier, f"Tier {tier}")
    print(f"    {label}: {cnt}")

print("\n  Shortlist preview:")
print(shortlist[["pert_iname", "Q1_raw_cs", "Q2_raw_cs",
                  "evidence_tier", "moa", "matched_hubs"]].to_string(index=False))

# ---------------------------------------------------------------------------
# Step 8: Bubble plot
# ---------------------------------------------------------------------------
print("\n" + "=" * 70)
print("STEP 8: Drug shortlist visualization")
print("=" * 70)

if len(shortlist) == 0:
    print("  No drugs in shortlist -- skipping plot.")
else:
    tier_colors = {1: "#d62728", 2: "#ff7f0e", 3: "#1f77b4"}
    tier_labels = {1: "Tier 1: top50 both + hub",
                   2: "Tier 2: top50 both",
                   3: "Tier 3: top100 both + hub"}

    fig, ax = plt.subplots(figsize=(10, 8))

    # Reference lines
    ax.axvline(-0.3, color="gray", linestyle=":", linewidth=1, alpha=0.6)
    ax.axhline(-0.3, color="gray", linestyle=":", linewidth=1, alpha=0.6)
    xlim_all = shortlist[["Q1_raw_cs", "Q2_raw_cs"]].values.flatten()
    xlim_all = xlim_all[~np.isnan(xlim_all)]
    margin = 0.05
    xmin, xmax = xlim_all.min() - margin, xlim_all.max() + margin
    diag = np.linspace(xmin, xmax, 100)
    ax.plot(diag, diag, color="lightgray", linewidth=1, linestyle="--", zorder=0)

    for tier in sorted(shortlist["evidence_tier"].unique()):
        sub = shortlist[shortlist["evidence_tier"] == tier]
        sizes = 80 + sub["n_hubs_targeted"] * 120  # min size 80
        ax.scatter(
            sub["Q1_raw_cs"], sub["Q2_raw_cs"],
            s=sizes, c=tier_colors[tier],
            alpha=0.85, edgecolors="white", linewidths=0.7,
            label=tier_labels[tier], zorder=3,
        )

    # Labels
    for _, r in shortlist.iterrows():
        if pd.notna(r["Q1_raw_cs"]) and pd.notna(r["Q2_raw_cs"]):
            ax.annotate(
                r["pert_iname"],
                xy=(r["Q1_raw_cs"], r["Q2_raw_cs"]),
                xytext=(4, 3), textcoords="offset points",
                fontsize=7, alpha=0.9,
            )

    ax.set_xlabel("Query 1 raw_cs (Inflammation Signature)", fontsize=11)
    ax.set_ylabel("Query 2 raw_cs (Multi-omic Signature)", fontsize=11)
    ax.set_title("ALS CMap Drug Candidates -- Cross-Query Shortlist\n"
                 "(bubble size = N hub genes targeted)", fontsize=12)
    ax.legend(fontsize=9, loc="upper left")
    ax.text(0.99, 0.01,
            "Grey dotted lines: raw_cs = -0.3 threshold",
            transform=ax.transAxes, fontsize=7,
            ha="right", va="bottom", color="gray")
    sns.despine(ax=ax)
    plt.tight_layout()
    out_plot = os.path.join(OUT_DIR, "final_drug_shortlist_plot.png")
    plt.savefig(out_plot, dpi=150)
    plt.close()
    print(f"  Saved: {out_plot}")

# ---------------------------------------------------------------------------
# Step 9: Summary report
# ---------------------------------------------------------------------------
print("\n" + "=" * 70)
print("STEP 9: Summary report")
print("=" * 70)

q1_moa_top10 = q1_moa[q1_moa["count_top50"] > 0].head(10)
q2_moa_top10 = q2_moa[q2_moa["count_top50"] > 0].head(10)

overlap100_hub_drugs = overlap100_hub[overlap100_hub["hub_target"] == True]

summary = f"""
CMap L1000 Analysis Summary -- ALS Drug Repurposing
=====================================================
Date: {pd.Timestamp.today().strftime('%Y-%m-%d')}

GCT FILES USED
--------------
  Query 1: {Q1_GCT}
  Query 2: {Q2_GCT}

COMPOUND COUNTS AFTER FILTERING
---------------------------------
  Filter step                   Query 1       Query 2
  Total rows in GCT:        {q1_stats['n_raw']:>10,}    {q2_stats['n_raw']:>10,}
  After trt_cp:             {q1_stats['n_cp']:>10,}    {q2_stats['n_cp']:>10,}
  After is_exemplar_sig==1: {q1_stats['n_ex']:>10,}    {q2_stats['n_ex']:>10,}
  After qc_pass==1:         {q1_stats['n_qc']:>10,}    {q2_stats['n_qc']:>10,}
  After tas >= 0.2:         {q1_stats['n_tas']:>10,}    {q2_stats['n_tas']:>10,}
  After moa known (named):  {q1_stats['n_named']:>10,}    {q2_stats['n_named']:>10,}

SCORE RANGES (all compounds, post-filter)
------------------------------------------
  Query 1 raw_cs:   {q1_all['raw_cs'].min():.4f} to {q1_all['raw_cs'].max():.4f}
  Query 1 norm_cs:  {q1_all['norm_cs'].min():.4f} to {q1_all['norm_cs'].max():.4f}
  Query 2 raw_cs:   {q2_all['raw_cs'].min():.4f} to {q2_all['raw_cs'].max():.4f}
  Query 2 norm_cs:  {q2_all['norm_cs'].min():.4f} to {q2_all['norm_cs'].max():.4f}

  Query 1: N raw_cs < -0.3: {q1_sn}   N raw_cs < -0.2: {q1_mn}   N raw_cs > +0.3: {q1_sp}
  Query 2: N raw_cs < -0.3: {q2_sn}   N raw_cs < -0.2: {q2_mn}   N raw_cs > +0.3: {q2_sp}

CROSS-QUERY OVERLAP (named drugs, most negative raw_cs)
--------------------------------------------------------
  Top 25 overlap:  {len(overlap25)} drugs
  Top 50 overlap:  {len(overlap50)} drugs
  Top 100 overlap: {len(overlap100)} drugs

  Top 50 overlap drugs targeting hub genes: {overlap50_hub['hub_target'].sum()}
  Top 100 overlap drugs targeting hub genes: {len(overlap100_hub_drugs)}

TOP 10 ENRICHED MOA CLASSES -- QUERY 1
----------------------------------------
{q1_moa_top10[['moa_term','count_top50','count_all','enrichment_ratio']].to_string(index=False)}

TOP 10 ENRICHED MOA CLASSES -- QUERY 2
----------------------------------------
{q2_moa_top10[['moa_term','count_top50','count_all','enrichment_ratio']].to_string(index=False)}

FINAL DRUG SHORTLIST ({len(shortlist)} drugs)
{'='*(28+len(str(len(shortlist))))}
{shortlist[['pert_iname','Q1_raw_cs','Q2_raw_cs','mean_cs','evidence_tier',
            'moa','matched_hubs','n_hubs_targeted']].to_string(index=False)}

TIER BREAKDOWN
--------------
{chr(10).join(f"  Tier {t}: {c}" for t, c in tier_counts.items())}

IMPORTANT NOTES
---------------
All shortlisted drugs require manual literature review on PubMed for
existing ALS, motor neuron disease, or neurodegeneration evidence before
final reporting. CMap profiles were generated in cancer cell lines --
connectivity scores indicate transcriptional reversal in those contexts
and require biological validation in motor neuron models.

Evidence tiers:
  Tier 1 = top 50 of BOTH queries AND targets a STRING hub gene (highest confidence)
  Tier 2 = top 50 of BOTH queries (strong cross-query signal)
  Tier 3 = top 100 of BOTH queries AND targets a STRING hub gene (hub-supported)

als_literature_flag = CHECK for all entries -- perform manual PubMed search.

Submit drugs for further validation:
  - iNeuron or iPSC-derived motor neuron transcriptomic validation
  - ALS mouse model literature check (SOD1, TDP-43, FUS)
  - Check FDA approval status and CNS penetrance
"""

with open(os.path.join(OUT_DIR, "cmap_analysis_summary.txt"), "w") as fh:
    fh.write(summary)

print(summary)
print(f"\nAll outputs written to: {OUT_DIR}/")
print("Pipeline complete.")
