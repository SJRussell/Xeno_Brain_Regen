#!/usr/bin/env python3
"""Aggregate STAR ReadsPerGene.out.tab files into a single gene x sample count matrix.

STAR ReadsPerGene format (one tab per sample):
    column 1: gene_id
    column 2: unstranded counts
    column 3: forward-stranded counts (--quantMode GeneCounts assumes -s 1)
    column 4: reverse-stranded counts (-s 2, our libtype)

First 4 rows are summary (N_unmapped, N_multimapping, N_noFeature, N_ambiguous).

Output:
    counts/star_gene_counts.tsv          - gene_id rows, sample columns, reverse-stranded
    counts/star_gene_summary.tsv         - per-sample summary metrics
"""
from __future__ import annotations
import sys
from pathlib import Path
import pandas as pd

STAR_DIR = Path("/project/def-srussell/nirosha_xeno/star")
COUNTS_DIR = Path("/project/def-srussell/nirosha_xeno/counts")
COUNTS_DIR.mkdir(parents=True, exist_ok=True)

SAMPLE_DIRS = sorted(p for p in STAR_DIR.iterdir() if p.is_dir() and (p / "ReadsPerGene.out.tab").exists())
print(f"Found {len(SAMPLE_DIRS)} STAR sample dirs")

# Read each ReadsPerGene.out.tab
counts_dict = {}
summary_rows = []
gene_ids_ref = None

for sd in SAMPLE_DIRS:
    sample = sd.name.replace("_R1_001", "").replace("5461_mix-", "")
    df = pd.read_csv(sd / "ReadsPerGene.out.tab", sep="\t", header=None,
                     names=["gene_id", "unstranded", "fwd", "rev"])
    # Top 4 rows are summary
    summary = df.iloc[:4].set_index("gene_id")
    body = df.iloc[4:].copy().reset_index(drop=True)

    if gene_ids_ref is None:
        gene_ids_ref = list(body["gene_id"])
    else:
        if list(body["gene_id"]) != gene_ids_ref:
            print(f"ERROR: gene order mismatch in {sample}", file=sys.stderr)
            sys.exit(1)

    # Reverse-stranded = our libtype
    counts_dict[sample] = body["rev"].values

    # Total counted (sum of body col rev)
    total_counted = int(body["rev"].sum())
    n_unmap = int(summary.loc["N_unmapped", "rev"])
    n_multi = int(summary.loc["N_multimapping", "rev"])
    n_nofeat = int(summary.loc["N_noFeature", "rev"])
    n_ambig = int(summary.loc["N_ambiguous", "rev"])
    n_input = n_unmap + n_multi + n_nofeat + n_ambig + total_counted

    summary_rows.append({
        "sample": sample,
        "input_reads": n_input,
        "counted_reads": total_counted,
        "unmapped": n_unmap,
        "multimapping": n_multi,
        "no_feature": n_nofeat,
        "ambiguous": n_ambig,
        "pct_counted": 100.0 * total_counted / n_input if n_input else 0.0,
    })

counts_df = pd.DataFrame(counts_dict, index=gene_ids_ref)
counts_df.index.name = "gene_id"
counts_df.to_csv(COUNTS_DIR / "star_gene_counts.tsv", sep="\t")
print(f"Counts matrix: {counts_df.shape[0]} genes × {counts_df.shape[1]} samples")
print(f"Wrote {COUNTS_DIR / 'star_gene_counts.tsv'}")

summary_df = pd.DataFrame(summary_rows).sort_values("sample").reset_index(drop=True)
summary_df.to_csv(COUNTS_DIR / "star_gene_summary.tsv", sep="\t", index=False)
print(f"Wrote {COUNTS_DIR / 'star_gene_summary.tsv'}")

# Quick QC: total counted reads per sample
print("\nPer-sample summary (counted reads, % counted):")
for _, row in summary_df.iterrows():
    print(f"  {row['sample']:<55s}  {int(row['counted_reads']):>10}  {row['pct_counted']:5.1f}%")

# Gene-level QC: how many genes have >0 reads in any sample
n_expressed = (counts_df > 0).any(axis=1).sum()
n_5cpm = (counts_df > 5).sum(axis=1)
n_genes_robust = (n_5cpm >= 3).sum()
print(f"\nGenes detected (>0 in any sample):     {n_expressed:,} / {counts_df.shape[0]:,}")
print(f"Genes robustly expressed (>5 in ≥3 sm): {n_genes_robust:,}")
