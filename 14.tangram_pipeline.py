#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Tangram mapping pipeline (paths sanitized, ready for GitHub).

This script maps single-cell annotations to spatial spots using Tangram.
All file paths are configurable via CLI flags or environment variables.
No absolute or personal paths are embedded in the code.

Env vars (optional, override defaults):
  - TANGRAM_ST_DIR  : directory containing spatial .h5ad files
  - TANGRAM_SC_FILE : path to single-cell reference .h5ad
  - TANGRAM_OUT_DIR : output directory

Example:
  python tangram_pipeline.py \
    --st-dir ./st_h5ad \
    --st-file sample1.h5ad \
    --sc-file ./scRef.h5ad \
    --out-dir ./outputs \
    --cluster-label Ann_Level2
"""

from __future__ import annotations

import argparse
import os
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import squidpy as sq  # noqa: F401 (import kept if downstream users extend)
import matplotlib as mpl  # noqa: F401
import matplotlib.pyplot as plt  # noqa: F401
import seaborn as sns  # noqa: F401
import skimage  # noqa: F401
import tangram as tg
from anndata import AnnData


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Tangram mapping pipeline (path-sanitized).")

    parser.add_argument(
        "--st-dir",
        type=str,
        default=os.getenv("TANGRAM_ST_DIR", "./st_h5ad"),
        help="Directory containing spatial .h5ad files (default: ./st_h5ad or $TANGRAM_ST_DIR).",
    )
    parser.add_argument(
        "--st-file",
        type=str,
        required=True,
        help="Spatial h5ad file name within --st-dir (e.g., sample1.h5ad).",
    )
    parser.add_argument(
        "--sc-file",
        type=str,
        default=os.getenv("TANGRAM_SC_FILE", "./scRef.h5ad"),
        help="Path to single-cell reference .h5ad (default: ./scRef.h5ad or $TANGRAM_SC_FILE).",
    )
    parser.add_argument(
        "--out-dir",
        type=str,
        default=os.getenv("TANGRAM_OUT_DIR", "./outputs"),
        help="Directory to write outputs (default: ./outputs or $TANGRAM_OUT_DIR).",
    )
    parser.add_argument(
        "--cluster-label",
        type=str,
        default="Ann_Level2",
        help="Obs column in single-cell data used as cluster labels (default: Ann_Level2).",
    )
    parser.add_argument(
        "--min-features",
        type=int,
        default=10,
        help="Filter: minimum detected features per spot/cell (default: 10).",
    )
    parser.add_argument(
        "--min-counts",
        type=int,
        default=300,
        help="Filter: minimum total counts per spot/cell (default: 300).",
    )
    parser.add_argument(
        "--top-n-markers",
        type=int,
        default=100,
        help="Number of top-ranked genes per cluster to aggregate as marker set (default: 100).",
    )
    return parser.parse_args()


def ensure_outdir(path: Path) -> None:
    """Create output directory if it does not exist."""
    path.mkdir(parents=True, exist_ok=True)


def safe_filter_by_qc(adata: AnnData, min_features: int, min_counts: int) -> AnnData:
    """
    Filter AnnData by QC metrics using common column names with safe fallbacks.

    Tries typical Scanpy names first, then falls back to RNA-seq style names
    if present. If none exist, skips filtering with a warning.
    """
    candidates_features = ["n_genes_by_counts", "nFeature_RNA"]
    candidates_counts = ["total_counts", "nCount_RNA"]

    feat_col = next((c for c in candidates_features if c in adata.obs), None)
    cnt_col = next((c for c in candidates_counts if c in adata.obs), None)

    if feat_col is None or cnt_col is None:
        print(
            f"[WARN] QC columns not found (features={feat_col}, counts={cnt_col}). "
            "Skipping QC filtering."
        )
        return adata

    mask = (adata.obs[feat_col] > min_features) & (adata.obs[cnt_col] > min_counts)
    print(f"[INFO] QC filter kept {mask.sum()}/{adata.n_obs} observations.")
    return adata[mask].copy()


def compute_marker_set(adata_sc: AnnData, groupby: str, top_n: int) -> list[str]:
    """
    Rank genes per cluster and return a unique marker set.

    Uses Scanpy's rank_genes_groups (non-raw by default to avoid hidden layers).
    """
    if groupby not in adata_sc.obs:
        raise KeyError(f"'{groupby}' not found in single-cell .obs columns.")

    sc.tl.rank_genes_groups(adata_sc, groupby=groupby, use_raw=False)
    names_df = pd.DataFrame(adata_sc.uns["rank_genes_groups"]["names"]).iloc[:top_n, :]
    markers = list(np.unique(names_df.melt().value.values))
    print(f"[INFO] Aggregated {len(markers)} unique marker genes from top {top_n} per cluster.")
    return markers


def run_tangram(
    adata_sc: AnnData,
    adata_st: AnnData,
    markers: list[str],
    cluster_label: str,
) -> tuple[AnnData, AnnData]:
    """
    Prepare data and run Tangram cluster-mode mapping + gene projection.

    Returns
    -------
    ad_map : AnnData
        Mapping object returned by Tangram (cell-to-spot mapping).
    ad_ge : AnnData
        Projected gene expression in spatial coordinates.
    """
    # Harmonize gene spaces and subset to marker set
    tg.pp_adatas(adata_sc, adata_st, genes=markers)

    # Map clusters instead of individual cells (robust and faster)
    ad_map = tg.map_cells_to_space(
        adata_sc,
        adata_st,
        mode="clusters",
        cluster_label=cluster_label,
    )

    # Project cluster annotations to spatial spots
    tg.project_cell_annotations(ad_map, adata_st, annotation=cluster_label)

    # Project genes to spatial space
    ad_ge = tg.project_genes(
        adata_map=ad_map,
        adata_sc=adata_sc,
        cluster_label=cluster_label,
    )
    return ad_map, ad_ge


def main() -> None:
    args = parse_args()

    # Resolve paths (all relative by default; no absolute personal paths here)
    st_dir = Path(args.st_dir)
    st_file = st_dir / args.st_file
    sc_file = Path(args.sc_file)
    out_dir = Path(args.out_dir)
    ensure_outdir(out_dir)

    sc.logging.print_header()

    # ---- Load data
    print(f"[INFO] Reading spatial data: {st_file}")
    adata_st = sc.read(st_file)

    print(f"[INFO] Reading single-cell reference: {sc_file}")
    adata_sc = sc.read(sc_file)

    # Optional: seed a placeholder 'cell_count' if downstream expects it
    if "cell_count" not in adata_st.obs:
        adata_st.obs["cell_count"] = 20

    # ---- QC filtering (safe, optional if columns exist)
    adata_st = safe_filter_by_qc(adata_st, args.min_features, args.min_counts)

    # ---- Inspect cluster label distribution (non-fatal side effect)
    if args.cluster_label in adata_sc.obs:
        print(adata_sc.obs[args.cluster_label].value_counts())
    else:
        raise KeyError(f"'{args.cluster_label}' not found in adata_sc.obs")

    # ---- Marker computation
    markers = compute_marker_set(
        adata_sc=adata_sc,
        groupby=args.cluster_label,
        top_n=args.top_n_markers,
    )

    # ---- Tangram
    ad_map, ad_ge = run_tangram(
        adata_sc=adata_sc,
        adata_st=adata_st,
        markers=markers,
        cluster_label=args.cluster_label,
    )

    # ---- Outputs (all within out_dir; no personal paths)
    # Names mimic your originals but kept generic and safe
    prefix = Path(args.st_file).stem  # use input file stem as prefix
    out_ad_ge = out_dir / f"{prefix}.ad_ge_clusters.h5ad"
    out_ad_st = out_dir / f"{prefix}.adata_st_clusters.h5ad"
    out_ad_map = out_dir / f"{prefix}.ad_map_clusters.h5ad"
    out_ct_pred = out_dir / f"{prefix}.tangram_ct_pred.csv"

    print(f"[INFO] Writing: {out_ad_ge}")
    ad_ge.write_h5ad(out_ad_ge)

    print(f"[INFO] Writing: {out_ad_st}")
    adata_st.write_h5ad(out_ad_st)

    print(f"[INFO] Writing: {out_ad_map}")
    ad_map.write_h5ad(out_ad_map)

    # tangram_ct_pred lives in .obsm after project_cell_annotations
    if "tangram_ct_pred" in adata_st.obsm:
        print(f"[INFO] Writing: {out_ct_pred}")
        adata_st.obsm["tangram_ct_pred"].to_csv(out_ct_pred)
    else:
        print("[WARN] 'tangram_ct_pred' not found in adata_st.obsm; skipping CSV export.")

    print("[DONE] Tangram pipeline completed successfully.")


if __name__ == "__main__":
    main()
