#!/usr/bin/env python
"""Fit the reference or downsampled K562 + TNF-alpha trajectories with ChronoPack (r9 data).

Recreates the synthetic downsampling experiment of the original Fig3 (formerly raw tradeSeq in
TradeSeqFits.SingleCell.SyntheticDownSample.ipynb) with cp.tl.fit_gene_dynamics:

  K562_trunc  reference:   time tags 1-10  (0-90 min, 10-min sampling), 8 knots
  K562_ds     downsampled: time tags 1,4,7,10 (0/30/60/90 min),          4 knots

Each run fits the 2,022 K562a genes of the r9 gene-trajectory table (pc_df.new.csv) on all cells
at its time points. prepare_inputs=True writes that subset's tradeSeq inputs to
output/chronopack/<run>/; the backend recomputes pseudotime ((tag - 1) / 9) and size factors
(full library size / mean) within the subset, as the original run did. The fitted
SingleCellExperiment is copied next to the inputs for the prediction step.
"""
from __future__ import annotations

import argparse
import json
import os
import shutil
import time
from pathlib import Path

import anndata as ad
import pandas as pd

import chronopack as cp

FIG_DIR = Path(__file__).resolve().parents[1]
ADATA = FIG_DIR / "input" / "adata_prefiltered.1e4-variant.fixed.h5ad"
PC_DF = FIG_DIR / "input" / "pc_df.new.csv"
OUT = FIG_DIR / "output" / "chronopack"
RUNS = {
    "K562_trunc": {"time_tags": list(range(1, 11)), "nknots": 8},
    "K562_ds": {"time_tags": [1, 4, 7, 10], "nknots": 4},
}


def k562a_genes() -> list[str]:
    pc = pd.read_csv(PC_DF, index_col=0)
    return sorted(pc.loc[pc["cell_type"] == "K562a", "gene"].astype(str).unique())


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--run", choices=list(RUNS), required=True)
    parser.add_argument("--workers", type=int, default=int(os.environ.get("SLURM_CPUS_PER_TASK", "4")))
    args = parser.parse_args()
    spec = RUNS[args.run]
    run_dir = OUT / args.run
    run_dir.mkdir(parents=True, exist_ok=True)

    adata = ad.read_h5ad(ADATA, backed="r")
    tags = pd.to_numeric(adata.obs["time"])
    sub = adata[tags.isin(spec["time_tags"]).to_numpy()].to_memory()
    genes = k562a_genes()
    missing = sorted(set(genes) - set(sub.var_names))
    if missing:
        raise ValueError(f"{len(missing)} K562a genes are not in the AnnData: {missing[:10]}")
    print(f"{args.run}: {sub.n_obs} cells at time tags {spec['time_tags']}, {len(genes)} genes, "
          f"{spec['nknots']} knots, {args.workers} workers", flush=True)
    print(sub.obs["time"].value_counts().sort_index().to_string(), flush=True)

    start = time.time()
    fit = cp.tl.fit_gene_dynamics(
        sub,
        output_dir=run_dir,
        prepare_inputs=True,
        genes=genes,
        final_nknots=spec["nknots"],
        workers=args.workers,
    )
    sce_path = run_dir / f"sce.1trajectory.nknots{spec['nknots']}.rds"
    shutil.copy2(fit.sidecar_paths["sce"], sce_path)
    fit.association.to_csv(run_dir / f"asso_df.1trajectory.nknots{spec['nknots']}.csv")
    (run_dir / "fit_metadata.json").write_text(json.dumps({
        "run": args.run, "time_tags": spec["time_tags"], "nknots": spec["nknots"], "n_cells": int(sub.n_obs),
        "n_genes": len(fit.selected_genes), "workers": args.workers, "minutes": round((time.time() - start) / 60, 1),
        "adata": str(ADATA.resolve()),
    }, indent=2))
    print(f"Fitted {len(fit.selected_genes)} genes in {(time.time() - start) / 60:.1f} min -> {sce_path}")
    print(f"FIT_OK {args.run}")


if __name__ == "__main__":
    main()
