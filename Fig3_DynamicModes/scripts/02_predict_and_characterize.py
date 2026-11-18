#!/usr/bin/env python
"""Curves, observed means and onset/peak times for the two ChronoPack downsampling fits.

For each run (K562_trunc, K562_ds) fitted by 01_fit_downsampling.py:
- cp.tl.predict_gene_dynamics: fitted mean and 95% band for every gene at minutes 0-90
  (91 points; the run's own cell_libsize_full.csv gives the reference offset, as in the fit);
- cp.tl.observed_gene_dynamics: normalized counts of the panel genes, averaged per time point
  (count scale, then log1p), as in the original figure.
Both runs' curves are then characterized together with cp.tl.compute_gene_trajectory_descriptors
(package defaults; var['time'] = (1..91)/91, as in Fig3 step 02). The descriptor times are the
grid index (1..91); 1 is subtracted so t_onset / t_peak are true minutes (0-90).

Writes output/chronopack/{curves.csv, observed_means.csv, trajectory_descriptors.csv}.
"""
from __future__ import annotations

from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

import chronopack as cp

FIG_DIR = Path(__file__).resolve().parents[1]
OUT = FIG_DIR / "output" / "chronopack"
RUNS = {"K562_trunc": 8, "K562_ds": 4}  # run -> knots
PANEL_GENES = ["NFKBIA", "CXCL1", "ATF3", "TNFAIP2", "IRF1", "HES1", "BIRC3", "CD69"]
TIME_MAX, N_POINTS = 90.0, 91


def main() -> None:
    curves, means = [], []
    for run, nknots in RUNS.items():
        run_dir = OUT / run
        genes = pd.read_csv(run_dir / f"asso_df.1trajectory.nknots{nknots}.csv", index_col=0).index.astype(str).tolist()
        c = cp.tl.predict_gene_dynamics(run_dir, genes, n_points=N_POINTS, time_max=TIME_MAX,
                                        sce_file=f"sce.1trajectory.nknots{nknots}.rds")
        curves.append(c.assign(run=run, t=c["t"].round(6)))
        obs = cp.tl.observed_gene_dynamics(run_dir, PANEL_GENES, time_max=TIME_MAX)
        m = (obs.assign(count=np.expm1(obs["expr_log1p"]), t=obs["t"].round(6))
             .groupby(["gene", "t"], as_index=False)["count"].mean()
             .assign(mean_log1p=lambda d: np.log1p(d["count"]), run=run).drop(columns="count"))
        means.append(m)
        print(f"{run}: {len(genes)} genes predicted; observed means at t = {sorted(m['t'].unique())}")
    curves = pd.concat(curves, ignore_index=True)
    means = pd.concat(means, ignore_index=True)
    curves.to_csv(OUT / "curves.csv", index=False)
    means.to_csv(OUT / "observed_means.csv", index=False)

    wide = curves.assign(row=curves["gene"] + "_" + curves["run"]).pivot(index="row", columns="t", values="mu_log1p")
    if wide.shape[1] != N_POINTS or wide.isna().any().any():
        raise ValueError(f"Expected a complete {N_POINTS}-point grid, got {wide.shape} with NaNs={wide.isna().sum().sum()}")
    traj = ad.AnnData(
        X=wide.to_numpy(float),
        obs=pd.DataFrame(index=wide.index.rename("gene_celltype")),
        var=pd.DataFrame({"time": np.arange(1, N_POINTS + 1) / N_POINTS}, index=[f"t{int(t)}" for t in wide.columns]),
    )
    cp.tl.compute_gene_trajectory_descriptors(traj)
    desc = traj.obs.copy()
    desc.insert(0, "cell_type", [i.rsplit("_", 2)[-2] + "_" + i.rsplit("_", 1)[-1] for i in desc.index])
    desc.insert(0, "gene", [i[: -len(ct) - 1] for i, ct in zip(desc.index, desc["cell_type"])])
    desc["t_onset"] = desc["t_onset"] - 1  # grid index (1..91) -> minutes (0..90)
    desc["t_peak"] = desc["t_peak"] - 1
    desc.to_csv(OUT / "trajectory_descriptors.csv")
    print(desc.groupby("cell_type")[["t_onset", "t_peak"]].describe().round(1).T.to_string())
    print("PREDICT_OK")


if __name__ == "__main__":
    main()
