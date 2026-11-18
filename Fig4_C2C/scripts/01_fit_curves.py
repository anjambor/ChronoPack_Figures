#!/usr/bin/env python
"""ChronoPack NB-GAM curves for the Fig4 trajectory panels (r9 data).

- THP1 ligands (the top 10 of ChronoPack, LIANA, NicheNet and TraSig; 30 ligands): one
  ChronoPack quick fit, cp.tl.fit_gene_dynamics(input_dir=<r9 THP1 tradeSeq inputs>, genes=...,
  final_nknots=8, tag="full"). It uses the same cells, offsets and 8 knots as the r9 final fit.
  tag="full" reads the unfiltered count matrix, so ligands removed from the r9 tradeSeq inputs
  (low-count or technical/tag-bias genes) can be fitted too. Ligands that are in the stored r9
  final fit (sce.1trajectory.nknots8.rds) are checked against it.
- K562 NFKBIA (the response curve): the stored r9 K562 final fit.

Curves: cp.tl.predict_gene_dynamics, log1p of the fitted mean count (and 95% band) at the fit's
reference library size, on a 1-min grid over 0-165 min.

Writes output/chronopack/fitted_curves.csv and ligand_fit_summary.csv.
"""
from __future__ import annotations

import os

import numpy as np
import pandas as pd

import chronopack as cp

from fig4_v2_common import (
    CHRONOPACK_DIR,
    CURVES_CSV,
    FIT_SUMMARY_CSV,
    K562_FIT_DIR,
    N_POINTS,
    NKNOTS,
    RECEIVER,
    RESPONSE_GENE,
    SENDER,
    STORED_SCE,
    THP1_FIT_DIR,
    TIME_MAX,
    curve_ligands,
)


def gene_list(fit_dir, name: str) -> set[str]:
    return set(pd.read_csv(fit_dir / name)["gene"].astype(str))


def reference_libsize(fit_dir) -> float:
    libsize = pd.read_csv(fit_dir / "cell_libsize_full.csv")["libsize_full"].to_numpy(float)
    return float(np.exp(np.mean(np.log(libsize))))


def total_counts(fit_dir, genes: list[str]) -> pd.Series:
    chunks = [c[c.index.isin(genes)] for c in pd.read_csv(fit_dir / "counts_genes_x_cells.full.csv.gz", index_col=0, chunksize=3000)]
    counts = pd.concat(chunks)
    return pd.Series({g: (float(counts.loc[g].sum()) if g in counts.index else np.nan) for g in genes})


def main() -> None:
    CHRONOPACK_DIR.mkdir(parents=True, exist_ok=True)
    ligands = curve_ligands()
    print(f"{len(ligands)} ligands: {', '.join(ligands)}")
    workers = int(os.environ.get("SLURM_CPUS_PER_TASK", "2"))

    fit = cp.tl.fit_gene_dynamics(input_dir=THP1_FIT_DIR, genes=ligands, final_nknots=NKNOTS, tag="full", workers=workers, verbose=False)
    thp1 = cp.tl.predict_gene_dynamics(fit, ligands, n_points=N_POINTS, time_max=TIME_MAX)
    thp1.insert(0, "cell_type", SENDER)
    k562 = cp.tl.predict_gene_dynamics(K562_FIT_DIR, [RESPONSE_GENE], n_points=N_POINTS, time_max=TIME_MAX, sce_file=STORED_SCE)
    k562.insert(0, "cell_type", RECEIVER)
    curves = pd.concat([thp1, k562], ignore_index=True)
    curves["t"] = curves["t"].round(6)
    curves.to_csv(CURVES_CSV, index=False)

    # Quick fit vs the stored r9 final fit, for ligands that are in it.
    stored_sig = set(pd.read_csv(THP1_FIT_DIR / "sig_df.1trajectory.nknots4.csv", index_col=0).index.astype(str))
    in_stored = [g for g in ligands if g in stored_sig]
    stored = cp.tl.predict_gene_dynamics(THP1_FIT_DIR, in_stored, n_points=N_POINTS, time_max=TIME_MAX, sce_file=STORED_SCE)
    diff = {g: float(np.max(np.abs(thp1.loc[thp1.gene == g, "mu_log1p"].to_numpy() - stored.loc[stored.gene == g, "mu_log1p"].to_numpy())))
            for g in in_stored}

    tested = gene_list(THP1_FIT_DIR, "subset_genes.csv")
    asso = fit.association
    r9_asso = pd.read_csv(THP1_FIT_DIR / "asso_df.1trajectory.nknots4.csv", index_col=0)
    counts = total_counts(THP1_FIT_DIR, ligands)
    summary = []
    for g in ligands:
        c = thp1[thp1.gene == g]
        summary.append({
            "gene": g,
            "cell_type": SENDER,
            "thp1_total_counts": counts[g],
            "in_r9_tradeseq_inputs": g in tested,
            "in_r9_final_fit": g in stored_sig,
            "r9_screen_pvalue": r9_asso["pvalue_1"].get(g, np.nan),
            "quick_fit_pvalue": asso["pvalue_1"].get(g, np.nan) if asso is not None and "pvalue_1" in asso else np.nan,
            "max_abs_diff_vs_r9_final_fit": diff.get(g, np.nan),
            "fit_peak_min": float(c.loc[c["mu_log1p"].idxmax(), "t"]),
            "fit_max_mu_log1p": float(c["mu_log1p"].max()),
        })
    summary = pd.DataFrame(summary)
    summary.to_csv(FIT_SUMMARY_CSV, index=False)
    with pd.option_context("display.width", 200):
        print(summary.to_string(index=False))
    print(f"Reference library size: {SENDER} {reference_libsize(THP1_FIT_DIR):.0f}, {RECEIVER} {reference_libsize(K562_FIT_DIR):.0f} UMIs")
    worst = max(diff.values()) if diff else 0.0
    print(f"Quick fit vs r9 final fit ({len(diff)} ligands): max |diff mu_log1p| = {worst:.3g}")
    if worst > 1e-6:
        raise SystemExit("Quick fit does not reproduce the stored r9 final fit.")
    print(f"Wrote {CURVES_CSV}\nWrote {FIT_SUMMARY_CSV}\nFIT_OK")


if __name__ == "__main__":
    main()
