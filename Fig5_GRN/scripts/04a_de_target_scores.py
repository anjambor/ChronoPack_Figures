#!/usr/bin/env python
"""Conventional differential-expression target scores, for the ChronoPack "TF prior x DE" ablation.

For each context, uses exactly the cells, genes and counts of the r9 tradeSeq fit that ChronoPack's
target scores come from (input/chronopack/<context>/tradeseq_fit: counts_genes_x_cells.csv.gz,
cell_times.csv.gz, cell_libsize_full.csv). Counts are normalized to 10,000 per cell (full library size)
and log1p-transformed; each later time point is compared with the first (time tag 1) by a two-sided
Wilcoxon rank-sum test (scanpy rank_genes_groups). A gene's DE p-value is its smallest p-value over
the time points, Bonferroni-corrected for the number of comparisons (capped at 1); its DE statistic
(used, like tradeSeq's Wald statistic, only to order equal scores) is the largest |z|.

The p-values then go through ChronoPack's own target scoring (cp.tl.grn.score_targets: BUM local FDR),
so the ablation replaces only the temporal test. Writes output/tf_balanced_global/de_scores/<context>/
de_association.csv and target_scores.scores.csv.
"""
from __future__ import annotations

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc

import chronopack as cp

from fig5_paths import CHRONOPACK_INPUT_DIR, CONTEXTS, OUTPUT_DIR

OUT_DIR = OUTPUT_DIR / "tf_balanced_global" / "de_scores"
PROTEIN_CODING = CHRONOPACK_INPUT_DIR / "protein_coding_genes.csv"


def load_fit_inputs(fit_dir) -> ad.AnnData:
    counts = pd.read_csv(fit_dir / "counts_genes_x_cells.csv.gz", index_col=0)
    times = pd.read_csv(fit_dir / "cell_times.csv.gz", index_col="cell_id")["time"]
    libsize = pd.read_csv(fit_dir / "cell_libsize_full.csv", index_col="cell_id")["libsize_full"]
    cells = counts.columns.astype(str)
    adata = ad.AnnData(X=counts.T.to_numpy(dtype=np.float32), obs=pd.DataFrame(index=cells), var=pd.DataFrame(index=counts.index.astype(str)))
    adata.obs["time"] = pd.Categorical(times.reindex(cells).astype(int).astype(str), categories=[str(t) for t in sorted(times.unique())])
    lib = libsize.reindex(cells).to_numpy(dtype=float)
    adata.X = np.log1p(adata.X / lib[:, None] * 1e4)
    return adata


def main() -> None:
    for context in CONTEXTS:
        fit_dir = CHRONOPACK_INPUT_DIR / context / "tradeseq_fit"
        out = OUT_DIR / context
        out.mkdir(parents=True, exist_ok=True)
        adata = load_fit_inputs(fit_dir)
        groups = list(adata.obs["time"].cat.categories)
        ref, later = groups[0], groups[1:]
        sc.tl.rank_genes_groups(adata, groupby="time", groups=later, reference=ref, method="wilcoxon", n_genes=adata.n_vars)
        res = adata.uns["rank_genes_groups"]
        pvals, zs = [], []
        for g in later:
            names = pd.Index(res["names"][g])
            pvals.append(pd.Series(res["pvals"][g], index=names).reindex(adata.var_names))
            zs.append(pd.Series(res["scores"][g], index=names).reindex(adata.var_names).abs())
        p = pd.concat(pvals, axis=1).min(axis=1) * len(later)
        table = pd.DataFrame({"pvalue_1": p.clip(upper=1.0), "waldStat_1": pd.concat(zs, axis=1).max(axis=1)})
        table.index.name = "gene"
        assoc = out / "de_association.csv"
        table.to_csv(assoc)
        scores = cp.tl.grn.score_targets(assoc, valid_targets=PROTEIN_CODING)
        scores.write(out, prefix="target_scores")
        print(f"{context}: {adata.n_obs} cells, {adata.n_vars} genes, {len(later)} comparisons vs time {ref}; "
              f"DE p < 0.05: {(table.pvalue_1 < 0.05).sum()}; scored targets: {len(scores.scores)}")
        print(scores.summary.to_string(index=False))


if __name__ == "__main__":
    main()
