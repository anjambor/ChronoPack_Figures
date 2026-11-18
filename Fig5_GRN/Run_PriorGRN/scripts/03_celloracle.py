#!/usr/bin/env python
"""CellOracle GRN with the hg38 promoter base GRN (run inside resources/containers/celloracle_0.18.0.sif).

Follows the CellOracle GRN tutorial (0.18): import_anndata_as_raw_count (X = r9 cpm layer; CellOracle's log1p of it is
the r9 log1p layer), import_TF_data(base GRN), perform_PCA, knn_imputation with the tutorial's n_comps (explained-
variance elbow, at most 50) and k = 2.5% of cells, get_links(alpha=10, bagging_number=20). One GRN unit = all cells,
the same single network per dataset as the other methods (the tutorial fits one GRN per cluster).
EdgeWeight = coef_abs, CellOracle's own edge ranking (filter_links keeps the top edges by coef_abs after p < 0.001);
every base-GRN-supported edge that CellOracle fits is ranked, with the full Links table (coef_mean, p) kept.
Usage: python scripts/03_celloracle.py <dataset> <n_jobs>
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import anndata as ad
import celloracle as co
import numpy as np
import pandas as pd
import scanpy as sc

sys.path.insert(0, str(Path(__file__).resolve().parent))
import common as c  # noqa: E402

ds, n_jobs = sys.argv[1], int(sys.argv[2])
d = c.NETWORK_DIR / ds / "CELLORACLE"
d.mkdir(parents=True, exist_ok=True)
adata = ad.read_h5ad(c.INPUTS_DIR / ds / "expr.h5ad")
adata.obs["grn_unit"] = pd.Categorical(["all"] * adata.n_obs)
emb = sc.pp.pca(np.asarray(adata.layers["log1p"].todense()), n_comps=2, random_state=c.SEED)  # colours only
adata.obsm["X_pca2"] = emb

oracle = co.Oracle()
oracle.import_anndata_as_raw_count(adata=adata, cluster_column_name="grn_unit", embedding_name="X_pca2")
assert np.allclose(np.asarray(oracle.adata.layers["normalized_count"].todense()),
                   np.asarray(adata.layers["log1p"].todense()), atol=1e-5), "CellOracle log1p != r9 log1p"
oracle.import_TF_data(TF_info_matrix=pd.read_parquet(c.BASE_GRN))
oracle.perform_PCA()
n_comps = int(np.where(np.diff(np.diff(np.cumsum(oracle.pca.explained_variance_ratio_)) > 0.002))[0][0])
n_comps = min(n_comps, 50)
k = int(0.025 * oracle.adata.n_obs)
oracle.knn_imputation(n_pca_dims=n_comps, k=k, balanced=True, b_sight=k * 8, b_maxl=k * 4, n_jobs=n_jobs)
links = oracle.get_links(cluster_name_for_GRN_unit="grn_unit", alpha=10, bagging_number=20, verbose_level=1,
                         n_jobs=n_jobs)
L = links.links_dict["all"].copy()
L.to_csv(d / "links_all.csv", index=False)
e = L.rename(columns={"source": "Gene1", "target": "Gene2", "coef_abs": "EdgeWeight"})
c.write_ranked_edges(e, d / "rankedEdges.csv")
s = {"dataset": ds, "n_cells": int(oracle.adata.n_obs), "n_genes": int(oracle.adata.n_vars), "n_pca_dims": n_comps,
     "k": k, "edges": int(len(L)), "edges_p<0.001": int((L.p < 0.001).sum()), "tfs": int(L.source.nunique()),
     "targets": int(L.target.nunique()), "fig5_tf_edges": {t: int((L.source == t).sum()) for t in c.FIG5_TFS}}
(d / "summary.json").write_text(json.dumps(s, indent=2))
print(json.dumps(s))
print("CELLORACLE_OK")
