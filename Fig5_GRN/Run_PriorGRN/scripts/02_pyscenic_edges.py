#!/usr/bin/env python
"""pySCENIC regulons -> rankedEdges.csv (run inside resources/containers/pyscenic_0.12.1.sif).

Edges are the TF -> target pairs of the activating regulons that `pyscenic ctx` keeps after cisTarget pruning
(df2regulons merges a TF's enriched motifs into one regulon); EdgeWeight is the GRNBoost2 importance that ctx carries
into the regulon. Pairs that cisTarget prunes are not predictions of pySCENIC and are left out.
Usage: python scripts/02_pyscenic_edges.py <dataset>
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import pandas as pd
from pyscenic.transform import df2regulons
from pyscenic.utils import load_motifs

sys.path.insert(0, str(Path(__file__).resolve().parent))
import common as c  # noqa: E402

ds = sys.argv[1]
d = c.NETWORK_DIR / ds / "PYSCENIC"
adj = pd.read_csv(d / "adjacencies.tsv", sep="\t")
regulons = df2regulons(load_motifs(d / "regulons.csv"))
rows = [(r.transcription_factor, g, w) for r in regulons for g, w in r.gene2weight.items()]
e = pd.DataFrame(rows, columns=["Gene1", "Gene2", "EdgeWeight"])
e = e[e.Gene1 != e.Gene2].groupby(["Gene1", "Gene2"], as_index=False).EdgeWeight.max()
# regulon weights are the GRNBoost2 importances of the same pairs
chk = e.merge(adj.rename(columns={"TF": "Gene1", "target": "Gene2"}), on=["Gene1", "Gene2"], how="left")
assert chk.importance.notna().all() and (chk.importance - chk.EdgeWeight).abs().max() < 1e-9, "weights != importances"
c.write_ranked_edges(e, d / "rankedEdges.csv")
s = {"dataset": ds, "grnboost2_pairs": int(len(adj)), "grnboost2_tfs": int(adj.TF.nunique()),
     "regulons": len(regulons), "regulon_tfs": sorted({r.transcription_factor for r in regulons}),
     "edges": int(len(e)), "targets": int(e.Gene2.nunique()),
     "fig5_tf_edges": {t: int((e.Gene1 == t).sum()) for t in c.FIG5_TFS}}
(d / "summary.json").write_text(json.dumps(s, indent=2))
print(json.dumps({k: v for k, v in s.items() if k != "regulon_tfs"}))
print("PYSCENIC_EDGES_OK")
