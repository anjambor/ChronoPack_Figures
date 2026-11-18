"""Where Run_BEELINE gets Fig5's scoring code (steps 05, 06, 10).

Scoring code: Fig5_GRN/scripts (Run_BEELINE lives inside Fig5_GRN since 2026-11-16), imported read-only with FIG5_METHOD_SET=expression and the full Fig5 universe.
Its benchmark_grn_systems.evaluate_one has the EPR top-k fix (2026-10-22): TP_at_k counts the hits at ranks 1..k.

History: until 2026-11-16 these steps imported Fig5_GRN_v3 (now ../old/Fig5_GRN_v3; same scoring code) and also scored
the metacell BEELINE run of Fig5 v1 (../old/Fig5_GRN), checked against v1's tables. The metacell round and those
checks were dropped (user request, 2026-11-16) so that Run_BEELINE depends on no folder in old/. v4 differs from v3
only in ChronoPack's TF rankings (ChronoPack_THP1 no longer scores CREB1).

Importing this module puts v4's script directories first on sys.path.
"""
from __future__ import annotations

import os
import sys
from pathlib import Path

V4 = Path(__file__).resolve().parents[2]  # Fig5_GRN (this folder is Fig5_GRN/Run_BEELINE)
os.environ["FIG5_METHOD_SET"] = "expression"
os.environ.pop("FIG5_UNIVERSE", None)  # Fig5's full universe (19,338 protein-coding targets)
for p in (V4 / "scripts/capsule_lib/02_visualize_publication_figures",
          V4 / "scripts/capsule_lib/01_generate_predictions_and_metrics", V4 / "scripts"):
    if str(p) not in sys.path:
        sys.path.insert(0, str(p))

V4_PER_TF = (V4 / "output" / "benchmark_outputs" / "tables" / "publication_figures"
             / "system_broad_filtered_publication_v2_per_tf_metrics_all_methods.csv")
