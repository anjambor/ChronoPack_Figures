# ChronoPack_figures

Scripts, inputs and outputs for the ChronoPack manuscript figures. Each `Fig*/` folder has
`input/`, `scripts/` and `output/`, a `run_figN.sh` driver, a Slurm wrapper in
`scripts/slurm/`, and a README describing provenance, differences from the original figures,
and external dependencies.

| Folder | Figure | Data |
|---|---|---|
| `Fig3_DynamicModes/` | both Fig3 figures: K562a broad/fine dynamic modes (v1's step, unchanged) and the downsampling timing agreement rebuilt with ChronoPack fits (reference 0–90 min vs 30-min downsampled) | ChronoPack r9 K562 + TNF-alpha cells and gene-trajectory table |
| `Fig4_C2C/` | Fig4 (THP1 -> K562 ligand prioritization vs LIANA / CellChat / NicheNet / TraSig) with the ligand trajectory panels drawn from ChronoPack NB-GAM fits | ChronoPack r9 THP1 / K562 tradeSeq inputs + v1's rankings (copied) |
| `FigS_C2C_ligand_specificity/` (renamed from `Fig4_C2C_ligand_specificity/` 2026-11-16) | ChronoPack vs best-comparator ligand rank (rebuild of the capsule's `chronopack_specific_ligand_recovery`) | Fig4_C2C's r9 ChronoPack ligand ranks + comparator tables (copied) |
| `Fig5_GRN/` | TF -> target GRN benchmark (panels 5b-5e) with the BEELINE comparators from the r9 single-cell P < 0.05 run, in two method sets (expression-only and with the motif-prior methods pySCENIC, CellOracle, Inferelator), ChronoPack's TF rankings from the current ChronoPack (`rank_regulators` defaults since commit e0c89fb, pinned OmniPath snapshot); also the with_prior shared-universe (`*_v2`) and Perturb-seq reference figures | ChronoPack r9 GRN inference + `Run_BEELINE` / `Run_PriorGRN` predictions (symlinked) + Replogle / Yao predictions (copied) |
| `Fig5_GRN/Run_BEELINE/` | not a figure (moved into Fig5_GRN 2026-11-16, whose BEELINE comparators it produces): BEELINE comparators (GENIE3 / GRNBOOST2 / PPCOR / SINCERITIES) rerun on the r9 single cells (K562 and THP1), TFs + 1000 / TFs + 2000, three gene selections | ChronoPack r9 K562 + TNF-alpha and THP1 + LPS cells |
| `Fig5_GRN/Run_PriorGRN/` | not a figure (moved into Fig5_GRN 2026-11-16): motif-prior GRN comparators (pySCENIC, CellOracle, Inferelator) on the Run_BEELINE P < 0.05 datasets | same cells and genes as Run_BEELINE |
| `Fig6_CrossContext/` | ORA / covariate / motif dot plot; regression-screen panels 6d/6e from the v4 kinetics tables (renamed families), plus the family-only CV R² figure | the eight `Fig6_Reproduction/reference/` tables of the user-supplied reproducibility package (copied) |
| `FigS_ClusteringMetrics/` | gene-trajectory clustering metrics (silhouette, bootstrap ARI / Jaccard) from chrono_seq TESTING.ipynb, recomputed on r9 | ChronoPack r9 gene clustering (symlink) + TESTING.ipynb's saved results (copied) |

The current folders were renamed on 2026-11-17 (Fig3_DynamicModes_v2 -> Fig3_DynamicModes, Fig4_C2C_v2 -> Fig4_C2C,
Fig5_GRN_v4 -> Fig5_GRN, Fig6_CrossContext_v2 -> Fig6_CrossContext); their output and script names keep their
version tags (e.g. `fig4_v2_*`, `run_fig5_v4.sh`). The superseded versions (the earlier Fig3_DynamicModes, Fig4_C2C,
Fig5_GRN v1-v3, Fig6_CrossContext, plus Fig6_Reproduction and C2C_information_baseline) were deleted on 2026-11-16;
their tracked files are in git history (commit 00ed526, under `old/`). No current folder depends on them.

Run each figure from its folder on the cluster, e.g. `cd Fig4_C2C && sbatch scripts/slurm/fig4_v2_pipeline.sbatch`,
or open its notebook (e.g. `Fig3_DynamicModes/Fig3_v2.ipynb`, `Fig6_CrossContext/Fig6.ipynb`, chronopack kernel), which explains and shows each step and
by default displays the existing outputs (`RERUN = True` regenerates them).

## Not in git

- Symlinks into `../ChronoPack/` (r9 data) are committed as links; they resolve only on the cluster
  with ChronoPack checked out next to this folder.
- Files over 50 MB (copied adata, synthetic-run tradeSeq fits, BEELINE predictions, NicheNet matrix,
  large output tables) are listed in `.gitignore` and stay on the cluster under
  `/new-stg/home/alex/tmp/ChronoPack_figures/`.
- Logs, matplotlib caches and Fig4's `output/verify_run_mode/` are not tracked.
