# Fig4_C2C_v2: ligand trajectories from ChronoPack NB-GAM fits

> **Folder renamed 2026-11-17:** `Fig4_C2C_v2/` -> `Fig4_C2C/` (user request). Below, `Fig4_C2C_v2` means this folder; the
> earlier `Fig4_C2C` it was compared with (and other superseded versions) were deleted on 2026-11-16 and are in git
> history (commit 00ed526, `old/`). Script and output names keep their version tags.

This is `../Fig4_C2C`'s THP1 -> K562 ligand prioritization figure with one change. In v1, the
right-hand trajectory panels show each method's top 10 ligands as per-time-point means of the r9
`log1p` layer. In v2 they show **ChronoPack NB-GAM (tradeSeq) fits** on the r9 data. Everything
else is v1's:
- the rankings, read from v1's tables (copied to `input/v1/`);
- panel A (TNF/TNFR rank);
- panel B (early log2FC, 45 vs 0 min, from the `log1p` layer);
- the manifest.

Step 03 checks that these are identical to v1's.

Output: `output/figures/fig4_v2_composite_temporal_prioritization.{pdf,png,svg}`, plus `_manifest.csv`
and `_positive_control_ranks.csv` (both identical to v1's).

Panel A tie blocks (2026-11-16, user request): a method's positive-control rank can sit in a block of rows tied on
its score, where the cited rank is only input order. The bar now shows the block's mid-rank and the label the block
("rank a-b (tie)"): LIANA 18-19 (TNF->TNFRSF1B tied with C3->IFITM1 on magnitude rank; bar 18.5), CellChat 23-50
(28 pairs at lr_probs = 0; bar 36.5). ChronoPack 1, NicheNet 2 and TraSig 3 are not tied. `_positive_control_ranks.csv`
keeps v1's `rank` column and adds `rank_min`, `rank_max`, `tie_n`, `mid_rank`; step 03 compares v1's columns.

## Curves

| Curve | Fit |
|---|---|
| THP1 ligands: 30, the top 10s of ChronoPack, LIANA, NicheNet and TraSig (CellChat is not drawn, as in v1) | One quick fit: `cp.tl.fit_gene_dynamics(input_dir=input/thp1_tradeseq, genes=..., final_nknots=8, tag="full")` |
| K562 NFKBIA | The stored r9 K562 final fit, `sce.1trajectory.nknots8.rds` |

- The quick fit uses the r9 THP1 fit's cells, offsets and 8 knots. For the 11 ligands that are in
  the stored r9 THP1 final fit, it reproduces that fit: max |Δ log1p| = 1e-14.
- `tag="full"` reads the unfiltered counts, so it also fits 6 ligands that are not in the r9 THP1
  tradeSeq inputs. They are marked † in the legends:
  - 4 low-count genes (< 10 counts in 787 cells): TNFSF12 (2), IL1A (3), TNFSF10 (3), CCL18 (9). Their
    fits are essentially zero, as are their v1 means.
  - APOE (4,810 counts) and TIMP1 (1,451), removed by `filter_technical_genes` (tag-bias / technical
    exclusion in `03_PrepareForTradeSeq`).
- Curves come from `cp.tl.predict_gene_dynamics`: log1p of the fitted mean count at each fit's
  reference library size (the geometric mean of library sizes: THP1 2,528, K562 3,306 UMIs), on a
  1-min grid over 0–165 min.
- **This scale differs from v1's.** It is lower than v1's normalized `log1p` layer, and the THP1
  and K562 reference sizes differ by ~30%. The panels are meant to compare timing, not levels.
- The 95% bands are in `output/chronopack/fitted_curves.csv` but not drawn: 10 bands per panel
  would be unreadable.

## Results (2026-10-05) vs v1

From `output/comparison_to_v1.txt`:
- **Dynamic ligands:** the fits track v1's means closely. Over the 12 time points, r = 0.96–0.99 for
  ChronoPack's top 10 (TNF, CXCL1/2/3/8, CCL3/3L3/4/5, IL1B) and for NFKBIA.
  - Fitted peaks are smoothed and sometimes later than v1's 15-min grid maxima. For example, NFKBIA
    peaks at 114 vs 105 min, and TNF at 85 vs 90 min.
  - None of ChronoPack's top 10 peaks inside 0–45 min, in either version.
- **Flat ligands** (most of LIANA's and TraSig's) have low r (median over all curves 0.62). For
  them, the fit removes time-point noise, so the curves become nearly flat lines and their peak
  minute is not meaningful.

## Run

From this directory, on a compute node:

```bash
sbatch scripts/slurm/fig4_v2_pipeline.sbatch            # ~2 min; reuses fitted curves if present
REFIT=1 sbatch scripts/slurm/fig4_v2_pipeline.sbatch    # refit too (~1 min more)
```

Logs go to `output/logs/`; a successful run ends with `FIG4_V2_PIPELINE_OK`. A `REFIT=1` run
reproduced `fitted_curves.csv` and the manifest byte-for-byte.

| Step | Script | What |
|---|---|---|
| 01 | `01_fit_curves.py` | Quick fit of the THP1 ligands, prediction of the THP1 curves and K562 NFKBIA, and the check against the stored r9 fit. Writes `output/chronopack/{fitted_curves,ligand_fit_summary}.csv` and fails if the check differs by > 1e-6. |
| 02 | `02_plot_fig4_v2.py` | v1's `07_plot_fig4.py`, with the trajectory panels reading step 01's curves |
| 03 | `03_compare_to_v1.py` | Writes `output/comparison_to_v1.txt` and `output/chronopack/fitted_vs_observed_curves.csv`. Fails if the manifest or ranks differ from v1's. |

Shared paths, constants and the top-10 selection (as in v1) are in `scripts/fig4_v2_common.py`.

**Notebook:** `Fig4_v2.ipynb` (chronopack kernel) walks through the steps with inline results. It
shows the existing outputs by default; `RERUN = True` or `FIG_RERUN=1` regenerates them. To
execute it on a compute node: `sbatch scripts/slurm/execute_notebook.sbatch`.

## Inputs (`input/`)

- **Symlinks to ChronoPack r9** (`.../full_pipeline_analyses_20260930_r9/data/thp1-k562-lps.hg38.rep/`):
  - `adata_full.chk1.h5ad` (`output/`): for panel B and the manifest, as in v1.
  - `thp1_tradeseq`, `k562_tradeseq` (`intmd/sc-{THP1,K562}_to_tradeSeq.fixed`): the r9 tradeSeq
    inputs and fits.
- **`v1/`** (copied from `../Fig4_C2C/output/`, 2026-10-05):
  - ChronoPack ligand scores and interactions;
  - the harmonized LIANA, CellChat, NicheNet and TraSig results;
  - v1's manifest, positive-control ranks and figure.

  To rebuild v2 after v1 changes, recopy these files.

## External dependencies

`~/miniconda3/envs/chronopack` (chronopack, and rpy2 with its R and tradeSeq) runs all steps.
`~/miniconda3/envs/tradeseq_ce` provides nbconvert for executing the notebook.
