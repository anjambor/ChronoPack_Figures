# Fig6_CrossContext_v2: dynamic-mode enrichment and regression-screen panels (v4 kinetics tables)

> **Folder renamed 2026-11-17:** `Fig6_CrossContext_v2/` -> `Fig6_CrossContext/` (user request). Below, `Fig6_CrossContext_v2` means this folder; the
> earlier `Fig6_CrossContext` it was compared with (and other superseded versions) were deleted on 2026-11-16 and are in git
> history (commit 00ed526, `old/`). Script and output names keep their version tags.

Same figures as `Fig6_CrossContext/`, with two changes:

- **Inputs:** the regression-screen panels use the v4 kinetics tables in `Fig6_Reproduction/reference/kinetics/`
  instead of the Run3 tables, with the renamed feature families.
- **New figure:** the per-family (family-only) CV R² figure.

All eight inputs are copied unchanged from `Fig6_Reproduction/reference/`. The five dot-plot tables are
byte-identical to v1's inputs.

## Outputs (`output/`)

- `fig6_ora_cov_pls_motifs_transposed_v2_pctl_fg4_q05.{pdf,png}` (+ `_points.csv`): the transposed dot plot of
  ORA terms, covariates and PLS motifs per dynamic mode × context. Pixel-identical to v1.
- `fig6d_drop_one_family_delta_r2.{pdf,png}` (+ `_values.csv`): the CV R² lost when each feature family is dropped
  from the elastic-net model.
- `fig6e_individual_partial_r2.{pdf,png}` (+ `_values.csv`): per-feature OLS partial R² (ring = FDR < 0.05).
  The values are identical to v1; only a family label changed.
- `fig6de_composite.{pdf,png}`: 6d + 6e on one 260 × 90 mm page.
- `fig6_full_model_cv_r2.{pdf,png}` (+ `_values.csv`): the elastic-net full-model CV R² per outcome (mean ± SD over
  10 CV repeats). Diamonds show the cell-type-only model. This is the publication version of
  `refs/figureA_full_model_cv_r2.png`.
- **New:** `fig6_family_only_cv_r2_gain.{pdf,png}`. A family × outcome heatmap in the 6d style: the CV R² of each
  family-only model minus the cell-type-only CV R², i.e. what a family adds over cell type when used alone.
- **New:** `fig6_family_only_cv_r2_raw.{pdf,png}`. The family-only CV R² itself, with "Cell type only" as the top
  row. Both family-only figures share `fig6_family_only_cv_r2_values.csv`, which holds the mean, SD, min and max
  over repeats, the cell-type-only value and the gain.
- `comparison_to_v1.txt`: the checks below.

The family-only models include cell type. The v4 screen partials cell type out of the outcome and the
features, fits the elastic net on the residuals, and adds the cell-type part back for prediction
(`Fig6_Reproduction/docs/methods.md`). So raw family-only R² is dominated by cell type on the timing outcomes
(t_peak: 0.28 for cell type alone, 0.28–0.30 for every family-only model). That is why the gain version is the
main figure.

The two figures answer different questions, read together:

- **Family-only gain:** the signal a family carries on its own.
- **6d:** the signal it carries that the other families do not.

## Feature families (v4)

The five families are disjoint. Their membership is in
`Fig6_Reproduction/data/processed/kinetics/feature_family_map.csv`.

| Family (label) | Run3 / v1 name | Features | Figures |
|---|---|---|---|
| `promoter_architecture` (Promoter architecture) | `promoter_architecture` | 71 promoter (PLS) motif/architecture | 6d, 6e, family-only |
| `enhancer_architecture` (Enhancer architecture) | `matched_enhancer_motifs` | 132: 66 pELS-matched + 66 dELS-matched motifs | 6d, family-only |
| `enhancer_burden` (Enhancer burden) | `distal_cre` (6d label "Enhancer CRE burden") | 4: n_pELS, n_dELS, total pELS/dELS bp | 6d, 6e, family-only |
| `rna_stability` (RNA stability) | `rna_stability` | 2: median log10 decay rate, is_ARE | 6d, 6e, family-only |
| `gene_length` (Gene length) | `gene_architecture` | 1: gene_length_kb | 6d, 6e, family-only |

- **Labels:** `plot_fig6de.py` checks 6d's row labels against the table's `family_label` column.
- **Dropped:** v4 has no `all_enhancer_motifs`, separate pELS/dELS-matched sets, or general-pELS motifs. The
  general-pELS motifs are left out of the model.

## Checks (`scripts/compare_to_v1.py`; result: ALL CHECKS PASSED, 2026-10-23)

- **Inputs:** every file matches the sha256 in `Fig6_Reproduction/reference/manifest.json` (copy:
  `input/v1/Fig6_Reproduction_reference_manifest.json`).
- **Dot plot:** PNG pixel-identical to v1; points table identical.
- **6e:**
  - The OLS table is byte-identical to v1's Run3 table, and the plotted values are identical.
  - The family column maps onto v1's through the renames.
  - The PNG differs from v1 only in the renamed label ("Enhancer CRE burden" → "Enhancer burden").
- **Full model:**
  - v4's full model equals Run3's full model with the general-pELS motifs dropped (to 1e-16).
  - Mean CV R² changes by −0.003 to +0.014 (largest: max induction rate 0.100 → 0.114).
  - The cell-type-only diamonds are unchanged.
- **6d:**
  - Drop-one ΔR² changes by at most 0.011 (median 0.001), because it is now measured against that smaller full model.
  - Largest change: promoter architecture on peak log2FC, 0.053 → 0.064.
- **Family-only:**
  - The four families that Run3 also has match it exactly.
  - Enhancer architecture's family-only model is new in v4.

## Run

From this directory (~1 min):

```bash
sbatch scripts/slurm/fig6_pipeline.sbatch
```

Logs go to `output/logs/`, and a successful run ends with `FIG6_V2_PIPELINE_OK`.

| Step | Script | What |
|---|---|---|
| 1 | `plot_fig6_dotplot.py --variant v2 --cutoffs q05` | ORA + covariate + PLS-motif dot plot |
| 2 | `plot_fig6de.py --panel d e de` | 6d, 6e and the 6de composite |
| 3 | `plot_fig6_full_model_r2.py` | full-model CV R² bar chart with cell-type-only markers |
| 4 | `plot_fig6_family_only_r2.py` | family-only CV R² heatmaps (gain over cell type; raw) |
| 5 | `compare_to_v1.py` | `output/comparison_to_v1.txt` |

**Changes from v1's scripts:**

- Input and output names (no `_run3` suffix).
- The 6d/6e family keys and labels.
- `draw_d` takes row and colour-bar labels, so the family-only heatmaps reuse it. 6d itself is drawn as before.
- `compare_to_original.py` (the Fig7 pixel comparison) is replaced by `compare_to_v1.py`.

## Notebook

`Fig6.ipynb` (chronopack kernel) is v1's notebook adapted to v2. It adds sections for the full-model and
family-only figures and runs `compare_to_v1.py` at the end.

The committed notebook was executed in place with the default settings on 2026-10-23. A `FIG_RERUN=1`
execution, written to a scratch copy, redrew every figure inside the chronopack kernel. All PNGs and tables
were identical to the pipeline outputs, and the checks passed.

- `sbatch scripts/slurm/execute_notebook.sbatch` executes it in place.
- `FIG_RERUN=1 sbatch scripts/slurm/execute_notebook.sbatch` also redraws every figure inside the kernel.

## Inputs (`input/`, from `Fig6_Reproduction/reference/`)

| File | Reference subfolder | Used by |
|---|---|---|
| `ORA_results_v2_newbg_newfg4.csv`, `curated_ordered_grouped_specific_terms.csv` | `ora/`, `curated/` | dot plot, ORA rows |
| `logistic_regression_screen_results_v3_b3_pctl_fg4.csv` | `logistic/` | dot plot, covariate rows |
| `motif_summary_S10_CRE-matched_with_counts.csv`, `tata_box_positional_S10_summary_with_counterfactual_q.csv` | `motifs/`, `tata/` | dot plot, PLS motif rows |
| `cv_summary_results.csv` (replaces `ENRegressionScreen.cv_summary_results.Run3.csv`) | `kinetics/` | 6d, full-model, family-only |
| `canonical9_ols_results.csv` (replaces `ENRegressionScreen.OLS.individual_feature_regression_results.Run3.csv`) | `kinetics/` | 6e |
| `cv_celltype_only_results.csv` | `kinetics/` | full-model diamonds, family-only gain and raw top row |

`refs/figureA_full_model_cv_r2.png` is the reference for the full-model figure, copied from v1.

## External dependencies

All inputs are inside this folder. Since 2026-11-16 (when `Fig6_CrossContext` and `Fig6_Reproduction` moved to
`../old/`), `compare_to_v1.py` reads only copies in `input/v1/`: v1's outputs it compares against (`output/`), v1's
Run3 tables (`input/`) and `Fig6_Reproduction_reference_manifest.json`. Rechecked 2026-11-16: ALL CHECKS PASSED,
report unchanged.

Software outside the folder: `~/miniconda3/envs/dynamic_mode_motifs` (Python: matplotlib 3.10.8, pandas, numpy,
Pillow).
