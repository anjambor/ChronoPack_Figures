# Fig3_DynamicModes_v2: downsampling timing agreement rebuilt with ChronoPack on r9 data

> **Folder renamed 2026-11-17:** `Fig3_DynamicModes_v2/` -> `Fig3_DynamicModes/` (user request). Below, `Fig3_DynamicModes_v2` means this folder; the
> earlier `Fig3_DynamicModes` it was compared with (and other superseded versions) were deleted on 2026-11-16 and are in git
> history (commit 00ed526, `old/`). Script and output names keep their version tags.

Fig3_DynamicModes' timing composite takes its mode assignments from ChronoPack r9, but its onset/peak
times, curves and ribbons come from an old synthetic downsampling run. That run was made with raw
tradeSeq on pre-r9 data (`tmp/old/synthetic_ds/TradeSeqFits.SingleCell.SyntheticDownSample.ipynb`).
v2 repeats that experiment with ChronoPack's own API on the r9 K562 + TNF-alpha cells, rebuilding the
timing composite. Since 2026-11-16 (user request) v2 also builds Figure 3's broad/fine dynamic-modes figure with v1's
step 01 copied unchanged (step 00), so v2 holds both Figure 3 figures and v1 has moved to `../old/Fig3_DynamicModes`.

Outputs:
- `output/dynamic_modes_v3_K562a_broad_fine_examples_v4.{pdf,png}` (+ `_fine_cluster_examples_v4.csv`,
  `_selected_fine_cluster_examples_v4.csv`): step 00, from the r9 gene-trajectory tables (`input/pc_df.new.csv`,
  `input/gene_trajectory_fit.csv`, `input/gene_trajectory_time_grid.csv`; symlinks). Step 04 requires its tables to
  equal v1's (copied with v1's figure and its check against the original Fig3 to `input/v1/modes/`); 2026-11-16:
  tables identical, PNG 0 pixels differ.
- `output/publication_timing_gene_composite_..._q005_no_ybreak_vector_v2.{pdf,png}` plus
  `_mode_summary.csv` (n and mean ± SE absolute offset per mode) and `_genes.csv`.

## Design (matches the original run, except as noted)

| Fit | Cells (r9 K562a) | Time points | Knots |
|---|---|---|---|
| reference `K562_trunc` | 5,674 | tags 1–10 = 0–90 min, 10-min sampling | 8 |
| downsampled `K562_ds` | 2,582 | tags 1, 4, 7, 10 = 0/30/60/90 min | 4 |

- **Genes:** the 2,022 K562a genes of r9 `pc_df.new.csv`. The original run fitted 1,676 pre-r9 genes.
- **Pseudotime and size factors:** recomputed within each subset by the ChronoPack backend, as in
  the original run: (tag − 1) / 9, and full library size / mean.
- **Curves:** 91 points over 0–90 min at each fit's geometric-mean library size, with 95% bands.
- **Onset/peak:** `cp.tl.compute_gene_trajectory_descriptors` with package defaults (horizontal
  baseline over the first 10 min, z-difference 0.2, 30 consecutive minutes). Times are reported in
  **true minutes**; v1 reported the grid index, i.e. minute + 1, which cancels in the offsets.
- **Figure:** the same script as Fig3. It uses the same r9 modes and filters (|peak log2FC| > 0.5,
  p < 0.05, q < 0.05), the same 8 panel genes, and the same layout and styling.

## Run

From this directory, on a compute node:

```bash
RUN=K562_trunc sbatch scripts/slurm/fit.sbatch   # reference fit, ~20 min (4 workers, 90G)
RUN=K562_ds    sbatch scripts/slurm/fit.sbatch   # downsampled fit, ~4 min
sbatch scripts/slurm/fig3_v2_pipeline.sbatch     # steps 02-04 (fits missing models itself), ~15 min
```

Logs go to `output/logs/` (`FIT_OK <run>`, `FIG3_V2_PIPELINE_OK`). `REFIT=1 ./run_fig3_v2.sh` forces
refits.

| Step | Script | What |
|---|---|---|
| 01 | `01_fit_downsampling.py --run <run>` | Subsets the r9 AnnData and runs `cp.tl.fit_gene_dynamics(sub, prepare_inputs=True, genes=..., final_nknots=k)`. The run's tradeSeq inputs and `sce.1trajectory.nknots{k}.rds` go to `output/chronopack/<run>/`. |
| 02 | `02_predict_and_characterize.py` | `cp.tl.predict_gene_dynamics` (curves and bands), `cp.tl.observed_gene_dynamics` (per-time-point means of the panel genes), `cp.tl.compute_gene_trajectory_descriptors`. Writes `output/chronopack/{curves,observed_means,trajectory_descriptors}.csv`. |
| 03 | `03_make_timing_composite.R` (tradeseq_ce) | Fig3's figure script, reading the step-02 tables instead of the old run's `.rds`/gamList files |
| 04 | `04_compare_to_v1.py` | `output/comparison_to_v1.txt` |

**Notebook:** `Fig3_v2.ipynb` (chronopack kernel) walks through the same steps with explanations and
inline results. It shows the existing outputs by default; `RERUN = True` or `FIG_RERUN=1` regenerates
them. To execute it on a compute node: `sbatch scripts/slurm/execute_notebook.sbatch`.

## Inputs (`input/`)

- **Symlinks to ChronoPack r9:**
  - `adata_prefiltered.1e4-variant.fixed.h5ad`: K562 + TNF-alpha cells, `layers['counts']`, time tags,
    `total_counts` (= r9 `cell_libsize_full.csv`).
  - `pc_df.new.csv`: modes, filters and the gene list.
  - `r9_full_fit/` (r9 `sc-K562-1e4_to_tradeSeq.fixed`): r9's full 0–110 min fit, used only for the
    sanity check in step 04.
- **`v1/`** (copied from v1's `output/`, now `../old/Fig3_DynamicModes/output/`): v1's mode summary, gene table, figure, and its
  ChronoPack descriptors of the old run, which are identical to `desc.rds`.

## Results (2026-10-05) vs v1

| Mode | Peak offset (min), v2 | Peak offset (min), v1 | Onset offset (min), v2 | Onset offset (min), v1 |
|---|---|---|---|---|
| Impulse up 1+2 | 12.93 ± 1.67 (n=129) | 10.83 ± 1.65 (n=100) | 5.74 ± 0.78 | 4.44 ± 0.68 |
| Trend up | 4.83 ± 1.34 (n=109) | 4.40 ± 1.34 (n=85) | 6.59 ± 1.00 | 8.08 ± 1.25 |
| Trend down 1+2 | 4.27 ± 0.98 (n=179) | 3.90 ± 1.19 (n=118) | 6.25 ± 0.75 | 5.43 ± 0.83 |
| Transient down | 17.06 ± 3.89 (n=17) | 17.29 ± 4.65 (n=14) | 6.12 ± 1.89 | 4.29 ± 1.55 |

- **v2 has 434 genes vs v1's 317.** All 434 pass the same r9 filters, but 117 were not among the
  1,676 genes fitted in the old run, so v1 silently dropped them. All 317 v1 genes are in v2.
- **On the 317 shared genes, timing barely changes:**

  | Fit | Measure | Pearson r | Mean difference |
  |---|---|---|---|
  | Reference | onset | 0.97 | 0.35 min |
  | Reference | peak | 1.00 | 0.05 min |
  | Downsampled | onset | 0.95 | 0.8 min |
  | Downsampled | peak | 0.99 | 1.0 min |

  Restricted to these genes, v2's bars are close to v1's (e.g. impulse-up peak 11.4 vs 10.8 min,
  onset 4.7 vs 4.4). Most of the change comes from the 117 added genes; for example, the added
  impulse-up genes have a 18.1 min peak offset.
- **Sanity check:** v2's 0–90 min reference fit agrees with r9's full 0–110 min fit over 0–90 min.
  Across 434 genes, the per-gene curve correlation has median 0.995, and the max |Δ log1p| has
  median 0.005.

## External dependencies

All inputs are inside this folder or `../ChronoPack/`. Software outside both:
- `~/miniconda3/envs/chronopack` (chronopack, rpy2 and its R with tradeSeq): steps 01, 02 and 04;
- `~/miniconda3/envs/tradeseq_ce` (R ggplot2/patchwork): step 03, and nbconvert for the notebook.
