# FigS_C2C_ligand_specificity: ChronoPack vs best-comparator ligand ranks

This folder rebuilds the capsule figure `c2c_downsample/benchmark_outputs/figures/chronopack_specific_ligand_recovery`
(2026-09-28, from `c2c_downsample/scripts/plot_chronopack_specific_ligands.py`) on Fig. 4's r9 rankings. It covers the
same THP1 -> K562 data. The figure asks which ligands ChronoPack ranks well that every comparator ranks poorly.

Nothing in `Fig4_C2C/`, `Fig4_C2C/` or the `chronopack` package is changed.

**Output:** `output/figures/chronopack_specific_ligand_recovery_fig4_runs.{pdf,png,svg}`
- **(a):** a square scatter of ChronoPack's ligand rank percentile (x) against the best comparator percentile (y). Both
  axes are inverted, so the top right is good for both.
- **(b):** a bar chart of the 14 ligands with the largest ChronoPack advantage. Each bar is labelled with ChronoPack's
  rank (#) and the method that gives the best comparator rank.

**Figure legend (draft):**

> **Figure S#. Ligands prioritized by ChronoPack but not by existing cell–cell communication methods.**
> Ligand rankings for THP1-to-K562 signalling from ChronoPack and from the four comparator analyses in Fig. 4: LIANA
> rank aggregate and CellChat (both at the post-stimulation static window), NicheNet (lagged mid–late window) and
> TraSig. A comparator's ligand rank is the rank of the ligand's best-scoring ligand–receptor pair, with tied scores
> given their mean rank. Each rank is expressed as a percentile of the ligands that the method scores, so lower values
> are better.
> **(a)** ChronoPack ligand rank percentile against the best percentile that any comparator assigns to the same
> ligand, for all 49 ChronoPack ligands. Both axes are inverted, so ligands ranked highly by both approaches lie at
> the top right. The dashed line marks equal ranking. Teal points are ligands that ChronoPack ranks higher than every
> comparator, with point size proportional to the difference. The 14 ligands with the largest difference, and TNF,
> are labelled.
> **(b)** The 14 ligands with the largest ChronoPack advantage, defined as the best comparator percentile minus the
> ChronoPack percentile. Each label gives ChronoPack's rank in parentheses, and the text beside each bar names the
> comparator that gives the ligand its best rank.

## Definitions

| | |
|---|---|
| ChronoPack rank | Fig. 4's ligand ranking: `zscore_composite_rank` in `chronopack_c2c.ligand_scores.csv`, from `cp.tl.c2c.infer` on the r9 THP1 fit. It has 49 ligands. Percentile = rank / 49. |
| Comparator rank | Computed within each method/run. A ligand's score is the score of its best pair. Ligands are ranked by that score, and tied scores share their mean rank. Percentile = ligand rank / number of ligands that run scores. |
| Best comparator | The minimum percentile over the 4 comparator runs below. A ligand that no run scores gets 1.05, as in the capsule. In practice, every one of ChronoPack's 49 ligands is scored by at least one run (NicheNet or TraSig). |
| Specificity | Best comparator percentile minus ChronoPack percentile. A positive value means ChronoPack ranks the ligand higher. |

Comparator runs: the 4 runs drawn in Fig. 4, one per method, read from Fig4_C2C's harmonized tables (nothing is
rerun):
- LIANA (rank_aggregate) `static_post_baseline`;
- CellChat (LIANA's implementation) `static_post_baseline`;
- NicheNet `lagged_mid_late`;
- TraSig `official_experimental_time_pseudotime`.

The harmonized tables also hold LIANA's other windows and scoring methods, 32 runs in all. A version taking the best
rank over all 32 (the capsule's definition) gave the same top 14 except VEGFA; it was dropped on 2026-10-21.

### Differences from the capsule figure

- **ChronoPack:** the capsule ranked 73 ligands from its re-implementation on the joint 2-trajectory fit. This folder
  uses Fig. 4's 49 ligands from `cp.tl.c2c.infer` on the THP1-only fit.
  - 28 capsule ligands are not in the r9 set: GDF11, INHBA, CCL15, ANGPTL4, PTN, LTA, TNFSF12 and others, mostly genes
    with fewer than 10 THP1 counts.
  - 4 ligands are new in r9: CXCL12, EDN1, HGF and POMC.
- **Comparator runs:** only Fig. 4's 4 runs. The capsule took the best over every comparator run, including LIANA's
  5 scoring methods in 6 windows each and 4 NicheNet windows.
- **Ligand-level comparator ranks:** the capsule divided a ligand's best *pair* rank by the largest such rank.
- **Ties:** tied comparator scores share their mean rank. The capsule used the harmonized tables' arbitrary order
  within ties, for example CellChat's runs of zero-score ties of 17–21 ligands.
- **Layout:** square scatter, panels labelled (a) and (b), and no figure title, footnote or axis-direction notes. (The
  capsule's x-axis read "better left", but the axis is inverted, so better ranks were on the right.)

## Results (2026-10-21)

From the step 01 log and `output/comparison_to_capsule.txt`.

ChronoPack ranks 18 of its 49 ligands higher than any of the 4 comparator runs does.

| Ligand | ChronoPack rank (of 49) | Best comparator | Comparator rank | Specificity |
|---|---|---|---|---|
| HBEGF | 13 | TraSig | 60 / 69 | 0.60 |
| IL10 | 12 | NicheNet | 58 / 78 | 0.50 |
| BMP2 | 35 | NicheNet | 77 / 78 | 0.27 |
| VEGFA | 17 | TraSig | 41 / 69 | 0.25 |
| CCL2 | 11 | TraSig | 31 / 69 | 0.23 |
| SEMA3C | 14 | TraSig | 34 / 69 | 0.21 |
| CXCL2 | 5 | TraSig | 20 / 69 | 0.19 |
| CSF1 | 24 | TraSig | 46 / 69 | 0.18 |
| HGF | 21 | TraSig | 40 / 69 | 0.15 |
| TNF | 1 | NicheNet | 2 / 78 | 0.005 |

- **Capsule vs r9:** the capsule's top 14 and r9's share 7 ligands: HBEGF, IL10, CCL2, CSF1, CXCL2, CXCL1 and CXCL3.
  - The capsule's leader, GDF11 (0.72), and INHBA, CCL15 and ANGPTL4 are not r9 ChronoPack ligands.
  - On the 45 shared ligands, Spearman ρ is 0.997 for ChronoPack's percentile, 0.72 for the best comparator percentile
    and 0.84 for specificity.

**Caveat:** each method's percentile is relative to its own candidate set, which ranges from 33 to 78 ligands.
LIANA and CellChat do not score many of ChronoPack's ligands, so for those ligands the best comparator is NicheNet or
TraSig.

## Run

From this directory:

```bash
sbatch scripts/slurm/ligand_specificity.sbatch    # ~3 min, mostly imports
```

Logs go to `output/logs/`. A successful run ends with `LIGAND_SPECIFICITY_PIPELINE_OK`.

| Step | Script | What it does |
|---|---|---|
| 01 | `01_ligand_specificity.py` | Checks that ChronoPack's composite-score order equals `zscore_composite_rank` (no ties), and that each comparator run's harmonized rank follows its score. Writes `output/tables/comparator_ligand_ranks.csv` (per run and ligand), `ligand_specificity_fig4_runs.csv` and the figure. |
| 02 | `02_compare_to_capsule.py` | Writes `output/comparison_to_capsule.txt`: ligand sets, agreement on shared ligands, top-14 overlap, and the largest changes. |

`scripts/common.py` holds the paths, the comparator runs and the tie-aware ranking helper. The helper matches
`C2C_information_baseline/scripts/common.tie_ranks`.

## Inputs (`input/`)

- `fig4/`, copied from `../Fig4_C2C/output/` on 2026-10-21:
  - `chronopack/chronopack_c2c.{ligand_scores,interactions}.csv`;
  - `benchmark/results/harmonized_{liana_results,liana_individual_methods,nichenet_results,trasig_results}.csv`.

  To follow changes in Fig4_C2C, recopy these files.
- `capsule/`, copied from `c2c_downsample/`: the original figure's table and PNG
  (`benchmark_outputs/figures/chronopack_specific_ligand_recovery.{csv,png}`) and its script
  (`scripts/plot_chronopack_specific_ligands.py`). Step 02 reads the table; the PNG and script are kept for reference.

## External dependencies

`~/miniconda3/envs/chronopack` (pandas, matplotlib and scipy) runs both steps. You can override it with `CP_PYTHON`.
