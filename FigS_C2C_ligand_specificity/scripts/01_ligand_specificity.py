#!/usr/bin/env python
"""ChronoPack ligand rank vs the best rank any comparator gives the same ligand (THP1 -> K562, Fig. 4 data).

Rebuild of the capsule figure c2c_downsample/benchmark_outputs/figures/chronopack_specific_ligand_recovery
(script copied to input/capsule/plot_chronopack_specific_ligands.py) on Fig. 4's r9 rankings:

- ChronoPack: Fig. 4's ligand ranking (chronopack_c2c.ligand_scores.csv, zscore_composite_rank; 49 ligands from
  cp.tl.c2c.infer on the r9 THP1 fit). Percentile = rank / 49.
- Comparators: the 4 runs drawn in Fig. 4 (common.FIG4_RUNS: LIANA and LIANA CellChat static_post_baseline,
  NicheNet lagged_mid_late, TraSig), read from Fig. 4's harmonized tables. Within each run a ligand's score is that
  of its best pair; ligands are ranked by it with ties at their mean rank, and percentile = ligand rank / number of
  ligands the run scores.
- Best comparator percentile = minimum over the runs; a ligand no run scores gets 1.05 (as in the capsule).
- Specificity = best comparator percentile - ChronoPack percentile.

Differences from the capsule definition: the capsule took the best rank over every comparator run (incl. LIANA's
5 scoring methods in 6 windows each and 4 NicheNet windows); comparator ranks are ligand-level (the capsule divided the ligand's best
*pair* rank by the largest such rank) and tied scores share their mean rank (the capsule used the harmonized
tables' arbitrary order within ties, e.g. CellChat's zero-score ties). The capsule's ChronoPack ranks (73 ligands,
joint 2-trajectory re-implementation) are replaced by Fig. 4's (49 ligands).

Writes output/tables/{comparator_ligand_ranks,ligand_specificity_fig4_runs}.csv and
output/figures/chronopack_specific_ligand_recovery_fig4_runs.{pdf,png,svg}.
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import common as c  # noqa: E402

import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

TEAL, GREY = "#218380", "#9aa0a6"
METHOD_LABEL = {"LIANA": "LIANA", "LIANA_CellChat": "CellChat", "LIANA_CellPhoneDB": "CellPhoneDB",
                "LIANA_Connectome": "Connectome", "LIANA_log2FC": "log2FC", "NicheNet": "NicheNet", "TraSig": "TraSig"}
LABEL_OFFSET = {"CCL3": ((0, 7), "center")}  # (offset points, ha) for labels that would collide


def chronopack_ranks() -> pd.DataFrame:
    lig = pd.read_csv(c.CHRONOPACK_LIGANDS)
    r = c.tie_ranks(lig, "composite_score", "ligand")
    assert (r["tie_size"] == 1).all(), "tied ChronoPack composite scores"
    assert (r["rank_min"] == r["zscore_composite_rank"]).all(), "composite_score order != zscore_composite_rank"
    out = r[["ligand", "best_pair", "zscore_composite_rank", "composite_score"]].rename(
        columns={"zscore_composite_rank": "chronopack_ligand_rank", "best_pair": "chronopack_best_pair"})
    out["chronopack_rank_pct"] = out["chronopack_ligand_rank"] / len(out)
    return out


def comparator_ligand_ranks() -> pd.DataFrame:
    rows = pd.concat([pd.read_csv(f) for f in c.COMPARATOR_TABLES], ignore_index=True, sort=False)
    rows = rows.loc[~rows["method"].astype(str).str.startswith("ChronoPack")]
    rows = rows.dropna(subset=["method", "run_id", "ligand", "rank", "score"])
    rows = rows.loc[pd.Series(list(zip(rows["method"], rows["run_id"])), index=rows.index).isin(c.FIG4_RUNS)]
    assert {tuple(x) for x in rows[["method", "run_id"]].drop_duplicates().to_numpy()} == c.FIG4_RUNS
    parts = []
    for (method, run), g in rows.groupby(["method", "run_id"], sort=True):
        g = g.sort_values("rank", kind="mergesort")
        # the harmonized rank must order scores (higher = better) so a ligand's best pair has its best score
        assert (np.diff(g["score"].to_numpy(float)) <= 1e-12).all(), f"{method}/{run}: rank does not follow score"
        best = (g.drop_duplicates("ligand", keep="first")[["ligand", "interaction_name", "score", "rank"]]
                .rename(columns={"interaction_name": "best_pair", "rank": "best_pair_rank"}))
        r = c.tie_ranks(best, "score", "ligand")
        r["method"], r["run_id"], r["n_ligands"] = method, run, len(r)
        r["rank_pct"] = r["rank_mean"] / len(r)
        parts.append(r)
    cols = ["method", "run_id", "ligand", "best_pair", "best_pair_rank", "score", "rank_min", "rank_mean",
            "rank_max", "tie_size", "n_ligands", "rank_pct"]
    return pd.concat(parts, ignore_index=True)[cols]


def specificity_table(cp: pd.DataFrame, comp: pd.DataFrame) -> pd.DataFrame:
    n_runs = comp[["method", "run_id"]].drop_duplicates().shape[0]
    best = (comp.sort_values(["ligand", "rank_pct", "method", "run_id"], kind="mergesort")
            .drop_duplicates("ligand", keep="first")
            .rename(columns={"method": "best_other_method", "run_id": "best_other_run_id",
                             "best_pair": "best_other_pair", "rank_mean": "best_other_ligand_rank",
                             "tie_size": "best_other_tie_size", "n_ligands": "best_other_n_ligands",
                             "rank_pct": "best_other_rank_pct"})
            [["ligand", "best_other_method", "best_other_run_id", "best_other_pair", "best_other_ligand_rank",
              "best_other_tie_size", "best_other_n_ligands", "best_other_rank_pct"]])
    n_scored = comp.groupby("ligand").size().rename("n_comparator_runs_scoring")
    t = cp.merge(best, on="ligand", how="left").merge(n_scored, on="ligand", how="left")
    t["n_comparator_runs_scoring"] = t["n_comparator_runs_scoring"].fillna(0).astype(int)
    t["n_comparator_runs"] = n_runs
    t["best_other_rank_pct_filled"] = t["best_other_rank_pct"].fillna(c.MISSING_PCT)
    t["best_other_label"] = np.where(
        t["best_other_method"].notna(),
        t["best_other_method"].map(METHOD_LABEL).fillna(t["best_other_method"]) + " / " + t["best_other_run_id"].astype(str),
        "not scored by any comparator run")
    t["chronopack_specificity"] = t["best_other_rank_pct_filled"] - t["chronopack_rank_pct"]
    return t.sort_values(["chronopack_specificity", "chronopack_ligand_rank"], ascending=[False, True]).reset_index(drop=True)


def plot(t: pd.DataFrame) -> None:
    top = t.head(c.N_LABEL)
    label = set(top["ligand"]) | {"TNF"}
    # the square scatter is height-limited at this size, so both panels have the same height
    fig, (ax, axb) = plt.subplots(1, 2, figsize=(10.2, 4.4), dpi=220,
                                  gridspec_kw={"width_ratios": [1.15, 1], "wspace": 0.42})

    d = t.sort_values("chronopack_ligand_rank")
    ax.scatter(d["chronopack_rank_pct"], d["best_other_rank_pct_filled"],
               c=np.where(d["chronopack_specificity"] > 0, TEAL, GREY),
               s=np.clip(78 * (d["chronopack_specificity"].clip(lower=0) + 0.12), 18, 82),
               alpha=0.88, edgecolor="white", linewidth=0.45, zorder=3)
    ax.plot([0, 1.05], [0, 1.05], color="0.58", lw=1.0, ls="--")
    if d["best_other_method"].isna().any():
        ax.axhline(c.MISSING_PCT, color="0.75", lw=0.8, ls=":")
        ax.text(1.06, c.MISSING_PCT - 0.012, "not scored by any comparator", fontsize=6.8, color="0.45", va="bottom")
    for r in d[d["ligand"].isin(label)].itertuples(index=False):
        dx, ha = LABEL_OFFSET.get(r.ligand, ((4, 3), "left"))
        ax.annotate(r.ligand, (r.chronopack_rank_pct, r.best_other_rank_pct_filled), xytext=dx, ha=ha,
                    textcoords="offset points", fontsize=7.2, color="#1f2933")
    ax.set_xlim(0, 1.07)
    ax.set_ylim(0, 1.09)
    ax.invert_xaxis()
    ax.invert_yaxis()
    ax.set_box_aspect(1)
    ax.set_xlabel("ChronoPack ligand rank percentile", fontsize=9)
    ax.set_ylabel("Best comparator ligand rank percentile", fontsize=9)
    ax.set_title("ChronoPack-specific ligand recovery", loc="left", fontsize=11)
    ax.grid(True, color="#e9ecef", lw=0.75)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)

    b = top.sort_values("chronopack_specificity")
    y = np.arange(len(b))
    axb.barh(y, b["chronopack_specificity"], color=TEAL, alpha=0.9)
    axb.set_yticks(y)
    axb.set_yticklabels([f"{l} (#{int(r)})" for l, r in zip(b["ligand"], b["chronopack_ligand_rank"])], fontsize=8)
    for yi, (v, method) in enumerate(zip(b["chronopack_specificity"], b["best_other_method"])):
        txt = METHOD_LABEL.get(method, method) if isinstance(method, str) else "none"
        axb.text(v + 0.01, yi, txt, va="center", fontsize=6.3, color="#4b5563")
    axb.set_xlim(0, max(1.0, b["chronopack_specificity"].max() + 0.2))
    axb.set_xlabel("Best comparator percentile - ChronoPack percentile", fontsize=9)
    axb.set_title("Largest ChronoPack advantage", loc="left", fontsize=11)
    axb.grid(True, axis="x", color="#e9ecef", lw=0.75)
    for s in ("top", "right", "left"):
        axb.spines[s].set_visible(False)
    axb.tick_params(axis="y", length=0)
    axb.set_axisbelow(True)

    # panel letters at the same height, left of each panel's tick labels
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    inv = fig.transFigure.inverted()
    boxes = [inv.transform(a.get_tightbbox(renderer)) for a in (ax, axb)]
    y = max(b[1, 1] for b in boxes)
    for b, letter in zip(boxes, ("(a)", "(b)")):
        fig.text(b[0, 0], y, letter, fontsize=12, fontweight="bold", va="bottom", ha="left")
    c.FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    for ext in ("png", "pdf", "svg"):
        f = c.FIGURE_DIR / f"chronopack_specific_ligand_recovery_fig4_runs.{ext}"
        fig.savefig(f, bbox_inches="tight")
        print("wrote", f)
    plt.close(fig)


def main() -> None:
    c.TABLE_DIR.mkdir(parents=True, exist_ok=True)
    cp = chronopack_ranks()
    comp = comparator_ligand_ranks()
    comp.to_csv(c.TABLE_DIR / "comparator_ligand_ranks.csv", index=False)
    print(f"ChronoPack: {len(cp)} ligands; comparators: {comp.groupby(['method', 'run_id']).ngroups} runs")
    print(comp.groupby(["method", "run_id"]).agg(n_ligands=("ligand", "size"), max_tie=("tie_size", "max"))
          .to_string())
    t = specificity_table(cp, comp)
    t.to_csv(c.TABLE_DIR / "ligand_specificity_fig4_runs.csv", index=False)
    with pd.option_context("display.width", 220):
        print(f"\n== top {c.N_LABEL} by specificity "
              f"({int((t['chronopack_specificity'] > 0).sum())}/{len(t)} ligands > 0; "
              f"{int(t['best_other_method'].isna().sum())} not scored by any comparator)")
        print(t.head(c.N_LABEL)[["ligand", "chronopack_ligand_rank", "chronopack_rank_pct", "best_other_label",
                                 "best_other_ligand_rank", "best_other_n_ligands", "best_other_rank_pct_filled",
                                 "chronopack_specificity"]].round(3).to_string(index=False))
        print(t.loc[t["ligand"] == "TNF", ["ligand", "chronopack_ligand_rank", "best_other_label",
                                            "best_other_ligand_rank", "best_other_rank_pct_filled",
                                            "chronopack_specificity"]].round(3).to_string(index=False))
    plot(t)
    print("LIGAND_SPECIFICITY_OK")


if __name__ == "__main__":
    main()
