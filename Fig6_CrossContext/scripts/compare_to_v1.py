#!/usr/bin/env python3
"""Check Fig6_CrossContext_v2 against its inputs' source and against Fig6_CrossContext (v1, Run3 tables).

1. Every file in input/ must match Fig6_Reproduction/reference/manifest.json (sha256).
2. The dot plot uses the same five tables as v1, so its PNG must be pixel-identical and its points table identical.
3. 6e: the OLS table is byte-identical to v1's Run3 table, so the plotted values must be identical; only the family
   names changed (distal_cre -> enhancer_burden, gene_architecture -> gene_length). The PNG differs only in the
   renamed 6e family label ("Enhancer CRE burden" -> "Enhancer burden"), reported as a pixel count.
4. 6d and the full-model figure: v4 values against Run3, matched through the family renames (reported).
5. Family-only models: the four families that Run3 also has must match it; enhancer architecture is new.
Writes output/comparison_to_v1.txt; exit code 1 if a check fails.
"""
from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from PIL import Image

ROOT = Path(__file__).resolve().parents[1]
OUT, INPUT = ROOT / "output", ROOT / "input"
# v1's outputs / inputs and Fig6_Reproduction's manifest, copied into input/v1/ on 2026-11-16 when Fig6_CrossContext
# and Fig6_Reproduction moved to ../old/ (no dependencies on old/)
V1 = INPUT / "v1"
V1_OUT, V1_INPUT = V1 / "output", V1 / "input"
MANIFEST = V1 / "Fig6_Reproduction_reference_manifest.json"
# v4 family -> Run3 family (Fig6_Reproduction/data/processed/kinetics/feature_family_map.csv, source_family).
RUN3_FAMILY = {
    "promoter_architecture": "promoter_architecture",
    "enhancer_architecture": "matched_enhancer_motifs",
    "enhancer_burden": "distal_cre",
    "rna_stability": "rna_stability",
    "gene_length": "gene_architecture",
}
lines: list[str] = []
failures: list[str] = []


def say(text: str = "") -> None:
    print(text)
    lines.append(text)


def check(ok: bool, label: str) -> None:
    say(f"  [{'OK' if ok else 'FAIL'}] {label}")
    if not ok:
        failures.append(label)


def png_diff(a: Path, b: Path) -> int | str:
    x, y = np.asarray(Image.open(a).convert("RGB")), np.asarray(Image.open(b).convert("RGB"))
    return f"size {x.shape} vs {y.shape}" if x.shape != y.shape else int(np.any(x != y, axis=2).sum())


def same_table(x: pd.DataFrame, y: pd.DataFrame) -> bool:
    if x.shape != y.shape or list(x.columns) != list(y.columns):
        return False
    return all(np.allclose(x[c], y[c], rtol=1e-12, atol=1e-12, equal_nan=True)
               if pd.api.types.is_numeric_dtype(x[c]) and pd.api.types.is_numeric_dtype(y[c])
               else x[c].astype(str).equals(y[c].astype(str)) for c in x.columns)


def inputs() -> None:
    say("== 1. inputs vs Fig6_Reproduction/reference/manifest.json")
    targets = {Path(t["file"]).name: t["sha256"] for t in json.loads(MANIFEST.read_text())["targets"]}
    for f in sorted(INPUT.glob("*.csv")):
        check(hashlib.sha256(f.read_bytes()).hexdigest() == targets.get(f.name), f"{f.name} sha256 matches the reference manifest")
    check(sorted(p.name for p in INPUT.glob("*.csv")) == sorted(targets), "input/ holds exactly the eight reference tables")


def dotplot() -> None:
    say("\n== 2. dot plot vs v1")
    stem = "fig6_ora_cov_pls_motifs_transposed_v2_pctl_fg4_q05"
    check(png_diff(OUT / f"{stem}.png", V1_OUT / f"{stem}.png") == 0, f"{stem}.png pixel-identical to v1")
    check(same_table(pd.read_csv(OUT / f"{stem}_points.csv"), pd.read_csv(V1_OUT / f"{stem}_points.csv")), f"{stem}_points.csv identical to v1")


def panel_e() -> None:
    say("\n== 3. 6e vs v1")
    check((INPUT / "canonical9_ols_results.csv").read_bytes()
          == (V1_INPUT / "ENRegressionScreen.OLS.individual_feature_regression_results.Run3.csv").read_bytes(),
          "canonical9_ols_results.csv byte-identical to v1's Run3 OLS table")
    a = pd.read_csv(OUT / "fig6e_individual_partial_r2_values.csv")
    b = pd.read_csv(V1_OUT / "fig6e_individual_partial_r2_run3_values.csv")
    check(same_table(a.drop(columns="family"), b.drop(columns="family")), "6e values identical to v1 (family column aside)")
    check(a["family"].map(RUN3_FAMILY).equals(b["family"]), "6e family column = v1's under the family renames")
    for stem, v1 in [("fig6e_individual_partial_r2", "fig6e_individual_partial_r2_run3"), ("fig6de_composite", "fig6de_composite_run3")]:
        n = png_diff(OUT / f"{stem}.png", V1_OUT / f"{v1}.png")
        say(f"    {stem}.png vs v1: {n} pixels differ (6e family label renamed; the composite also contains 6d)")


def regression() -> None:
    say("\n== 4. 6d and full-model CV R^2: v4 vs Run3")
    v4 = pd.read_csv(INPUT / "cv_summary_results.csv")
    r3 = pd.read_csv(V1_INPUT / "ENRegressionScreen.cv_summary_results.Run3.csv")
    d4 = pd.read_csv(OUT / "fig6d_drop_one_family_delta_r2_values.csv")
    d1 = pd.read_csv(V1_OUT / "fig6d_drop_one_family_delta_r2_run3_values.csv")
    m = d4.assign(family=d4.family.map(RUN3_FAMILY)).merge(d1, on=["family", "outcome"], suffixes=("_v4", "_run3"))
    check(len(m) == len(d4) == len(d1), f"6d: all {len(d4)} cells matched to v1 through the family renames")
    m["delta"] = m.delta_r2_v4 - m.delta_r2_run3
    top = m.reindex(m.delta.abs().sort_values(ascending=False).index).head(5)
    say(f"    6d ΔR² change: max |Δ| {m.delta.abs().max():.4f}, median {m.delta.abs().median():.4f}; largest:")
    for r in top.itertuples():
        say(f"      {r.family:24s} {r.outcome:20s} {r.delta_r2_run3:+.4f} -> {r.delta_r2_v4:+.4f}")
    f4 = pd.read_csv(OUT / "fig6_full_model_cv_r2_values.csv").set_index("outcome")
    f1 = pd.read_csv(V1_OUT / "fig6_full_model_cv_r2_run3_values.csv").set_index("outcome")
    check(np.allclose(f4.r2_celltype_only, f1.r2_celltype_only.reindex(f4.index), rtol=0, atol=1e-14),
          "full-model figure: cell-type-only diamonds unchanged")
    say("    full-model CV R² (mean): Run3 -> v4")
    for o in f4.index:
        say(f"      {o:20s} {f1.loc[o, 'r2_full_mean']:.4f} -> {f4.loc[o, 'r2_full_mean']:.4f}  ({f4.loc[o, 'r2_full_mean'] - f1.loc[o, 'r2_full_mean']:+.4f})")
    # v4's full model leaves out the general-pELS motifs; Run3's model without them is full minus that drop-one delta.
    g = r3[(r3.level == "drop_one_delta") & (r3.family == "pels_motif_general")].set_index("outcome")
    full4 = v4[v4.level == "full_model"].set_index("outcome").r2_mean
    err = (full4 - (g.r2_full_mean - g.r2_mean).reindex(full4.index)).abs().max()
    check(err < 1e-12, f"v4 full model = Run3's full model without the general-pELS motifs (max |Δ| {err:.1e})")


def family_only() -> None:
    say("\n== 5. family-only models: v4 vs Run3")
    v4 = pd.read_csv(INPUT / "cv_summary_results.csv")
    r3 = pd.read_csv(V1_INPUT / "ENRegressionScreen.cv_summary_results.Run3.csv")
    a = v4[v4.level == "family_only"].assign(run3=lambda d: d.family.map(RUN3_FAMILY))
    b = r3[r3.level == "family_only"]
    m = a.merge(b, left_on=["run3", "outcome"], right_on=["family", "outcome"], suffixes=("_v4", "_run3"))
    cols = ["r2_mean", "r2_sd", "r2_min", "r2_max"]
    err = max(float((m[f"{c}_v4"] - m[f"{c}_run3"]).abs().max()) for c in cols)
    check(sorted(m.family_v4.unique()) == ["enhancer_burden", "gene_length", "promoter_architecture", "rna_stability"]
          and err < 1e-12, f"4 families also in Run3 match it ({len(m)} rows, max |Δ| {err:.1e})")
    say(f"    new in v4: {sorted(set(a.family) - set(m.family_v4))} (Run3 had no family-only model for matched_enhancer_motifs)")


def main() -> int:
    inputs()
    dotplot()
    panel_e()
    regression()
    family_only()
    say(f"\n{'ALL CHECKS PASSED' if not failures else 'FAILED: ' + '; '.join(failures)}")
    (OUT / "comparison_to_v1.txt").write_text("\n".join(lines) + "\n")
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
