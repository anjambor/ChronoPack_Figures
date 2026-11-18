#!/usr/bin/env python
"""Compare the r9 rebuild (Fig. 4 comparator runs) with the capsule figure's table (input/capsule/chronopack_specific_ligand_recovery.csv).

Reports the ligand sets (capsule 73 vs r9 49 ChronoPack ligands), rank-percentile agreement on shared ligands, the
top-14 overlap and the ligands whose specificity changes most. Writes output/comparison_to_capsule.txt.
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import common as c  # noqa: E402

import pandas as pd  # noqa: E402
from scipy.stats import spearmanr  # noqa: E402


def main() -> None:
    old = pd.read_csv(c.CAPSULE_TABLE)
    new = pd.read_csv(c.TABLE_DIR / "ligand_specificity_fig4_runs.csv")
    lines = [f"Capsule: {len(old)} ChronoPack ligands; r9 rebuild: {len(new)}"]
    only_old, only_new = sorted(set(old.ligand) - set(new.ligand)), sorted(set(new.ligand) - set(old.ligand))
    lines.append(f"Only in capsule ({len(only_old)}): {', '.join(only_old)}")
    lines.append(f"Only in r9 ({len(only_new)}): {', '.join(only_new) or '-'}")
    m = old.merge(new, on="ligand", suffixes=("_capsule", "_r9"))
    lines.append(f"\nShared ligands: {len(m)}")
    for col in ["chronopack_rank_pct", "best_other_rank_pct_filled", "chronopack_specificity"]:
        rho = spearmanr(m[f"{col}_capsule"], m[f"{col}_r9"]).statistic
        lines.append(f"  Spearman {col}: {rho:.3f}")
    n = c.N_LABEL
    top_old, top_new = (list(t.sort_values("chronopack_specificity", ascending=False).head(n).ligand)
                        for t in (old, new))
    lines.append(f"\nTop {n} by specificity")
    lines.append(f"  capsule: {', '.join(top_old)}")
    lines.append(f"  r9:      {', '.join(top_new)}")
    lines.append(f"  capsule & r9: {len(set(top_old) & set(top_new))}; "
                 f"capsule top {n} not among r9 ligands: {', '.join(l for l in top_old if l in only_old) or '-'}")
    lines.append(f"\nLigands with specificity > 0: capsule {(old.chronopack_specificity > 0).sum()}/{len(old)}, "
                 f"r9 {(new.chronopack_specificity > 0).sum()}/{len(new)}")
    cols = ["ligand", "chronopack_ligand_rank_capsule", "chronopack_ligand_rank_r9", "best_other_label_capsule",
            "best_other_label_r9", "chronopack_specificity_capsule", "chronopack_specificity_r9"]
    m["delta"] = m["chronopack_specificity_r9"] - m["chronopack_specificity_capsule"]
    lines.append("\nShared ligands, largest specificity change (r9 - capsule):")
    with pd.option_context("display.width", 250):
        lines.append(m.reindex(m["delta"].abs().sort_values(ascending=False).index).head(12)[cols + ["delta"]]
                     .round(3).to_string(index=False))
        lines.append("\nShared ligands in either top list:")
        lines.append(m[m.ligand.isin(set(top_old) | set(top_new))].sort_values("chronopack_specificity_r9",
                                                                               ascending=False)[cols]
                     .round(3).to_string(index=False))
    text = "\n".join(lines)
    (c.OUTPUT_DIR / "comparison_to_capsule.txt").write_text(text + "\n")
    print(text)
    print("COMPARE_OK")


if __name__ == "__main__":
    main()
