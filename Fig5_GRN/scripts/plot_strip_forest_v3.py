#!/usr/bin/env python
"""5b EPR / AUPRCr strip-forests, v3: the v2 plots (plot_publication_v4_epr_strip_forest.py) with
TF labels placed to avoid clashes.

Medians, IQR/range bars, axes, ticks and the clipped-value annotations are drawn exactly as in v2.
Dots start at v2's positions (row + per-TF offset); after the layout is fixed, dots in the same row
whose x positions overlap each other or the median diamond (e.g. several TFs at EPR = 0) are spread
vertically in TF order within the row (a median stays at the row centre with the dots split around
it). Each label then gets a spot beside its dot (right as in v2, shifted up/down, left, or a little
farther out) that avoids dots, medians, the legend, clipped-value annotations and other labels; a thin
leader line joins a label to its dot when it is not directly beside it. See place_labels for the
search. Spots that still overlap or cross are reported in the run output.

Settings differ per benchmark (dicts below): K562 ChIP separates dots only as far as needed; TRRUST spaces
dots one label height apart. Both keep v2's figure size (Fig5_GRN_v3, user request; TRRUST was drawn at
10.6 x 6.8 in with greedy placement only), let labels reach slightly into the gap between rows and refine
the placement.

Writes publication_strip_forest_{benchmark}_per_tf_EPR_AUPRCr_all_methods_all_labeled_v3.pdf to 5b/.
"""
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import FormatStrFormatter
import numpy as np
import pandas as pd

import benchmark_grn_systems as b
import plot_publication_v3_trrust_strip_forest as strip
import plot_publication_v4_epr_strip_forest as v2
from plot_publication_v2_bounded_violins import format_outlier_value

OUT_DIR = v2.OUT_DIR
LABEL_FONTSIZE = 7.6
# Labels stay within row_y +/- ROW_BAND (data units; rows are 0.86 apart). K562 ChIP (original height)
# lets labels reach slightly into the gap toward the next row when that space is free.
ROW_BAND = {"K562_ChIP": 0.50, "TRRUST": 0.50}  # Fig5_GRN_v3: TRRUST 0.43 -> 0.50 at v2's height
# Dot spreading per benchmark: band (spread dots stay within row_y +/- band, data units) and target spacing
# ("label": one label height, so labels can sit level with their dots; "dot": just clear of each other).
DOT_SPREAD = {"K562_ChIP": (0.30, "dot"), "TRRUST": (0.36, "label")}
# Refine the greedy label placement by repeated best-response moves (see place_labels).
REFINE_LABELS = {"K562_ChIP": True, "TRRUST": True}  # Fig5_GRN_v3: TRRUST refined too at v2's height
DOT_MARGIN_PT = 1.0  # clearance around dots and medians, points
LABEL_GAP_PT = 2.4  # gap between a dot's edge and its label (v2: ~2% of the x range)
LABEL_CLEARANCE_PT = 1.5  # minimum horizontal space between two labels (so names never run together)
LEADER_COLOR = "#8c8c8c"
# Figure size per benchmark. K562 ChIP keeps v2's size; TRRUST is taller, giving rows room to stack labels.
FIGSIZE = {"K562_ChIP": (10.6, 4.65), "TRRUST": (10.6, 4.65)}  # Fig5_GRN_v3: TRRUST at v2's height (user request; was 6.8)


def segment_hits_box(p0, p1, box) -> bool:
    """True if the segment p0-p1 passes through the (x0, y0, x1, y1) box (Liang-Barsky clip)."""
    (x0, y0), (x1, y1) = p0, p1
    dx, dy = x1 - x0, y1 - y0
    t0, t1 = 0.0, 1.0
    for p, q in ((-dx, x0 - box[0]), (dx, box[2] - x0), (-dy, y0 - box[1]), (dy, box[3] - y0)):
        if p == 0:
            if q < 0:
                return False
        else:
            t = q / p
            if p < 0:
                t0 = max(t0, t)
            else:
                t1 = min(t1, t)
            if t0 > t1:
                return False
    return True


def segments_cross(a0, a1, b0, b1) -> bool:
    def orient(p, q, r):
        return np.sign((q[0] - p[0]) * (r[1] - p[1]) - (q[1] - p[1]) * (r[0] - p[0]))
    return (orient(a0, a1, b0) != orient(a0, a1, b1)) and (orient(b0, b1, a0) != orient(b0, b1, a1))


def method_points(data: pd.DataFrame, benchmark: str, metric: str, xmin: float, xmax: float, y_positions: dict) -> list[dict]:
    """One record per drawn dot, with the same position, size and clipping as v2."""
    offsets = strip.tf_offsets(strip.ordered_tfs_for_benchmark(benchmark))
    out = []
    for method in strip.METHOD_ORDER:
        rows = data[data["method"].eq(method)]
        for _, row in rows.iterrows():
            true_x = float(row[metric])
            if not np.isfinite(true_x):
                continue
            out.append({
                "method": method, "TF": row["TF"], "true_x": true_x,
                "x": float(np.clip(true_x, xmin, xmax)), "y": y_positions[method] + offsets.get(row["TF"], 0.0),
                "row_y": y_positions[method], "high_clip": true_x > xmax, "low_clip": true_x < xmin,
                "size": 64 if "ChronoPack" in method else 50,
            })
    return out


def draw_metric(ax, data, benchmark, metric, title, xmin, xmax, y_positions, clipped_rows, summary_rows) -> list[dict]:
    """v2's plot_metric without the plain point labels (clipped-value annotations are kept)."""
    x_range = xmax - xmin
    edge_pad = x_range * 0.018
    ax.axvline(1, color="#9a9a9a", lw=0.9, ls="--", alpha=0.8, zorder=0)
    ax.grid(axis="x", color="#e7e7e7", linewidth=0.65)
    ax.set_axisbelow(True)
    medians = []
    for method in strip.METHOD_ORDER:
        rows = data[data["method"].eq(method)]
        if rows.empty:
            continue
        row_y = y_positions[method]
        values = rows[metric].replace([np.inf, -np.inf], np.nan).dropna()
        clipped_values = values.clip(lower=xmin, upper=xmax)
        if values.empty:
            continue
        q1, q3 = np.percentile(values, [25, 75])
        ax.hlines(row_y, max(xmin, float(q1)), min(xmax, float(q3)), color="#2a2a2a", linewidth=1.35, alpha=0.22, zorder=1)
        ax.hlines(row_y, float(clipped_values.min()), float(clipped_values.max()), color="#2a2a2a", linewidth=0.85, alpha=0.13, zorder=1)
        median_val = float(values.median())
        median_plot = float(np.clip(median_val, xmin, xmax))
        summary_rows.append({"benchmark": benchmark, "method": method, "metric": title, "statistic": "median",
                             "value": median_val, "plotted_value": median_plot, "n_TFs": int(values.shape[0])})
        ax.scatter(median_plot, row_y, marker="D", s=82, facecolor="white", edgecolor="#111111", linewidth=1.2, zorder=6)
        medians.append({"x": median_plot, "y": row_y, "size": 82})

    points = method_points(data, benchmark, metric, xmin, xmax, y_positions)
    for p in points:
        color = strip.softened_color(b.METHOD_COLORS[p["method"]], mix=0.10 if "ChronoPack" in p["method"] else 0.22)
        marker = ">" if p["high_clip"] else "<" if p["low_clip"] else "o"
        p["artist"] = ax.scatter(p["x"], p["y"], s=p["size"], marker=marker, color=color,
                   edgecolor="#111111" if p["high_clip"] or p["low_clip"] else "white",
                   linewidth=0.65, clip_on=False, zorder=4)
        if p["high_clip"] or p["low_clip"]:
            clipped_rows.append({"benchmark": benchmark, "method": p["method"], "TF": p["TF"], "metric": title,
                                 "true_value": p["true_x"], "plotted_value": p["x"]})
        if p["high_clip"]:
            p["annotation"] = ax.annotate(
                f"{p['TF']} {format_outlier_value(p['true_x'])}", xy=(xmax, p["y"]), xytext=(xmax + edge_pad, p["y"]),
                ha="left", va="center", fontsize=8.2, color="#222222",
                arrowprops={"arrowstyle": "-", "color": "#777777", "lw": 0.55}, clip_on=False, zorder=8)

    ax.set_title("")
    ax.set_xlabel(title, fontsize=12.2)
    ax.set_xlim(xmin, xmax)
    if metric in {"AUPRC_ratio", "EPR"} and xmin < 1.0:
        ax.set_xticks(np.arange(1.0, xmax + 0.001, 1.0))
    ticks = v2.X_TICKS.get((benchmark, metric))
    if ticks is not None:
        ax.set_xticks(ticks)
        if not all(float(t).is_integer() for t in ticks):  # Fig5: integer ticks (shared-universe TRRUST) keep plain labels
            ax.xaxis.set_major_formatter(FormatStrFormatter("%.1f"))
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(False)
    ax.tick_params(axis="y", length=0)
    ax.tick_params(axis="x", labelsize=10.5, bottom=True, length=3.5, width=0.8, color=ax.spines["bottom"].get_edgecolor())
    return points, medians


def spread_dots(fig, ax, points: list[dict], medians: list[dict], dot_band: float, spacing: str) -> None:
    """Spread dots of one row that overlap in x (with each other or the row's median diamond) vertically.

    Dots keep their TF order; spacing is one label height where the DOT_BAND allows. A median diamond in
    a cluster stays fixed at the row centre, with the cluster's dots split above and below it.
    """
    renderer = fig.canvas.get_renderer()
    px_per_pt = fig.dpi / 72.0
    to_disp, to_data = ax.transData.transform, ax.transData.inverted().transform
    probe = ax.text(0, 0, "RELA", fontsize=LABEL_FONTSIZE)
    h = probe.get_window_extent(renderer).height
    probe.remove()
    if spacing == "dot":
        h = (np.sqrt(max(p["size"] for p in points)) + DOT_MARGIN_PT) * px_per_pt
    rows: dict[float, list[dict]] = {}
    for p in points:
        if not (p["high_clip"] or p["low_clip"]):
            rows.setdefault(p["row_y"], []).append(p)
    median_by_row = {m["y"]: m for m in medians}
    for row_y, pts in rows.items():
        _, y_a = to_disp((0, row_y - dot_band))
        _, y_b = to_disp((0, row_y + dot_band))
        band_lo, band_hi = min(y_a, y_b), max(y_a, y_b)
        members = list(pts)
        if row_y in median_by_row:
            m = median_by_row[row_y]
            members.append({"x": m["x"], "y": row_y, "size": m["size"], "fixed": True})
        members.sort(key=lambda p: p["x"])
        clusters, current = [], [members[0]]
        for p in members[1:]:
            q = current[-1]
            reach = ((np.sqrt(p["size"]) + np.sqrt(q["size"])) / 2 + DOT_MARGIN_PT) * px_per_pt
            if to_disp((p["x"], 0))[0] - to_disp((q["x"], 0))[0] < reach:
                current.append(p)
            else:
                clusters.append(current)
                current = [p]
        clusters.append(current)
        for cl in clusters:
            dots = [p for p in cl if not p.get("fixed")]
            fixed = [p for p in cl if p.get("fixed")]
            if len(cl) < 2 or not dots:
                continue
            dots.sort(key=lambda p: p["y"])  # TF order; data y grows downwards on screen
            if fixed:
                k = len(dots) // 2  # median in the middle of the stack, at the row centre
                order = dots[:k] + fixed + dots[k:]
                y_k = to_disp((0, row_y))[1]
                above, below = k, len(order) - 1 - k
                step = h
                if above:
                    step = min(step, (band_hi - y_k) / above)
                if below:
                    step = min(step, (y_k - band_lo) / below)
                new_disp = [y_k + (k - i) * step for i in range(len(order))]
            else:
                order = dots
                step = min(h, (band_hi - band_lo) / (len(order) - 1))
                centre = float(np.mean([to_disp((0, p["y"]))[1] for p in order]))
                half = step * (len(order) - 1) / 2
                centre = min(max(centre, band_lo + half), band_hi - half)
                new_disp = [centre + half - i * step for i in range(len(order))]
            for p, yd in zip(order, new_disp):
                if p.get("fixed"):
                    continue
                p["y"] = float(to_data((to_disp((p["x"], 0))[0], yd))[1])
                p["artist"].set_offsets([[p["x"], p["y"]]])


def overlap_area(a, b) -> float:
    w = min(a[2], b[2]) - max(a[0], b[0])
    h = min(a[3], b[3]) - max(a[1], b[1])
    return w * h if w > 0 and h > 0 else 0.0


def place_labels(fig, ax, points: list[dict], medians: list[dict], extra_obstacles: list, refine: bool = False,
                 row_band: float = 0.43) -> list[str]:
    """Label placement in display (pixel) coordinates; returns warnings for unresolved overlaps.

    Each label has candidate spots beside its dot (right, shifted up/down, left, a little farther out),
    inside the axes and its row band. Labels are first placed one at a time (crowded dots first) at
    the cheapest spot that is clear of everything placed so far. With ``refine``, every label is then
    revisited and moved to its best spot given all the others, repeating until nothing improves; the
    score penalises overlaps (heavily), crossings of a leader with a label, dot or other leader
    (strongly) and distance from the dot (lightly).
    """
    renderer = fig.canvas.get_renderer()
    px_per_pt = fig.dpi / 72.0
    to_disp = ax.transData.transform
    to_data = ax.transData.inverted().transform
    ax_box = ax.get_window_extent(renderer)
    gap = LABEL_GAP_PT * px_per_pt
    clear = LABEL_CLEARANCE_PT * px_per_pt

    def square(x, y, size_pt2, margin_pt):
        cx, cy = to_disp((x, y))
        r = (np.sqrt(size_pt2) / 2 + margin_pt) * px_per_pt
        return (cx - r, cy - r, cx + r, cy + r)

    dot_boxes = {id(p): square(p["x"], p["y"], p["size"], DOT_MARGIN_PT) for p in points}
    fixed = list(dot_boxes.values()) + [square(m["x"], m["y"], m["size"], DOT_MARGIN_PT) for m in medians] + list(extra_obstacles)

    labelled = [p for p in points if not p["high_clip"] and not p["low_clip"]]
    disp = {id(p): to_disp((p["x"], p["y"])) for p in labelled}
    texts, sizes = {}, {}
    for p in labelled:
        t = ax.text(0, 0, p["TF"], ha="left", va="center", fontsize=LABEL_FONTSIZE, color="#2d2d2d", clip_on=False, zorder=5)
        bb = t.get_window_extent(renderer)
        texts[id(p)], sizes[id(p)] = t, (bb.width, bb.height)
    h_ref = float(np.median([hh for _w, hh in sizes.values()])) if sizes else 1.0

    def candidates_for(p):
        w, h = sizes[id(p)]
        cx, cy = disp[id(p)]
        r = np.sqrt(p["size"]) / 2 * px_per_pt
        _, a = to_disp((0, p["row_y"] + row_band))  # y axis is inverted: larger data y is lower
        _, b_ = to_disp((0, p["row_y"] - row_band))
        band_lo, band_hi = min(a, b_), max(a, b_)
        out = []
        dys = [0.0] + [sg * k * 0.25 * h for k in range(1, 17) for sg in (1, -1)]
        for side in ("right", "left"):
            for dx_extra in (0.0, 0.35 * w, 0.7 * w, w + 2 * gap):
                for dy in dys:
                    x0 = cx + r + gap + dx_extra if side == "right" else cx - r - gap - dx_extra - w
                    box = (x0, cy + dy - h / 2, x0 + w, cy + dy + h / 2)
                    if box[0] < ax_box.x0 or box[2] > ax_box.x1 + 2 * px_per_pt or box[1] < band_lo or box[3] > band_hi:
                        continue
                    tx, ty = min(max(cx, box[0]), box[2]), min(max(cy, box[1]), box[3])
                    d = np.hypot(tx - cx, ty - cy)
                    seg = None if d <= r + gap + 1 else ((cx + (tx - cx) * r / d, cy + (ty - cy) * r / d), (tx, ty))
                    cost = abs(dy) / h + dx_extra / w * 1.2 + (0.6 if side == "left" else 0.0)
                    out.append((cost, box, seg))
        if not out:  # nothing fits the band/axes: v2's default spot
            box = (cx + r + gap, cy - h / 2, cx + r + gap + w, cy + h / 2)
            out.append((0.0, box, None))
        out.sort(key=lambda c: c[0])
        return out

    cands = {id(p): candidates_for(p) for p in labelled}
    own_dot = {id(p): dot_boxes[id(p)] for p in labelled}
    W_OVERLAP, W_CROSS = 1.0, 40.0  # per px^2 of overlap; per crossing

    def clash(pid, cand, others):
        """Overlap area + weighted crossings of one candidate against fixed obstacles and other labels."""
        _cost, box, seg = cand
        ov = sum(overlap_area(box, o) for o in fixed)
        crosses = 0
        padded = (box[0] - clear, box[1], box[2] + clear, box[3])
        for qid, (_c2, box2, seg2) in others:
            if qid == pid:
                continue
            ov += overlap_area(padded, box2)
            if seg2 is not None and segment_hits_box(seg2[0], seg2[1], box):
                crosses += 1
            if seg is not None:
                if segment_hits_box(seg[0], seg[1], box2):
                    crosses += 1
                if seg2 is not None and segments_cross(seg[0], seg[1], seg2[0], seg2[1]):
                    crosses += 1
        if seg is not None:
            crosses += sum(segment_hits_box(seg[0], seg[1], o) for o in fixed if o is not own_dot[pid])
        return ov, crosses

    def score(pid, cand, others):
        ov, crosses = clash(pid, cand, others)
        return W_OVERLAP * ov + W_CROSS * crosses + cand[0], ov

    # 1. Greedy: crowded dots first, each at its cheapest spot clear of everything placed so far.
    crowd = {id(p): sum(np.hypot(*(disp[id(p)] - disp[id(q)])) < 2.5 * h_ref for q in labelled) for p in labelled}
    order = sorted(labelled, key=lambda p: (-crowd[id(p)], p["x"], p["y"]))
    chosen: dict[int, tuple] = {}
    for p in order:
        placed = list(chosen.items())
        best, best_key = None, None
        for cand in cands[id(p)]:
            ov, crosses = clash(id(p), cand, placed)
            key = (ov > 0, crosses > 0, ov, crosses, cand[0])
            if best_key is None or key < best_key:
                best, best_key = cand, key
            if ov == 0 and crosses == 0:
                break
        chosen[id(p)] = best

    # 2. Optional refinement: move each label to its best spot given all others until stable.
    if refine:
        for _sweep in range(30):
            improved = False
            for p in order:
                pid = id(p)
                others = list(chosen.items())
                current, _ = score(pid, chosen[pid], others)
                for cand in cands[pid]:
                    sc, _ = score(pid, cand, others)
                    if sc < current - 1e-9:
                        chosen[pid], current, improved = cand, sc, True
            if not improved:
                break

    warnings = []
    for p in labelled:
        _cost, box, seg = chosen[id(p)]
        ov, crosses = clash(id(p), chosen[id(p)], list(chosen.items()))
        if ov > 0 or crosses > 0:
            warnings.append(f"{p['method']} {p['TF']}: overlap {ov / px_per_pt ** 2:.0f} pt^2, {crosses} leader crossing(s)")
        texts[id(p)].set_position(to_data((box[0], (box[1] + box[3]) / 2)))
        if seg is not None:
            (x_a, y_a), (x_b, y_b) = to_data(seg[0]), to_data(seg[1])
            ax.add_line(Line2D([x_a, x_b], [y_a, y_b], color=LEADER_COLOR, lw=0.5, zorder=3.5, clip_on=False))
    return warnings


def plot_benchmark(all_data: pd.DataFrame, benchmark: str) -> Path:
    config = v2.BENCHMARK_CONFIGS[benchmark]
    data = all_data[all_data["benchmark"].eq(benchmark) & all_data["method"].isin(b.METHOD_ORDER)].copy()
    y_values, fig_height, (k562_label_y, thp1_label_y) = b.strip_rows(b.METHOD_ORDER, FIGSIZE[benchmark][1])  # Fig5_GRN_v3
    y_positions = dict(zip(b.METHOD_ORDER, y_values))

    fig, axes = plt.subplots(1, 2, figsize=(FIGSIZE[benchmark][0], fig_height), sharey=True)
    clipped_rows, summary_rows, drawn = [], [], []
    for ax, (metric, title, xmin, xmax, _threshold) in zip(axes, config["metrics"]):
        drawn.append(draw_metric(ax, data, benchmark, metric, title, xmin, xmax, y_positions, clipped_rows, summary_rows))
    for ax in axes:
        ax.set_yticks([y_positions[m] for m in b.METHOD_ORDER])
        ax.set_yticklabels([strip.method_label(m) for m in b.METHOD_ORDER], fontsize=11.0)
        ax.set_ylim(max(y_values) + 0.55, -0.62)
    axes[0].text(-0.30, k562_label_y, "K562", transform=axes[0].transAxes, ha="center", va="center", rotation=90, fontsize=9.8, color="#555555")
    axes[0].text(-0.30, thp1_label_y, "THP1", transform=axes[0].transAxes, ha="center", va="center", rotation=90, fontsize=9.8, color="#555555")
    legend = axes[1].legend(
        handles=[Line2D([0], [0], marker="D", linestyle="", markerfacecolor="white", markeredgecolor="#111111", markersize=7.8, label="Median")],
        frameon=False, loc="lower right", fontsize=10.0, handletextpad=0.4)
    fig.suptitle(config["label"], fontsize=16.2, y=0.985)
    fig.tight_layout(rect=[0.055, 0.02, 1, 0.945], w_pad=2.0)

    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    warnings = []
    for ax, (points, medians) in zip(axes, drawn):
        spread_dots(fig, ax, points, medians, *DOT_SPREAD[benchmark])
    for i, (ax, (points, medians)) in enumerate(zip(axes, drawn)):
        extra = [tuple(p["annotation"].get_window_extent(renderer).extents) for p in points if "annotation" in p]
        if i == 1:
            extra.append(tuple(legend.get_window_extent(renderer).extents))
        warnings += place_labels(fig, ax, points, medians, extra, refine=REFINE_LABELS[benchmark],
                                 row_band=ROW_BAND[benchmark])

    out = OUT_DIR / f"publication_strip_forest_{benchmark}_per_tf_EPR_AUPRCr_all_methods_all_labeled_v3.pdf"
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, bbox_inches="tight")
    if PREVIEW_DIR is not None:
        fig.savefig(Path(PREVIEW_DIR) / out.with_suffix(".png").name, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(out)
    for w in warnings:
        print("  unresolved:", w)
    return out


PREVIEW_DIR = None  # set to a directory to also write PNG previews there


def main() -> None:
    strip.OUT_DIR = OUT_DIR
    strip.ordered_tfs_for_benchmark = v2.ordered_tfs_for_benchmark
    data = v2.filter_rows(pd.read_csv(v2.METRICS_TABLE))
    for benchmark in ["TRRUST", "K562_ChIP"]:
        plot_benchmark(data, benchmark)


if __name__ == "__main__":
    main()
