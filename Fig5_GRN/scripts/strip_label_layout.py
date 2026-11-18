"""Sparse TF labels for the composite's 5b strip-forests, in the style of ref/FigureComposite.v6.png.

The reference labels were placed by hand: each label sits right beside its dot (right, a diagonal, left, above or
below; never at the end of a leader line), dots were nudged up or down within their row to make room, and dots
with no clear spot (inner dots of a stack at the same value, crowded clusters) were left unlabelled.

layout() reproduces that as one optimisation per axes. Each free dot (clipped dots keep their position and
annotation) chooses a vertical offset within its row (ROW_BAND data units of the row centre) and a label slot
(SLOTS, or no label). The cost of a layout adds:
  * HARD for any overlap: label-label (with clearance), label-dot, label-median/legend/annotation, dot-dot
    (centres closer than DOT_TOUCH x the sum of radii), dot-median; and for a label outside the axes or farther
    than LABEL_BAND from its row centre;
  * AMBIGUOUS for a label that is not clearly nearer its own dot than every other TF dot (so any label shown
    names the dot next to it; median diamonds are only obstacles);
  * the cost of leaving a label out (HIDE_COST; more for ChronoPack rows and RELA / NFKB1 / JUN, the TFs of 5d),
    so a label is dropped only when it cannot be placed cleanly;
  * LINE_COST for a label lying across a row's grey range line (the reference puts such labels on a diagonal);
  * small preferences: slot order (SLOTS: right first, as in the reference) and dots near their row centre.
Dots with identical values in a row (e.g. several TFs at EPR = 0) are not optimised: they form a fixed compact
stack centred on the row (STACK_SPACING x the dot diameter apart, in TF order; a median at the same value is drawn
over the stack's centre), as in the reference; only their labels are placed.
Simulated annealing over single-dot moves (seeded, so deterministic), then exhaustive best-response sweeps.
Geometry uses each label's ink box (TextPath extents), so labels can sit as close to dots as in the reference.
"""
from __future__ import annotations

import math

import numpy as np
from matplotlib.font_manager import FontProperties
from matplotlib.textpath import TextPath

# (name, preference cost); the reference mostly labels to the right or on a right diagonal
SLOTS = [("R", 0.0), ("UR", 0.25), ("DR", 0.25), ("L", 0.6), ("UL", 0.75), ("DL", 0.75), ("U", 1.0), ("D", 1.0)]
HIDDEN = len(SLOTS)
HARD = 1000.0
AMBIGUOUS = 30.0
HIDE_COST = 10.0
HIDE_BONUS = {"chrono": 4.0, "key_tf": 3.0}
KEY_TFS = {"RELA", "NFKB1", "JUN"}
DISPLACEMENT = 1.5  # cost of a dot at the edge of its band (quadratic in the offset)
ROW_BAND = 0.36  # dots stay within row +/- this (data units; rows are 0.86 apart)
LABEL_BAND = 0.52  # label centres stay within row +/- this
N_OFFSETS = 17
DOT_TOUCH = 0.6  # dot centres at least 0.6 x the sum of radii apart (the reference's stacks overlap by ~40%)
LINE_COST = 3.0
STACK_SPACING = 0.5  # tied dots: centre spacing / dot diameter (the reference's stacks overlap by about half)
GAP_PT = 1.6  # dot edge to label ink
LABEL_CLEAR_PT = (2.2, 1.4)  # label-label clearance (x, y)
DOT_CLEAR_PT = 0.8  # label-dot clearance
AMBIG_MARGIN_PT = 1.5
SWEEPS_PER_POINT = 2500


def ink_box_pt(text: str, fontsize: float, family: str) -> tuple[float, float, float, float]:
    """Ink extents (points) of text drawn with ha='left', va='baseline' at the origin."""
    e = TextPath((0, 0), text, size=fontsize, prop=FontProperties(family=family)).get_extents()
    return e.x0, e.y0, e.x1, e.y1


def _overlap(a, b) -> float:
    w = min(a[2], b[2]) - max(a[0], b[0])
    h = min(a[3], b[3]) - max(a[1], b[1])
    return w * h if w > 0 and h > 0 else 0.0


def _dist_to_box(x, y, box) -> float:
    dx = max(box[0] - x, 0.0, x - box[2])
    dy = max(box[1] - y, 0.0, y - box[3])
    return math.hypot(dx, dy)


def layout(fig, ax, points: list[dict], medians: list[dict], extra_obstacles: list, *, fontsize: float,
           family: str, color: str = "#2d2d2d", seed: int = 0) -> list[dict]:
    """Move free dots within their rows and add the labels that fit; returns one record per dot."""
    renderer = fig.canvas.get_renderer()
    px = fig.dpi / 72.0
    T = ax.transData
    inv = T.inverted()
    axb = ax.get_window_extent(renderer)
    gap, clr_x, clr_y = GAP_PT * px, LABEL_CLEAR_PT[0] * px, LABEL_CLEAR_PT[1] * px
    dclr, amb_m = DOT_CLEAR_PT * px, AMBIG_MARGIN_PT * px

    free = [p for p in points if not (p["high_clip"] or p["low_clip"])]
    n = len(free)
    cx = np.array([T.transform((p["x"], p["row_y"]))[0] for p in free])
    yc = np.array([T.transform((p["x"], p["row_y"]))[1] for p in free])
    r = np.array([math.sqrt(p["size"]) / 2 * px for p in free])

    def band_px(row_y, half):
        return abs(T.transform((0, row_y + half))[1] - T.transform((0, row_y))[1])

    offs = [np.linspace(-1, 1, N_OFFSETS) * band_px(p["row_y"], ROW_BAND) for p in free]
    lband = np.array([band_px(p["row_y"], LABEL_BAND) for p in free])
    ink = [tuple(v * px for v in ink_box_pt(p["TF"], fontsize, family)) for p in free]
    hide = np.array([HIDE_COST + (HIDE_BONUS["chrono"] if "ChronoPack" in p["method"] else 0.0)
                     + (HIDE_BONUS["key_tf"] if p["TF"] in KEY_TFS else 0.0) for p in free])

    # fixed obstacles: medians (diamonds as circles for dots, squares for labels), clipped dots, legend, annotations
    fixed_circles, fixed_boxes, fixed_dots = [], list(extra_obstacles), []  # fixed_dots: clipped dots (ambiguity)
    for m in medians:
        mx, my = T.transform((m["x"], m["y"]))
        mr = math.sqrt(m["size"]) / 2 * px
        fixed_circles.append((mx, my, mr))
        fixed_boxes.append((mx - mr, my - mr, mx + mr, my + mr))
    lines = []  # each row's range line (draw_metric: clipped min to max at the row centre), as thin boxes
    for row_y in {p["row_y"] for p in points}:
        xs = [T.transform((p["x"], row_y))[0] for p in points if p["row_y"] == row_y]
        ly = T.transform((0, row_y))[1]
        lines.append((min(xs), ly - 0.7 * px, max(xs), ly + 0.7 * px))
    band_max = np.array([o[-1] for o in offs])
    stacks: dict[tuple, list[int]] = {}  # tied dots: fixed compact stacks
    stack_id = np.full(n, -1)
    for i, p in enumerate(free):
        stacks.setdefault((p["row_y"], round(cx[i] * 2) / 2), []).append(i)
    for (row_y, x), members in stacks.items():
        if len(members) < 2:
            continue
        members.sort(key=lambda i: free[i]["y"])  # TF order, top of the row first
        sp = STACK_SPACING * 2 * r[members].max()
        for j, i in enumerate(members):
            offs[i] = np.array([sp * ((len(members) - 1) / 2 - j)])
            stack_id[i] = members[0]

    for p in points:
        if p["high_clip"] or p["low_clip"]:
            qx, qy = T.transform((p["x"], p["y"]))
            qr = math.sqrt(p["size"]) / 2 * px
            fixed_circles.append((qx, qy, qr))
            fixed_dots.append((qx, qy, qr))
            fixed_boxes.append((qx - qr, qy - qr, qx + qr, qy + qr))

    def label_box(i, y, s):
        x0i, y0i, x1i, y1i = ink[i]
        w, h = x1i - x0i, y1i - y0i
        c, ri = cx[i], r[i]
        name = SLOTS[s][0]
        if name == "R":
            bx, by = c + ri + gap, y - h / 2
        elif name == "L":
            bx, by = c - ri - gap - w, y - h / 2
        elif name == "UR":
            bx, by = c + 0.8 * ri + 0.5 * gap, y + 0.45 * ri
        elif name == "DR":
            bx, by = c + 0.8 * ri + 0.5 * gap, y - 0.45 * ri - h
        elif name == "UL":
            bx, by = c - 0.8 * ri - 0.5 * gap - w, y + 0.45 * ri
        elif name == "DL":
            bx, by = c - 0.8 * ri - 0.5 * gap - w, y - 0.45 * ri - h
        elif name == "U":
            bx, by = c - w / 2, y + ri + 0.5 * gap
        else:  # D
            bx, by = c - w / 2, y - ri - 0.5 * gap - h
        return (bx, by, bx + w, by + h)

    # neighbours: dots close enough in x to interact, in the same or an adjacent row
    max_w = max((b[2] - b[0] for b in ink), default=0.0)
    reach = 2 * max_w + 4 * r.max() + 4 * gap if n else 0.0
    row_ids = sorted({p["row_y"] for p in free})
    row_index = np.array([row_ids.index(p["row_y"]) for p in free])
    nbrs = [[j for j in range(n) if j != i and abs(cx[i] - cx[j]) < reach and abs(row_index[i] - row_index[j]) <= 1]
            for i in range(n)]

    def unary(i, k, s):
        y = yc[i] + offs[i][k]
        c = DISPLACEMENT * min(1.0, abs(offs[i][k]) / max(band_max[i], 1e-9)) ** 2 if len(offs[i]) > 1 else 0.0
        for fx, fy, fr in (fixed_circles if stack_id[i] < 0 else ()):  # a stack may sit under its median
            d = math.hypot(cx[i] - fx, y - fy)
            lim = DOT_TOUCH * (r[i] + fr)
            if d < lim:
                c += HARD * (0.2 + (lim - d) / lim)
        if s == HIDDEN:
            return c + hide[i], None
        box = label_box(i, y, s)
        c += SLOTS[s][1]
        if box[0] < axb.x0 + 1 or box[2] > axb.x1 + 2 * px or box[1] < axb.y0 + 1 or box[3] > axb.y1 - 1:
            c += HARD
        if abs((box[1] + box[3]) / 2 - yc[i]) > lband[i]:
            c += HARD
        h2 = (ink[i][3] - ink[i][1]) ** 2
        for fb in fixed_boxes:
            ov = _overlap(box, (fb[0] - dclr, fb[1] - dclr, fb[2] + dclr, fb[3] + dclr))
            if ov:
                c += HARD * (0.2 + ov / h2)
        c += LINE_COST * sum(_overlap(box, ln) > 0 for ln in lines)
        own = _dist_to_box(cx[i], y, box) - r[i]
        for fx, fy, fr in fixed_dots:  # a diamond is never mistaken for a TF's dot
            if _dist_to_box(fx, fy, box) - fr < own + amb_m:
                c += AMBIGUOUS
        return c, box

    def pair(i, yi, bi, j, yj, bj):
        c = 0.0
        d = math.hypot(cx[i] - cx[j], yi - yj)
        lim = DOT_TOUCH * (r[i] + r[j])
        if d < lim and (stack_id[i] < 0 or stack_id[i] != stack_id[j]):
            c += HARD * (0.2 + (lim - d) / lim)
        if bi is not None:
            if _dist_to_box(cx[j], yj, bi) < r[j] + dclr:
                c += HARD
            elif _dist_to_box(cx[j], yj, bi) - r[j] < _dist_to_box(cx[i], yi, bi) - r[i] + amb_m:
                c += AMBIGUOUS
        if bj is not None:
            if _dist_to_box(cx[i], yi, bj) < r[i] + dclr:
                c += HARD
            elif _dist_to_box(cx[i], yi, bj) - r[i] < _dist_to_box(cx[j], yj, bj) - r[j] + amb_m:
                c += AMBIGUOUS
        if bi is not None and bj is not None:
            ov = _overlap((bi[0] - clr_x, bi[1] - clr_y, bi[2] + clr_x, bi[3] + clr_y), bj)
            if ov:
                c += HARD * (0.2 + ov / (bi[3] - bi[1]) ** 2)
        return c

    rng = np.random.default_rng(seed)
    k = np.array([int(np.argmin(np.abs(o - (T.transform((0, p["y"]))[1] - yc[i])))) for i, (p, o) in enumerate(zip(free, offs))])
    s = np.zeros(n, dtype=int)
    cache = [unary(i, k[i], s[i]) for i in range(n)]  # (unary cost, label box)

    def local(i, ki, si, ui=None):
        u, bi = ui if ui is not None else unary(i, ki, si)
        yi = yc[i] + offs[i][ki]
        return u + sum(pair(i, yi, bi, j, yc[j] + offs[j][k[j]], cache[j][1]) for j in nbrs[i]), (u, bi)

    if n:
        iters = SWEEPS_PER_POINT * n
        temps = np.geomspace(40.0, 0.05, iters)
        for t in temps:
            i = int(rng.integers(n))
            move = rng.random()
            ki, si = k[i], s[i]
            if move < 0.45:
                si = int(rng.integers(HIDDEN + 1))
            elif move < 0.8:
                ki = int(np.clip(ki + rng.choice([-3, -2, -1, 1, 2, 3]), 0, len(offs[i]) - 1))
            else:
                ki, si = int(rng.integers(len(offs[i]))), int(rng.integers(HIDDEN + 1))
            if ki == k[i] and si == s[i]:
                continue
            old, _ = local(i, k[i], s[i], cache[i])
            new, u_new = local(i, ki, si)
            if new <= old or rng.random() < math.exp(-(new - old) / t):
                k[i], s[i], cache[i] = ki, si, u_new
        for _sweep in range(6):  # exhaustive best response
            improved = False
            for i in range(n):
                best, best_state = local(i, k[i], s[i], cache[i])[0], None
                for ki in range(len(offs[i])):
                    for si in range(HIDDEN + 1):
                        c, u = local(i, ki, si)
                        if c < best - 1e-9:
                            best, best_state = c, (ki, si, u)
                if best_state is not None:
                    k[i], s[i], cache[i] = best_state
                    improved = True
            if not improved:
                break

    records = []
    for i, p in enumerate(free):
        y = yc[i] + offs[i][k[i]]
        p["y"] = float(inv.transform((cx[i], y))[1])
        p["artist"].set_offsets([[p["x"], p["y"]]])
        cost, _ = local(i, k[i], s[i], cache[i])
        box = cache[i][1]
        if box is not None:
            ax_x, ax_y = inv.transform((box[0] - ink[i][0], box[1] - ink[i][1]))
            ax.text(ax_x, ax_y, p["TF"], ha="left", va="baseline", fontsize=fontsize, color=color, clip_on=False, zorder=5)
        records.append({"method": p["method"], "TF": p["TF"], "value": p["true_x"], "labelled": box is not None,
                        "slot": SLOTS[s[i]][0] if s[i] != HIDDEN else "", "dot_offset_rows": p["y"] - p["row_y"],
                        "residual_cost": round(cost - (hide[i] if box is None else 0.0), 2)})
    for p in points:
        if p["high_clip"] or p["low_clip"]:
            records.append({"method": p["method"], "TF": p["TF"], "value": p["true_x"], "labelled": "annotation" in p,
                            "slot": "clipped", "dot_offset_rows": p["y"] - p["row_y"], "residual_cost": 0.0})
    return records
