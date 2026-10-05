"""Rung 3 of the calibration ladder: predict each control's first failure from σ₁(x) in the certified-positive
regime only (x ≤ x_pos), by a least-squares line through ln σ₁ against x, and by the last two points."""
import json, sys, glob
import mpmath as mp
mp.mp.dps = 20
S = "data/research_r1"
XPOS = {"ftstar": 5.31217, "dh": 10.80490, "z1": 15.79984}
BRACKET = {"ftstar": "(5.31217, 5.36250)", "dh": "(10.805, 30.745)", "z1": "(15.800, 19.844)"}


def lsq(xs, ys):
    n = len(xs)
    mx, my = sum(xs) / n, sum(ys) / n
    a = sum((x - mx) * (y - my) for x, y in zip(xs, ys)) / sum((x - mx) ** 2 for x in xs)
    return a, my - a * mx


for fn in ("ftstar", "dh", "z1"):
    for sector in ("odd", "even") if fn == "ftstar" else ("even", "odd"):
        rows = []
        for f in sorted(glob.glob(f"{S}/ctrl_{fn}_{sector}*.json")):
            rows += json.load(open(f))
        pts = sorted({(float(r["x"]), float(r["sigma1"])) for r in rows if r.get("sigma1") not in (None, "0.0")})
        if not pts:
            continue
        cross = None
        for (x1, s1), (x2, s2) in zip(pts, pts[1:]):
            if s1 < 1 <= s2:
                cross = (x1, x2)
        pos = [(x, s) for x, s in pts if x <= XPOS[fn] + 1e-6]
        line = f"{fn:<6} {sector:<4} σ₁ crossing in grid: {cross};  certified bracket {BRACKET[fn]}"
        if len(pos) >= 2:
            a, b = lsq([x for x, _ in pos], [mp.log(s) for _, s in pos])
            x_ls = -b / a
            (xa, sa), (xb, sb) = pos[-2], pos[-1]
            a2 = (mp.log(sb) - mp.log(sa)) / (xb - xa)
            x_2 = xb - mp.log(sb) / a2
            line += f";  from x ≤ {XPOS[fn]} ({len(pos)} pts): LS line → x_c = {mp.nstr(x_ls, 5)} (slope {mp.nstr(a, 4)}), last-two → {mp.nstr(x_2, 5)}"
        print(line)
        print("        σ₁(x):", ", ".join(f"{x:g}:{s:.4g}" for x, s in pts))
