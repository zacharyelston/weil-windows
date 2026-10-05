"""Extrapolation-protocol prediction for the planted 11a1 cases: LS line through ln σ₁ against x over the
cheap regime x/N ≤ 3 (five points), compared with the exact crossing."""
import json, sys
import mpmath as mp
mp.mp.dps = 20
S = "data/research_r1"
N = 11


def lsq(xs, ys):
    n = len(xs)
    mx, my = sum(xs) / n, sum(ys) / n
    a = sum((x - mx) * (y - my) for x, y in zip(xs, ys)) / sum((x - mx) ** 2 for x in xs)
    return a, my - a * mx


for sector in ("even", "odd"):
    try:
        d = json.load(open(f"{S}/plant_11a1_{sector}.json"))
    except Exception as e:
        print(sector, "missing", e)
        continue
    rows = d["rows"] if isinstance(d, dict) else d
    cross = d.get("crossings", {}) if isinstance(d, dict) else {}
    for key in rows[0]["cases"]:
        pts = [(float(r["x"]), float(r["cases"][key]["sigma1"])) for r in rows if r["cases"][key]["sigma1"] is not None]
        cheap = [(x, s) for x, s in pts if x <= 3 * N + 1e-9]
        a, b = lsq([x for x, _ in cheap], [mp.log(s) for _, s in cheap])
        x_ls = -b / a
        (xa, sa), (xb, sb) = cheap[-2], cheap[-1]
        a2 = (mp.log(sb) - mp.log(sa)) / (xb - xa)
        x_2 = xb - mp.log(sb) / a2
        c = cross.get(key)
        print(f"{sector:<4} {key:<10}: cheap-regime σ₁ = {[(x, round(s, 5)) for x, s in cheap]}; LS → x_c = {mp.nstr(x_ls, 5)} (x/N {mp.nstr(x_ls / N, 4)}), "
              f"last-two → {mp.nstr(x_2, 5)} (x/N {mp.nstr(x_2 / N, 4)}); exact crossing: {c}")
