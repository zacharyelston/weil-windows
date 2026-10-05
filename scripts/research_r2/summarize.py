#!/usr/bin/env python3
"""Tables for docs/RESEARCH_R2_SAMPLING.md from data/research_r2/*.json (every quoted R2 number is printed here).

Usage: .venv/bin/python scripts/research_r2/summarize.py
"""
import glob
import json
import math
import os

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..")
DATA = os.path.join(ROOT, "data", "research_r2")


def load(name):
    p = os.path.join(DATA, name)
    return json.load(open(p)) if os.path.exists(p) else []


def fmt(v, d=3):
    return f"{v:.{d}g}"


def slope(xs, ys):
    xb, yb = sum(xs) / len(xs), sum(ys) / len(ys)
    return sum((a - xb) * (b - yb) for a, b in zip(xs, ys)) / sum((a - xb) ** 2 for a in xs)


def bound_table(rows, label, rate):
    """rate = 1 (e^{-c}) or 2 (e^{-2c}): report measured/bound range and the fitted c-power p of bound*e^{rate c}."""
    print(f"\n### {label}\n")
    print("| object | sector | n | x/q range | measured/bound (min .. max) | p in bound ~ c^p e^{-%dc} |" % rate)
    print("|---|---|---|---|---|---|")
    allr = []
    keys = []
    for r in rows:
        k = (r["D"], r["parity"])
        if k not in keys:
            keys.append(k)
    for D, par in keys:
        rs = [r for r in rows if r["D"] == D and r["parity"] == par and r.get("ratio_measured_over_bound")]
        if not rs:
            continue
        rat = [float(r["ratio_measured_over_bound"]) for r in rs]
        allr += rat
        cs = [float(r["c"]) for r in rs]
        lb = [math.log(float(r["bound"])) + rate * float(r["c"]) for r in rs]
        p = slope([math.log(c) for c in cs], lb) if len(rs) > 2 else float("nan")
        n = rs[0].get("n_index")
        print(f"| {'zeta' if D == 1 else 'chi_' + str(D)} | {par} | {n} | {rs[0]['xq']}..{rs[-1]['xq']} | {fmt(min(rat))} .. {fmt(max(rat))} | {p:.2f} |")
    if allr:
        print(f"\nAll {len(allr)} points: measured/bound in [{fmt(min(allr))}, {fmt(max(allr))}]; every ratio <= 1: {all(v <= 1 for v in allr)}")


def trial_table(rows):
    print("\n### Sharp E-map prolate trial functions (Rayleigh quotient in our basis, same N as the scan)\n")
    print("| object | sector | n | x/q | trial/measured | R_trial | R_measured |")
    print("|---|---|---|---|---|---|---|")
    keys = []
    for r in rows:
        k = (r["D"], r["parity"])
        if k not in keys:
            keys.append(k)
    allr = []
    for D, par in keys:
        rs = [r for r in rows if r["D"] == D and r["parity"] == par]
        tm = [float(r["trial_over_measured"]) for r in rs]
        allr += tm
        print(f"| {'zeta' if D == 1 else 'chi_' + str(D)} | {par} | {rs[0]['n_index']} | {rs[0]['xq']}..{rs[-1]['xq']} | "
              f"{fmt(min(tm))} .. {fmt(max(tm))} | {fmt(float(rs[0]['R_trial']))} .. {fmt(float(rs[-1]['R_trial']))} | "
              f"{fmt(float(rs[0]['R_measured']))} .. {fmt(float(rs[-1]['R_measured']))} |")
    if allr:
        srt = sorted(allr)
        med = srt[len(srt) // 2]
        print(f"\nAll {len(allr)} points: trial/measured in [{fmt(min(allr))}, {fmt(max(allr))}], median {fmt(med)}, "
              f"{sum(1 for v in allr if v <= 1.5)} of {len(allr)} below 1.5")


def leakage_tables(rows):
    zf = [r for r in rows if r.get("object") == "zeta fit points"]
    d1 = [r for r in rows if "n_rule" in r and r.get("object") != "zeta fit points"]
    g2 = [r for r in rows if "weight" in r]
    print("\n### Degree 1: Fuchs asymptotic against exact leakage\n")
    print("| object | sector | n | R_F range | R_F drift (x/q 3 -> last) | R_ex range | R_ex drift | n*_F | n*_ex |")
    print("|---|---|---|---|---|---|---|---|---|")
    for r in d1:
        rf = [float(v) for v in r["R_fuchs"]]
        re_ = [float(v) for v in r["R_exact"]]
        print(f"| {r['object']} | {r['sector']} | {r['n_rule']} | {fmt(min(rf))}..{fmt(max(rf))} | {r['R_fuchs_last_over_first_from3']} | "
              f"{fmt(min(re_))}..{fmt(max(re_))} | {r['R_exact_last_over_first_from3']} | {r['n_star_fuchs']} | {r['n_star_exact']} |")
    if d1:
        ch = [r for r in d1 if r["object"] != "D=1"]
        dev = [float(r["n_star_exact"]) - r["n_rule"] for r in ch]
        devF = [float(r["n_star_fuchs"]) - r["n_rule"] for r in ch]
        print(f"\nCharacters: n*_F - n in [{min(devF):.2f}, {max(devF):.2f}]; n*_ex - n in [{min(dev):.2f}, {max(dev):.2f}]")
    if zf:
        print("\n### zeta's fit points (best upper bounds and certified lower bounds) against the exact leakage\n")
        print("| sector | n | x | R_F (upper bound) | R_ex (upper bound) | R_ex (certified lower bound) |")
        print("|---|---|---|---|---|---|")
        for r in zf:
            print(f"| {r['sector']} | {r['n_rule']} | {', '.join(r['x'])} | {', '.join(r['R_fuchs_ub'])} | {', '.join(r['R_exact_ub'])} | "
                  f"{', '.join(v if v else '-' for v in r['R_exact_cert'])} |")
    print("\n### GL(2): continuous index against the exact Hankel-prolate leakage (nu = k - 1, c = 4 pi v)\n")
    print("| object | N | k | eps | sector | m | n*-k measured | n*-k Hankel prolate | shift |")
    print("|---|---|---|---|---|---|---|---|---|")
    for r in g2:
        print(f"| {r['object']} | {r['level']} | {r['weight']} | {r['eps']:+d} | {r['sector']} | {r['m']} | {r['n_star_minus_k_measured']} | "
              f"{r['n_star_minus_k_prolate']} | {r['shift']} |")
    if g2:
        lv = [float(r["shift"]) for r in g2 if r["level"] > 1]
        l1 = [float(r["shift"]) for r in g2 if r["level"] == 1]
        raw1 = [float(r["n_star_minus_k_measured"]) for r in g2 if r["level"] == 1]
        print(f"\nshift (measured - Hankel prolate): level N > 1 in [{min(lv):.2f}, {max(lv):.2f}] (mean {sum(lv)/len(lv):.2f}); "
              f"level 1 in [{min(l1):.2f}, {max(l1):.2f}] (mean {sum(l1)/len(l1):.2f}); all |shift| <= {max(abs(v) for v in lv + l1):.2f}")
        for m in (0, 1):
            for lab, sel in (("N > 1", lambda r: r["level"] > 1), ("1", lambda r: r["level"] == 1)):
                vs = [float(v) for r in g2 if r["m"] == m and sel(r) for v in r["R_exact"]]
                print(f"GL(2) R_ex, level {lab}, m = {m}: [{min(vs):.3g}, {max(vs):.3g}]")
        if d1:
            vs = [float(v) for r in d1 for v in r["R_exact"]]
            print(f"degree-1 R_ex (all objects incl. zeta): [{min(vs):.3g}, {max(vs):.3g}]")


def gap_table(rows):
    if not rows:
        return
    print("\n### Root-number rule: ln(lambda_odd/lambda_even) at v = 10 against eps * ln(l_{nu,1}/l_{nu,0})\n")
    print("| object | N | k | eps | measured | Hankel-prolate prediction |")
    print("|---|---|---|---|---|---|")
    for r in rows:
        print(f"| {r['object']} | {r['level']} | {r['weight']} | {r['eps']:+d} | {r['measured_ln_odd_over_even']} | {r['prolate_predicted']} |")
    ok = sum(1 for r in rows if (float(r["measured_ln_odd_over_even"]) > 0) == (r["eps"] > 0))
    d = [abs(float(r["measured_ln_odd_over_even"]) - float(r["prolate_predicted"])) for r in rows]
    print(f"\nsign = eps: {ok}/{len(rows)}; |measured - predicted| in [{min(d):.2f}, {max(d):.2f}]")


def hankel_table(rows, label):
    if not rows:
        return
    print(f"\n### {label}\n")
    print("| nu | k | m | n*-k (exact Hankel prolate) | asymptotic |")
    print("|---|---|---|---|---|")
    for r in rows:
        print(f"| {r['nu']} | {r['k']} | {r['m']} | {r['n_star_minus_k']} | {r['asymptotic_n_minus_k']} |")


def main():
    bound_table(load("hermite_bound.json"), "Theorem 3 (Hermite trial functions; rate e^{-c})", 1)
    bound_table(load("kb_bound.json") + load("kb_bound_b.json"), "Theorem 4 (Kaiser-Bessel trial functions; rate e^{-2c})", 2)
    trial_table(load("prolate_trial_a.json") + load("prolate_trial_b.json") + load("prolate_trial_c.json"))
    leakage_tables(load("exact_leakage.json"))
    gap_table(load("root_number_gap.json"))
    hankel_table(load("hankel_index.json"), "Hankel-prolate continuous index on v = 6..10")
    hankel_table(load("hankel_index_v14_18.json"), "Hankel-prolate continuous index on v = 14..18 (prediction grid)")


if __name__ == "__main__":
    main()
