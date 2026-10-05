#!/usr/bin/env python3
"""Regenerate every number quoted in docs/DECAY_LAW.md, docs/ANTIPERIODIC.md and the reviewer notes of
docs/LIT_TANGENTS.md from the committed data, and compare it with the stated value.

This is the audit manifest: claim IDs match docs/BRIEF_AUDIT_DECAY_LAW.md. It checks transcription and
arithmetic only (data → doc). Whether the data are right is the auditor's job.

Usage:  .venv/bin/python scripts/decay_claims_check.py        (exit code 1 if any claim fails)
"""

import glob
import json
import math
import os
import re
import sys

import mpmath as mp

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import decay_law_mp as dl  # noqa: E402

PI4 = 4 * math.pi
RESULTS = []


def claim(cid, text, stated, got, tol, rel=False):
    """Record a claim. stated/got are floats or tuples of floats; tol is absolute (or relative if rel)."""
    s = stated if isinstance(stated, tuple) else (stated,)
    g = got if isinstance(got, tuple) else (got,)
    ok = len(s) == len(g) and all((abs(a - b) <= tol * abs(a)) if rel else (abs(a - b) <= tol) for a, b in zip(s, g))
    RESULTS.append((cid, ok, text, s, g))


def fmt(v):
    return ", ".join(f"{x:.6g}" for x in v)


def load(path):
    return json.load(open(path))


def scanhi(D):
    return load(f"data/connes/decay/scanhi_D{D}.json")


def shifted(rows, par, xq):
    r = next(r for r in rows if r["parity"] == par and float(r["xq"]) == xq)
    return math.log(float(r["lambda_n2"])) + PI4 * xq


def ell(n, c):
    return 0.5 * 4 * math.sqrt(math.pi) * 8 ** n * c ** (n + 0.5) * math.exp(-2 * c) / math.factorial(n)


def rel_diffs(path, pattern):
    return [float(m) for m in re.findall(pattern + r".*?rel diff ([0-9.e+-]+)", open(path).read())]


# ---------------------------------------------------------------- ANTIPERIODIC.md
def antiperiodic():
    txt = open("data/connes/antiperiodic/check.run.txt").read()
    d0 = [float(v) for v in re.findall(r"shift 0 vs build_form.*?max \|Δ\| = ([0-9.e+-]+)", txt)]
    claim("A1", "shift 0 reproduces build_form exactly (max |Δ|, even and odd)", (0.0, 0.0), tuple(d0), 0)
    anti_even = [float(v) for v in re.findall(r"antiperiodic even ω.*?rel diff ([0-9.e+-]+)", txt)]
    per_odd = [float(v) for v in re.findall(r"periodic     odd  ω.*?rel diff ([0-9.e+-]+)", txt)]
    off_odd = [float(v) for v in re.findall(r"antiperiodic odd  c = .*?rel diff ([0-9.e+-]+)", txt)]
    off_even = [float(v) for v in re.findall(r"antiperiodic even c = .*?rel diff ([0-9.e+-]+)", txt)]
    claim("A2a", "antiperiodic even diagonal checks: range of rel diff", (4.5e-7, 4.8e-5), (min(anti_even), max(anti_even)), 0.1, rel=True)
    claim("A2b", "antiperiodic odd b0+b1 off-diagonal checks (x = 5, 13)", (1.6e-8, 5.3e-9), tuple(off_odd), 0.1, rel=True)
    claim("A2c", "antiperiodic even b1−2b2+b3 checks (x = 5, 13)", (2.5e-6, 1.3e-5), tuple(off_even), 0.1, rel=True)
    claim("A2d", "periodic odd controls: range of rel diff", (1.0e-6, 2.0e-5), (min(per_odd), max(per_odd)), 0.1, rel=True)

    stated = {  # ratio λ_anti/λ_per at N = 16, 32, 48, 64, 100, 140, then anti48/per140
        ("0.8", "even"): (0.91, 1.03, 0.94, 0.96, 1.01, 1.04, 1.06), ("0.8", "odd"): (1.00, 1.08, 1.12, 1.06, 0.98, 0.96, 1.16),
        ("1.0", "even"): (1.34, 0.98, 1.03, 0.96, 1.00, 0.95, 1.20), ("1.0", "odd"): (1.23, 0.99, 1.04, 1.01, 1.01, 1.06, 1.31),
        ("1.19", "even"): (0.34, 0.66, 1.10, 1.03, 0.95, 0.98, 2.28), ("1.19", "odd"): (0.40, 0.76, 1.15, 1.03, 0.99, 1.02, 2.35),
        ("1.3", "even"): (0.32, 0.70, 0.62, 0.79, 0.99, 1.00, 1.6e3), ("1.3", "odd"): (0.32, 0.78, 0.63, 0.82, 1.03, 0.98, 8.5e2),
        ("1.495", "even"): (0.25, 0.19, 0.81, 0.21, 0.69, 1.00, 1.8e21), ("1.495", "odd"): (0.29, 0.21, 0.87, 0.21, 0.72, 1.03, 4.3e20)}
    for (a, par), st in stated.items():
        rows = load(f"data/connes/antiperiodic/sweep_a{a}.json")
        t = {(r["basis"], r["n"]): float(r["lambda_min"]) for r in rows if r["parity"] == par}
        got = tuple(t[("antiperiodic", n)] / t[("periodic", n)] for n in (16, 32, 48, 64, 100, 140)) + (t[("antiperiodic", 48)] / t[("periodic", 140)],)
        claim(f"A3[{a},{par}]", f"ratio table row a = {a} {par}", st, got, 0.051, rel=True)
    rows = load("data/connes/antiperiodic/sweep_a1.495.json")
    t = {r["n"]: float(r["lambda_min"]) for r in rows if r["parity"] == "even" and r["basis"] == "periodic"}
    dec = (math.log10(t[48]) - math.log10(t[64])) / 16
    claim("A4", "a = 1.495 even, periodic: decades per mode between N = 48 and 64", 0.58, dec, 0.01)
    claim("A4b", "factor 4.8 at N = 64 (anti vs periodic) is worth ≈ 1.2 modes", 1.2, math.log10(4.8) / dec, 0.05)
    claim("A5", "sweep size: 5 a × 6 N × 2 sectors × 2 bases runs", 120,
          sum(len(load(f)) for f in glob.glob("data/connes/antiperiodic/sweep_a*.json")), 0)


# ---------------------------------------------------------------- DECAY_LAW.md: ζ
def zeta_fit():
    mp.mp.dps = 30
    pts = dl.zeta_points()
    for par, (al, ga, be, rs) in {"even": (12.5556, 4.99, 15.55, 0.035), "odd": (12.549, 7.20, 18.85, 0.040)}.items():
        P = pts[par]
        c, r = dl._lsq([[-x, mp.log(x), 1] for x, _, _ in P], [mp.log(l) for _, l, _ in P])
        claim(f"Z1[{par}]", f"ζ fit {par}: α, γ, β, max resid", (al, ga, be, rs), (float(c[0]), float(c[1]), float(c[2]), float(r)), 0.006)
        claim(f"Z1b[{par}]", f"ζ fit {par}: α/4π", {"even": 0.99915, "odd": 0.99859}[par], float(c[0]) / PI4, 0.0001)
    lr = [mp.log(o[1] / e[1]) for e, o in zip(pts["even"], pts["odd"])]
    c, _ = dl._lsq([[mp.log(e[0]), 1] for e in pts["even"]], lr)
    claim("Z2", "λ_odd/λ_even ∝ x^p", 2.29, float(c[0]), 0.006)


# ---------------------------------------------------------------- DECAY_LAW.md: validation of L(χ) forms
def validation():
    odd, even = [], []
    for f in glob.glob("data/connes/decay/checks/check_D*.run.txt"):
        odd += rel_diffs(f, r"odd (?:k|c)")
        even += rel_diffs(f, r"even c")
    claim("V1", "L(χ) odd-sector zero checks: count, min, max rel diff (6 characters × 3)", (18, 3.9e-6, 5.0e-4), (len(odd), min(odd), max(odd)), 0.1, rel=True)
    claim("V2", "L(χ) even-sector zero checks: count, min, max rel diff (6 characters × 2)", (12, 1.4e-9, 4.3e-7), (len(even), min(even), max(even)), 0.1, rel=True)


# ---------------------------------------------------------------- DECAY_LAW.md: conductor scan
def conductor():
    mp.mp.dps = 30
    conv = []
    for D in (1, 5, 8, -3, -4, -7, -20):
        for r in scanhi(D):
            conv.append(abs(math.log(float(r["lambda_n1"]) / float(r["lambda_n2"]))))
    claim("C0", "convergence |ln(λ_5k*/λ_9k*)| max over all scan points", 0.17, max(conv), 0.005)
    table = {  # D: (even linear, even log, γ, odd linear, odd log, γ) on x/q ∈ [3, 8]
        1: (0.920, 1.002, 5.32, 0.885, 1.017, 8.57), 5: (0.994, 0.994, -0.02, 0.962, 0.987, 1.66),
        8: (0.989, 0.999, 0.61, 0.958, 0.994, 2.32), -3: (0.974, 1.004, 1.93, 0.940, 1.000, 3.90),
        -4: (0.975, 1.003, 1.82, 0.943, 0.998, 3.55), -7: (0.968, 0.987, 1.26, 0.940, 0.984, 2.83),
        -20: (0.975, 0.996, 1.36, 0.941, 0.996, 3.57)}
    for D, st in table.items():
        got = []
        for par in ("even", "odd"):
            R = sorted((r for r in scanhi(D) if r["parity"] == par and 3 <= float(r["xq"]) <= 8), key=lambda r: float(r["xq"]))
            xq = [mp.mpf(r["xq"]) for r in R]
            y = [mp.log(mp.mpf(r["lambda_n2"])) for r in R]
            c1, _ = dl._lsq([[u, 1] for u in xq], y)
            c2, _ = dl._lsq([[u, mp.log(u), 1] for u in xq], y)
            got += [float(-c1[0] / (4 * mp.pi)), float(-c2[0] / (4 * mp.pi)), float(c2[1])]
        claim(f"C1[{D}]", f"P-DL1 table row D = {D}", st, tuple(got), 0.0051)
    classes = {"even χ": (5, 8), "odd χ": (-3, -4, -7, -20)}
    stated = {("even χ", "even"): ((5.3, 5.7), (5.8, 6.6), (6.3, 6.9)), ("even χ", "odd"): ((12.2, 12.4), (14.7, 14.9), (16.8, 16.9)),
              ("odd χ", "even"): ((8.3, 9.3), (9.8, 11.2), (11.0, 12.2)), ("odd χ", "odd"): ((15.1, 15.9), (18.9, 19.3), (21.5, 22.1))}
    for (cls, par), st in stated.items():
        for (lo, hi), xq in zip(st, (2.0, 5.0, 10.0)):
            v = [shifted(scanhi(D), par, xq) for D in classes[cls]]
            claim(f"C2[{cls},{par},{xq:g}]", f"class band ln λ + 4π x/q, {cls} {par} at x/q = {xq:g}", (lo, hi), (min(v), max(v)), 0.051)
    for par, st in {"even": (18.5, 23.7, 27.1), "odd": (22.6, 30.7, 35.6)}.items():
        claim(f"C2[ζ,{par}]", f"ζ ln λ + 4π x at x = 2, 5, 10 ({par})", st, tuple(shifted(scanhi(1), par, xq) for xq in (2.0, 5.0, 10.0)), 0.051)
    gap = [shifted(scanhi(1), "even", xq) - shifted(scanhi(D), "even", xq) for D in (5, 8) for xq in (2.0, 10.0)]
    claim("C3", "ζ above even χ at equal x/q: e^{13} to e^{21} (range of the ln gap at x/q = 2 and 10)", (12.9, 20.8), (min(gap), max(gap)), 0.06)


def prolate_index():
    assigned = {5: (0, 2), 8: (0, 2), -3: (1, 3), -4: (1, 3), -7: (1, 3), -20: (1, 3), 1: (4, 6)}
    stated = {5: ((23, 35), (3.8, 6.2)), 8: ((15, 20), (3.0, 5.4)), -3: ((8.6, 14), (2.7, 5.6)), -4: ((6.7, 9.8), (2.7, 4.6)),
              -7: ((3.1, 4.6), (1.7, 3.7)), -20: ((3.3, 4.1), (3.9, 6.4)), 1: ((2.1, 7.9), (0.37, 4.6))}
    flat_stated = {5: (0, 2), 8: (0, 2), -3: (1, 3), -4: (1, 3), -7: (1, 3), -20: (1, 3), 1: (5, 7)}
    for D in assigned:
        got_r, got_f = [], []
        for par, n in zip(("even", "odd"), assigned[D]):
            R = sorted((r for r in scanhi(D) if r["parity"] == par), key=lambda r: float(r["xq"]))
            rat = lambda m: [float(r["lambda_n2"]) / ell(m, 2 * math.pi * float(r["xq"])) for r in R]
            v = rat(n)
            got_r += [v[0], v[-1]]
            got_f.append(min(range(8), key=lambda m: abs(math.log(rat(m)[-1] / rat(m)[1]))))
        st = stated[D]
        claim(f"P1[{D}]", f"λ/ℓ_n first → last (x/q = 2 → 10), D = {D}", (st[0][0], st[0][1], st[1][0], st[1][1]), tuple(got_r), 0.051, rel=True)
        claim(f"P2[{D}]", f"flattest prolate index (even, odd), D = {D}", flat_stated[D], tuple(got_f), 0)


def dh():
    rows = []
    for f in glob.glob("data/connes/decay/scanhi_dh_*.json"):
        rows += load(f)
    t = {(r["parity"], float(r["xq"])): float(r["lambda_n2"]) for r in rows}
    stated = {("even", 2): 8.58e-8, ("even", 3): 4.42e-13, ("even", 4): 4.21e-18, ("even", 5): 1.42e-23, ("even", 6): 9.39e-29,
              ("even", 6.25): -8.33e-30, ("even", 6.5): -1.90e-23, ("odd", 2): 6.25e-5, ("odd", 3): 1.23e-9, ("odd", 4): 1.35e-14,
              ("odd", 5): 1.02e-19, ("odd", 6): 6.21e-25, ("odd", 6.25): 4.56e-27, ("odd", 6.5): -2.98e-26}
    for k, v in stated.items():
        claim(f"D1[{k[0]},{k[1]:g}]", f"DH λ_min {k[0]} at x/q = {k[1]:g}", v, t[k], 0.006, rel=True)
    # DH inside the odd-χ class band (both sectors) for x/q ∈ {2, 3, 4, 5, 6}
    inside = []
    for par in ("even", "odd"):
        for xq in (2.0, 3.0, 4.0, 5.0, 6.0):
            band = [shifted(scanhi(D), par, xq) for D in (-3, -4, -7, -20)]
            v = math.log(t[(par, xq)]) + PI4 * xq
            inside.append(1.0 if min(band) <= v <= max(band) else 0.0)
    claim("D2", "DH inside the odd-χ band at all 10 (sector, x/q ≤ 6) points", 10.0, sum(inside), 0)


def tower():
    rows = load("data/connes/decay/tower.json")
    g = {(r["a"], r["parity"], r["n"]): r for r in rows}
    for par, st1, st2 in (("even", (1.03, 1.09, 1.28), (1.19, 2.48, 9.99)), ("odd", (1.00, 1.03, 1.07), (1.21, 1.54, 3.02))):
        claim(f"T1[{par}]", f"λ(ζ_K)/λ(L(χ−20)) at a = 1.1, 1.3, 1.495 ({par}, N = 96)", st1,
              tuple(float(g[(a, par, 96)]["ratio_zetaK_over_chi20"]) for a in ("1.1", "1.3", "1.495")), 0.0051)
        claim(f"T2[{par}]", f"λ(ζ_H)/λ(ζ_K) at a = 1.1, 1.3, 1.495 ({par}, N = 96)", st2,
              tuple(float(g[(a, par, 96)]["ratio_zetaH_over_zetaK"]) for a in ("1.1", "1.3", "1.495")), 0.0051)
    md = max(float(r["max_diff_vs_control"]) for r in rows)
    claim("T3", "ζ + L(χ−20) matches control_normalisation's ζ_K (max entry difference ≤ 1e-59)", 1.0, 1.0 if md <= 1e-59 else 0.0, 0)
    txt = open("data/connes/decay/tower_control.txt").read()
    line = next(l for l in txt.splitlines() if l.startswith("a=1.495 even"))
    vals = [float(v) for v in re.findall(r": ([0-9.]+)(?:;|$)", line.split("λ(ζ_K·extra)/λ(ζ_K):")[1])]
    claim("T4", "control boosts at e^{2.99} even: cover, (−7,13), (−4,−3), (5,8), (−3,8), χ−4 alone, χ5 alone",
          (9.9, 23.1, 9.5, 6.6, 2.95, 7.7, 5.1), (vals[0], vals[4], vals[2], vals[3], vals[1], vals[5], vals[6]), 0.051)


def pdl2():
    d = load("data/connes/decay/pdl2.json")
    rels = [float(v) for r in d["rows"] for v in r["R_over_R10"].values()]
    claim("Q1", "P-DL2: R/R(10) range over x/q = 12, 14, 16, all 12 cases", (0.79, 1.20), (min(rels), max(rels)), 0.006)
    claim("Q2", "P-DL2: flattest index (3 vs 16) ≠ assigned, out of 12", 0.0, float(d["flattest_mismatch"]), 0)
    conv = [abs(math.log(float(r["lambda_n1"]) / float(r["lambda_n2"]))) for f in glob.glob("data/connes/decay/scanext_D*.json") for r in load(f)]
    claim("Q3", "P-DL2 extension convergence: point count, max |ln(λ_5k*/λ_9k*)|", (36, 0.144), (len(conv), max(conv)), 0.0011)
    prec = [float(r["rel_change"]) for f in glob.glob("data/connes/decay/precision/precext_*.json") for r in load(f)]
    claim("Q4", "P-DL2 extension: points precision-rechecked at +40 digits, and all stable to < 1e-13", (36, 1.0), (len(prec), 1.0 if max(prec) < 1e-13 else 0.0), 0)


def products():
    for t, (lo2, hi8) in {"1_-20": (1.54, 1.0e18), "-3_-4": (285, 6.9e24), "5_8": (2.2e3, 2.0e28)}.items():
        rows = load(f"data/connes/decay/product_{t}.json")
        r2 = min(float(r["rho"]) for r in rows if float(r["xq"]) == 2)
        r8 = min(float(r["rho"]) for r in rows if float(r["xq"]) == 8)
        claim(f"R1[{t}]", f"P-DL3 {t}: min ρ over sectors at x/q_max = 2 and at 8", (lo2, hi8), (r2, r8), 0.051, rel=True)
    mp.mp.dps = 30
    stated = {("1_-20", "even"): 2.080, ("1_-20", "odd"): 2.288, ("-3_-4", "even"): 2.182, ("-3_-4", "odd"): 1.788,
              ("5_8", "even"): 2.040, ("5_8", "odd"): 2.052}
    for (t, par), st in stated.items():
        R = sorted((r for r in load(f"data/connes/decay/productv_{t}.json") if r["parity"] == par), key=lambda r: float(r["v"]))
        v = [mp.mpf(r["v"]) for r in R]
        c, _ = dl._lsq([[-u, mp.log(u), 1] for u in v], [mp.log(mp.mpf(r["lambda_n2"])) for r in R])
        claim(f"R2[{t},{par}]", f"P-DL5 log-corrected α/4π, {t} {par}", st, float(c[0] / (4 * mp.pi)), 0.0006)
    conv = [abs(math.log(float(r["lambda_n1"]) / float(r["lambda_n2"]))) for f in glob.glob("data/connes/decay/productv_*.json") for r in load(f)]
    claim("R3", "P-DL5 convergence max |Δ ln λ| (≤ 0.20 stated, ≤ 0.3 required)", 0.20, max(conv), 0.005)


def gl2():
    mp.mp.dps = 30
    import gl2_form_mp  # noqa: F401  (ensures the GL(2) code imports)
    for lab, st in {"11a1": {"even": (0.2105, 2.006), "odd": (0.2010, 1.916)}, "37b1": {"even": (0.2048, 1.953), "odd": (0.2089, 1.991)}}.items():
        rows = load(f"data/connes/gl2/scan_{lab}.json")
        for par, (a_st, b_st) in st.items():
            S = sorted((r for r in rows if r["parity"] == par), key=lambda r: float(r["xq"]))
            y = [mp.log(mp.mpf(r["lambda_n2"])) for r in S]
            got = []
            for v in ([mp.mpf(r["xq"]) for r in S], [mp.sqrt(mp.mpf(r["xq"])) for r in S]):
                c, _ = dl._lsq([[-u, mp.log(u), 1] for u in v], y)
                got.append(float(c[0] / (4 * mp.pi)))
            claim(f"G1[{lab},{par}]", f"P-DL4 {lab} {par}: α/4π for A′ (x/N) and B′ (√(x/N))", (a_st, b_st), tuple(got), 0.0006)
    st6 = {("15a1", "even"): 2.138, ("15a1", "odd"): 2.161, ("19a1", "even"): 2.134, ("19a1", "odd"): 2.157,
           ("delta", "even"): 2.102, ("delta", "odd"): 2.169}
    for (obj, par), st in st6.items():
        S = sorted((r for r in load(f"data/connes/gl2/scan6_{obj}.json") if r["parity"] == par), key=lambda r: float(r["v"]))
        v = [mp.mpf(r["v"]) for r in S]
        c, _ = dl._lsq([[-u, mp.log(u), 1] for u in v], [mp.log(mp.mpf(r["lambda_n2"])) for r in S])
        claim(f"G2[{obj},{par}]", f"P-DL6a {obj} {par}: log-corrected α/4π", st, float(c[0] / (4 * mp.pi)), 0.0006)
    flat = {("15a1", "even"): 1, ("15a1", "odd"): 3, ("19a1", "even"): 1, ("19a1", "odd"): 4, ("delta", "even"): 14, ("delta", "odd"): 16}
    for (obj, par), n_st in flat.items():
        S = sorted((r for r in load(f"data/connes/gl2/scan6_{obj}.json") if r["parity"] == par), key=lambda r: float(r["v"]))
        R = lambda m, r: float(r["lambda_n2"]) / ell(m, 4 * math.pi * float(r["v"]))
        got = min(range(21), key=lambda m: abs(math.log(R(m, S[-1]) / R(m, S[1]))))
        claim(f"G3[{obj},{par}]", f"P-DL6b {obj} {par}: flattest n (c = 4πv, metric v5 vs v2)", float(n_st), float(got), 0)


def pdl8():
    d = load("data/connes/decay/pdl8.json")
    a = [float(r["alpha_over_4pi"]) for r in d["rows"]]
    claim("E1", "P-DL8: index hits out of 12, and α/4π range", (11, 0.984, 1.026), (float(d["index_hits"]), min(a), max(a)), 0.0006)
    gates = [float(m) for f in glob.glob("data/connes/decay/checks/check8_D*.run.txt")
             for m in re.findall(r"GATE PASS: 5 tests, max absolute error ([0-9.e+-]+)", open(f).read())]
    claim("E2", "P-DL8 gates: forms passed, and worst absolute error", (6, 9.9e-5), (float(len(gates)), max(gates)), 0.011, rel=True)
    conv = [abs(math.log(float(r["lambda_n1"]) / float(r["lambda_n2"]))) for f in glob.glob("data/connes/decay/scan8_D*.json") for r in load(f)]
    claim("E3", "P-DL8 convergence max |Δ ln λ|", 0.156, max(conv), 0.0006)


def pdl9_10():
    mp.mp.dps = 30
    eps = {"delta": 1, "f16": 1, "f18": -1, "f20": 1, "f22": -1, "f26": -1}
    al, mism = [], 0
    for o, e in eps.items():
        R = load(f"data/connes/gl2/scan9_{o}.json")
        R = R if isinstance(R, list) else R["rows"]
        lam = {(float(r["v"]), r["parity"]): mp.mpf(r["lambda_n2"]) for r in R}
        for par in ("even", "odd"):
            vs = sorted(v for v, p in lam if p == par)
            c, _ = dl._lsq([[-mp.mpf(v), mp.log(v), 1] for v in vs], [mp.log(lam[(v, par)]) for v in vs])
            al.append(float(c[0] / (4 * mp.pi)))
        mism += sum((1 if lam[(v, "odd")] > lam[(v, "even")] else -1) != e for v in sorted({v for v, _ in lam}))
    claim("H1", "P-DL9 level-1: α/4π range over 12 cases, and ε-rule mismatches over 30 points", (1.952, 2.074, 0.0), (min(al), max(al), float(mism)), 0.0011)
    R = load("data/connes/gl2/scan9_37a1.json")
    R = R if isinstance(R, list) else R["rows"]
    lam = {(float(r["v"]), r["parity"]): mp.mpf(r["lambda_n2"]) for r in R}
    a37 = []
    for par in ("even", "odd"):
        vs = sorted(v for v, p in lam if p == par)
        c, _ = dl._lsq([[-mp.mpf(v), mp.log(v), 1] for v in vs], [mp.log(lam[(v, par)]) for v in vs])
        a37.append(float(c[0] / (4 * mp.pi)))
    neg = sum(lam[(v, "odd")] < lam[(v, "even")] for v in sorted({v for v, _ in lam}))
    claim("H1b", "P-DL9 37a1: α/4π (even, odd) and points with λ_odd < λ_even (of 5)", (1.889, 1.850, 5.0), (a37[0], a37[1], float(neg)), 0.0011)
    al3 = []
    for t in ("-3_-4_5", "-4_5_8", "1_-3_-4"):
        R = load(f"data/connes/decay/productv3_{t}.json")
        for par in ("even", "odd"):
            S = sorted((r for r in R if r["parity"] == par), key=lambda r: float(r["v"]))
            c, _ = dl._lsq([[-mp.mpf(r["v"]), mp.log(mp.mpf(r["v"])), 1] for r in S], [mp.log(mp.mpf(r["lambda_n2"])) for r in S])
            al3.append(float(c[0] / (4 * mp.pi)))
    claim("H2", "P-DL10 degree 3: log-corrected α/4π range over 6 cases", (2.743, 3.203), (min(al3), max(al3)), 0.0011)


def pdl11():
    d = load("data/connes/gl2/pdl11_fit.json")
    al = [float(r["alpha_over_4pi"]) for r in d["11a"]]
    claim("K1", "P-DL11a α/4π range over 10 cases", (1.950, 2.062), (min(al), max(al)), 0.0011)
    off = {(c["obj"], c["parity"]): float(c["n_minus_k"]) for c in d["11b"]}
    ev = [off[(o, "even")] for o in ("g4", "g6", "g8", "11a1")]
    od = [off[(o, "odd")] for o in ("g4", "g6", "g8", "11a1")]
    claim("K2", "P-DL11 post hoc: ε = +1 level-N offsets, even range and odd range", (-0.71, -0.59, 1.62, 1.88), (min(ev), max(ev), min(od), max(od)), 0.006)
    claim("K4", "P-DL11 post hoc: 37a1 offsets (even, odd)", (1.32, -0.67), (off[("37a1", "even")], off[("37a1", "odd")]), 0.006)
    claim("K3", "P-DL11b within ±1.2 (of 10), off by > 2", (0.0, 1.0), (float(d["11b_within"]), float(d["11b_off_by_more_than_2"])), 0)


def pdl12():
    d = load("data/connes/gl2/pdl12_fit.json")
    within = sum(1 for r in d["12a"] if abs(float(r["n_minus_k"]) - float(r["hypothesis"])) <= 0.5)
    al = [float(r["alpha_over_4pi"]) for r in d["12b"]]
    claim("P12", "P-DL12: index cases within ±0.5 (of 8); α/4π range", (8.0, 1.881, 2.06), (float(within), min(al), max(al)), 0.0011)
    conv = [abs(math.log(float(r["lambda_n1"]) / float(r["lambda_n2"]))) for f in glob.glob("data/connes/gl2/scan12_*.json")
            for r in (lambda x: x if isinstance(x, list) else x["rows"])(load(f))]
    claim("P12b", "P-DL12 convergence: points, max |Δ ln λ|", (40.0, 0.067), (float(len(conv)), max(conv)), 0.0011)


def square_torus():
    mp.mp.dps = 30
    rows = load("data/connes/decay/square_torus.json")
    zs = [float(r["Zsq_n96"]) for r in rows]
    z1e = {round(float(r["x"]), 2): float(r["Z1_n96"]) for r in rows if r["parity"] == "even"}
    z1o = {round(float(r["x"]), 2): float(r["Z1_n96"]) for r in rows if r["parity"] == "odd"}
    claim("ST1", "P-ST1: Z_sq min over 14 points > 0; Z1 even sign flips between 13.46 and 19.89; odd between 19.89 and 24.53",
          (1.0, 1.0, 1.0), (1.0 if min(zs) > 0 else 0.0, 1.0 if (z1e[13.46] > 0 > z1e[19.89]) else 0.0,
                            1.0 if (z1o[19.89] > 0 > z1o[24.53]) else 0.0), 0)
    a2s = sorted({r["a2"] for r in rows}, key=float)
    by = {(r["a2"], r["parity"]): mp.mpf(r["Zsq_n96"]) for r in rows}
    got = []
    for par in ("even", "odd"):
        v = [mp.sqrt(mp.exp(mp.mpf(a)) / 4) for a in a2s]
        c, _ = dl._lsq([[-u, mp.log(u), 1] for u in v], [mp.log(by[(a, par)]) for a in a2s])
        got.append(float(c[0] / (4 * mp.pi)))
    claim("ST2", "P-ST2: Z_sq log-corrected α/4π (even, odd) in v = √(x/4)", (2.061, 2.133), tuple(got), 0.0011)
    claim("ST3", "P-ST3: Z_sq odd > even at every support (count of 7)", 7.0, float(sum(by[(a, "odd")] > by[(a, "even")] for a in a2s)), 0)
    conv = [abs(math.log(float(r["Zsq_n64"]) / float(r["Zsq_n96"]))) for r in rows]
    claim("ST4", "P-ST Z_sq convergence max |ln(λ64/λ96)|", 0.074, max(conv), 0.0011)


def r2():
    vals = []
    for f in ("hermite_bound", "hermite_matrix_check", "kb_bound", "kb_bound_b"):
        d = load(f"data/research_r2/{f}.json")
        rows = d if isinstance(d, list) else next(v for v in d.values() if isinstance(v, list) and v and isinstance(v[0], dict))
        vals += [float(r["ratio_measured_over_bound"]) for r in rows if r.get("ratio_measured_over_bound") not in (None, "")]
    claim("R2a", "R2: points tabulated and max Ritz/bound ratio (consistency only, NOT a validation: both are upper bounds)", (334.0, 0.0293), (float(len(vals)), max(vals)), 0.0006)
    cb = load("data/research_r2/bounds_at_certified_supports.json")["rows"]
    claim("R2b", "R2 validation: bounds above the certified lower bound at the certified supports (count, violations)", (16.0, 0.0),
          (float(len(cb)), float(sum(1 for r in cb if float(r["bound"]) < float(r["certified_lower_bound"])))), 0)


def r389():
    mp.mp.dps = 30
    rows = load("data/connes/gl2/scanR389.json")
    lam = {(float(r["v"]), r["parity"]): mp.mpf(r["lambda_n2"]) for r in rows}
    vs = sorted({v for v, _ in lam})
    lr = [float(mp.log(lam[(v, "odd")] / lam[(v, "even")])) for v in vs]
    claim("R389a", "P-R389: points with odd below even (of 5), and ln(λ_odd/λ_even) range", (5.0, -8.755, -7.381), (float(sum(x < 0 for x in lr)), min(lr), max(lr)), 0.002)
    al = []
    for par in ("even", "odd"):
        c, _ = dl._lsq([[-mp.mpf(v), mp.log(v), 1] for v in vs], [mp.log(lam[(v, par)]) for v in vs])
        al.append(float(c[0] / (4 * mp.pi)))
    claim("R389b", "P-R389 rate: log-corrected α/4π (even, odd)", (1.953, 1.937), tuple(al), 0.0011)


def r2_rigorous():
    def rows(f):
        d = load(f"data/research_r2_rigorous/{f}.json")
        return d if isinstance(d, list) else next(v for v in d.values() if isinstance(v, list) and v and isinstance(v[0], dict))
    def key(r):
        return (r["construction"].lower()[:2], int(r["D"]), r["parity"], round(float(r.get("support_2a") or 0), 4), round(float(r.get("xq") or -1), 4))
    bounds = {key(r): float(r["bound_upper"]) for f in ("hermite_cert", "hermite_cert_supports", "kb_cert", "kb_cert_supports") for r in rows(f)}
    rq = [(float(r.get("rq_hi") or r["rq"]), bounds[key(r)]) for f in ("rq_check_grid", "rq_check_supports") for r in rows(f)]
    claim("R2C1", "R2 Theorem C validation (a): trial RQ ≤ B_cert (cases, failures)", (386.0, 0.0), (float(len(rq)), float(sum(a > b for a, b in rq))), 0)
    cert = {(1.6, "even"): 1.0277e-17, (1.6, "odd"): 9.118e-15, (2.38, "even"): 6.81311646951e-48, (2.38, "odd"): 3.9196686160e-44,
            (2.6, "even"): 5.76344479222e-62, (2.6, "odd"): 5.36011304048e-58, (2.99, "even"): 3.50114217641e-96, (2.99, "odd"): 8.25626494001e-92}
    sup = [(float(r["bound_upper"]), cert[(round(float(r["support_2a"]), 2), r["parity"])]) for f in ("hermite_cert_supports", "kb_cert_supports") for r in rows(f)]
    claim("R2C2", "R2 Theorem C validation (b): B_cert ≥ certified value at the four supports (cases, failures)", (16.0, 0.0), (float(len(sup)), float(sum(b < c for b, c in sup))), 0)


def lit_notes():
    log = open("data/connes/zhu_window_a1495_n180.run.txt").read()
    e180 = float(re.search(r"N=180 even: eps0=([0-9.e+-]+)", log).group(1))
    o180 = float(re.search(r"N=180 odd: eps0=([0-9.e+-]+)", log).group(1))
    e140 = float(re.search(r"N=140 even: eps0=([0-9.e+-]+)", log).group(1))
    claim("L1", "ζ at e^{2.99}: even N = 180, N = 140, odd N = 180", (6.27e-96, 8.27e-96, 1.48e-91), (e180, e140, o180), 0.005, rel=True)
    claim("L2", "even N = 180 upper bound is below the literal leading term of Zhu's form, 2.26e-95 (0.59% in −ln λ; not an exclusion)", 1.0, 1.0 if e180 < 2.26e-95 else 0.0, 0)
    claim("L2b", "−ln λ gap to Zhu's leading term at e^{2.99}, in % (review round 1)", 0.59,
          100 * (math.log(e180) / math.log(2.26e-95) - 1), 0.01)
    claim("L3", "odd N = 180 is inside the kill window [3, 25]e-92 and above the bracket [5.7, 13]e-92", 1.0,
          1.0 if (3e-92 <= o180 <= 25e-92 and o180 > 13e-92) else 0.0, 0)


DOC_STRINGS = [  # the manifest's stated values must also be what the docs say
    ("docs/DECAY_LAW.md", "α = 12.5556 (α/4π = 0.99915), γ = 4.99, β = 15.55, max residual 0.035"),
    ("docs/DECAY_LAW.md", "α = 12.549 (α/4π = 0.99859), γ = 7.20, β = 18.85, max residual 0.040"),
    ("docs/DECAY_LAW.md", "x^{2.29}"),
    ("docs/DECAY_LAW.md", "relative error 3.9e-6 to 5.0e-4"),
    ("docs/DECAY_LAW.md", "relative error 1.4e-9 to 4.3e-7"),
    ("docs/DECAY_LAW.md", "|ln(λ_5k*/λ_9k*)| ≤ 0.17"),
    ("docs/DECAY_LAW.md", "| 1 (ζ) | 1 | even + pole | 0.920 | 1.002 | 5.32 | 0.885 | 1.017 | 8.57 |"),
    ("docs/DECAY_LAW.md", "| −20 | 20 | odd | 0.975 | 0.996 | 1.36 | 0.941 | 0.996 | 3.57 |"),
    ("docs/DECAY_LAW.md", "e^{13} (x/q = 2) to e^{21} (x/q = 10)"),
    ("docs/DECAY_LAW.md", "| 6.25 | 31.25 | **−8.33e-30** |"),
    ("docs/DECAY_LAW.md", "| 6.5 | 32.5 | **−1.90e-23** | — | | **−2.98e-26** |"),
    ("docs/ANTIPERIODIC.md", "all 120 runs"),
    ("docs/ANTIPERIODIC.md", "| 1.495 | even | 0.25 | 0.19 | 0.81 | 0.21 | 0.69 | 1.00 | 1.8e21 |"),
    ("docs/ANTIPERIODIC.md", "periodic odd controls: 1e-6 to 2e-5"),
    ("docs/LIT_TANGENTS.md", "The N = 180 upper bound is 1.48·10⁻⁹¹"),
    ("docs/LIT_TANGENTS.md", "**not** an exclusion of the conjecture"),
    ("docs/DECAY_LAW.md", "and R runs from 1.55 to 38.5"),
    ("docs/DECAY_LAW.md", "up to N = 838 (D = −20, x/q = 16)"),
    ("docs/DECAY_LAW.md", "λ falls from 1e-14–1.5e-9 to 3e-85–1e-77"),
    ("docs/DECAY_LAW.md", "Every R/R(10) lies in [0.79, 1.20]"),
    ("docs/DECAY_LAW.md", "sign(ln λ_odd/λ_even) = ε at all 35 grid points"),
    ("docs/LIT_TANGENTS.md", "λ(ζ_H)/λ(ζ_K) = 1.19, 2.48, **9.99** (even) and 1.21, 1.54, **3.02** (odd)"),
]


def docs_text():
    for path, sub in DOC_STRINGS:
        claim(f"DOC[{os.path.basename(path)}]", f"doc contains: {sub[:70]}", 1.0, 1.0 if sub in open(path).read() else 0.0, 0)


def main():
    os.chdir(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
    for fn in (antiperiodic, zeta_fit, validation, conductor, prolate_index, dh, tower, pdl2, products, gl2, pdl8, pdl9_10, pdl11, pdl12, square_torus, r2, r389, r2_rigorous, lit_notes, docs_text):
        fn()
    bad = 0
    for cid, ok, text, s, g in RESULTS:
        bad += not ok
        print(f"{'PASS' if ok else 'FAIL'}  {cid:<22} {text}\n      stated: {fmt(s)}\n      data:   {fmt(g)}")
    print(f"\n{len(RESULTS) - bad}/{len(RESULTS)} claims match the committed data")
    sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()
