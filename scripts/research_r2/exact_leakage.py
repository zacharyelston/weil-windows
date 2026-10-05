#!/usr/bin/env python3
"""Measured window minima against the exact prolate leakage instead of the Fuchs asymptotic.

Degree 1 (scanhi/scanext/scan8): R_F = lambda / l_n(c) (l_n = Fuchs/2, our convention) and
R_ex = lambda / ((1 - lambda_n(c))/2) with the exact 1D eigenvalue (hankel_prolate.leakage at nu = -1/2, +1/2), c = 2 pi x/q.
GL(2) (scan9, scan11, scan12): R_ex = lambda / ((1 - lambda_{nu,m}(c))/2) with nu = k - 1, c = 4 pi v, m = 0 in the sector
(-1)^s = epsilon and m = 1 in the other; the continuous index n* (P-DL11 metric) from the exact leakage is also reported
as the shift Delta = n*_measured - n*_prolate (0 if the exact generalised prolate accounts for the whole index).

Usage: .venv/bin/python scripts/research_r2/exact_leakage.py --json data/research_r2/exact_leakage.json
"""
import argparse
import glob
import json
import os
import sys

import mpmath as mp

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from hankel_prolate import leakage  # noqa: E402
from decay_law_mp import _ell  # noqa: E402


def lsq_slope(xs, ys):
    xb, yb = mp.fsum(xs) / len(xs), mp.fsum(ys) / len(ys)
    return mp.fsum((a - xb) * (b - yb) for a, b in zip(xs, ys)) / mp.fsum((a - xb) ** 2 for a in xs)


def degree1(rows_out):
    files = sorted(glob.glob("data/connes/decay/scanhi_D*.json") + glob.glob("data/connes/decay/scan8_D*.json"))
    for f in files:
        D = int(f.split("_D")[-1].split(".json")[0])
        rows = json.load(open(f))
        ext = f.replace("scanhi_", "scanext_")
        if "scanhi_" in f and os.path.exists(ext):
            rows += json.load(open(ext))
        kappa = 0 if D > 0 else 1
        pole = 1 if D == 1 else 0
        for s, par in enumerate(("even", "odd")):
            n = kappa + 2 * s + 4 * pole
            pts = []
            for r in rows:
                if r["parity"] != par:
                    continue
                xq = mp.mpf(r["xq"])
                c = 2 * mp.pi * xq
                mp.mp.dps = int(2 * float(c) / 2.3) + 40
                nu = mp.mpf(-1) / 2 if n % 2 == 0 else mp.mpf(1) / 2
                lk = leakage(nu, c, n // 2)
                lam = mp.mpf(r["lambda_n2"])
                pts.append((xq, c, lam, lam / _ell(n, c), lam / (lk / 2), lk))
            mp.mp.dps = 30
            pts.sort()
            sel = [p for p in pts if p[0] >= 3]
            xs = [mp.log(p[1]) for p in sel]
            nF = lsq_slope(xs, [mp.log(p[2]) + 2 * p[1] for p in sel]) - mp.mpf(1) / 2
            nE = lsq_slope(xs, [mp.log(p[2]) - mp.log(p[5]) for p in sel]) + n
            rows_out.append({"object": f"D={D}", "sector": par, "n_rule": n, "xq": [mp.nstr(p[0], 4) for p in pts],
                             "R_fuchs": [mp.nstr(p[3], 5) for p in pts], "R_exact": [mp.nstr(p[4], 5) for p in pts],
                             "R_fuchs_last_over_first_from3": mp.nstr(sel[-1][3] / sel[0][3], 4),
                             "R_exact_last_over_first_from3": mp.nstr(sel[-1][4] / sel[0][4], 4),
                             "n_star_fuchs": mp.nstr(nF, 4), "n_star_exact": mp.nstr(nE, 4)})
            print(f"D={D} {par} n={n}: R_F {mp.nstr(pts[0][3],4)}..{mp.nstr(pts[-1][3],4)} (x/q>=3 drift {mp.nstr(sel[-1][3]/sel[0][3],4)})  "
                  f"R_ex {mp.nstr(pts[0][4],4)}..{mp.nstr(pts[-1][4],4)} (drift {mp.nstr(sel[-1][4]/sel[0][4],4)})  n*_F={mp.nstr(nF,4)} n*_ex={mp.nstr(nE,4)}", flush=True)


EPS = {"delta": 1, "f16": 1, "f20": 1, "f18": -1, "f22": -1, "f26": -1, "g4": 1, "g6": 1, "g8": 1,
       "11a1": 1, "37a1": -1, "14a1": 1, "17a1": 1, "43a1": -1, "53a1": -1}


def gl2(rows_out):
    files = sorted(glob.glob("data/connes/gl2/scan9_*.json") + glob.glob("data/connes/gl2/scan11_*.json") + glob.glob("data/connes/gl2/scan12_*.json"))
    for f in files:
        rows = json.load(open(f))
        obj = rows[0]["obj"]
        if obj == "37a1" and "scan9_" in f:
            continue  # the v = 2..4 grid; the asymptotic grid is scan11
        eps = EPS[obj]
        k = int(rows[0]["weight"])
        nu = mp.mpf(k - 1)
        for s, par in enumerate(("even", "odd")):
            m = 0 if (1 if par == "even" else -1) == eps else 1
            pts = []
            for r in rows:
                if r["parity"] != par:
                    continue
                v = mp.mpf(r["v"])
                c = 4 * mp.pi * v
                mp.mp.dps = int(2 * float(c) / 2.3) + 40
                lk = leakage(nu, c, m)
                lam = mp.mpf(r["lambda_n2"])
                pts.append((v, c, lam, lam / (lk / 2), lk))
            mp.mp.dps = 30
            pts.sort()
            xs = [mp.log(p[1]) for p in pts]
            nmeas = lsq_slope(xs, [mp.log(p[2]) + 2 * p[1] for p in pts]) - mp.mpf(1) / 2
            nprol = lsq_slope(xs, [mp.log(p[4]) + 2 * p[1] for p in pts]) - mp.mpf(1) / 2
            rows_out.append({"object": obj, "level": rows[0]["N"], "weight": k, "eps": eps, "sector": par, "m": m,
                             "v": [mp.nstr(p[0], 4) for p in pts], "R_exact": [mp.nstr(p[3], 5) for p in pts],
                             "n_star_minus_k_measured": mp.nstr(nmeas - k, 4), "n_star_minus_k_prolate": mp.nstr(nprol - k, 4),
                             "shift": mp.nstr(nmeas - nprol, 4)})
            print(f"{obj} N={rows[0]['N']} k={k} eps={eps} {par} m={m}: n*-k measured {mp.nstr(nmeas-k,4)} prolate {mp.nstr(nprol-k,4)} "
                  f"shift {mp.nstr(nmeas-nprol,4)}; R_ex {mp.nstr(pts[0][3],4)}..{mp.nstr(pts[-1][3],4)}", flush=True)


CERT = {"even": {"4.95303242439511": "1.0277e-17", "10.8049028639313": "6.81311646951e-48", "13.4637380350017": "5.76344479222e-62",
                 "19.8856824915647": "3.50114217641e-96"},
        "odd": {"4.95303242439511": "9.118e-15", "10.8049028639313": "3.9196686160e-44", "13.4637380350017": "5.36011304048e-58",
                "19.8856824915647": "8.25626494001e-92"}}


def zeta_fit(rows_out):
    """zeta's best upper bounds (fit_zeta.json) and certified lower bounds (docs/CERTIFICATE_238.md) against exact leakage."""
    fz = json.load(open("data/connes/decay/fit_zeta.json"))
    for par, n in (("even", 4), ("odd", 6)):
        pts = []
        for xs, lam, src in fz[par]["points"]:
            c = 2 * mp.pi * mp.mpf(xs)
            mp.mp.dps = int(2 * float(c) / 2.3) + 40
            lk = leakage(mp.mpf(-1) / 2, c, n // 2)
            ub = mp.mpf(lam)
            cert = mp.mpf(CERT[par][xs]) if xs in CERT[par] else None
            pts.append((mp.mpf(xs), ub / _ell(n, c), ub / (lk / 2), cert / (lk / 2) if cert else None))
        mp.mp.dps = 30
        rows_out.append({"object": "zeta fit points", "sector": par, "n_rule": n, "x": [mp.nstr(p[0], 6) for p in pts],
                         "R_fuchs_ub": [mp.nstr(p[1], 4) for p in pts], "R_exact_ub": [mp.nstr(p[2], 4) for p in pts],
                         "R_exact_cert": [mp.nstr(p[3], 4) if p[3] else None for p in pts]})
        print(f"zeta fit {par} n={n}: x = {[mp.nstr(p[0],4) for p in pts]}; R_F(ub) {[mp.nstr(p[1],3) for p in pts]}; "
              f"R_ex(ub) {[mp.nstr(p[2],3) for p in pts]}; R_ex(cert) {[mp.nstr(p[3],3) if p[3] else None for p in pts]}", flush=True)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--json")
    args = ap.parse_args()
    out = []
    degree1(out)
    zeta_fit(out)
    gl2(out)
    if args.json:
        with open(args.json, "w") as fh:
            json.dump(out, fh, indent=1)


if __name__ == "__main__":
    main()
