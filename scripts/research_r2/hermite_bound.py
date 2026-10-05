#!/usr/bin/env python3
"""Theorem 2 with Hermite trial functions (Riemann's k for zeta's even sector): an unconditional, explicit upper
bound for the window minimum, evaluated on the decay-law grid, against the measured minima.

For each (D, x, sector): psi = Hermite trial of index n = kappa + 2s + 4*pole (r2lib.hermite_choice),
f = E-map window function, and
    ||f_win||^2 = 2 int_0^{L/2} f^2 du                (sector symmetry f(-u) = sigma f(u), checked)
    B = 2 int_{L/2}^inf |f| e^{u/2} du
    A = 2 (|f(L/2)| e^{L/4} + int_{L/2}^inf |f'| e^{u/2} du)
    bound = zero_sum_bound(A, B) / ||f_win||^2         (>= lambda_min in that sector, unconditionally)
Also reported: the trial function's Rayleigh quotient in rh2's basis (matrix from the explicit formula; also an
unconditional upper bound) as a consistency check that Q(f_win) is as small as the bound says.

Usage: .venv/bin/python scripts/research_r2/hermite_bound.py --d 1,5,-4 --xq 2,3,4,6 --json data/research_r2/hermite_bound.json
"""
import argparse
import json
import os
import sys
import time

import mpmath as mp

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
import r2lib  # noqa: E402
from decay_law_mp import zeros_side, lam_min_acb  # noqa: E402
from progress import Job  # noqa: E402


def measured(D, xq, parity):
    import glob
    for f in [g for g in glob.glob(f"data/connes/decay/scan*_D{D}.json") if g.split("/")[-1].split("_")[0] in ("scanhi", "scanext", "scan8")]:
        for r in json.load(open(f)):
            if abs(float(r["xq"]) - float(xq)) < 1e-9 and r["parity"] == parity:
                return mp.mpf(r["lambda_n2"]), r["n2"], f
    return None, None, None


def project(win, L, n, parity, nodes):
    """Coefficients of f_win in rh2's orthonormal basis by Gauss-Legendre on [-L/2, L/2] (nodes on [-1,1])."""
    panels = n // 3 + 6                      # composite Gauss-Legendre, >= 3 nodes per oscillation of the top mode
    width = L / panels
    vals = []
    for p in range(panels):
        a = -L / 2 + p * width
        for t, wt in nodes:
            u = a + (t + 1) * width / 2
            vals.append((u, wt * width / 2, win.f(u)))
    if parity == "even":
        c = [mp.fsum(wt * fv for _, wt, fv in vals) / mp.sqrt(L)]
        c += [mp.sqrt(2 / L) * mp.fsum(wt * fv * mp.cos(2 * mp.pi * k * u / L) for u, wt, fv in vals) for k in range(1, n + 1)]
    else:
        c = [mp.sqrt(2 / L) * mp.fsum(wt * fv * mp.sin(2 * mp.pi * k * u / L) for u, wt, fv in vals) for k in range(1, n + 1)]
    return c, mp.fsum(wt * fv * fv for _, wt, fv in vals)


def gl_nodes(m):
    """Gauss-Legendre nodes/weights on [-1, 1] (Golub-Welsch, scripts/k_lambda_mp.py)."""
    from k_lambda_mp import gl_nodes as gln
    return gln(m)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--d", default="1,5,8,-3,-4,-7,-20")
    ap.add_argument("--xq", default="2,3,4,5,6,7,8,10")
    ap.add_argument("--dps", type=int, default=40)
    ap.add_argument("--matrix", action="store_true", help="also the Rayleigh quotient in rh2's basis (slow)")
    ap.add_argument("--json")
    args = ap.parse_args()
    mp.mp.dps = args.dps
    Ds = [int(v) for v in args.d.split(",")]
    xqs = [mp.mpf(v) for v in args.xq.split(",")]
    job = Job("r2 hermite bound", total=len(Ds) * len(xqs) * 2, args=vars(args))
    rows = []
    nodes = gl_nodes(20) if args.matrix else None
    for D in Ds:
        q = 1 if D == 1 else abs(D)
        kappa = 0 if D > 0 else 1
        pole = D == 1
        for xq in xqs:
            x = xq * q
            L = mp.log(x)
            c = 2 * mp.pi * xq
            for s, parity in enumerate(("even", "odd")):
                t0 = time.time()
                psi = r2lib.hermite_choice(kappa, s, pole)
                win = r2lib.Window(psi, D, x)
                # sector check at two interior points
                sig = [win.f(-u) / win.f(u) for u in (L / 7, L / 3)]
                sigma = 1 if parity == "even" else -1
                assert all(abs(v - sigma) < mp.mpf(10) ** (-args.dps // 2) for v in sig), (D, xq, parity, sig)
                half = L / 2
                brk = [mp.mpf(0), half / 4, half / 2, 3 * half / 4, half]
                norm2 = 2 * mp.quad(lambda u: win.f(u) ** 2, brk)
                # out-of-window pieces: f decays like e^{-pi w^2/q}, so integrate to where y = w/sqrt q reaches ycut
                umax = mp.log(win.ycut * win.sq)
                ub = [half + (umax - half) * j / 40 for j in range(41)]
                B = 2 * mp.quad(lambda u: abs(win.f(u)) * mp.exp(u / 2), ub)
                A = 2 * (abs(win.f(half)) * mp.exp(L / 4) + mp.quad(lambda u: abs(win.fprime(u)) * mp.exp(u / 2), ub))
                if D == 1:
                    Z = r2lib.zero_sum_bound(A, B, r2lib.zeta_Nplus, zeta=True, t0=mp.mpf("14.134"))
                else:
                    Z = r2lib.zero_sum_bound(A, B, lambda t: r2lib.dirichlet_Nplus(t, q, kappa), zeta=False)
                bound = Z / norm2
                lam, n2, src = measured(D, xq, parity)
                row = {"D": D, "q": q, "xq": mp.nstr(xq, 6), "x": mp.nstr(x, 10), "parity": parity, "n_index": psi.n,
                       "c": mp.nstr(c, 10), "norm2": mp.nstr(norm2, 12), "A": mp.nstr(A, 12), "B": mp.nstr(B, 12),
                       "Tc": mp.nstr(A / B, 8), "zero_sum_bound": mp.nstr(Z, 12), "bound": mp.nstr(bound, 12),
                       "measured": mp.nstr(lam, 12) if lam is not None else None,
                       "ratio_measured_over_bound": mp.nstr(lam / bound, 6) if lam is not None else None,
                       "ln_bound_plus_c": mp.nstr(mp.log(bound) + c, 8),
                       "ln_measured_plus_2c": mp.nstr(mp.log(lam) + 2 * c, 8) if lam is not None else None}
                if args.matrix and lam is not None:
                    M = zeros_side(D, x, n2, parity)
                    cv, _ = project(win, L, n2, parity, nodes)
                    cm = mp.matrix(cv)
                    rq = (cm.T * M * cm)[0] / mp.fsum(v * v for v in cv)
                    row["rayleigh_matrix"] = mp.nstr(rq, 12)
                    row["matrix_n"] = n2
                rows.append(row)
                job.step(f"D={D} x/q={mp.nstr(xq,4)} {parity} n={psi.n}: bound={mp.nstr(bound,6)} measured={mp.nstr(lam,6) if lam else None} "
                         f"ratio={row['ratio_measured_over_bound']} Tc={mp.nstr(A/B,5)} ({time.time()-t0:.0f}s)"
                         + (f" RQ={row['rayleigh_matrix']}" if 'rayleigh_matrix' in row else ""))
                if args.json:
                    with open(args.json, "w") as fh:
                        json.dump(rows, fh, indent=1)
    job.done()


if __name__ == "__main__":
    main()
