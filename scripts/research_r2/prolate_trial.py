#!/usr/bin/env python3
"""Sharp E-map trial functions (Connes' k_lambda generalised) on the decay-law grid: how close does the
prolate construction come to the measured window minimum?

For L(s, chi_D) (D = 1: zeta) and sector s, with lam = sqrt(x/q), c = 2 pi lam^2 = 2 pi x/q:
    phi = time-limited prolate xi_n on [-lam, lam], n = kappa + 2s           (no pole)
    phi = beta0(n+4) xi_n - beta0(n) xi_{n+4}  (vanishing integral)        (zeta, n = 2s)
    G(w) = w^{1/2} sum_{n <= lam sqrt(q)/w} chi(n) phi(n w / sqrt q),  f(u) = G(e^u),  u in [-L/2, L/2]
    f_s = (f(u) + sigma f(-u))/2  (sector projection; f is in the sector up to leakage).
The Rayleigh quotient of f_s in rh2's N-mode basis (N = the scan's n2), with rh2's explicit-formula matrix
(decay_law_mp.zeros_side), is an unconditional upper bound for the window minimum; it is compared with the scan's
lambda at the same N and with the Fuchs-Slepian leakage l_n(c).

Usage: .venv/bin/python scripts/research_r2/prolate_trial.py --d 5,-4 --xq 2,3,4 --json data/research_r2/prolate_trial.json
"""
import argparse
import glob
import json
import os
import sys
import time

import mpmath as mp

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from decay_law_mp import _ell, kronecker, zeros_side  # noqa: E402
from k_lambda_mp import gl_nodes, prolate_even  # noqa: E402
from progress import Job  # noqa: E402


def legendre_vals(t, kmax):
    p = [mp.mpf(1), t]
    for k in range(1, kmax):
        p.append(((2 * k + 1) * t * p[k] - k * p[k - 1]) / (k + 1))
    return p


class Phi:
    def __init__(self, c, n, lam, pole):
        parity = n % 2
        if pole:
            ks, (ba, bb) = prolate_even(c, which=(n // 2, n // 2 + 2), parity=parity)
            self.coef = [bb[0] * a - ba[0] * b for a, b in zip(ba, bb)]
        else:
            ks, (b,) = prolate_even(c, which=(n // 2,), parity=parity)
            self.coef = b
        self.ks, self.lam = ks, mp.mpf(lam)
        self.norm = [mp.sqrt(k + mp.mpf(1) / 2) for k in ks]

    def __call__(self, y):
        t = y / self.lam
        if abs(t) > 1:
            return mp.mpf(0)
        P = legendre_vals(t, self.ks[-1])
        return mp.fsum(b * nk * P[k] for b, nk, k in zip(self.coef, self.norm, self.ks))


def measured(D, xq, parity):
    for f in [g for g in glob.glob(f"data/connes/decay/scan*_D{D}.json") if g.split("/")[-1].split("_")[0] in ("scanhi", "scanext", "scan8")]:
        for r in json.load(open(f)):
            if abs(float(r["xq"]) - float(xq)) < 1e-9 and r["parity"] == parity:
                return mp.mpf(r["lambda_n2"]), r["n2"], r["dps"]
    return None, None, None


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--d", default="1,5,8,-3,-4,-7,-20")
    ap.add_argument("--xq", default="2,3,4,6")
    ap.add_argument("--json")
    args = ap.parse_args()
    Ds = [int(v) for v in args.d.split(",")]
    xqs = [mp.mpf(v) for v in args.xq.split(",")]
    job = Job("r2 prolate trial", total=len(Ds) * len(xqs) * 2, args=vars(args))
    rows = []
    for D in Ds:
        q = 1 if D == 1 else abs(D)
        kappa = 0 if D > 0 else 1
        pole = D == 1
        chi = (lambda n: 1) if D == 1 else (lambda n, D=D: kronecker(D, n))
        for xq in xqs:
            for s, parity in enumerate(("even", "odd")):
                lam_m, n2, dps = measured(D, xq, parity)
                if lam_m is None:
                    continue
                t0 = time.time()
                mp.mp.dps = dps
                x = xq * q
                L = mp.log(x)
                c = 2 * mp.pi * xq
                lam = mp.sqrt(xq)
                n = kappa + 2 * s
                phi = Phi(c, n, lam, pole)
                sq = mp.sqrt(q)
                top = lam * sq                     # = sqrt(x): G(w) = 0 for w > sqrt(x)

                def G(w):
                    nmax = int(mp.floor(top / w + mp.mpf(10) ** (-dps // 2)))
                    return mp.sqrt(w) * mp.fsum(chi(k) * phi(k * w / sq) for k in range(1, nmax + 1) if chi(k))

                # composite Gauss-Legendre on [0, L/2], split where a term switches on (w = sqrt(x)/k)
                bps = sorted({mp.mpf(0), L / 2} | {mp.log(top / k) for k in range(1, int(top) + 1) if 0 < mp.log(top / k) < L / 2})
                gl = gl_nodes(16)
                hmax = L / (n2 / 2 + 8)
                nodes = []
                for a, b in zip(bps[:-1], bps[1:]):
                    p = int(mp.ceil((b - a) / hmax))
                    h = (b - a) / p
                    for i in range(p):
                        lo = a + i * h
                        nodes += [(lo + (t + 1) * h / 2, wt * h / 2) for t, wt in gl]
                sigma = 1 if parity == "even" else -1
                vals = [(u, wt, (G(mp.exp(u)) + sigma * G(mp.exp(-u))) / 2) for u, wt in nodes]
                # f_s on [-L/2, 0] is sigma times the mirror image; integrate over [0, L/2] and double
                if parity == "even":
                    cvec = [2 * mp.fsum(wt * fv for _, wt, fv in vals) / mp.sqrt(L)]
                    cvec += [2 * mp.sqrt(2 / L) * mp.fsum(wt * fv * mp.cos(2 * mp.pi * k * u / L) for u, wt, fv in vals) for k in range(1, n2 + 1)]
                else:
                    cvec = [2 * mp.sqrt(2 / L) * mp.fsum(wt * fv * mp.sin(2 * mp.pi * k * u / L) for u, wt, fv in vals) for k in range(1, n2 + 1)]
                norm2 = 2 * mp.fsum(wt * fv * fv for _, wt, fv in vals)
                cn2 = mp.fsum(v * v for v in cvec)
                M = zeros_side(D, x, n2, parity)
                cm = mp.matrix(cvec)
                rq = (cm.T * M * cm)[0] / cn2
                nidx = n + 4 * pole
                ell = _ell(nidx, c)
                row = {"D": D, "q": q, "xq": mp.nstr(xq, 6), "x": mp.nstr(x, 10), "parity": parity, "n_index": nidx, "N": n2,
                       "dps": dps, "c": mp.nstr(c, 10), "rayleigh": mp.nstr(rq, 12), "measured": mp.nstr(lam_m, 12),
                       "trial_over_measured": mp.nstr(rq / lam_m, 6), "outside_span": mp.nstr(1 - cn2 / norm2, 4),
                       "R_trial": mp.nstr(rq / ell, 6), "R_measured": mp.nstr(lam_m / ell, 6)}
                rows.append(row)
                job.step(f"D={D} x/q={mp.nstr(xq,4)} {parity} n={nidx} N={n2}: RQ={mp.nstr(rq,6)} measured={mp.nstr(lam_m,6)} "
                         f"trial/measured={row['trial_over_measured']} R_trial={row['R_trial']} R_meas={row['R_measured']} ({time.time()-t0:.0f}s)")
                if args.json:
                    with open(args.json, "w") as fh:
                        json.dump(rows, fh, indent=1)
    job.done()


if __name__ == "__main__":
    main()
