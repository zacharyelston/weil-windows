#!/usr/bin/env python3
"""Generalised (Hankel) prolate leakage 1 - lambda_{nu,m}(c) and the continuous prolate index it implies.

For a weight-k GL(2) newform the Voronoi formula replaces Poisson: the self-reciprocal transform is the Hankel
transform of order nu = k - 1 in the square-root variable, and the E-map construction of
docs/RESEARCH_R2_SAMPLING.md section 5 leads to the finite Hankel transform (Slepian 1964, "PSWF IV")
    (H phi)(r) = int_0^1 J_nu(c r rho) sqrt(c r rho) phi(rho) d rho,   lambda = c gamma^2   (concentration),
with c = 4 pi sqrt(x/N) (= d T*).  nu = -1/2, +1/2 reproduce the 1D even/odd prolates (n = 2m, 2m+1).

Method: the commuting operator -d/dr (1-r^2) d/dr + (nu^2 - 1/4)/r^2 + c^2 r^2 is tridiagonal in the orthonormal basis
T_k(r) = 2^{(nu+2)/2} r^{nu+1/2} p_k(1 - 2r^2), p_k the orthonormal Jacobi polynomials for (1-x)^nu (Jacobi matrix J),
with c = 0 eigenvalues (nu+2k+1/2)(nu+2k+3/2) and r^2 = (I - J)/2. The eigenvalue gamma follows from r -> 0:
    gamma = c^{nu+1/2} d_0 / (2^nu Gamma(nu+1) sqrt(2(nu+1)) sum_k d_k sqrt(2(2k+nu+1)) binom(k+nu, k)).
The continuous index n* (P-DL11 metric): slope of ln(1 - lambda) + 2c against ln c over c = 4 pi v, v = 6..10, minus 1/2.

Usage: .venv/bin/python scripts/research_r2/hankel_prolate.py --check
       .venv/bin/python scripts/research_r2/hankel_prolate.py --nu 1,3,5,7,11,15,17,19,21,25 --json data/research_r2/hankel_index.json
"""
import argparse
import json
import os
import sys
import time

import mpmath as mp

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from progress import Job  # noqa: E402


def jacobi_rec(nu, K):
    """Orthonormal Jacobi recurrence for weight (1-x)^nu (alpha = nu, beta = 0): x p_k = a_{k+1} p_{k+1} + b_k p_k + a_k p_{k-1}."""
    al = mp.mpf(nu)
    b = [-al / (al + 2)]
    for k in range(1, K):
        s = 2 * k + al
        b.append(-al * al / (s * (s + 2)))
    a = [mp.mpf(0)]
    for k in range(1, K):
        s = 2 * k + al
        a.append(mp.sqrt(4 * k * k * (k + al) ** 2 / (s * s * (s + 1) * (s - 1))))
    return a, b


def tridiag(nu, c, K):
    a, b = jacobi_rec(nu, K)
    nu = mp.mpf(nu)
    diag = [(nu + 2 * k + mp.mpf(1) / 2) * (nu + 2 * k + mp.mpf(3) / 2) + c * c * (1 - b[k]) / 2 for k in range(K)]
    off = [-c * c * a[k + 1] / 2 for k in range(K - 1)]
    return diag, off


def sturm_count(diag, off, x):
    """Number of eigenvalues < x of the symmetric tridiagonal matrix."""
    cnt, d = 0, mp.mpf(1)
    for i in range(len(diag)):
        d = diag[i] - x - (off[i - 1] ** 2 / d if i > 0 else 0)
        if d == 0:
            d = mp.mpf(10) ** (-mp.mp.dps)
        if d < 0:
            cnt += 1
    return cnt


def eig_m(diag, off, m):
    """m-th smallest eigenvalue by bisection and its eigenvector by inverse iteration (Thomas algorithm)."""
    lo = min(diag) - 2 * max(abs(v) for v in off) - 1
    hi = max(diag) + 2 * max(abs(v) for v in off) + 1
    for _ in range(int(mp.mp.prec) + 20):
        mid = (lo + hi) / 2
        if sturm_count(diag, off, mid) > m:
            hi = mid
        else:
            lo = mid
    lamb = (lo + hi) / 2
    K = len(diag)
    v = [mp.mpf(1)] * K
    shift = lamb + (hi - lo)
    for _ in range(4):
        # solve (T - shift) y = v
        cp, dp = [mp.mpf(0)] * K, [mp.mpf(0)] * K
        bb = [diag[i] - shift for i in range(K)]
        cp[0] = off[0] / bb[0]
        dp[0] = v[0] / bb[0]
        for i in range(1, K):
            den = bb[i] - off[i - 1] * cp[i - 1]
            cp[i] = off[i] / den if i < K - 1 else 0
            dp[i] = (v[i] - off[i - 1] * dp[i - 1]) / den
        y = [mp.mpf(0)] * K
        y[-1] = dp[-1]
        for i in range(K - 2, -1, -1):
            y[i] = dp[i] - cp[i] * y[i + 1]
        nrm = mp.sqrt(mp.fsum(t * t for t in y))
        v = [t / nrm for t in y]
    if v[0] < 0:
        v = [-t for t in v]
    return lamb, v


def leakage(nu, c, m, K=None):
    """1 - lambda_{nu,m}(c)."""
    c = mp.mpf(c)
    nu = mp.mpf(nu)
    if K is None:
        K = int(c + 3 * m + 100)
    diag, off = tridiag(nu, c, K)
    _, d = eig_m(diag, off, m)
    num = c ** (nu + mp.mpf(1) / 2) * d[0] / (2 ** nu * mp.gamma(nu + 1) * mp.sqrt(2 * (nu + 1)))
    den = mp.fsum(d[k] * mp.sqrt(2 * (2 * k + nu + 1)) * mp.binomial(k + nu, k) for k in range(K))
    gam = num / den
    lam = c * gam * gam
    return 1 - lam


def fuchs(n, c):
    return 4 * mp.sqrt(mp.pi) * mp.mpf(8) ** n * c ** (n + mp.mpf(1) / 2) * mp.exp(-2 * c) / mp.factorial(n)


def check():
    sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
    from prolate_ratio_mp import prolate_even_lambdas
    for c in (10, 25):
        mp.mp.dps = int(2 * c / 2.3) + 40
        lam1d = prolate_even_lambdas(c, count=3)
        for m in range(3):
            v = leakage(-0.5, c, m)
            print(f"c={c} nu=-1/2 m={m}: 1-lambda = {mp.nstr(v, 12)}   1D Bouwkamp n={2*m}: {mp.nstr(1 - lam1d[m], 12)}   Fuchs {mp.nstr(fuchs(2*m, mp.mpf(c)), 6)}")
        for m in range(2):
            v = leakage(0.5, c, m)
            print(f"c={c} nu=+1/2 m={m}: 1-lambda = {mp.nstr(v, 12)}   Fuchs n={2*m+1}: {mp.nstr(fuchs(2*m+1, mp.mpf(c)), 6)}")
        v1, v2 = leakage(3, c, 0), leakage(3, c, 0, K=int(c / 2 + 140))
        print(f"c={c} nu=3 m=0 basis check: {mp.nstr(v1, 15)} vs {mp.nstr(v2, 15)}")


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--check", action="store_true")
    ap.add_argument("--nu", default="1,3,5,7,11,15,17,19,21,25")
    ap.add_argument("--v", default="6,7,8,9,10")
    ap.add_argument("--m", default="0,1")
    ap.add_argument("--json")
    args = ap.parse_args()
    if args.check:
        check()
        return
    nus = [mp.mpf(s) for s in args.nu.split(",")]
    vs = [mp.mpf(s) for s in args.v.split(",")]
    ms = [int(s) for s in args.m.split(",")]
    job = Job("r2 hankel prolate", total=len(nus) * len(ms), args=vars(args))
    rows = []
    for nu in nus:
        for m in ms:
            t0 = time.time()
            pts = []
            for v in vs:
                c = 4 * mp.pi * v
                mp.mp.dps = int(2 * float(c) / 2.3) + 40
                lk = leakage(nu, c, m)
                pts.append((mp.log(c), mp.log(lk) + 2 * c, lk, c))
            mp.mp.dps = 30
            xs = [p[0] for p in pts]
            ys = [p[1] for p in pts]
            xb, yb = mp.fsum(xs) / len(xs), mp.fsum(ys) / len(ys)
            slope = mp.fsum((a - xb) * (b - yb) for a, b in zip(xs, ys)) / mp.fsum((a - xb) ** 2 for a in xs)
            nstar = slope - mp.mpf(1) / 2
            k = nu + 1
            row = {"nu": mp.nstr(nu, 4), "k": mp.nstr(k, 4), "m": m, "v": [mp.nstr(v, 4) for v in vs],
                   "c": [mp.nstr(p[3], 8) for p in pts], "leakage": [mp.nstr(p[2], 10) for p in pts],
                   "n_star": mp.nstr(nstar, 6), "n_star_minus_k": mp.nstr(nstar - k, 6),
                   "asymptotic_n_minus_k": mp.nstr(2 * m - mp.mpf(1) / 2, 4)}
            rows.append(row)
            job.step(f"nu={mp.nstr(nu,3)} (k={mp.nstr(k,3)}) m={m}: n*-k = {row['n_star_minus_k']} (asymptotic {row['asymptotic_n_minus_k']}) ({time.time()-t0:.0f}s)")
            if args.json:
                with open(args.json, "w") as fh:
                    json.dump(rows, fh, indent=1)
    job.done()


if __name__ == "__main__":
    main()
