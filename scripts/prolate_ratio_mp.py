#!/usr/bin/env python3
"""Prolate ratio test: ε(λ) / (1 − χ₂(λ)) for ζ and Davenport–Heilbronn.

Connes (arXiv:2602.04022, §6.4, Fig. 1) observes that the smallest eigenvalue
ε(λ) of the Weil form QW_λ tracks 1 − χ₂(λ), a purely archimedean quantity:
χ₂ is the eigenvalue of the time- and band-limited Fourier transform on
[−λ, λ] for the prolate function h_{4,λ}. In Slepian's notation χ₂² = λ₄(c)
with bandwidth c = 2πλ² = 2πx.

λ₄(c) is computed by Bouwkamp's method: even prolate functions ψ_n are the
eigenvectors of the prolate operator −d/dx(1−x²)d/dx + c²x² in the normalised
Legendre basis (a tridiagonal matrix), and for even n the Fourier eigenvalue is
μ_n = ∫ψ_n / ψ_n(0) = √2 β₀ / ψ_n(0), with λ_n = (c/2π) μ_n².

ε(λ) comes from scripts/connes_letter_mp.py (ζ: pole-free form; DH: full form).

Usage: .venv/bin/python scripts/prolate_ratio_mp.py --function zeta --x 3,5,7,10,13 --n 100 --dps 150
"""

import argparse
import json
import os
import sys
import time

import mpmath as mp

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import connes_letter_mp as cl  # noqa: E402


def prolate_even_lambdas(c, count=3, k_max=None):
    """Slepian energy eigenvalues λ_0, λ_2, λ_4, … (even n) at bandwidth c."""
    c = mp.mpf(c)
    if k_max is None:
        k_max = int(2 * c + 80)
    ks = list(range(0, k_max + 1, 2))
    m = len(ks)
    T = mp.matrix(m, m)
    for i, k in enumerate(ks):
        T[i, i] = k * (k + 1) + c * c * (2 * k * (k + 1) - 1) / ((2 * k + 3) * (2 * k - 1))
        if i + 1 < m:
            off = c * c * (k + 2) * (k + 1) / ((2 * k + 3) * mp.sqrt((2 * k + 1) * (2 * k + 5)))
            T[i, i + 1] = off
            T[i + 1, i] = off
    vals, vecs = mp.eigsy(T)
    order = sorted(range(m), key=lambda i: vals[i])

    # P̄_k(0) = sqrt(k + 1/2) P_k(0), P_k(0) = (−1)^{k/2} (k−1)!!/k!! for even k.
    pbar0 = []
    p0 = mp.mpf(1)
    for k in ks:
        if k > 0:
            p0 = -p0 * (k - 1) / k
        pbar0.append(mp.sqrt(k + mp.mpf(1) / 2) * p0)

    out = []
    for j in order[:count]:
        beta = [vecs[i, j] for i in range(m)]
        psi0 = mp.fsum(b * p for b, p in zip(beta, pbar0))
        mu = mp.sqrt(2) * beta[0] / psi0
        out.append(c / (2 * mp.pi) * mu * mu)
    return out


PROLATE_DPS = None


def one_minus_chi2(x):
    """1 − χ₂(√x) with χ₂ = √λ₄(2πx), at PROLATE_DPS digits if set."""
    with mp.workdps(PROLATE_DPS or mp.mp.dps):
        lam = prolate_even_lambdas(2 * mp.pi * x, count=3)
        l4 = lam[2]
        omc = (1 - l4) / (1 + mp.sqrt(l4))
    return +omc, [+v for v in lam]


def connes_asymptotic(x):
    """Connes §6.4: 1 − χ₂ ~ (2^14/3) √2 π^5 exp(−4π e^L + 9L/2), e^L = x."""
    return mp.mpf(2) ** 14 / 3 * mp.sqrt(2) * mp.pi ** 5 * mp.exp(-4 * mp.pi * x + mp.mpf(9) / 2 * mp.log(x))


def smallest_eigenvalue(function, x, n):
    L, E, pp = cl.build_form(mp.mpf(x), n, function)
    n1 = n + 1
    if function == "zeta":
        v = mp.matrix(n1, 1)
        v[0] = cl.int_cos_cosh(0, L) / mp.sqrt(L)
        for k in range(1, n1):
            v[k] = mp.sqrt(2 / L) * cl.int_cos_cosh(2 * mp.pi * k / L, L)
        Q = E + 2 * v * v.T
    else:
        Q = E
    vals, _ = mp.eigsy(Q)
    return min(vals[i] for i in range(n1))


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--function", choices=["zeta", "dh", "prolate"], default="zeta")
    ap.add_argument("--x", default="3,5,7,10,13", help="values of x = λ²")
    ap.add_argument("--n", default="100", help="basis sizes N (comma list; the last is used for the ratio)")
    ap.add_argument("--dps", type=int, default=150)
    ap.add_argument("--prolate-dps", type=int, help="working precision for 1 − χ₂ (default: --dps)")
    ap.add_argument("--json")
    args = ap.parse_args()
    mp.mp.dps = args.dps
    global PROLATE_DPS
    PROLATE_DPS = args.prolate_dps
    xs = [mp.mpf(v) for v in args.x.split(",")]
    ns = [int(v) for v in args.n.split(",")]
    from progress import Job

    job = Job(f"prolate ratio {args.function}", total=len(xs) * (1 + (len(ns) if args.function != "prolate" else 0)), args=vars(args))
    rows = []
    for x in xs:
        t0 = time.time()
        omc, lam = one_minus_chi2(x)
        job.step(f"x={mp.nstr(x, 4)} 1-chi2={mp.nstr(omc, 6)}")
        asym = connes_asymptotic(x)
        row = {"x": float(x), "one_minus_chi2": mp.nstr(omc, 12), "asymptotic": mp.nstr(asym, 12),
               "lambda0_2_4": [mp.nstr(1 - l, 6) for l in lam]}
        line = f"x={mp.nstr(x, 4):>5}  1−χ₂ = {mp.nstr(omc, 6):>12}  (Connes asymptotic {mp.nstr(asym, 6)}, ratio {mp.nstr(omc / asym, 6)})"
        if args.function != "prolate":
            eps = []
            for n in ns:
                e = smallest_eigenvalue(args.function, x, n)
                eps.append((n, e))
                job.step(f"x={mp.nstr(x, 4)} N={n} eps0={mp.nstr(e, 6)}")
            e_last = eps[-1][1]
            row["eps"] = {str(n): mp.nstr(e, 12) for n, e in eps}
            row["ratio"] = mp.nstr(e_last / omc, 8)
            line += "   ε₀: " + ", ".join(f"N={n}: {mp.nstr(e, 6)}" for n, e in eps) + f"   ε/(1−χ₂) = {mp.nstr(e_last / omc, 6)}"
        print(line + f"   [{time.time() - t0:.0f}s]", flush=True)
        job.result(line.strip())
        rows.append(row)
    if args.json:
        with open(args.json, "w") as f:
            json.dump({"function": args.function, "dps": args.dps, "rows": rows}, f, indent=2)
    job.done()


if __name__ == "__main__":
    main()
