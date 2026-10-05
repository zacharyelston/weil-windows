#!/usr/bin/env python3
"""Connes' form QW_λ for the discriminant −20 pair (degree 2, conductor 20).

Step 3 of the Epstein control (scripts/epstein_mp.py: steps 1–2).
  Z₁  : Epstein zeta of x² + 5y² (no Euler product; off-line zeros from 15.67i up)
  ζ_K : Dedekind zeta of Q(√−5) = ζ(s) L(s, χ₋₂₀) (Euler product)
Both have Λ(s) = ½ 20^{s/2} Γ_R(s) Γ_R(s+1) F(s) = Λ(1 − s) and a simple pole at
s = 1, so

  Q_arch = Q_arch[Γ_R(s)] + Q_arch[Γ_R(s+1)] + log 20 · I
         = Q_ζ,arch + Q_DH,arch + log 4 · I      (Q_DH,arch already holds log 5 · I),
  Q      = Q_arch − Σ_{n ≤ x} c_n/√n · (prime-side term),   −F'/F = Σ c_n n^{−s},
  QW_λ   = Q + 2 v vᵀ (pole term dropped, as for ζ).

`validate` checks the explicit formula on a C³ test function against the census
zeros (on- and off-line); `sweep` reports ε₀ and the number of negative
eigenvalues of QW_λ for both functions.

Usage:
  .venv/bin/python scripts/epstein_connes_mp.py validate --zeros z1_census_60.json
  .venv/bin/python scripts/epstein_connes_mp.py sweep --x 2,3,4,5,7,10 --n 48 --dps 50
"""

import argparse
import json
import math
import os
import sys

import mpmath as mp

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import connes_letter_mp as cl  # noqa: E402
import epstein_mp as ep  # noqa: E402
from progress import Job  # noqa: E402


def coefficients(function, n_max):
    if function == "Z1":
        direct, _ = ep.coefficients(n_max)
        return direct
    # ζ_K: a_n = Σ_{d | n} χ₋₂₀(d)
    return [mp.mpf(0)] + [mp.mpf(sum(ep.CHI_M20[d % 20] for d in range(1, n + 1) if n % d == 0)) for n in range(1, n_max + 1)]


def log_derivative(a, n_max):
    c = [mp.mpf(0)] * (n_max + 1)
    for k in range(2, n_max + 1):
        c[k] = a[k] * mp.log(k)
    for d in range(2, n_max + 1):
        if c[d] == 0:
            continue
        for m in range(2 * d, n_max + 1, d):
            c[m] -= c[d] * a[m // d]
    return c


def pole_vector(x, n):
    L = mp.log(x)
    return mp.matrix([cl.int_cos_cosh(2 * mp.pi * k / L, L) * (1 / mp.sqrt(L) if k == 0 else mp.sqrt(2 / L)) for k in range(n + 1)])


def forms(function, x, n):
    """(E, QW): E = −Σ_v W_v (with pole), QW = E + 2vvᵀ (Connes' pole-free form)."""
    x = mp.mpf(x)
    n_max = int(mp.floor(x))
    a = coefficients(function, max(n_max, 2))
    c = log_derivative(a, max(n_max, 2))
    terms = [(mp.log(k), c[k] / mp.sqrt(k)) for k in range(2, n_max + 1) if c[k] != 0]
    _, Qz, _ = cl.build_form(x, n, "zeta", terms=[])
    _, Qd, _ = cl.build_form(x, n, "dh", terms=[])
    _, Qp, _ = cl.build_form(x, n, "zeta", terms=terms, include_arch=False) if terms else (None, mp.zeros(n + 1, n + 1), None)
    E = Qz + Qd + mp.log(4) * mp.eye(n + 1) + Qp
    v = pole_vector(x, n)
    return E, E + 2 * v * v.T, v, terms


def fourier_even_complex(c, z, L):
    n = len(c) - 1
    sinc = lambda w: L / 2 if abs(w) < mp.mpf(10) ** (-mp.mp.dps + 5) else mp.sin(w * L / 2) / w
    val = c[0] * 2 * sinc(z) / mp.sqrt(L)
    for k in range(1, n + 1):
        om = 2 * mp.pi * k / L
        val += c[k] * mp.sqrt(2 / L) * (sinc(z - om) + sinc(z + om))
    return val


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    v = sub.add_parser("validate")
    v.add_argument("--zeros", required=True, help="census JSON for Z1 from epstein_mp.py")
    v.add_argument("--x", type=float, default=13)
    v.add_argument("--dps", type=int, default=30)
    s = sub.add_parser("sweep")
    s.add_argument("--x", default="2,3,4,5,7,10")
    s.add_argument("--n", type=int, default=48)
    s.add_argument("--dps", type=int, default=50)
    s.add_argument("--json")
    args = ap.parse_args()
    mp.mp.dps = args.dps

    if args.cmd == "validate":
        x = mp.mpf(args.x)
        L = mp.log(x)
        nb = 8
        E, _, vpole, _ = forms("Z1", x, nb)
        # φ(u) = (1 + cos(2πu/L))² = 3/2 + 2 cos θ + ½ cos 2θ: vanishes to 4th order at ±L/2.
        c = [mp.mpf(0)] * (nb + 1)
        c[0] = mp.mpf(3) / 2 * mp.sqrt(L)
        c[1] = 2 * mp.sqrt(L / 2)
        c[2] = mp.mpf(1) / 2 * mp.sqrt(L / 2)
        cv = mp.matrix(c)
        q = (cv.T * E * cv)[0, 0]
        with open(args.zeros) as f:
            d = json.load(f)
        on = mp.fsum(fourier_even_complex(c, mp.mpf(g), L) ** 2 for g in d["on_line"])
        off = mp.fsum(2 * mp.re(fourier_even_complex(c, mp.mpc(z["gamma"], -(z["beta"] - 0.5)), L) ** 2) for z in d["off_line"])
        pole = 2 * mp.fsum(cc * vv for cc, vv in zip(c, vpole)) ** 2
        rhs = 2 * (on + off) - pole
        print(f"Z1 explicit formula at x = {args.x}: QW(φ) = {mp.nstr(q, 12)};  2Σ_ρ − 2φ̂(i/2)² = {mp.nstr(rhs, 12)}  "
              f"(on {mp.nstr(on, 6)}, off {mp.nstr(off, 6)}, pole {mp.nstr(pole, 6)});  rel. diff {mp.nstr(abs(q - rhs) / abs(q), 3)}")
        tail = fourier_even_complex(c, mp.mpf(d["t2"]), L) ** 2
        print(f"  zeros used up to γ = {d['t2']}; φ̂(γ_max)² = {mp.nstr(tail, 3)} sets the truncation scale")
        return

    rows = []
    xs_list = args.x.split(",")
    job = Job("epstein connes sweep", total=2 * len(xs_list), args=vars(args))
    print(f"{'x':>5} {'function':>8}  {'ε₀':>12} {'ε₁':>12} {'n₋':>3}   prime side (n: c_n/√n)")
    for xs in xs_list:
        for fn in ["zeta_K", "Z1"]:
            E, Q, _, terms = forms(fn, xs, args.n)
            vals = sorted(mp.eigsy(Q, eigvals_only=True))
            tol = mp.mpf(10) ** (-args.dps + 15) * max(abs(t) for t in vals)
            nneg = sum(1 for t in vals if t < -tol)
            ts = " ".join(f"{int(mp.nint(mp.exp(l)))}:{mp.nstr(a, 3)}" for l, a in terms)
            print(f"{xs:>5} {fn:>8}  {mp.nstr(vals[0], 5):>12} {mp.nstr(vals[1], 5):>12} {nneg:>3}   {ts}", flush=True)
            rows.append({"x": xs, "function": fn, "n": args.n, "dps": args.dps, "eps": [mp.nstr(t, 15) for t in vals[:4]], "n_negative": nneg})
            job.result(f"x={xs} {fn} eps0={mp.nstr(vals[0], 6)} n_neg={nneg}")
            job.step(f"x={xs} {fn}")
        if args.json:
            with open(args.json, "w") as f:
                json.dump(rows, f, indent=2)
    job.done()


if __name__ == "__main__":
    main()
