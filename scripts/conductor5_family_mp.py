#!/usr/bin/env python3
"""Even conductor-5 family switching the Euler product off with the functional
equation fixed (docs/LITERATURE_PASS.md, Q5 and recommended next shot).

  F_t(s) = (1 − t) L(s, χ₅) + t (1 + √5·5^{−s}) ζ(s),   χ₅ = (·/5),
  Λ_t(s) = (5/π)^{s/2} Γ(s/2) F_t(s) = Λ_t(1 − s).

t = 0 (L(s, χ₅)) and t = 1 ((1 + √5·5^{−s})ζ) have Euler products; 0 < t < 1 does
not (a₂a₃ ≠ a₆). For t > 0 there is a simple pole at s = 1. The choice
t* = L(¾, χ₅)/(L(¾, χ₅) − G(¾)), G = (1 + √5·5^{−s})ζ, puts real zeros at
s = ¾ and s = ¼, i.e. γ = ∓i/4: an off-line pair that only ODD test functions
see (pair term −2φ̂(i/4)² for odd φ, +2φ̂(i/4)² for even φ).

Subcommands:
  check     functional equation for several t; F_{t*}(¾) = F_{t*}(¼) = 0
  census    zeros of F_t on and off the line (γ ∈ [t1, t2]) plus real zeros
  validate  explicit formula in both parity sectors against the census zeros
  sweep     smallest eigenvalues of the zeros-side form, even and odd sectors

Usage:
  .venv/bin/python scripts/conductor5_family_mp.py check
  .venv/bin/python scripts/conductor5_family_mp.py census --t tstar --t2 40 --json data/conductor5/census_tstar_40.json
  .venv/bin/python scripts/conductor5_family_mp.py validate --zeros data/conductor5/census_tstar_40.json --x 5
  .venv/bin/python scripts/conductor5_family_mp.py sweep --t 0,tstar,1 --x 2,3,5,7,10 --n 32 --dps 40
"""

import argparse
import json
import os
import sys

import mpmath as mp

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import connes_letter_mp as cl  # noqa: E402
import epstein_mp as ep  # noqa: E402  (generic argument-principle and Newton helpers)
from progress import Job  # noqa: E402

CHI5 = [0, 1, -1, -1, 1]


def L5(s):
    return mp.dirichlet(s, CHI5)


def G(s):
    return (1 + mp.sqrt(5) * mp.power(5, -s)) * mp.zeta(s)


def tstar(sigma=mp.mpf(3) / 4):
    return L5(sigma) / (L5(sigma) - G(sigma))


def parse_t(v):
    return tstar() if v == "tstar" else mp.mpf(v)


def F(t):
    return lambda s: (1 - t) * L5(s) + t * G(s)


def completed(t, s):
    return mp.exp(s / 2 * mp.log(5 / mp.pi) + mp.loggamma(s / 2)) * F(t)(s)


def theta(t_im):
    s = mp.mpc(0.5, t_im)
    return mp.im(s / 2 * mp.log(5 / mp.pi) + mp.loggamma(s / 2))


def hardy(f, t_im):
    return mp.re(mp.exp(1j * theta(t_im)) * f(mp.mpc(0.5, t_im)))


def coefficients(t, n_max):
    r5 = mp.sqrt(5)
    return [mp.mpf(0)] + [(1 - t) * CHI5[n % 5] + t * (1 + (r5 if n % 5 == 0 else 0)) for n in range(1, n_max + 1)]


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


def forms(t, x, n, parity):
    """Zeros-side form of F_t on the given parity sector (window [λ⁻¹, λ], x = λ²)."""
    x = mp.mpf(x)
    n_max = int(mp.floor(x))
    a = coefficients(t, max(n_max, 2))
    c = log_derivative(a, max(n_max, 2))
    terms = [(mp.log(k), c[k] / mp.sqrt(k)) for k in range(2, n_max + 1) if c[k] != 0]
    L, Qa, _ = cl.build_form(x, n, "zeta", terms=[], parity=parity)
    dim = Qa.rows
    E = Qa + mp.log(5) * mp.eye(dim)
    if terms:
        _, Qp, _ = cl.build_form(x, n, "zeta", terms=terms, include_arch=False, parity=parity)
        E = E + Qp
    if t == 0:
        return L, E, E  # entire: no pole correction
    if parity == "even":
        v = mp.matrix([cl.int_cos_cosh(2 * mp.pi * k / L, L) * (1 / mp.sqrt(L) if k == 0 else mp.sqrt(2 / L)) for k in range(n + 1)])
        return L, E, E + 2 * v * v.T
    w = cl.odd_pole_vector(L, n)
    return L, E, E - 2 * w * w.T


def gamma_values(zeros):
    """All γ = (ρ − 1/2)/i for the zero set closed under ρ ↦ 1 − ρ, ρ ↦ ρ̄."""
    out = []
    for g in zeros["on_line"]:
        out += [mp.mpf(g), -mp.mpf(g)]
    for z in zeros["off_line"]:
        d, g = mp.mpf(z["beta"]) - mp.mpf(0.5), mp.mpf(z["gamma"])
        out += [mp.mpc(g, -d), mp.mpc(g, d), mp.mpc(-g, -d), mp.mpc(-g, d)]
    for b in zeros.get("real", []):
        d = mp.mpf(b) - mp.mpf(0.5)
        if d > 0:
            out += [mp.mpc(0, -d), mp.mpc(0, d)]
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("check")
    c = sub.add_parser("census")
    c.add_argument("--t", default="tstar")
    c.add_argument("--t1", type=float, default=0.5)
    c.add_argument("--t2", type=float, default=40)
    c.add_argument("--dps", type=int, default=20)
    c.add_argument("--json")
    v = sub.add_parser("validate")
    v.add_argument("--t", default="tstar")
    v.add_argument("--zeros", required=True)
    v.add_argument("--x", type=float, default=5)
    v.add_argument("--dps", type=int, default=30)
    s = sub.add_parser("sweep")
    s.add_argument("--t", default="0,tstar,1")
    s.add_argument("--x", default="2,3,5,7,10")
    s.add_argument("--n", type=int, default=32)
    s.add_argument("--dps", type=int, default=40)
    s.add_argument("--json")
    args = ap.parse_args()

    if args.cmd == "check":
        mp.mp.dps = 30
        ts = tstar()
        print(f"t* = {mp.nstr(ts, 20)}   (L(¾,χ₅) = {mp.nstr(L5(mp.mpf(0.75)), 12)}, G(¾) = {mp.nstr(G(mp.mpf(0.75)), 12)})")
        for t in [mp.mpf(0), ts, mp.mpf("0.5"), mp.mpf(1)]:
            errs = [abs(completed(t, s_) - completed(t, 1 - s_)) / abs(completed(t, s_)) for s_ in [mp.mpc(0.3, 6), mp.mpc(0.9, 23), mp.mpc(-0.4, 41)]]
            print(f"  t = {mp.nstr(t, 8):<12} |Λ(s) − Λ(1−s)|/|Λ(s)| = " + "  ".join(mp.nstr(e, 3) for e in errs))
        f = F(ts)
        print(f"  F_t*(3/4) = {mp.nstr(f(mp.mpf(0.75)), 3)},  F_t*(1/4) = {mp.nstr(f(mp.mpf(0.25)), 3)}")
        grid = [mp.mpf(k) / 200 for k in range(1, 200)]
        signs = [mp.sign(f(g)) for g in grid]
        changes = [mp.nstr(grid[i], 3) for i in range(1, len(grid)) if signs[i] != signs[i - 1]]
        print(f"  sign changes of F_t* on (0, 1): {changes}")
        a = coefficients(ts, 12)
        print(f"  a_n(t*), n = 1..12: " + " ".join(mp.nstr(v_, 4) for v_ in a[1:]) + f"   a₂a₃ − a₆ = {mp.nstr(a[2] * a[3] - a[6], 4)}")
        return

    if args.cmd == "census":
        mp.mp.dps = args.dps
        t = parse_t(args.t)
        f = F(t)
        job = Job("conductor5 census", total=4, args=vars(args))
        z = lambda ti: hardy(f, ti)
        on = []
        a, za = mp.mpf(args.t1), z(args.t1)
        n_steps, k = int((args.t2 - args.t1) / 0.02) + 1, 0
        while a < args.t2:
            k += 1
            job.sub(f"on-line scan at γ={mp.nstr(a, 5)} ({len(on)} zeros)", k, n_steps, every=max(1, n_steps // 20))
            b = a + mp.mpf("0.02")
            zb = z(b)
            if za * zb < 0:
                on.append(mp.findroot(z, (a, b), solver="bisect", tol=mp.mpf(10) ** (-mp.mp.dps + 4)))
            a, za = b, zb
        job.step(f"on-line zeros: {len(on)}")
        off = ep.off_line_zeros(f, args.t1, args.t2, job=job)
        job.step(f"off-line pairs: {len(off)}")
        total = ep.zero_count_rect(f, -2, 3, args.t1, args.t2)
        job.step(f"argument principle total: {total}")
        grid = [(mp.mpf(k) + mp.mpf("0.5")) / 400 for k in range(0, 400)]
        real = [mp.findroot(f, (grid[i - 1], grid[i]), solver="bisect") for i in range(1, len(grid)) if mp.sign(f(grid[i])) != mp.sign(f(grid[i - 1]))]
        job.step(f"real zeros: {len(real)}")
        job.result(f"total={total} on={len(on)} off_pairs={len(off)} real={[mp.nstr(r, 8) for r in real]}")
        job.done()
        print(f"F_t, t = {mp.nstr(t, 10)}: γ ∈ [{args.t1}, {args.t2}]: argument principle {total}; on the line {len(on)}; off-line pairs {len(off)} → {len(on) + 2 * len(off)}")
        print(f"  real zeros in (0, 1): {[mp.nstr(r, 12) for r in real]}")
        for zz, r in off:
            print(f"  off-line β = {mp.nstr(zz.real, 8)}  γ = {mp.nstr(zz.imag, 10)}")
        if args.json:
            os.makedirs(os.path.dirname(os.path.abspath(args.json)), exist_ok=True)
            with open(args.json, "w") as fh:
                json.dump({"t": mp.nstr(t, 25), "t1": args.t1, "t2": args.t2, "total": total, "on_line": [float(g) for g in on],
                           "off_line": [{"beta": float(zz.real), "gamma": float(zz.imag)} for zz, _ in off],
                           "real": [float(r) for r in real]}, fh, indent=2)
        return

    if args.cmd == "validate":
        mp.mp.dps = args.dps
        t = parse_t(args.t)
        with open(args.zeros) as fh:
            zs = json.load(fh)
        gammas = gamma_values(zs)
        x = mp.mpf(args.x)
        for parity in ["even", "odd"]:
            nb = 8
            L, E, QW = forms(t, x, nb, parity)
            if parity == "even":
                coeff = [mp.mpf(3) / 2 * mp.sqrt(L), 2 * mp.sqrt(L / 2), mp.sqrt(L / 2) / 2] + [mp.mpf(0)] * (nb - 2)
                ft = lambda z: cl.fourier_even(coeff, z, L) if not isinstance(z, mp.mpc) else ep_fourier_even_complex(coeff, z, L)
            else:
                coeff = [mp.sqrt(L / 2), mp.sqrt(L / 2) / 2] + [mp.mpf(0)] * (nb - 2)
                ft = lambda z: cl.fourier_odd(coeff, z, L)
            cv = mp.matrix(coeff)
            q = (cv.T * QW * cv)[0, 0]
            zsum = mp.fsum(ft(g) * mp.conj(ft(mp.conj(g))) for g in gammas)
            real_pair = mp.fsum(ft(g) * mp.conj(ft(mp.conj(g))) for g in gammas if mp.re(g) == 0)
            print(f"  {parity:<4} sector, x = {args.x}: zeros-side form = {mp.nstr(q, 12)};  Σ over census zeros = {mp.nstr(mp.re(zsum), 12)}  "
                  f"(|Im| {mp.nstr(abs(mp.im(zsum)), 2)}; real pair contributes {mp.nstr(mp.re(real_pair), 6)})")
        return

    mp.mp.dps = args.dps
    rows = []
    t_list, x_list = args.t.split(","), args.x.split(",")
    job = Job("conductor5 sweep", total=2 * len(t_list) * len(x_list), args=vars(args))
    print(f"{'t':>12} {'x':>5}   {'ε₀ even':>12} {'ε₀ odd':>12}  n₋(even, odd)")
    for tv in t_list:
        t = parse_t(tv)
        for xs in x_list:
            res = {}
            for parity in ["even", "odd"]:
                _, _, QW = forms(t, xs, args.n, parity)
                vals = sorted(mp.eigsy(QW, eigvals_only=True))
                tol = mp.mpf(10) ** (-args.dps + 12) * max(abs(u) for u in vals)
                res[parity] = (vals[0], sum(1 for u in vals if u < -tol))
                job.step(f"t={tv} x={xs} parity={parity} eps0={mp.nstr(vals[0], 5)}")
            print(f"{mp.nstr(t, 8):>12} {xs:>5}   {mp.nstr(res['even'][0], 5):>12} {mp.nstr(res['odd'][0], 5):>12}  ({res['even'][1]}, {res['odd'][1]})", flush=True)
            job.result(f"t={tv} x={xs} eps_even={mp.nstr(res['even'][0], 6)} eps_odd={mp.nstr(res['odd'][0], 6)} n_neg=({res['even'][1]},{res['odd'][1]})")
            rows.append({"t": mp.nstr(t, 20), "x": xs, "n": args.n, "dps": args.dps,
                         "eps_even": mp.nstr(res["even"][0], 15), "eps_odd": mp.nstr(res["odd"][0], 15),
                         "neg_even": res["even"][1], "neg_odd": res["odd"][1]})
            if args.json:
                os.makedirs(os.path.dirname(os.path.abspath(args.json)), exist_ok=True)
                with open(args.json, "w") as fh:
                    json.dump(rows, fh, indent=2)
    job.done()


def ep_fourier_even_complex(c, z, L):
    n = len(c) - 1
    sinc = lambda w: L / 2 if abs(w) < mp.mpf(10) ** (-mp.mp.dps + 5) else mp.sin(w * L / 2) / w
    val = c[0] * 2 * sinc(z) / mp.sqrt(L)
    for k in range(1, n + 1):
        om = 2 * mp.pi * k / L
        val += c[k] * mp.sqrt(2 / L) * (sinc(z - om) + sinc(z + om))
    return val


if __name__ == "__main__":
    main()
