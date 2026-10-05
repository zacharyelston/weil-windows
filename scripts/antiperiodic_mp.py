#!/usr/bin/env python3
"""The "720°" basis: Weil's form on half-integer (antiperiodic) frequencies.

rh2's bases on the window [−L/2, L/2] use the periodic frequencies ω_k = 2πk/L: cos(ω_k u), with zero
slope at the ends, and sin(ω_k u), which vanishes there. The antiperiodic frequencies
ω_k = 2π(k + ½)/L close only after two trips across the window:
  - even sector: cos(ω_k u) vanishes at u = ±L/2 (Dirichlet), like ζ's near-minimiser;
  - odd sector:  sin(ω_k u) has zero slope at the ends (Neumann).

The entries of Q in the exponential basis e^{iωu}/√L depend on the frequencies only through transforms
of the archimedean and prime kernels at ω_j and ω_k, and through the boundary phase
e^{i(ω_j − ω_k)L/2}. That phase is still (−1)^{j−k} on the shifted lattice, so connes_letter_mp's
closed forms carry over unchanged (`build_form_general`, with shift = 0 or ½).

Commands:
  check   shift = 0 reproduces connes_letter_mp.build_form; shift = ½ matches Σ_ρ |F(γ)|² over ζ zeros.
  sweep   λ_min of the zeros-side form against N, periodic versus antiperiodic, both sectors, several a.

Usage:
  .venv/bin/python scripts/antiperiodic_mp.py check
  .venv/bin/python scripts/antiperiodic_mp.py sweep --a 0.8,1.0,1.19,1.3,1.495 --n 16,32,48,64,100 --json data/connes/antiperiodic_sweep.json
"""

import argparse
import json
import os
import sys
import time

import mpmath as mp

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import connes_letter_mp as cl  # noqa: E402
from progress import Job  # noqa: E402


def build_form_general(lambda_sq, n, parity="even", shift=0, terms=None):
    """ζ's geometric-side form E on the basis √(2/L) cos(ω_k u) (even) or √(2/L) sin(ω_k u) (odd),
    ω_k = 2π(k + shift)/L. For shift = 0, the even basis also has b_0 = 1/√L and the odd basis starts
    at k = 1, as in connes_letter_mp.build_form. Returns (L, E, omegas)."""
    L = mp.log(lambda_sq)
    quarter = mp.mpf(1) / 4
    tail_terms = int(mp.mp.dps * mp.log(10) / (2 * L)) + 10
    shift = mp.mpf(shift)
    # The archimedean closed forms truncate ∫_0^∞ at x = L, where e^{−zL} = e^{−aL}·e^{iωL}. On both lattices
    # sin(ωL) = 0 and cos(ωL) = cL = +1 (periodic) or −1 (antiperiodic), so every E = e^{−aL} becomes cL·E.
    assert shift in (0, mp.mpf(1) / 2), "closed forms need sin(ωL) = 0"
    cL = 1 if shift == 0 else -1
    if shift == 0:
        ks = list(range(0, n + 1)) if parity == "even" else list(range(1, n + 1))
    else:
        ks = list(range(0, n + 1))
    om = [2 * mp.pi * (k + shift) / L for k in ks]

    def j_term(m, w):
        a = 2 * m + mp.mpf(1) / 2
        E = cL * mp.exp(-a * L)
        d = a * a + w * w
        c = L * a * (1 - E) / d + (((1 - E) + a * L * E) * d - 2 * a * a * (1 - E)) / (d * d)
        return 2 * c - 2 * L * (1 - mp.exp(-(a + mp.mpf(1) / 2) * L)) / (a + mp.mpf(1) / 2)

    def I_of(w):          # odd in w
        if w == 0:
            return mp.mpf(0)
        z = mp.mpc(quarter, abs(w) / 2)
        main = mp.im(mp.digamma(z)) / 2
        tail = cL * mp.fsum(abs(w) * mp.exp(-(2 * m + mp.mpf(1) / 2) * L) / ((2 * m + mp.mpf(1) / 2) ** 2 + w * w) for m in range(tail_terms))
        return (main - tail) * (1 if w > 0 else -1)

    def J_of(w):          # even in w
        w = abs(w)
        z = mp.mpc(quarter, w / 2)
        main = L * (mp.digamma(mp.mpf(1) / 2) - mp.re(mp.digamma(z))) - mp.re(mp.psi(1, z)) / 2

        def t0(m):
            a = 2 * m + mp.mpf(1) / 2
            d = a * a + w * w
            return 2 * L * a / d + 2 * (w * w - a * a) / (d * d) - 2 * L / (a + mp.mpf(1) / 2)

        return main + mp.fsum(j_term(m, w) - t0(m) for m in range(tail_terms))

    pp = list(terms) if terms is not None else [(mp.log(nn), lp / mp.sqrt(nn)) for nn, lp in cl.prime_powers(int(lambda_sq))]
    Ivals = [I_of(w) for w in om]
    Jvals = [J_of(w) for w in om]
    Pvals = [mp.fsum(a * mp.sin(w * ln) for ln, a in pp) for w in om]          # odd in w
    Dvals = [mp.fsum(a * 2 * (L - ln) * mp.cos(w * ln) for ln, a in pp) for w in om]   # even in w
    diag_const = mp.log(4 * mp.pi) + mp.euler + mp.log(mp.tanh(L / 2))

    # Exponential basis indexed by (i, s): frequency s·om[i], s = ±1.
    def q(i, s, j, t):
        w1, w2 = s * om[i], t * om[j]
        if i == j and s == t:
            return (-(diag_const * L + Jvals[i]) - Dvals[i]) / L
        delta = w1 - w2
        nint = int(mp.nint(delta * L / (2 * mp.pi)))
        sign = 1 if nint % 2 == 0 else -1
        I1, I2 = s * Ivals[i], t * Ivals[j]
        P1, P2 = s * Pvals[i], t * Pvals[j]
        return -2 * sign * ((I2 - I1) + (P2 - P1)) / delta / L

    m = len(ks)
    E = mp.matrix(m, m)
    s2 = mp.sqrt(2)
    for a in range(m):
        for b in range(a, m):
            if parity == "even":
                if shift == 0 and ks[a] == 0 and ks[b] == 0:
                    v = q(a, 1, b, 1)
                elif shift == 0 and ks[a] == 0:
                    v = (q(a, 1, b, 1) + q(a, 1, b, -1)) / s2
                else:
                    v = (q(a, 1, b, 1) + q(a, 1, b, -1) + q(a, -1, b, 1) + q(a, -1, b, -1)) / 2
            else:
                v = (q(a, 1, b, 1) - q(a, 1, b, -1) - q(a, -1, b, 1) + q(a, -1, b, -1)) / 2
            E[a, b] = v
            E[b, a] = v
    return L, E, om


def zeros_side(lambda_sq, n, parity, shift):
    """Zeros-side form: E + 2vvᵀ (even, v_k = ∫ b_k cosh(u/2)) or E − 2wwᵀ (odd, w_k = ∫ b_k sinh(u/2))."""
    L, E, om = build_form_general(lambda_sq, n, parity, shift)
    if parity == "even":
        v = mp.matrix([cl.int_cos_cosh(w, L) * (1 / mp.sqrt(L) if w == 0 else mp.sqrt(2 / L)) for w in om])
        return L, E + 2 * v * v.T, om
    w = mp.matrix([mp.sqrt(2 / L) * cl.int_sin_sinh(x, L) for x in om])
    return L, E - 2 * w * w.T, om


def F_hat(om, parity, L, t):
    """F(t) = ∫ b(u) e^{itu} du for b = √(2/L) cos(ωu) (real) or √(2/L) sin(ωu) (imaginary part returned)."""
    def S(x):
        return L / 2 if abs(x) < mp.mpf(10) ** (-mp.mp.dps + 5) else mp.sin(x * L / 2) / x
    if parity == "even":
        return mp.sqrt(2 / L) * (S(t - om) + S(t + om))
    return mp.sqrt(2 / L) * (S(t - om) - S(t + om))


def cmd_check(args):
    mp.mp.dps = 40
    job = Job("antiperiodic check", args=vars(args))
    x = mp.mpf(13)
    for parity in ["even", "odd"]:
        _, E0, _ = build_form_general(x, 12, parity, 0)
        _, Eref, _ = cl.build_form(x, 12, "zeta", parity=parity)
        d = max(abs(E0[i, j] - Eref[i, j]) for i in range(E0.rows) for j in range(E0.cols))
        job.result(f"shift 0 vs build_form, {parity}, x = 13, N = 12: max |Δ| = {mp.nstr(d, 3)}")
    # independent: zeros-side Σ_ρ |F(γ)|² for single antiperiodic basis functions (they vanish or are flat at
    # the edges, so F ~ 1/t² or 1/t and the zero sum converges)
    nz = args.zeros
    gam = [mp.im(mp.zetazero(k)) for k in range(1, nz + 1)]
    job.log(f"{nz} zeros, γ_max = {mp.nstr(gam[-1], 8)}")
    for xs in ["5", "13"]:
        x = mp.mpf(xs)
        # Sharp checks are the edge-vanishing functions (F ~ 1/t², tail ~ 1/T²): antiperiodic even, periodic odd.
        # The others jump at the edges (F ~ 1/t), so their truncated zero sums converge only like log T/T.
        for shift_name, shift in [("periodic", 0), ("antiperiodic", mp.mpf(1) / 2)]:
            for parity in ["even", "odd"]:
                L, Q, om = zeros_side(x, 6, parity, shift)
                for i in (1, 3):
                    Fz = mp.fsum(2 * F_hat(om[i], parity, L, g) ** 2 for g in gam)
                    sharp = (shift_name, parity) in [("antiperiodic", "even"), ("periodic", "odd")]
                    job.result(f"x = {xs}, {shift_name:<12} {parity:<4} ω = {mp.nstr(om[i] * L / (2 * mp.pi), 3)}·2π/L: Q_ii = {mp.nstr(Q[i, i], 14)};  "
                               f"2Σ_ρ F(γ)² ({nz} zeros) = {mp.nstr(Fz, 14)};  rel diff {mp.nstr(abs(Q[i, i] - Fz) / abs(Q[i, i]), 3)}"
                               + ("  [sharp]" if sharp else ""))
        # Off-diagonal entries: edge-vanishing combinations. Antiperiodic odd b_0 + b_1 (sin(π(2k+1)/2) = (−1)^k
        # cancels at the edges); antiperiodic even b_1 − 2b_2 + b_3 (generic, already edge-vanishing).
        for parity, c in [("odd", [1, 1, 0, 0]), ("even", [0, 1, -2, 1])]:
            L, Q, om = zeros_side(x, 6, parity, mp.mpf(1) / 2)
            cv = mp.matrix(c + [0] * (Q.rows - len(c)))
            qv = (cv.T * Q * cv)[0]
            Fz = mp.fsum(2 * mp.fsum(ci * F_hat(om[i], parity, L, g) for i, ci in enumerate(c)) ** 2 for g in gam)
            job.result(f"x = {xs}, antiperiodic {parity:<4} c = {c}: cᵀQc = {mp.nstr(qv, 14)};  2Σ_ρ F(γ)² = {mp.nstr(Fz, 14)};  "
                       f"rel diff {mp.nstr(abs(qv - Fz) / abs(qv), 3)}  [sharp]")
    job.done()


def cmd_sweep(args):
    rows = []
    a_list, n_list = args.a.split(","), [int(s) for s in args.n.split(",")]
    job = Job("antiperiodic sweep", total=len(a_list) * len(n_list) * 4, args=vars(args))
    for a_s in a_list:
        a = mp.mpf(a_s)
        mp.mp.dps = args.dps if float(a_s) < 1.4 else args.dps_high
        x = mp.exp(2 * a)
        for parity in ["even", "odd"]:
            for shift_name, shift in [("periodic", 0), ("antiperiodic", mp.mpf(1) / 2)]:
                vals = []
                for n in n_list:
                    t0 = time.time()
                    _, Q, _ = zeros_side(x, n, parity, shift)
                    lam = min(mp.eigsy(Q, eigvals_only=True))
                    vals.append(lam)
                    job.step(f"a={a_s} {parity} {shift_name} N={n} λ={mp.nstr(lam, 6)} ({time.time() - t0:.0f}s)")
                    rows.append({"a": a_s, "x": mp.nstr(x, 12), "parity": parity, "basis": shift_name, "n": n,
                                 "lambda_min": mp.nstr(lam, 20), "dps": mp.mp.dps})
                job.result(f"a={a_s} {parity:<4} {shift_name:<12} " + "  ".join(f"N={n}:{mp.nstr(v, 6)}" for n, v in zip(n_list, vals)))
                if args.json:
                    os.makedirs(os.path.dirname(os.path.abspath(args.json)), exist_ok=True)
                    with open(args.json, "w") as fh:
                        json.dump(rows, fh, indent=2)
    job.done()


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    c = sub.add_parser("check")
    c.add_argument("--zeros", type=int, default=200)
    s = sub.add_parser("sweep")
    s.add_argument("--a", default="0.8,1.0,1.19,1.3,1.495")
    s.add_argument("--n", default="16,32,48,64,100")
    s.add_argument("--dps", type=int, default=120)
    s.add_argument("--dps-high", type=int, default=220)
    s.add_argument("--json")
    args = ap.parse_args()
    {"check": cmd_check, "sweep": cmd_sweep}[args.cmd](args)


if __name__ == "__main__":
    main()
