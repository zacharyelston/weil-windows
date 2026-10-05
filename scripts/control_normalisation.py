#!/usr/bin/env python3
"""Normalisation checks: the Q of docs/CONTROL_CERTIFICATES.md against our zeros-side matrices.

For each (function, sector, x):
  (a) the T₀ = 1/√ℓ entry (even) or the b_1 = √(2/ℓ) sin(2πu/ℓ) entry (odd): rigorous u-space value
      (control_uspace.py) against our matrix entry;
  (b) one test function vanishing to 4th order at ±ℓ/2,
        even φ = (1 + cos θ)² = 3/2 + 2 cos θ + ½ cos 2θ,  odd ψ = sin θ (1 + cos θ)² = 5/4 sin θ + sin 2θ + 1/4 sin 3θ,
      θ = 2πu/ℓ: rigorous u-space value, rigorous frequency-side value (acb.integral of Ψ|F|² on [0, T] plus a
      tail bound from |F(t)| ≤ ‖f⁗‖₁/t⁴), and our cᵀ M c.
native matrices (mpmath, --dps digits):
  ftstar: conductor5_family_mp.forms(t*, x, n, parity)[2]   (E + 2vvᵀ even, E − 2wwᵀ odd)
  dh:     connes_letter_mp.build_form(x, n, "dh", parity)    (no pole)
  z1:     even: epstein_connes_mp.forms("Z1", x, n)[1]; odd: assembled the same way with parity="odd", minus 2wwᵀ
  zeta:   connes_letter_mp.build_form(x, n) + 2vvᵀ (even) / − 2wwᵀ (odd)

Usage: .venv/bin/python scripts/control_normalisation.py --cases ftstar:even:7,ftstar:odd:7,dh:even:32,z1:even:20 \
         --json data/controls/normalisation.json
"""

import argparse
import json
import os
import sys
import time

import mpmath as mp
from flint import arb, ctx

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import conductor5_family_mp as c5  # noqa: E402
import connes_letter_mp as cl  # noqa: E402
import epstein_connes_mp as ec  # noqa: E402
from control_cert_lib import Control  # noqa: E402
from control_uspace import Q_uspace, TrigPoly, freq_side  # noqa: E402
from progress import Job  # noqa: E402


def native_matrix(function, x, n, parity):
    x = mp.mpf(x)
    L = mp.log(x)
    if function == "ftstar":
        return c5.forms(c5.tstar(), x, n, parity)[2]
    if function == "dh":
        return cl.build_form(x, n, "dh", parity=parity)[1]
    if function == "zeta":
        _, E, _ = cl.build_form(x, n, "zeta", parity=parity)
        if parity == "even":
            v = mp.matrix([cl.int_cos_cosh(2 * mp.pi * k / L, L) * (1 / mp.sqrt(L) if k == 0 else mp.sqrt(2 / L)) for k in range(n + 1)])
            return E + 2 * v * v.T
        w = cl.odd_pole_vector(L, n)
        return E - 2 * w * w.T
    if function in ("z1", "zetak"):
        name = "Z1" if function == "z1" else "zeta_K"
        if parity == "even":
            return ec.forms(name, x, n)[1]
        n_max = int(mp.floor(x))
        a = ec.coefficients(name, max(n_max, 2))
        c = ec.log_derivative(a, max(n_max, 2))
        terms = [(mp.log(k), c[k] / mp.sqrt(k)) for k in range(2, n_max + 1) if c[k] != 0]
        _, Qz, _ = cl.build_form(x, n, "zeta", terms=[], parity="odd")
        _, Qd, _ = cl.build_form(x, n, "dh", terms=[], parity="odd")
        E = Qz + Qd + mp.log(4) * mp.eye(n)
        if terms:
            E = E + cl.build_form(x, n, "zeta", terms=terms, include_arch=False, parity="odd")[1]
        w = cl.odd_pole_vector(L, n)
        return E - 2 * w * w.T
    raise ValueError(function)


def test_amplitudes(parity):
    """Unnormalised amplitudes (A_0, A_1, …) of the 4th-order-vanishing test function."""
    if parity == "even":
        return [mp.mpf(3) / 2, mp.mpf(2), mp.mpf(1) / 2]
    return [mp.mpf(0), mp.mpf(5) / 4, mp.mpf(1), mp.mpf(1) / 4]


def to_native_coeffs(A, L, parity):
    """native orthonormal coefficients from unnormalised amplitudes (even includes x_0; odd starts at k = 1)."""
    if parity == "even":
        return [A[0] * mp.sqrt(L)] + [A[k] * mp.sqrt(L / 2) for k in range(1, len(A))]
    return [A[k] * mp.sqrt(L / 2) for k in range(1, len(A))]


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--cases", default="ftstar:even:7,ftstar:odd:7,dh:even:32,z1:even:20")
    ap.add_argument("--n", type=int, default=8, help="native basis size (exact for these finite trigonometric polynomials)")
    ap.add_argument("--dps", type=int, default=40)
    ap.add_argument("--prec", type=int, default=192)
    ap.add_argument("--tfreq", type=int, default=2000)
    ap.add_argument("--json")
    args = ap.parse_args()
    mp.mp.dps = args.dps
    cases = [c.split(":") for c in args.cases.split(",")]
    job = Job("control normalisation", total=len(cases), args=vars(args))
    rows = []
    for fn, parity, xs in cases:
        t0 = time.time()
        ctx.prec = args.prec
        ctrl = Control(fn)
        if xs.startswith("e^"):
            ell_a, x_int, x_mp = arb(xs[2:]), None, mp.exp(mp.mpf(xs[2:]))
        else:
            ell_a, x_int, x_mp = arb(xs).log(), (int(xs) if xs.isdigit() else None), mp.mpf(xs)
        L = mp.log(x_mp)
        M = native_matrix(fn, x_mp, args.n, parity)
        # (a) T0 / b1 entry
        dim = args.n + 1 if parity == "even" else args.n
        e0 = [arb(1)] + [arb(0)] * (dim - 1)
        ent = Q_uspace(ctrl, ell_a, parity, e0, args.prec, x_int)
        native_ent = M[0, 0]
        # (b) test function
        A = test_amplitudes(parity)
        cvec = to_native_coeffs(A, L, parity) + [mp.mpf(0)] * (dim - (len(A) if parity == "even" else len(A) - 1))
        cv = mp.matrix(cvec)
        native_q = (cv.T * M * cv)[0, 0]
        norm_mp = sum(c * c for c in cvec)
        xs_arb = [arb(mp.nstr(c, args.dps)) for c in cvec]
        uq = Q_uspace(ctrl, ell_a, parity, xs_arb, args.prec, x_int)
        tp = uq["tp"]
        Dk = ell_a * sum((abs(tp.A[k]) * tp.om[k] ** 4 for k in range(1, tp.n + 1)), arb(0))
        ctx.prec = 128
        fq, tail = freq_side(ctrl, tp, args.tfreq, 128, deriv_bound=Dk, kbound=4, x_int=x_int)
        fq_enc = fq + arb(0, tail.upper())
        ctx.prec = args.prec
        diff_native = abs(uq["Q"] - arb(mp.nstr(native_q, args.dps)))
        overlap = fq_enc.overlaps(uq["Q"])
        row = {"function": fn, "parity": parity, "x": xs, "native_n": args.n, "native_dps": args.dps,
               "entry": {"uspace": ent["Q"].str(30, radius=True), "native": mp.nstr(native_ent, 30),
                         "abs_diff": mp.nstr(abs(mp.mpf(ent["Q"].mid().str(40, radius=False)) - native_ent), 3)},
               "test_function": {"uspace": uq["Q"].str(30, radius=True), "freq_side": fq_enc.str(20, radius=True),
                                 "freq_tail_bound": tail.str(3), "T": args.tfreq, "native": mp.nstr(native_q, 30),
                                 "norm2_uspace": uq["norm2"].str(20), "norm2_native": mp.nstr(norm_mp, 20),
                                 "uspace_vs_native_abs_diff": diff_native.upper().str(3), "freq_overlaps_uspace": bool(overlap)},
               "prime_n": uq["n_terms"], "elapsed_s": round(time.time() - t0, 1)}
        rows.append(row)
        job.result(f"{fn} {parity} x={xs}: entry uspace {ent['Q'].str(20)} vs native {mp.nstr(native_ent, 20)}; "
                   f"test fn uspace {uq['Q'].str(20)}, freq {fq_enc.str(12)} (tail <= {tail.str(2)}), native {mp.nstr(native_q, 20)}; "
                   f"|uspace - native| <= {diff_native.upper().str(3)}; freq overlaps uspace: {overlap}")
        job.step(f"{fn} {parity} x={xs}")
        if args.json:
            os.makedirs(os.path.dirname(os.path.abspath(args.json)), exist_ok=True)
            with open(args.json, "w") as fh:
                json.dump(rows, fh, indent=2)
    job.done()


if __name__ == "__main__":
    main()
