#!/usr/bin/env python3
"""Rigorous (Arb) evaluation of Q(f) for a trigonometric polynomial f on [−ℓ/2, ℓ/2], in u-space.

f = Σ x_k b_k in rh2's orthonormal bases on [−ℓ/2, ℓ/2], ℓ = log x, ω_k = 2πk/ℓ:
  even: b_0 = 1/√ℓ, b_k = √(2/ℓ) cos(ω_k u)       (connes_letter_mp / epstein_connes_mp / conductor5 even)
  odd:  b_k = √(2/ℓ) sin(ω_k u), k = 1 … n          (conductor5 odd)

Identity used (docs/CONTROL_CERTIFICATES.md §1.2). For Re κ > 0, from ψ(z) = −γ + Σ_m [1/(m+1) − 1/(m+z)],
  Re ψ(κ + it/2) = ψ(κ) + ∫_0^∞ 2 e^{−2κu} (1 − cos tu)/(1 − e^{−2u}) du,
and Tonelli (non-negative integrand) with (1/π)∫_0^∞ |F|² (1 − cos tu) dt = g(0) − g(u) give
  (1/π)∫_0^∞ Re ψ(κ + it/2)|F|² dt = g(0)[ψ(κ) + Σ_{m≥0} e^{−2(κ+m)ℓ}/(κ+m)] + ∫_0^ℓ (g(0) − g(u))/u · φ_κ(u) du,
  φ_κ(u) = 2u e^{−2κu}/(1 − e^{−2u}) = e^{(1−2κ)u} u/sinh u   (g = 0 on [ℓ, ∞)).
Hence
  Q(f) = pole(f) + g(0)[K + Σ_κ (ψ(κ) + Σ_m e^{−2(κ+m)ℓ}/(κ+m))] + ∫_0^ℓ (g(0) − g(u))/u · Σ_κ φ_κ(u) du
         − Σ_{log n < ℓ} (2c_n/√n) g(log n).
For a trigonometric polynomial g has the closed form (§1.2)
  g(u) = (ℓ − u) c(u) + Σ_m σ_m sin(ω_m u),  c(u) = A_0² + ½ Σ_k A_k² cos(ω_k u)   (A = unnormalised amplitudes),
so the integrand is entire apart from the poles of u/sinh u at ±iπ, ±2iπ, … and acb.integral encloses it.

Subcommands:
  witness   Q(f)/‖f‖² for coefficients in a JSON file ({"coefficients": [...], "x": ..., "parity": ...})
  entry     Q(b_j) (diagonal entry) or the T₀ entry
  check     u-space vs frequency-side enclosures on a C³ test function (vanishing to 4th order at ±ℓ/2)
"""

import argparse
import json
import os
import sys
import time

from flint import acb, arb, ctx

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from control_cert_lib import Control, q  # noqa: E402


def amplitudes(ell, parity, xs):
    """Unnormalised amplitudes: even f = A_0 + Σ_{k≥1} A_k cos ω_k u; odd f = Σ_{k≥1} A_k sin ω_k u.
    xs: rh2 coefficients (even: x_0 … x_n; odd: x_1 … x_n). Returns A with A[0] (0 for odd)."""
    r2 = (2 / ell).sqrt()
    if parity == "even":
        return [xs[0] / ell.sqrt()] + [x * r2 for x in xs[1:]]
    return [arb(0)] + [x * r2 for x in xs]


class TrigPoly:
    def __init__(self, ell, parity, xs):
        self.ell = ell
        self.parity = parity
        self.A = amplitudes(ell, parity, xs)
        n = len(self.A) - 1
        self.n = n
        tp = 2 * arb.pi()
        self.om = [tp * k / ell for k in range(n + 1)]
        om, A = self.om, self.A
        sig = [arb(0)] * (n + 1)
        for m in range(1, n + 1):
            s = arb(0)
            for k in range(1, n + 1):
                sgn = 1 if (k + m) % 2 == 0 else -1
                if parity == "even":
                    term = -1 / (om[k] + om[m])
                    if k != m:
                        term += 1 / (om[k] - om[m])
                    s += sgn * A[k] * term / 2
                else:
                    term = 1 / (om[k] + om[m])
                    if k != m:
                        term += 1 / (om[k] - om[m])
                    s += sgn * A[k] * term
            if parity == "even":
                s += -A[0] * (1 if m % 2 == 0 else -1) / om[m]
                sig[m] = 2 * A[m] * s
            else:
                sig[m] = A[m] * s
        self.sig = sig
        self.c0 = A[0] * A[0] if parity == "even" else arb(0)
        self.half_sq = [A[k] * A[k] / 2 for k in range(n + 1)]

    def c(self, u):
        return self.c0 + sum((self.half_sq[k] * (self.om[k] * u).cos() for k in range(1, self.n + 1)), arb(0) if isinstance(u, arb) else acb(0))

    def g(self, u):
        """g(u) = ∫ f(v) f(v + u) dv for real 0 ≤ u ≤ ℓ."""
        v = (self.ell - u) * self.c(u)
        for m in range(1, self.n + 1):
            v += self.sig[m] * (self.om[m] * u).sin()
        return v

    def norm2(self):
        return self.ell * (self.c0 + sum(self.half_sq[1:], arb(0)))

    def dg_over_u(self, z):
        """(g(0) − g(z))/z as an entire function (acb z)."""
        ell = self.ell
        tot = acb(0)
        for k in range(1, self.n + 1):
            w2 = self.om[k] / 2
            x = w2 * z
            tot += ell * (2 * self.half_sq[k]) * x.sin() * w2 * x.sinc()
        tot += self.c(z)
        for m in range(1, self.n + 1):
            tot -= self.sig[m] * self.om[m] * (self.om[m] * z).sinc()
        return tot

    def moment(self, s):
        """∫ f(u) e^{su} du for real s ≠ 0."""
        sh = 2 * (s * self.ell / 2).sinh()
        tot = arb(0)
        for k in range(0, self.n + 1):
            if k == 0:
                if self.parity == "even":
                    tot += self.A[0] * sh / s
                continue
            sgn = 1 if k % 2 == 0 else -1
            den = s * s + self.om[k] * self.om[k]
            tot += self.A[k] * sgn * sh * ((s if self.parity == "even" else -self.om[k]) / den)
        return tot

    def fourier(self, z):
        """F(z) = ∫ f e^{izu} du for acb z; for odd f returns F(z)/i (real on the real axis)."""
        h = self.ell / 2
        tot = acb(0)
        if self.parity == "even":
            tot += self.A[0] * self.ell * (z * h).sinc()
            for k in range(1, self.n + 1):
                tot += self.A[k] * h * (((z - self.om[k]) * h).sinc() + ((z + self.om[k]) * h).sinc())
        else:
            for k in range(1, self.n + 1):
                tot += self.A[k] * h * (((z - self.om[k]) * h).sinc() - ((z + self.om[k]) * h).sinc())
        return tot


def phi_sum(ctrl, z):
    """Σ_κ e^{(1−2κ)z} z/sinh z (acb z)."""
    inv = 1 / (acb(0, 1) * z).sinc()         # sinc(iz) = sinh(z)/z
    return sum(((1 - 2 * arb(k)) * z).exp() * inv for k in ctrl.kappas)


def const_term(ctrl, ell, prec):
    """K + Σ_κ [ψ(κ) + Σ_{m≥0} e^{−2(κ+m)ℓ}/(κ+m)], with a rigorous geometric tail."""
    tot = ctrl.K
    r = (-2 * ell).exp()
    M = int(prec * 0.7 / float(ell.lower())) + 5
    for k in ctrl.kappas:
        kk = arb(k)
        tot += kk.digamma()
        s = arb(0)
        for m in range(M):
            s += (-2 * (kk + m) * ell).exp() / (kk + m)
        tail = (-2 * (kk + M) * ell).exp() / (kk + M) / (1 - r)
        s += arb(0, tail.upper())
        tot += s
    return tot


def Q_uspace(ctrl, ell, parity, xs, prec=256, x_int=None):
    ctx.prec = prec
    tp = TrigPoly(ell, parity, xs)
    g0 = tp.norm2()
    C0 = const_term(ctrl, ell, prec)
    comb = ctrl.comb(ell, x_int)
    prime = sum((cn * tp.g(ln) for _, cn, ln in comb), arb(0))
    pole = arb(0)
    if ctrl.pole:
        pole = 2 * tp.moment(q(1, 2)) * tp.moment(q(-1, 2))

    def integrand(z, analytic):
        return tp.dg_over_u(z) * phi_sum(ctrl, z)

    b_mid = arb(ell.mid())
    t0 = time.time()
    I = acb.integral(integrand, 0, b_mid, rel_tol=arb(2) ** (-prec + 20), abs_tol=arb(2) ** (-prec + 20)).real
    # endpoint uncertainty: ∫ over [mid, ℓ] lies in h(ball ℓ)·[−rad, rad]
    h_end = integrand(acb(ell), False).real
    I += h_end * arb(0, ell.rad())
    Q = pole + g0 * C0 + I - prime
    out = {"Q": Q, "norm2": g0, "ratio": Q / g0, "pole": pole, "const": C0, "integral": I, "prime": prime,
           "n_terms": [n for n, _, _ in comb], "t_integral_s": time.time() - t0, "tp": tp}
    return out


def freq_side(ctrl, tp, T, prec=256, deriv_bound=None, kbound=4, x_int=None):
    """(1/π)∫_0^T Ψ |F|² dt + pole by acb.integral, plus a rigorous tail bound for t ≥ T using
    |F(t)| ≤ D/t^k, D = ‖f^{(k)}‖_1 (f ∈ C^{k−1} vanishing to order k at ±ℓ/2), Ψ ≤ Σ_κ log(t/2) + |K| + A + 2.
    Returns (enclosure, tail_bound)."""
    ctx.prec = prec
    comb = ctrl.comb(tp.ell, x_int)

    def Psi_c(z):
        v = acb(ctrl.K)
        for k in ctrl.kappas:
            kk = arb(k)
            v += ((acb(kk) + acb(0, 1) * z / 2).digamma() + (acb(kk) - acb(0, 1) * z / 2).digamma()) / 2
        for _, cn, ln in comb:
            v -= cn * (z * ln).cos()
        return v

    def integrand(z, analytic):
        F = tp.fourier(z)
        return Psi_c(z) * F * F

    sgn = 1 if tp.parity == "even" else -1      # |F|² = F² (even) or (F/i)² (odd)
    I = acb.integral(integrand, 0, T, rel_tol=arb(2) ** (-prec + 30), abs_tol=arb(2) ** (-prec + 30)).real / arb.pi()
    pole = 2 * tp.moment(q(1, 2)) * tp.moment(q(-1, 2)) if ctrl.pole else arb(0)
    tail = None
    if deriv_bound is not None:
        A = sum((abs(cn) for _, cn, _ in comb), arb(0))
        B = abs(ctrl.K) + A + 2
        nk = len(ctrl.kappas)
        Tt = arb(T)
        p = 2 * kbound
        # ∫_T^∞ (nk log(t/2) + B) D²/t^p dt = D² [nk (log(T/2) + 1/(p−1)) + B] / ((p−1) T^{p−1})
        tail = deriv_bound ** 2 * (nk * ((Tt / 2).log() + arb(1) / (p - 1)) + B) / ((p - 1) * Tt ** (p - 1)) / arb.pi()
    return pole + I, tail


def load_coeffs(path):
    with open(path) as fh:
        d = json.load(fh)
    return d


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    w = sub.add_parser("witness")
    w.add_argument("--function", required=True)
    w.add_argument("--json", required=True, help="witness JSON (coefficients as decimal strings, x, parity)")
    w.add_argument("--parity")
    w.add_argument("--x")
    w.add_argument("--prec", type=int, default=256)
    w.add_argument("--out")
    e = sub.add_parser("entry")
    e.add_argument("--function", required=True)
    e.add_argument("--x", required=True, help="x = e^ℓ, decimal or 'e^r' for ℓ = 2a with a = r/2")
    e.add_argument("--parity", default="even")
    e.add_argument("--k", type=int, default=0)
    e.add_argument("--prec", type=int, default=256)
    args = ap.parse_args()
    ctx.prec = args.prec

    def parse_ell(xs):
        if xs.startswith("e^"):
            return arb(xs[2:]), None
        return arb(xs).log(), (int(xs) if xs.isdigit() else None)

    if args.cmd == "entry":
        ctrl = Control(args.function)
        ell, x_int = parse_ell(args.x)
        n = max(args.k, 1)
        xs = [arb(0)] * (n + 1 if args.parity == "even" else n)
        idx = args.k if args.parity == "even" else args.k - 1
        xs[idx] = arb(1)
        r = Q_uspace(ctrl, ell, args.parity, xs, args.prec, x_int)
        print(f"{args.function} {args.parity} x={args.x} entry k={args.k}: Q = {r['Q'].str(30, radius=True)}  (norm² {r['norm2'].str(10)}) n={r['n_terms']}")
        return

    d = load_coeffs(args.json)
    ctrl = Control(args.function)
    parity = args.parity or d.get("parity", "odd")
    x = args.x or str(d["x"])
    ell, x_int = parse_ell(x)
    xs = [arb(c) for c in d["coefficients"]]
    t0 = time.time()
    r = Q_uspace(ctrl, ell, parity, xs, args.prec, x_int)
    print(f"{args.function} witness x={x} parity={parity} N={len(xs)} prec={args.prec}: "
          f"Q/‖f‖² = {r['ratio'].str(25, radius=True)}  ({time.time() - t0:.1f}s)")
    print(f"  Q = {r['Q'].str(25, radius=True)}; ‖f‖² = {r['norm2'].str(25, radius=True)}; upper bound of Q/‖f‖²: {r['ratio'].upper().str(12)}")
    print(f"  pieces: pole {r['pole'].str(20)}, const·g(0) {(r['const'] * r['norm2']).str(20)}, integral {r['integral'].str(20)}, prime {r['prime'].str(20)}; n = {r['n_terms']}")
    if args.out:
        out = {"function": args.function, "x": x, "parity": parity, "n_coeffs": len(xs), "prec_bits": args.prec,
               "Q": r["Q"].str(40, radius=True), "norm2": r["norm2"].str(40, radius=True),
               "ratio": r["ratio"].str(40, radius=True), "ratio_upper": r["ratio"].upper().str(20),
               "ratio_lower": r["ratio"].lower().str(20),
               "pieces": {k: r[k].str(30, radius=True) for k in ("pole", "const", "integral", "prime")},
               "prime_n": r["n_terms"], "source": args.json}
        os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
        with open(args.out, "w") as fh:
            json.dump(out, fh, indent=2)


if __name__ == "__main__":
    main()
