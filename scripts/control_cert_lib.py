#!/usr/bin/env python3
"""Shared definitions for the control-function certificates (docs/CONTROL_CERTIFICATES.md).

Every function F here has real Dirichlet coefficients a_n (a_1 = 1), a completed form
  Λ(s) = (archimedean factor) F(s) = Λ(1 − s),
and −F'/F = Σ c_n n^{−s}. With f supported in [−a, a] (our window [−ℓ/2, ℓ/2], ℓ = 2a = log x),
F(t) = ∫ f e^{itu} du and g = f ⋆ f̃ (g(y) = ∫ f(v) f(v + y) dv, supported in [−ℓ, ℓ]):

  Q(f) = pole(f) + (1/π) ∫_0^∞ Ψ(t) |F(t)|² dt,
  Ψ(t) = Σ_κ Re ψ(κ + it/2) + K − Σ_{log n < ℓ} (2 c_n/√n) cos(t log n),
  pole(f) = 2 F(i/2) F(−i/2) when F has a simple pole at s = 1, else 0.

  function  κ           K                    pole   a_n
  zeta      1/4         −log π               yes    1
  ftstar    1/4         log(5/π)             yes    (1 − t*) χ₅(n) + t* (1 + √5 [5 | n])
  dh        3/4         log(5/π)             no     Re((1 − iκ_DH) χ(n)),  χ(2) = i mod 5
  z1        1/4, 3/4    log 20 − 2 log π     yes    r(n)/2, r(n) = #{x² + 5y² = n}
  zetak     1/4, 3/4    log 20 − 2 log π     yes    Σ_{d|n} χ₋₂₀(d)   (ζ_K, K = ℚ(√−5): Z₁'s Euler partner)

All quantities are python-flint (Arb) balls at the current ctx.prec.
"""

import math

from flint import acb, arb, ctx, fmpq, fmpz

CHI5 = [0, 1, -1, -1, 1]                     # (·/5), real, even
DH_CHI = {0: (0, 0), 1: (1, 0), 2: (0, 1), 3: (0, -1), 4: (-1, 0)}   # χ mod 5 with χ(2) = i (odd)

FUNCTIONS = ("zeta", "ftstar", "dh", "z1", "zetak")

CHI_M4 = [0, 1, 0, -1]
CHI_M20 = [CHI_M4[n % 4] * CHI5[n % 5] for n in range(20)]   # (−20/n) = (−4/n)(n/5), multiplicative


def q(p, r=1):
    return arb(fmpq(p, r))


def is_zero_exact(x):
    return x.is_exact() and x.is_zero()


def tstar():
    """t* = L(¾, χ₅)/(L(¾, χ₅) − G(¾)), G(s) = (1 + √5·5^{−s}) ζ(s); Hurwitz-zeta form of L."""
    s = q(3, 4)
    L5 = sum((CHI5[k] * arb.zeta(s, q(k, 5)) for k in range(1, 5)), arb(0)) * arb(5) ** (-s)
    G = (1 + arb(5).sqrt() * arb(5) ** (-s)) * arb.zeta(s)
    return L5 / (L5 - G)


def kappa_dh():
    r5 = arb(5).sqrt()
    return ((10 - 2 * r5).sqrt() - 2) / (r5 - 1)


def z1_r(nmax):
    r = [0] * (nmax + 1)
    lim = math.isqrt(nmax) + 1
    for x in range(-lim, lim + 1):
        for y in range(-lim, lim + 1):
            v = x * x + 5 * y * y
            if 0 < v <= nmax:
                r[v] += 1
    return r


def coefficients(name, nmax):
    """a_0 … a_nmax as balls (a_0 = 0)."""
    if name == "zeta":
        return [arb(0)] + [arb(1)] * nmax
    if name == "ftstar":
        t = tstar()
        r5 = arb(5).sqrt()
        return [arb(0)] + [(1 - t) * CHI5[n % 5] + t * (1 + (r5 if n % 5 == 0 else 0)) for n in range(1, nmax + 1)]
    if name == "dh":
        k = kappa_dh()
        return [arb(0)] + [DH_CHI[n % 5][0] + k * DH_CHI[n % 5][1] for n in range(1, nmax + 1)]
    if name == "z1":
        r = z1_r(nmax)
        return [arb(0)] + [q(r[n], 2) for n in range(1, nmax + 1)]
    if name == "zetak":
        return [arb(0)] + [arb(sum(CHI_M20[d % 20] for d in range(1, n + 1) if n % d == 0)) for n in range(1, nmax + 1)]
    raise ValueError(name)


def von_mangoldt(n):
    for p in range(2, n + 1):
        if n % p == 0:
            m = n
            while m % p == 0:
                m //= p
            return arb(p).log() if m == 1 else arb(0)
    return arb(0)


def log_derivative(name, nmax):
    """c_0 … c_nmax with −F'/F = Σ c_n n^{−s}: a_n log n = Σ_{d | n} c_d a_{n/d} (a_1 = 1, c_1 = 0)."""
    if name == "zeta":
        return [arb(0), arb(0)] + [von_mangoldt(n) for n in range(2, nmax + 1)]
    if name == "zetak":
        # ζ_K = ζ · L(χ₋₂₀): −ζ_K′/ζ_K = Σ Λ(n)(1 + χ₋₂₀(n)) n^{−s}, exactly (zero off prime powers)
        return [arb(0), arb(0)] + [von_mangoldt(n) * (1 + CHI_M20[n % 20]) for n in range(2, nmax + 1)]
    a = coefficients(name, nmax)
    c = [arb(0)] * (nmax + 1)
    for n in range(2, nmax + 1):
        v = a[n] * arb(n).log()
        for d in range(2, n):
            if n % d == 0 and not is_zero_exact(c[d]):
                v -= c[d] * a[n // d]
        c[n] = v
    return c


class Control:
    """One control function: archimedean data, pole flag and prime comb for a window of full width ℓ."""

    def __init__(self, name):
        if name not in FUNCTIONS:
            raise ValueError(name)
        self.name = name
        self.pole = name != "dh"
        self.kappas = [fmpq(1, 4)] if name in ("zeta", "ftstar") else ([fmpq(3, 4)] if name == "dh" else [fmpq(1, 4), fmpq(3, 4)])  # z1, zetak: degree 2
        logpi = arb.pi().log()
        self.K = {"zeta": -logpi, "ftstar": arb(5).log() - logpi, "dh": arb(5).log() - logpi,
                  "z1": arb(20).log() - 2 * logpi, "zetak": arb(20).log() - 2 * logpi}[name]

    def comb(self, ell, x_int=None):
        """[(n, 2c_n/√n, log n)] for 2 ≤ n < e^ℓ (certified), dropping exact zeros. If e^ℓ = x_int is an
        integer, n = x_int is dropped: log n = ℓ and g(ℓ) = 0 for every f ∈ L²[−ℓ/2, ℓ/2]."""
        nmax = int(math.floor(float(ell.exp().upper()))) + 1
        c = log_derivative(self.name, max(nmax, 2))
        out = []
        for n in range(2, nmax + 1):
            if x_int is not None and n >= x_int:
                continue
            ln = arb(n).log()
            if ln >= ell:
                continue
            if not ln < ell:
                raise ValueError(f"cannot decide log {n} < ell")
            if is_zero_exact(c[n]):
                continue
            out.append((n, 2 * c[n] / arb(n).sqrt(), ln))
        return out

    def psi_arch(self, t):
        """Σ_κ Re ψ(κ + it/2) + K for a real ball t."""
        v = self.K
        for k in self.kappas:
            v += acb(arb(k), t / 2).digamma().real
        return v

    def Psi(self, t, comb):
        v = self.psi_arch(t)
        for _, cn, ln in comb:
            v -= cn * (t * ln).cos()
        return v

    def envelope_beta(self, T, A):
        """β* = Σ_κ [log(T/2) − 1/T] + K − A; Ψ ≥ β* on [T, ∞) (Binet; needs T ≥ 24κ²/(12κ − 1))."""
        for k in self.kappas:
            kk = arb(k)
            tk = 24 * kk * kk / (12 * kk - 1)
            assert T >= tk, f"T = {T} below the envelope threshold {tk}"
        T = arb(T)
        return len(self.kappas) * ((T / 2).log() - 1 / T) + self.K - A
