#!/usr/bin/env python3
"""Zeros-side Weil form for a primitive degree-2 L-function (an elliptic curve E/ℚ), and P-DL4 (docs/DECAY_LAW.md).

General archimedean block. For one factor Γ_R(s + μ), μ ≥ 0 real, put σ = ¼ + μ/2 and, as in connes_letter_mp,
G(w) = ∫ F(v + w) F(v) dv (support [−L, L]), h = |F̂|². With ψ(z) = −γ + ∫_0^∞ (e^{−t} − e^{−zt})/(1 − e^{−t}) dt
(t = 2w), Weil's archimedean term is

  (1/2π)∫ h(t) [Re ψ(σ + it/2) − log π] dt
      = [ψ(σ) − log π] G(0) + ∫_0^L K_σ(w) (G(0) − G(w)) dw + G(0) ∫_L^∞ K_σ,   K_σ(w) = 2 Σ_m e^{−a_m w},

a_m = 2m + ½ + μ: the μ = 0 kernel 2K(w) of connes_letter_mp with every exponent shifted by μ. Keeping our
subtraction 2L e^{−w/2} inside J_k (J_k = ∫_0^L [2(L − w)cos(ωw) − 2L e^{−w/2}] K_μ(w) dw, K_μ = K_σ/2), the
ω-independent part collects into

  −diag_const(μ) = ψ(σ) − log π + ∫_0^∞ K_σ (1 − e^{−w/2}) + ∫_L^∞ K_σ e^{−w/2}
                 = ψ((1 + μ)/2) − log π + Σ_m 2 e^{−(2m+1+μ)L}/(2m + 1 + μ),

using ∫_0^∞ K_σ(1 − e^{−w/2}) dw = ψ(σ + ¼) − ψ(σ). At μ = 0 this is −(log 4π + γ + log tanh(L/2)), our diag_const.
The closed forms of I_k, J_k carry over with z_μ = (½ + μ)/2 + iω/2:
  I_k = ½ Im ψ(z_μ) − Σ_m ω e^{−aL}/(a² + ω²),
  J_k = L(ψ((1 + μ)/2) − Re ψ(z_μ)) − ½ Re ψ'(z_μ) + Σ_m [exact − E-free] (decays like e^{−2mL}).

L(E, s) (analytic normalisation, L(E, s) = Σ a_n n^{−1/2} n^{−s}): Λ(s) = N^{s/2} Γ_C(s + ½) L(E, s) = εΛ(1 − s),
Γ_C(s + ½) = Γ_R(s + ½) Γ_R(s + 3/2). Zeros-side form = block(½) + block(3/2) + log N · I − Σ_{p^k ≤ x} prime terms,
with prime terms (log p^k, (α_p^k + β_p^k) log p / p^{k/2}) in our convention; no pole.

Commands:
  validate   hard validations 1, 2 (block vs build_form at μ = 0 and μ = 1 + log 5) and assembly cross-checks
  curves     a_p by point counting, curve verification (discriminant, conductor, torsion, eta product for 11a1)
  afe        smoothed approximate functional equation checks (two splitting parameters; Λ(s) = εΛ(1 − s))
  zeros      zeros of Λ(½ + it) from sign changes of the real Z(t); count check by the argument principle
  check      hard validation 3: form vs 2 Σ_{γ>0} F(γ)² at x = 13, both sectors
  scan       P-DL4: λ_min at x/N ∈ grid, both sectors, N_basis = 5k*, 9k*, precision-stable
  fit        P-DL4: fit ln λ = −α v + γ ln v + β for the three registered variables, survive/kill per curve and sector
"""

import argparse
import json
import os
import random
import sys
import time

import mpmath as mp

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import connes_letter_mp as cl  # noqa: E402
import decay_law_mp as dl  # noqa: E402
from progress import Job  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "..", "data", "connes", "gl2")


# ---------------------------------------------------------------------------------------------------------------
# Archimedean block for Γ_R(s + μ) and assembly on our bases
# ---------------------------------------------------------------------------------------------------------------

def arch_vectors(L, n, mu):
    """(I_k, J_k for k = 0..n, diag_const) for one factor Γ_R(s + μ) on the window of length L."""
    mu = mp.mpf(mu)
    half = mp.mpf(1) / 2
    tail_terms = int(mp.mp.dps * mp.log(10) / (2 * L)) + 10
    sig = (half + mu) / 2           # Re z_μ
    b = (1 + mu) / 2                # ψ(b) from Σ 1/(a + ½), a + ½ = 2(m + b)
    om = lambda k: 2 * mp.pi * k / L

    def j_term(a, w):
        E = mp.exp(-a * L)
        d = a * a + w * w
        c = L * a * (1 - E) / d + (((1 - E) + a * L * E) * d - 2 * a * a * (1 - E)) / (d * d)
        return 2 * c - 2 * L * (1 - mp.exp(-(a + half) * L)) / (a + half)

    def t0(a, w):
        d = a * a + w * w
        return 2 * L * a / d + 2 * (w * w - a * a) / (d * d) - 2 * L / (a + half)

    I, J = [], []
    for k in range(n + 1):
        w = om(k)
        z = mp.mpc(sig, w / 2)
        if k == 0:
            I.append(mp.mpf(0))
        else:
            tail = mp.fsum(w * mp.exp(-(2 * m + half + mu) * L) / ((2 * m + half + mu) ** 2 + w * w) for m in range(tail_terms))
            I.append(mp.im(mp.digamma(z)) / 2 - tail)
        main = L * (mp.digamma(b) - mp.re(mp.digamma(z))) - mp.re(mp.psi(1, z)) / 2
        J.append(main + mp.fsum(j_term(2 * m + half + mu, w) - t0(2 * m + half + mu, w) for m in range(tail_terms)))
    diag_const = mp.log(mp.pi) - mp.digamma(b) - mp.fsum(2 * mp.exp(-(2 * m + 1 + mu) * L) / (2 * m + 1 + mu) for m in range(tail_terms))
    return I, J, diag_const


def assemble(L, n, parity, I, J, diag_const, terms, extra_diag=0):
    """Same matrix as connes_letter_mp.build_form, from the arch vectors (I, J, diag_const), the prime list
    terms = [(log n, coefficient)] and a constant added to the diagonal (the conductor term log N).
    In the e_k basis: q(k, k) = (−(diag_const L + J_k) − D_k)/L + extra_diag and, for j ≠ k,
    q(j, k) = (−1)^{j−k} (A_k − A_j)/(2π(j − k)) with A_k = −2(I_k + P_k) (odd in k)."""
    om = lambda k: 2 * mp.pi * k / L
    P = [mp.fsum(a * mp.sin(om(k) * ln) for ln, a in terms) for k in range(n + 1)]
    D = [mp.fsum(a * 2 * (L - ln) * mp.cos(om(k) * ln) for ln, a in terms) for k in range(n + 1)]
    A = [-2 * (I[k] + P[k]) for k in range(n + 1)]
    d = [(-(diag_const * L + J[k]) - D[k]) / L + extra_diag for k in range(n + 1)]
    twopi = 2 * mp.pi

    def q(j, k):
        if j == k:
            return d[abs(k)]
        Ak = A[k] if k >= 0 else -A[-k]
        Aj = A[j] if j >= 0 else -A[-j]
        sign = 1 if (j - k) % 2 == 0 else -1
        return sign * (Ak - Aj) / (twopi * (j - k))

    if parity == "odd":
        M = mp.matrix(n, n)
        for a in range(1, n + 1):
            for b in range(a, n + 1):
                v = (q(a, b) - q(a, -b) - q(-a, b) + q(-a, -b)) / 2
                M[a - 1, b - 1] = v
                M[b - 1, a - 1] = v
        return M
    s2 = mp.sqrt(2)
    M = mp.matrix(n + 1, n + 1)
    for a in range(n + 1):
        for b in range(a, n + 1):
            if a == 0 and b == 0:
                v = q(0, 0)
            elif a == 0:
                v = (q(0, b) + q(0, -b)) / s2
            else:
                v = (q(a, b) + q(a, -b) + q(-a, b) + q(-a, -b)) / 2
            M[a, b] = v
            M[b, a] = v
    return M


def form(x, n, parity, mus, log_cond=0, terms=()):
    """Zeros-side-style form: Σ_μ block(μ) + log_cond · I − primes(terms)."""
    x = mp.mpf(x)
    L = mp.log(x)
    I = [mp.mpf(0)] * (n + 1)
    J = [mp.mpf(0)] * (n + 1)
    dc = mp.mpf(0)
    for mu in mus:
        Ii, Ji, di = arch_vectors(L, n, mu)
        I = [u + v for u, v in zip(I, Ii)]
        J = [u + v for u, v in zip(J, Ji)]
        dc += di
    return assemble(L, n, parity, I, J, dc, list(terms), extra_diag=log_cond)


def maxdiff(A, B):
    return max(abs(A[i, j] - B[i, j]) for i in range(A.rows) for j in range(A.cols))


def maxabs(A):
    return max(abs(A[i, j]) for i in range(A.rows) for j in range(A.cols))


def cmd_validate(args):
    """Hard validations 1 and 2, plus cross-checks of the prime/conductor path against decay_law_mp.zeros_side."""
    mp.mp.dps = args.dps
    job = Job("gl2 validate", args=vars(args))
    out = []
    ok = True
    tol = mp.mpf(10) ** (-(args.dps - 10))
    for xs in args.x.split(","):
        x = mp.mpf(xs)
        for parity in ("even", "odd"):
            n = args.n
            ref0 = cl.build_form(x, n, "zeta", terms=[], parity=parity)[1]
            M0 = form(x, n, parity, [0])
            d0 = maxdiff(M0, ref0)
            ref1 = cl.build_form(x, n, "dh", terms=[], parity=parity)[1]
            M1 = form(x, n, parity, [1], log_cond=mp.log(5))
            d1 = maxdiff(M1, ref1)
            # prime path: full ζ (Γ_R(s), Λ(n)/√n), even χ_5 and odd χ_{−4} against zeros_side (pole removed for ζ)
            pp = [(mp.log(nn), lp / mp.sqrt(nn)) for nn, lp in cl.prime_powers(int(x))]
            dz = maxdiff(form(x, n, parity, [0], terms=pp), cl.build_form(x, n, "zeta", parity=parity)[1])
            chi = lambda D: [(mp.log(nn), dl.kronecker(D, nn) * lp / mp.sqrt(nn)) for nn, lp in cl.prime_powers(int(x)) if dl.kronecker(D, nn)]
            d5 = maxdiff(form(x, n, parity, [0], log_cond=mp.log(5), terms=chi(5)), dl.zeros_side(5, x, n, parity))
            d4 = maxdiff(form(x, n, parity, [1], log_cond=mp.log(4), terms=chi(-4)), dl.zeros_side(-4, x, n, parity))
            row = {"x": xs, "parity": parity, "n": n, "dps": args.dps, "max_entry": mp.nstr(maxabs(ref0), 5),
                   "V1_mu0_vs_zeta": mp.nstr(d0, 3), "V2_mu1_log5_vs_dh": mp.nstr(d1, 3),
                   "zeta_with_primes": mp.nstr(dz, 3), "chi5_vs_zeros_side": mp.nstr(d5, 3), "chi-4_vs_zeros_side": mp.nstr(d4, 3)}
            passed = all(v < tol for v in (d0, d1, dz, d5, d4))
            ok &= passed
            row["pass"] = passed
            out.append(row)
            job.result(f"x={xs} {parity} n={n} dps={args.dps}: max|block(0) − build_form(zeta)| = {mp.nstr(d0, 3)};  "
                       f"max|block(1) + log5 − build_form(dh)| = {mp.nstr(d1, 3)};  ζ with primes {mp.nstr(dz, 3)};  "
                       f"χ5 {mp.nstr(d5, 3)};  χ−4 {mp.nstr(d4, 3)}  (max entry {mp.nstr(maxabs(ref0), 4)}) → {'PASS' if passed else 'FAIL'}")
    job.result(f"hard validations 1 and 2: {'PASS' if ok else 'FAIL'} (tolerance {mp.nstr(tol, 2)})")
    if args.json:
        with open(args.json, "w") as fh:
            json.dump({"rows": out, "pass": ok, "tolerance": mp.nstr(tol, 3)}, fh, indent=1)
    job.done()
    if not ok:
        sys.exit(1)


# ---------------------------------------------------------------------------------------------------------------
# Elliptic curves: a_p by point counting, a_n, Euler-product prime terms
# ---------------------------------------------------------------------------------------------------------------

# Cremona labels, minimal Weierstrass models [a1, a2, a3, a4, a6]. 11a1, 37b1: P-DL4. 15a1, 19a1: P-DL6.
CURVES = {
    "11a1": {"ainv": [0, -1, 1, -10, -20], "N": 11},
    "37b1": {"ainv": [0, 1, 1, -23, -50], "N": 37},
    "15a1": {"ainv": [1, 1, 1, -10, -10], "N": 15},
    "19a1": {"ainv": [0, 1, 1, -9, -15], "N": 19},
    "37a1": {"ainv": [0, 0, 1, -1, 0], "N": 37, "eps": -1},          # P-DL9: rank 1, ε = −1 (verified by the AFE)
    # P-DL12 (models from memory; conductor, ε and rank verified by objects6 / afe-both). eps = None until the AFE decides.
    "14a1": {"ainv": [1, 0, 1, 4, -6], "N": 14, "eps": 1},
    "17a1": {"ainv": [1, -1, 1, -1, -14], "N": 17, "eps": 1},
    "43a1": {"ainv": [0, 1, 1, 0, 0], "N": 43, "eps": -1},
    "53a1": {"ainv": [1, -1, 1, 0, 0], "N": 53, "eps": -1},
    # P-R389: the smallest-conductor rank-2 curve; ε = +1 with a double central zero (m0 = 2), verified by afe-both.
    "389a1": {"ainv": [0, 1, 1, -2, 0], "N": 389, "eps": 1, "m0": 2},
}
KNOWN_AP = {"37a1": {2: -2, 3: -3, 5: -2}}                           # values supplied with the P-DL9 task, checked

# Ramanujan's Δ (P-DL6): weight 12, level 1, L(Δ, s) = Σ τ(n) n^{−11/2} n^{−s}, Γ_C(s + 11/2), ε = +1.
DELTA = "delta"


# P-DL7: the level-1 normalised Hecke eigenforms f_k = Δ·E_{k−12} (S_k(1) is one-dimensional), root number i^k.
EIGEN = {"f16": 16, "f18": 18, "f20": 20, "f22": 22, "f26": 26}


# P-DL11: newforms that are single eta products, in one-dimensional spaces S_k(Γ₀(N)).
# "eps" is the root number found by the AFE (both signs tested; see objects11 / afe_checks11.json).
ETA_NEWFORMS = {
    "g4": {"k": 4, "N": 5, "factors": ((1, 4), (5, 4)), "eps": 1},      # η(τ)⁴η(5τ)⁴
    "g6": {"k": 6, "N": 3, "factors": ((1, 6), (3, 6)), "eps": 1},      # η(τ)⁶η(3τ)⁶
    "g8": {"k": 8, "N": 2, "factors": ((1, 8), (2, 8)), "eps": 1},      # η(τ)⁸η(2τ)⁸
}


def N_of(label):
    if label in ETA_NEWFORMS:
        return ETA_NEWFORMS[label]["N"]
    return 1 if (label == DELTA or label in EIGEN) else CURVES[label]["N"]


def weight(label):
    if label in ETA_NEWFORMS:
        return ETA_NEWFORMS[label]["k"]
    return 12 if label == DELTA else EIGEN.get(label, 2)


def m0_of(label):
    """Order of the central zero. Defaults to the value forced by ε for rank ≤ 1 (0 if ε = +1, 1 if ε = −1); a curve
    may override it ("m0" in CURVES), as 389a1 (ε = +1, m0 = 2) does."""
    if label in CURVES and "m0" in CURVES[label]:
        return CURVES[label]["m0"]
    return 0 if eps_of(label) == 1 else 1


def eps_of(label):
    """Root number: +1 for the curves (verified) and Δ; i^k for level-1 forms of weight k."""
    if label in EIGEN:
        return 1 if EIGEN[label] % 4 == 0 else -1
    if label in CURVES:
        e = CURVES[label].get("eps", 1)
        if e is None:
            raise ValueError(f"root number of {label} not yet determined (run afe-both)")
        return e
    if label in ETA_NEWFORMS:
        e = ETA_NEWFORMS[label]["eps"]
        if e is None:
            raise ValueError(f"root number of {label} not yet determined (run afe with --eps-test)")
        return e
    return 1


def _sigma_table(e, nmax):
    s = [0] * (nmax + 1)
    for d in range(1, nmax + 1):
        de = d ** e
        for m in range(d, nmax + 1, d):
            s[m] += de
    return s


_ETA_CACHE = {}


def eta_product_fast(nmax, factors):
    """a_0..a_nmax of q Π_n Π_{(m, e)} (1 − q^{mn})^e (Σ m e = 24), from Euler's pentagonal series
    Π(1 − q^n) = Σ_j (−1)^j q^{j(3j−1)/2} (j ∈ ℤ) and FLINT polynomial powers; the naive product is a cross-check."""
    key = (nmax, factors)
    if key in _ETA_CACHE:
        return _ETA_CACHE[key]
    from flint import fmpz_poly
    L = nmax                                    # need powers q^0..q^{nmax−1} of the product
    base = [0] * L
    j = 0
    while True:
        done = True
        for jj in (j, -j) if j else (0,):
            e = jj * (3 * jj - 1) // 2
            if e < L:
                base[e] += (-1) ** abs(jj)
                done = False
        if done and j > 0:
            break
        j += 1
    trunc = lambda Q: fmpz_poly([int(c) for c in Q.coeffs()[:L]])
    P = fmpz_poly([1])
    for m, e in factors:
        sub = [0] * L
        for i, c in enumerate(base):
            if c and i * m < L:
                sub[i * m] = c
        Pm = fmpz_poly(sub)
        for _ in range(e):
            P = trunc(P * Pm)
    c = [int(v) for v in P.coeffs()] + [0] * (L + 1)
    out = [0] + c[:nmax]                        # a_n = coefficient of q^{n−1} in the product
    _ETA_CACHE[key] = out
    return out


_EIGEN_CACHE = {}


def eigen_table(k, nmax):
    """a_n, n = 0..nmax, of f_k = Δ·E_4^a E_6^b with 4a + 6b = k − 12 (E_4 = 1 + 240Σσ_3 q^n, E_6 = 1 − 504Σσ_5 q^n)."""
    key = (k, nmax)
    if key in _EIGEN_CACHE:
        return _EIGEN_CACHE[key]
    from flint import fmpz_poly
    a, b = {4: (1, 0), 6: (0, 1), 8: (2, 0), 10: (1, 1), 14: (2, 1)}[k - 12]
    s3, s5 = _sigma_table(3, nmax), _sigma_table(5, nmax)
    E4 = fmpz_poly([1] + [240 * s3[n] for n in range(1, nmax + 1)])
    E6 = fmpz_poly([1] + [-504 * s5[n] for n in range(1, nmax + 1)])
    tau = tau_table(nmax + 1)
    P = fmpz_poly([0] + tau[1:nmax + 1])
    trunc = lambda Q: fmpz_poly([int(c) for c in Q.coeffs()[:nmax + 1]])
    for _ in range(a):
        P = trunc(P * E4)
    for _ in range(b):
        P = trunc(P * E6)
    c = [int(v) for v in P.coeffs()] + [0] * (nmax + 1)
    _EIGEN_CACHE[key] = c[:nmax + 1]
    return _EIGEN_CACHE[key]


def mus_of(label):
    """Γ_C(s + (k−1)/2) = Γ_R(s + (k−1)/2) Γ_R(s + (k+1)/2): μ = ½, 3/2 for curves, 11/2, 13/2 for Δ."""
    k = weight(label)
    return [mp.mpf(k - 1) / 2, mp.mpf(k + 1) / 2]


_TAU_CACHE = {}


def tau_table(nmax):
    """τ(n), n = 0..nmax, from Δ = q Π (1 − q^n)^24, computed as q·(Π(1 − q^n)³)^8 with Jacobi's identity
    Π(1 − q^n)³ = Σ_m (−1)^m (2m + 1) q^{m(m+1)/2} (FLINT integer polynomials)."""
    if nmax in _TAU_CACHE:
        return _TAU_CACHE[nmax]
    from flint import fmpz_poly
    coeffs = [0] * nmax
    m = 0
    while m * (m + 1) // 2 < nmax:
        coeffs[m * (m + 1) // 2] += (-1) ** m * (2 * m + 1)
        m += 1
    A = fmpz_poly(coeffs)
    trunc = lambda P: fmpz_poly([int(c) for c in P.coeffs()[:nmax]])
    A2 = trunc(A * A)
    A4 = trunc(A2 * A2)
    A8 = trunc(A4 * A4)
    c = [int(v) for v in A8.coeffs()] + [0] * nmax
    tau = [0] + c[:nmax]                       # τ(n) = coefficient of q^{n−1} in Π(1 − q^n)^24
    _TAU_CACHE[nmax] = tau
    return tau


def tau_naive(nmax):
    """τ(n) by multiplying out Π_{n<nmax} (1 − q^n)^24 term by term (independent of Jacobi's identity; slow)."""
    c = [0] * nmax
    c[0] = 1
    for n in range(1, nmax):
        for _ in range(24):
            for i in range(nmax - 1, n - 1, -1):
                c[i] -= c[i - n]
    return [0] + c[:nmax]


def primes_upto(n):
    sieve = bytearray([1]) * (n + 1)
    sieve[0:2] = b"\x00\x00"
    for p in range(2, int(n ** 0.5) + 1):
        if sieve[p]:
            sieve[p * p::p] = bytearray(len(sieve[p * p::p]))
    return [p for p in range(n + 1) if sieve[p]]


def invariants(ainv):
    a1, a2, a3, a4, a6 = ainv
    b2 = a1 * a1 + 4 * a2
    b4 = 2 * a4 + a1 * a3
    b6 = a3 * a3 + 4 * a6
    b8 = a1 * a1 * a6 + 4 * a2 * a6 - a1 * a3 * a4 + a2 * a3 * a3 - a4 * a4
    c4 = b2 * b2 - 24 * b4
    c6 = -b2 ** 3 + 36 * b2 * b4 - 216 * b6
    disc = -b2 * b2 * b8 - 8 * b4 ** 3 - 27 * b6 * b6 + 9 * b2 * b4 * b6
    return {"b2": b2, "b4": b4, "b6": b6, "b8": b8, "c4": c4, "c6": c6, "disc": disc}


def count_points(ainv, p):
    """#Ẽ(F_p): all projective points of the reduced Weierstrass cubic (∞ and any singular point included).
    Then a_p = p + 1 − #Ẽ(F_p) at good p and at multiplicative/additive p alike (model minimal at p)."""
    a1, a2, a3, a4, a6 = ainv
    cnt = 1
    if p == 2:
        for x in range(2):
            for y in range(2):
                if (y * y + a1 * x * y + a3 * y - (x ** 3 + a2 * x * x + a4 * x + a6)) % 2 == 0:
                    cnt += 1
        return cnt
    for x in range(p):
        f = x ** 3 + a2 * x * x + a4 * x + a6
        b = a1 * x + a3
        disc = (b * b + 4 * f) % p
        if disc == 0:
            cnt += 1
        elif pow(disc, (p - 1) // 2, p) == 1:
            cnt += 2
    return cnt


_AP_CACHE = {}


def ap_table(label, pmax):
    """{p: a_p}: point counting for curves; τ(p) for Δ."""
    key = (label, pmax)
    if key not in _AP_CACHE:
        if label == DELTA:
            tau = tau_table(max(pmax, 2))
            _AP_CACHE[key] = {p: tau[p] for p in primes_upto(pmax)}
        elif label in EIGEN:
            co = eigen_table(EIGEN[label], max(pmax, 2))
            _AP_CACHE[key] = {p: co[p] for p in primes_upto(pmax)}
        elif label in ETA_NEWFORMS:
            co = eta_product_fast(max(pmax, 2) + 1, ETA_NEWFORMS[label]["factors"])
            _AP_CACHE[key] = {p: co[p] for p in primes_upto(pmax)}
        else:
            ainv = CURVES[label]["ainv"]
            _AP_CACHE[key] = {p: p + 1 - count_points(ainv, p) for p in primes_upto(pmax)}
    return _AP_CACHE[key]


def an_table(label, nmax):
    """Arithmetic a_n, n = 0..nmax (a_0 = 0). Curves: a_{p^k} = a_p a_{p^{k−1}} − p a_{p^{k−2}} at good p, a_p^k at
    p | N. Δ: τ(n) straight from the q-expansion (hecke_table rebuilds it from τ(p) as a check)."""
    if label == DELTA:
        return tau_table(nmax)
    if label in EIGEN:
        return eigen_table(EIGEN[label], nmax)
    if label in ETA_NEWFORMS:
        return eta_product_fast(nmax + 1, ETA_NEWFORMS[label]["factors"])[:nmax + 1]
    return hecke_table(label, nmax)


def hecke_table(label, nmax):
    """a_n from a_p by multiplicativity and the Hecke recursion a_{p^k} = a_p a_{p^{k−1}} − χ(p) p^{w−1} a_{p^{k−2}}
    (w the weight; χ(p) = 0 at p | N, where a_{p^k} = a_p^k)."""
    N = N_of(label)
    pw = weight(label) - 1
    ap = ap_table(label, nmax)
    a = [0] * (nmax + 1)
    a[1] = 1
    spf = list(range(nmax + 1))
    for p in primes_upto(int(nmax ** 0.5) + 1):
        for m in range(p * p, nmax + 1, p):
            if spf[m] == m:
                spf[m] = p
    for n in range(2, nmax + 1):
        p = spf[n]
        m, k = n, 0
        while m % p == 0:
            m //= p
            k += 1
        pk = n // m
        if m > 1:
            a[n] = a[pk] * a[m]
            continue
        # n = p^k
        if k == 1:
            a[n] = ap[p]
        elif N % p == 0:
            a[n] = ap[p] * a[n // p]
        else:
            a[n] = ap[p] * a[n // p] - p ** pw * a[n // (p * p)]
    return a


def euler_terms(label, x):
    """Prime terms [(log p^k, (α_p^k + β_p^k) log p / p^{k/2})] for p^k ≤ x (our convention; analytic normalisation:
    α_p + β_p = a_p/p^{(w−1)/2}, α_p β_p = 1 at good p; α_p = a_p/√p, β_p = 0 at p | N; w = weight)."""
    N = N_of(label)
    x = int(mp.floor(x))
    ap = ap_table(label, max(x, 2))
    half = mp.mpf(weight(label) - 1) / 2
    out = []
    for p in primes_upto(x):
        s1 = mp.mpf(ap[p]) / (mp.sqrt(p) if weight(label) == 2 else mp.mpf(p) ** half)
        prev2, prev1 = mp.mpf(2), s1          # s_0, s_1 at good p
        pk, k = p, 1
        while pk <= x:
            if N % p == 0:
                sk = s1 ** k
            elif k == 1:
                sk = s1
            else:
                sk = s1 * prev1 - prev2
                prev2, prev1 = prev1, sk
            out.append((mp.log(pk), sk * mp.log(p) / mp.sqrt(pk)))
            pk *= p
            k += 1
    return out


def log_derivative_coeffs(label, nmax):
    """Λ_E(n) of −L'/L(E, s) = Σ Λ_E(n) n^{−s} (analytic normalisation) from the Dirichlet coefficients b_n = a_n/√n
    by c_n = b_n log n − Σ_{d | n, 1 < d < n} c_d b_{n/d}: an independent route to the prime-power coefficients."""
    a = an_table(label, nmax)
    half = mp.mpf(weight(label) - 1) / 2
    b = [mp.mpf(0)] + [mp.mpf(a[n]) / (mp.sqrt(n) if weight(label) == 2 else mp.mpf(n) ** half) for n in range(1, nmax + 1)]
    c = [mp.mpf(0)] * (nmax + 1)
    for n in range(2, nmax + 1):
        c[n] = b[n] * mp.log(n)
    for d in range(2, nmax + 1):
        if c[d] == 0:
            continue
        for m in range(2 * d, nmax + 1, d):
            c[m] -= c[d] * b[m // d]
    return c


def zeros_side_E(label, x, n, parity):
    """Zeros-side form: block(μ1) + block(μ2) + log N · I − Σ_{p^k ≤ x} prime terms; (μ1, μ2) = (½, 3/2) for curves,
    (11/2, 13/2) for Δ (N = 1)."""
    return form(x, n, parity, mus_of(label), log_cond=mp.log(N_of(label)), terms=euler_terms(label, x))


def eta_product(nmax, factors):
    """Coefficients a_1..a_nmax of q Π_{n≥1} Π_{(m, e) in factors} (1 − q^{mn})^e (Σ m e = 24, so the leading power is q)."""
    c = [0] * nmax                     # coefficients of Π (series in q, index = power), length nmax (powers 0..nmax−1)
    c[0] = 1
    for n in range(1, nmax):
        for mult, reps in factors:
            e = mult * n
            if e >= nmax:
                continue
            for _ in range(reps):
                for i in range(nmax - 1, e - 1, -1):
                    c[i] -= c[i - e]
    return [0] + c[:nmax]              # a_n = coefficient of q^n = c[n − 1]


def eta_product_11(nmax):
    """q Π (1 − q^n)² (1 − q^{11n})², the weight-2 newform of level 11."""
    return eta_product(nmax, ((1, 2), (11, 2)))


# Weight-2 newforms that are eta products (independent of point counting): level 11 and level 15.
ETA_FORMS = {"11a1": ((1, 2), (11, 2)), "15a1": ((1, 1), (3, 1), (5, 1), (15, 1)), "14a1": ((1, 1), (2, 1), (7, 1), (14, 1))}


def cmd_curves(args):
    """Curve verification: discriminant and conductor, a_p, torsion congruence, eta product (11a1), log-derivative."""
    mp.mp.dps = 30
    job = Job("gl2 curves", args=vars(args))
    report = {}
    for label, cv in ((lb, CURVES[lb]) for lb in ("11a1", "37b1")):   # P-DL4 curves (prime N)
        N, ainv = cv["N"], cv["ainv"]
        inv = invariants(ainv)
        D = inv["disc"]
        m, k = abs(D), 0
        while m % N == 0:
            m //= N
            k += 1
        cond_ok = (m == 1 and inv["c4"] % N != 0 and k < 12)   # Δ = ±N^k, c4 prime to N: multiplicative at N only
        ap = ap_table(label, args.pmax)
        tors = 5 if label == "11a1" else 3
        tors_ok = all((p + 1 - ap[p]) % tors == 0 for p in ap if N % p and p != tors)
        hasse_ok = all(ap[p] ** 2 <= 4 * p for p in ap if N % p)
        job.result(f"{label} {ainv}: Δ = {D} = {'−' if D < 0 else ''}{N}^{k}, c4 = {inv['c4']} (N ∤ c4: {inv['c4'] % N != 0}), "
                   f"j = {mp.nstr(mp.mpf(inv['c4']) ** 3 / D, 12)} → conductor {N}: {cond_ok};  a_{N} = {ap[N]} "
                   f"({'split' if ap[N] == 1 else 'non-split'} multiplicative);  #Ẽ(F_p) ≡ 0 mod {tors} for good p ≤ {args.pmax}: {tors_ok};  Hasse: {hasse_ok}")
        job.result(f"{label} a_p, p ≤ 50: " + ", ".join(f"a_{p}={ap[p]}" for p in ap if p <= 50))
        a = an_table(label, args.nmax)
        rep = {"ainv": ainv, "N": N, "disc": D, "disc_exponent": k, "c4": inv["c4"], "conductor_ok": cond_ok,
               "torsion_congruence_mod": tors, "torsion_ok": tors_ok, "hasse_ok": hasse_ok,
               "a_p_upto_50": {p: ap[p] for p in ap if p <= 50}}
        if label == "11a1":
            known = {2: -2, 3: -1, 5: 1, 7: -2}
            known_ok = all(ap[p] == v for p, v in known.items())
            eta = eta_product_11(args.nmax + 1)
            bad = [n for n in range(1, args.nmax + 1) if eta[n] != a[n]]
            job.result(f"11a1: a_2, a_3, a_5, a_7 = {[ap[p] for p in known]} vs known {list(known.values())}: {known_ok};  "
                       f"a_n (point counting + Hecke) vs η(τ)²η(11τ)² for n ≤ {args.nmax}: {len(bad)} mismatches")
            rep.update({"known_ap_ok": known_ok, "eta_mismatches": len(bad), "eta_nmax": args.nmax})
        # prime-power coefficients two ways
        c = log_derivative_coeffs(label, args.lognmax)
        et = {int(mp.nint(mp.exp(ln))): v for ln, v in euler_terms(label, args.lognmax)}
        dev_pp = max(abs(c[n] / mp.sqrt(n) - et[n]) for n in et)
        dev_other = max(abs(c[n]) for n in range(2, args.lognmax + 1) if n not in et)
        job.result(f"{label}: Euler terms vs Dirichlet log-derivative for n ≤ {args.lognmax}: max |Δ| on prime powers {mp.nstr(dev_pp, 3)}, "
                   f"max |c_n| off prime powers {mp.nstr(dev_other, 3)}")
        rep.update({"logderiv_dev_primepowers": mp.nstr(dev_pp, 3), "logderiv_max_off_primepowers": mp.nstr(dev_other, 3)})
        report[label] = rep
    if args.json:
        with open(args.json, "w") as fh:
            json.dump(report, fh, indent=1)
    job.done()


def _factor(n):
    out, p = {}, 2
    n = abs(n)
    while p * p <= n:
        while n % p == 0:
            out[p] = out.get(p, 0) + 1
            n //= p
        p += 1
    if n > 1:
        out[n] = out.get(n, 0) + 1
    return out


def _logderiv_check(label, nmax):
    c = log_derivative_coeffs(label, nmax)
    et = {int(mp.nint(mp.exp(ln))): v for ln, v in euler_terms(label, nmax)}
    # Euler terms are native entries Λ(n)/√n; the recursion gives Λ(n) itself.
    dev_pp = max(abs(c[n] / mp.sqrt(n) - et[n]) for n in et)
    dev_other = max(abs(c[n]) for n in range(2, nmax + 1) if n not in et)
    return dev_pp, dev_other


def cmd_objects6(args):
    """P-DL6 object verification. Curves: conductor (Tate: Δ = ±Π_{p|N} p^{k_p}, k_p < 12, c4 prime to N ⇒ N
    squarefree multiplicative), a_p by point counting, gcd of #Ẽ(F_p) (torsion), Hasse, eta product (15a1), Euler
    terms vs log-derivative. Δ: τ(n) two ways (Jacobi route vs naive product), known values, τ(n) ≡ σ_11(n) mod 691,
    Hecke multiplicativity from τ(p) alone, Deligne's bound, Euler terms vs log-derivative."""
    import math
    mp.mp.dps = 30
    job = Job("gl2 objects6", args=vars(args))
    report = {}
    for label in args.objects.split(","):
        if label == DELTA:
            tau = tau_table(args.nmax)
            naive = tau_naive(args.naive)
            two_ways = sum(1 for n in range(1, args.naive + 1) if naive[n] != tau[n])
            known = [1, -24, 252, -1472, 4830, -6048, -16744, 84480, -113643, -115920]
            known_ok = tau[1:11] == known
            sig11 = lambda n: sum(d ** 11 for d in range(1, n + 1) if n % d == 0)
            cong_bad = sum(1 for n in range(1, args.cong + 1) if (tau[n] - sig11(n)) % 691)
            hk = hecke_table(DELTA, args.nmax)
            hecke_bad = sum(1 for n in range(1, args.nmax + 1) if hk[n] != tau[n])
            ap = ap_table(DELTA, args.nmax)
            deligne = all(ap[p] ** 2 <= 4 * p ** 11 for p in ap)
            dpp, doth = _logderiv_check(DELTA, args.lognmax)
            job.result(f"Δ: τ(1..10) = {tau[1:11]} vs known: {known_ok};  Jacobi route vs naive Π(1 − q^n)^24, n ≤ {args.naive}: "
                       f"{two_ways} mismatches;  τ(n) ≢ σ_11(n) mod 691 for n ≤ {args.cong}: {cong_bad};  Hecke (from τ(p), p^11) "
                       f"vs q-expansion, n ≤ {args.nmax}: {hecke_bad} mismatches;  Deligne |τ(p)| ≤ 2p^{{11/2}}, p ≤ {args.nmax}: {deligne};  "
                       f"Euler terms vs log-derivative (n ≤ {args.lognmax}): {mp.nstr(dpp, 3)}, off prime powers {mp.nstr(doth, 3)}")
            report[label] = {"tau_1_10": tau[1:11], "known_ok": known_ok, "jacobi_vs_naive_mismatches": two_ways, "naive_nmax": args.naive,
                             "ramanujan_691_failures": cong_bad, "cong_nmax": args.cong, "hecke_mismatches": hecke_bad,
                             "deligne_ok": deligne, "logderiv_dev_primepowers": mp.nstr(dpp, 3), "logderiv_max_off_primepowers": mp.nstr(doth, 3),
                             "weight": 12, "N": 1, "mus": ["11/2", "13/2"]}
            continue
        cv = CURVES[label]
        N, ainv = cv["N"], cv["ainv"]
        inv = invariants(ainv)
        D = inv["disc"]
        fN, fD = _factor(N), _factor(D)
        cond_ok = (all(e == 1 for e in fN.values()) and set(fD) == set(fN) and all(fD[p] < 12 for p in fD)
                   and all(inv["c4"] % p for p in fN))
        ap = ap_table(label, args.pmax)
        good = [p for p in ap if N % p]
        g = g_odd = 0
        for p in good:
            g = math.gcd(g, p + 1 - ap[p])
            if p > 2:
                g_odd = math.gcd(g_odd, p + 1 - ap[p])
        hasse_ok = all(ap[p] ** 2 <= 4 * p for p in good)
        bad = {p: ap[p] for p in fN}
        job.result(f"{label} {ainv}: Δ = {D} = {'−' if D < 0 else ''}{' · '.join(f'{p}^{e}' for p, e in sorted(fD.items()))}, c4 = {inv['c4']} "
                   f"(prime to N: {all(inv['c4'] % p for p in fN)}), j = {mp.nstr(mp.mpf(inv['c4']) ** 3 / D, 12)} → conductor {N}: {cond_ok};  "
                   f"a_p at p | N: {bad} ({', '.join('split' if v == 1 else 'non-split' for v in bad.values())} multiplicative);  "
                   f"gcd of #Ẽ(F_p) over good p ≤ {args.pmax}: {g} (odd p: {g_odd});  Hasse: {hasse_ok}")
        job.result(f"{label} a_p, p ≤ 50: " + ", ".join(f"a_{p}={ap[p]}" for p in ap if p <= 50))
        known_ok = None
        if label in KNOWN_AP:
            known_ok = all(ap[p] == v for p, v in KNOWN_AP[label].items())
            job.result(f"{label}: a_p for p in {list(KNOWN_AP[label])} = {[ap[p] for p in KNOWN_AP[label]]} vs supplied {list(KNOWN_AP[label].values())}: {known_ok}")
        rep = {"known_ap_ok": known_ok, "ainv": ainv, "N": N, "disc": D, "disc_factor": {str(p): e for p, e in fD.items()}, "c4": inv["c4"], "conductor_ok": cond_ok,
               "bad_ap": {str(p): v for p, v in bad.items()}, "gcd_points": g, "gcd_points_odd_p": g_odd, "hasse_ok": hasse_ok,
               "a_p_upto_50": {p: ap[p] for p in ap if p <= 50}}
        if label in ETA_FORMS:
            a = an_table(label, args.nmax)
            eta = eta_product(args.nmax + 1, ETA_FORMS[label])
            nbad = sum(1 for n in range(1, args.nmax + 1) if eta[n] != a[n])
            job.result(f"{label}: a_n (point counting + Hecke) vs eta product {ETA_FORMS[label]}, n ≤ {args.nmax}: {nbad} mismatches")
            rep.update({"eta_factors": ETA_FORMS[label], "eta_mismatches": nbad, "eta_nmax": args.nmax})
        dpp, doth = _logderiv_check(label, args.lognmax)
        job.result(f"{label}: Euler terms vs Dirichlet log-derivative for n ≤ {args.lognmax}: max |Δ| on prime powers {mp.nstr(dpp, 3)}, "
                   f"max |c_n| off prime powers {mp.nstr(doth, 3)}")
        rep.update({"logderiv_dev_primepowers": mp.nstr(dpp, 3), "logderiv_max_off_primepowers": mp.nstr(doth, 3)})
        report[label] = rep
    if args.json:
        with open(args.json, "w") as fh:
            json.dump(report, fh, indent=1)
    job.done()


# ---------------------------------------------------------------------------------------------------------------
# Λ(E, s) by the smoothed approximate functional equation
# ---------------------------------------------------------------------------------------------------------------

def afe_arith(label, w, t=1, eps=1, an=None, extra_digits=0):
    """Λ_a(w) = (√N/2π)^w Γ(w) Σ a_n n^{−w} (arithmetic normalisation of weight k: w = s + (k−1)/2,
    Λ_a(w) = εΛ_a(k − w)), from Λ_a(w) = Σ_n a_n [c_n^{−w} Γ(w, c_n t) + ε c_n^{w−k} Γ(k − w, c_n/t)], c_n = 2πn/√N,
    for any splitting t > 0 (f(1/y) = ε y^k f(y) for f(y) = Σ a_n e^{−2πny/√N}). Uses the caller's mp.dps.
    k = 2 keeps the P-DL4 term count exactly; for k > 2 the count also covers |a_n| ~ n^{(k−1)/2} and Re(k − w)."""
    N = N_of(label)
    k = weight(label)
    sN = mp.sqrt(N)
    tmin = min(mp.mpf(t), 1 / mp.mpf(t))
    X = (mp.mp.dps + extra_digits) * mp.log(10) + 20 + 3 * abs(mp.re(w))
    if k != 2:
        X = (mp.mp.dps + extra_digits) * mp.log(10) + 20 + 3 * max(abs(mp.re(w)), abs(k - mp.re(w)))
        n0 = int(X * sN / (2 * mp.pi * tmin)) + 2
        X += (mp.mpf(k + 1) / 2) * mp.log(n0)
    nmax = int(X * sN / (2 * mp.pi * tmin)) + 2
    if an is None or len(an) <= nmax:
        an = an_table(label, nmax)
    tot = mp.mpc(0)
    for n in range(1, nmax + 1):
        if an[n] == 0:
            continue
        c = 2 * mp.pi * n / sN
        tot += an[n] * (c ** (-w) * mp.gammainc(w, c * t) + eps * c ** (w - k) * mp.gammainc(k - w, c / t))
    return tot


def _work_dps(tau, digits=30):
    return digits + int(0.7 * abs(float(tau))) + 10


def hardy_Z(label, tau, an=None, digits=30):
    """Z(τ) = e^{iθ(τ)} L(½ + iτ), real for ε = +1: Λ_a(k/2 + iτ)/((√N/2π)^{k/2}|Γ(k/2 + iτ)|). With t = 1 the AFE is
    2 Re Σ a_n c_n^{−w} Γ(w, c_n), w = k/2 + iτ (k − w = w̄ on the line)."""
    N = N_of(label)
    k = weight(label)
    with mp.workdps(_work_dps(tau, digits)):
        tau = mp.mpf(tau)
        w = mp.mpc(mp.mpf(k) / 2, tau)
        sN = mp.sqrt(N)
        X = mp.mp.dps * mp.log(10) + 23
        if k != 2:
            X += 3 * mp.mpf(k) / 2 + (mp.mpf(k + 1) / 2) * mp.log(int(X * sN / (2 * mp.pi)) + 2)
        nmax = int(X * sN / (2 * mp.pi)) + 2
        if an is None or len(an) <= nmax:
            an = an_table(label, nmax)
        tot = mp.mpc(0)
        for n in range(1, nmax + 1):
            if an[n]:
                c = 2 * mp.pi * n / sN
                tot += an[n] * c ** (-w) * mp.gammainc(w, c)
        if k == 2 and eps_of(label) == 1:
            val = 2 * mp.re(tot) * 2 * mp.pi / (sN * abs(mp.gamma(w)))        # P-DL4 expression, unchanged
        elif k == 2:
            val = 2 * mp.im(tot) * 2 * mp.pi / (sN * abs(mp.gamma(w)))        # ε = −1: Z = Λ/(i·norm), odd in τ
        elif eps_of(label) == 1:
            val = 2 * mp.re(tot) / ((sN / (2 * mp.pi)) ** (mp.mpf(k) / 2) * abs(mp.gamma(w)))
        else:
            # ε = −1: c^{w−k} = conj(c^{−w}) on the line, so Λ = Σ a_n (X_n − X̄_n) = 2i Im Σ X_n; Z = Λ/(i·norm), odd in τ
            val = 2 * mp.im(tot) / ((sN / (2 * mp.pi)) ** (mp.mpf(k) / 2) * abs(mp.gamma(w)))
    return +val


def theta_E(label, T):
    N = N_of(label)
    return T * mp.log(mp.sqrt(N) / (2 * mp.pi)) + mp.im(mp.loggamma(mp.mpc(mp.mpf(weight(label)) / 2, T)))


def count_zeros_argument(label, T, sigma_hi=None, step=mp.mpf("0.1")):
    """N(T) = (θ(T) + arg L_a(k/2 + iT))/π, arg by continuous variation along w = σ + iT from σ = sigma_hi
    (default k/2 + 3.5, where |L_a − 1| < 1) down to σ = k/2. The result must be an integer; returns
    (value, max |Δarg| per step)."""
    N = N_of(label)
    k = weight(label)
    lo = mp.mpf(k) / 2
    if sigma_hi is None:
        sigma_hi = lo + mp.mpf("3.5")
    with mp.workdps(_work_dps(T, 25)):
        T = mp.mpf(T)
        sN = mp.sqrt(N)
        an = an_table(label, 2000)

        def La(sig):
            w = mp.mpc(sig, T)
            return afe_arith(label, w, an=an, eps=eps_of(label)) / ((sN / (2 * mp.pi)) ** w * mp.gamma(w))

        sig = mp.mpf(sigma_hi)
        v = La(sig)
        arg, maxstep = mp.arg(v), mp.mpf(0)
        while sig > lo:
            s2 = max(lo, sig - step)
            v2 = La(s2)
            d = mp.arg(v2 / v)
            if abs(d) > mp.mpf("0.5"):
                # refine this step
                sub = s2
                h = (sig - s2) / 8
                vv, cur = v, sig
                for _ in range(8):
                    cur -= h
                    vn = La(cur)
                    dd = mp.arg(vn / vv)
                    maxstep = max(maxstep, abs(dd))
                    arg += dd
                    vv = vn
                v, sig = vv, sub
                continue
            maxstep = max(maxstep, abs(d))
            arg += d
            v, sig = v2, s2
        return (theta_E(label, T) + arg) / mp.pi, maxstep


def cmd_afe(args):
    """AFE checks (weight k): (a) two splitting parameters; (b) Λ_a(w) = εΛ_a(k − w) with t ≠ 1; (c) vs the Dirichlet
    series at Re w = k/2 + 4 and k/2 + 3.5; (d) the wrong sign ε = −1 fails (b); (e) the central value L(½) ≠ 0."""
    job = Job("gl2 afe", args=vars(args))
    rng = random.Random(args.seed)
    report = {}
    for label in args.curves.split(","):
        N = N_of(label)
        k = weight(label)
        h = mp.mpf(k) / 2
        rep = {"split": [], "fe": [], "fe_wrong_sign": [], "dirichlet": []}
        mp.mp.dps = args.dps
        pts = [mp.mpc(rng.uniform(float(h) - 2.5, float(h) + 2.5), rng.uniform(-25, 25)) for _ in range(args.npts)]
        for w in pts:
            e = eps_of(label)
            with mp.workdps(_work_dps(mp.im(w), args.dps)):
                v1 = afe_arith(label, w, t=1, eps=e)
                v2 = afe_arith(label, w, t=mp.mpf("1.3"), eps=e)
                v3 = afe_arith(label, w, t=mp.mpf("0.75"), eps=e)
                fe = afe_arith(label, k - w, t=mp.mpf("1.3"), eps=e)
                wrong1 = afe_arith(label, w, t=mp.mpf("1.3"), eps=-e)
                wrong2 = afe_arith(label, k - w, t=mp.mpf("1.3"), eps=-e)
                r_split = max(abs(v2 - v1), abs(v3 - v1)) / abs(v1)
                r_fe = abs(fe - e * v2) / abs(v2)
                r_wrong = abs(wrong2 + e * wrong1) / abs(wrong1)    # the other sign would need Λ(k − w) = −εΛ(w)
            rep["split"].append(mp.nstr(r_split, 3))
            rep["fe"].append(mp.nstr(r_fe, 3))
            rep["fe_wrong_sign"].append(mp.nstr(r_wrong, 3))
            job.result(f"{label} w = {mp.nstr(w, 6)} (ε = {e:+d}): t = 1 vs 1.3, 0.75: rel {mp.nstr(r_split, 3)};  Λ({k} − w) vs εΛ(w) (t = 1.3): rel "
                       f"{mp.nstr(r_fe, 3)};  with the other sign ε = {-e:+d}: rel {mp.nstr(r_wrong, 3)}")
        # (c) Dirichlet series in its region of absolute convergence, independent of the functional equation
        mp.mp.dps = 30
        nds = 30000
        an = an_table(label, nds)
        for w in (mp.mpc(h + 4, 3), mp.mpc(h + 4, -11), mp.mpc(h + mp.mpf("3.5"), 20)):
            with mp.workdps(_work_dps(mp.im(w), 30)):
                ds = mp.fsum(an[n] * mp.power(n, -w) for n in range(1, nds + 1) if an[n])
                tail = mp.fsum(2 * mp.mpf(n) ** ((mp.mpf(k) - 1) / 2 - mp.re(w)) for n in range(nds + 1, nds + 2))   # first omitted term
                lam_ds = (mp.sqrt(N) / (2 * mp.pi)) ** w * mp.gamma(w) * ds
                lam = afe_arith(label, w, eps=eps_of(label))
                r = abs(lam - lam_ds) / abs(lam)
            rep["dirichlet"].append({"w": mp.nstr(w, 6), "rel": mp.nstr(r, 3)})
            job.result(f"{label} w = {mp.nstr(w, 6)}: AFE vs Dirichlet series (n ≤ {nds}, first omitted term ~{mp.nstr(tail, 2)}): rel {mp.nstr(r, 3)}")
        # (e) central value L_a(k/2) = L(½) (analytic)
        mp.mp.dps = 40
        if k == 2 and eps_of(label) == 1:
            L1 = afe_arith(label, mp.mpf(1)) / ((mp.sqrt(N) / (2 * mp.pi)))      # Γ(1) = 1
            L1b = afe_arith(label, mp.mpf(1), t=mp.mpf("1.3")) / ((mp.sqrt(N) / (2 * mp.pi)))
        else:
            e = eps_of(label)
            nrm = (mp.sqrt(N) / (2 * mp.pi)) ** h * mp.gamma(h)
            L1 = afe_arith(label, h, eps=e) / nrm
            L1b = afe_arith(label, h, t=mp.mpf("1.3"), eps=e) / nrm
            if e == -1:
                # L(½) = 0 is forced; the derivative (the central zero is simple iff L'(½) ≠ 0), two splitting parameters
                dd = mp.mpf("1e-12")
                nrm_d = lambda u: (mp.sqrt(N) / (2 * mp.pi)) ** u * mp.gamma(u)
                for tt, key in ((1, "Lprime_half"), (mp.mpf("1.3"), "Lprime_half_t1.3")):
                    der = (afe_arith(label, h + dd, t=tt, eps=e) / nrm_d(h + dd) - afe_arith(label, h - dd, t=tt, eps=e) / nrm_d(h - dd)) / (2 * dd)
                    rep[key] = mp.nstr(mp.re(der), 15)
                job.result(f"{label}: ε = −1, so L(½) = {mp.nstr(mp.re(L1), 5)} (forced 0); L'(½) [analytic] = {rep['Lprime_half']} "
                           f"(t = 1.3: {rep['Lprime_half_t1.3']})")
        rep["L_half"] = mp.nstr(mp.re(L1), 20)
        rep["L_half_t1.3"] = mp.nstr(mp.re(L1b), 20)
        job.result(f"{label}: L(½) [analytic] = L_a({mp.nstr(h, 3)}) [arithmetic] = {mp.nstr(mp.re(L1), 20)} (t = 1.3: {mp.nstr(mp.re(L1b), 20)})")
        report[label] = rep
    if args.json:
        with open(args.json, "w") as fh:
            json.dump(report, fh, indent=1)
    job.done()


def cmd_zeros(args):
    """Ordinates 0 < γ ≤ T of L(E, s) from sign changes of Z_E on a grid of step h, refined by findroot; counts checked
    against N(T) from the argument principle at several heights. With --tmin/--chunk, one τ-range of a split run."""
    label = args.curve
    job = Job(f"gl2 zeros {label}", total=None, args=vars(args))
    mp.mp.dps = 30
    an = an_table(label, 4000)
    h = mp.mpf(args.step)
    T = mp.mpf(args.tmax)
    zs = []
    t0 = mp.mpf(args.tmin) if args.tmin else h
    t, zt = t0, hardy_Z(label, t0, an)
    vals = [(t, zt)]
    while t < T:
        t2 = t + h
        z2 = hardy_Z(label, t2, an)
        vals.append((t2, z2))
        if zt * z2 < 0:
            # Z is accurate to ~30 digits (it raises its own precision by 0.7|τ| for the AFE cancellation), so the
            # root is refined at 30 digits; findroot's residual test is then on Z's real accuracy.
            try:
                r = mp.findroot(lambda u: hardy_Z(label, u, an), (t, t2), solver="anderson")
            except ValueError:
                a, b, fa = t, t2, zt
                while b - a > mp.mpf("1e-25"):
                    m = (a + b) / 2
                    fm = hardy_Z(label, m, an)
                    if fm * fa > 0:
                        a, fa = m, fm
                    else:
                        b = m
                r = (a + b) / 2
                job.warn(f"findroot failed near {mp.nstr(t, 8)}; bisection to 1e-25")
            zs.append(r)
            job.log(f"zero {len(zs)}: {mp.nstr(r, 15)}")
            if args.json:
                with open(args.json + ".partial", "w") as fh:
                    json.dump([mp.nstr(g, 25) for g in zs], fh)
        t, zt = t2, z2
        if int(t / h) % 100 == 0:
            job.sub(f"zeros of L({label})", int(t), int(T), every=1)
    job.result(f"{label}: {len(zs)} zeros in ({mp.nstr(t0 if args.tmin else 0, 6)}, {args.tmax}] with step {args.step}; "
               f"first {[mp.nstr(g, 10) for g in zs[:4]]}")
    if args.chunk:
        # one τ-chunk of a split run: save the zeros only; `zeros-merge` does the central check and the counts
        with open(args.json, "w") as fh:
            json.dump({"curve": label, "tmin": args.tmin, "tmax": args.tmax, "step": args.step,
                       "zeros": [mp.nstr(g, 25) for g in zs]}, fh, indent=1)
        job.done()
        return
    _finish_zeros(label, zs, args, job, an, T)


def cmd_zeros_merge(args):
    """Merge τ-chunks from `zeros --chunk` (consecutive [tmin, tmax] ranges on the same step), then run the central-zero
    check and the argument-principle counts exactly as `zeros` does, and write the usual zeros_<label>.json."""
    label = args.curve
    job = Job(f"gl2 zeros-merge {label}", total=None, args=vars(args))
    mp.mp.dps = 30
    parts = sorted((json.load(open(f)) for f in args.parts.split(",")), key=lambda d: float(d["tmin"] or 0))
    for a, b in zip(parts, parts[1:]):
        assert float(a["tmax"]) == float(b["tmin"]), "chunks must be contiguous"
        assert a["step"] == b["step"]
    raw = sorted(mp.mpf(g) for d in parts for g in d["zeros"])
    # A chunk's grid may overrun its tmax by one step (accumulated rounding of t + h), so a zero just above a chunk
    # boundary can be found by both neighbours: identical to ~1e-25, removed here and counted.
    zs = [g for i, g in enumerate(raw) if i == 0 or g - raw[i - 1] > mp.mpf("1e-12")]
    if len(zs) != len(raw):
        job.log(f"removed {len(raw) - len(zs)} duplicate zero(s) at chunk boundaries")
    T = mp.mpf(parts[-1]["tmax"])
    zs = [g for g in zs if g <= T]
    args.tmax, args.step = float(T), parts[0]["step"]
    job.result(f"{label}: merged {len(parts)} chunks → {len(zs)} zeros in (0, {args.tmax}]; first {[mp.nstr(g, 10) for g in zs[:4]]}")
    an = an_table(label, 4000)
    _finish_zeros(label, zs, args, job, an, T)


def _finish_zeros(label, zs, args, job, an, T):
    # Central zero. The argument principle on the symmetric rectangle gives (θ + arg L)/π = N(T) + m0/2, m0 the
    # order at s = ½ (0 for ε = +1 with L(½) ≠ 0; odd for ε = −1). Directly: Z(τ) ~ Z^{(m0)}(0) τ^{m0}/m0! near 0.
    m0_expected = m0_of(label)
    central = {}
    if m0_expected:
        z1, z2 = hardy_Z(label, mp.mpf("1e-4"), an), hardy_Z(label, mp.mpf("2e-4"), an)
        central = {"Z(1e-4)": mp.nstr(z1, 10), "Z(2e-4)": mp.nstr(z2, 10), "ratio": mp.nstr(z2 / z1, 10),
                   "expected_ratio": 2 ** m0_expected, "Zprime0_est": mp.nstr(z1 / mp.mpf("1e-4"), 10)}
        job.result(f"{label}: central zero: Z(2e-4)/Z(1e-4) = {mp.nstr(z2 / z1, 10)} (2^m0 = {2 ** m0_expected} expected: "
                   f"2 for a simple zero, 4 for a double, 8 for a triple)")
    # argument-principle counts at heights midway between found zeros
    checks = []
    for Tc in [float(v) for v in args.count_at.split(",")]:
        if Tc > float(T):
            continue
        below = [g for g in zs if g <= Tc]
        above = [g for g in zs if g > Tc]
        if below and above:
            Tm = (below[-1] + above[0]) / 2
        else:
            Tm = mp.mpf(Tc)
        nT, mx = count_zeros_argument(label, Tm)
        cnt = sum(1 for g in zs if g <= Tm)
        ok = abs(nT - cnt - mp.mpf(m0_expected) / 2) < mp.mpf("1e-6")
        rec = {"T": mp.nstr(Tm, 10), "N_argument": mp.nstr(nT, 12), "found": cnt, "max_darg": mp.nstr(mx, 3), "ok": ok}
        if m0_expected:
            rec["central_multiplicity_from_count"] = mp.nstr(2 * (nT - cnt), 8)
        checks.append(rec)
        job.result(f"{label}: (θ + arg L)/π at T = {mp.nstr(Tm, 8)}: {mp.nstr(nT, 12)}; found {cnt} zeros in (0, T]"
                   + (f", so m0 = 2·({mp.nstr(nT, 10)} − {cnt}) = {mp.nstr(2 * (nT - cnt), 8)}" if m0_expected else "")
                   + f" → {'OK' if ok else 'MISMATCH'} (max |Δarg| per step {mp.nstr(mx, 2)})")
    if args.json:
        with open(args.json, "w") as fh:
            json.dump({"curve": label, "N": N_of(label), "weight": weight(label), "eps": eps_of(label), "tmax": args.tmax,
                       "step": args.step, "central_multiplicity": m0_expected if all(c["ok"] for c in checks) else None,
                       "central_check": central, "zeros": [mp.nstr(g, 25) for g in zs], "count_checks": checks}, fh, indent=1)
    job.done()


# ---------------------------------------------------------------------------------------------------------------
# Hard validation 3: the form against the zero sum
# ---------------------------------------------------------------------------------------------------------------

def tail_estimate(label, Ffun, T, L):
    """Expected truncation error 2∫_T^∞ F(t)² dN(t), with the smooth density dN/dt = θ_E'(t)/π =
    (log(√N/2π) + Re ψ(1 + it))/π, integrated over [T, 40T] in pieces of length 2π/L (the tail beyond is
    O(40⁻³) of this)."""
    N = N_of(label)
    h = mp.mpf(weight(label)) / 2
    with mp.workdps(20):
        dens = lambda t: (mp.log(mp.sqrt(N) / (2 * mp.pi)) + mp.re(mp.digamma(mp.mpc(h, t)))) / mp.pi
        pts = mp.linspace(mp.mpf(T), 40 * mp.mpf(T), int(39 * T * L / (2 * mp.pi)) + 2)
        return 2 * mp.quad(lambda t: Ffun(t) ** 2 * dens(t), pts)


def cmd_check(args):
    """Q(c) = cᵀ M c against 2 Σ_{0<γ≤T} F_c(γ)² for edge-vanishing test functions (as decay_law_mp check):
    odd: single sin modes k = 2, 4 and b_2 − b_3; even: b_1 + b_2 and b_2 − b_4."""
    job = Job("gl2 check", args=vars(args))
    out = []
    for label in args.curves.split(","):
        zd = json.load(open(os.path.join(OUT, f"zeros_{label}{args.zeros_tag}.json")))
        assert all(c["ok"] for c in zd["count_checks"]), f"zero count check failed for {label}"
        assert (zd.get("central_multiplicity") or 0) == m0_of(label), f"central multiplicity not established for {label}"
        for xs in args.x.split(","):
            mp.mp.dps = 30
            x = mp.mpf(xs)
            L = mp.log(x)
            gam = [mp.mpf(g) for g in zd["zeros"]]
            T = mp.mpf(zd["tmax"])
            cases = []
            Mo = zeros_side_E(label, x, 6, "odd")
            for k in (2, 4):
                c = [0] * 6
                c[k - 1] = 1
                cases.append(("odd", c, Mo))
            cases.append(("odd", [1 if i == 1 else -1 if i == 2 else 0 for i in range(6)], Mo))
            Me = zeros_side_E(label, x, 6, "even")
            cases += [("even", [0, 1, 1, 0, 0, 0, 0], Me), ("even", [0, 0, 1, 0, -1, 0, 0], Me)]
            for par, c, M in cases:
                if par == "odd":
                    Ff = lambda t, c=c: mp.fsum(ci * dl.F_per_odd(i + 1, L, t) for i, ci in enumerate(c) if ci)
                else:
                    Ff = lambda t, c=c: dl.F_per_even(c, L, t)
                cv = mp.matrix(c)
                q = (cv.T * M * cv)[0]
                fz = mp.fsum(2 * Ff(g) ** 2 for g in gam) + (zd.get("central_multiplicity") or 0) * Ff(mp.mpf(0)) ** 2
                tail = tail_estimate(label, Ff, T, L)
                rel = abs(q - fz) / abs(q)
                rel_t = abs(q - fz - tail) / abs(q)
                out.append({"curve": label, "x": xs, "parity": par, "c": c, "Q": mp.nstr(q, 15), "zero_sum": mp.nstr(fz, 15),
                            "T": mp.nstr(T, 6), "n_zeros": len(gam), "abs_diff": mp.nstr(q - fz, 4), "rel_diff": mp.nstr(rel, 4),
                            "expected_tail": mp.nstr(tail, 4), "rel_diff_after_tail": mp.nstr(rel_t, 4)})
                job.result(f"{label} x={xs} {par} c={c}: cᵀQc = {mp.nstr(q, 12)};  2Σ_{{γ≤{mp.nstr(T, 4)}}} F(γ)² = {mp.nstr(fz, 12)};  "
                           f"Q − Σ = {mp.nstr(q - fz, 3)} (rel {mp.nstr(rel, 3)});  expected tail {mp.nstr(tail, 3)};  "
                           f"rel after tail {mp.nstr(rel_t, 3)}")
    if args.gate:
        # P-DL6 zero-sum gate (docs/GL2_DECAY.md §9, fixed before any data): at x = 13, every test has relative error
        # ≤ 1e-4 and (Q − Σ)/tail ∈ [0.8, 1.25].
        gate = {}
        for label in args.curves.split(","):
            rs = [r for r in out if r["curve"] == label and r["x"] == "13"]
            ok = bool(rs) and all(mp.mpf(r["rel_diff"]) <= mp.mpf("1e-4") and
                                  mp.mpf("0.8") <= mp.mpf(r["abs_diff"]) / mp.mpf(r["expected_tail"]) <= mp.mpf("1.25") for r in rs)
            gate[label] = ok
            ratios = [mp.nstr(mp.mpf(r["abs_diff"]) / mp.mpf(r["expected_tail"]), 3) for r in rs]
            job.result(f"P-DL6 zero-sum gate {label} (x = 13): max rel {mp.nstr(max(mp.mpf(r['rel_diff']) for r in rs), 3)}, "
                       f"(Q − Σ)/tail = {ratios} → {'PASS' if ok else 'FAIL'}")
        out = {"rows": out, "gate": gate}
    if args.gate7:
        # P-DL7 gate (docs/DECAY_LAW.md, 6fd392d): at x = 13 with zeros to T = 240, every test has (Q − Σ)/tail ∈ [0.8, 1.25]
        # and |Q − Σ| ≤ 2e-4.
        gate = {}
        for label in args.curves.split(","):
            rs = [r for r in out if r["curve"] == label and r["x"] == "13"]
            T_ok = bool(rs) and all(mp.mpf(r["T"]) >= 240 for r in rs)
            ok = T_ok and all(abs(mp.mpf(r["abs_diff"])) <= mp.mpf("2e-4") and
                              mp.mpf("0.8") <= mp.mpf(r["abs_diff"]) / mp.mpf(r["expected_tail"]) <= mp.mpf("1.25") for r in rs)
            gate[label] = ok
            ratios = [mp.nstr(mp.mpf(r["abs_diff"]) / mp.mpf(r["expected_tail"]), 3) for r in rs]
            job.result(f"P-DL7 gate {label} (x = 13, T = {rs[0]['T'] if rs else '?'}): max |Q − Σ| = "
                       f"{mp.nstr(max(abs(mp.mpf(r['abs_diff'])) for r in rs), 3)}, (Q − Σ)/tail = {ratios} → {'PASS' if ok else 'FAIL'}")
        out = {"rows": out, "gate7": gate}
    if args.json:
        with open(args.json, "w") as fh:
            json.dump(out, fh, indent=1)
    job.done()


# ---------------------------------------------------------------------------------------------------------------
# P-DL4 scan and fit
# ---------------------------------------------------------------------------------------------------------------

def lam_min_inv_checked(M, tol=mp.mpf("1e-14"), maxit=40):
    """Inverse iteration as decay_law_mp.lam_min_inv (FLINT ball solve, midpoints), but iterated to a relative
    Rayleigh-quotient change ≤ tol (lam_min_inv stops silently after 6). If that is not reached it returns NaN, which
    stable_lams can never accept, so the round is retried at higher precision."""
    from flint import arb, arb_mat, ctx
    ctx.prec = int(mp.mp.dps * 3.33) + 20
    n = M.rows
    A = arb_mat([[arb(mp.nstr(M[i, j], mp.mp.dps + 5)) for j in range(n)] for i in range(n)])
    v = arb_mat([[1] for _ in range(n)])
    prev = None
    for it in range(maxit):
        w = A.solve(v)
        nrm = sum((w[i, 0] ** 2 for i in range(n)), arb(0)).sqrt()
        v = arb_mat([[w[i, 0] / nrm] for i in range(n)])
        lam = mp.mpf((v.transpose() * A * v)[0, 0].mid().str(mp.mp.dps, radius=False))
        if prev is not None and abs(lam - prev) <= abs(lam) * tol:
            return lam
        prev = lam
    print(f"WARN inverse iteration not converged in {maxit} steps (n = {n}, dps {mp.mp.dps}); returning NaN", flush=True)
    return mp.nan


def cmd_scan(args):
    """P-DL4 as registered: x/N ∈ grid, both sectors, N_basis = int(f k*) + 8 for f = 5, 9 (k* = xL/N, the
    convention of decay_law_mp scan/product), every λ precision-stable (decay_law_mp.stable_lams: rebuild at
    dps + 40 until two consecutive minima agree to 1e-8; the lower-precision value is kept)."""
    label = args.curve
    N = CURVES[label]["N"]
    xqs = [mp.mpf(s) for s in args.xq.split(",")]
    rows = []
    path = args.json or os.path.join(OUT, f"scan_{label}.json")
    if os.path.exists(path) and args.resume:
        rows = json.load(open(path))
    done = {(r["xq"], r["parity"]) for r in rows}
    job = Job(f"gl2 scan {label}", total=len(xqs) * 2, args=vars(args))
    for xq in xqs:
        x = xq * N
        L = mp.log(x)
        kstar = float(x * L / N)
        for parity in ("even", "odd"):
            if (mp.nstr(xq, 10), parity) in done:
                job.step(f"x/N={mp.nstr(xq, 4)} {parity}: already done")
                continue
            ns = [max(args.nmin, int(f * kstar) + 8) for f in (args.f1, args.f2)]
            res = {}
            for n in ns:
                t0 = time.time()

                def build(n=n):
                    return {"E": zeros_side_E(label, x, n, parity)}

                dps0 = int(40 + 1.5 * 4 * float(mp.pi) * float(xq) / 2.303)
                lam, at = dl.stable_lams(build, dps0, method=lam_min_inv_checked)
                rec = {"n": n, "lambda": lam["E"], "dps": at["E"], "secs": time.time() - t0}
                if n <= args.acb_max:
                    mp.mp.dps = at["E"]
                    la = dl.lam_min_acb(build()["E"])
                    rec["acb"] = la
                    rec["acb_rel"] = abs(la - lam["E"]) / abs(lam["E"])
                res[n] = rec
                job.log(f"x/N={mp.nstr(xq, 4)} x={mp.nstr(x, 6)} {parity} N_basis={n}: λ = {mp.nstr(lam['E'], 12)} stable at dps {at['E']} "
                        f"({rec['secs']:.0f}s)" + (f"; acb {mp.nstr(rec['acb'], 12)} (rel {mp.nstr(rec['acb_rel'], 2)})" if "acb" in rec else ""))
            r1, r2 = res[ns[0]], res[ns[1]]
            conv = mp.log(r1["lambda"] / r2["lambda"])
            job.step(f"x/N={mp.nstr(xq, 4)} x={mp.nstr(x, 6)} {parity}: N={ns[0]}: {mp.nstr(r1['lambda'], 8)}  N={ns[1]}: {mp.nstr(r2['lambda'], 8)}  "
                     f"ln(λ1/λ2) = {mp.nstr(conv, 3)}  (dps {r1['dps']}, {r2['dps']})")
            row = {"curve": label, "N": N, "xq": mp.nstr(xq, 10), "x": mp.nstr(x, 12), "parity": parity, "kstar": kstar,
                   "n1": ns[0], "lambda_n1": mp.nstr(r1["lambda"], 15), "dps_n1": r1["dps"],
                   "n2": ns[1], "lambda_n2": mp.nstr(r2["lambda"], 15), "dps_n2": r2["dps"],
                   "ln_l1_over_l2": mp.nstr(conv, 6), "secs": round(r1["secs"] + r2["secs"], 1)}
            for key, r in (("n1", r1), ("n2", r2)):
                if "acb" in r:
                    row[f"acb_{key}"] = mp.nstr(r["acb"], 15)
                    row[f"acb_rel_{key}"] = mp.nstr(r["acb_rel"], 3)
            rows.append(row)
            with open(path, "w") as fh:
                json.dump(rows, fh, indent=1)
    job.done()


FALLBACKS = []


def lam_min_inv_or_full(M):
    """Inverse iteration (lam_min_inv_checked); if it does not converge, which happens when λ_2/λ_1 is close to 1
    (pre-asymptotic points where λ_min = O(1)), FLINT's full eigensolver (decay_law_mp.lam_min_acb, as used for
    DH). Each fallback is recorded in FALLBACKS."""
    v = lam_min_inv_checked(M)
    if mp.isnan(v):
        FALLBACKS.append((M.rows, mp.mp.dps))
        v = dl.lam_min_acb(M)
    return v


def cmd_scan6(args):
    """P-DL6 as registered: v = √(x/N) on the grid, both sectors, N_basis = ⌊f·k⌋ + 8 for f = 5, 9 with k = 2vL
    (nmin 24, n2 ≥ n1 + 8, as decay_law_mp product-v), every λ precision-stable (stable_lams, dps0 = 40 + 1.5·8πv/ln 10
    as product-v), inverse iteration to 1e-14 (NaN on non-convergence, never accepted); FLINT's full eigensolver
    as a PSD check for N_basis ≤ acb_max."""
    label = args.obj
    N = N_of(label)
    vs = [mp.mpf(s) for s in args.v.split(",")]
    path = args.json or os.path.join(OUT, f"scan6_{label}.json")
    rows = json.load(open(path)) if (args.resume and os.path.exists(path)) else []
    done = {(r["v"], r["parity"]) for r in rows}
    job = Job(f"gl2 scan6 {label}", total=len(vs) * 2, args=vars(args))
    for v in vs:
        x = v * v * N
        L = mp.log(x)
        kk = float(2 * v * L)
        for parity in ("even", "odd"):
            if (mp.nstr(v, 10), parity) in done:
                job.step(f"v={mp.nstr(v, 4)} {parity}: already done")
                continue
            ns = [max(args.nmin, int(f * kk) + 8) for f in (args.f1, args.f2)]
            ns[1] = max(ns[1], ns[0] + 8)
            res = {}
            for n in ns:
                t0 = time.time()

                def build(n=n):
                    return {"L": zeros_side_E(label, x, n, parity)}

                dps0 = int(40 + 1.5 * 8 * float(mp.pi) * float(v) / 2.303)
                del FALLBACKS[:]
                method = lam_min_inv_or_full if args.method == "inv+full" else lam_min_inv_checked
                lam, at = dl.stable_lams(build, dps0, method=method)
                rec = {"n": n, "lambda": lam["L"], "dps": at["L"], "secs": time.time() - t0, "fallbacks": len(FALLBACKS)}
                if n <= args.acb_max:
                    mp.mp.dps = at["L"]
                    la = dl.lam_min_acb(build()["L"])
                    rec["acb"] = la
                    rec["acb_rel"] = abs(la - lam["L"]) / abs(lam["L"])
                res[n] = rec
                job.log(f"v={mp.nstr(v, 4)} x={mp.nstr(x, 8)} {parity} N_basis={n}: λ = {mp.nstr(lam['L'], 12)} stable at dps {at['L']} "
                        f"({rec['secs']:.0f}s)" + (f"; acb {mp.nstr(rec['acb'], 12)} (rel {mp.nstr(rec['acb_rel'], 2)})" if "acb" in rec else ""))
            r1, r2 = res[ns[0]], res[ns[1]]
            conv = mp.log(r1["lambda"] / r2["lambda"])
            job.step(f"v={mp.nstr(v, 4)} x={mp.nstr(x, 8)} {parity}: N={ns[0]}: {mp.nstr(r1['lambda'], 8)}  N={ns[1]}: {mp.nstr(r2['lambda'], 8)}  "
                     f"ln(λ1/λ2) = {mp.nstr(conv, 3)}  (dps {r1['dps']}, {r2['dps']})")
            row = {"obj": label, "N": N, "weight": weight(label), "v": mp.nstr(v, 10), "x": mp.nstr(x, 15), "parity": parity, "k": kk,
                   "n1": ns[0], "lambda_n1": mp.nstr(r1["lambda"], 15), "dps_n1": r1["dps"],
                   "n2": ns[1], "lambda_n2": mp.nstr(r2["lambda"], 15), "dps_n2": r2["dps"],
                   "ln_l1_over_l2": mp.nstr(conv, 6), "secs": round(r1["secs"] + r2["secs"], 1),
                   "method": args.method, "full_eigensolver_fallbacks_n1": r1["fallbacks"], "full_eigensolver_fallbacks_n2": r2["fallbacks"]}
            for key, r in (("n1", r1), ("n2", r2)):
                if "acb" in r:
                    row[f"acb_{key}"] = mp.nstr(r["acb"], 15)
                    row[f"acb_rel_{key}"] = mp.nstr(r["acb_rel"], 3)
            rows.append(row)
            with open(path, "w") as fh:
                json.dump(rows, fh, indent=1)
    job.done()


PDL6_PRED = {"15a1": (2, 4), "19a1": (2, 4), DELTA: (12, 14)}


def cmd_fit6(args):
    """P-DL6 as registered (docs/DECAY_LAW.md, f5697aa).
    6a: ln λ = −αv + γ ln v + β over the five grid points (λ at the 9k basis); prediction α/4π ∈ [1.7, 2.3] for every
        object and sector; kill if any α/4π is outside [1.4, 2.6].
    6b: c = 4πv, n ∈ 0..20, flattest = argmin_n |ln R(v₅)/R(v₂)|, R = λ/ℓ_n(c); predicted n = 2/4 (curves), 12/14 (Δ);
        success iff ≥ 4 of 6 exact and all 6 within ±1; kill if Δ has flattest n ≤ 8 in either sector, or any curve
        case is off by ≥ 2."""
    mp.mp.dps = 30
    four_pi = 4 * mp.pi
    out = {"6a": [], "6b": [], "6a_5k": []}
    for label in args.objects.split(","):
        rows = json.load(open(os.path.join(args.dir, f"scan6_{label}.json")))
        for si, par in enumerate(("even", "odd")):
            R = sorted((r for r in rows if r["parity"] == par), key=lambda r: float(r["v"]))
            v = [mp.mpf(r["v"]) for r in R]
            for key, dest in (("lambda_n2", "6a"), ("lambda_n1", "6a_5k")):
                lam = [mp.mpf(r[key]) for r in R]
                y = [mp.log(l) for l in lam]
                c, res = dl._lsq([[-vi, mp.log(vi), 1] for vi in v], y)
                cl_, rl = dl._lsq([[-vi, 1] for vi in v], y)
                a4 = c[0] / four_pi
                status = "in [1.7, 2.3]" if mp.mpf("1.7") <= a4 <= mp.mpf("2.3") else ("in kill-free band" if mp.mpf("1.4") <= a4 <= mp.mpf("2.6") else "KILL")
                rec = {"obj": label, "parity": par, "basis": "9k" if key == "lambda_n2" else "5k", "alpha_over_4pi": mp.nstr(a4, 6),
                       "gamma": mp.nstr(c[1], 6), "beta": mp.nstr(c[2], 6), "max_resid": mp.nstr(res, 4),
                       "linear_alpha_over_4pi": mp.nstr(cl_[0] / four_pi, 6), "linear_max_resid": mp.nstr(rl, 4), "status": status,
                       "conv_ln_l1_l2": [r["ln_l1_over_l2"] for r in R]}
                out[dest].append(rec)
                print(f"6a [{rec['basis']}] {label:<5} {par:<4}: α/4π = {mp.nstr(a4, 5):>7}  γ = {mp.nstr(c[1], 4):>7}  β = {mp.nstr(c[2], 5):>8}  "
                      f"max resid {mp.nstr(res, 3):>6};  linear α/4π = {mp.nstr(cl_[0] / four_pi, 5)}  → {status}")
            lam = [mp.mpf(r["lambda_n2"]) for r in R]
            ratios = {m: [l / dl._ell(m, four_pi * vi) for l, vi in zip(lam, v)] for m in range(21)}
            metric = {m: abs(mp.log(r[4] / r[1])) for m, r in ratios.items()}
            flat = min(metric, key=lambda m: metric[m])
            pred = PDL6_PRED[label][si]
            second = sorted(metric, key=lambda m: metric[m])[1]
            rec = {"obj": label, "parity": par, "predicted_n": pred, "flattest_n": flat, "off_by": flat - pred,
                   "metric": {m: mp.nstr(metric[m], 4) for m in metric}, "runner_up": second,
                   "R_flattest": [mp.nstr(t, 4) for t in ratios[flat]], "R_predicted": [mp.nstr(t, 4) for t in ratios[pred]]}
            out["6b"].append(rec)
            print(f"6b {label:<5} {par:<4}: flattest n = {flat} (metric {mp.nstr(metric[flat], 3)}; runner-up {second}: {mp.nstr(metric[second], 3)}), "
                  f"predicted {pred} (metric {mp.nstr(metric[pred], 3)}) → off by {flat - pred};  R(flattest) = {rec['R_flattest']}")
    a = [mp.mpf(r["alpha_over_4pi"]) for r in out["6a"]]
    a_hold = all(mp.mpf("1.7") <= t <= mp.mpf("2.3") for t in a)
    a_kill = any(not (mp.mpf("1.4") <= t <= mp.mpf("2.6")) for t in a)
    out["verdict_6a"] = "KILLED" if a_kill else ("HOLDS" if a_hold else "NOT KILLED, prediction partly missed")
    exact = sum(1 for r in out["6b"] if r["off_by"] == 0)
    within1 = all(abs(r["off_by"]) <= 1 for r in out["6b"])
    kill_b = any(r["obj"] == DELTA and r["flattest_n"] <= 8 for r in out["6b"]) or any(r["obj"] != DELTA and abs(r["off_by"]) >= 2 for r in out["6b"])
    out["verdict_6b"] = "KILLED" if kill_b else ("HOLDS" if exact >= 4 and within1 else "NOT KILLED, prediction partly missed")
    out["6b_exact"] = exact
    out["6b_all_within_1"] = within1
    print(f"P-DL6a: {out['verdict_6a']} (α/4π range {mp.nstr(min(a), 4)}–{mp.nstr(max(a), 4)});  "
          f"P-DL6b: {out['verdict_6b']} ({exact}/6 exact, all within ±1: {within1})")
    if args.json:
        with open(args.json, "w") as fh:
            json.dump(out, fh, indent=1)


def cmd_precision(args):
    """Recompute every saved 9k* minimum at its accepted dps + extra (stable_lams compared dps against dps + 40;
    this is a further, independent rise) and report the relative change."""
    label = args.curve
    src = args.file or os.path.join(OUT, f"scan_{label}.json")
    rows = json.load(open(src))
    dest = args.out or os.path.join(OUT, f"precision_{label}.json")
    out = []
    job = Job(f"gl2 precision {label}", total=len(rows), args=vars(args))
    for r in rows:
        mp.mp.dps = r["dps_n2"] + args.extra
        M = zeros_side_E(label, mp.mpf(r["x"]), r["n2"], r["parity"])
        v = lam_min_inv_or_full(M) if r.get("method") == "inv+full" else lam_min_inv_checked(M)
        old = mp.mpf(r["lambda_n2"])
        rel = abs(v - old) / abs(v)
        grid = f"x/N={r['xq']}" if "xq" in r else f"v={r['v']}"
        job.step(f"{label} {grid} {r['parity']} N={r['n2']}: dps {r['dps_n2']} → {mp.mp.dps}: {mp.nstr(old, 12)} → {mp.nstr(v, 12)}  rel {mp.nstr(rel, 3)}")
        out.append({"curve": label, **({"xq": r["xq"]} if "xq" in r else {"v": r["v"]}), "parity": r["parity"], "n2": r["n2"],
                    "dps": r["dps_n2"], "dps_recheck": mp.mp.dps,
                    "lambda_saved": r["lambda_n2"], "lambda_recheck": mp.nstr(v, 15), "rel_change": mp.nstr(rel, 4)})
        with open(dest, "w") as fh:
            json.dump(out, fh, indent=1)
    job.done()


HYPOTHESES = [("A'", "x/N"), ("C'", "x/sqrt(N)"), ("B'", "sqrt(x/N)")]


def hyp_var(h, xq, N):
    return {"A'": xq, "C'": xq * mp.sqrt(N), "B'": mp.sqrt(xq)}[h]


def cmd_fit(args):
    """P-DL4 rule: per curve and sector, fit ln λ = −α v + γ ln v + β over the five grid points (λ at N_basis = 9k*);
    a hypothesis survives there iff α/4π ∈ [0.8, 1.2]. Also the same fit at 5k* (robustness, not the registered
    value) and the exploratory flattest prolate index n ∈ 0..9 of λ/ℓ_n(2πv)."""
    mp.mp.dps = 30
    four_pi = 4 * mp.pi
    out = {"fits": [], "prolate": []}
    survive = {}
    for label in args.curves.split(","):
        N = CURVES[label]["N"]
        rows = json.load(open(os.path.join(args.dir, f"scan_{label}.json")))
        for par in ("even", "odd"):
            R = sorted((r for r in rows if r["parity"] == par), key=lambda r: float(r["xq"]))
            xq = [mp.mpf(r["xq"]) for r in R]
            for key in ("lambda_n2", "lambda_n1"):
                y = [mp.log(mp.mpf(r[key])) for r in R]
                for h, vname in HYPOTHESES:
                    v = [hyp_var(h, u, N) for u in xq]
                    c, res = dl._lsq([[-vi, mp.log(vi), 1] for vi in v], y)
                    a4 = c[0] / four_pi
                    ok = mp.mpf("0.8") <= a4 <= mp.mpf("1.2")
                    resid = [float(-c[0] * vi + c[1] * mp.log(vi) + c[2] - yi) for vi, yi in zip(v, y)]
                    rec = {"curve": label, "parity": par, "basis": "9k*" if key == "lambda_n2" else "5k*", "hypothesis": h, "v": vname,
                           "alpha": mp.nstr(c[0], 8), "alpha_over_4pi": mp.nstr(a4, 6), "gamma": mp.nstr(c[1], 6), "beta": mp.nstr(c[2], 6),
                           "max_resid": mp.nstr(res, 4), "resid": [round(t, 4) for t in resid], "survives": bool(ok)}
                    out["fits"].append(rec)
                    if key == "lambda_n2":
                        survive[(label, par, h)] = ok
                    print(f"{label} {par:<4} [{rec['basis']}] {h:<2} v = {vname:<10}: α/4π = {mp.nstr(a4, 5):>8}  γ = {mp.nstr(c[1], 4):>7}  "
                          f"β = {mp.nstr(c[2], 5):>8}  max resid {mp.nstr(res, 3):>8}  → {'survives' if ok else 'killed'}")
    # exploratory: flattest prolate index under each scaling (the surviving one is flagged)
    for label in args.curves.split(","):
        N = CURVES[label]["N"]
        rows = json.load(open(os.path.join(args.dir, f"scan_{label}.json")))
        for par in ("even", "odd"):
            R = sorted((r for r in rows if r["parity"] == par), key=lambda r: float(r["xq"]))
            xq = [mp.mpf(r["xq"]) for r in R]
            lam = [mp.mpf(r["lambda_n2"]) for r in R]
            for h, vname in HYPOTHESES:
                v = [hyp_var(h, u, N) for u in xq]
                ratios = {m: [l / dl._ell(m, 2 * mp.pi * vi) for l, vi in zip(lam, v)] for m in range(10)}
                metric = {m: abs(mp.log(r[-1] / r[0])) for m, r in ratios.items()}            # x/N = 12 against 2
                metric4 = {m: abs(mp.log(r[-1] / r[1])) for m, r in ratios.items()}           # x/N = 12 against 4
                flat = min(metric, key=lambda m: metric[m])
                flat4 = min(metric4, key=lambda m: metric4[m])
                rec = {"curve": label, "parity": par, "hypothesis": h, "survives": survive[(label, par, h)], "flattest_n": flat,
                       "flattest_n_4_12": flat4, "metric_ln_R12_over_R2": {m: mp.nstr(metric[m], 4) for m in metric},
                       "metric_ln_R12_over_R4": {m: mp.nstr(metric4[m], 4) for m in metric4},
                       "R_flattest": [mp.nstr(t, 4) for t in ratios[flat]]}
                out["prolate"].append(rec)
                print(f"prolate {label} {par:<4} {h:<2} ({'survives' if rec['survives'] else 'killed'}): flattest n = {flat} "
                      f"(|ln R(12)/R(2)| = {mp.nstr(metric[flat], 3)}; 4 vs 12: n = {flat4}), R = {rec['R_flattest']}")
    out["verdict"] = {f"{k[0]} {k[1]} {k[2]}": ("survives" if v else "killed") for k, v in survive.items()}
    if args.json:
        with open(args.json, "w") as fh:
            json.dump(out, fh, indent=1)


def cmd_posthoc(args):
    """NOT REGISTERED (found after the P-DL4 verdict): local slopes of ln λ in √(x/N) and in x/N, the collapse of the
    two curves at equal x/N, the B′-variable fit with the rate fixed at 8π, and the flattest prolate index for
    c = 4π√(x/N) (twice B′'s c)."""
    mp.mp.dps = 30
    four_pi = 4 * mp.pi
    out = {"note": "post hoc, not registered", "rows": []}
    lam = {}
    for label in args.curves.split(","):
        rows = json.load(open(os.path.join(args.dir, f"scan_{label}.json")))
        for par in ("even", "odd"):
            R = sorted((r for r in rows if r["parity"] == par), key=lambda r: float(r["xq"]))
            xq = [mp.mpf(r["xq"]) for r in R]
            ls = [mp.mpf(r["lambda_n2"]) for r in R]
            lam[(label, par)] = dict(zip([r["xq"] for r in R], ls))
            v = [mp.sqrt(u) for u in xq]
            y = [mp.log(l) for l in ls]
            s_sqrt = [-(y[i + 1] - y[i]) / (v[i + 1] - v[i]) / four_pi for i in range(len(v) - 1)]
            s_lin = [-(y[i + 1] - y[i]) / (xq[i + 1] - xq[i]) / four_pi for i in range(len(v) - 1)]
            c2, r2 = dl._lsq([[mp.log(vi), 1] for vi in v], [yi + 8 * mp.pi * vi for yi, vi in zip(y, v)])
            rat = {m: [l / dl._ell(m, 4 * mp.pi * vi) for l, vi in zip(ls, v)] for m in range(10)}
            met = {m: abs(mp.log(r[-1] / r[0])) for m, r in rat.items()}
            met4 = {m: abs(mp.log(r[-1] / r[1])) for m, r in rat.items()}
            f2, f4 = min(met, key=met.get), min(met4, key=met4.get)
            rec = {"curve": label, "parity": par, "local_slope_over_4pi_sqrt": [mp.nstr(s, 4) for s in s_sqrt],
                   "local_slope_over_4pi_linear": [mp.nstr(s, 4) for s in s_lin],
                   "fixed_8pi_sqrt_fit": {"gamma": mp.nstr(c2[0], 4), "beta": mp.nstr(c2[1], 5), "max_resid": mp.nstr(r2, 3)},
                   "ln_lambda_plus_8pi_sqrt": [round(float(yi + 8 * mp.pi * vi), 3) for yi, vi in zip(y, v)],
                   "prolate_c_4pi_sqrt": {"flattest_2_12": f2, "metric": mp.nstr(met[f2], 3), "flattest_4_12": f4,
                                          "R": [mp.nstr(t, 4) for t in rat[f2]]}}
            out["rows"].append(rec)
            print(f"{label} {par:<4}: local slope/4π in √(x/N) {rec['local_slope_over_4pi_sqrt']};  in x/N {rec['local_slope_over_4pi_linear']}")
            print(f"           ln λ + 8π√(x/N) = {rec['ln_lambda_plus_8pi_sqrt']};  fixed-8π fit γ = {mp.nstr(c2[0], 4)}, β = {mp.nstr(c2[1], 5)}, "
                  f"resid {mp.nstr(r2, 3)};  c = 4π√(x/N): flattest n = {f2} (2 vs 12, metric {mp.nstr(met[f2], 3)}), {f4} (4 vs 12); "
                  f"R = {rec['prolate_c_4pi_sqrt']['R']}")
    labels = args.curves.split(",")
    if len(labels) == 2:
        for par in ("even", "odd"):
            a, b = lam[(labels[0], par)], lam[(labels[1], par)]
            d = {k: mp.nstr(mp.log(a[k] / b[k]), 4) for k in a}
            out[f"collapse_{par}_ln_{labels[0]}_over_{labels[1]}"] = d
            print(f"collapse {par}: ln(λ_{labels[0]}/λ_{labels[1]}) at equal x/N = {d}")
    if args.json:
        with open(args.json, "w") as fh:
            json.dump(out, fh, indent=1)


# ---------------------------------------------------------------------------------------------------------------
# P-DL7: level-1 eigenforms f_k, k ∈ {16, 18, 20, 22, 26}
# ---------------------------------------------------------------------------------------------------------------

def eisenstein_direct(m, nmax):
    """E_m = 1 − (2m/B_m) Σ σ_{m−1}(n) q^n from its Bernoulli coefficient (exact), independent of the E_4, E_6 products."""
    from fractions import Fraction
    p_, q_ = mp.bernfrac(m)
    c = Fraction(-2 * m) / Fraction(int(p_), int(q_))
    assert c.denominator == 1
    sig = _sigma_table(m - 1, nmax)
    return int(c), [1] + [int(c) * sig[n] for n in range(1, nmax + 1)]


def delta_times_direct(k, nmax):
    """Δ·E_{k−12} with E_{k−12} from eisenstein_direct (M_{k−12}(1) is one-dimensional, so this must equal eigen_table)."""
    from flint import fmpz_poly
    c1, E = eisenstein_direct(k - 12, nmax)
    tau = tau_table(nmax + 1)
    P = fmpz_poly([0] + tau[1:nmax + 1]) * fmpz_poly(E)
    co = [int(v) for v in P.coeffs()] + [0] * (nmax + 1)
    return c1, co[:nmax + 1]


def bernoulli_numerator_primes(k):
    """Numerator of B_k/k (= −ζ(1−k) up to sign) and its prime factors, computed exactly (mpmath bernfrac)."""
    from fractions import Fraction
    p_, q_ = mp.bernfrac(k)
    val = Fraction(int(p_), int(q_)) / k
    num = abs(val.numerator)
    return num, sorted(_factor(num))


def cmd_objects7(args):
    """P-DL7 form verification: a_1 = 1; the q-expansion two ways (E_4, E_6 products vs E_{k−12} from its Bernoulli
    coefficient −2m/B_m), a_2 = −24 + c_1(E_{k−12}); Hecke multiplicativity (rebuild a_n from a_p with p^{k−1});
    Deligne's bound; the Eisenstein congruence a_n ≡ σ_{k−1}(n) mod every prime dividing the numerator of B_k/k
    (Ramanujan's 691 for k = 12); Euler terms vs the log-derivative; root number i^k."""
    mp.mp.dps = 30
    job = Job("gl2 objects7", args=vars(args))
    report = {}
    for label in args.objects.split(","):
        k = EIGEN[label]
        a = eigen_table(k, args.nmax)
        hk = hecke_table(label, args.nmax)
        hecke_bad = sum(1 for n in range(1, args.nmax + 1) if hk[n] != a[n])
        ap = ap_table(label, args.nmax)
        deligne = all(ap[p] ** 2 <= 4 * p ** (k - 1) for p in ap)
        num, primes = bernoulli_numerator_primes(k)
        sig = _sigma_table(k - 1, args.cong)
        cong_bad = {q: sum(1 for n in range(1, args.cong + 1) if (a[n] - sig[n]) % q) for q in primes}
        dpp, doth = _logderiv_check(label, args.lognmax)
        eps = eps_of(label)
        c1, direct = delta_times_direct(k, args.nmax)
        direct_bad = sum(1 for n in range(1, args.nmax + 1) if direct[n] != a[n])
        a2_hand = -24 + c1                      # coefficient of q² in (q − 24q² + …)(1 + c1 q + …)
        job.result(f"{label} (k = {k}, ε = i^k = {eps:+d}): a_1..a_5 = {a[1:6]};  a_2 = {a[2]} (= −24 + c_1(E_{k - 12}) = −24 + {c1} = {a2_hand}: "
                   f"{a[2] == a2_hand});  Δ·E_{k - 12} with E_{k - 12} from −2m/B_m vs the E_4, E_6 product, n ≤ {args.nmax}: {direct_bad} mismatches;  "
                   f"Hecke (from a_p, p^{k - 1}) vs q-expansion, n ≤ {args.nmax}: {hecke_bad} mismatches;  Deligne, p ≤ {args.nmax}: {deligne};  "
                   f"numerator of B_{k}/{k} = {num} = {' · '.join(map(str, primes))}: a_n ≢ σ_{k - 1}(n) for n ≤ {args.cong}: {cong_bad};  "
                   f"Euler terms vs log-derivative (n ≤ {args.lognmax}): {mp.nstr(dpp, 3)}, off prime powers {mp.nstr(doth, 3)}")
        report[label] = {"k": k, "eps": eps, "a_1_10": a[1:11], "a2_hand": a2_hand, "a2_ok": a[2] == a2_hand,
                         "eisenstein_c1": c1, "direct_eisenstein_mismatches": direct_bad,
                         "hecke_mismatches": hecke_bad, "hecke_nmax": args.nmax, "deligne_ok": deligne,
                         "bernoulli_numerator": num, "congruence_primes": primes, "congruence_failures": cong_bad, "cong_nmax": args.cong,
                         "logderiv_dev_primepowers": mp.nstr(dpp, 3), "logderiv_max_off_primepowers": mp.nstr(doth, 3),
                         "mus": [f"{k - 1}/2", f"{k + 1}/2"]}
    if args.json:
        with open(args.json, "w") as fh:
            json.dump(report, fh, indent=1)
    job.done()


def cmd_central7(args):
    """Non-gating check of the central zero (ε = −1). The five gate tests all have F(0) = 0 (no b_0 component), so
    they cannot see γ = 0. The even combination c = (1, 1/√2, 0, …) vanishes at the window edges (b_0 + b_1/√2 = 0
    at u = ±L/2) and has F(0) = √L ≠ 0: compare Q with Σ_{γ>0} 2F² + m0·F(0)² and with the sum without the central term."""
    job = Job("gl2 central7", args=vars(args))
    out = []
    for label in args.objects.split(","):
        zd = json.load(open(os.path.join(OUT, f"zeros_{label}.json")))
        assert all(c["ok"] for c in zd["count_checks"])
        m0 = zd.get("central_multiplicity") or 0
        mp.mp.dps = 30
        x = mp.mpf(args.x)
        L = mp.log(x)
        c = [mp.mpf(1), 1 / mp.sqrt(2), 0, 0, 0, 0, 0]
        Ff = lambda t: dl.F_per_even(c, L, t)
        M = zeros_side_E(label, x, 6, "even")
        cv = mp.matrix(c)
        q = (cv.T * M * cv)[0]
        gam = [mp.mpf(g) for g in zd["zeros"]]
        s_pos = mp.fsum(2 * Ff(g) ** 2 for g in gam)
        f0 = Ff(mp.mpf(0))
        tail = tail_estimate(label, Ff, mp.mpf(zd["tmax"]), L)
        with_c = s_pos + m0 * f0 ** 2
        rec = {"obj": label, "x": args.x, "eps": eps_of(label), "m0": m0, "Q": mp.nstr(q, 15), "F0_sq": mp.nstr(f0 ** 2, 10),
               "sum_pos": mp.nstr(s_pos, 15), "sum_with_central": mp.nstr(with_c, 15), "expected_tail": mp.nstr(tail, 4),
               "diff_with_central": mp.nstr(q - with_c, 4), "ratio_with_central": mp.nstr((q - with_c) / tail, 4),
               "diff_without_central": mp.nstr(q - s_pos, 6)}
        out.append(rec)
        job.result(f"{label} (ε = {eps_of(label):+d}, m0 = {m0}) x = {args.x}, c = (1, 1/√2, 0…): Q = {mp.nstr(q, 12)};  F(0)² = {mp.nstr(f0 ** 2, 8)};  "
                   f"Q − (Σ_γ>0 + m0F(0)²) = {mp.nstr(q - with_c, 3)} (tail {mp.nstr(tail, 3)}, ratio {mp.nstr((q - with_c) / tail, 3)});  "
                   f"without the central term: Q − Σ_γ>0 = {mp.nstr(q - s_pos, 6)}")
    if args.json:
        with open(args.json, "w") as fh:
            json.dump(out, fh, indent=1)
    job.done()


def cmd_gate9(args):
    """P-DL9 gate for a fresh object (docs/DECAY_LAW.md, 40f46ee: "the P-DL7 gate, with the central zero included"):
    at x = 13 with zeros to T ≥ 240, the five standard tests (check7-style JSON) and the central-zero test
    (central7-style JSON) must all have |Q − Σ| ≤ 2e-4 and (Q − Σ)/tail ∈ [0.8, 1.25]."""
    mp.mp.dps = 30
    label = args.obj
    chk = json.load(open(args.check))
    rows = chk["rows"] if isinstance(chk, dict) else chk
    tests = []
    for r in rows:
        if r["curve"] == label and r["x"] == "13":
            tests.append({"test": f"{r['parity']} c={r['c']}", "T": r["T"], "diff": mp.mpf(r["abs_diff"]), "tail": mp.mpf(r["expected_tail"])})
    for r in json.load(open(args.central)):
        if r["obj"] == label and r["x"] == "13":
            tests.append({"test": "even c=(1, 1/√2, 0…) with the central zero", "T": None, "diff": mp.mpf(r["diff_with_central"]),
                          "tail": mp.mpf(r["expected_tail"]), "F0_sq": r["F0_sq"], "diff_without_central": r["diff_without_central"]})
    T_ok = all(t["T"] is None or mp.mpf(t["T"]) >= 240 for t in tests)
    out = []
    ok = T_ok and len(tests) == 6
    for t in tests:
        ratio = t["diff"] / t["tail"]
        good = abs(t["diff"]) <= mp.mpf("2e-4") and mp.mpf("0.8") <= ratio <= mp.mpf("1.25")
        ok &= good
        out.append({**{k: (mp.nstr(v, 6) if isinstance(v, mp.mpf) else v) for k, v in t.items()}, "ratio": mp.nstr(ratio, 4), "pass": bool(good)})
        print(f"{label} {t['test']}: |Q − Σ| = {mp.nstr(abs(t['diff']), 3)}, (Q − Σ)/tail = {mp.nstr(ratio, 3)} → {'pass' if good else 'FAIL'}")
    print(f"P-DL9 gate {label}: {len(tests)} tests, T ≥ 240: {T_ok} → {'PASS' if ok else 'FAIL'}")
    if args.json:
        with open(args.json, "w") as fh:
            json.dump({"obj": label, "tests": out, "pass": bool(ok)}, fh, indent=1)


def cmd_fit9(args):
    """P-DL9 as registered (docs/DECAY_LAW.md, 40f46ee), for any list of objects with scan9_<obj>.json.
    9a: ln λ = −αv + γ ln v + β (9k values); prediction α/4π ∈ [1.7, 2.3] for every object and sector; kill outside
        [1.4, 2.6]; the linear slope is reported alongside.
    9b: sign(ln(λ_odd/λ_even)) = ε at every grid point, for every object; kill on any mismatch at a point where both
        λ < 1e-6."""
    mp.mp.dps = 30
    four_pi = 4 * mp.pi
    out = {"9a": [], "9a_5k": [], "9b": []}
    for label in args.objects.split(","):
        rows = json.load(open(os.path.join(args.dir, f"scan9_{label}.json")))
        e = eps_of(label)
        for par in ("even", "odd"):
            R = sorted((r for r in rows if r["parity"] == par), key=lambda r: float(r["v"]))
            v = [mp.mpf(r["v"]) for r in R]
            for key, dest in (("lambda_n2", "9a"), ("lambda_n1", "9a_5k")):
                y = [mp.log(mp.mpf(r[key])) for r in R]
                c, res = dl._lsq([[-vi, mp.log(vi), 1] for vi in v], y)
                cl_, rl = dl._lsq([[-vi, 1] for vi in v], y)
                a4 = c[0] / four_pi
                status = ("in [1.7, 2.3]" if mp.mpf("1.7") <= a4 <= mp.mpf("2.3") else
                          "in kill-free band" if mp.mpf("1.4") <= a4 <= mp.mpf("2.6") else "KILL")
                rec = {"obj": label, "eps": e, "parity": par, "basis": "9k" if key == "lambda_n2" else "5k", "v": [r["v"] for r in R],
                       "alpha_over_4pi": mp.nstr(a4, 6), "gamma": mp.nstr(c[1], 6), "beta": mp.nstr(c[2], 6), "max_resid": mp.nstr(res, 4),
                       "linear_alpha_over_4pi": mp.nstr(cl_[0] / four_pi, 6), "linear_max_resid": mp.nstr(rl, 4), "status": status,
                       "conv_ln_l1_l2": [r["ln_l1_over_l2"] for r in R]}
                out[dest].append(rec)
                print(f"9a [{rec['basis']}] {label:<5} {par:<4}: α/4π = {mp.nstr(a4, 5):>7}  γ = {mp.nstr(c[1], 4):>7}  β = {mp.nstr(c[2], 5):>8}  "
                      f"max resid {mp.nstr(res, 3):>6};  linear α/4π = {mp.nstr(cl_[0] / four_pi, 5)}  → {status}")
        ev = {r["v"]: mp.mpf(r["lambda_n2"]) for r in rows if r["parity"] == "even"}
        od = {r["v"]: mp.mpf(r["lambda_n2"]) for r in rows if r["parity"] == "odd"}
        for vv in sorted(ev, key=float):
            lr = mp.log(od[vv] / ev[vv])
            sign = 1 if lr > 0 else -1
            small = ev[vv] < mp.mpf("1e-6") and od[vv] < mp.mpf("1e-6")
            rec = {"obj": label, "eps": e, "v": vv, "ln_odd_over_even": mp.nstr(lr, 6), "sign": sign, "match": sign == e,
                   "both_below_1e-6": bool(small), "lambda_even": mp.nstr(ev[vv], 6), "lambda_odd": mp.nstr(od[vv], 6)}
            out["9b"].append(rec)
            print(f"9b {label:<5} ε = {e:+d} v = {vv:<4}: ln(λ_odd/λ_even) = {mp.nstr(lr, 4):>7}  sign {sign:+d}  "
                  f"{'match' if sign == e else 'MISMATCH'}{'' if small else '  (a λ ≥ 1e-6: cannot kill)'}")
    a4s = [mp.mpf(r["alpha_over_4pi"]) for r in out["9a"]]
    out["verdict_9a"] = ("KILLED" if any(not (mp.mpf("1.4") <= t <= mp.mpf("2.6")) for t in a4s) else
                         "HOLDS" if all(mp.mpf("1.7") <= t <= mp.mpf("2.3") for t in a4s) else "NOT KILLED, prediction partly missed")
    mism = [r for r in out["9b"] if not r["match"]]
    out["verdict_9b"] = ("KILLED" if any(r["both_below_1e-6"] for r in mism) else
                         "HOLDS" if not mism else "NOT KILLED, mismatches only where a λ ≥ 1e-6")
    print(f"P-DL9a ({args.objects}): {out['verdict_9a']} (α/4π {mp.nstr(min(a4s), 4)}–{mp.nstr(max(a4s), 4)});  "
          f"P-DL9b: {out['verdict_9b']} ({len(out['9b']) - len(mism)}/{len(out['9b'])} grid points match)")
    if args.json:
        with open(args.json, "w") as fh:
            json.dump(out, fh, indent=1)


def cmd_objects11(args):
    """P-DL11 object verification (eta-product newforms g4, g6, g8), and the root number from the AFE with BOTH signs:
    q-expansion two ways (FLINT pentagonal route vs the naive product, n ≤ naive); a_1 = 1; a_mn = a_m a_n for coprime
    m, n with mn ≤ nmax; a_{p²} = a_p² − p^{k−1} at good p; a_{p^r} = a_p^r at p | N; the full Hecke rebuild from a_p;
    Deligne |a_p| ≤ 2p^{(k−1)/2} at good p and |a_N| = N^{(k−2)/2} at the bad prime; Euler terms vs the log-derivative;
    the Atkin–Lehner prediction ε = i^k w_N with w_N = −a_N/N^{k/2−1}; then, for ε = +1 and ε = −1 separately, the AFE
    split test (t = 1, 1.3, 0.75) and Λ(k − w) = εΛ(w) at random w. The sign that passes both is recorded."""
    import math
    mp.mp.dps = 30
    job = Job("gl2 objects11", args=vars(args))
    rng = random.Random(args.seed)
    report = {}
    for label in args.objects.split(","):
        f = ETA_NEWFORMS[label]
        k, N = f["k"], f["N"]
        a = eta_product_fast(args.nmax + 1, f["factors"])
        naive = eta_product(args.naive + 1, f["factors"])
        two_ways = sum(1 for n in range(1, args.naive + 1) if naive[n] != a[n])
        mult_bad = mult_n = 0
        for m in range(2, args.nmax + 1):
            for n in range(m + 1, args.nmax // m + 1):
                if math.gcd(m, n) == 1:
                    mult_n += 1
                    mult_bad += a[m * n] != a[m] * a[n]
        good = [p_ for p_ in primes_upto(int(args.nmax ** 0.5)) if N % p_]
        sq_bad = sum(1 for p_ in good if a[p_ * p_] != a[p_] ** 2 - p_ ** (k - 1))
        bad_pow = []
        for p_ in _factor(N):
            r, q_ = 1, p_
            while q_ <= args.nmax:
                bad_pow.append(a[q_] == a[p_] ** r)
                r, q_ = r + 1, q_ * p_
        hk = hecke_table(label, args.nmax)
        hecke_bad = sum(1 for n in range(1, args.nmax + 1) if hk[n] != a[n])
        ap = {p_: a[p_] for p_ in primes_upto(args.nmax)}
        deligne = all(ap[p_] ** 2 <= 4 * p_ ** (k - 1) for p_ in ap if N % p_)
        bad_abs = all(abs(ap[p_]) == p_ ** ((k - 2) // 2) for p_ in _factor(N))
        wN = -a[N] // N ** (k // 2 - 1)
        eps_pred = (1 if k % 4 == 0 else -1) * wN
        dpp, doth = _logderiv_check_eta(label, args.lognmax)
        job.result(f"{label} (k = {k}, N = {N}): a_1..a_8 = {a[1:9]};  fast vs naive product, n ≤ {args.naive}: {two_ways} mismatches;  "
                   f"a_mn = a_m a_n for {mult_n} coprime pairs, mn ≤ {args.nmax}: {mult_bad} failures;  a_p² = a_p² − p^{k - 1} at good p "
                   f"(p² ≤ {args.nmax}): {sq_bad} failures;  a_(N^r) = a_N^r: {all(bad_pow)} ({len(bad_pow)} powers);  Hecke rebuild vs "
                   f"q-expansion, n ≤ {args.nmax}: {hecke_bad} mismatches;  Deligne (good p ≤ {args.nmax}): {deligne};  |a_{N}| = {N}^{(k - 2) // 2}: "
                   f"{bad_abs} (a_{N} = {a[N]});  Euler terms vs log-derivative: {mp.nstr(dpp, 3)}, off prime powers {mp.nstr(doth, 3)};  "
                   f"Atkin–Lehner w_{N} = {wN:+d} → predicted ε = i^{k}·w_{N} = {eps_pred:+d}")
        # AFE with both signs
        sgn = {}
        h = mp.mpf(k) / 2
        pts = [mp.mpc(rng.uniform(float(h) - 2.5, float(h) + 2.5), rng.uniform(-25, 25)) for _ in range(args.npts)]
        for e in (1, -1):
            sp, fe = [], []
            for w in pts:
                with mp.workdps(_work_dps(mp.im(w), 30)):
                    v1 = afe_arith(label, w, t=1, eps=e)
                    v2 = afe_arith(label, w, t=mp.mpf("1.3"), eps=e)
                    v3 = afe_arith(label, w, t=mp.mpf("0.75"), eps=e)
                    vf = afe_arith(label, k - w, t=mp.mpf("1.3"), eps=e)
                    sp.append(max(abs(v2 - v1), abs(v3 - v1)) / abs(v1))
                    fe.append(abs(vf - e * v2) / abs(v2))
            sgn[e] = {"split_max": max(sp), "split_min": min(sp), "fe_max": max(fe), "fe_min": min(fe)}
            job.result(f"{label} ε = {e:+d}: split rel {mp.nstr(min(sp), 2)}–{mp.nstr(max(sp), 2)};  Λ(k − w) vs εΛ(w): rel "
                       f"{mp.nstr(min(fe), 2)}–{mp.nstr(max(fe), 2)}")
        passes = [e for e in (1, -1) if sgn[e]["split_max"] < mp.mpf("1e-30") and sgn[e]["fe_max"] < mp.mpf("1e-30")]
        fails = [e for e in (1, -1) if sgn[e]["split_min"] > mp.mpf("1e-6") or sgn[e]["fe_min"] > mp.mpf("1e-6")]
        eps = passes[0] if (len(passes) == 1 and len(fails) == 1) else None
        job.result(f"{label}: ε = {eps:+d} (the other sign fails); Atkin–Lehner prediction {eps_pred:+d}: {eps == eps_pred}" if eps
                   else f"{label}: ROOT NUMBER UNDETERMINED (passes {passes}, fails {fails})")
        # Dirichlet series and central value with the AFE's sign
        dirich = []
        if eps:
            mp.mp.dps = 30
            nds = 30000
            an = an_table(label, nds)
            for w in (mp.mpc(h + 4, 3), mp.mpc(h + mp.mpf("3.5"), 20)):
                with mp.workdps(_work_dps(mp.im(w), 30)):
                    ds = mp.fsum(an[n] * mp.power(n, -w) for n in range(1, nds + 1) if an[n])
                    lam = afe_arith(label, w, eps=eps)
                    lam_ds = (mp.sqrt(N) / (2 * mp.pi)) ** w * mp.gamma(w) * ds
                    dirich.append(mp.nstr(abs(lam - lam_ds) / abs(lam), 3))
            mp.mp.dps = 40
            nrm = lambda u: (mp.sqrt(N) / (2 * mp.pi)) ** u * mp.gamma(u)
            Lc = [afe_arith(label, h, t=tt, eps=eps) / nrm(h) for tt in (1, mp.mpf("1.3"))]
            cen = {"L_half": [mp.nstr(mp.re(v), 15) for v in Lc]}
            if eps == -1:
                dd = mp.mpf("1e-12")
                cen["Lprime_half"] = [mp.nstr(mp.re((afe_arith(label, h + dd, t=tt, eps=eps) / nrm(h + dd)
                                                     - afe_arith(label, h - dd, t=tt, eps=eps) / nrm(h - dd)) / (2 * dd)), 15)
                                      for tt in (1, mp.mpf("1.3"))]
            job.result(f"{label}: AFE vs Dirichlet series (Re w = k/2 + 4, k/2 + 3.5): {dirich};  central values {cen}")
        else:
            cen = {}
        report[label] = {"k": k, "N": N, "a_1_10": a[1:11], "fast_vs_naive_mismatches": two_ways, "naive_nmax": args.naive,
                         "coprime_pairs": mult_n, "multiplicativity_failures": mult_bad, "ap2_failures_good_p": sq_bad,
                         "bad_prime_powers_ok": all(bad_pow), "hecke_mismatches": hecke_bad, "hecke_nmax": args.nmax,
                         "deligne_ok": deligne, "bad_prime_abs_ok": bad_abs, "a_N": a[N], "atkin_lehner_w": wN,
                         "eps_predicted": eps_pred, "eps_afe": eps,
                         "afe_sign_test": {str(e): {kk: mp.nstr(vv, 3) for kk, vv in d.items()} for e, d in sgn.items()},
                         "dirichlet_rel": dirich, "central": cen,
                         "logderiv_dev_primepowers": mp.nstr(dpp, 3), "logderiv_max_off_primepowers": mp.nstr(doth, 3)}
    if args.json:
        with open(args.json, "w") as fh:
            json.dump(report, fh, indent=1)
    job.done()


def _logderiv_check_eta(label, nmax):
    """_logderiv_check without eps_of (the root number is not needed for the coefficients)."""
    return _logderiv_check(label, nmax)


DELTA11 = {(1, "even"): mp.mpf("0.7"), (1, "odd"): mp.mpf("3.1"), (-1, "even"): mp.mpf("3.8"), (-1, "odd"): mp.mpf("1.3")}


def continuous_index(v, lam, nmax=80, extrapolate=False):
    """The registered continuous index (docs/DECAY_LAW.md, 3778445): the n at which the least-squares slope of
    ln(λ/ℓ_n(c)) against ln c over the grid (c = 4πv) crosses zero, linearly interpolated between integers.
    Returns (n*, slopes, S − ½) where S is the slope of ln λ + 2c against ln c; since ln ℓ_n(c) = const_n +
    (n + ½) ln c − 2c, the slope is exactly S − n − ½ and n* = S − ½ (a cross-check, not the method)."""
    c = [4 * mp.pi * vi for vi in v]
    xs = [mp.log(ci) for ci in c]
    xm = mp.fsum(xs) / len(xs)

    def slope(ys):
        ym = mp.fsum(ys) / len(ys)
        return mp.fsum((x - xm) * (y - ym) for x, y in zip(xs, ys)) / mp.fsum((x - xm) ** 2 for x in xs)

    sl = [slope([mp.log(l / dl._ell(n, ci)) for l, ci in zip(lam, c)]) for n in range(nmax + 1)]
    nstar = None
    for n in range(nmax):
        if sl[n] == 0:
            nstar = mp.mpf(n)
            break
        if sl[n] > 0 > sl[n + 1]:
            nstar = n + sl[n] / (sl[n] - sl[n + 1])
            break
    if nstar is None and extrapolate:
        # the slope is exactly linear in n, so a crossing outside 0..nmax is extrapolated linearly (and flagged by callers)
        nstar = sl[0] / (sl[0] - sl[1]) if sl[0] < 0 else nmax + sl[nmax] / (sl[nmax - 1] - sl[nmax])
    S = slope([mp.log(l) + 2 * ci for l, ci in zip(lam, c)])
    return nstar, sl, S - mp.mpf(1) / 2


def cmd_fit11(args):
    """P-DL11 as registered (docs/DECAY_LAW.md, 3778445), on scan11_<obj>.json for g4, g6, g8, 11a1, 37a1.
    11a: ln λ = −αv + γ ln v + β (9k values); prediction α/4π ∈ [1.7, 2.3]; kill outside [1.4, 2.6]; linear slope alongside.
    11b: continuous index n* (continuous_index); n* − k within ±1.2 of δ(ε, s), δ = 0.7, 3.1, 3.8, 1.3 for (+1, even),
         (+1, odd), (−1, even), (−1, odd); holds if ≥ 8 of the 10 (object, sector) cases are within ±1.2; kill if more
         than 2 cases are off by more than 2.
    11c: sign ln(λ_odd/λ_even) = ε at every point; kill on any mismatch where both λ < 1e-6."""
    mp.mp.dps = 30
    four_pi = 4 * mp.pi
    out = {"11a": [], "11a_5k": [], "11b": [], "11c": []}
    for label in args.objects.split(","):
        rows = json.load(open(os.path.join(args.dir, f"scan11_{label}.json")))
        e, k = eps_of(label), weight(label)
        for par in ("even", "odd"):
            R = sorted((r for r in rows if r["parity"] == par), key=lambda r: float(r["v"]))
            v = [mp.mpf(r["v"]) for r in R]
            for key, dest in (("lambda_n2", "11a"), ("lambda_n1", "11a_5k")):
                y = [mp.log(mp.mpf(r[key])) for r in R]
                c, res = dl._lsq([[-vi, mp.log(vi), 1] for vi in v], y)
                cl_, rl = dl._lsq([[-vi, 1] for vi in v], y)
                a4 = c[0] / four_pi
                status = ("in [1.7, 2.3]" if mp.mpf("1.7") <= a4 <= mp.mpf("2.3") else
                          "in kill-free band" if mp.mpf("1.4") <= a4 <= mp.mpf("2.6") else "KILL")
                out[dest].append({"obj": label, "k": k, "eps": e, "parity": par, "basis": "9k" if key == "lambda_n2" else "5k",
                                  "v": [r["v"] for r in R], "alpha_over_4pi": mp.nstr(a4, 6), "gamma": mp.nstr(c[1], 6),
                                  "beta": mp.nstr(c[2], 6), "max_resid": mp.nstr(res, 4), "linear_alpha_over_4pi": mp.nstr(cl_[0] / four_pi, 6),
                                  "status": status, "conv_ln_l1_l2": [r.get("ln_l1_over_l2") for r in R]})
                print(f"11a [{'9k' if key == 'lambda_n2' else '5k'}] {label:<5} {par:<4}: α/4π = {mp.nstr(a4, 5):>7}  γ = {mp.nstr(c[1], 4):>7}  "
                      f"max resid {mp.nstr(res, 3):>6};  linear α/4π = {mp.nstr(cl_[0] / four_pi, 5)}  → {status}")
            lam = [mp.mpf(r["lambda_n2"]) for r in R]
            nstar, sl, check = continuous_index(v, lam)
            d = nstar - k if nstar is not None else None
            delta = DELTA11[(e, par)]
            off = abs(d - delta) if d is not None else mp.inf
            out["11b"].append({"obj": label, "k": k, "eps": e, "parity": par, "n_star": mp.nstr(nstar, 6) if nstar is not None else None,
                               "S_minus_half": mp.nstr(check, 6), "n_minus_k": mp.nstr(d, 5) if d is not None else None, "delta": mp.nstr(delta, 3),
                               "off": mp.nstr(off, 4), "within_1.2": bool(off <= mp.mpf("1.2")), "off_by_more_than_2": bool(off > 2)})
            print(f"11b {label:<5} k = {k:>2} ε = {e:+d} {par:<4}: n* = {mp.nstr(nstar, 5)} (S − ½ = {mp.nstr(check, 5)}), n* − k = {mp.nstr(d, 4)}, "
                  f"δ = {mp.nstr(delta, 2)}, |off| = {mp.nstr(off, 3)}  → {'within ±1.2' if off <= mp.mpf('1.2') else 'off > 2' if off > 2 else 'off (1.2, 2]'}")
        ev = {r["v"]: mp.mpf(r["lambda_n2"]) for r in rows if r["parity"] == "even"}
        od = {r["v"]: mp.mpf(r["lambda_n2"]) for r in rows if r["parity"] == "odd"}
        for vv in sorted(ev, key=float):
            lr = mp.log(od[vv] / ev[vv])
            sign = 1 if lr > 0 else -1
            small = ev[vv] < mp.mpf("1e-6") and od[vv] < mp.mpf("1e-6")
            out["11c"].append({"obj": label, "eps": e, "v": vv, "ln_odd_over_even": mp.nstr(lr, 6), "sign": sign, "match": sign == e,
                               "both_below_1e-6": bool(small)})
            print(f"11c {label:<5} ε = {e:+d} v = {vv:<5}: ln(λ_odd/λ_even) = {mp.nstr(lr, 4):>7}  {'match' if sign == e else 'MISMATCH'}"
                  f"{'' if small else '  (a λ ≥ 1e-6: cannot kill)'}")
    a4s = [mp.mpf(r["alpha_over_4pi"]) for r in out["11a"]]
    out["verdict_11a"] = ("KILLED" if any(not (mp.mpf("1.4") <= t <= mp.mpf("2.6")) for t in a4s) else
                          "HOLDS" if all(mp.mpf("1.7") <= t <= mp.mpf("2.3") for t in a4s) else "NOT KILLED, prediction partly missed")
    within = sum(1 for r in out["11b"] if r["within_1.2"])
    far = sum(1 for r in out["11b"] if r["off_by_more_than_2"])
    out["verdict_11b"] = "KILLED" if far > 2 else ("HOLDS" if within >= 8 else "NOT KILLED, prediction partly missed")
    out["11b_within"], out["11b_off_by_more_than_2"], out["11b_cases"] = within, far, len(out["11b"])
    mism = [r for r in out["11c"] if not r["match"]]
    out["verdict_11c"] = ("KILLED" if any(r["both_below_1e-6"] for r in mism) else
                          "HOLDS" if not mism else "NOT KILLED, mismatches only where a λ ≥ 1e-6")
    print(f"P-DL11a: {out['verdict_11a']} (α/4π {mp.nstr(min(a4s), 4)}–{mp.nstr(max(a4s), 4)});  P-DL11b: {out['verdict_11b']} "
          f"({within}/{len(out['11b'])} within ±1.2, {far} off by > 2);  P-DL11c: {out['verdict_11c']} "
          f"({len(out['11c']) - len(mism)}/{len(out['11c'])} points match)")
    if args.json:
        with open(args.json, "w") as fh:
            json.dump(out, fh, indent=1)


def afe_sign_and_rank(label, rng, npts, job):
    """Root number from the AFE with BOTH signs (split test t = 1, 1.3, 0.75 and Λ(k − w) = εΛ(w) at random w; the
    passing sign must pass both to < 1e-30 and the other must fail > 1e-6), then the Dirichlet-series check and the
    central value with that sign: rank 0 if ε = +1 and L(½) ≠ 0; rank 1 if ε = −1 (L(½) = 0 forced) and L'(½) ≠ 0."""
    k, N = weight(label), N_of(label)
    h = mp.mpf(k) / 2
    mp.mp.dps = 30
    pts = [mp.mpc(rng.uniform(float(h) - 2.5, float(h) + 2.5), rng.uniform(-25, 25)) for _ in range(npts)]
    sgn = {}
    for e in (1, -1):
        sp, fe = [], []
        for w in pts:
            with mp.workdps(_work_dps(mp.im(w), 30)):
                v1 = afe_arith(label, w, t=1, eps=e)
                v2 = afe_arith(label, w, t=mp.mpf("1.3"), eps=e)
                v3 = afe_arith(label, w, t=mp.mpf("0.75"), eps=e)
                vf = afe_arith(label, k - w, t=mp.mpf("1.3"), eps=e)
                sp.append(max(abs(v2 - v1), abs(v3 - v1)) / abs(v1))
                fe.append(abs(vf - e * v2) / abs(v2))
        sgn[e] = {"split_max": max(sp), "split_min": min(sp), "fe_max": max(fe), "fe_min": min(fe)}
        job.result(f"{label} ε = {e:+d}: split rel {mp.nstr(min(sp), 2)}–{mp.nstr(max(sp), 2)};  Λ(k − w) vs εΛ(w): rel "
                   f"{mp.nstr(min(fe), 2)}–{mp.nstr(max(fe), 2)}")
    passes = [e for e in (1, -1) if sgn[e]["split_max"] < mp.mpf("1e-30") and sgn[e]["fe_max"] < mp.mpf("1e-30")]
    fails = [e for e in (1, -1) if sgn[e]["split_min"] > mp.mpf("1e-6") or sgn[e]["fe_min"] > mp.mpf("1e-6")]
    eps = passes[0] if (len(passes) == 1 and len(fails) == 1) else None
    out = {"afe_sign_test": {str(e): {kk: mp.nstr(vv, 3) for kk, vv in d.items()} for e, d in sgn.items()}, "eps_afe": eps}
    if eps is None:
        job.result(f"{label}: ROOT NUMBER UNDETERMINED (passes {passes}, fails {fails})")
        return out
    nds = 30000
    an = an_table(label, nds)
    dirich = []
    for w in (mp.mpc(h + 4, 3), mp.mpc(h + mp.mpf("3.5"), 20)):
        with mp.workdps(_work_dps(mp.im(w), 30)):
            ds = mp.fsum(an[n] * mp.power(n, -w) for n in range(1, nds + 1) if an[n])
            lam = afe_arith(label, w, eps=eps)
            dirich.append(mp.nstr(abs(lam - (mp.sqrt(N) / (2 * mp.pi)) ** w * mp.gamma(w) * ds) / abs(lam), 3))
    mp.mp.dps = 40
    nrm = lambda u: (mp.sqrt(N) / (2 * mp.pi)) ** u * mp.gamma(u)
    Lc = [mp.re(afe_arith(label, h, t=tt, eps=eps) / nrm(h)) for tt in (1, mp.mpf("1.3"))]
    dd = mp.mpf("1e-12")
    Ld = [mp.re((afe_arith(label, h + dd, t=tt, eps=eps) / nrm(h + dd) - afe_arith(label, h - dd, t=tt, eps=eps) / nrm(h - dd)) / (2 * dd))
          for tt in (1, mp.mpf("1.3"))]
    d2 = mp.mpf("1e-6")
    L2 = [mp.re((afe_arith(label, h + d2, t=tt, eps=eps) / nrm(h + d2) - 2 * afe_arith(label, h, t=tt, eps=eps) / nrm(h)
                 + afe_arith(label, h - d2, t=tt, eps=eps) / nrm(h - d2)) / d2 ** 2) for tt in (1, mp.mpf("1.3"))]
    out["Lsecond_half"] = [mp.nstr(v, 12) for v in L2]
    if eps == 1:
        if all(abs(v) > mp.mpf("1e-10") for v in Lc) and abs(Lc[0] - Lc[1]) < mp.mpf("1e-25"):
            rank = 0
        elif (all(abs(v) < mp.mpf("1e-30") for v in Lc) and all(abs(v) < mp.mpf("1e-12") for v in Ld)
              and all(abs(v) > mp.mpf("1e-6") for v in L2) and abs(L2[0] - L2[1]) < mp.mpf("1e-6") * abs(L2[0])):
            rank = 2
        else:
            rank = None
    else:
        rank = 1 if (all(abs(v) < mp.mpf("1e-30") for v in Lc) and all(abs(v) > mp.mpf("1e-10") for v in Ld)
                     and abs(Ld[0] - Ld[1]) < mp.mpf("1e-12")) else None
    out.update({"dirichlet_rel": dirich, "L_half": [mp.nstr(v, 15) for v in Lc], "Lprime_half": [mp.nstr(v, 15) for v in Ld],
                "analytic_rank": rank})
    job.result(f"{label}: ε = {eps:+d} (the other sign fails);  AFE vs Dirichlet series: {dirich};  L(½) = {[mp.nstr(v, 12) for v in Lc]};  "
               f"L'(½) = {[mp.nstr(v, 12) for v in Ld]}  → analytic rank {rank}")
    return out


def cmd_afe_both(args):
    """Root number by the AFE with both signs, and the analytic rank (afe_sign_and_rank), for any objects."""
    job = Job("gl2 afe-both", args=vars(args))
    rng = random.Random(args.seed)
    report = {label: afe_sign_and_rank(label, rng, args.npts, job) for label in args.objects.split(",")}
    if args.json:
        with open(args.json, "w") as fh:
            json.dump(report, fh, indent=1)
    job.done()


HYP12 = {(1, "even"): mp.mpf("-0.64"), (1, "odd"): mp.mpf("1.81"), (-1, "even"): mp.mpf("1.32"), (-1, "odd"): mp.mpf("-0.67")}


def cmd_fit12(args):
    """P-DL12 as registered (docs/DECAY_LAW.md, b22392c), on scan12_<obj>.json for 14a1, 17a1, 43a1, 53a1; the P-DL11
    metric (continuous_index, c = 4πv).
    12a: n* − k within ±0.5 of the hypothesis (ε = +1: −0.64 even, +1.81 odd; ε = −1: +1.32 even, −0.67 odd) in at least
         6 of the 8 cases; kill if more than 2 cases are off by more than 1.0.
    12b: log-corrected α/4π ∈ [1.7, 2.3]; kill outside [1.4, 2.6]; linear slope alongside.
    12c: sign ln(λ_odd/λ_even) = ε at every point; kill on any mismatch with both λ < 1e-6."""
    mp.mp.dps = 30
    four_pi = 4 * mp.pi
    out = {"12a": [], "12b": [], "12b_5k": [], "12c": []}
    for label in args.objects.split(","):
        rows = json.load(open(os.path.join(args.dir, f"scan12_{label}.json")))
        e, k = eps_of(label), weight(label)
        for par in ("even", "odd"):
            R = sorted((r for r in rows if r["parity"] == par), key=lambda r: float(r["v"]))
            v = [mp.mpf(r["v"]) for r in R]
            for key, dest in (("lambda_n2", "12b"), ("lambda_n1", "12b_5k")):
                y = [mp.log(mp.mpf(r[key])) for r in R]
                c, res = dl._lsq([[-vi, mp.log(vi), 1] for vi in v], y)
                cl_, rl = dl._lsq([[-vi, 1] for vi in v], y)
                a4 = c[0] / four_pi
                status = ("in [1.7, 2.3]" if mp.mpf("1.7") <= a4 <= mp.mpf("2.3") else
                          "in kill-free band" if mp.mpf("1.4") <= a4 <= mp.mpf("2.6") else "KILL")
                out[dest].append({"obj": label, "k": k, "eps": e, "parity": par, "basis": "9k" if key == "lambda_n2" else "5k",
                                  "v": [r["v"] for r in R], "alpha_over_4pi": mp.nstr(a4, 6), "gamma": mp.nstr(c[1], 6),
                                  "beta": mp.nstr(c[2], 6), "max_resid": mp.nstr(res, 4), "linear_alpha_over_4pi": mp.nstr(cl_[0] / four_pi, 6),
                                  "status": status, "conv_ln_l1_l2": [r.get("ln_l1_over_l2") for r in R]})
                print(f"12b [{'9k' if key == 'lambda_n2' else '5k'}] {label:<5} {par:<4}: α/4π = {mp.nstr(a4, 5):>7}  γ = {mp.nstr(c[1], 4):>7}  "
                      f"max resid {mp.nstr(res, 3):>6};  linear α/4π = {mp.nstr(cl_[0] / four_pi, 5)}  → {status}")
            lam = [mp.mpf(r["lambda_n2"]) for r in R]
            nstar, sl, check = continuous_index(v, lam, extrapolate=True)
            extrap = not (0 <= nstar <= 80)
            d = nstar - k if nstar is not None else None
            hyp = HYP12[(e, par)]
            off = abs(d - hyp) if d is not None else mp.inf
            out["12a"].append({"obj": label, "k": k, "eps": e, "parity": par, "n_star": mp.nstr(nstar, 6) if nstar is not None else None,
                               "S_minus_half": mp.nstr(check, 6), "n_minus_k": mp.nstr(d, 5) if d is not None else None,
                               "hypothesis": mp.nstr(hyp, 3), "off": mp.nstr(off, 4), "within_0.5": bool(off <= mp.mpf("0.5")),
                               "off_by_more_than_1": bool(off > 1), "extrapolated_outside_0_80": bool(extrap)})
            print(f"12a {label:<5} k = {k} ε = {e:+d} {par:<4}: n* = {mp.nstr(nstar, 5)} (S − ½ = {mp.nstr(check, 5)}), n* − k = {mp.nstr(d, 4)}, "
                  f"hypothesis {mp.nstr(hyp, 3)}, |off| = {mp.nstr(off, 3)}  → "
                  f"{'within ±0.5' if off <= mp.mpf('0.5') else 'off > 1.0' if off > 1 else 'off (0.5, 1.0]'}"
                  f"{'  (crossing extrapolated outside n = 0..80)' if extrap else ''}")
        ev = {r["v"]: mp.mpf(r["lambda_n2"]) for r in rows if r["parity"] == "even"}
        od = {r["v"]: mp.mpf(r["lambda_n2"]) for r in rows if r["parity"] == "odd"}
        for vv in sorted(ev, key=float):
            lr = mp.log(od[vv] / ev[vv])
            sign = 1 if lr > 0 else -1
            small = ev[vv] < mp.mpf("1e-6") and od[vv] < mp.mpf("1e-6")
            out["12c"].append({"obj": label, "eps": e, "v": vv, "ln_odd_over_even": mp.nstr(lr, 6), "sign": sign, "match": sign == e,
                               "both_below_1e-6": bool(small)})
            print(f"12c {label:<5} ε = {e:+d} v = {vv:<5}: ln(λ_odd/λ_even) = {mp.nstr(lr, 4):>7}  {'match' if sign == e else 'MISMATCH'}"
                  f"{'' if small else '  (a λ ≥ 1e-6: cannot kill)'}")
    within = sum(1 for r in out["12a"] if r["within_0.5"])
    far = sum(1 for r in out["12a"] if r["off_by_more_than_1"])
    out["verdict_12a"] = "KILLED" if far > 2 else ("HOLDS" if within >= 6 else "NOT KILLED, prediction partly missed")
    out["12a_within"], out["12a_off_by_more_than_1"], out["12a_cases"] = within, far, len(out["12a"])
    a4s = [mp.mpf(r["alpha_over_4pi"]) for r in out["12b"]]
    out["verdict_12b"] = ("KILLED" if any(not (mp.mpf("1.4") <= t <= mp.mpf("2.6")) for t in a4s) else
                          "HOLDS" if all(mp.mpf("1.7") <= t <= mp.mpf("2.3") for t in a4s) else "NOT KILLED, prediction partly missed")
    mism = [r for r in out["12c"] if not r["match"]]
    out["verdict_12c"] = ("KILLED" if any(r["both_below_1e-6"] for r in mism) else
                          "HOLDS" if not mism else "NOT KILLED, mismatches only where a λ ≥ 1e-6")
    print(f"P-DL12a: {out['verdict_12a']} ({within}/{len(out['12a'])} within ±0.5, {far} off by > 1.0);  P-DL12b: {out['verdict_12b']} "
          f"(α/4π {mp.nstr(min(a4s), 4)}–{mp.nstr(max(a4s), 4)});  P-DL12c: {out['verdict_12c']} ({len(out['12c']) - len(mism)}/{len(out['12c'])} match)")
    if args.json:
        with open(args.json, "w") as fh:
            json.dump(out, fh, indent=1)


def cmd_fit7(args):
    """P-DL7 as registered (docs/DECAY_LAW.md, 6fd392d).
    7a: ln λ = −αv + γ ln v + β (9k values); prediction α/4π ∈ [1.7, 2.3] for every form and sector; kill outside [1.4, 2.6].
    Index: c = 4πv, n ∈ 0..40, flattest = argmin_n |ln R(v₅)/R(v₂)|.
    7b: n_even(k) = a·k + b over k ∈ {12, 16, 18, 20, 22, 26}, Δ's k = 12 value from pdl6_fit.json; prediction
        a ∈ [0.8, 1.2]; kill outside [0.6, 1.5].
    7c: n_odd − n_even ∈ {1, 2, 3} for every k; kill if ≤ 0 or ≥ 5 for any k."""
    mp.mp.dps = 30
    four_pi = 4 * mp.pi
    out = {"7a": [], "7a_5k": [], "index": []}
    flat = {}
    for label in args.objects.split(","):
        k = EIGEN[label]
        rows = json.load(open(os.path.join(args.dir, f"scan7_{label}.json")))
        for par in ("even", "odd"):
            R = sorted((r for r in rows if r["parity"] == par), key=lambda r: float(r["v"]))
            v = [mp.mpf(r["v"]) for r in R]
            for key, dest in (("lambda_n2", "7a"), ("lambda_n1", "7a_5k")):
                y = [mp.log(mp.mpf(r[key])) for r in R]
                c, res = dl._lsq([[-vi, mp.log(vi), 1] for vi in v], y)
                cl_, rl = dl._lsq([[-vi, 1] for vi in v], y)
                a4 = c[0] / four_pi
                status = ("in [1.7, 2.3]" if mp.mpf("1.7") <= a4 <= mp.mpf("2.3") else
                          "in kill-free band" if mp.mpf("1.4") <= a4 <= mp.mpf("2.6") else "KILL")
                rec = {"obj": label, "k": k, "parity": par, "basis": "9k" if key == "lambda_n2" else "5k", "alpha_over_4pi": mp.nstr(a4, 6),
                       "gamma": mp.nstr(c[1], 6), "beta": mp.nstr(c[2], 6), "max_resid": mp.nstr(res, 4),
                       "linear_alpha_over_4pi": mp.nstr(cl_[0] / four_pi, 6), "status": status, "conv_ln_l1_l2": [r["ln_l1_over_l2"] for r in R]}
                out[dest].append(rec)
                print(f"7a [{rec['basis']}] {label:<4} {par:<4}: α/4π = {mp.nstr(a4, 5):>7}  γ = {mp.nstr(c[1], 4):>7}  β = {mp.nstr(c[2], 5):>8}  "
                      f"max resid {mp.nstr(res, 3):>6};  linear α/4π = {mp.nstr(cl_[0] / four_pi, 5)}  → {status}")
            lam = [mp.mpf(r["lambda_n2"]) for r in R]
            ratios = {m: [l / dl._ell(m, four_pi * vi) for l, vi in zip(lam, v)] for m in range(41)}
            metric = {m: abs(mp.log(r[4] / r[1])) for m, r in ratios.items()}
            fm = min(metric, key=lambda m: metric[m])
            second = sorted(metric, key=lambda m: metric[m])[1]
            flat[(k, par)] = fm
            out["index"].append({"obj": label, "k": k, "parity": par, "flattest_n": fm, "runner_up": second,
                                 "metric": {m: mp.nstr(metric[m], 4) for m in metric}, "R_flattest": [mp.nstr(t, 4) for t in ratios[fm]]})
            print(f"index {label:<4} {par:<4}: flattest n = {fm} (metric {mp.nstr(metric[fm], 3)}; runner-up {second}: {mp.nstr(metric[second], 3)})")
    p6 = json.load(open(os.path.join(args.dir, "pdl6_fit.json")))
    d6 = {r["parity"]: r["flattest_n"] for r in p6["6b"] if r["obj"] == DELTA}
    ks = [12] + [EIGEN[lb] for lb in args.objects.split(",")]
    ne = [d6["even"]] + [flat[(kk, "even")] for kk in ks[1:]]
    c, res = dl._lsq([[kk, 1] for kk in ks], [mp.mpf(t) for t in ne])
    a_ = c[0]
    vb = "KILLED" if not (mp.mpf("0.6") <= a_ <= mp.mpf("1.5")) else ("HOLDS" if mp.mpf("0.8") <= a_ <= mp.mpf("1.2") else "NOT KILLED, prediction missed")
    out["7b"] = {"k": ks, "n_even": ne, "a": mp.nstr(a_, 6), "b": mp.nstr(c[1], 6), "max_resid": mp.nstr(res, 4), "delta_from": "pdl6_fit.json", "verdict": vb}
    offs = {kk: flat[(kk, "odd")] - flat[(kk, "even")] for kk in ks[1:]}
    offs_all = dict(offs)
    offs_all[12] = d6["odd"] - d6["even"]
    vc = ("KILLED" if any(o <= 0 or o >= 5 for o in offs.values()) else
          "HOLDS" if all(o in (1, 2, 3) for o in offs.values()) else "NOT KILLED, prediction partly missed")
    out["7c"] = {"offsets": {str(kk): o for kk, o in sorted(offs.items())}, "delta_offset_reference": offs_all[12], "verdict": vc}
    a4s = [mp.mpf(r["alpha_over_4pi"]) for r in out["7a"]]
    va = ("KILLED" if any(not (mp.mpf("1.4") <= t <= mp.mpf("2.6")) for t in a4s) else
          "HOLDS" if all(mp.mpf("1.7") <= t <= mp.mpf("2.3") for t in a4s) else "NOT KILLED, prediction partly missed")
    out["verdict_7a"] = va
    print(f"P-DL7a: {va} (α/4π {mp.nstr(min(a4s), 4)}–{mp.nstr(max(a4s), 4)})")
    print(f"P-DL7b: n_even = {dict(zip(ks, ne))}; fit a = {mp.nstr(a_, 4)}, b = {mp.nstr(c[1], 4)}, max resid {mp.nstr(res, 3)} → {vb}")
    print(f"P-DL7c: offsets n_odd − n_even = {offs} (Δ, k = 12, for reference: {offs_all[12]}) → {vc}")
    if args.json:
        with open(args.json, "w") as fh:
            json.dump(out, fh, indent=1)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    v = sub.add_parser("validate")
    v.add_argument("--x", default="13,40")
    v.add_argument("--n", type=int, default=12)
    v.add_argument("--dps", type=int, default=50)
    v.add_argument("--json")
    c = sub.add_parser("curves")
    c.add_argument("--pmax", type=int, default=2000)
    c.add_argument("--nmax", type=int, default=2000)
    c.add_argument("--lognmax", type=int, default=500)
    c.add_argument("--json")
    a = sub.add_parser("afe")
    a.add_argument("--curves", default="11a1,37b1")
    a.add_argument("--npts", type=int, default=5)
    a.add_argument("--dps", type=int, default=30)
    a.add_argument("--seed", type=int, default=20261002)
    a.add_argument("--json")
    z = sub.add_parser("zeros")
    z.add_argument("--curve", required=True)
    z.add_argument("--tmax", type=float, default=120)
    z.add_argument("--step", default="0.05")
    z.add_argument("--count-at", default="20,40,60,80,100,120")
    z.add_argument("--json")
    z.add_argument("--tmin", type=float, default=0, help="start of the τ-range (split runs)")
    z.add_argument("--chunk", action="store_true", help="save the zeros of this τ-range only (merge with zeros-merge)")
    zm = sub.add_parser("zeros-merge")
    zm.add_argument("--curve", required=True)
    zm.add_argument("--parts", required=True, help="comma-separated chunk JSONs")
    zm.add_argument("--count-at", default="40,80,120,160,200,240")
    zm.add_argument("--json", required=True)
    k = sub.add_parser("check")
    k.add_argument("--curves", default="11a1,37b1")
    k.add_argument("--x", default="13")
    k.add_argument("--gate", action="store_true", help="apply the P-DL6 zero-sum gate at x = 13")
    k.add_argument("--gate7", action="store_true", help="apply the P-DL7 zero-sum gate at x = 13 (T ≥ 240)")
    k.add_argument("--zeros-tag", default="", help="read zeros_<label><tag>.json (e.g. _T240 for a diagnostic)")
    k.add_argument("--json")
    s = sub.add_parser("scan")
    s.add_argument("--curve", required=True)
    s.add_argument("--xq", default="2,4,6,8,12")
    s.add_argument("--f1", type=float, default=5)
    s.add_argument("--f2", type=float, default=9)
    s.add_argument("--nmin", type=int, default=24)
    s.add_argument("--acb-max", type=int, default=120, help="also run FLINT's full eigensolver (PSD check) for N_basis ≤ this")
    s.add_argument("--resume", action="store_true")
    s.add_argument("--json")
    ph = sub.add_parser("posthoc")
    ph.add_argument("--curves", default="11a1,37b1")
    ph.add_argument("--dir", default=OUT)
    ph.add_argument("--json")
    pr = sub.add_parser("precision")
    pr.add_argument("--curve", required=True)
    pr.add_argument("--extra", type=int, default=80)
    pr.add_argument("--file", help="scan JSON (default: P-DL4 scan_<curve>.json)")
    pr.add_argument("--out", help="output JSON (default: precision_<curve>.json)")
    o6 = sub.add_parser("objects6")
    o6.add_argument("--objects", default="15a1,19a1,delta")
    o6.add_argument("--pmax", type=int, default=2000)
    o6.add_argument("--nmax", type=int, default=2000)
    o6.add_argument("--naive", type=int, default=300)
    o6.add_argument("--cong", type=int, default=1000)
    o6.add_argument("--lognmax", type=int, default=500)
    o6.add_argument("--json")
    s6 = sub.add_parser("scan6")
    s6.add_argument("--obj", required=True)
    s6.add_argument("--v", default=None, help="default: 1.5,2,2.5,3,3.5 (curves), 2,2.5,3,3.5,4 (delta)")
    s6.add_argument("--f1", type=float, default=5)
    s6.add_argument("--f2", type=float, default=9)
    s6.add_argument("--nmin", type=int, default=24)
    s6.add_argument("--acb-max", type=int, default=120)
    s6.add_argument("--resume", action="store_true")
    s6.add_argument("--method", choices=["inv", "inv+full"], default="inv",
                    help="inv: inverse iteration only (P-DL6); inv+full: fall back to FLINT's full eigensolver when it does not converge")
    s6.add_argument("--json")
    o7 = sub.add_parser("objects7")
    o7.add_argument("--objects", default="f16,f18,f20,f22,f26")
    o7.add_argument("--nmax", type=int, default=2000)
    o7.add_argument("--cong", type=int, default=1000)
    o7.add_argument("--lognmax", type=int, default=500)
    o7.add_argument("--json")
    c7 = sub.add_parser("central7")
    c7.add_argument("--objects", default="f16,f18,f20,f22,f26")
    c7.add_argument("--x", default="13")
    c7.add_argument("--json")
    g9 = sub.add_parser("gate9")
    g9.add_argument("--obj", required=True)
    g9.add_argument("--check", required=True)
    g9.add_argument("--central", required=True)
    g9.add_argument("--json")
    f9 = sub.add_parser("fit9")
    f9.add_argument("--objects", default="delta,f16,f18,f20,f22,f26,37a1")
    f9.add_argument("--dir", default=OUT)
    f9.add_argument("--json")
    o11 = sub.add_parser("objects11")
    o11.add_argument("--objects", default="g4,g6,g8")
    o11.add_argument("--nmax", type=int, default=2000)
    o11.add_argument("--naive", type=int, default=300)
    o11.add_argument("--lognmax", type=int, default=500)
    o11.add_argument("--npts", type=int, default=5)
    o11.add_argument("--seed", type=int, default=20261002)
    o11.add_argument("--json")
    f11 = sub.add_parser("fit11")
    f11.add_argument("--objects", default="g4,g6,g8,11a1,37a1")
    f11.add_argument("--dir", default=OUT)
    f11.add_argument("--json")
    ab = sub.add_parser("afe-both")
    ab.add_argument("--objects", required=True)
    ab.add_argument("--npts", type=int, default=5)
    ab.add_argument("--seed", type=int, default=20261002)
    ab.add_argument("--json")
    f12 = sub.add_parser("fit12")
    f12.add_argument("--objects", default="14a1,17a1,43a1,53a1")
    f12.add_argument("--dir", default=OUT)
    f12.add_argument("--json")
    f7 = sub.add_parser("fit7")
    f7.add_argument("--objects", default="f16,f18,f20,f22,f26")
    f7.add_argument("--dir", default=OUT)
    f7.add_argument("--json")
    f6 = sub.add_parser("fit6")
    f6.add_argument("--objects", default="15a1,19a1,delta")
    f6.add_argument("--dir", default=OUT)
    f6.add_argument("--json")
    f = sub.add_parser("fit")
    f.add_argument("--curves", default="11a1,37b1")
    f.add_argument("--dir", default=OUT)
    f.add_argument("--json")
    args = ap.parse_args()
    if args.cmd == "scan6" and args.v is None:
        args.v = "2,2.5,3,3.5,4" if weight(args.obj) > 2 else "1.5,2,2.5,3,3.5"
    {"validate": cmd_validate, "curves": cmd_curves, "afe": cmd_afe, "zeros": cmd_zeros, "check": cmd_check, "scan": cmd_scan,
     "fit": cmd_fit, "precision": cmd_precision, "posthoc": cmd_posthoc, "objects6": cmd_objects6, "scan6": cmd_scan6,
     "fit6": cmd_fit6, "objects7": cmd_objects7, "fit7": cmd_fit7, "central7": cmd_central7, "zeros-merge": cmd_zeros_merge,
     "gate9": cmd_gate9, "fit9": cmd_fit9, "objects11": cmd_objects11, "fit11": cmd_fit11, "afe-both": cmd_afe_both, "fit12": cmd_fit12}[args.cmd](args)


if __name__ == "__main__":
    main()
