#!/usr/bin/env python3
"""Positive-side certificate for the control functions: Zhu's Theorem 1.1 reduction, generalised, in Arb.

Generalises scripts/zhu_certificate_iv.py (ζ, mpmath.iv) to
  Ψ(t) = Σ_κ Re ψ(κ + it/2) + K − Σ_{log n < 2a} (2c_n/√n) cos(t log n)        (control_cert_lib.Control)
on f supported in [−a, a] (our x = e^{2a}). With A = Σ 2|c_n|/√n and the Binet envelope
  Re ψ(κ + it/2) ≥ log(t/2) − 2κ/t² − 1/(12κt) ≥ log(t/2) − 1/t      (t ≥ 24κ²/(12κ − 1)),
Ψ ≥ β* := Σ_κ [log(T♯/2) − 1/T♯] + K − A on [T♯, ∞), so (Theorem 1.1)
  Q(f) ≥ R(f) = pole + (1/π)∫_0^{T♯} (Ψ − β)|F|² + β ‖f‖²     (β ≤ β*).
R is evaluated on the Legendre basis T_n(u) = P̄_n(u/a)/√a of one parity; in sign-adjusted coordinates
  A = s·2 p pᵀ + Σ_q c_q S(t_q) S(t_q)ᵀ + β I,  S_n(t) = √(2a(2n+1)) j_n(a t),  c_q = w_q (Ψ(t_q) − β)/π,
  p_n = ±√(2a(2n+1)) i_n(a/2) (absent without a pole),  s = +1 (even), −1 (odd).

Rigour (all in Arb ball arithmetic):
  * Gauss–Legendre nodes and weights: arb.legendre_p_root enclosures of the exact rule; quadrature error
    of the exact rule from the Chebyshev-coefficient bound on Bernstein ellipses (bound M below);
  * Ψ(t_q) via acb.digamma; j_n(a t_q) by the ratio continued fraction (n > a t_q; Pincherle) and upward
    recurrence below; i_n by arb.bessel_i;
  * Gram sums by arb_mat products (balls);
  * λ_min(A) by a floating Cholesky factor L̃ (gmpy2) and the ball residual A − μI − L̃L̃ᵀ;
  * coupling to Legendre degrees above the block (ε_B, ε_D) by the series bounds of the audit.

Usage:
  .venv/bin/python scripts/control_certificate.py --function zeta --sector even --a 0.8 --tsharp 200 --nmodes 200 \
      --workers 2 --json data/controls/cert_zeta_even_a0.8.json
"""

import argparse
import json
import math
import multiprocessing as mproc
import os
import sys
import time
from fractions import Fraction

import gmpy2
from flint import acb, arb, arb_mat, ctx, fmpq

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from control_cert_lib import Control, q  # noqa: E402
from progress import Job  # noqa: E402


def a_arb(a_str):
    fr = Fraction(a_str)
    return arb(fmpq(fr.numerator, fr.denominator))


def pack(x):
    m, e = x.mid().man_exp()
    r, re_ = x.rad().man_exp()
    return (int(m), int(e), int(r), int(re_))


def unpack(t):
    return arb((t[0], t[1]), (t[2], t[3]))


# ---------------------------------------------------------------- spherical Bessel enclosures

def sph_j(x, nmax, extra=60):
    """Enclosures of j_0(x) … j_nmax(x) for a positive ball x.
    For n ≥ nx > x (so (2n+1)/x > 2) the true ratios r_n = j_n/j_{n−1} lie in (0, 1]: by Pincherle's theorem
    r_n is the value of 1/((2n+1)/x − 1/((2n+3)/x − …)), and with every partial denominator ≥ 2 each
    approximant lies in (0, 1]. The exact identity r_n = 1/((2n+1)/x − r_{n+1}) then maps the ball [0, 1]
    at n = ntop into enclosures of all r_n, nx ≤ n ≤ ntop. Below nx: upward recurrence from j_0, j_1."""
    xm = float(x.mid())
    nx = max(1, int(xm) + 1)
    assert arb(nx) > x
    ntop = max(nmax, 2 * int(xm) + 1) + extra
    r = arb(0.5, 0.5)
    rs = [None] * (ntop + 1)
    for n in range(ntop, nx - 1, -1):
        r = 1 / ((2 * n + 1) / x - r)
        rs[n] = r
    j = [None] * (max(nmax, nx) + 1)
    sx, cx = x.sin(), x.cos()
    j[0] = sx / x
    if nx - 1 >= 1:
        j[1] = (j[0] - cx) / x
        for n in range(1, nx - 1):
            j[n + 1] = (2 * n + 1) / x * j[n] - j[n - 1]
    prev = j[nx - 1]
    for n in range(nx, nmax + 1):
        prev = prev * rs[n]
        j[n] = prev
    return j[: nmax + 1]


def sph_i(z, n):
    """i_n(z) = √(π/(2z)) I_{n+1/2}(z)."""
    return (arb.pi() / (2 * z)).sqrt() * z.bessel_i(arb(n) + q(1, 2))


# ---------------------------------------------------------------- node worker

_G = {}


def _init_worker(cfg):
    ctx.prec = cfg["prec"]
    ctx.threads = 1
    _G.update(cfg)
    _G["ctrl"] = Control(cfg["function"])
    _G["a"] = a_arb(cfg["a"])
    _G["comb"] = _G["ctrl"].comb(2 * _G["a"])
    _G["beta"] = unpack(cfg["beta"])
    _G["inv_pi"] = 1 / arb.pi()
    _G["norm"] = [(2 * _G["a"] * (2 * n + 1)).sqrt() for n in cfg["degrees"]]
    # Ball arithmetic overestimates the upward recurrence below the turning point by up to ~0.7 bits per
    # unit of a·t, so nodes and Bessel values are computed with that many extra bits and rounded after.
    _G["prec_hi"] = cfg["prec"] + int(0.9 * float(_G["a"].upper()) * cfg["tsharp"]) + 64
    ctx.prec = _G["prec_hi"]
    _G["gl"] = [arb.legendre_p_root(cfg["gl"], i, weight=True) for i in range(cfg["gl"])]
    _G["h"] = arb(fmpq(cfg["tsharp"], cfg["panels"]))
    _G["a_hi"] = a_arb(cfg["a"])
    ctx.prec = cfg["prec"]


def _panel_task(rng):
    """Gram contribution of panels rng[0] … rng[1]−1 (upper triangle, packed) and Ψ range on the nodes."""
    k0, k1 = rng
    degrees = _G["degrees"]
    nmax = max(degrees)
    a, beta, h, ctrl, comb = _G["a"], _G["beta"], _G["h"], _G["ctrl"], _G["comb"]
    N = len(degrees)
    gram = arb_mat(N, N)
    psi_lo, psi_hi = None, None
    cols_S, cols_C = [], []

    def flush():
        nonlocal gram
        if not cols_S:
            return
        m = len(cols_S)
        Sm = arb_mat([[cols_S[c][i] for c in range(m)] for i in range(N)])
        Cm = arb_mat([[cols_C[c][i] for c in range(m)] for i in range(N)])
        gram = gram + Cm * Sm.transpose()
        cols_S.clear()
        cols_C.clear()

    P, Phi = _G["prec"], _G["prec_hi"]
    for k in range(k0, k1):
        for xi, wi in _G["gl"]:
            ctx.prec = Phi
            t_hi = (k + q(1, 2)) * h + xi * h / 2
            j_hi = sph_j(_G["a_hi"] * t_hi, nmax)
            w_hi = wi * h / 2
            ctx.prec = P
            t, w = +t_hi, +w_hi
            j = [+v for v in j_hi]
            psi = ctrl.Psi(t, comb)
            lo, hi = float(psi.lower()), float(psi.upper())
            psi_lo = lo if psi_lo is None else min(psi_lo, lo)
            psi_hi = hi if psi_hi is None else max(psi_hi, hi)
            c = w * (psi - beta) * _G["inv_pi"]
            S = [_G["norm"][i] * j[n] for i, n in enumerate(degrees)]
            cols_S.append(S)
            cols_C.append([c * s for s in S])
            if len(cols_S) >= 256:
                flush()
    flush()
    packed = [[pack(gram[i, j]) for j in range(i, N)] for i in range(N)]
    return packed, psi_lo, psi_hi, k1 - k0


# ---------------------------------------------------------------- analytic bounds

def ellipse_M(ctrl, a, tsharp, beta, comb, nmax, rho, h):
    """Bound on |(Ψ(z) − β) S_m(z) S_n(z)/π| on the panels' Bernstein ellipses (semi-minor b, semi-major a_e).
    Ψ(z) = Σ_κ ½[ψ(κ + iz/2) + ψ(κ − iz/2)] + K − Σ c_n' cos(z log n); with z = t + is, |s| ≤ b < 2κ:
      w = κ ± iz/2 has Re w ≥ κ − b/2;  ψ(w) = ψ(1 + w) − 1/w;
      |ψ(v)| ≤ log|v| + π/2 + 1/(2 Re v) + 1/(12 (Re v)²)   (Re v ≥ 1; Binet);
      |cos(z log n)| ≤ cosh(b log n);  |j_n(a z)| ≤ e^{a|Im z|}  ⇒  |S_n| ≤ √(2a(2n+1)) e^{ab}."""
    b = h * (rho - 1 / rho) / 4
    a_e = h * (rho + 1 / rho) / 4
    G = abs(ctrl.K) + abs(beta)
    for k in ctrl.kappas:
        kk = arb(k)
        assert b < 2 * kk
        re_w = kk - b / 2
        v_re = 1 + re_w
        v_abs = 1 + kk + b / 2 + (tsharp + a_e) / 2
        G += 1 / re_w + v_abs.log() + arb.pi() / 2 + 1 / (2 * v_re) + 1 / (12 * v_re * v_re)
    for _, cn, ln in comb:
        G += abs(cn) * (b * ln).cosh()
    S2 = 2 * a * (2 * nmax + 1) * (2 * a * b).exp()
    return G * S2 / arb.pi(), G, b, a_e


def coupling_bounds(a, degrees, tsharp, G, pole_abs):
    """ε_B (block ↔ Legendre degrees ≥ n0) and ε_D (tail block − β I), as in zhu_certificate_iv.py.
    For t ∈ [0, T♯]: |S_n| ≤ d_n = √(2a(2n+1)) (aT♯)^n/(2n+1)!!, |S_m| ≤ √(2a(2m+1)),
    |p_n| ≤ e_n = √(2a(2n+1)) (a/2)^n/(2n+1)!! · exp((a/2)²/(2(2n+3))), |Ψ − β| ≤ G."""
    T = arb(tsharp)
    n0 = max(degrees) + 2

    def dfact(n):
        v = arb(1)
        for k in range(1, 2 * n + 2, 2):
            v *= k
        return v

    def d(n):
        return (2 * a * (2 * n + 1)).sqrt() * (a * T) ** n / dfact(n)

    def e(n):
        z = a / 2
        return (2 * a * (2 * n + 1)).sqrt() * z ** n / dfact(n) * (z * z / (2 * (2 * n + 3))).exp()

    def q_of(n, y):
        return (arb(2 * n + 5) / (2 * n + 1)).sqrt() * y ** 2 / ((2 * n + 3) * (2 * n + 5))

    n1 = n0
    while not q_of(n1, a * T) < q(1, 2):
        n1 += 2
    q_d, q_e = q_of(n1, a * T), q_of(n1, a / 2)
    assert q_e < q(1, 2)
    d_list = [d(n) for n in range(n0, n1, 2)]
    sum_d = sum(d_list, arb(0)) + d(n1) / (1 - q_d)
    sum_e = sum((e(n) for n in range(n0, n1, 2)), arb(0)) + e(n1) / (1 - q_e)
    d_max = max(d_list + [d(n1)], key=lambda v: float(v.upper()))
    e0 = e(n0)
    sum_p = sum(pole_abs, arb(0))
    max_p = max(pole_abs, key=lambda v: float(v.upper())) if pole_abs else arb(0)
    sum_sqrt = sum(((2 * a * (2 * m + 1)).sqrt() for m in degrees), arb(0))
    max_sqrt = (2 * a * (2 * max(degrees) + 1)).sqrt()
    c = G * T / arb.pi()
    B_inf = 2 * e0 * sum_p + c * d_max * sum_sqrt
    B_one = 2 * sum_e * max_p + c * sum_d * max_sqrt
    eps_B = (B_inf * B_one).sqrt()
    eps_D = 2 * e0 * sum_e + c * d_max * sum_d
    return {"n0": n0, "n1": n1, "d_n0": d(n0), "e_n0": e0, "d_max": d_max, "eps_B": eps_B, "eps_D": eps_D}


# ---------------------------------------------------------------- floating Cholesky and ball residual

def to_mpfr(x):
    m, e = x.mid().man_exp()
    return gmpy2.mpfr(int(m)) * gmpy2.mpfr(2) ** int(e) if int(m) != 0 else gmpy2.mpfr(0)


def mpfr_to_arb(v):
    if v == 0:
        return arb(0)
    m, e = v.as_integer_ratio()
    # exact: denominator is a power of two
    k = e.bit_length() - 1
    assert e == 1 << k
    return arb((int(m), -k))


def fs(v, d=15):
    """Decimal string of an mpfr (d significant digits)."""
    return mpfr_to_arb(v).str(d, radius=False)


def cholesky(Mf):
    """Floating Cholesky of a list-of-lists of mpfr (current gmpy2 precision). Returns L or None."""
    N = len(Mf)
    L = [[gmpy2.mpfr(0)] * N for _ in range(N)]
    for j in range(N):
        Lj = L[j]
        s = Mf[j][j] - gmpy2.fsum([Lj[k] * Lj[k] for k in range(j)]) if j else Mf[j][j]
        if s <= 0:
            return None, j
        d = gmpy2.sqrt(s)
        Lj[j] = d
        for i in range(j + 1, N):
            Li = L[i]
            s = Mf[i][j] - gmpy2.fsum([Li[k] * Lj[k] for k in range(j)]) if j else Mf[i][j]
            Li[j] = s / d
    return L, None


def tri_solve(L, b, k=None):
    """Solve L Lᵀ x = b on the leading k×k block."""
    N = k or len(L)
    y = [gmpy2.mpfr(0)] * N
    for i in range(N):
        y[i] = (b[i] - gmpy2.fsum([L[i][j] * y[j] for j in range(i)])) / L[i][i]
    x = [gmpy2.mpfr(0)] * N
    for i in range(N - 1, -1, -1):
        x[i] = (y[i] - gmpy2.fsum([L[j][i] * x[j] for j in range(i + 1, N)])) / L[i][i]
    return x


def smallest_eig(L, Mf, k=None, iters=40, tol_digits=30):
    """λ_min of the leading k×k block of Mf by inverse iteration with its Cholesky factor."""
    N = k or len(L)
    v = [gmpy2.mpfr(1) / (1 + i) for i in range(N)]
    lam = None
    for _ in range(iters):
        x = tri_solve(L, v, N)
        nrm = gmpy2.sqrt(gmpy2.fsum([c * c for c in x]))
        v = [c / nrm for c in x]
        Mv = [gmpy2.fsum([Mf[i][j] * v[j] for j in range(N)]) for i in range(N)]
        new = gmpy2.fsum([v[i] * Mv[i] for i in range(N)])
        if lam is not None and abs(new - lam) <= abs(new) * gmpy2.mpfr(10) ** (-tol_digits):
            lam = new
            break
        lam = new
    res = gmpy2.sqrt(gmpy2.fsum([(Mv[i] - lam * v[i]) ** 2 for i in range(N)]))
    return lam, res, v


def residual_bound(A, L, mu, prec):
    """Upper bound on ‖A_true − μI − L̃L̃ᵀ‖_∞ over every symmetric A_true in the ball matrix A."""
    old = ctx.prec
    ctx.prec = prec
    N = A.nrows()
    Lm = arb_mat([[mpfr_to_arb(L[i][j]) for j in range(N)] for i in range(N)])
    R = A - Lm * Lm.transpose()
    mu_a = mpfr_to_arb(mu)
    worst = arb(0)
    for i in range(N):
        row = arb(0)
        for j in range(N):
            rij = R[i, j] - (mu_a if i == j else 0)
            row += arb(abs(rij).upper())
        if row.upper() > worst.upper():
            worst = arb(row.upper())
    ctx.prec = old
    return worst


# ---------------------------------------------------------------- main

def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--function", choices=["zeta", "ftstar", "dh", "z1", "zetak"], required=True)
    ap.add_argument("--sector", choices=["even", "odd"], default="even")
    ap.add_argument("--a", default="0.8", help="half-width a of the support [−a, a] (exact decimal); x = e^{2a}")
    ap.add_argument("--tsharp", type=int, default=200)
    ap.add_argument("--nmodes", type=int, default=200)
    ap.add_argument("--panels", type=int, default=None, help="default 4·T♯ (width 1/4)")
    ap.add_argument("--gl", type=int, default=32)
    ap.add_argument("--rho", default="4.2360679774997896964091736687", help="Bernstein ellipse parameter (2+√5 ⇒ b = h)")
    ap.add_argument("--prec", type=int, default=256)
    ap.add_argument("--chol-prec", type=int, default=256, help="bits for the floating Cholesky (gmpy2)")
    ap.add_argument("--mu", default="auto", help="Cholesky shift(s), comma separated; 'auto' = λ₁(1 − 1e-9)")
    ap.add_argument("--workers", type=int, default=2)
    ap.add_argument("--tasks", type=int, default=None)
    ap.add_argument("--blocks", default="")
    ap.add_argument("--json")
    args = ap.parse_args()

    ctx.prec = args.prec
    ctx.threads = min(2, args.workers)
    gmpy2.get_context().precision = args.chol_prec
    panels = args.panels or 4 * args.tsharp
    first = 0 if args.sector == "even" else 1
    degrees = [first + 2 * i for i in range(args.nmodes)]
    sign = 1 if args.sector == "even" else -1
    tag = f"control cert {args.function} {args.sector} a={args.a} T{args.tsharp} N{args.nmodes}"
    job = Job(tag, total=None, args=vars(args))
    t0 = time.time()

    ctrl = Control(args.function)
    a = a_arb(args.a)
    ell = 2 * a
    comb = ctrl.comb(ell)
    A_comb = sum((abs(cn) for _, cn, _ in comb), arb(0))
    A_up = arb(A_comb.upper())
    beta_star = ctrl.envelope_beta(args.tsharp, A_up)
    lo = beta_star.lower()
    beta_t = arb((int((lo * arb(2) ** 190).floor().unique_fmpz()), -190))
    assert beta_t <= beta_star
    T1 = None
    job.log(f"{args.function} kappas={[str(k) for k in ctrl.kappas]} K={ctrl.K.str(15)} pole={ctrl.pole}; a={args.a}, x=e^(2a)={ell.exp().str(15)}")
    job.log(f"comb n = {[n for n, _, _ in comb]}; 2c_n/sqrt(n) = {[cn.str(8) for _, cn, _ in comb]}")
    job.log(f"A = {A_comb.str(20)}; beta* = {beta_star.str(20)}; beta~ = {beta_t.str(25)}")

    # quadrature error of the exact Gauss rule on every panel (Chebyshev bound, even k only)
    rho = arb(args.rho)
    h = arb(fmpq(args.tsharp, panels))
    M, Gbound, b_semi, a_semi = ellipse_M(ctrl, a, args.tsharp, beta_t, comb, max(degrees), rho, h)
    n = args.gl
    eps_Q = arb(args.tsharp) * M * (2 + arb(2) / (4 * n * n - 1)) * rho ** (-2 * n) / (1 - rho ** (-2))
    c_err = args.nmodes * eps_Q
    job.log(f"GL{n} exact rule; ellipse rho={args.rho[:8]}, b={b_semi.str(6)}; |Psi-beta| <= G = {Gbound.str(6)}; M = {M.str(6)}; "
            f"eps_Q = {eps_Q.str(4)}; c_err = N*eps_Q = {c_err.str(4)}")

    pole = []
    if ctrl.pole:
        for d in degrees:
            sg = (-1) ** (d // 2) if args.sector == "even" else (-1) ** ((d - 1) // 2)
            pole.append(sg * (2 * a * (2 * d + 1)).sqrt() * sph_i(a / 2, d))
        job.log(f"pole p_0 = {pole[0].str(20)}; |p| = {sum((p * p for p in pole), arb(0)).sqrt().str(15)}")
    pole_abs = [abs(p) for p in pole]
    coup = coupling_bounds(a, degrees, args.tsharp, Gbound, pole_abs)
    job.log(f"coupling: n0 = {coup['n0']}, d_n0 = {coup['d_n0'].str(4)}, eps_B = {coup['eps_B'].str(4)}, eps_D = {coup['eps_D'].str(4)}")

    # node sums
    ntask = args.tasks or max(4 * args.workers, min(64, panels // 16))
    bounds = [round(i * panels / ntask) for i in range(ntask + 1)]
    ranges = [(bounds[i], bounds[i + 1]) for i in range(ntask) if bounds[i + 1] > bounds[i]]
    job.total = len(ranges)
    cfg = {"function": args.function, "a": args.a, "degrees": degrees, "beta": pack(beta_t), "prec": args.prec,
           "gl": args.gl, "tsharp": args.tsharp, "panels": panels}
    N = len(degrees)
    gram = [[arb(0)] * N for _ in range(N)]
    psi_lo, psi_hi = None, None
    ctxm = mproc.get_context("fork")
    with ctxm.Pool(args.workers, initializer=_init_worker, initargs=(cfg,)) as pool:
        for packed, plo, phi, cnt in pool.imap_unordered(_panel_task, ranges):
            for i in range(N):
                row = packed[i]
                gi = gram[i]
                for jj, tpl in enumerate(row):
                    gi[i + jj] = gi[i + jj] + unpack(tpl)
            psi_lo = plo if psi_lo is None else min(psi_lo, plo)
            psi_hi = phi if psi_hi is None else max(psi_hi, phi)
            job.step(f"{cnt} panels")
    job.log(f"nodes done in {time.time() - t0:.0f}s; Psi on nodes in [{psi_lo:.6f}, {psi_hi:.6f}]")

    A = arb_mat(N, N)
    max_rad = 0.0
    for i in range(N):
        for j in range(i, N):
            v = gram[i][j]
            if pole:
                v = v + sign * 2 * pole[i] * pole[j]
            if i == j:
                v = v + beta_t
            A[i, j] = v
            A[j, i] = v
            max_rad = max(max_rad, float(v.rad()))
    job.log(f"A assembled; max entry radius {max_rad:.3e}")

    Mf = [[to_mpfr(A[i, j]) for j in range(N)] for i in range(N)]
    t1 = time.time()
    L0, fail = cholesky(Mf)
    eigs = {}
    block_rows = []
    if L0 is None:
        job.result(f"midpoint matrix not positive definite (Cholesky breaks at {fail}); no certificate")
        lam1 = None
    else:
        lam1, res1, _ = smallest_eig(L0, Mf)
        eigs = {"lambda_1": float(lam1), "lambda_1_str": fs(lam1, 15), "residual": fs(res1, 3)}
        job.result(f"lambda_1(A_mid) = {fs(lam1, 15)} (residual {fs(res1, 2)}; {time.time() - t1:.0f}s)")
        for kb in [int(x) for x in args.blocks.split(",") if x]:
            if kb >= N:
                continue
            lb, _, _ = smallest_eig(L0, Mf, k=kb)
            block_rows.append({"k": kb, "max_degree": degrees[kb - 1], "lambda_1": fs(lb, 12)})
            job.result(f"leading block k={kb} (degree <= {degrees[kb - 1]}): lambda_1 = {fs(lb, 12)}")

    certs = []
    best_cert = None
    mus = [s for s in args.mu.split(",") if s and s != "auto"]
    if "auto" in args.mu and lam1 is not None and lam1 > 0:
        mus.append(fs(lam1 * (1 - gmpy2.mpfr("1e-9")), 12))
    def chol_at(mu):
        return cholesky([[Mf[i][j] - (mu if i == j else 0) for j in range(N)] for i in range(N)])

    if "auto" in args.mu and lam1 is not None and lam1 > 0 and chol_at(gmpy2.mpfr(mus[-1]))[0] is None:
        # inverse iteration had not converged (an eigenvalue cluster near β): step the shift down
        base = gmpy2.mpfr(mus[-1])
        for fac in ["1e-6", "1e-4", "1e-3", "1e-2", "0.1", "0.5"]:
            trial = base * (1 - gmpy2.mpfr(fac))
            if chol_at(trial)[0] is not None:
                break
        job.log(f"auto shift {mus[-1]} fails; using mu = {fs(trial, 12)} = auto·(1 − {fac})")
        mus[-1] = fs(trial, 12)
    for mu_s in mus:
        mu = gmpy2.mpfr(mu_s)
        L, fail = chol_at(mu)
        if L is None:
            job.result(f"Cholesky at mu = {mu_s}: FAILED at {fail}")
            certs.append({"mu": mu_s, "ok": False})
            continue
        r = residual_bound(A, L, mu, 2 * args.chol_prec + 128)
        lamA = mpfr_to_arb(mu) - r - c_err
        l1, l2 = lamA.lower(), (beta_t - coup["eps_D"]).lower()      # exact lower endpoints
        lam_cert_lo = ((l1 if l1 < l2 else l2) - coup["eps_B"]).lower()
        ok = lam_cert_lo > 0
        if ok and (best_cert is None or lam_cert_lo > best_cert):
            best_cert = lam_cert_lo
        job.result(f"Cholesky at mu = {mu_s}: residual_inf <= {r.str(4)}; lambda_min(A) >= {lamA.lower().str(15)}; "
                   f"certified Q >= {lam_cert_lo.str(15)} ||f||^2  ({'POSITIVE' if ok else 'not positive'})")
        certs.append({"mu": mu_s, "ok": True, "residual_inf_upper": r.str(6), "lambda_A_lower": lamA.lower().str(20),
                      "Q_lower": lam_cert_lo.str(20), "positive": bool(ok)})

    out = {
        "function": args.function, "sector": args.sector, "a": args.a, "x": ell.exp().str(20),
        "kappas": [str(k) for k in ctrl.kappas], "K": ctrl.K.str(20), "pole": ctrl.pole,
        "comb": [{"n": n_, "2c_n/sqrt(n)": cn.str(20)} for n_, cn, _ in comb],
        "tsharp": args.tsharp, "nmodes": args.nmodes, "degrees": [degrees[0], degrees[-1]], "panels": panels, "gl": args.gl,
        "prec_bits": args.prec, "chol_prec_bits": args.chol_prec,
        "A": A_comb.str(20), "beta_star": beta_star.str(20), "beta_used": beta_t.str(30),
        "rho": args.rho, "ellipse_b": b_semi.str(8), "G_bound": Gbound.str(8), "M_bound": M.str(8),
        "eps_Q": eps_Q.str(5), "c_err": c_err.str(5),
        "pole_p0": pole[0].str(25) if pole else None,
        "coupling": {k: (v.str(5) if not isinstance(v, int) else v) for k, v in coup.items()},
        "psi_on_nodes": [psi_lo, psi_hi], "max_entry_radius": max_rad, "A00": A[0, 0].str(25, radius=True),
        "eig_mid": eigs, "certificates": certs, "leading_blocks": block_rows,
        "elapsed_s": round(time.time() - t0, 1),
    }
    out["Q_lower_best"] = best_cert.str(20) if best_cert is not None else None
    out["certified_positive"] = best_cert is not None
    job.result(f"elapsed {time.time() - t0:.0f}s")
    if args.json:
        os.makedirs(os.path.dirname(os.path.abspath(args.json)), exist_ok=True)
        with open(args.json, "w") as fh:
            json.dump(out, fh, indent=2)
    job.done()


if __name__ == "__main__":
    main()
