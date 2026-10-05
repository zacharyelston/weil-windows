#!/usr/bin/env python3
"""Theorem 2 with Kaiser-Bessel trial functions: an unconditional upper bound with the conjectured rate e^{-2c}.

phi is time-limited to [-lam, lam], lam = sqrt(x/q) (self-dual variable), built from the Kaiser-Bessel window of
order m, w(t) = (1 - t^2)^{m/2} I_m(beta sqrt(1 - t^2)) on [-1, 1]:
    kappa = 0 (even chi, zeta):  phi(y) = p(y/lam) w(y/lam), p(t) = 1 + a t^2 (a = 0 unless the pole needs psi(0) = 0)
    kappa = 1 (odd chi):         phi(y) = (y/lam) w(y/lam)
Closed forms (Lewitt 1990 for w; derivatives by d/domega K_nu = -omega K_{nu+1}):
    int_{-1}^{1} w(t) e^{-i omega t} dt = sqrt(2 pi) beta^m K_{m+1/2}(omega),
    K_nu(omega) = I_nu(z)/z^nu, z = sqrt(beta^2 - omega^2)   (|omega| < beta),
                = J_nu(s)/s^nu, s = sqrt(omega^2 - beta^2)   (|omega| > beta).
psi = (phi + conj(mu) phi_hat)/2 has psi_hat = mu psi with mu = (-i)^n, n = kappa + 2s (s = sector), so the E-map window
function is exactly in sector s. Out of band (|y| > lam) psi = conj(mu) phi_hat/2 is bounded by explicit envelopes:
    |J_{k+1/2}(s)| <= sqrt(2/(pi s)) sum_{j<=k} (k+j)!/(j!(k-j)! 2^j) s^{-j}   (Hankel's terminating expansion),
    |J_nu(s)/s^nu| <= 1/(2^nu Gamma(nu+1))                                    (DLMF 10.14.4, nu >= -1/2).
B and A of Theorem 2 are bounded by sums over n of integrals of these envelopes; ||f_win||^2 is computed by
quadrature. bound = zero_sum_bound(A, B)/||f_win||^2 >= lambda_min in the sector, unconditionally.

Usage: .venv/bin/python scripts/research_r2/kb_bound.py --d 1,5,-4 --xq 2,4,8 --json data/research_r2/kb_bound.json
"""
import argparse
import glob
import json
import os
import sys
import time

import mpmath as mp

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
import r2lib  # noqa: E402
from decay_law_mp import kronecker, zeros_side  # noqa: E402
from progress import Job  # noqa: E402


def K(nu, omega, beta):
    """K_nu(omega) for half-integer nu = k + 1/2: an entire function of omega^2,
    2^{-nu} sum_j ((beta^2 - omega^2)/4)^j / (j! Gamma(nu + j + 1)); elementary closed forms away from the band edge
    (I_{k+1/2}(z) = sqrt(2z/pi) i_k(z), J_{k+1/2}(s) = sqrt(2s/pi) j_k(s), spherical Bessel recurrences, stable for z, s > 2k+6)."""
    k = int(nu - mp.mpf(1) / 2)
    z2 = beta * beta - omega * omega
    if abs(z2) < (2 * k + 8) ** 2:
        tot, term, j = mp.mpf(0), 1 / mp.gamma(nu + 1), 0
        while True:
            tot += term
            j += 1
            term = term * (z2 / 4) / (j * (nu + j))
            if abs(term) < mp.eps * abs(tot) and j > 3:
                return tot / 2 ** nu
    if z2 > 0:
        z = mp.sqrt(z2)
        sh, ch = mp.sinh(z), mp.cosh(z)
        a, b = sh / z, (z * ch - sh) / (z * z)
        for i in range(1, k):
            a, b = b, a - (2 * i + 1) / z * b
        ik = a if k == 0 else b
        return mp.sqrt(2 / mp.pi) * ik / z ** k
    s_ = mp.sqrt(-z2)
    sn, cs = mp.sin(s_), mp.cos(s_)
    a, b = sn / s_, sn / (s_ * s_) - cs / s_
    for i in range(1, k):
        a, b = b, (2 * i + 1) / s_ * b - a
    jk = a if k == 0 else b
    return mp.sqrt(2 / mp.pi) * jk / s_ ** k


_HANK = {}


def _hank(k):
    if k not in _HANK:
        _HANK[k] = [mp.factorial(k + j) / (mp.factorial(j) * mp.factorial(k - j) * 2 ** j) for j in range(k + 1)]
    return _HANK[k]


def K_env(nu, omega, beta):
    """|K_nu(omega)| <= this for omega >= beta (nu half-integer >= 1/2)."""
    s2 = omega * omega - beta * beta
    cap = 1 / (2 ** nu * mp.gamma(nu + 1))
    if s2 <= 0:
        return cap
    s = mp.sqrt(s2)
    k = int(nu - mp.mpf(1) / 2)
    co = _hank(k)
    hank, sp = mp.mpf(0), mp.mpf(1)
    for j in range(k + 1):
        hank += co[j] * sp
        sp /= s
    return min(cap, mp.sqrt(2 / mp.pi) * s ** (-nu - mp.mpf(1) / 2) * hank)


class KBPsi:
    """psi on the line from a Kaiser-Bessel phi; phi_hat as sum_i c_i omega^{p_i} K_{nu + j_i}(omega) (times a phase)."""

    def __init__(self, m, beta, lam, kappa, s, pole):
        self.m, self.beta, self.lam, self.kappa = m, mp.mpf(beta), mp.mpf(lam), kappa
        self.n = kappa + 2 * s
        self.mu = (-1j) ** self.n
        self.nu = m + mp.mpf(1) / 2
        pref = self.lam * mp.sqrt(2 * mp.pi) * self.beta ** m
        if kappa == 0:
            base = [(1, 0, 0)]                     # K_nu
            t2 = [(1, 0, 1), (-1, 2, 2)]            # FT of t^2 w: K_{nu+1} - omega^2 K_{nu+2}
            self.a = mp.mpf(0)
            if pole:
                # psi(0) = 0: I_m(beta) + conj(mu) pref (K_nu(0) + a K_{nu+1}(0)) = 0 (phi(0) = w(0) = I_m(beta) ... times 0^m
                # limit: (1-t^2)^{m/2} I_m(beta sqrt(1-t^2)) at t = 0 is I_m(beta))
                cm = complex(self.mu).conjugate()
                cm = mp.mpf(cm.real)
                k0, k1 = K(self.nu, 0, self.beta), K(self.nu + 1, 0, self.beta)
                self.a = -(mp.besseli(m, self.beta) + cm * pref * k0) / (cm * pref * k1)
            self.terms = [(c, p, j) for c, p, j in base] + [(self.a * c, p, j) for c, p, j in t2]
            self.phase = 1
        else:
            self.terms = [(1, 1, 1)]                # FT of t w: -i omega K_{nu+1}
            self.phase = -1j
            self.a = mp.mpf(0)
        self.pref = pref
        # conj(mu) * phase must be real: psi real
        cp = complex(self.mu).conjugate() * self.phase
        assert abs(cp.imag) < 1e-12
        self.cp = mp.mpf(cp.real)

    def phi(self, y):
        t = y / self.lam
        if abs(t) >= 1:
            return mp.mpf(0)
        r = mp.sqrt(1 - t * t)
        w = r ** self.m * mp.besseli(self.m, self.beta * r)
        return (1 + self.a * t * t) * w if self.kappa == 0 else t * w

    def phihat_real(self, y):
        """phi_hat(y) / phase (real)."""
        om = 2 * mp.pi * self.lam * y
        return self.pref * mp.fsum(c * om ** p * K(self.nu + j, om, self.beta) for c, p, j in self.terms)

    def __call__(self, y):
        return (self.phi(y) + self.cp * self.phihat_real(y)) / 2

    def tail_env(self, y):
        """|psi(y)| <= this for |y| >= lam."""
        om = 2 * mp.pi * self.lam * abs(y)
        return self.pref / 2 * mp.fsum(abs(c) * om ** p * K_env(self.nu + j, om, self.beta) for c, p, j in self.terms)

    def tail_deriv_env(self, y):
        """|psi'(y)| <= this for |y| >= lam (d/dy = 2 pi lam d/domega; d/domega om^p K_{nu+j} = p om^{p-1} K_{nu+j} - om^{p+1} K_{nu+j+1})."""
        om = 2 * mp.pi * self.lam * abs(y)
        tot = mp.mpf(0)
        for c, p, j in self.terms:
            if p:
                tot += abs(c) * p * om ** (p - 1) * K_env(self.nu + j, om, self.beta)
            tot += abs(c) * om ** (p + 1) * K_env(self.nu + j + 1, om, self.beta)
        return self.pref / 2 * 2 * mp.pi * self.lam * tot


def gl_composite(fn, bps, nodes, hmax=mp.mpf("0.08")):
    tot = mp.mpf(0)
    for a, b in zip(bps[:-1], bps[1:]):
        p = int(mp.ceil((b - a) / hmax))
        h = (b - a) / p
        for i in range(p):
            lo = a + i * h
            tot += mp.fsum(wt * fn(lo + (t + 1) * h / 2) for t, wt in nodes) * h / 2
    return tot


def measured(D, xq, parity):
    for f in [g for g in glob.glob(f"data/connes/decay/scan*_D{D}.json") if g.split("/")[-1].split("_")[0] in ("scanhi", "scanext", "scan8")]:
        for r in json.load(open(f)):
            if abs(float(r["xq"]) - float(xq)) < 1e-9 and r["parity"] == parity:
                return mp.mpf(r["lambda_n2"]), r["n2"], r["dps"]
    return None, None, None


def window_G(psi, D, q, w, nmax_extra=60):
    """G(w) = w^{1/2} sum_n chi(n) psi(n w / sqrt q); out-of-band terms decay like y^{-m-1}."""
    sq = mp.sqrt(q)
    nin = int(psi.lam * sq / w) + 1
    tot = mp.mpf(0)
    for n in range(1, nin + nmax_extra):
        a = 1 if D == 1 else kronecker(D, n)
        if a:
            tot += a * psi(n * w / sq)
    return mp.sqrt(w) * tot


def bounds_AB(psi, D, q, x):
    """Envelope bounds for A and B (Theorem 2), in the variable y = n w / sqrt q (env0 >= |psi|, env1 >= |psi'| out of band):
       B = 2 int_{sqrt x}^inf |S(w)| dw <= 2 sqrt q sum_n (1/n) int_{n lam}^inf env0(y) dy
       |f(L/2)| e^{L/4} = sqrt x |S(sqrt x)| <= sqrt x sum_n env0(n lam)
       int |f'| e^{u/2} du <= int (|S|/2 + w |S'|) dw <= B/4 + sqrt q sum_n (1/n) int_{n lam}^inf y env1(y) dy.
    The n-sums are explicit to N0 = 40; beyond, the monotonicity of env*y^{m+1} (asserted on a geometric grid) bounds them."""
    lam = psi.lam
    sq = mp.sqrt(q)
    N0 = 40
    grid = [N0 * lam * mp.mpf(2) ** (j / mp.mpf(4)) for j in range(0, 80)]
    for env in (psi.tail_env, psi.tail_deriv_env):
        g = [env(y) * y ** (psi.m + 1) for y in grid]
        assert all(b <= a * (1 + mp.mpf("1e-12")) for a, b in zip(g, g[1:])), "env*y^(m+1) not nonincreasing"
    I1 = lambda Y: mp.quad(psi.tail_env, [Y, 2 * Y, 10 * Y, 100 * Y, mp.inf])
    I2 = lambda Y: mp.quad(lambda y: y * psi.tail_deriv_env(y), [Y, 2 * Y, 10 * Y, 100 * Y, mp.inf])
    s1 =mp.fsum(I1(n * lam) / n for n in range(1, N0 + 1))
    s2 = mp.fsum(I2(n * lam) / n for n in range(1, N0 + 1))
    s0 = mp.fsum(psi.tail_env(n * lam) for n in range(1, N0 + 1))
    m = psi.m
    # beyond N0 the envelopes decay like y^{-m-1} with env(y) y^{m+1} decreasing (checked): I1(Y) Y^m and
    # I2(Y) Y^{m-1} and env(Y) Y^{m+1} are then nonincreasing, which gives the tail sums below.
    t1 = I1(N0 * lam) * N0 ** m * mp.zeta(m + 1, N0 + 1)
    t2 = I2(N0 * lam) * N0 ** (m - 1) * mp.zeta(m, N0 + 1)
    t0 = psi.tail_env(N0 * lam) * N0 ** (m + 1) * mp.zeta(m + 1, N0 + 1)
    B = 2 * sq * (s1 + t1)
    edge = mp.sqrt(x) * (s0 + t0)
    Aint = B / 4 + sq * (s2 + t2)
    A = 2 * (edge + Aint)
    return A, B


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--d", default="1,5,8,-3,-4,-7,-20")
    ap.add_argument("--xq", default="2,3,4,5,6,7,8,10")
    ap.add_argument("--m", type=int, default=2)
    ap.add_argument("--dbeta", default="0,1,2,3", help="beta = c - dbeta; the best bound is kept")
    ap.add_argument("--dps", type=int, default=25)
    ap.add_argument("--matrix", action="store_true")
    ap.add_argument("--json")
    args = ap.parse_args()
    mp.mp.dps = args.dps
    global GL16
    from k_lambda_mp import gl_nodes
    GL16 = gl_nodes(16)
    Ds = [int(v) for v in args.d.split(",")]
    xqs = [mp.mpf(v) for v in args.xq.split(",")]
    job = Job("r2 kb bound", total=len(Ds) * len(xqs) * 2, args=vars(args))
    rows = []
    for D in Ds:
        q = 1 if D == 1 else abs(D)
        kappa = 0 if D > 0 else 1
        pole = D == 1
        for xq in xqs:
            x = xq * q
            L = mp.log(x)
            c = 2 * mp.pi * xq
            lam = mp.sqrt(xq)
            for s, parity in enumerate(("even", "odd")):
                t0 = time.time()
                best = None
                for db in [mp.mpf(v) for v in args.dbeta.split(",")]:
                    psi = KBPsi(args.m, c - db, lam, kappa, s, pole)
                    sigma = 1 if parity == "even" else -1
                    for w in (mp.mpf("1.3"), mp.mpf("2.1")):
                        g1, g2 = window_G(psi, D, q, w, 3000), window_G(psi, D, q, 1 / w, 3000)
                        assert abs(g2 - sigma * g1) <= mp.mpf("1e-5") * abs(g1), (D, xq, parity, g1, g2)
                    # ||f_win||^2 = 2 int_0^{L/2} G(e^u)^2 du, breakpoints where terms switch band (w = lam sqrt q / n)
                    bps = sorted({mp.mpf(0), L / 2} | {mp.log(lam * mp.sqrt(q) / n) for n in range(1, int(lam * mp.sqrt(q)) + 1)
                                                     if 0 < mp.log(lam * mp.sqrt(q) / n) < L / 2})
                    norm2 = 2 * gl_composite(lambda u: window_G(psi, D, q, mp.exp(u)) ** 2, bps, GL16)
                    A, B = bounds_AB(psi, D, q, x)
                    if D == 1:
                        Z = r2lib.zero_sum_bound(A, B, r2lib.zeta_Nplus, zeta=True, t0=mp.mpf("14.134"))
                    else:
                        Z = r2lib.zero_sum_bound(A, B, lambda t: r2lib.dirichlet_Nplus(t, q, kappa), zeta=False)
                    bound = Z / norm2
                    if best is None or bound < best[0]:
                        best = (bound, db, A, B, Z, norm2, psi)
                bound, db, A, B, Z, norm2, psi = best
                lam_m, n2, dps_m = measured(D, xq, parity)
                row = {"D": D, "q": q, "xq": mp.nstr(xq, 6), "x": mp.nstr(x, 10), "parity": parity, "n_index": psi.n + 4 * pole,
                       "m": args.m, "beta": mp.nstr(c - db, 10), "c": mp.nstr(c, 10), "a_t2": mp.nstr(psi.a, 8),
                       "norm2": mp.nstr(norm2, 10), "A": mp.nstr(A, 10), "B": mp.nstr(B, 10), "Tc": mp.nstr(A / B, 6),
                       "zero_sum_bound": mp.nstr(Z, 10), "bound": mp.nstr(bound, 10),
                       "ln_bound_plus_2c": mp.nstr(mp.log(bound) + 2 * c, 8),
                       "measured": mp.nstr(lam_m, 12) if lam_m is not None else None,
                       "ratio_measured_over_bound": mp.nstr(lam_m / bound, 6) if lam_m is not None else None}
                rows.append(row)
                job.step(f"D={D} x/q={mp.nstr(xq,4)} {parity}: beta=c-{mp.nstr(db,2)} bound={mp.nstr(bound,6)} "
                         f"measured={mp.nstr(lam_m,6) if lam_m else None} ratio={row['ratio_measured_over_bound']} "
                         f"ln(bound)+2c={row['ln_bound_plus_2c']} ({time.time()-t0:.0f}s)")
                if args.json:
                    with open(args.json, "w") as fh:
                        json.dump(rows, fh, indent=1)
    job.done()


if __name__ == "__main__":
    main()
