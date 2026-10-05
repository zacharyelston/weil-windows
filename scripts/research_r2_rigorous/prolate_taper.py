#!/usr/bin/env python3
"""Step 4 experiment [N, floating point]: the sharp prolate E-map trial with a short C^2 edge taper.

The sharp trial (R2 section 3.6) uses phi = time-limited prolate on [-lam, lam], which jumps at the edge, so psi decays
like 1/y, (A2) fails and Theorem 2's B diverges. Here phi = P(t) tau(t), t = y/lam, with
  P = alpha xi_n + xi_{n+4} (zeta, n = 2s; alpha fixed by psi(0) = 0), xi_k the prolates at c = 2 pi lam^2 (Legendre
      coefficients from k_lambda_mp.prolate_even),
  tau(t) = S((1 - |t|)/eps) on 1 - eps < |t| < 1 (S(v) = 10v^3 - 15v^4 + 6v^5, C^2), 1 inside, eps = kappa/c^2
      (the prolate's own edge scale: near t = 1 it behaves like I_0(c sqrt(2(1-t))), which varies on 1 - t ~ 1/c^2).
phi is C^2 with phi''' of bounded variation, so psi = O(y^-4), (A1)-(A4) hold and Theorem 2 applies (section 9 of
docs/R2_THEOREMS.md). This script computes, in floating point only:
  - the trial's Rayleigh quotient in rh2's basis (zeros_side), as in rq_check.py;
  - a floating estimate of Theorem 2's bound: B, A from sampled |Phi(omega)| (Lemma H form), edge sum, norm, zero
    integral (Bellotti-Wong closed form). Not certified: sampling on a grid, extrapolated omega^-4 tail.
Phi(omega) = int_{-1}^{1} phi(t) cos(omega t) dt is evaluated as the Legendre-Bessel series of P minus the transform of
P (1 - tau) on the two short edge intervals (Gauss-Legendre), at dps ~ c/ln 10 + 40 (the out-of-band value is e^-c
relative to the terms of the series).

Usage: .venv/bin/python scripts/research_r2_rigorous/prolate_taper.py --supports 1.6,2.38 --kappa 2,8 --json out.json
"""
import argparse
import json
import math
import os
import sys
import time

import mpmath as mp

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
import cert_lib as C  # noqa: E402
from flint import arb  # noqa: E402
from decay_law_mp import zeros_side  # noqa: E402
from k_lambda_mp import gl_nodes, prolate_even  # noqa: E402
from progress import Job  # noqa: E402

SUPPORT_N = {"1.6": 60, "2.38": 100, "2.6": 120, "2.99": 180}


def smooth(v):
    return v ** 3 * (10 - 15 * v + 6 * v * v)


class Taper:
    def __init__(self, c, s, kappa_eps, dps):
        self.c = mp.mpf(c)
        self.eps = mp.mpf(kappa_eps) / self.c ** 2
        ks, (bn, bn4) = prolate_even(self.c, which=(s, s + 2), parity=0)
        self.ks = ks
        self.K = ks[-1]
        self.mu_bar = (-1) ** s                      # mu = (-i)^{2s}
        self.gl = gl_nodes(40)
        self.bn, self.bn4 = bn, bn4
        self.beta = None
        self.alpha = mp.mpf(0)
        # alpha from psi(0) = 0: psi(0) = (phi(0) + mu_bar lam Phi(0))/2, linear in alpha (lam supplied later)

    def set_lam(self, lam, om_max):
        self.lam = mp.mpf(lam)
        self.om_max = om_max

        def parts(b):
            self.beta = b
            self.edge_nodes(self.om_max)
            p0 = self.P(mp.mpf(0))
            return p0 + self.mu_bar * self.lam * self.Phi(mp.mpf(0))
        u = parts(self.bn)
        v = parts(self.bn4)
        self.alpha = -v / u
        self.beta = [self.alpha * a + b for a, b in zip(self.bn, self.bn4)]
        self.edge_nodes(self.om_max)
        assert abs(self.psi(mp.mpf(0))) < mp.mpf(10) ** (-mp.mp.dps // 2)

    def P(self, t):
        return mp.fsum(b * mp.sqrt(k + mp.mpf(1) / 2) * mp.legendre(k, t) for b, k in zip(self.beta, self.ks))

    def Pvals(self, t):
        """P(t) via the Legendre recurrence (all k at once)."""
        p_prev, p = mp.mpf(1), t
        vals = {0: p_prev, 1: p}
        for k in range(1, self.K):
            p_prev, p = p, ((2 * k + 1) * t * p - k * p_prev) / (k + 1)
            vals[k + 1] = p
        return mp.fsum(b * mp.sqrt(k + mp.mpf(1) / 2) * vals[k] for b, k in zip(self.beta, self.ks))

    def tau(self, t):
        a = 1 - abs(t)
        if a >= self.eps:
            return mp.mpf(1)
        if a <= 0:
            return mp.mpf(0)
        return smooth(a / self.eps)

    def phi(self, t):
        if abs(t) >= 1:
            return mp.mpf(0)
        return self.Pvals(t) * self.tau(t)

    def sph_j(self, om, kmax):
        """j_0..j_kmax at om > 0: forward recurrence if om > kmax + 10, else Miller backward normalised by j0 or j1."""
        if om > kmax + 10:
            j = [mp.sin(om) / om, mp.sin(om) / om ** 2 - mp.cos(om) / om]
            for k in range(1, kmax):
                j.append((2 * k + 1) / om * j[k] - j[k - 1])
            return j
        ks = int(kmax + om + 40)
        a, b = mp.mpf(0), mp.mpf(10) ** -50
        out = [None] * (kmax + 1)
        for k in range(ks, 0, -1):
            a, b = b, (2 * k + 1) / om * b - a           # b = j_{k-1}, a = j_k
            if k - 1 <= kmax:
                out[k - 1] = b
            if k <= kmax:
                out[k] = a
        j0, j1 = mp.sin(om) / om, mp.sin(om) / om ** 2 - mp.cos(om) / om
        sc = j0 / out[0] if abs(j0) > abs(j1) else j1 / out[1]
        return [v * sc for v in out]

    def edge_nodes(self, om_max):
        """Fixed Gauss-Legendre nodes on [1 - eps, 1] fine enough for cos(om t), om <= om_max; P (1 - tau) cached."""
        a0 = 1 - self.eps
        npan = max(2, int(om_max * self.eps / 3) + 2)
        h = self.eps / npan
        self._en = []
        for p in range(npan):
            lo = a0 + p * h
            for t, wt in self.gl:
                tt = lo + (t + 1) * h / 2
                self._en.append((tt, wt * h / 2 * self.Pvals(tt) * (1 - self.tau(tt))))

    def edge_int(self, om, weight):
        """2 int_{1-eps}^{1} weight(t) P(t)(1 - tau(t)) dt on the cached nodes."""
        return 2 * mp.fsum(wv * weight(tt) for tt, wv in self._en)

    def Phi(self, om):
        if om == 0:
            return mp.sqrt(2) * self.beta[0] * (1 if self.ks[0] == 0 else 0) - self.edge_int(om, lambda t: 1)
        j = self.sph_j(om, self.K)
        ser = mp.fsum(b * mp.sqrt(k + mp.mpf(1) / 2) * 2 * (-1) ** (k // 2) * j[k] for b, k in zip(self.beta, self.ks))
        return ser - self.edge_int(om, lambda t: mp.cos(om * t))

    def Phi1(self, om):
        """int t phi(t) sin(om t) dt = -Phi'(om)."""
        j = self.sph_j(om, self.K + 1)
        tot = mp.mpf(0)
        for b, k in zip(self.beta, self.ks):
            nk = mp.sqrt(k + mp.mpf(1) / 2) / (2 * k + 1)
            tot += b * nk * ((k + 1) * 2 * (-1) ** (k // 2) * j[k + 1] + (k * 2 * (-1) ** ((k - 2) // 2) * j[k - 1] if k >= 1 else 0))
        return tot - self.edge_int(om, lambda t: t * mp.sin(om * t))

    def psi(self, y):
        t = y / self.lam
        return (self.phi(t) + self.mu_bar * self.lam * self.Phi(2 * mp.pi * self.lam * abs(y))) / 2

    def J3(self):
        """sum of |jumps of phi'''| (at +-1 and +-(1 - eps)), for the omega^-4 tail estimate."""
        return 2 * 60 / self.eps ** 3 * (abs(self.Pvals(mp.mpf(1))) + abs(self.Pvals(1 - self.eps)))


def run_point(sup, s, kappa, h_om=None):
    x = mp.e ** mp.mpf(sup)
    c = 2 * mp.pi * x
    lam = mp.sqrt(x)
    L = mp.log(x)
    dps = int(float(c) / math.log(10)) + 40
    mp.mp.dps = dps
    x = mp.e ** mp.mpf(sup)
    c = 2 * mp.pi * x
    lam = mp.sqrt(x)
    L = mp.log(x)
    t0 = time.time()
    tr = Taper(c, s, kappa, dps)
    Om = 20 / tr.eps + 10 * c
    tr.set_lam(lam, Om)
    eps = tr.eps
    h = h_om or (2 * mp.pi / 8)
    # sampled out-of-band integrals (trapezoid) of |Phi|(1+log) and om|Phi1|(1+log)
    n = int((Om - c) / h) + 1
    I0 = I1 = mp.mpf(0)
    prev = None
    for i in range(n + 1):
        om = c + i * h
        lg = 1 + mp.log(om / c)
        f0 = abs(tr.Phi(om)) * lg
        f1 = om * abs(tr.Phi1(om)) * lg
        w = h / 2 if i in (0, n) else h
        I0 += w * f0
        I1 += w * f1
    J3 = tr.J3()
    OmE = c + n * h
    tail0 = J3 * (OmE ** -3 / 3 * (1 + mp.log(OmE / c)) + OmE ** -3 / 9)
    tail1 = J3 * (OmE ** -2 / 2 * (1 + mp.log(OmE / c)) + OmE ** -2 / 4) * 4      # om * 4 J3/om^5
    B = (I0 + tail0) / (2 * mp.pi)
    Aint = B / 4 + (I1 + tail1) / (4 * mp.pi)
    # edge: S(sqrt x) = sum_n mu_bar lam Phi(n c)/2
    Se = mp.mpf(0)
    nmax = int(Om / c) + 1
    for k in range(1, nmax + 1):
        Se += tr.mu_bar * lam * tr.Phi(k * c) / 2
    edge = mp.sqrt(x) * abs(Se)
    A = 2 * (edge + Aint)
    # norm: 2 int_1^{sqrt x} S^2 dw, S = sum_n psi(n w) (in band) + first 5 out-of-band terms
    rx = mp.sqrt(x)
    K = int(mp.floor(rx))
    gl = gl_nodes(30)
    nrm = mp.mpf(0)
    for k in range(1, K + 1):
        a0 = max(mp.mpf(1), rx / (k + 1))
        b0 = rx / k
        if b0 <= a0:
            continue
        npan = int((b0 - a0) * 4) + 2
        hh = (b0 - a0) / npan
        for p in range(npan):
            lo = a0 + p * hh
            for t, wt in gl:
                w = lo + (t + 1) * hh / 2
                S = mp.fsum(tr.psi(m * w) for m in range(1, k + 6))
                nrm += wt * hh / 2 * S * S
    norm2 = 2 * nrm
    Z, Tc, TL = C.zero_sum_bound(arb(mp.nstr(A, 30)), arb(mp.nstr(B, 30)), 1)
    bound = mp.mpf(Z.mid().str(30, radius=False)) / norm2
    # Rayleigh quotient in rh2's basis
    N = SUPPORT_N[sup]
    hmode = L / (2 * N + 8)
    nodes, vals = [], []
    for k in range(1, K + 1):
        a0 = mp.log(max(mp.mpf(1), rx / (k + 1)))
        b0 = mp.log(rx / k)
        if b0 <= a0:
            continue
        npan = int((b0 - a0) / hmode) + 1
        hh = (b0 - a0) / npan
        for p in range(npan):
            lo = a0 + p * hh
            for t, wt in gl:
                u = lo + (t + 1) * hh / 2
                w = mp.e ** u
                S = mp.fsum(tr.psi(m * w) for m in range(1, k + 6))
                nodes.append((u, wt * hh / 2))
                vals.append(mp.e ** (u / 2) * S)
    parity = ("even", "odd")[s]
    if parity == "even":
        cv = [2 * mp.fsum(wt * v for (u, wt), v in zip(nodes, vals)) / mp.sqrt(L)]
        cv += [2 * mp.sqrt(2 / L) * mp.fsum(wt * v * mp.cos(2 * mp.pi * kk * u / L) for (u, wt), v in zip(nodes, vals)) for kk in range(1, N + 1)]
    else:
        cv = [2 * mp.sqrt(2 / L) * mp.fsum(wt * v * mp.sin(2 * mp.pi * kk * u / L) for (u, wt), v in zip(nodes, vals)) for kk in range(1, N + 1)]
    M = zeros_side(1, x, N, parity)
    cm = mp.matrix(cv)
    rq = (cm.T * M * cm)[0] / mp.fsum(v * v for v in cv)
    return {"support_2a": sup, "parity": parity, "kappa": kappa, "eps": mp.nstr(eps, 6), "c": mp.nstr(c, 8), "dps": dps,
            "alpha": mp.nstr(tr.alpha, 10), "Omega": mp.nstr(OmE, 8), "samples": n + 1,
            "B_est": mp.nstr(B, 8), "A_est": mp.nstr(A, 8), "edge": mp.nstr(edge, 6), "tail0_share": mp.nstr(tail0 / (I0 + tail0), 3),
            "norm2": mp.nstr(norm2, 10), "Tc": Tc.str(6, radius=False), "bound_est": mp.nstr(bound, 6), "N": N,
            "rq": mp.nstr(rq, 8), "bound_over_rq": mp.nstr(bound / rq, 4), "seconds": round(time.time() - t0, 1)}


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--supports", default="1.6,2.38")
    ap.add_argument("--kappa", default="2,8")
    ap.add_argument("--json", required=True)
    ap.add_argument("--resume", action="store_true", help="keep the rows already in --json and skip their tasks")
    args = ap.parse_args()
    C.setprec()
    tasks = [(sup, s, float(k)) for sup in args.supports.split(",") for s in (0, 1) for k in args.kappa.split(",")]
    rows = json.load(open(args.json)) if args.resume and os.path.exists(args.json) else []
    done = {(r["support_2a"], r["parity"], float(r["kappa"])) for r in rows}
    tasks = [t for t in tasks if (t[0], ("even", "odd")[t[1]], t[2]) not in done]
    job = Job("r2 prolate taper", total=len(tasks), args=vars(args))
    for sup, s, k in tasks:
        r = run_point(sup, s, k)
        rows.append(r)
        job.step(f"support {sup} {r['parity']} kappa={k}: bound_est={r['bound_est']} RQ={r['rq']} B={r['B_est']} A={r['A_est']} "
                 f"Tc={r['Tc']} norm2={r['norm2']} ({r['seconds']}s)")
        with open(args.json, "w") as fh:
            json.dump(rows, fh, indent=1)
    job.done()


if __name__ == "__main__":
    main()
