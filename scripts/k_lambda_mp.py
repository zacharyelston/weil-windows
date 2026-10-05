#!/usr/bin/env python3
"""Connes' prolate approximant k_λ versus the true minimiser θ_x (arXiv:2602.04022 §6.4–6.6).

k_λ(v) = E(h_λ)(v) = v^{1/2} Σ_{n ≥ 1} h_λ(n v),   v ∈ [λ⁻¹, λ],  x = λ²,
where h_λ is the combination of the even prolate functions h_{0,λ}, h_{4,λ}
(bandwidth c = 2πλ², supported on [−λ, λ]) with vanishing integral. Only
n ≤ λ² contribute. Connes' remaining step (§6.6) is that k_λ is a
sufficiently good approximation of θ_x, the minimiser of QW_λ (here: η from
scripts/connes_letter_mp.py, pole-free form, which already lives on
[λ⁻¹, λ] in the additive variable w = log v ∈ [−L/2, L/2], L = log x).

Reported:
  * odd fraction ‖k_odd‖/‖k‖ (k_λ is only approximately invariant under v ↦ 1/v),
  * sin²∠(k_even, θ_x) and the share of k_even outside the N-term basis,
  * QW_λ(P k)/‖P k‖² for the basis projection, against ε₀,
  * zeros of k̂_λ against the zeta zeros (Fact 6.4 predicts slow convergence).

Usage: .venv/bin/python scripts/k_lambda_mp.py --x 7,13 --n 100 --dps 80
"""

import argparse
import json
import os
import sys
import time

import mpmath as mp

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import connes_letter_mp as cl  # noqa: E402


def prolate_even(c, which=(0, 2), k_max=None, parity=0):
    """Prolate functions of the given parity at bandwidth c (the `which`-th of
    that parity, ascending), as Legendre coefficient vectors on k ≡ parity (mod 2)
    in the normalised Legendre basis."""
    c = mp.mpf(c)
    if k_max is None:
        k_max = int(2 * c + 80)
    ks = list(range(parity, k_max + 1, 2))
    m = len(ks)
    T = mp.matrix(m, m)
    for i, k in enumerate(ks):
        T[i, i] = k * (k + 1) + c * c * (2 * k * (k + 1) - 1) / ((2 * k + 3) * (2 * k - 1))
        if i + 1 < m:
            off = c * c * (k + 2) * (k + 1) / ((2 * k + 3) * mp.sqrt((2 * k + 1) * (2 * k + 5)))
            T[i, i + 1] = off
            T[i + 1, i] = off
    vals, vecs = mp.eigsy(T)
    order = sorted(range(m), key=lambda i: vals[i])
    return ks, [[vecs[i, order[j]] for i in range(m)] for j in which]


def legendre_even_series(beta, ks, t):
    """Σ β_i P̄_{k_i}(t), P̄_k = sqrt(k + 1/2) P_k (any parity)."""
    k_max = ks[-1]
    p_prev, p = mp.mpf(1), t  # P_0, P_1
    vals = {0: p_prev, 1: p}
    for k in range(1, k_max):
        p_next = ((2 * k + 1) * t * p - k * p_prev) / (k + 1)
        p_prev, p = p, p_next
        vals[k + 1] = p
    return mp.fsum(b * mp.sqrt(k + mp.mpf(1) / 2) * vals[k] for b, k in zip(beta, ks))


class KLambda:
    def __init__(self, x):
        self.x = mp.mpf(x)
        self.lam = mp.sqrt(self.x)
        self.L = mp.log(self.x)
        ks, (b0, b4) = prolate_even(2 * mp.pi * self.x, which=(0, 2))
        self.ks = ks
        # ∫_{−1}^{1} ψ = √2 β₀; vanishing integral: h = β₀⁽⁴⁾ ψ₀ − β₀⁽⁰⁾ ψ₄.
        self.coef = [b4[0] * a - b0[0] * b for a, b in zip(b0, b4)]
        self.n_max = int(mp.floor(self.x))

    def h(self, y):
        """h_λ(y) on [−λ, λ], zero outside."""
        if abs(y) > self.lam:
            return mp.mpf(0)
        return legendre_even_series(self.coef, self.ks, y / self.lam)

    def k(self, w):
        """k_λ(e^w), w ∈ [−L/2, L/2]."""
        v = mp.exp(w)
        n_top = int(mp.floor(self.lam / v + mp.mpf(10) ** (-mp.mp.dps // 2)))
        return mp.sqrt(v) * mp.fsum(self.h(n * v) for n in range(1, n_top + 1))

    def breakpoints(self):
        """w where a term switches on: v = λ/n, and the mirror images."""
        pts = set()
        for n in range(1, self.n_max + 1):
            w = mp.log(self.lam / n)
            if abs(w) < self.L / 2:
                pts.add(w)
                pts.add(-w)
        pts.add(-self.L / 2)
        pts.add(self.L / 2)
        return sorted(pts)


class KLambdaDH(KLambda):
    """Davenport–Heilbronn analogue: k(u) = u^{1/2} Σ a_n h(n u) with a_n the DH
    coefficients and h the first odd prolate function on [−T, T], T = λ/√5
    (bandwidth c = 2πT² = 2πx/5, matching the odd Γ-factor and conductor 5).
    The window [u₀/λ, u₀λ] is centred at the symmetry point u₀ = 1/√5; in the
    shifted variable w = log(u/u₀) it is [−L/2, L/2] with the same breakpoints."""

    def __init__(self, x):
        self.x = mp.mpf(x)
        self.lam = mp.sqrt(self.x)
        self.L = mp.log(self.x)
        self.u0 = 1 / mp.sqrt(5)
        self.T = self.lam / mp.sqrt(5)
        ks, (b1,) = prolate_even(2 * mp.pi * self.T ** 2, which=(0,), parity=1)
        self.ks = ks
        self.coef = b1
        self.n_max = int(mp.floor(self.x))
        r5 = mp.sqrt(5)
        kappa = (mp.sqrt(10 - 2 * r5) - 2) / (r5 - 1)
        chi = {1: (1, 0), 2: (0, 1), 3: (0, -1), 4: (-1, 0), 0: (0, 0)}
        self.a = {n: chi[n % 5][0] + kappa * chi[n % 5][1] for n in range(1, self.n_max + 2)}

    def h(self, y):
        if abs(y) > self.T:
            return mp.mpf(0)
        return legendre_even_series(self.coef, self.ks, y / self.T)

    def k(self, w):
        u = self.u0 * mp.exp(w)
        n_top = int(mp.floor(self.T / u + mp.mpf(10) ** (-mp.mp.dps // 2)))
        return mp.sqrt(u) * mp.fsum(self.a[n] * self.h(n * u) for n in range(1, n_top + 1))


def gl_nodes(n):
    """Gauss–Legendre nodes/weights on [−1, 1] at current precision (Golub–Welsch)."""
    J = mp.matrix(n, n)
    for k in range(1, n):
        b = k / mp.sqrt(4 * k * k - 1)
        J[k, k - 1] = b
        J[k - 1, k] = b
    vals, vecs = mp.eigsy(J)
    return sorted(((vals[i], 2 * vecs[0, i] ** 2) for i in range(n)), key=lambda p: p[0])


def quadrature_grid(kl, panel, gl):
    """Nodes (w, weight) on [−L/2, L/2], panels split at the breakpoints."""
    nodes = []
    bp = kl.breakpoints()
    for a, b in zip(bp[:-1], bp[1:]):
        m = max(1, int(mp.ceil((b - a) / panel)))
        h = (b - a) / m
        for p in range(m):
            lo = a + p * h
            for t, wt in gl:
                nodes.append((lo + (t + 1) * h / 2, wt * h / 2))
    return nodes


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--function", choices=["zeta", "dh"], default="zeta")
    ap.add_argument("--x", default="7,13")
    ap.add_argument("--n", type=int, default=100)
    ap.add_argument("--dps", type=int, default=80)
    ap.add_argument("--panel", type=float, default=0.04, help="max panel width in w")
    ap.add_argument("--gl", type=int, default=40, help="Gauss–Legendre nodes per panel")
    ap.add_argument("--zeros", type=int, default=10)
    ap.add_argument("--json")
    args = ap.parse_args()
    mp.mp.dps = args.dps
    gl = gl_nodes(args.gl)
    from progress import Job

    x_list = args.x.split(",")
    job = Job(f"k_lambda {args.function}", total=4 * len(x_list), args=vars(args))
    out = []
    for xs in x_list:
        t0 = time.time()
        x = mp.mpf(xs)
        dh = args.function == "dh"
        kl = KLambdaDH(x) if dh else KLambda(x)
        L, E, _ = cl.build_form(x, args.n, args.function)
        n1 = args.n + 1
        if dh:
            Q = E  # DH has no pole term
        else:
            v = mp.matrix(n1, 1)
            v[0] = cl.int_cos_cosh(0, L) / mp.sqrt(L)
            for j in range(1, n1):
                v[j] = mp.sqrt(2 / L) * cl.int_cos_cosh(2 * mp.pi * j / L, L)
            Q = E + 2 * v * v.T  # pole-free form (Connes' QW_λ)
        vals, vecs = mp.eigsy(Q)
        order = sorted(range(n1), key=lambda i: vals[i])
        eps0 = vals[order[0]]
        eta = [vecs[r, order[0]] for r in range(n1)]
        job.step(f"x={xs} form + eigensolve, eps0={mp.nstr(eps0, 6)}")

        nodes = quadrature_grid(kl, mp.mpf(args.panel), gl)
        kv = []
        for i_node, (w, wt) in enumerate(nodes):
            kv.append((w, wt, kl.k(w)))
            job.sub(f"x={xs} evaluating k_λ at quadrature nodes", i_node + 1, len(nodes), every=max(1, len(nodes) // 20))
        job.step(f"x={xs} k_λ evaluated at {len(nodes)} nodes")
        kmap = {}
        # k_even(w) = (k(w) + k(−w))/2: the grid is symmetric, so pair node i with its mirror.
        ws = [w for w, _, _ in kv]
        kvals = [k for _, _, k in kv]
        n_nodes = len(kv)
        k_even = [(kvals[i] + kvals[n_nodes - 1 - i]) / 2 for i in range(n_nodes)]
        k_odd = [(kvals[i] - kvals[n_nodes - 1 - i]) / 2 for i in range(n_nodes)]
        wts = [wt for _, wt, _ in kv]
        norm2 = lambda f: mp.fsum(wt * fv * fv for wt, fv in zip(wts, f))
        nk, ne, no = norm2(kvals), norm2(k_even), norm2(k_odd)

        # Project k_even on the orthonormal cosine basis b_0 = 1/√L, b_j = √(2/L) cos(ω_j w).
        def basis(j, w):
            return 1 / mp.sqrt(L) if j == 0 else mp.sqrt(2 / L) * mp.cos(2 * mp.pi * j * w / L)

        ck = [mp.fsum(wt * ke * basis(j, w) for w, wt, ke in zip(ws, wts, k_even)) for j in range(n1)]
        cnorm2 = mp.fsum(c * c for c in ck)
        outside = 1 - cnorm2 / ne
        dot = mp.fsum(a * b for a, b in zip(ck, eta))
        cos2 = dot * dot / (ne * mp.fsum(e * e for e in eta))
        cvec = mp.matrix(ck)
        qw_proj = (cvec.T * Q * cvec)[0, 0] / cnorm2
        # Overlap within the basis only (excludes the part of k outside the span).
        cos2_in = dot * dot / (cnorm2 * mp.fsum(e * e for e in eta))
        job.step(f"x={xs} projection, sin2={mp.nstr(1 - cos2, 6)}")

        # Zeros of k̂_even(t) = ∫ k_even(w) cos(tw) dw near the zeta zeros.
        def khat(t):
            return mp.fsum(wt * ke * mp.cos(t * w) for w, wt, ke in zip(ws, wts, k_even))

        zrows = []
        for j in range(1, (0 if dh else args.zeros) + 1):
            g = mp.zetazero(j).imag
            try:
                z = mp.findroot(khat, g)
                zrows.append((j, g, abs(z - g)))
            except Exception:
                zrows.append((j, g, None))
        print(
            f"x = {xs}: ε₀ = {mp.nstr(eps0, 6)};  ‖k_odd‖/‖k‖ = {mp.nstr(mp.sqrt(no / nk), 4)};"
            f"  sin²∠(k_even, θ_x) = {mp.nstr(1 - cos2, 6)}  (within basis: {mp.nstr(1 - cos2_in, 6)}; outside span: {mp.nstr(outside, 4)})"
        )
        print(f"    QW(P k)/‖P k‖² = {mp.nstr(qw_proj, 6)}   (nodes {n_nodes}, {time.time() - t0:.0f}s)")
        print("    |zero(k̂_λ) − γ_j|: " + "  ".join(f"{j}:{mp.nstr(d, 3) if d is not None else '—'}" for j, g, d in zrows))
        job.step(f"x={xs} zeros of k̂_λ ({len(zrows)})")
        job.result(f"x={xs} eps0={mp.nstr(eps0, 6)} QW(Pk)={mp.nstr(qw_proj, 6)} sin2={mp.nstr(1 - cos2, 6)}")
        out.append({
            "function": args.function, "x": float(x), "n": args.n, "dps": args.dps, "eps0": mp.nstr(eps0, 15),
            "odd_fraction": mp.nstr(mp.sqrt(no / nk), 10), "sin2_angle": mp.nstr(1 - cos2, 15),
            "sin2_angle_in_basis": mp.nstr(1 - cos2_in, 15), "outside_span": mp.nstr(outside, 10),
            "qw_projection": mp.nstr(qw_proj, 15),
            "zeros": [{"j": j, "diff": mp.nstr(d, 10) if d is not None else None} for j, g, d in zrows],
        })
    if args.json:
        with open(args.json, "w") as f:
            json.dump(out, f, indent=2)
    job.done()


if __name__ == "__main__":
    main()
