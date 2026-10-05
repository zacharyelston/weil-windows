#!/usr/bin/env python3
"""How fast does the window minimum λ_min(x) of Weil's form decay, and what sets the rate?

For ζ the zeros-side minima fit ln λ_min ≈ −4πx + γ log x + β (x = e^{2a}): the Slepian prolate rate e^{−2c} at
c = 2πx. This script measures the same quantity for real primitive Dirichlet L-functions L(s, χ_D), conductor
q = |D|, to see how the rate depends on q. Two candidate laws:
  zero counting  F (type L/2) can vanish on the zeros only below T* = 2πx/q, where the zero density
                 (1/2π)log(qT/2π) reaches L/2π:  ln λ ≈ −4πx/q.
  self-dual      the conductor rescales the window variable by √q:  ln λ ≈ −4πx/√q.

The zeros-side form is additive over Euler factors, so ζ_K = ζ·L(χ_{−20}) is the sum of two forms here.
Blocks: χ even uses ζ's archimedean form (Γ_R(s)); χ odd uses the DH block (Γ_R(s+1), built for conductor 5).
The conductor enters as + log q on the diagonal, and only ζ carries the pole vector.

Commands:
  check  L(χ_5) (even) and L(χ_{−4}) (odd): zeros-side diagonals against Σ over zeros found from Hardy's Z.
  scan   λ_min over an x-grid for each D, both sectors, at two basis sizes (convergence).
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


def kronecker(D, n):
    """Kronecker symbol (D/n) for n ≥ 0."""
    if n == 0:
        return 1 if abs(D) == 1 else 0
    if n == 1:
        return 1
    res = 1
    while n % 2 == 0:
        n //= 2
        if D % 2 == 0:
            return 0
        res *= 1 if D % 8 in (1, 7) else -1
    a, m = D % n, n
    while a:
        while a % 2 == 0:
            a //= 2
            if m % 8 in (3, 5):
                res = -res
        a, m = m, a
        if a % 4 == 3 and m % 4 == 3:
            res = -res
        a %= m
    return res if m == 1 else 0


def zeros_side(D, x, n, parity):
    """rh2's zeros-side matrix for L(s, χ_D) (D = 1: ζ) on the periodic basis (even: k = 0..n, odd: k = 1..n)."""
    x = mp.mpf(x)
    L = mp.log(x)
    if D == "dh":           # Davenport–Heilbronn: conductor 5, Γ_R(s+1), no Euler product, off-line zeros
        return cl.build_form(x, n, "dh", parity=parity)[1]
    q = abs(D)
    pp = cl.prime_powers(int(x))
    if D == 1:
        _, E, _ = cl.build_form(x, n, "zeta", parity=parity)
        if parity == "even":
            v = mp.matrix([cl.int_cos_cosh(2 * mp.pi * k / L, L) * (1 / mp.sqrt(L) if k == 0 else mp.sqrt(2 / L)) for k in range(n + 1)])
            return E + 2 * v * v.T
        w = cl.odd_pole_vector(L, n)
        return E - 2 * w * w.T
    terms = []
    for nn, lp in pp:
        c = kronecker(D, nn)
        if c:
            terms.append((mp.log(nn), c * lp / mp.sqrt(nn)))
    if D > 0:      # even character: Γ_R(s), conductor log q added
        _, E, _ = cl.build_form(x, n, "zeta", terms=terms, parity=parity)
        return E + mp.log(q) * mp.eye(E.rows)
    _, E, _ = cl.build_form(x, n, "dh", terms=terms, parity=parity)       # Γ_R(s+1), log 5 built in
    return E + (mp.log(q) - mp.log(5)) * mp.eye(E.rows)


def lam_min_inv(M, iters=6):
    """Smallest-|λ| eigenvalue by inverse iteration with FLINT's solver (midpoints). The forms here are positive and
    the gap λ_2/λ_1 is huge, so 3–4 iterations converge; we stop when the Rayleigh quotient is stable to 1e-12."""
    from flint import arb, arb_mat, ctx
    ctx.prec = int(mp.mp.dps * 3.33) + 20
    n = M.rows
    A = arb_mat([[arb(mp.nstr(M[i, j], mp.mp.dps + 5)) for j in range(n)] for i in range(n)])
    v = arb_mat([[1] for _ in range(n)])
    prev = None
    for _ in range(iters):
        w = A.solve(v)
        nrm = sum((w[i, 0] ** 2 for i in range(n)), arb(0)).sqrt()
        v = arb_mat([[w[i, 0] / nrm] for i in range(n)])
        lam = mp.mpf((v.transpose() * A * v)[0, 0].mid().str(mp.mp.dps, radius=False))
        if prev is not None and abs(lam - prev) <= abs(lam) * mp.mpf("1e-12"):
            break
        prev = lam
    return lam


def lam_min_acb(M):
    """Smallest eigenvalue from FLINT's full eigensolver (midpoints): safe when negative eigenvalues may exist."""
    from flint import acb_mat, arb, arb_mat, ctx
    ctx.prec = int(mp.mp.dps * 3.33) + 20
    A = arb_mat([[arb(mp.nstr(M[i, j], mp.mp.dps + 5)) for j in range(M.cols)] for i in range(M.rows)])
    vals = [mp.mpf(z.real.mid().str(mp.mp.dps, radius=False)) for z in acb_mat(A).eig(nonstop=True)]
    if all(mp.isfinite(v) for v in vals):
        return min(vals)
    # FLINT could not isolate a cluster (indeterminate balls): min() over NaNs is order-dependent, so fall back
    return min(mp.eigsy(M, eigvals_only=True))


def hardy_zeros(D, tmax, step=mp.mpf("0.1"), job=None):
    """Ordinates 0 < γ ≤ tmax of L(s, χ_D) from sign changes of Λ(½ + it) (real for real primitive χ, root number 1)."""
    q = abs(D)
    chi = [kronecker(D, r) for r in range(q)]
    kappa = 0 if D > 0 else 1

    def Z(t):
        s = mp.mpc(mp.mpf(1) / 2, t)
        g = (mp.mpf(q) / mp.pi) ** (s / 2) * mp.gamma((s + kappa) / 2)
        return mp.re(g * mp.dirichlet(s, chi))

    out, t, zt = [], step, Z(step)
    while t < tmax:
        t2 = t + step
        z2 = Z(t2)
        if zt * z2 < 0:
            out.append(mp.findroot(Z, (t, t2), solver="anderson"))
        t, zt = t2, z2
        if job is not None:
            job.sub(f"zeros of L(χ_{D})", int(t), int(tmax), every=10)
    return out


def F_anti_even(i, L, t):
    """F(t) for the antiperiodic even basis function √(2/L) cos(2π(i+½)u/L) (vanishes at the edges: F ~ 1/t²)."""
    w = 2 * mp.pi * (i + mp.mpf(1) / 2) / L
    S = lambda y: L / 2 if abs(y) < mp.mpf(10) ** (-mp.mp.dps + 5) else mp.sin(y * L / 2) / y
    return mp.sqrt(2 / L) * (S(t - w) + S(t + w))


def F_per_odd(k, L, t):
    """F(t)/i for √(2/L) sin(2πku/L) (vanishes at the edges)."""
    w = 2 * mp.pi * k / L
    S = lambda y: L / 2 if abs(y) < mp.mpf(10) ** (-mp.mp.dps + 5) else mp.sin(y * L / 2) / y
    return mp.sqrt(2 / L) * (S(t - w) - S(t + w))


def F_per_even(c, L, t):
    """F(t) for Σ c_k b_k on the periodic even basis b_0 = 1/√L, b_k = √(2/L) cos(2πku/L)."""
    S = lambda y: L / 2 if abs(y) < mp.mpf(10) ** (-mp.mp.dps + 5) else mp.sin(y * L / 2) / y
    v = c[0] * 2 * S(t) / mp.sqrt(L)
    for k in range(1, len(c)):
        w = 2 * mp.pi * k / L
        v += c[k] * mp.sqrt(2 / L) * (S(t - w) + S(t + w))
    return v


def cmd_check(args):
    """Both sectors on edge-vanishing test functions, so the truncated zero sum converges like 1/T² or faster.
    Odd: single sin modes and b_2 − b_3. Even: b_1 + b_2 and b_2 − b_4 (cos(2πk·½) = (−1)^k cancels at the edges)."""
    mp.mp.dps = 25
    job = Job("decay law check", args=vars(args))
    absdiffs = []
    x = mp.mpf(args.x)
    L = mp.log(x)
    for D in [int(d) for d in args.d.split(",")]:
        t0 = time.time()
        gam = hardy_zeros(D, args.tmax, job=job)
        job.log(f"D = {D}: {len(gam)} zeros up to {args.tmax} ({time.time() - t0:.0f}s); first {[mp.nstr(g, 8) for g in gam[:3]]}")
        M = zeros_side(D, x, 6, "odd")
        for k in (2, 4):
            Fz = mp.fsum(2 * F_per_odd(k, L, g) ** 2 for g in gam)
            Q = M[k - 1, k - 1]
            absdiffs.append(abs(Q - Fz))
            job.result(f"D = {D}, x = {args.x}, odd k = {k}: Q_kk = {mp.nstr(Q, 12)};  2Σ_ρ F(γ)² = {mp.nstr(Fz, 12)};  rel diff {mp.nstr(abs(Q - Fz) / abs(Q), 3)}")
        c = [0, 1, -1, 0, 0, 0]
        cv = mp.matrix(c)
        qv = (cv.T * M * cv)[0]
        Fz = mp.fsum(2 * mp.fsum(ci * F_per_odd(i + 1, L, g) for i, ci in enumerate(c)) ** 2 for g in gam)
        absdiffs.append(abs(qv - Fz))
        job.result(f"D = {D}, x = {args.x}, odd c = {c}: cᵀQc = {mp.nstr(qv, 12)};  2Σ_ρ F(γ)² = {mp.nstr(Fz, 12)};  rel diff {mp.nstr(abs(qv - Fz) / abs(qv), 3)}")
        Me = zeros_side(D, x, 6, "even")
        for c in ([0, 1, 1, 0, 0, 0, 0], [0, 0, 1, 0, -1, 0, 0]):
            cv = mp.matrix(c)
            qv = (cv.T * Me * cv)[0]
            Fz = mp.fsum(2 * F_per_even(c, L, g) ** 2 for g in gam)
            absdiffs.append(abs(qv - Fz))
            job.result(f"D = {D}, x = {args.x}, even c = {c}: cᵀQc = {mp.nstr(qv, 12)};  2Σ_ρ F(γ)² = {mp.nstr(Fz, 12)};  rel diff {mp.nstr(abs(qv - Fz) / abs(qv), 3)}")
    job.done()
    if args.gate_abs is not None:
        worst = max(absdiffs)
        ok = worst <= mp.mpf(args.gate_abs)
        print(f"GATE {'PASS' if ok else 'FAIL'}: {len(absdiffs)} tests, max absolute error {mp.nstr(worst, 3)} (limit {args.gate_abs})", flush=True)
        sys.exit(0 if ok else 3)


def cmd_scan(args):
    """λ_min at x = (x/q)·q on a grid in x/q, so that every conductor is compared at the same predicted rate."""
    Ds = [s if s == "dh" else int(s) for s in args.d.split(",")]
    xqs = [mp.mpf(s) for s in args.xq.split(",")]
    rows = []
    job = Job("decay law scan", total=len(Ds) * len(xqs) * 2, args=vars(args))
    for D in Ds:
        q = 5 if D == "dh" else abs(D)
        for xq in xqs:
            x = xq * q
            L = mp.log(x)
            kmax = float(x * L / q)            # k at the crossover frequency T* = 2πx/q
            mp.mp.dps = int(40 + 1.5 * 4 * float(mp.pi) * float(xq) / 2.303)
            for parity in ["even", "odd"]:
                lam = []
                for f in (args.f1, args.f2):
                    n = max(args.nmin, int(f * kmax) + 8)
                    t0 = time.time()
                    M = zeros_side(D, x, n, parity)
                    lm = lam_min_acb(M) if args.acb else lam_min_inv(M) if args.inv else min(mp.eigsy(M, eigvals_only=True))
                    lam.append((n, lm, time.time() - t0))
                (n1, l1, _), (n2, l2, s2) = lam
                job.step(f"D={D} x/q={mp.nstr(xq, 4)} x={mp.nstr(x, 6)} {parity} N={n1}:{mp.nstr(l1, 8)} N={n2}:{mp.nstr(l2, 8)} ({s2:.0f}s, dps {mp.mp.dps})")
                rows.append({"D": D, "q": q, "xq": mp.nstr(xq, 10), "x": mp.nstr(x, 12), "parity": parity, "n1": n1,
                             "lambda_n1": mp.nstr(l1, 15), "n2": n2, "lambda_n2": mp.nstr(l2, 15), "dps": mp.mp.dps})
                if args.json:
                    with open(args.json, "w") as fh:
                        json.dump(rows, fh, indent=1)
    job.done()


def stable_lams(build, dps0, step=40, rtol=mp.mpf("1e-8"), maxdps=600, method=None):
    """Eigenvalue minima that are stable under a precision increase. `build()` returns {label: matrix} at the
    current mp.dps. Each round rebuilds at dps + step; a label is accepted once its minima at two consecutive
    precisions agree to rtol (relative), and the lower-precision value is kept. Returns ({label: λ}, {label: dps})."""
    dps, prev, acc, at = dps0, {}, {}, {}
    while True:
        mp.mp.dps = dps
        mats = build()
        for k, M in mats.items():
            if k in acc:
                continue
            v = (method or lam_min_acb)(M)
            if k in prev and abs(v - prev[k]) <= rtol * abs(v):
                acc[k], at[k] = prev[k], dps - step
            else:
                prev[k] = v
        if len(acc) == len(mats):
            return acc, at
        if dps + step > maxdps:
            raise RuntimeError(f"eigenvalues not stable by {maxdps} digits: {sorted(set(mats) - set(acc))}")
        dps += step


def cmd_tower(args):
    """Lit-review Test 1: the quadratic tower ζ ⊂ ζ_K = ζ·L(χ_{−20}) ⊂ ζ_H = ζ_K·L(χ_{−4})L(χ_5), H = ℚ(i, √5).
    Zeros-side forms add over Euler factors, so every member is a sum of the four degree-1 blocks. Every saved
    eigenvalue is precision-stable (`stable_lams`): ζ alone needs ~130+ digits at a = 1.495 (λ ~ 1e-74 at N = 48)."""
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from control_normalisation import rh2_matrix
    rows = []
    job = Job("tower test", args=vars(args))
    for a_s in args.a.split(","):
        for parity in ["even", "odd"]:
            for n in [int(v) for v in args.n.split(",")]:
                def build():
                    x = mp.exp(2 * mp.mpf(a_s))
                    Z = {D: zeros_side(D, x, n, parity) for D in (1, -20, -4, 5)}
                    return {"zeta": Z[1], "L(chi-20)": Z[-20], "zeta_K": Z[1] + Z[-20], "L(chi-4)": Z[-4], "L(chi5)": Z[5],
                            "L_K(psi)": Z[-4] + Z[5], "zeta_H": Z[1] + Z[-20] + Z[-4] + Z[5]}
                lam, at = stable_lams(build, args.dps)
                mp.mp.dps = max(at.values()) + 40
                x = mp.exp(2 * mp.mpf(a_s))
                ref = rh2_matrix("zetak", x, n, parity)
                zk = zeros_side(1, x, n, parity) + zeros_side(-20, x, n, parity)
                d = max(abs(zk[i, j] - ref[i, j]) for i in range(ref.rows) for j in range(ref.cols))
                r1, r2 = lam["zeta_K"] / lam["L(chi-20)"], lam["zeta_H"] / lam["zeta_K"]
                job.result(f"a={a_s} x={mp.nstr(x, 6)} {parity} N={n}: |ζ_K − control| = {mp.nstr(d, 2)} (dps {mp.mp.dps});  " +
                           "  ".join(f"{k}={mp.nstr(v, 6)}@{at[k]}" for k, v in lam.items()) +
                           f";  ζ_K/L(χ−20) = {mp.nstr(r1, 5)}  ζ_H/ζ_K = {mp.nstr(r2, 5)}")
                rows.append({"a": a_s, "x": mp.nstr(x, 12), "parity": parity, "n": n, "max_diff_vs_control": mp.nstr(d, 3),
                             "control_dps": mp.mp.dps, **{k: mp.nstr(v, 15) for k, v in lam.items()},
                             "stable_dps": at, "ratio_zetaK_over_chi20": mp.nstr(r1, 8), "ratio_zetaH_over_zetaK": mp.nstr(r2, 8)})
                if args.json:
                    with open(args.json, "w") as fh:
                        json.dump(rows, fh, indent=1)
    job.done()


def cmd_precision(args):
    """Recompute saved scan minima (the n2 value) at dps + extra with the same eigen-method; report the relative change."""
    import glob
    rows_out = []
    job = Job("precision recheck", args=vars(args))
    for f in sorted(glob.glob(args.glob)):
        for r in json.load(open(f)):
            if float(r["xq"]) < args.xq_min:
                continue
            D = r["D"] if r["D"] == "dh" else int(r["D"])
            mp.mp.dps = r["dps"] + args.extra
            M = zeros_side(D, mp.mpf(r["x"]), r["n2"], r["parity"])
            v = lam_min_acb(M) if args.acb else lam_min_inv(M)
            old = mp.mpf(r["lambda_n2"])
            rel = abs(v - old) / abs(v)
            job.result(f"D={D} x/q={r['xq']} {r['parity']} N={r['n2']}: dps {r['dps']} → {mp.mp.dps}: {mp.nstr(old, 12)} → {mp.nstr(v, 12)}  rel change {mp.nstr(rel, 3)}")
            rows_out.append({"D": r["D"], "xq": r["xq"], "parity": r["parity"], "n2": r["n2"], "dps": r["dps"], "dps_recheck": mp.mp.dps,
                             "lambda_saved": r["lambda_n2"], "lambda_recheck": mp.nstr(v, 15), "rel_change": mp.nstr(rel, 4)})
            if args.json:
                with open(args.json, "w") as fh:
                    json.dump(rows_out, fh, indent=1)
    job.done()


def _ell(n, c):
    """Fuchs–Slepian leakage ½·4√π 8ⁿ c^{n+½} e^{−2c}/n! (the ½ from 1 − χ ≈ (1 − λ)/2)."""
    return 2 * mp.sqrt(mp.pi) * mp.mpf(8) ** n * c ** (n + mp.mpf(1) / 2) * mp.exp(-2 * c) / mp.factorial(n)


def cmd_pdl2(args):
    """Evaluate P-DL2 (docs/DECAY_LAW.md) exactly as registered:
    R = λ/ℓ_n(2πx/q) with n = κ + 2s must stay in [0.5, 2]× R(10) at x/q = 12, 14, 16 (kill outside [0.25, 4]);
    the assigned n must stay the flattest (metric |ln R(16)/R(3)|) in all 12 cases (kill if > 2 differ)."""
    mp.mp.dps = 30
    assigned = {5: (0, 2), 8: (0, 2), -3: (1, 3), -4: (1, 3), -7: (1, 3), -20: (1, 3)}
    within, killed, flat_bad, out = True, False, 0, []
    for D, ns in assigned.items():
        rows = json.load(open(f"data/connes/decay/scanhi_D{D}.json")) + json.load(open(f"data/connes/decay/scanext_D{D}.json"))
        for par, n in zip(("even", "odd"), ns):
            lam = {float(r["xq"]): mp.mpf(r["lambda_n2"]) for r in rows if r["parity"] == par}
            R = lambda m, u: lam[u] / _ell(m, 2 * mp.pi * u)
            rel = {u: R(n, u) / R(n, 10.0) for u in (12.0, 14.0, 16.0)}
            ok = all(mp.mpf("0.5") <= v <= 2 for v in rel.values())
            kill = any(not (mp.mpf("0.25") <= v <= 4) for v in rel.values())
            flat = min(range(8), key=lambda m: abs(mp.log(R(m, 16.0) / R(m, 3.0))))
            within &= ok
            killed |= kill
            flat_bad += flat != n
            print(f"D={D:>4} {par:<4} n={n}: R(10) = {mp.nstr(R(n, 10.0), 4)};  R/R(10) at 12, 14, 16 = "
                  f"{', '.join(mp.nstr(v, 4) for v in rel.values())}  [{'in' if ok else 'OUT'} 0.5–2]  flattest(3 vs 16) = {flat}"
                  f"{'' if flat == n else '  ≠ assigned'}")
            out.append({"D": D, "parity": par, "n": n, "R10": mp.nstr(R(n, 10.0), 8),
                        "R_over_R10": {str(u): mp.nstr(v, 8) for u, v in rel.items()}, "flattest_3_16": flat})
    verdict = ("KILLED" if killed or flat_bad > 2 else "HOLDS" if within and flat_bad == 0 else "NOT KILLED, prediction partly missed")
    print(f"P-DL2: R within [0.5, 2]: {within};  any R outside [0.25, 4]: {killed};  flattest ≠ assigned: {flat_bad}/12  →  {verdict}")
    if args.json:
        with open(args.json, "w") as fh:
            json.dump({"rows": out, "within": within, "killed": killed, "flattest_mismatch": flat_bad, "verdict": verdict}, fh, indent=1)


def cmd_product(args):
    """P-DL3 (docs/DECAY_LAW.md): ρ = λ(Σ factor forms)/λ(largest-conductor factor) on a grid in x/q_max.
    Products are PSD sums, so inverse iteration is used, inside the precision-stability loop."""
    Ds = [int(v) for v in args.factors.split(",")]
    qmax = max(abs(D) for D in Ds)
    worst = next(D for D in Ds if abs(D) == qmax)
    rows = []
    job = Job(f"product {args.factors}", total=len(args.xq.split(",")) * 2, args=vars(args))
    for xq_s in args.xq.split(","):
        xq = mp.mpf(xq_s)
        for parity in ["even", "odd"]:
            res = {}
            x = xq * qmax
            ns = [max(args.nmin, int(f * float(x * mp.log(x) / qmax)) + 8) for f in (args.f1, args.f2)]
            ns[1] = max(ns[1], ns[0] + 8)          # the convergence check needs two distinct sizes
            for n in ns:
                def build():
                    Z = {D: zeros_side(D, xq * qmax, n, parity) for D in Ds}
                    tot = Z[Ds[0]]
                    for D in Ds[1:]:
                        tot = tot + Z[D]
                    return {"product": tot, "worst": Z[worst]}
                lam, at = stable_lams(build, int(40 + 1.5 * 4 * float(mp.pi) * float(xq) / 2.303), method=lam_min_inv)
                res[n] = (lam, at)
            (n1, (l1, a1)), (n2, (l2, a2)) = sorted(res.items())
            rho = l2["product"] / l2["worst"]
            job.step(f"x/q_max={xq_s} x={mp.nstr(xq * qmax, 6)} {parity}: N={n2} λ(product)={mp.nstr(l2['product'], 6)} "
                     f"λ(L(χ{worst}))={mp.nstr(l2['worst'], 6)}  ρ={mp.nstr(rho, 5)}  (N={n1}: ρ={mp.nstr(l1['product'] / l1['worst'], 5)}; dps {a2})")
            rows.append({"factors": args.factors, "worst": worst, "xq": xq_s, "x": mp.nstr(xq * qmax, 12), "parity": parity,
                         "n1": n1, "product_n1": mp.nstr(l1["product"], 15), "worst_n1": mp.nstr(l1["worst"], 15),
                         "n2": n2, "product_n2": mp.nstr(l2["product"], 15), "worst_n2": mp.nstr(l2["worst"], 15),
                         "rho": mp.nstr(rho, 10), "stable_dps": a2})
            if args.json:
                with open(args.json, "w") as fh:
                    json.dump(rows, fh, indent=1)
    job.done()


def cmd_product_v(args):
    """P-DL5/P-DL10 (docs/DECAY_LAW.md): the product's own λ on a grid v = (x/Q)^{1/d}, Q = Π|D|, d = --degree
    (default 2), N_basis = f·vL + 8. Only the product's minimum is computed."""
    Ds = [int(v) for v in args.factors.split(",")]
    Q = 1
    for D in Ds:
        Q *= abs(D)
    rows = []
    job = Job(f"product-v {args.factors}", total=len(args.v.split(",")) * 2, args=vars(args))
    for v_s in args.v.split(","):
        v = mp.mpf(v_s)
        x = Q * v ** args.degree
        for parity in ["even", "odd"]:
            ns = [max(args.nmin, int(f * float(v * mp.log(x))) + 8) for f in (args.f1, args.f2)]
            ns[1] = max(ns[1], ns[0] + 8)
            lam = {}
            for n in ns:
                def build():
                    tot = None
                    for D in Ds:
                        Z = zeros_side(D, x, n, parity)
                        tot = Z if tot is None else tot + Z
                    return {"product": tot}
                l, at = stable_lams(build, int(40 + 1.5 * 4 * args.degree * float(mp.pi) * float(v) / 2.303), method=lam_min_inv)
                lam[n] = (l["product"], at["product"])
            (n1, (l1, _)), (n2, (l2, a2)) = sorted(lam.items())
            job.step(f"v={v_s} x={mp.nstr(x, 6)} {parity}: N={n1}:{mp.nstr(l1, 8)}  N={n2}:{mp.nstr(l2, 8)}  (dps {a2})")
            rows.append({"factors": args.factors, "Q": Q, "degree": args.degree, "v": v_s, "x": mp.nstr(x, 12), "parity": parity, "n1": n1,
                         "lambda_n1": mp.nstr(l1, 15), "n2": n2, "lambda_n2": mp.nstr(l2, 15), "stable_dps": a2})
            if args.json:
                with open(args.json, "w") as fh:
                    json.dump(rows, fh, indent=1)
    job.done()


def cmd_index_test(args):
    """P-DL8 (docs/DECAY_LAW.md): for each D, both sectors, the flattest prolate index n ∈ 0..7 under
    |ln R(10)/R(3)| with c = 2πx/q (prediction n = κ + 2s), and the log-corrected rate fit over x/q ∈ [3, 8]."""
    mp.mp.dps = 30
    rows_out, hits, rate_ok, rate_kill = [], 0, True, False
    Ds = [int(v) for v in args.ds.split(",")]
    for D in Ds:
        rows = json.load(open(args.scan_glob.format(D=D)))
        kappa = 0 if D > 0 else 1
        for s_i, par in enumerate(("even", "odd")):
            lam = {float(r["xq"]): mp.mpf(r["lambda_n2"]) for r in rows if r["parity"] == par}
            R = lambda m, u: lam[u] / _ell(m, 2 * mp.pi * u)
            flat = min(range(8), key=lambda m: abs(mp.log(R(m, 10.0) / R(m, 3.0))))
            pred = kappa + 2 * s_i
            hits += flat == pred
            xs = sorted(u for u in lam if 3 <= u <= 8)
            c, res = _lsq([[mp.mpf(u), mp.log(u), 1] for u in xs], [mp.log(lam[u]) for u in xs])
            a = -c[0] / (4 * mp.pi)
            rate_ok &= mp.mpf("0.95") <= a <= mp.mpf("1.05")
            rate_kill |= not (mp.mpf("0.9") <= a <= mp.mpf("1.1"))
            print(f"D={D:>4} {par:<4}: flattest n = {flat} (predicted {pred}){'' if flat == pred else '  ✗'};  R(3) = {mp.nstr(R(pred, 3.0), 4)}, "
                  f"R(10) = {mp.nstr(R(pred, 10.0), 4)};  log-corrected α/4π = {mp.nstr(a, 5)} (γ = {mp.nstr(c[1], 3)}, resid {mp.nstr(res, 2)})")
            rows_out.append({"D": D, "parity": par, "flattest": flat, "predicted": pred, "alpha_over_4pi": mp.nstr(a, 8),
                             "R3": mp.nstr(R(pred, 3.0), 6), "R10": mp.nstr(R(pred, 10.0), 6)})
    n = 2 * len(Ds)
    va = "HOLDS" if hits >= n - 1 else ("KILLED" if hits <= n - 3 else "NOT KILLED, missed")
    vb = "KILLED" if rate_kill else ("HOLDS" if rate_ok else "NOT KILLED, missed")
    print(f"P-DL8a: {hits}/{n} flattest = κ + 2s → {va};  P-DL8b: → {vb}")
    if args.json:
        with open(args.json, "w") as fh:
            json.dump({"rows": rows_out, "index_hits": hits, "verdict_a": va, "verdict_b": vb}, fh, indent=1)


def _lsq(A, y):
    A, y = mp.matrix(A), mp.matrix(y)
    c = mp.lu_solve(A.T * A, A.T * y)
    return c, max(abs(v) for v in A * c - y)


def cmd_fit(args):
    """Per (D, sector): slope of ln λ in x/q (linear, and with a γ ln(x/q) term) over x/q ≥ --from; collapse against ζ."""
    import glob
    mp.mp.dps = 30
    data = {}
    for f in sorted(glob.glob(args.glob)):
        for r in json.load(open(f)):
            data.setdefault((r["D"], r["parity"]), []).append(r)
    four_pi = 4 * mp.pi
    out = []
    for (D, par), R in sorted(data.items()):
        R.sort(key=lambda r: float(r["xq"]))
        xq = [mp.mpf(r["xq"]) for r in R]
        y = [mp.log(mp.mpf(r["lambda_n2"])) for r in R]
        conv = [float(mp.log(mp.mpf(r["lambda_n1"]) / mp.mpf(r["lambda_n2"]))) for r in R]
        sel = [i for i, v in enumerate(xq) if args.from_xq <= v <= args.to_xq]
        c1, r1 = _lsq([[xq[i], 1] for i in sel], [y[i] for i in sel])
        c2, r2 = _lsq([[xq[i], mp.log(xq[i]), 1] for i in sel], [y[i] for i in sel])
        line = (f"D={D:>4} {par:<4} linear slope/4π = {mp.nstr(-c1[0] / four_pi, 4)} (resid {mp.nstr(r1, 2)});  "
                f"with log: slope/4π = {mp.nstr(-c2[0] / four_pi, 4)}, γ = {mp.nstr(c2[1], 3)}, β = {mp.nstr(c2[2], 4)} (resid {mp.nstr(r2, 2)});  "
                f"ln λ+4π·x/q = {[round(float(v + four_pi * u), 2) for u, v in zip(xq, y)]};  conv ln(λ1/λ2) = {[round(v, 2) for v in conv]}")
        print(line)
        out.append({"D": D, "parity": par, "slope_linear_over_4pi": mp.nstr(-c1[0] / four_pi, 6),
                    "slope_log_over_4pi": mp.nstr(-c2[0] / four_pi, 6), "gamma": mp.nstr(c2[1], 6), "beta": mp.nstr(c2[2], 6)})
    if args.json:
        with open(args.json, "w") as fh:
            json.dump(out, fh, indent=1)


ZETA_SOURCES = """Best committed upper bound λ_N per (a, sector) for ζ, used by `fit-zeta`:
  a = 0.8    even: cosine N = 800, data/connes/zhu_window_ext.json
             odd:  Legendre 200 modes, data/connes/zhu_legendre_upper_odd.json (below sine N = 800)
  a = 1.0, 1.19, 1.3: min over periodic and antiperiodic bases at N = 140, data/connes/antiperiodic/sweep_a*.json
  a = 1.495  N = 180 (both sectors), data/connes/zhu_window_a1495_n180.run.txt (zhu_window_mp.py log, 260 digits)"""


def zeta_points():
    """{parity: [(x, λ, source)]} from the committed files listed in ZETA_SOURCES."""
    import re
    pts = {"even": [], "odd": []}
    x08 = mp.exp(mp.mpf("1.6"))
    ext = json.load(open("data/connes/zhu_window_ext.json"))
    pts["even"].append((x08, mp.mpf(next(r for r in ext if r["a"] == "0.8" and r["n"] == 800 and r["parity"] == "even")["eps"][0]), "cosine N=800"))
    leg = json.load(open("data/connes/zhu_legendre_upper_odd.json"))
    pts["odd"].append((x08, mp.mpf(leg["blocks"][-1]["lambda_U"]), "Legendre 200"))
    for a in ("1.0", "1.19", "1.3"):
        rows = json.load(open(f"data/connes/antiperiodic/sweep_a{a}.json"))
        for par in ("even", "odd"):
            r = min((r for r in rows if r["parity"] == par and r["n"] == 140), key=lambda r: mp.mpf(r["lambda_min"]))
            pts[par].append((mp.exp(2 * mp.mpf(a)), mp.mpf(r["lambda_min"]), f"{r['basis']} N=140"))
    log = open("data/connes/zhu_window_a1495_n180.run.txt").read()
    for par in ("even", "odd"):
        m = re.search(rf"RESULT N=180 {par}: eps0=([0-9.e+-]+)", log)
        pts[par].append((mp.exp(2 * mp.mpf("1.495")), mp.mpf(m.group(1)), "N=180"))
    return pts


def cmd_fit_zeta(args):
    """ln λ_ζ(x) = −αx + γ ln x + β over the committed best upper bounds (ZETA_SOURCES)."""
    mp.mp.dps = 30
    print(ZETA_SOURCES)
    pts = zeta_points()
    out = {}
    for par in ("even", "odd"):
        P = pts[par]
        for x, lam, src in P:
            print(f"  {par:<4} x = {mp.nstr(x, 8):>10}  λ = {mp.nstr(lam, 10):>18}  ({src})")
        c, r = _lsq([[-x, mp.log(x), 1] for x, _, _ in P], [mp.log(l) for _, l, _ in P])
        c2, r2 = _lsq([[mp.log(x), 1] for x, _, _ in P], [mp.log(l) + 4 * mp.pi * x for x, l, _ in P])
        print(f"{par}: α = {mp.nstr(c[0], 6)} (α/4π = {mp.nstr(c[0] / (4 * mp.pi), 6)}), γ = {mp.nstr(c[1], 4)}, β = {mp.nstr(c[2], 5)}, max resid {mp.nstr(r, 3)};"
              f"  with α := 4π: γ = {mp.nstr(c2[0], 4)}, β = {mp.nstr(c2[1], 5)}, max resid {mp.nstr(r2, 3)}")
        out[par] = {"alpha": mp.nstr(c[0], 10), "alpha_over_4pi": mp.nstr(c[0] / (4 * mp.pi), 10), "gamma": mp.nstr(c[1], 8),
                    "beta": mp.nstr(c[2], 8), "max_resid": mp.nstr(r, 4), "points": [[mp.nstr(x, 15), mp.nstr(l, 15), s] for x, l, s in P]}
    lr = [mp.log(o[1] / e[1]) for e, o in zip(pts["even"], pts["odd"])]
    c, r = _lsq([[mp.log(e[0]), 1] for e in pts["even"]], lr)
    print(f"ln(λ_odd/λ_even) = {mp.nstr(c[0], 4)} ln x + {mp.nstr(c[1], 4)} (max resid {mp.nstr(r, 3)})")
    out["odd_even_power"] = mp.nstr(c[0], 8)
    if args.json:
        with open(args.json, "w") as fh:
            json.dump(out, fh, indent=1)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    c = sub.add_parser("check")
    c.add_argument("--x", default="13")
    c.add_argument("--tmax", type=float, default=120)
    c.add_argument("--d", default="5,-4,-20")
    c.add_argument("--gate-abs", default=None, help="exit 3 unless every absolute error ≤ this (P-DL8 gate)")
    s = sub.add_parser("scan")
    s.add_argument("--d", default="1,-3,-4,5,-7,8,-20")
    s.add_argument("--xq", default="2,3,4,5,6,8")
    s.add_argument("--f1", type=float, default=1.5)
    s.add_argument("--f2", type=float, default=2.5)
    s.add_argument("--nmin", type=int, default=24)
    s.add_argument("--inv", action="store_true", help="inverse iteration (FLINT) instead of mp.eigsy")
    s.add_argument("--acb", action="store_true", help="FLINT full eigensolver (needed when λ_min may be negative)")
    s.add_argument("--json")
    t = sub.add_parser("tower")
    t.add_argument("--a", default="1.1,1.3,1.495")
    t.add_argument("--n", default="48,96")
    t.add_argument("--dps", type=int, default=60)
    t.add_argument("--json")
    pr = sub.add_parser("precision")
    pr.add_argument("--glob", required=True)
    pr.add_argument("--xq-min", type=float, default=10)
    pr.add_argument("--extra", type=int, default=40)
    pr.add_argument("--acb", action="store_true")
    pr.add_argument("--json")
    pp = sub.add_parser("product")
    pp.add_argument("--factors", required=True, help="comma-separated D values, e.g. 1,-20")
    pp.add_argument("--xq", default="2,4,6,8", help="x/q_max grid")
    pp.add_argument("--f1", type=float, default=5)
    pp.add_argument("--f2", type=float, default=9)
    pp.add_argument("--nmin", type=int, default=24)
    pp.add_argument("--json")
    pv = sub.add_parser("product-v")
    pv.add_argument("--factors", required=True)
    pv.add_argument("--v", default="1.5,2,2.5,3,3.5")
    pv.add_argument("--f1", type=float, default=5)
    pv.add_argument("--f2", type=float, default=9)
    pv.add_argument("--nmin", type=int, default=24)
    pv.add_argument("--degree", type=int, default=2, help="d in v = (x/Q)^{1/d}")
    pv.add_argument("--json")
    it = sub.add_parser("index-test")
    it.add_argument("--ds", required=True)
    it.add_argument("--scan-glob", default="data/connes/decay/scan8_D{D}.json")
    it.add_argument("--json")
    p2 = sub.add_parser("pdl2")
    p2.add_argument("--json")
    fz = sub.add_parser("fit-zeta")
    fz.add_argument("--json")
    f = sub.add_parser("fit")
    f.add_argument("--glob", default="data/connes/decay/scanhi_D*.json")
    f.add_argument("--from-xq", type=float, default=3)
    f.add_argument("--to-xq", type=float, default=1e9)
    f.add_argument("--json")
    args = ap.parse_args()
    {"check": cmd_check, "scan": cmd_scan, "fit": cmd_fit, "tower": cmd_tower, "fit-zeta": cmd_fit_zeta, "precision": cmd_precision, "pdl2": cmd_pdl2, "product": cmd_product, "product-v": cmd_product_v, "index-test": cmd_index_test}[args.cmd](args)


if __name__ == "__main__":
    main()
