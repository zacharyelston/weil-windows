#!/usr/bin/env python3
"""Discriminant −20: an Euler product and a non-Euler-product function with the
same functional equation (degree 2, conductor 20).

  ζ_K(s) = ζ(s) L(s, χ₋₂₀)                       Dedekind zeta of Q(√−5): Euler product
  G(s)   = L(s, χ₋₄) L(s, χ₅)                     genus-character product: Euler product
  Z₁(s)  = ½[ζ_K(s) + G(s)] = Σ r₁(n)/2 · n^{−s}  Epstein zeta of x² + 5y²: no Euler product
  Z₂(s)  = ½[ζ_K(s) − G(s)]                       Epstein zeta of 2x² + 2xy + 3y²

All four satisfy Λ(s) = 20^{s/2} (2π)^{−s} Γ(s) F(s) = Λ(1 − s).

Step 1 (`check`): functional equations, coefficients against representation
counts. Step 2 (`census`): zeros of Z₁ on and off the critical line.

Usage:
  .venv/bin/python scripts/epstein_mp.py check
  .venv/bin/python scripts/epstein_mp.py census --t1 1 --t2 120
"""

import argparse
import json
import math
import sys

import mpmath as mp

from progress import Job

CHI_M4 = [0, 1, 0, -1]           # (−4/n)
CHI_5 = [0, 1, -1, -1, 1]        # (n/5)
CHI_M20 = [CHI_M4[n % 4] * CHI_5[n % 5] for n in range(20)]  # (−20/n)


def L(s, chi):
    return mp.dirichlet(s, chi)


def zeta_K(s):
    return mp.zeta(s) * L(s, CHI_M20)


def genus(s):
    return L(s, CHI_M4) * L(s, CHI_5)


def Z1(s):
    return (zeta_K(s) + genus(s)) / 2


def Z2(s):
    return (zeta_K(s) - genus(s)) / 2


FUNCS = {"zeta_K": zeta_K, "genus": genus, "Z1": Z1, "Z2": Z2}


def gamma_factor_log(s):
    return s / 2 * mp.log(20) - s * mp.log(2 * mp.pi) + mp.loggamma(s)


def completed(f, s):
    return mp.exp(gamma_factor_log(s)) * f(s)


def theta(t):
    return mp.im(gamma_factor_log(mp.mpc(0.5, t)))


def hardy(f, t):
    v = mp.exp(1j * theta(t)) * f(mp.mpc(0.5, t))
    return v


def coefficients(n_max):
    """a_n for Z₁ from representation counts r₁(n)/2, and from the L-function convolution."""
    r = [0] * (n_max + 1)
    lim = int(math.isqrt(n_max)) + 1
    for x in range(-lim, lim + 1):
        for y in range(-lim, lim + 1):
            v = x * x + 5 * y * y
            if 0 < v <= n_max:
                r[v] += 1
    direct = [mp.mpf(r[n]) / 2 for n in range(n_max + 1)]
    conv = [mp.mpf(0)] * (n_max + 1)
    for n in range(1, n_max + 1):
        dk = sum(CHI_M20[d % 20] for d in range(1, n + 1) if n % d == 0)
        gn = sum(CHI_M4[d % 4] * CHI_5[(n // d) % 5] for d in range(1, n + 1) if n % d == 0)
        conv[n] = mp.mpf(dk + gn) / 2
    return direct, conv


def zero_count_rect(f, s0, s1, t0, t1):
    corners = [mp.mpc(s0, t0), mp.mpc(s1, t0), mp.mpc(s1, t1), mp.mpc(s0, t1), mp.mpc(s0, t0)]
    total = mp.mpf(0)
    for a, b in zip(corners[:-1], corners[1:]):
        length = abs(b - a)
        u, h = mp.mpf(0), min(mp.mpf("0.05"), length)
        fa = f(a)
        while u < length - mp.mpf("1e-12"):
            step = min(h, length - u)
            fb = f(a + (b - a) * ((u + step) / length))
            d = mp.arg(fb / fa)
            if abs(d) > mp.pi / 8 and step > mp.mpf("1e-6"):
                h = step / 2
                continue
            total += d
            u += step
            fa = fb
            if abs(d) < mp.pi / 32:
                h = min(h * 1.5, mp.mpf("0.25"))
    return int(mp.nint(total / (2 * mp.pi)))


def newton(f, s, iters=60):
    start = s
    for _ in range(iters):
        h = mp.mpf("1e-8")
        fv = f(s)
        df = (f(s + h) - f(s - h)) / (2 * h)
        step = fv / df
        if abs(step) > 1:
            step /= abs(step)
        s -= step
        if abs(s - start) > 5:
            return s, mp.inf
        if abs(step) < mp.mpf(10) ** (-mp.mp.dps + 4):
            break
    return s, abs(f(s))


def off_line_zeros(f, t1, t2, box=2, gap=mp.mpf("0.002"), s_max=3, job=None):
    out = []

    def search(s0, s1, a, b, depth):
        n = zero_count_rect(f, s0, s1, a, b)
        if n == 0:
            return
        if n == 1:
            for fs, ft in [(0.5, 0.5), (0.25, 0.25), (0.75, 0.25), (0.25, 0.75), (0.75, 0.75)]:
                z, r = newton(f, mp.mpc(s0 + fs * (s1 - s0), a + ft * (b - a)))
                if s0 - 1e-9 <= z.real <= s1 + 1e-9 and a - 1e-9 <= z.imag <= b + 1e-9 and r < 1e-8:
                    out.append((z, r))
                    return
        if depth > 30:
            print(f"warning: unresolved box [{s0},{s1}]x[{a},{b}] n={n}", file=sys.stderr)
            return
        if (s1 - s0) > (b - a):
            m = (s0 + s1) / 2
            search(s0, m, a, b, depth + 1)
            search(m, s1, a, b, depth + 1)
        else:
            m = (a + b) / 2
            search(s0, s1, a, m, depth + 1)
            search(s0, s1, m, b, depth + 1)

    t = mp.mpf(t1)
    n_boxes = int(mp.ceil((t2 - t1) / box))
    k = 0
    while t < t2:
        search(mp.mpf(0.5) + gap, mp.mpf(s_max), t, min(t + box, mp.mpf(t2)), 0)
        t += box
        k += 1
        if job:
            job.sub(f"off-line boxes up to γ={mp.nstr(t, 5)} ({len(out)} found)", k, n_boxes, every=max(1, n_boxes // 20))
    return sorted(out, key=lambda p: p[0].imag)


def on_line_zeros(f, t1, t2, step=0.02, job=None):
    z = lambda t: mp.re(hardy(f, t))
    out = []
    a, za = mp.mpf(t1), z(t1)
    n_steps = int((t2 - t1) / step) + 1
    k = 0
    while a < t2:
        k += 1
        if job:
            job.sub(f"on-line scan at γ={mp.nstr(a, 5)} ({len(out)} zeros)", k, n_steps, every=max(1, n_steps // 20))
        b = a + step
        zb = z(b)
        if za * zb < 0:
            out.append(mp.findroot(z, (a, b), solver="bisect", tol=mp.mpf(10) ** (-mp.mp.dps + 4)))
        a, za = b, zb
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("check")
    c = sub.add_parser("census")
    c.add_argument("--t1", type=float, default=1)
    c.add_argument("--t2", type=float, default=120)
    c.add_argument("--function", default="Z1", choices=list(FUNCS))
    c.add_argument("--dps", type=int, default=20)
    c.add_argument("--json")
    args = ap.parse_args()

    if args.cmd == "check":
        mp.mp.dps = 30
        print("functional equation |Λ(s) − Λ(1−s)| / |Λ(s)|:")
        for name, f in FUNCS.items():
            errs = []
            for s in [mp.mpc(0.3, 7), mp.mpc(0.8, 31), mp.mpc(-0.4, 55)]:
                a, b = completed(f, s), completed(f, 1 - s)
                errs.append(abs(a - b) / abs(a))
            print(f"  {name:<7} " + "  ".join(mp.nstr(e, 3) for e in errs))
        print("Hardy function Im/|Re| on the line (should vanish):")
        for name in ["zeta_K", "Z1"]:
            vals = [hardy(FUNCS[name], t) for t in [10, 33.3, 77.7]]
            print(f"  {name:<7} " + "  ".join(mp.nstr(abs(mp.im(v)) / max(abs(mp.re(v)), mp.mpf(1e-30)), 3) for v in vals))
        direct, conv = coefficients(200)
        bad = [n for n in range(1, 201) if direct[n] != conv[n]]
        print(f"coefficients of Z₁: r₁(n)/2 vs ½(Σχ₋₂₀(d) + Σχ₋₄(d)χ₅(n/d)) for n ≤ 200: {'all agree' if not bad else 'mismatch at ' + str(bad[:10])}")
        print("  a_n, n = 1..30: " + " ".join(str(int(direct[n])) for n in range(1, 31)))
        sv = mp.mpc(2, 3)
        dsum = mp.fsum(direct[n] * mp.power(n, -sv) for n in range(1, 201))
        print(f"  Z₁(2+3i) = {mp.nstr(Z1(sv), 10)}; partial Dirichlet sum to 200 = {mp.nstr(dsum, 10)}")
        return

    mp.mp.dps = args.dps
    f = FUNCS[args.function]
    job = Job(f"epstein census {args.function}", total=3, args=vars(args))
    on = on_line_zeros(f, args.t1, args.t2, job=job)
    job.step(f"on-line zeros: {len(on)}")
    off = off_line_zeros(f, args.t1, args.t2, job=job)
    job.step(f"off-line pairs: {len(off)}")
    total = zero_count_rect(f, -2, 3, args.t1, args.t2)
    job.step(f"argument principle total: {total}")
    job.done()
    print(f"{args.function}: zeros with γ in [{args.t1}, {args.t2}]")
    print(f"  argument principle on [−2, 3]: {total}")
    print(f"  on the line: {len(on)};  off the line: {len(off)} pairs → accounted {len(on) + 2 * len(off)} of {total}")
    for z, r in off:
        print(f"    β = {mp.nstr(z.real, 8)}  γ = {mp.nstr(z.imag, 10)}   |f| = {mp.nstr(r, 2)}")
    if args.json:
        with open(args.json, "w") as fh:
            json.dump({"function": args.function, "t1": args.t1, "t2": args.t2, "total": total,
                       "on_line": [float(z) for z in on],
                       "off_line": [{"beta": float(z.real), "gamma": float(z.imag)} for z, _ in off]}, fh, indent=2)


if __name__ == "__main__":
    main()
