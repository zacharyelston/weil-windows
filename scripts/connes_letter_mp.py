#!/usr/bin/env python3
"""High-precision reproduction of Connes' "Letter to Riemann" (arXiv:2602.04022, §5).

Weil's quadratic form QW_λ on even test functions supported in [λ⁻¹, λ]
(additive variable u ∈ [−L/2, L/2], L = 2 log λ), restricted to Weil's subspace
φ̂(±i/2) = 0. The minimiser η of QW_λ has a Fourier transform whose real zeros
approximate the zeta zeros. Same construction as src/rh/connes_letter.rs, in
mpmath at `--dps` digits.

Exact pieces (ω_k = 2πk/L, ω_k L = 2πk, K(w) = e^{w/2}/(e^w − e^{−w}) = Σ_m e^{−a_m w},
a_m = 2m + 1/2):
  I_k = ∫_0^L sin(ω_k w) K(w) dw = Σ_m ω (1 − e^{−aL}) / (a² + ω²)
  J_k = ∫_0^L [2(L − w) cos(ω_k w) − 2L e^{−w/2}] K(w) dw
      = Σ_m [2 C(a) − 2L (1 − e^{−(a+1/2)L})/(a + 1/2)],
  C(a) = ∫_0^L (L − w) cos(ωw) e^{−aw} dw
       = L a (1−E)/(a²+ω²) + [((1−E) + aLE)(a²+ω²) − 2a²(1−E)]/(a²+ω²)²,  E = e^{−aL}.

Usage: .venv/bin/python scripts/connes_letter_mp.py --lambda-sq 13 --n 100 --dps 120 --zeros 50
"""

import argparse
import json
import time

import mpmath as mp

CHECK_SERIES = False


def prime_powers(n_max):
    """(n, log p) for prime powers 2 ≤ n ≤ n_max."""
    out = []
    for p in range(2, n_max + 1):
        if all(p % q for q in range(2, int(p ** 0.5) + 1)):
            pk = p
            while pk <= n_max:
                out.append((pk, mp.log(p)))
                pk *= p
    return out


def dh_log_derivative(n_max):
    """c_n of −f'/f for the Davenport–Heilbronn function (a_1 = 1)."""
    r5 = mp.sqrt(5)
    kappa = (mp.sqrt(10 - 2 * r5) - 2) / (r5 - 1)
    chi_re_im = {1: (1, 0), 2: (0, 1), 3: (0, -1), 4: (-1, 0), 0: (0, 0)}

    def a(n):
        re, im = chi_re_im[n % 5]
        return re + kappa * im  # Re((1 − iκ)(re + i·im))

    c = [mp.mpf(0)] * (n_max + 1)
    av = [mp.mpf(0)] + [a(k) for k in range(1, n_max + 1)]
    for k in range(2, n_max + 1):
        c[k] = av[k] * mp.log(k)
    for d in range(2, n_max + 1):
        if c[d] == 0:
            continue
        for m in range(2 * d, n_max + 1, d):
            c[m] -= c[d] * av[m // d]
    return c


def build_form(lambda_sq, n, function="zeta", terms=None, include_arch=True, parity="even", progress_job=None):
    """QW_λ on the even cosine basis. `terms` overrides the prime-side list
    [(log n, coefficient)], and include_arch=False drops the archimedean part
    (ψ, log 4π + γ, conductor and sech terms), so that the form splits as
    Q = Q_arch + Σ_n Q_n. parity="odd" returns the odd sector in the basis
    √(2/L) sin(ω_k u), k = 1..n (an n × n matrix)."""
    L = mp.log(lambda_sq)
    om = lambda k: 2 * mp.pi * k / L

    quarter = mp.mpf(1) / 4
    # Tail terms carry E = e^{−aL}, a = 2m + 1/2: geometric, summed directly.
    tail_terms = int(mp.mp.dps * mp.log(10) / (2 * L)) + 10

    def i_k_series(k):
        w = om(k)
        if k == 0:
            return mp.mpf(0)
        return mp.nsum(lambda m: w * (1 - mp.exp(-(2 * m + mp.mpf(1) / 2) * L)) / ((2 * m + mp.mpf(1) / 2) ** 2 + w * w), [0, mp.inf])

    def j_term(m, w):
        a = 2 * m + mp.mpf(1) / 2
        E = mp.exp(-a * L)
        d = a * a + w * w
        c = L * a * (1 - E) / d + (((1 - E) + a * L * E) * d - 2 * a * a * (1 - E)) / (d * d)
        return 2 * c - 2 * L * (1 - mp.exp(-(a + mp.mpf(1) / 2) * L)) / (a + mp.mpf(1) / 2)

    def j_k_series(k):
        w = om(k)
        return mp.nsum(lambda m: j_term(m, w), [0, mp.inf])

    def i_k(k):
        # Σ_m ω/(a² + ω²) = Im ψ(1/4 + iω/2)/2, minus the E-tail.
        w = om(k)
        if k == 0:
            return mp.mpf(0)
        z = mp.mpc(quarter, w / 2)
        main = mp.im(mp.digamma(z)) / 2
        tail = mp.fsum(w * mp.exp(-(2 * m + mp.mpf(1) / 2) * L) / ((2 * m + mp.mpf(1) / 2) ** 2 + w * w) for m in range(tail_terms))
        return main - tail

    def j_k(k):
        # E-free part: L(ψ(1/2) − Re ψ(1/4 + iω/2)) − ½ Re ψ'(1/4 + iω/2); the rest decays like e^{−2mL}.
        w = om(k)
        z = mp.mpc(quarter, w / 2)
        main = L * (mp.digamma(mp.mpf(1) / 2) - mp.re(mp.digamma(z))) - mp.re(mp.psi(1, z)) / 2

        def t0(m):
            a = 2 * m + mp.mpf(1) / 2
            d = a * a + w * w
            return 2 * L * a / d + 2 * (w * w - a * a) / (d * d) - 2 * L / (a + mp.mpf(1) / 2)

        tail = mp.fsum(j_term(m, w) - t0(m) for m in range(tail_terms))
        return main + tail

    if terms is not None:
        pp = list(terms)
    elif function == "dh":
        cn = dh_log_derivative(int(lambda_sq))
        pp = [(mp.log(k), cn[k] / mp.sqrt(k)) for k in range(2, int(lambda_sq) + 1) if cn[k] != 0]
    else:
        pp = [(mp.log(nn), lp / mp.sqrt(nn)) for nn, lp in prime_powers(int(lambda_sq))]

    # Odd conductor-5 archimedean correction (DH): log 5 · F(0) + ½∫ sech(u/2)(F(u)+F(−u)),
    # sech(u/2) = 2 Σ_m (−1)^m e^{−(m+1/2)u}.
    def sech_sin(k):
        w = om(k)
        if k == 0:
            return mp.mpf(0)
        return 2 * mp.nsum(lambda m: (-1) ** int(m) * w * (1 - mp.exp(-(m + mp.mpf(1) / 2) * L)) / ((m + mp.mpf(1) / 2) ** 2 + w * w), [0, mp.inf])

    def sech_lcos(k):
        w = om(k)

        def term(m):
            b = m + mp.mpf(1) / 2
            E = mp.exp(-b * L)
            d = b * b + w * w
            return (-1) ** int(m) * (L * b * (1 - E) / d + (((1 - E) + b * L * E) * d - 2 * b * b * (1 - E)) / (d * d))

        return 2 * mp.nsum(term, [0, mp.inf])

    odd = function == "dh" and include_arch
    zero = [mp.mpf(0)] * (n + 1)

    def tracked(name, fn):
        vals = []
        for k in range(n + 1):
            vals.append(fn(k))
            if progress_job is not None:
                progress_job.sub(f"build_form {name}", k + 1, n + 1, every=max(1, (n + 1) // 10))
        return vals

    IS = tracked("sech·sin series", sech_sin) if odd else zero
    SD = tracked("sech·(L−u)cos series", sech_lcos) if odd else zero
    log_q = mp.log(5) if odd else mp.mpf(0)
    I = tracked("I_k", i_k) if include_arch else zero
    J = tracked("J_k", j_k) if include_arch else zero
    if CHECK_SERIES and include_arch:
        di = max(abs(I[k] - i_k_series(k)) for k in range(n + 1))
        dj = max(abs(J[k] - j_k_series(k)) for k in range(n + 1))
        print(f"closed form vs nsum: max |ΔI| = {mp.nstr(di, 3)}, max |ΔJ| = {mp.nstr(dj, 3)}")
    P = [mp.fsum(a * mp.sin(om(k) * ln) for ln, a in pp) for k in range(n + 1)]
    D = [mp.fsum(a * 2 * (L - ln) * mp.cos(om(k) * ln) for ln, a in pp) for k in range(n + 1)]
    diag_const = mp.log(4 * mp.pi) + mp.euler + mp.log(mp.tanh(L / 2)) if include_arch else mp.mpf(0)

    def get(v, k, odd):
        return v[k] if k >= 0 else (-v[-k] if odd else v[-k])

    def q(j, k):
        if j == k:
            return (-(diag_const * L + get(J, k, False)) - get(D, k, False) + log_q * L + get(SD, k, False)) / L
        sign = 1 if (j - k) % 2 == 0 else -1
        delta = om(j) - om(k)
        off = -2 * sign * ((get(I, k, True) - get(I, j, True)) + (get(P, k, True) - get(P, j, True)))
        off += sign * (get(IS, k, True) - get(IS, j, True))
        return off / delta / L

    if parity == "odd":
        # b_k = (e_k − e_{−k})/(i√2) = √(2/L) sin(ω_k u), k = 1..n.
        E = mp.matrix(n, n)
        for a in range(1, n + 1):
            for b in range(a, n + 1):
                v = (q(a, b) - q(a, -b) - q(-a, b) + q(-a, -b)) / 2
                E[a - 1, b - 1] = v
                E[b - 1, a - 1] = v
        return L, E, pp

    s2 = mp.sqrt(2)
    E = mp.matrix(n + 1, n + 1)
    for a in range(n + 1):
        for b in range(a, n + 1):
            if a == 0 and b == 0:
                v = q(0, 0)
            elif a == 0:
                v = (q(0, b) + q(0, -b)) / s2
            else:
                v = (q(a, b) + q(a, -b) + q(-a, b) + q(-a, -b)) / 2
            E[a, b] = v
            E[b, a] = v
    return L, E, pp


def int_cos_cosh(w, L):
    """∫_{−L/2}^{L/2} cos(ωu) cosh(u/2) du."""
    a = mp.mpf(1) / 2
    f = lambda u: (a * mp.sinh(a * u) * mp.cos(w * u) + w * mp.cosh(a * u) * mp.sin(w * u)) / (a * a + w * w)
    return f(L / 2) - f(-L / 2)


def int_sin_sinh(w, L):
    """∫_{−L/2}^{L/2} sin(ωu) sinh(u/2) du."""
    a = mp.mpf(1) / 2
    f = lambda u: (a * mp.cosh(a * u) * mp.sin(w * u) - w * mp.sinh(a * u) * mp.cos(w * u)) / (a * a + w * w)
    return 2 * (f(L / 2) - f(0))


def odd_pole_vector(L, n):
    """w_k = ĝ_k(i/2) = ∫ b_k(u) e^{u/2} du for the odd basis b_k = √(2/L) sin(ω_k u)."""
    return mp.matrix([mp.sqrt(2 / L) * int_sin_sinh(2 * mp.pi * k / L, L) for k in range(1, n + 1)])


def fourier_odd(c, z, L):
    """φ̂(z) = ∫ φ(u) e^{−izu} du for φ = Σ c_k √(2/L) sin(ω_k u), k = 1..len(c); z may be complex."""
    def S(x):
        return L / 2 if abs(x) < mp.mpf(10) ** (-mp.mp.dps + 5) else mp.sin(x * L / 2) / x
    val = mp.mpc(0)
    for k, ck in enumerate(c, start=1):
        w = 2 * mp.pi * k / L
        val += ck * (S(z - w) - S(z + w))
    return -1j * mp.sqrt(2 / L) * val


def fourier_even(c, t, L):
    n = len(c) - 1
    sinc = lambda x: L / 2 if x == 0 else mp.sin(x * L / 2) / x
    v = c[0] * 2 * sinc(t) / mp.sqrt(L)
    for k in range(1, n + 1):
        w = 2 * mp.pi * k / L
        v += c[k] * mp.sqrt(2 / L) * (sinc(t - w) + sinc(t + w))
    return v


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--lambda-sq", type=float, default=13)
    ap.add_argument("--n", type=int, default=60)
    ap.add_argument("--dps", type=int, default=100)
    ap.add_argument("--zeros", type=int, default=50)
    ap.add_argument("--function", choices=["zeta", "dh"], default="zeta")
    ap.add_argument("--variant", choices=["weil", "nopole", "full"], default="weil",
                    help="weil: restrict to φ̂(±i/2) = 0; nopole: drop the pole term, QW = Σ_ρ |φ̂(γ)|², all even φ")
    ap.add_argument("--check-series", action="store_true", help="compare closed-form archimedean terms with nsum")
    ap.add_argument("--json")
    args = ap.parse_args()
    mp.mp.dps = args.dps
    global CHECK_SERIES
    CHECK_SERIES = args.check_series
    from progress import Job

    n_zeros = args.zeros if args.function == "zeta" else 0
    job = Job(f"connes letter {args.function} x={args.lambda_sq}", total=2 + n_zeros, args=vars(args))
    t0 = time.time()
    lambda_sq = mp.mpf(args.lambda_sq)
    L, E, pp = build_form(lambda_sq, args.n, args.function, progress_job=job)
    job.step(f"form built (N={args.n})")
    n1 = args.n + 1
    print(f"λ² = {args.lambda_sq}, N = {args.n}, dps = {args.dps}; prime powers {[int(mp.nint(mp.exp(l))) for l, _ in pp]}; built in {time.time() - t0:.1f}s")

    # Constraint / pole vector v_k = b̂_k(i/2).
    v = mp.matrix(n1, 1)
    v[0] = int_cos_cosh(0, L) / mp.sqrt(L)
    for k in range(1, n1):
        v[k] = mp.sqrt(2 / L) * int_cos_cosh(2 * mp.pi * k / L, L)
    t1 = time.time()
    if args.variant == "weil":
        # Weil's subspace: Householder basis of v^⊥.
        nv = mp.norm(v)
        u = v.copy()
        u[0] += nv if v[0] >= 0 else -nv
        u = u / mp.norm(u)
        H = mp.eye(n1) - 2 * u * u.T
        B = mp.matrix([[H[r, c] for c in range(1, n1)] for r in range(n1)])
        Q = B.T * E * B
    elif args.variant == "nopole":
        # Drop the pole term: QW + 2 φ̂(i/2)² = Σ_ρ |φ̂(γ)|² on all even φ.
        B = mp.eye(n1)
        Q = E + 2 * v * v.T
    else:
        # The form as built (DH has no pole term).
        B = mp.eye(n1)
        Q = E
    vals, vecs = mp.eigsy(Q)
    order = sorted(range(len(vals)), key=lambda i: vals[i])
    print(f"variant {args.variant}; eigensolve {time.time() - t1:.1f}s; smallest eigenvalues: " + ", ".join(mp.nstr(vals[i], 6) for i in order[:4]))
    job.step("eigensolve done")
    job.result("smallest eigenvalues " + ", ".join(mp.nstr(vals[i], 6) for i in order[:4]))
    y = mp.matrix([vecs[r, order[0]] for r in range(Q.rows)])
    c = B * y
    print(f"η̂(i/2) = {mp.nstr(sum(c[k] * v[k] for k in range(n1)), 6)}")

    # Zeros of η̂ near each zeta zero (skipped for DH).
    rows = []
    for k in range(1, (args.zeros if args.function == "zeta" else 0) + 1):
        g = mp.zetazero(k).imag
        try:
            z = mp.findroot(lambda t: fourier_even(c, t, L), g)
            diff = abs(z - g)
        except Exception:
            z, diff = None, None
        rows.append((k, g, z, diff))
        job.step(f"zero {k}: |diff| = {mp.nstr(diff, 3) if diff is not None else 'no root'}")
    print("k   γ_k                      |zero(η̂) − γ_k|")
    for k, g, z, d in rows:
        print(f"{k:<3} {mp.nstr(g, 20):<24} {mp.nstr(d, 6) if d is not None else 'no root'}")
    if args.json:
        with open(args.json, "w") as f:
            json.dump(
                {
                    "lambda_sq": args.lambda_sq,
                    "n": args.n,
                    "dps": args.dps,
                    "variant": args.variant,
                    "eigenvalues": [mp.nstr(vals[i], 30) for i in order[:6]],
                    "zeros": [{"k": k, "gamma": mp.nstr(g, 40), "diff": mp.nstr(d, 10) if d is not None else None} for k, g, z, d in rows],
                },
                f,
                indent=2,
            )
    job.done()
    print(f"total {time.time() - t0:.1f}s")


if __name__ == "__main__":
    main()
