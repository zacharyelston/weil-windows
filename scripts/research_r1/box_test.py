"""The local-count regime of the reverse test, at any height, without a matrix.
Test function f = cos(γ₀ u) on [−a, a] (box window; the box maximises (∫w)²/∫w², so it is the best
single-window lower bound on κ(γ₀) = sup |⟨c,f⟩|²/Q(f) up to the shaping gain the matrix finds).
Q(f) = pole(f) + (1/π)∫₀^∞ Ψ(t)|F̂(t)|² dt,  Ψ = Re ψ(¼ + it/2) − log π − Σ_{log n<2a} (2Λ(n)/√n) cos(t log n),
F̂(t) = sin((t−γ₀)a)/(t−γ₀) + sin((t+γ₀)a)/(t+γ₀).
Numerics: unit panels on [max(0, γ₀−W), γ₀+W]; outside, sin² is replaced by its mean ½ and the comb by 0
(error O(1/(aW²)) relative), with the log tail integrated exactly.
R4_lower(δ) = 4⟨c,f⟩²/(Q(f) + 4⟨s,f⟩²) ≤ R4(γ₀, δ);  ⟨c,f⟩ = ∫ f cosh δu cos γ₀u, ⟨s,f⟩ = ∫ f sinh δu sin γ₀u.
Usage: python box_test.py --a 0.8 --bands -1,-0.5,0,0.5,1,2,3,4 --per-band 8 --delta 0,0.5 --json out.json
   (bands are offsets of log(γ₀/2π) from 4a; --gammas overrides with explicit heights)
"""
import argparse, json, os, sys, time
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(os.getcwd(), "scripts"))
import mpmath as mp
import connes_letter_mp as cl
from progress import Job


def Q_box(a, g0, pp, W=150):
    """Q(f) for f = cos(γ₀u)·1_[−a,a]. Returns (Q, arch, primes, pole, norm2)."""
    def F(t):
        tm, tp = t - g0, t + g0
        v = (a if abs(tm) < mp.mpf("1e-12") else mp.sin(tm * a) / tm)
        v += (a if abs(tp) < mp.mpf("1e-12") else mp.sin(tp * a) / tp)
        return v
    def Psi_arch(t):
        return mp.re(mp.digamma(mp.mpc(mp.mpf(1) / 4, t / 2))) - mp.log(mp.pi)
    def Psi_pr(t):   # Σ (2Λ(n)/√n) cos(t log n)
        return mp.fsum(2 * coef * mp.cos(t * ln) for ln, coef in pp)
    lo = max(mp.mpf(0), g0 - W)
    hi = g0 + W
    pts = [lo + k for k in range(int(mp.floor(hi - lo)) + 1)] + [hi]
    arch = mp.quad(lambda t: Psi_arch(t) * F(t) ** 2, pts) / mp.pi
    prim = mp.quad(lambda t: Psi_pr(t) * F(t) ** 2, pts) / mp.pi
    # averaged tails: |F̂|² → ½/τ² around ±γ₀ (the near peak), comb → 0, Ψ_arch → log(t/2π)
    tail = mp.mpf(0)
    if lo > 0:   # [0, lo]: τ = γ₀ − t ∈ [W, γ₀]; include both peaks' contributions
        tail += mp.quad(lambda t: mp.log(max(t, mp.mpf(1)) / (2 * mp.pi)) * (mp.mpf(1) / 2) * (1 / (g0 - t) ** 2 + 1 / (g0 + t) ** 2), [0, lo]) / mp.pi
    tail += mp.quad(lambda t: mp.log(t / (2 * mp.pi)) * (mp.mpf(1) / 2) * (1 / (t - g0) ** 2 + 1 / (t + g0) ** 2), [hi, 10 * hi, mp.inf]) / mp.pi
    arch += tail
    # pole term 2(∫ f cosh(u/2))²; ∫_{−a}^{a} cos(γ₀u) cosh(u/2) du in closed form
    h = mp.mpf(1) / 2
    pole_int = 2 * (h * mp.sinh(h * a) * mp.cos(g0 * a) + g0 * mp.cosh(h * a) * mp.sin(g0 * a)) / (h * h + g0 * g0)
    pole = 2 * pole_int ** 2
    norm2 = a + mp.sin(2 * g0 * a) / (2 * g0)
    return arch - prim + pole, arch, prim, pole, norm2, tail


def cs_box(a, g0, d):
    """⟨c,f⟩ = ∫ cos²(γ₀u) cosh(δu) du and ⟨s,f⟩ = ∫ cos(γ₀u) sin(γ₀u) sinh(δu) du over [−a, a], closed forms."""
    def Icc(nu):   # ∫ cos(νu) cosh(δu)
        if d == 0:
            return 2 * a if nu == 0 else 2 * mp.sin(nu * a) / nu
        return 2 * (d * mp.sinh(d * a) * mp.cos(nu * a) + nu * mp.cosh(d * a) * mp.sin(nu * a)) / (d * d + nu * nu)
    def Iss(nu):   # ∫ sin(νu) sinh(δu)
        if d == 0 or nu == 0:
            return mp.mpf(0)
        return 2 * (d * mp.cosh(d * a) * mp.sin(nu * a) - nu * mp.sinh(d * a) * mp.cos(nu * a)) / (d * d + nu * nu)
    cf = (Icc(0) + Icc(2 * g0)) / 2
    sf = Iss(2 * g0) / 2
    return cf, sf


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--a", required=True)
    ap.add_argument("--delta", default="0,0.5")
    ap.add_argument("--gammas", default=None)
    ap.add_argument("--bands", default="-1,-0.5,0,0.5,1,2,3,4")
    ap.add_argument("--per-band", type=int, default=8)
    ap.add_argument("--dps", type=int, default=15)
    ap.add_argument("--W", type=float, default=150)
    ap.add_argument("--json", required=True)
    args = ap.parse_args()
    mp.mp.dps = args.dps
    a = mp.mpf(args.a)
    x = mp.exp(2 * a)
    pp = [(mp.log(nn), lp_ / mp.sqrt(nn)) for nn, lp_ in cl.prime_powers(int(mp.floor(x)) + 1)]
    pp = [(ln, c) for ln, c in pp if ln < 2 * a]
    A_full = mp.fsum(2 * c for _, c in pp)
    A_box = mp.fsum(2 * c * (1 - ln / (2 * a)) for ln, c in pp)   # comb weighted by the box autocorrelation
    if args.gammas:
        gammas = [(None, mp.mpf(v)) for v in args.gammas.split(",")]
    else:
        gammas = []
        for b in args.bands.split(","):
            centre = 4 * a + mp.mpf(b)
            for k in range(args.per_band):
                # spread the samples over one e-fold around the band centre with an irrational stride
                off = (k * mp.sqrt(2)) % 1 - mp.mpf("0.5")
                gammas.append((b, 2 * mp.pi * mp.exp(centre + off)))
    deltas = [mp.mpf(v) for v in args.delta.split(",")]
    job = Job(f"r1 box test a={args.a}", total=len(gammas), args=vars(args))
    rows = []
    for band, g0 in gammas:
        t0 = time.time()
        Q, arch, prim, pole, norm2, tail = Q_box(a, g0, pp, W=mp.mpf(args.W))
        row = {"a": args.a, "band": band, "gamma0": mp.nstr(g0, 12), "Q": mp.nstr(Q, 8), "arch": mp.nstr(arch, 8), "primes": mp.nstr(prim, 8),
               "pole": mp.nstr(pole, 4), "norm2": mp.nstr(norm2, 8), "tail": mp.nstr(tail, 4), "log_g0_over_2pi": mp.nstr(mp.log(g0 / (2 * mp.pi)), 8),
               "A_full": mp.nstr(A_full, 6), "A_box": mp.nstr(A_box, 6), "rho_eff": mp.nstr(mp.log(g0 / (2 * mp.pi)) - Q / norm2, 6), "R4_lower": {}}
        for d in deltas:
            cf, sf = cs_box(a, g0, d)
            row["R4_lower"][mp.nstr(d, 4)] = mp.nstr(4 * cf * cf / (Q + 4 * sf * sf), 8)
        rows.append(row)
        job.result(f"a={args.a} band={band} γ₀={mp.nstr(g0, 8)}: Q/‖f‖²={mp.nstr(Q / norm2, 6)} (log(γ₀/2π)={row['log_g0_over_2pi']}, primes/‖f‖²={mp.nstr(prim / norm2, 5)}); R4_lower={row['R4_lower']} ({time.time()-t0:.0f}s)")
        job.step()
        with open(args.json, "w") as fh:
            json.dump({"rows": rows, "A_full": mp.nstr(A_full, 8), "A_box": mp.nstr(A_box, 8), "args": vars(args)}, fh, indent=1)
    job.done()


if __name__ == "__main__":
    main()
