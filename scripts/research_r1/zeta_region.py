"""ζ: the reverse (consistency) region at a certified support.
For each (a, sector, N): eigendecompose rh2's full zeros-side matrix Q_a; then
  (1) comb check: R2(γ_k) = 2 cᵀQ⁻¹c at the true zeros γ_k (δ = 0, simple pair) must be ≤ 1; between zeros it is huge;
  (2) for δ in a list: the profile R4(δ, γ₀) = 4 cᵀ(Q + 4ssᵀ)⁻¹c on a γ₀ grid and the last crossing T_excl(δ) of R4 = 1;
  (3) even sector: the real-pair test R2(δ, 0) = 2 cᵀQ⁻¹c (c_k = ∫ b_k cosh δu).
Usage: python zeta_region.py --jobs 1.495:even:120:150,1.495:odd:120:150 --json out.json
"""
import argparse, json, os, sys, time
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(os.getcwd(), "scripts"))
import mpmath as mp
import r1lib as R
from control_normalisation import rh2_matrix
from progress import Job


DELTAS = ["0", "0.01", "0.05", "0.1", "0.25", "0.4", "0.5"]
NZEROS = 50


def run_job(a_s, sector, n, dps, job, zeros_cache):
    mp.mp.dps = dps
    a = mp.mpf(a_s)
    x = mp.exp(2 * a)
    L = 2 * a
    Tstar = 2 * mp.pi * x
    t0 = time.time()
    Q = rh2_matrix("zeta", x, n, sector)
    spec = R.Spectral(Q)
    job.log(f"a={a_s} {sector} N={n} dps={dps}: built+eig in {time.time()-t0:.1f}s; λ_min = {mp.nstr(spec.lam_min, 8)}; T* = {mp.nstr(Tstar, 6)}")
    out = {"a": a_s, "x": mp.nstr(x, 15), "sector": sector, "n": n, "dps": dps, "lam_min": mp.nstr(spec.lam_min, 15), "Tstar": mp.nstr(Tstar, 10)}
    # (1) comb check at the true zeros
    if dps not in zeros_cache:
        zeros_cache[dps] = [mp.zetazero(k).imag for k in range(1, NZEROS + 1)]
    gam = zeros_cache[dps]
    comb = []
    for k, g in enumerate(gam, 1):
        c, s = R.pair_vectors(0, g, L, n, sector)
        r2 = 2 * spec.inv_form(c)
        comb.append((k, mp.nstr(g, 12), mp.nstr(r2, 10)))
    g1 = gam[0]
    shifted = []
    for eps in ("1e-5", "1e-10", "1e-20", "1e-30"):
        c, s = R.pair_vectors(0, g1 + mp.mpf(eps), L, n, sector)
        shifted.append((eps, mp.nstr(2 * spec.inv_form(c), 6)))
    # midpoints between zeros
    mids = []
    for k in range(len(gam) - 1):
        gm = (gam[k] + gam[k + 1]) / 2
        c, s = R.pair_vectors(0, gm, L, n, sector)
        mids.append((mp.nstr(gm, 8), mp.nstr(2 * spec.inv_form(c), 6)))
    out["comb_R2_at_zeros"] = comb
    out["comb_R2_max_at_zeros"] = mp.nstr(max(mp.mpf(r) for _, _, r in comb), 12)
    out["comb_R2_shifted_from_gamma1"] = shifted
    out["comb_R2_midpoints"] = mids
    job.log(f"  comb: max 2κ at zeros = {out['comb_R2_max_at_zeros']}; at γ₁+1e-20: {shifted[2][1]}; midpoints min/max {min(mp.mpf(v) for _, v in mids)} / {max(mp.mpf(v) for _, v in mids)}")
    # (2) profiles and T_excl
    deltas = DELTAS
    # box-test-function estimate: T_excl ≈ 2π e^{4a} = 2π x²; scan to 1.4× that (and at least 6 T*)
    # ... but the cosine basis only represents heights up to ω_N = 2πN/L; beyond that the matrix route is blind.
    omega_N = 2 * mp.pi * n / L
    gmax = omega_N          # scan the whole range the basis can represent; a crossing at the end means "> ω_N"
    out["omega_N"] = mp.nstr(omega_N, 8)
    h = Tstar / 30
    grid = [mp.mpf("0.5") + i * h for i in range(int(gmax / h) + 1)]
    profiles, texcl = {}, {}
    for d_s in deltas:
        d = mp.mpf(d_s)
        m = 4 if d != 0 else 2
        prof = []
        last_above = None
        for g in grid:
            c, s = R.pair_vectors(d, g, L, n, sector)
            r, cQc, sQs = R.reverse_test(spec, c, s, m if d != 0 else 4)   # δ = 0: test a DOUBLE on-line zero (m = 4)
            prof.append((float(g), float(mp.log10(r)) if r > 0 else None))
            if r > 1:
                last_above = g
        # bisect the last crossing
        if last_above is not None and last_above < grid[-1]:
            lo, hi = last_above, last_above + h
            for _ in range(30):
                mid = (lo + hi) / 2
                c, s = R.pair_vectors(d, mid, L, n, sector)
                r, _, _ = R.reverse_test(spec, c, s, 4)
                if r > 1:
                    lo = mid
                else:
                    hi = mid
            tex = (lo + hi) / 2
        else:
            tex = last_above
        profiles[d_s] = prof
        at_end = last_above is not None and last_above >= grid[-1]
        texcl[d_s] = (">" if at_end else "") + (mp.nstr(tex, 10) if tex is not None else "none")
        job.log(f"  δ={d_s}: T_excl = {texcl[d_s]}  (T_excl/T* = {mp.nstr(tex / Tstar, 6) if tex is not None else None}; basis limit ω_N = {mp.nstr(omega_N, 6)})")
    out["profiles_log10R4"] = profiles
    out["T_excl"] = texcl
    # (3) real pair, even sector
    if sector == "even":
        rp = []
        for d_s in ("0.05", "0.1", "0.25", "0.4", "0.49"):
            c, s = R.pair_vectors(mp.mpf(d_s), 0, L, n, "even")
            rp.append((d_s, mp.nstr(2 * spec.inv_form(c), 10)))
        out["real_pair_even_R2"] = rp
        job.log(f"  even real-pair R2 = {rp}")
    job.step(f"a={a_s} {sector} N={n}")
    return out


def main():
    global DELTAS, NZEROS
    ap = argparse.ArgumentParser()
    ap.add_argument("--jobs", required=True, help="a:sector:n:dps,...")
    ap.add_argument("--json", required=True)
    ap.add_argument("--deltas", default=",".join(DELTAS))
    ap.add_argument("--nzeros", type=int, default=NZEROS)
    args = ap.parse_args()
    DELTAS = args.deltas.split(",")
    NZEROS = args.nzeros
    jobs = [j.split(":") for j in args.jobs.split(",")]
    job = Job(f"r1 zeta region {args.jobs}", total=len(jobs), args=vars(args))
    zeros_cache = {}
    results = []
    for a_s, sector, n, dps in jobs:
        results.append(run_job(a_s, sector, int(n), int(dps), job, zeros_cache))
        with open(args.json, "w") as fh:
            json.dump(results, fh, indent=1)
    job.done()


if __name__ == "__main__":
    main()
