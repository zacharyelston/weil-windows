"""Pre-registration support: plant a quadruple ½ ± δ ± iγ₀ into the validated 11a1 form and compute
σ₁(x) = 4 sᵀ(P + 4ccᵀ)⁻¹ s with P = Q_11a1(x) (exact), λ_min of the planted form, and the bisected
crossing. Usage: python plant_gl2.py --cases 0.25:30,0.10:30,0.25:60 --xn 1,1.5,2,2.5,3,4,5,6,8 --sector even --n 64 --dps 60 --json out.json
"""
import argparse, json, os, sys, time
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(os.getcwd(), "scripts"))
import mpmath as mp
import r1lib as R
import gl2_form_mp as g2
from progress import Job

NCOND = 11


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cases", default="0.25:30,0.10:30,0.25:60")
    ap.add_argument("--xn", default="1,1.5,2,2.5,3,4,5,6,8")
    ap.add_argument("--sector", default="even")
    ap.add_argument("--label", default="11a1")
    ap.add_argument("--n", type=int, default=64)
    ap.add_argument("--dps", type=int, default=60)
    ap.add_argument("--json", required=True)
    args = ap.parse_args()
    mp.mp.dps = args.dps
    cases = [tuple(mp.mpf(v) for v in c.split(":")) for c in args.cases.split(",")]
    xns = [mp.mpf(v) for v in args.xn.split(",")]
    job = Job(f"r1 plant {args.label} {args.sector}", total=len(xns), args=vars(args))
    rows = []
    cache = {}

    def P_of(x):
        key = mp.nstr(x, 20)
        if key not in cache:
            cache[key] = g2.zeros_side_E(args.label, x, args.n, args.sector)
        return cache[key]

    for xn in xns:
        x = xn * NCOND
        L = mp.log(x)
        t0 = time.time()
        P = P_of(x)
        lamP = R.lam_min(P)
        spec = R.Spectral(P) if lamP > 0 else None
        row = {"x_over_N": mp.nstr(xn, 6), "x": mp.nstr(x, 12), "lam_min_P": mp.nstr(lamP, 10), "cases": {}}
        for d, g in cases:
            Tm, c, s, m = R.pair_matrix(d, g, L, args.n, args.sector)
            lamQ = R.lam_min(P + Tm)
            sig = R.forward_test(spec, c, s, m) if spec is not None else None
            row["cases"][f"{mp.nstr(d,4)}:{mp.nstr(g,6)}"] = {"lam_min_planted": mp.nstr(lamQ, 10), "sigma1": mp.nstr(sig, 10) if sig is not None else None,
                                                             "naive_min": mp.nstr(R.naive_min(c, s, m), 8)}
        rows.append(row)
        job.result(f"x/N={mp.nstr(xn,4)} x={mp.nstr(x,6)}: λ_min(P)={row['lam_min_P']}; " + "; ".join(f"{k}: λ={v['lam_min_planted']} σ1={v['sigma1']}" for k, v in row["cases"].items()) + f" ({time.time()-t0:.0f}s)")
        job.step()
        with open(args.json, "w") as fh:
            json.dump(rows, fh, indent=1)
    # bisect the crossing of each case between the last positive and first negative grid points
    cross = {}
    for d, g in cases:
        key = f"{mp.nstr(d,4)}:{mp.nstr(g,6)}"
        lams = [(xn, mp.mpf(r["cases"][key]["lam_min_planted"])) for xn, r in zip(xns, rows)]
        lo = hi = None
        for (x1, l1), (x2, l2) in zip(lams, lams[1:]):
            if l1 > 0 and l2 < 0:
                lo, hi = x1, x2
                break
        if lo is None:
            cross[key] = None
            continue
        for _ in range(16):
            mid = (lo + hi) / 2
            x = mid * NCOND
            P = P_of(x)
            Tm, c, s, m = R.pair_matrix(d, g, mp.log(x), args.n, args.sector)
            if R.lam_min(P + Tm) > 0:
                lo = mid
            else:
                hi = mid
        cross[key] = {"x_over_N_lo": mp.nstr(lo, 8), "x_over_N_hi": mp.nstr(hi, 8), "x_lo": mp.nstr(lo * NCOND, 8), "x_hi": mp.nstr(hi * NCOND, 8)}
        job.result(f"crossing {key}: x/N in ({mp.nstr(lo,7)}, {mp.nstr(hi,7)}], x in ({mp.nstr(lo*NCOND,7)}, {mp.nstr(hi*NCOND,7)}]")
    with open(args.json, "w") as fh:
        json.dump({"rows": rows, "crossings": cross, "args": vars(args)}, fh, indent=1)
    job.done()


if __name__ == "__main__":
    main()
