"""Controls: decompose Q = Σ_pairs T_i + P with the census off-line zeros; check P ⪰ 0; the forward
Sherman–Morrison quantity σ₁ for the first pair (P₁ = Q − T₁ holds the other pairs); the naive criterion.
Usage: python controls.py --function dh --sector even --x 10,13,... --n 80 --dps 70 --json out.json
"""
import argparse, json, os, sys, time
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(os.getcwd(), "scripts"))
import mpmath as mp
import r1lib as R
from control_normalisation import rh2_matrix
from progress import Job

# census off-line zeros (β > ½ representatives); docs/DAVENPORT_HEILBRONN.md, docs/EPSTEIN.md, docs/CONDUCTOR5_FAMILY.md
PAIRS = {
    "dh": [("0.8085171824566382", "85.69934848537758"), ("0.650830", "114.163343"), ("0.574356", "166.479306"),
           ("0.724258", "176.702461"), ("0.869531", "240.404672")],
    "z1": [("0.9329697", "15.66824953"), ("0.9376669", "29.98339524"), ("0.6969271", "36.37406369"),
           ("0.8231873", "44.00011318"), ("0.6349076", "46.75840866"), ("0.5974373", "52.74326641"),
           ("0.8111721", "58.50773733")],
    "ftstar": [("0.75", "0")],
}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--function", required=True)
    ap.add_argument("--sector", required=True)
    ap.add_argument("--x", required=True)
    ap.add_argument("--n", type=int, default=64)
    ap.add_argument("--dps", type=int, default=50)
    ap.add_argument("--json", required=True)
    args = ap.parse_args()
    mp.mp.dps = args.dps
    xs = args.x.split(",")
    job = Job(f"r1 controls {args.function} {args.sector}", total=len(xs), args=vars(args))
    rows = []
    for xs_ in xs:
        t0 = time.time()
        x = mp.exp(mp.mpf(xs_[2:])) if xs_.startswith("e^") else mp.mpf(xs_)
        L = mp.log(x)
        Q = rh2_matrix(args.function, x, args.n, args.sector)
        n = args.n
        Ts, cs, ss, ms = [], [], [], []
        for b_s, g_s in PAIRS[args.function]:
            d = mp.mpf(b_s) - mp.mpf(1) / 2
            Tm, c, s, m = R.pair_matrix(d, mp.mpf(g_s), L, n, args.sector)
            Ts.append(Tm); cs.append(c); ss.append(s); ms.append(m)
        P_all = Q
        for Tm in Ts:
            P_all = P_all - Tm
        P1 = Q - Ts[0]
        lamQ = R.lam_min(Q)
        lamP = R.lam_min(P_all)
        lamP1 = R.lam_min(P1)
        row = {"function": args.function, "sector": args.sector, "x": mp.nstr(x, 12), "n": n, "dps": args.dps,
               "lam_min_Q": mp.nstr(lamQ, 12), "lam_min_P_all": mp.nstr(lamP, 12), "lam_min_P1": mp.nstr(lamP1, 12),
               "naive_min_pair1": mp.nstr(R.naive_min(cs[0], ss[0], ms[0]), 8)}
        if lamP1 > 0:
            spec = R.Spectral(P1)
            if ms[0] == 2 and PAIRS[args.function][0][1] == "0":
                # real pair: even +2ccᵀ (never fails), odd −2ssᵀ: fails iff 2 sᵀP⁻¹s > 1
                if args.sector == "odd":
                    sigma = 2 * spec.inv_form(ss[0])
                else:
                    sigma = mp.mpf(0)
                    row["reverse_R2_even_realpair"] = mp.nstr(2 * spec.inv_form(cs[0]), 10)
            else:
                sigma = R.forward_test(spec, cs[0], ss[0], ms[0])
            row["sigma1"] = mp.nstr(sigma, 12)
            row["ln_sigma1"] = mp.nstr(mp.log(sigma), 10) if sigma > 0 else None
        # attribution on the minimiser of Q
        vals, vecs = mp.eigsy(Q)
        i0 = min(range(len(vals)), key=lambda i: vals[i])
        f = mp.matrix([vecs[r, i0] for r in range(n if args.sector == "odd" else n + 1)])
        row["on_minimiser"] = {"Q": mp.nstr((f.T * Q * f)[0], 8), "T_pair1": mp.nstr((f.T * Ts[0] * f)[0], 8),
                               "T_all_pairs": mp.nstr(mp.fsum((f.T * Tm * f)[0] for Tm in Ts), 8), "P_all": mp.nstr((f.T * P_all * f)[0], 8)}
        rows.append(row)
        job.result(f"{args.function} {args.sector} x={xs_}: λ_min(Q)={row['lam_min_Q']} λ_min(P_all)={row['lam_min_P_all']} "
                   f"λ_min(P1)={row['lam_min_P1']} σ1={row.get('sigma1')} naive={row['naive_min_pair1']} ({time.time()-t0:.0f}s)")
        job.step(xs_)
        with open(args.json, "w") as fh:
            json.dump(rows, fh, indent=1)
    job.done()


if __name__ == "__main__":
    main()
