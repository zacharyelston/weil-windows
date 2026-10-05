#!/usr/bin/env python3
"""P-ST (docs/DECAY_LAW.md): window minima of two torus spectral zeta functions with the same seam.
  Z_sq = 4 ζ L(χ_{-4})  (torus of x² + y², class number 1, Euler product): form = zeros_side(1) + zeros_side(-4)
  Z1   = ½[ζ_K + L(χ_{-4})L(χ_5)]  (torus of x² + 5y², class number 2, no Euler product): control_normalisation "z1"
Usage: .venv/bin/python scripts/square_torus_mp.py --json data/connes/decay/square_torus.json
"""
import argparse, json, os, sys
import mpmath as mp
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import decay_law_mp as dl                      # noqa: E402
from control_normalisation import native_matrix   # noqa: E402
from progress import Job                        # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--a2", default="1.6,2.0,2.38,2.6,2.99,3.2,3.5", help="supports 2a")
    ap.add_argument("--n", default="64,96")
    ap.add_argument("--json")
    args = ap.parse_args()
    rows = []
    job = Job("square torus P-ST", total=len(args.a2.split(",")) * 2, args=vars(args))
    for a2 in args.a2.split(","):
        x = mp.exp(mp.mpf(a2))
        for parity in ("even", "odd"):
            out = {}
            for n in [int(v) for v in args.n.split(",")]:
                def build(n=n):
                    return {"Zsq": dl.zeros_side(1, x, n, parity) + dl.zeros_side(-4, x, n, parity),
                            "Z1": native_matrix("z1", x, n, parity)}
                lam, at = dl.stable_lams(build, int(40 + 1.5 * 8 * float(mp.pi) * float(mp.sqrt(x / 4)) / 2.303), method=dl.lam_min_acb)
                out[n] = {k: (mp.nstr(v, 15), at[k]) for k, v in lam.items()}
            job.step(f"2a={a2} x={mp.nstr(x, 6)} {parity}: " + "  ".join(f"N={n}: Zsq={out[n]['Zsq'][0]} Z1={out[n]['Z1'][0]}" for n in out))
            rows.append({"a2": a2, "x": mp.nstr(x, 12), "parity": parity,
                         **{f"Zsq_n{n}": out[n]["Zsq"][0] for n in out}, **{f"Z1_n{n}": out[n]["Z1"][0] for n in out},
                         "stable_dps": {f"{k}_n{n}": out[n][k][1] for n in out for k in out[n]}})
            if args.json:
                with open(args.json, "w") as fh:
                    json.dump(rows, fh, indent=1)
    job.done()


if __name__ == "__main__":
    main()
