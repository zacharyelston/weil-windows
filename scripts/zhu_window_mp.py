#!/usr/bin/env python3
"""ζ's finite-basis minimum of the zeros-side Weil form at Zhu's window, x = e^{2a}.

Zhu (arXiv:2608.24827 v2) claims Q(f) ≥ 8.9·10⁻¹⁸ ‖f‖² for f supported in
[−0.8, 0.8], i.e. a = 0.8 and our x = e^{1.6} ≈ 4.953. Finite trigonometric
minima decrease to the full-space infimum (Connes–Consani, Cor. 2.4), so each
value printed here is an upper bound on that infimum in its parity sector:

  even sector: QW = E + 2vvᵀ   (v_k = b̂_k(i/2), cosine basis, k = 0..N)
  odd sector:  QW = E − 2wwᵀ   (w_k = ĝ_k(i/2), sine basis, k = 1..N)

Whether these match Zhu's Q depends on the normalisation map
(docs/BRIEF_ZHU_AUDIT.md, task 1); the numbers are reported as our own.

Usage: .venv/bin/python scripts/zhu_window_mp.py --a 0.8 --n 64,80,100 --dps 130 --json data/connes/zhu_window.json
"""

import argparse
import json
import os
import sys

import mpmath as mp

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import connes_letter_mp as cl  # noqa: E402
from progress import Job  # noqa: E402


def sector_minimum(x, n, parity):
    L, E, _ = cl.build_form(x, n, "zeta", parity=parity)
    if parity == "even":
        v = mp.matrix([cl.int_cos_cosh(2 * mp.pi * k / L, L) * (1 / mp.sqrt(L) if k == 0 else mp.sqrt(2 / L)) for k in range(n + 1)])
        Q = E + 2 * v * v.T
    else:
        w = cl.odd_pole_vector(L, n)
        Q = E - 2 * w * w.T
    vals = sorted(mp.eigsy(Q, eigvals_only=True))
    return vals[:3]


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--a", default="0.8", help="half-width of Zhu's window; x = e^{2a}")
    ap.add_argument("--n", default="64,80,100")
    ap.add_argument("--dps", type=int, default=130)
    ap.add_argument("--json")
    args = ap.parse_args()
    mp.mp.dps = args.dps
    a = mp.mpf(args.a)
    x = mp.exp(2 * a)
    n_list = [int(s) for s in args.n.split(",")]
    job = Job(f"zhu window a={args.a}", total=2 * len(n_list), args=vars(args))
    job.log(f"x = e^(2a) = {mp.nstr(x, 20)}; prime powers ≤ x: {[m for m, _ in cl.prime_powers(int(mp.floor(x)))]}")
    rows = []
    for n in n_list:
        for parity in ["even", "odd"]:
            vals = sector_minimum(x, n, parity)
            job.result(f"N={n} {parity}: eps0={mp.nstr(vals[0], 12)} eps1={mp.nstr(vals[1], 6)}")
            job.step(f"N={n} {parity}")
            rows.append({"a": args.a, "x": mp.nstr(x, 25), "n": n, "dps": args.dps, "parity": parity,
                         "eps": [mp.nstr(t, 20) for t in vals]})
            if args.json:
                os.makedirs(os.path.dirname(os.path.abspath(args.json)), exist_ok=True)
                with open(args.json, "w") as fh:
                    json.dump(rows, fh, indent=2)
    job.done()


if __name__ == "__main__":
    main()
