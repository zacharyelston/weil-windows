#!/usr/bin/env python3
"""Root-number rule: measured ln(lambda_odd/lambda_even) at v = 10 against the Hankel-prolate leakage ratio
ln[(1 - lambda_{nu,1}(c)) / (1 - lambda_{nu,0}(c))], nu = k - 1, c = 40 pi (the E-map ground state m = 0 sits in the
sector (-1)^s = epsilon, so the predicted sign of ln(lambda_odd/lambda_even) is epsilon).

Usage: .venv/bin/python scripts/research_r2/root_number_gap.py --json data/research_r2/root_number_gap.json
"""
import argparse
import glob
import json
import os
import sys

import mpmath as mp

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from hankel_prolate import leakage  # noqa: E402
from exact_leakage import EPS  # noqa: E402


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--json")
    args = ap.parse_args()
    c = 4 * mp.pi * 10
    mp.mp.dps = int(2 * float(c) / 2.3) + 40
    files = glob.glob("data/connes/gl2/scan9_*.json") + glob.glob("data/connes/gl2/scan11_*.json") + glob.glob("data/connes/gl2/scan12_*.json")
    out = []
    for f in sorted(files):
        rows = json.load(open(f))
        obj = rows[0]["obj"]
        d = {r["parity"]: mp.mpf(r["lambda_n2"]) for r in rows if abs(float(r["v"]) - 10) < 1e-9}
        if len(d) != 2:
            continue
        k = int(rows[0]["weight"])
        nu = mp.mpf(k - 1)
        meas = mp.log(d["odd"] / d["even"])
        pro = mp.log(leakage(nu, c, 1) / leakage(nu, c, 0))
        eps = EPS[obj]
        out.append({"object": obj, "level": rows[0]["N"], "weight": k, "eps": eps, "measured_ln_odd_over_even": mp.nstr(meas, 5),
                    "prolate_predicted": mp.nstr(eps * pro, 5)})
        print(f"{obj} N={rows[0]['N']} k={k} eps={eps:+d}: measured {mp.nstr(meas, 4)}  predicted eps*ln(l1/l0) = {mp.nstr(eps * pro, 4)}")
    if args.json:
        with open(args.json, "w") as fh:
            json.dump(out, fh, indent=1)


if __name__ == "__main__":
    main()
