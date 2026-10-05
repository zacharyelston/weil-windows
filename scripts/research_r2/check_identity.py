#!/usr/bin/env python3
"""Theorem 1 check: F_win(gamma) + F_out(gamma) = F_psi(gamma) = 0 at nontrivial zeros, for the Hermite E-map window
functions (zeta both sectors, L(chi_-4) even, L(chi_5) odd) at x = 4, 40 digits.

Usage: .venv/bin/python scripts/research_r2/check_identity.py --json data/research_r2/check_identity.json
"""
import argparse
import json
import os
import sys

import mpmath as mp

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
import r2lib  # noqa: E402
from decay_law_mp import hardy_zeros  # noqa: E402


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--json")
    args = ap.parse_args()
    mp.mp.dps = 40
    x = mp.mpf(4)
    L = mp.log(x)
    out = []
    for D, kappa, s, pole in [(1, 0, 0, True), (1, 0, 1, True), (-4, 1, 0, False), (5, 0, 1, False)]:
        psi = r2lib.hermite_choice(kappa, s, pole)
        win = r2lib.Window(psi, D, x)
        umax = mp.log(win.ycut * win.sq)
        trig = mp.cos if s == 0 else mp.sin
        Fw = lambda t: mp.quad(lambda u: win.f(u) * trig(t * u), [-L / 2, 0, L / 2])
        Fo = lambda t: 2 * mp.quad(lambda u: win.f(u) * trig(t * u), [L / 2, L / 2 + 1, umax])
        zs = [mp.zetazero(j).imag for j in (1, 2, 3)] if D == 1 else hardy_zeros(D, 25)[:3]
        for g in zs:
            a, b = Fw(g), Fo(g)
            row = {"D": D, "sector": "even" if s == 0 else "odd", "gamma": mp.nstr(g, 15), "F_win": mp.nstr(a, 10),
                   "F_out": mp.nstr(b, 10), "sum": mp.nstr(a + b, 3)}
            out.append(row)
            print(row, flush=True)
    worst = max(abs(mp.mpf(r["sum"])) for r in out)
    print(f"max |F_win + F_out| at zeros = {mp.nstr(worst, 3)}")
    if args.json:
        with open(args.json, "w") as fh:
            json.dump(out, fh, indent=1)


if __name__ == "__main__":
    main()
