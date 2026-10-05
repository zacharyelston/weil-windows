#!/usr/bin/env python3
"""Checks on the zero-count inputs of the certified bound (docs/R2_THEOREMS.md section 2). Each can fail.

1. N+(t) >= N(t) at every zero ordinate up to a height H (and just after it), with N(t) counted from zeros found
   independently: zeta via mpmath.zetazero (first 200 ordinates), L(s, chi) via sign changes of the real
   completed L-function on the line (decay_law_mp.hardy_zeros, even characters only).
2. The closed-form tail integrals int_T^inf N+(t) t^-3 dt (cert_lib) against mpmath quadrature of the same N+.
3. R2's zeta count (Trudgian 2014 journal constants, used in scripts/research_r2/r2lib.py) against the one used here
   (Bellotti-Wong 2025 Thm 1.1) on the range of T_c that occurs.

Usage: .venv/bin/python scripts/research_r2_rigorous/zero_count_check.py --json data/research_r2_rigorous/zero_count_check.json
"""
import argparse
import json
import os
import sys

import mpmath as mp

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
import cert_lib as C  # noqa: E402
from flint import arb  # noqa: E402


def r2_zeta_N(t):
    """R2's N+ (scripts/research_r2/r2lib.py zeta_Nplus, doubled): theta/pi + 1 + 0.112 log t + 0.278 loglog t + 2.510."""
    t = mp.mpf(t)
    return 2 * (mp.siegeltheta(t) / mp.pi + 1 + mp.mpf("0.112") * mp.log(t) + mp.mpf("0.278") * mp.log(mp.log(t)) + mp.mpf("2.510"))


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--json", required=True)
    args = ap.parse_args()
    C.setprec()
    mp.mp.dps = 30
    out = {}
    # ---- 1. zeta
    zs = [mp.im(mp.zetazero(k)) for k in range(1, 201)]
    worst = None
    fails = 0
    for k, g in enumerate(zs, start=1):
        for t in (g, g + mp.mpf("1e-9")):
            Nt = 2 * k                                   # |Im rho| <= t: k zeros above, k below
            ub = C.zeta_Nplus(arb(mp.nstr(t, 25)))
            slack = float(ub.lower()) - Nt
            if slack < 0:
                fails += 1
            worst = slack if worst is None or slack < worst else worst
    out["zeta"] = {"zeros_used": 200, "max_height": mp.nstr(zs[-1], 10), "failures": fails, "min_slack": worst}
    print(f"zeta: N+(t) >= N(t) at 400 heights up to {mp.nstr(zs[-1], 8)}: failures {fails}, min slack {worst:.3f}")
    # ---- 1b. even characters (hardy_zeros handles root number 1, chi even)
    from decay_law_mp import hardy_zeros
    out["dirichlet"] = {}
    for D in (5, 8, 12, 13, 17):
        q = abs(D)
        gs = hardy_zeros(D, mp.mpf(60), step=mp.mpf("0.05"))
        fails, worst = 0, None
        for k, g in enumerate(gs, start=1):
            for t in (g, g + mp.mpf("1e-9")):
                Nt = 2 * k
                if arb(mp.nstr(t, 25)) <= C.dirichlet_tlow(q):
                    continue
                ub = C.dirichlet_Nplus(arb(mp.nstr(t, 25)), q, 1)
                slack = float(ub.lower()) - Nt
                fails += slack < 0
                worst = slack if worst is None or slack < worst else worst
        out["dirichlet"][str(D)] = {"zeros_found_to_60": len(gs), "failures": int(fails), "min_slack": worst}
        print(f"chi_{D}: {len(gs)} zeros to height 60, failures {fails}, min slack {worst:.3f}")
    # ---- 2. closed-form tails vs quadrature
    rows = []
    for T in (14, 20, 40, 100, 400):
        cf = C.zeta_tail_integral(arb(T))
        num = mp.quad(lambda t: mp.mpf(C.zeta_Nplus(arb(mp.nstr(t, 30))).mid().str(30, radius=False)) / t ** 3, [T, 10 * T, 100 * T, mp.inf])
        rows.append({"L": "zeta", "T": T, "closed_form_upper": cf.upper().str(12, radius=False), "quad": mp.nstr(num, 12),
                     "closed_ge_quad": float(cf.upper()) >= float(num) * (1 - 1e-10)})
    for D, chim1 in ((5, 1), (-4, -1), (-20, -1), (13, 1)):
        q = abs(D)
        for T in (10, 30, 100, 400):
            if arb(T) < C.dirichlet_tlow(q):
                continue
            cf = C.dirichlet_tail_integral(arb(T), q, chim1)
            num = mp.quad(lambda t: mp.mpf(C.dirichlet_Nplus(arb(mp.nstr(t, 30)), q, chim1).mid().str(30, radius=False)) / t ** 3,
                          [T, 10 * T, 100 * T, mp.inf])
            rows.append({"L": f"chi_{D}", "T": T, "closed_form_upper": cf.upper().str(12, radius=False), "quad": mp.nstr(num, 12),
                         "closed_ge_quad": float(cf.upper()) >= float(num) * (1 - 1e-10),
                         "ratio": mp.nstr(mp.mpf(cf.upper().str(20, radius=False)) / num, 8)})
    out["tail_integrals"] = rows
    for r in rows:
        print(f"  {r['L']:>8} T={r['T']:>4}: closed {r['closed_form_upper']}  quad {r['quad']}  ok={r['closed_ge_quad']}")
    # ---- 3. R2's zeta count against Bellotti-Wong on 14..2000
    comp = []
    for T in (14.2, 20, 30, 50, 100, 300, 1000, 2000):
        comp.append({"T": T, "R2_trudgian_N": mp.nstr(r2_zeta_N(T), 8), "BW_N": C.zeta_Nplus(arb(T)).str(8, radius=False)})
    out["r2_vs_bw"] = comp
    for r in comp:
        print(f"  T={r['T']:>7}: R2 (Tru14) N+ = {r['R2_trudgian_N']}   BW N+ = {r['BW_N']}")
    with open(args.json, "w") as fh:
        json.dump(out, fh, indent=1)


if __name__ == "__main__":
    main()
