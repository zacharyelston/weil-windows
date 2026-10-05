#!/usr/bin/env python3
"""Proposition H (docs/R2_THEOREMS.md section 5.1): the Hermite bound has rate e^{-c} for ALL x >= x0 = 2q, with an
explicit certified constant:  lambda_s(x) <= K_H c^{d+1} e^{-c}  (d = deg P, c = 2 pi x/q).

Proof ingredients (all inequalities pointing up; lam = sqrt(x/q) >= lam0 = sqrt 2, Z0^2 = 2 pi lam0^2 = 4 pi):
  - ||f_win(x)||^2 is nondecreasing in x (S(w) does not depend on x; the window grows), so >= N0 := ||f_win(x0)||^2
    (certified lower endpoint from hermite_cert.json, x/q = 2);
  - G_k(Z) <= g_k Z^{k-1} e^{-Z^2/2}, g_k = 1/(1 - (k-1)/Z0^2) (k >= 2), 1 (k <= 1)   [Gamma(s,X) <= X^{s-1}e^{-X}/(1-(s-1)/X)];
  - I_R(Y) <= K_R Y^{deg R - 1} e^{-pi Y^2}, K_R = (2 pi)^{-1/2} sum_k |r_k| g_k (2 pi)^{(k-1)/2} lam0^{k - deg R};
  - sum_n n^j e^{-pi lam^2 (n^2 - 1)} <= Theta_j := the same at lam0 (decreasing in lam), certified with a tail bound;
  - B <= 2 sqrt q K_P Theta_{d-2} lam^{d-1} e^{-pi lam^2};
    A <= 2 sqrt q [C'_P Theta_d + K_P Theta_{d-2}/(2 lam0^2) + K_Q Theta_{dQ-2} lam0^{dQ - d - 2}...] lam^{d+1} e^{-pi lam^2}
    (dQ = d + 2; C'_P = sum_k |p_k| (2 pi)^{k/2} lam0^{k-d});
  - zero sum <= 2 A^2 I(t_low) (I decreasing, T_L >= t_low).
Check that can fail: K_H c^{d+1} e^{-c} >= B_cert(x) at every certified grid point (the chain above dominates the
certified one term by term).

Usage: .venv/bin/python scripts/research_r2_rigorous/hermite_rate.py --json data/research_r2_rigorous/hermite_rate.json
"""
import argparse
import json
import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import cert_lib as C  # noqa: E402
from flint import arb  # noqa: E402


def sup_str(b, digits=8):
    """Decimal string >= the upper endpoint of b (pad by 10^-(digits-1) relative before rounding)."""
    v = C.up(b)
    return (v * (1 + arb(10) ** -(digits - 1))).str(digits, radius=False)


def theta(j, lam0, nmax=30):
    """sum_{n>=1} n^j e^{-pi lam0^2 (n^2-1)} (upper bound ball): explicit to nmax, geometric tail."""
    a = arb.pi() * lam0 * lam0
    tot = sum((arb(n) ** j * (-a * (n * n - 1)).exp() for n in range(1, nmax + 1)), arb(0))
    t = C.gauss_tail(a.exp(), max(j, 0), a, nmax)        # n^j <= n^max(j,0)
    return tot + arb(0, t)


def KR(R, lam0, Z02):
    dR = len(R) - 1
    tot = arb(0)
    for k, r in enumerate(R):
        if not r:
            continue
        g = 1 / (1 - arb(k - 1) / Z02) if k >= 2 else arb(1)
        tot += abs(r) * g * (2 * arb.pi()) ** (arb(k - 1) / 2) * lam0 ** (k - dR)
    return tot / (2 * arb.pi()).sqrt()


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--json", required=True)
    args = ap.parse_args()
    C.setprec()
    rows = json.load(open("data/research_r2_rigorous/hermite_cert.json"))["rows"]
    lam0 = arb(2).sqrt()
    Z02 = 4 * arb.pi()
    out = []
    groups = {}
    for r in rows:
        groups.setdefault((r["D"], r["s"]), []).append(r)
    nfail = 0
    ncheck = 0
    for (D, s), rs in sorted(groups.items()):
        q = 1 if D == 1 else abs(D)
        kappa = 0 if D > 0 else 1
        tr = C.HermiteTrial(kappa, s, D == 1)
        P, Q = tr.P, tr.Q
        d, dQ = len(P) - 1, len(Q) - 1
        assert Q[-1] != 0 and dQ == d + 2
        r0 = next(r for r in rs if float(r["xq"]) == 2.0)
        N0 = C.lo(arb(r0["norm2"].strip("[]").split("+/-")[0]) - arb(r0["norm2"].strip("[]").split("+/-")[1]))
        assert N0 > 0
        Kp, Kq = KR(P, lam0, Z02), KR(Q, lam0, Z02)
        Cp = sum((abs(p) * (2 * arb.pi()) ** (arb(k) / 2) * lam0 ** (k - d) for k, p in enumerate(P) if p), arb(0))
        sq = arb(q).sqrt()
        th = {j: theta(j, lam0) for j in (d - 2, d, dQ - 2)}
        KA = 2 * sq * (Cp * th[d] + Kp * th[d - 2] / (2 * lam0 * lam0) + Kq * th[dQ - 2])
        if D == 1:
            I = C.zeta_tail_integral(C.ZETA_T0)
        else:
            I = C.dirichlet_tail_integral(C.dirichlet_tlow(q), q, 1 if D > 0 else -1)
        KH = 2 * KA * KA * I / (N0 * (2 * arb.pi()) ** (d + 1))
        checks = []
        for r in rs:
            c = 2 * arb.pi() * arb(r["xq"])
            prop = KH * c ** (d + 1) * (-c).exp()
            ok = C.up(prop) >= arb(r["bound_upper"])          # must hold
            ncheck += 1
            nfail += not ok
            checks.append({"xq": r["xq"], "prop_bound": sup_str(prop, 6), "B_cert": r["bound_upper"],
                           "ratio": (C.up(prop) / arb(r["bound_upper"])).str(4, radius=False), "ok": bool(ok)})
        row = {"D": D, "q": q, "s": s, "d": d, "x0": 2 * q, "N0": N0.str(10, radius=False), "K_A": sup_str(KA),
               "I_tlow": sup_str(I), "K_H": sup_str(KH), "checks": checks}
        out.append(row)
        rmin = min(float(ch["ratio"]) for ch in checks)
        print(f"D={D:>4} s={s} d={d}: lambda_s(x) <= {row['K_H']} c^{d + 1} e^-c for x >= {2 * q}; "
              f"prop/B_cert >= {rmin:.3g} over {len(checks)} grid points")
    print(f"check K_H c^(d+1) e^-c >= B_cert: {ncheck - nfail}/{ncheck}")
    with open(args.json, "w") as fh:
        json.dump({"rows": out, "check": [ncheck - nfail, ncheck]}, fh, indent=1)


if __name__ == "__main__":
    main()
