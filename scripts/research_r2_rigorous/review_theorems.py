#!/usr/bin/env python3
"""Reproduce the numerical evidence accompanying REVIEW_R2_THEOREMS.md.

This audits stored intervals and gives an exact finite-dimensional example.
It does not independently reproduce the certificate quadratures or zero census.
"""
import json
from fractions import Fraction
from pathlib import Path

from flint import arb
import cert_lib as C
import hermite_rate as H

ROOT = Path(__file__).resolve().parents[2]


def generate():
    C.setprec()
    all_trials = []
    for name in ("hermite_cert", "hermite_cert_supports", "kb_cert", "kb_cert_supports"):
        for row in json.loads((ROOT / f"data/research_r2_rigorous/{name}.json").read_text())["rows"]:
            for trial in row.get("candidates", [row]):
                all_trials.append((row, trial))
    bad_range = []
    bad_norm = []
    for row, trial in all_trials:
        if not arb(trial["norm2"]) > 0:
            bad_norm.append([row["D"], row["parity"], row["x"]])
        if row["D"] != 1:
            # The actual lower quotient endpoint, not a midpoint, must be in range.
            tc = C.up(arb(trial["A"])) / C.up(arb(trial["B"]))
            if not C.lo(tc) > arb(5) / 7:
                bad_range.append([row["D"], row["parity"], row["x"]])
    # Repair Proposition H's use of t_low without relying on a zero-free
    # theorem below BMOR's range: prove the analytic A_H/B_H is above t_low
    # for all lambda>=lambda0 using its increasing lambda^2 factor.
    uniform_cutoffs = []
    for row in json.loads((ROOT / "data/research_r2_rigorous/hermite_rate.json").read_text())["rows"]:
        D, s, q = row["D"], row["s"], row["q"]
        if D == 1:
            continue
        tr = C.HermiteTrial(int(D < 0), s, False)
        d = len(tr.P)-1
        lam0 = arb(2).sqrt()
        kp, kq = H.KR(tr.P, lam0, 4*arb.pi()), H.KR(tr.Q, lam0, 4*arb.pi())
        cp = sum((abs(p)*(2*arb.pi())**(arb(k)/2)*lam0**(k-d) for k,p in enumerate(tr.P) if p), arb(0))
        theta_low, theta_d = H.theta(d-2, lam0), H.theta(d, lam0)
        ka = 2*arb(q).sqrt()*(cp*theta_d + kp*theta_low/(2*lam0*lam0) + kq*theta_d)
        kb = 2*arb(q).sqrt()*kp*theta_low
        ratio = C.lo(ka*lam0*lam0/kb)
        threshold = C.up(C.dirichlet_tlow(q))
        uniform_cutoffs.append({"D": D, "s": s, "quotient_lower": ratio.str(25),
                                "threshold_upper": threshold.str(25), "passed": bool(ratio >= threshold)})
    # Q(v)=v_0^2+100v_1^2. Project (1,1) onto the second coordinate.
    rq = Fraction(1 + 100, 2)
    projected = Fraction(100)
    trial_upper = Fraction(60)
    assert rq <= trial_upper < projected
    return {
        "scope": "stored interval audit, not an independent rerun",
        "dirichlet_range_failures": len(bad_range),
        "nonpositive_norm_certificates": len(bad_norm),
        "bad_range_rows": bad_range,
        "bad_norm_rows": bad_norm,
        "hermite_uniform_cutoff_failures": sum(not r["passed"] for r in uniform_cutoffs),
        "hermite_uniform_cutoffs": uniform_cutoffs,
        "projection_example": {"diagonal": [1, 100], "vector": [1, 1],
                               "rq": float(rq), "valid_trial_upper": float(trial_upper),
                               "projected_rq": float(projected)},
        "mellin_counterexample": {"q": 5, "delta": 2, "imag_z": 0.75,
                                  "exponent_at_origin": -1.25},
        "sources": {"BW_second": [0.11200, 0.12567, 3.77417],
                    "Platt_S": 2.5167, "HSW_corrected_C3": 9.4925,
                    "BMOR_min_height": "5/7", "BMOR_zero_free_ell": 1.567},
    }


if __name__ == "__main__":
    out = ROOT / "data/research_r2_rigorous/review_theorems.json"
    out.write_text(json.dumps(generate(), indent=2) + "\n")
    print(out)
