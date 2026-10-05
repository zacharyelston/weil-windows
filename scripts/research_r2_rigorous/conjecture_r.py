#!/usr/bin/env python3
"""Pre-register R targets; compute falsifiable finite gates, including full norm upper.

The upper norm reuses committed certified main-integral balls. It reconstructs
their original dyadic domains, bounds the remainder upward, and pays for gaps.
It is not an independent quadrature rerun and is not an asymptotic norm proof.
"""
import argparse
from fractions import Fraction
import hashlib
import json
import math
from pathlib import Path
import subprocess

from flint import arb
import cert_lib as C

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "data/research_r2_rigorous"
DOC = "docs/R2_CONJECTURE_R.md"


def registration():
    hank = {}
    for k in range(2, 6):
        coeff = [Fraction(math.factorial(k+j), 2**j * math.factorial(j) * math.factorial(k-j)) for j in range(k+1)]
        odd_fact = math.prod(range(1, 2*k+2, 2))
        crossover_ratio = odd_fact * sum((a / 100**(k+1+j) for j, a in enumerate(coeff)), Fraction(0))
        assert sum(coeff) < 10000 and crossover_ratio < 1
        hank[str(k)] = {"sum": str(sum(coeff)), "hank100_over_cap": str(crossover_ratio)}
    return {
        "id": "R-REG-INITIAL",
        "status": "M1 and final rate are unproved targets; constants fixed before gates",
        "x0_over_q": 2,
        "C1_over_sqrt_q": 10**15, "p1": 8,
        "C2": "10^12 exp(4 pi q)", "p2": 8,
        "C": "20 * 10^42 * q * (1 + log(3q)) * exp(4 + 4 pi q)", "p": 24,
        "a_over_c": 10,
        "split_over_c": 100,
        "hank_coefficient_sum_upper": 10000,
        "cap_weight_upper": 6,
        "envelope_single_upper": "7 * 10^10 * c^5",
        "term_coefficient_sum_upper": "41 c",
        "J_upper": "3 * 10^12 * c^6",
        "edge_upper": "5 * 10^5 * sqrt(q) * c^(7/2)",
        "count_K": "5 * (1 + log(3q))",
        "fixed_hank_checks": hank,
        "formal_even_difference": {"constant": "1/8", "r": "-1/4", "zeta_minus_a_limit": "1/8"},
        "sheet_fit_range": [3.45, 11.5],
        "proof_constants": {"beta_min": 10, "I_seven_half_polynomial_lower": "2/5",
                            "split_cap_endpoint": 100, "p_max": 4, "k_min": 2, "k_max": 5},
    }


def proof_metadata():
    """Exact rational coefficients and constants used in the partial proof.

    Formal asymptotic coefficients are not explicit remainder estimates.
    This output does not alter the frozen registration or test its constants.
    """
    half = {}
    for k in (2, 3):
        half[str(k)] = [str(Fraction((-1)**j * math.factorial(k+j),
                                    2**j * math.factorial(j)*math.factorial(k-j))) for j in range(k+1)]
    corr = Fraction(4*2**2-1, 8)
    # r=2*pi*y^2, beta=b, c=b+2. Contributions from the exponent,
    # prefactor and first I2 / terminating I_(5/2) correction.
    W1 = [-corr, Fraction(1)-Fraction(3,4), -Fraction(1,8)]
    F1 = [Fraction(1)-Fraction(3), -Fraction(1)+Fraction(3,2), -Fraction(1,8)]
    diff = [a-b for a,b in zip(W1,F1)]
    assert diff[0] == Fraction(1,8) and diff[1] == -Fraction(1,4) and diff[2] == 0
    return {
        "status": "M2 proved for bounded-split exact envelope expressions; M1 and rate unproved",
        "stop_reason": "No explicit uniform normalized, pole-cancelled lattice remainder established from read sources",
        "half_integer_e_plus_polynomials": half,
        "formal_W1": [str(x) for x in W1], "formal_F1": [str(x) for x in F1],
        "formal_difference": [str(x) for x in diff], "formal_a_zeta_minus_limit": str(diff[0]),
        "scalar_bounds": {"I2_upper": 1, "I7half_lower": "2/5", "a_upper": 10,
                          "single_G_upper": "7e10", "J_upper": "3e12", "edge_upper": "5e5",
                          "count_K_factor": 5, "stieltjes_factor": 4,
                          "conditional_p": 2*registration()["p1"]+registration()["p2"]},
    }


def source_rows():
    for name in ("kb_cert", "kb_cert_supports"):
        for row in json.loads((DATA / (name + ".json")).read_text())["rows"]:
            selected = [t for t in row["candidates"] if t["dbeta"] == 2 and t["m"] == 2]
            assert len(selected) == 1
            yield name, row, selected[0]


def parameters(row):
    x = arb(row["support_2a"]).exp() if row["support_2a"] else arb(row["xq"]) * row["q"]
    c = 2 * arb.pi() * x / row["q"]
    trial = C.KBTrial(2, c-2, (x/row["q"]).sqrt(), int(row["D"] < 0), row["s"], row["D"] == 1)
    return x, c, trial


def norm_upper(row, cand, x, c, trial):
    """Return an upper endpoint for 2 int_1^sqrt(x) S(w)^2 dw.

    Original integration endpoints are reconstructed at original precision.
    Stored I intervals enclose the integral on those domains (trusting the
    original certificate). Their printed endpoint intervals must contain the
    reconstruction. A global bound on S pays for all omitted dyadic gaps.
    """
    rx, sq = x.sqrt(), arb(row["q"]).sqrt()
    covered, subtotal = arb(0), arb(0)
    previous_a = None
    K = math.floor(float(rx.mid()))
    if arb(K) > rx:
        K -= 1
    expected_bands = []
    for k in range(1, K+1):
        lb = rx/(k+1)
        assert lb > 1 or lb <= 1
        a = C.up(lb) if lb > 1 else arb(1)
        if C.lo(rx/k) > a:
            expected_bands.append(k)
    # At an integral sqrt(x), the final formal band has zero length and the
    # certificate correctly omits it. No positive-length band may be absent.
    assert [p["k"] for p in cand["pieces"]] == expected_bands
    for piece in cand["pieces"]:
        k, M = piece["k"], piece["M"]
        lb = rx / (k+1)
        assert lb > 1 or lb <= 1
        a = C.up(lb) if lb > 1 else arb(1)
        b = C.lo(rx/k)
        assert b > a and a >= 1 and b <= C.lo(rx)
        assert arb(piece["a"]).contains(a) and arb(piece["b"]).contains(b)
        if previous_a is not None:
            assert b <= previous_a
        previous_a = a
        length = b-a
        om = c*(k+M+1)*a/rx
        eps = C.up(trial.pref/2 * C.C_E(trial.env_terms(), trial.beta, om)
                   * (c*a/rx)**-3 / (2*arb(max(k+M, 1))**2))
        # Both stored and freshly computed epsilon upper endpoints are valid.
        # Taking their maximum absorbs textual rounding and reconstruction.
        eps = max(eps, C.up(arb(piece["eps_rest"])))
        I = C.up(arb(piece["int_Smain2"]))
        assert I >= 0
        subtotal += (I.sqrt() + eps*length.sqrt())**2
        covered += length
    # All reconstructed endpoints are exact dyadics. This deliberately covers
    # an extra sliver up to upper(sqrt(x)), in addition to the actual gaps.
    gap = C.up(rx)-1-covered
    assert gap >= 0
    N = math.ceil(float(C.up(rx)))+1
    assert arb(N) > C.up(rx)
    im = trial.beta.bessel_i(2)
    # |phi| <= (1+|a|) I2(beta); |phihat| <= 2 lambda sup|phi|.
    psi_sup = (1+abs(trial.a))*im*(1+2*trial.lam)/2
    h = c/rx
    tail = trial.pref/2 * C.C_E(trial.env_terms(), trial.beta, h*(N+1)) * h**-3/(2*arb(N)**2)
    S_sup = C.up(N*psi_sup + tail)
    gap_payment = gap*S_sup**2
    total = C.up(2*(subtotal+gap_payment))
    assert total >= C.lo(arb(cand["norm2"]))
    return total, C.up(gap), C.up(2*gap_payment)


def gates():
    reg = json.loads((DATA / "conjecture_r_registration.json").read_text())
    assert reg == registration(), "registration changed"
    committed = subprocess.check_output(["git", "show", f"HEAD:{DOC}"], cwd=ROOT, text=True)
    assert "R-REG-INITIAL" in committed, "commit registration before running gates"
    committed_reg = subprocess.check_output(["git", "show", "HEAD:data/research_r2_rigorous/conjecture_r_registration.json"], cwd=ROOT, text=True)
    assert json.loads(committed_reg) == reg
    out = {"registration": reg["id"], "registration_commit": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip(),
           "scope": "finite necessary tests; no proof of M1 or asymptotic rate",
           "norm_method": "stored certified integral upper endpoints, reconstructed dyadic intervals, plus remainder and gap bounds",
           "source_sha256": {name: hashlib.sha256((DATA/(name+".json")).read_bytes()).hexdigest() for name in ("kb_cert", "kb_cert_supports")},
           "cert_lib_sha256": hashlib.sha256((ROOT/"scripts/research_r2_rigorous/cert_lib.py").read_bytes()).hexdigest(),
           "rows": [], "stopped_on_failure": False}
    for name, row, cand in source_rows():
        assert row["prec"] == 256
        C.setprec(row["prec"])
        x, c, trial = parameters(row)
        assert x >= reg["x0_over_q"]*row["q"]
        q = arb(row["q"])
        C1 = arb(reg["C1_over_sqrt_q"])*q.sqrt()
        C2 = arb(10)**12 * (4*arb.pi()*q).exp()
        Crate = 20*arb(10)**42*q*(1+(3*q).log())*(4+4*arb.pi()*q).exp()
        H = C1*c**reg["p1"]
        proposed_norm = (2*trial.beta).exp()*c**(-reg["p2"])/C2
        proposed_rate = Crate*c**reg["p"]*(-2*c).exp()
        result = {"file": name, "D": row["D"], "s": row["s"], "x": x.str(25), "c": c.str(25),
                  "dbeta": cand["dbeta"], "m": cand["m"], "A_upper": C.up(arb(cand["A"])).str(25),
                  "B_upper": C.up(arb(cand["B"])).str(25), "proposed_AB_lower": C.lo(H).str(25),
                  "G1": bool(C.lo(H) >= C.up(arb(cand["A"])) and C.lo(H) >= C.up(arb(cand["B"]))) }
        if result["G1"]:
            nu, gap, payment = norm_upper(row, cand, x, c, trial)
            result.update({"actual_norm_upper": nu.str(25), "gap_length_upper": gap.str(25), "gap_norm_payment_upper": payment.str(25),
                           "proposed_norm_upper": C.up(proposed_norm).str(25), "G2": bool(C.up(proposed_norm) <= nu)})
        if result.get("G2"):
            result.update({"proposed_rate_lower": C.lo(proposed_rate).str(25), "candidate_bound_upper": cand["bound_upper"],
                           "G3": bool(C.lo(proposed_rate) >= arb(cand["bound_upper"]))})
        out["rows"].append(result)
        if not all(result.get(g, False) for g in ("G1", "G2", "G3")):
            out["stopped_on_failure"] = True
            break
    out["gates"] = {g: {"tested": sum(g in r for r in out["rows"]), "passed": sum(r.get(g, False) for r in out["rows"])} for g in ("G1", "G2", "G3")}
    (DATA/"conjecture_r_gates.json").write_text(json.dumps(out, indent=2)+"\n")
    print(json.dumps(out["gates"]))
    if out["stopped_on_failure"]:
        raise SystemExit("Gate failed or was inconclusive; stop, do not tune constants.")


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("mode", choices=("register", "gates", "proof"))
    args = ap.parse_args()
    if args.mode == "register":
        (DATA/"conjecture_r_registration.json").write_text(json.dumps(registration(), indent=2)+"\n")
        print(DATA/"conjecture_r_registration.json")
    elif args.mode == "gates":
        gates()
    else:
        (DATA/"conjecture_r_proof.json").write_text(json.dumps(proof_metadata(), indent=2)+"\n")
        print(DATA/"conjecture_r_proof.json")
