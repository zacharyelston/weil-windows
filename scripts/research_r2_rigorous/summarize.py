#!/usr/bin/env python3
"""Validation verdicts and the tables of docs/R2_THEOREMS.md, printed from the committed JSON only.

Checks that can fail (never two upper bounds against each other):
  (a) RQ(trial) <= B_cert   at every point where RQ was computed     (data/research_r2_rigorous/rq_check*.json)
  (b) B_cert >= certified lower bound at the four certified supports  (both sectors, both constructions)
  (c) RQ(trial) >= certified lower bound at the four supports        (an upper bound against a lower bound)
  (d) Theorem 1(i) symmetry ball encloses 0; psi(0) ball encloses 0 for zeta; zero-count checks.
Comparisons that are NOT validations (reported for orientation only): B_cert / R2's floating bound; B_cert / Ritz minimum.

Usage: .venv/bin/python scripts/research_r2_rigorous/summarize.py [--json data/research_r2_rigorous/summary.json]
"""
import argparse
import glob
import json
import math
import os

import mpmath as mp

mp.mp.dps = 40
D = "data/research_r2_rigorous"

# certified LOWER bounds (rh2 certificates; docs/CERTIFICATE_238.md, data/connes/*certificate*.json)
CERT = {
    "1.6": {"even": ("1.0276895590199815698e-17", "data/connes/zhu_certificate_even.json"),
            "odd": ("9.1183384453799998655e-15", "data/connes/zhu_certificate_odd.json")},
    "2.38": {"even": ("6.8131164696007119625e-48", "data/connes/grid_certificate_a119_even.json"),
             "odd": ("3.9196686160954221947e-44", "data/connes/grid_certificate_a119_odd.json")},
    "2.6": {"even": ("5.7634447922293648458e-62", "data/connes/grid_certificate_a130_even.json"),
            "odd": ("5.3601130404895638862e-58", "data/connes/grid_certificate_a130_odd.json")},
    "2.99": {"even": ("3.5011421764109500236e-96", "data/connes/grid_certificate_a1495_even.json"),
             "odd": ("8.2562649400193819720e-92", "data/connes/grid_certificate_a1495_odd.json")},
}
# rh2's best finite-basis UPPER bounds at the supports (data/connes/decay/fit_zeta.json) -- orientation only
RITZ = {("1.6", "even"): "1.63558583932233e-17", ("1.6", "odd"): "1.54916134936799e-14",
        ("2.38", "even"): "1.01419882049612e-47", ("2.38", "odd"): "5.64114073526898e-44",
        ("2.6", "even"): "9.12455269172982e-62", ("2.6", "odd"): "8.42951149251576e-58",
        ("2.99", "even"): "6.26578516781e-96", ("2.99", "odd"): "1.47929992461e-91"}
NAME = {1: "ζ", 5: "χ₅", 8: "χ₈", 12: "χ₁₂", 13: "χ₁₃", 17: "χ₁₇", -3: "χ₋₃", -4: "χ₋₄", -7: "χ₋₇", -8: "χ₋₈",
        -11: "χ₋₁₁", -15: "χ₋₁₅", -20: "χ₋₂₀"}


def lower_from_ball(s):
    """Lower endpoint of an Arb ball string '[mid +/- rad]' (for the certified lower bounds stored as balls)."""
    s = s.strip("[]")
    mid, rad = s.split("+/-")
    return mp.mpf(mid) - mp.mpf(rad)


def load(name):
    p = os.path.join(D, name)
    return json.load(open(p))["rows"] if os.path.exists(p) else []


def sci(v, d=3):
    return mp.nstr(mp.mpf(v), d, min_fixed=0, max_fixed=0)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--json", default=os.path.join(D, "summary.json"))
    args = ap.parse_args()
    out = {}
    # certified lower bounds, re-read from the certificate files where they are stored as balls
    certlow = {}
    for sup, d in CERT.items():
        for par, (val, f) in d.items():
            v = mp.mpf(val)
            if os.path.exists(f):
                j = json.load(open(f))
                if "certified_lower_bound" in j:
                    v = lower_from_ball(j["certified_lower_bound"])
                else:
                    v = mp.mpf(max(j["certificates"], key=lambda c: mp.mpf(c["Q_lower"]))["Q_lower"])
            certlow[(sup, par)] = v
    hs, ks = load("hermite_cert_supports.json"), load("kb_cert_supports.json")
    rq = []
    for f in sorted(glob.glob(os.path.join(D, "rq_check*.json"))):
        rq += json.load(open(f))["rows"]

    def rq_of(cons, Dv, par, xq=None, sup=None):
        hits = [r for r in rq if r["construction"] == cons and r["D"] == Dv and r["parity"] == par and
                (r.get("support_2a") == sup if sup else (r.get("xq") is not None and abs(float(r["xq"]) - float(xq)) < 1e-9))
                and r.get("n_out", 10) == 10]
        return hits[0] if hits else None

    # ---------------- table: zeta at the certified supports
    print("\n### Table S: ζ at the four certified supports (B_cert = certified upper bound; cert = certified lower bound)\n")
    print("| support 2a | x | sector | cert (lower) | Hermite B_cert | KB B_cert | KB β | KB RQ (a) | Hermite RQ (a) | KB B_cert/cert |")
    print("|---|---|---|---|---|---|---|---|---|---|")
    rows_s = []
    nb, nb_ok, nc, nc_ok = 0, 0, 0, 0
    for sup in ("1.6", "2.38", "2.6", "2.99"):
        for par in ("even", "odd"):
            h = next(r for r in hs if r["support_2a"] == sup and r["parity"] == par)
            k = next(r for r in ks if r["support_2a"] == sup and r["parity"] == par)
            cl = certlow[(sup, par)]
            hb, kb_ = mp.mpf(h["bound_upper"]), mp.mpf(k["bound_upper"])
            okb = [hb >= cl, kb_ >= cl]
            nb += 2
            nb_ok += sum(okb)
            rk, rh = rq_of("kb", 1, par, sup=sup), rq_of("hermite", 1, par, sup=sup)
            okc = []
            for r in (rk, rh):
                if r:
                    nc += 1
                    okc.append(mp.mpf(r["rq"]) >= cl)
                    nc_ok += okc[-1]
            beta = "c" if k["best_dbeta"] == 0 else f"c−{k['best_dbeta']}"
            x = math.exp(float(sup))
            rows_s.append({"support_2a": sup, "x": x, "parity": par, "cert_lower": mp.nstr(cl, 12), "hermite_B_cert": h["bound_upper"],
                           "kb_B_cert": k["bound_upper"], "kb_dbeta": k["best_dbeta"], "b_ok": okb,
                           "kb_rq": rk["rq"] if rk else None, "hermite_rq": rh["rq"] if rh else None, "c_ok": okc})
            print(f"| {sup} | {x:.4f} | {par} | {sci(cl, 5)} | {sci(hb, 4)} | {sci(kb_, 4)} | {beta} | "
                  f"{sci(rk['rq'], 3) if rk else '—'} | {sci(rh['rq'], 3) if rh else '—'} | {sci(kb_ / cl, 3)} |")
    out["supports"] = rows_s
    print(f"\n(b) B_cert >= certified lower bound: {nb_ok}/{nb}.   (c) RQ >= certified lower bound: {nc_ok}/{nc}.")
    out["check_b"] = [nb_ok, nb]
    out["check_c"] = [nc_ok, nc]
    # ---------------- (a): RQ <= B_cert everywhere computed
    na, na_ok, stable = 0, 0, 0
    worst = None
    fails = []
    for r in rq:
        if r.get("n_out", 10) != 10:
            continue
        na += 1
        ok = mp.mpf(r["rq"]) <= mp.mpf(r["bound_upper"])
        na_ok += ok
        stable += bool(r["stable"])
        ratio = mp.mpf(r["bound_upper"]) / mp.mpf(r["rq"]) if mp.mpf(r["rq"]) > 0 else None
        if ratio is not None and (worst is None or ratio < worst[0]):
            worst = (ratio, r)
        if not ok:
            fails.append(r)
    byc = {}
    for r in rq:
        if r.get("n_out", 10) != 10:
            continue
        cons = r["construction"]
        b = byc.setdefault(cons, [0, 0, None, None])
        b[0] += 1
        ok = mp.mpf(r["rq"]) <= mp.mpf(r["bound_upper"])
        b[1] += ok
        ratio = mp.mpf(r["bound_upper"]) / mp.mpf(r["rq"])
        b[2] = ratio if b[2] is None or ratio < b[2] else b[2]
        b[3] = ratio if b[3] is None or ratio > b[3] else b[3]
    print(f"\n(a) RQ(trial) <= B_cert: {na_ok}/{na} (precision-stable: {stable}/{na})")
    for cons, b in byc.items():
        print(f"    {cons}: {b[1]}/{b[0]}, B_cert/RQ in [{sci(b[2])}, {sci(b[3])}]")
    if worst:
        w = worst[1]
        print(f"    tightest: {w['construction']} D={w['D']} {w.get('xq') or 'e^' + str(w.get('support_2a'))} {w['parity']}: B_cert/RQ = {sci(worst[0])}")
    for r in fails:
        print(f"    FAIL: {r}")
    out["check_a"] = {"ok": na_ok, "n": na, "stable": stable, "by_construction": {k: [v[0], v[1], mp.nstr(v[2], 4), mp.nstr(v[3], 4)] for k, v in byc.items()},
                      "fails": fails}
    conv = [r for r in rq if r.get("n_out", 10) != 10]
    if conv:
        print("\n    n_out convergence (KB):")
        for r in conv:
            base = rq_of(r["construction"], r["D"], r["parity"], xq=r.get("xq"), sup=r.get("support_2a"))
            if base:
                print(f"      D={r['D']} {r.get('xq') or r.get('support_2a')} {r['parity']}: RQ(n_out=10) {base['rq']}  RQ(n_out={r['n_out']}) {r['rq']}")
    # ---------------- decomposition (diagnostic, NOT a validation): B_cert/RQ = Theorem 2's slack on that trial;
    # RQ/lambda = the trial's own excess, bracketed by RQ/Ritz <= RQ/lambda <= RQ/cert at the supports
    print("\n### Table D: where the gap B_cert/λ comes from (supports; diagnostic only)\n")
    print("| support | sector | construction | B_cert/RQ (Thm 2 slack) | RQ/Ritz ≤ RQ/λ ≤ RQ/cert (trial excess) |")
    print("|---|---|---|---|---|")
    dec = []
    for sup in ("1.6", "2.38", "2.6", "2.99"):
        for par in ("even", "odd"):
            for cons, rows_ in (("kb", ks), ("hermite", hs)):
                r = next(r for r in rows_ if r["support_2a"] == sup and r["parity"] == par)
                q_ = rq_of(cons, 1, par, sup=sup)
                if not q_:
                    continue
                rqv = mp.mpf(q_["rq"])
                sl = mp.mpf(r["bound_upper"]) / rqv
                lo_, hi_ = rqv / mp.mpf(RITZ[(sup, par)]), rqv / certlow[(sup, par)]
                dec.append({"support_2a": sup, "parity": par, "construction": cons, "slack": mp.nstr(sl, 4),
                            "trial_excess_low": mp.nstr(lo_, 4), "trial_excess_high": mp.nstr(hi_, 4)})
                print(f"| {sup} | {par} | {cons} | {sci(sl)} | {sci(lo_)} .. {sci(hi_)} |")
    out["decomposition_supports"] = dec
    gdec = {}
    for r in rq:
        if r.get("n_out", 10) != 10 or not r.get("xq"):
            continue
        cert = next((c for c in load("kb_cert.json" if r["construction"] == "kb" else "hermite_cert.json")
                     if c["D"] == r["D"] and c["parity"] == r["parity"] and abs(float(c["xq"]) - float(r["xq"])) < 1e-9), None)
        if not cert or not cert.get("measured_ritz"):
            continue
        g = gdec.setdefault(r["construction"], {"slack": [], "excess": []})
        g["slack"].append(mp.mpf(r["bound_upper"]) / mp.mpf(r["rq"]))
        g["excess"].append(mp.mpf(r["rq"]) / mp.mpf(cert["measured_ritz"]))
    for cons, g in gdec.items():
        print(f"\n    grid, {cons}: Theorem-2 slack B_cert/RQ in [{sci(min(g['slack']))}, {sci(max(g['slack']))}]; "
              f"trial excess RQ/Ritz (a lower bound for RQ/λ) in [{sci(min(g['excess']))}, {sci(max(g['excess']))}] over {len(g['slack'])} points")
        out[f"decomposition_grid_{cons}"] = {"slack": [mp.nstr(min(g["slack"]), 4), mp.nstr(max(g["slack"]), 4)],
                                            "excess_vs_ritz": [mp.nstr(min(g["excess"]), 4), mp.nstr(max(g["excess"]), 4)], "n": len(g["slack"])}
    # ---------------- character / grid tables
    for cons, fname in (("hermite", "hermite_cert.json"), ("kb", "kb_cert.json")):
        rows = load(fname)
        if not rows:
            continue
        print(f"\n### Table {cons}: certified bounds on R2's grid (ranges over x/q)\n")
        print("| object | sector | n | x/q | B_cert range | B_cert / R2 float | ln B_cert + " + ("c" if cons == "hermite" else "2c") + " range | p of C c^p fit (max resid) | sym ⊇ 0 |")
        print("|---|---|---|---|---|---|---|---|---|")
        groups = {}
        for r in rows:
            groups.setdefault((r["D"], r["parity"]), []).append(r)
        summ = []
        for (Dv, par), rs in sorted(groups.items(), key=lambda kv: (abs(kv[0][0]) if kv[0][0] != 1 else 0, kv[0][0] < 0, kv[0][1])):
            rs.sort(key=lambda r: float(r["xq"]))
            q = 1 if Dv == 1 else abs(Dv)
            bs = [mp.mpf(r["bound_upper"]) for r in rs]
            rat = [mp.mpf(r["bound_upper"]) / mp.mpf(r["r2_bound_float"]) for r in rs if r.get("r2_bound_float")]
            k = 1 if cons == "hermite" else 2
            lnr = [mp.log(b) + k * 2 * mp.pi * float(r["xq"]) for b, r in zip(bs, rs)]
            sym = all(r["sym_contains_0"] for r in rs)
            summ.append({"D": Dv, "parity": par, "n_points": len(rs), "B_min": mp.nstr(min(bs), 4), "B_max": mp.nstr(max(bs), 4),
                         "ratio_to_r2": [mp.nstr(min(rat), 4), mp.nstr(max(rat), 4)] if rat else None, "sym_all": sym})
            # least-squares fit ln B_cert + k c = alpha + p ln c (the shape C c^p e^{-k c}); max residual
            X = [mp.log(k * 2 * mp.pi * float(r["xq"]) / k) for r in rs]
            Y = lnr
            n_ = len(X)
            mx, my = sum(X) / n_, sum(Y) / n_
            pfit = sum((a - mx) * (b - my) for a, b in zip(X, Y)) / sum((a - mx) ** 2 for a in X)
            afit = my - pfit * mx
            res = max(abs(b - afit - pfit * a) for a, b in zip(X, Y))
            summ[-1].update({"p_fit": mp.nstr(pfit, 4), "max_resid": mp.nstr(res, 3)})
            print(f"| {NAME.get(Dv, Dv)} | {par} | {rs[0]['index']} | {rs[0]['xq']}..{rs[-1]['xq']} | {sci(max(bs))} .. {sci(min(bs))} | "
                  f"{mp.nstr(min(rat), 3)} .. {mp.nstr(max(rat), 3)} | {mp.nstr(min(lnr), 3)} .. {mp.nstr(max(lnr), 3)} | {mp.nstr(pfit, 3)} ({mp.nstr(res, 2)}) | {sym} |")
        out[f"grid_{cons}"] = summ
        allr = [mp.mpf(r["bound_upper"]) / mp.mpf(r["r2_bound_float"]) for r in rows if r.get("r2_bound_float")]
        print(f"\n{cons}: {len(rows)} points; B_cert / R2 floating bound in [{mp.nstr(min(allr), 4)}, {mp.nstr(max(allr), 4)}]; "
              f"symmetry ball encloses 0 at {sum(r['sym_contains_0'] for r in rows)}/{len(rows)}"
              + ("" if cons == "hermite" else f"; ψ(0) ∋ 0 for ζ at {sum(bool(r.get('psi0_contains_0')) for r in rows if r['D'] == 1)}/{sum(r['D'] == 1 for r in rows)}"))
        if cons == "hermite":
            print(f"   root-free certificate used for every out-of-band integral: {sum(r['root_free_certified'] for r in rows)}/{len(rows)}")
    # ---------------- (e) Proposition H consistency
    hr = os.path.join(D, "hermite_rate.json")
    if os.path.exists(hr):
        h = json.load(open(hr))
        print(f"\n(e) Proposition H bound K_H c^(d+1) e^-c >= B_cert at every Hermite grid point: {h['check'][0]}/{h['check'][1]}")
        out["check_e"] = h["check"]
    # ---------------- KB family scan at the supports (certified; any member is a valid bound)
    var = [(2, "kb_cert_supports_m2_var.json"), (3, "kb_cert_supports_m3.json")]
    if all(os.path.exists(os.path.join(D, f)) for _, f in var):
        print("\n### Table V: certified KB bounds at the supports over (m, β = c − Δβ) (all [C]; Theorem C's table uses m = 2, Δβ ∈ {0, 2})\n")
        print("| support | sector | m = 2: Δβ = 0 | 1 | 2 | 3 | 4 | m = 3: Δβ = 0 | 2 | 4 | best / table value |")
        print("|---|---|---|---|---|---|---|---|---|---|---|")
        tab = {(r["support_2a"], r["parity"]): r for r in ks}
        v2 = {(r["support_2a"], r["parity"]): r for r in json.load(open(os.path.join(D, var[0][1])))["rows"]}
        v3 = {(r["support_2a"], r["parity"]): r for r in json.load(open(os.path.join(D, var[1][1])))["rows"]}
        vout = []
        for key in sorted(v2, key=lambda k: (float(k[0]), k[1])):
            c2 = {c["dbeta"]: mp.mpf(c["bound_upper"]) for c in v2[key]["candidates"]}
            c3 = {c["dbeta"]: mp.mpf(c["bound_upper"]) for c in v3[key]["candidates"]}
            best = min(list(c2.values()) + list(c3.values()))
            ratio = best / mp.mpf(tab[key]["bound_upper"])
            vout.append({"support_2a": key[0], "parity": key[1], "best": mp.nstr(best, 6), "best_over_table": mp.nstr(ratio, 4)})
            print(f"| {key[0]} | {key[1]} | " + " | ".join(sci(c2[d]) for d in (0, 1, 2, 3, 4)) + " | " +
                  " | ".join(sci(c3[d]) for d in (0, 2, 4)) + f" | {mp.nstr(ratio, 3)} |")
        out["kb_family_supports"] = vout
    # ---------------- tapered prolate experiment [N]
    tp = os.path.join(D, "prolate_taper.json")
    if os.path.exists(tp):
        trows = json.load(open(tp))
        print("\n### Table P: tapered sharp prolate [N] against the certified KB bound (ζ)\n")
        print("| support | sector | κ (ε = κ/c²) | taper RQ | RQ/Ritz | taper Thm-2 estimate [N] | KB B_cert [C] | estimate/KB | taper slack (estimate/RQ) | KB slack (B_cert/RQ) | T_c taper |")
        print("|---|---|---|---|---|---|---|---|---|---|---|")
        pout = []
        for r in sorted(trows, key=lambda r: (float(r["support_2a"]), r["parity"], float(r["kappa"]))):
            k = (r["support_2a"], r["parity"])
            kbr = next(x for x in ks if x["support_2a"] == k[0] and x["parity"] == k[1])
            kq = rq_of("kb", 1, k[1], sup=k[0])
            be, rqv, kbb = mp.mpf(r["bound_est"]), mp.mpf(r["rq"]), mp.mpf(kbr["bound_upper"])
            pout.append({"support_2a": k[0], "parity": k[1], "kappa": r["kappa"], "estimate_over_kb": mp.nstr(be / kbb, 4)})
            print(f"| {k[0]} | {k[1]} | {int(float(r['kappa']))} | {sci(rqv)} | {mp.nstr(rqv / mp.mpf(RITZ[k]), 3)} | {sci(be)} | {sci(kbb)} | "
                  f"{mp.nstr(be / kbb, 3)} | {sci(be / rqv)} | {sci(kbb / mp.mpf(kq['rq'])) if kq else '—'} | {mp.nstr(mp.mpf(r['Tc']), 5)} |")
        out["taper"] = pout
    with open(args.json, "w") as fh:
        json.dump(out, fh, indent=1, default=str)


if __name__ == "__main__":
    main()
