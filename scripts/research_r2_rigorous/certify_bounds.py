#!/usr/bin/env python3
"""Certified (Arb) evaluation of R2's unconditional upper bound B_cert(x) >= lambda_s(x) (docs/R2_THEOREMS.md, Thms 2-4).

For each (L-function, x, sector) and construction (Hermite or Kaiser-Bessel trial), computes balls for A, B (upper
endpoints used), ||f_win||^2 (lower endpoint used) and the zero-count integral, and the final ball
    B_cert = 2 A^2 int_{T_L}^inf N+(t) t^-3 dt / ||f_win||^2,     reported upper endpoint = the certified bound.
Also records two checks that can fail: Theorem 1(i) (E(1/w) - sigma E(w) must enclose 0) and, for zeta, psi(0) = 0.

Usage:
  .venv/bin/python scripts/research_r2_rigorous/certify_bounds.py --construction hermite --d 1,5 --xq 2,4 --json out.json
  .venv/bin/python scripts/research_r2_rigorous/certify_bounds.py --construction kb --d 1,5 --xq 2,4 --dbeta 0,2 --json out.json
  .venv/bin/python scripts/research_r2_rigorous/certify_bounds.py --construction kb --zeta-supports --json out.json
"""
import argparse
import glob
import json
import multiprocessing as mproc
import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
import cert_lib as C  # noqa: E402
from flint import arb, ctx  # noqa: E402
from progress import Job  # noqa: E402

SUPPORTS = ["1.6", "2.38", "2.6", "2.99"]          # 2a for the four certified supports (x = e^{2a})


def fmt_up(b, digits=12):
    """Decimal string >= the upper endpoint of b (rounding error of str() < 1e-11 relative, so pad by 1e-11)."""
    v = C.up(b)
    v = v * (1 + arb(10) ** -11) if v > 0 else v * (1 - arb(10) ** -11)
    return v.str(digits, radius=False)


def fmt_down(b, digits=12):
    v = C.lo(b)
    v = v * (1 - arb(10) ** -11) if v > 0 else v * (1 + arb(10) ** -11)
    return v.str(digits, radius=False)


def bstr(b):
    return b.str(15, radius=True)


def r2_reference(construction, D, xq, parity):
    files = ["data/research_r2/hermite_bound.json"] if construction == "hermite" else \
        ["data/research_r2/kb_bound.json", "data/research_r2/kb_bound_b.json"]
    for f in files:
        if not os.path.exists(f):
            continue
        for r in json.load(open(f)):
            if r["D"] == D and abs(float(r["xq"]) - float(xq)) < 1e-9 and r["parity"] == parity:
                return r
    return None


def task_point(t):
    C.setprec(t["prec"])
    D, s, cons = t["D"], t["s"], t["construction"]
    q = 1 if D == 1 else abs(D)
    kappa = 0 if D > 0 else 1
    pole = D == 1
    parity = ("even", "odd")[s]
    if t.get("support"):
        x = (arb(t["support"])).exp()
        xq_label = None
    else:
        x = arb(t["xq"]) * q if float(t["xq"]).is_integer() else arb(t["xq"]) * q
        xq_label = t["xq"]
    t0 = time.time()
    row = {"construction": cons, "D": D, "q": q, "parity": parity, "s": s, "x": bstr(x),
           "xq": xq_label, "support_2a": t.get("support"), "prec": t["prec"]}
    sigma = 1 if s == 0 else -1
    if cons == "hermite":
        tr = C.HermiteTrial(kappa, s, pole)
        Q = C.hermite_quantities(tr, D, x, s)
        Z, Tc, TL = C.zero_sum_bound(Q["A"], Q["B"], D)
        bound = Z / C.lo(Q["norm2"])
        w0 = float(x.mid()) ** 0.25
        sym, e1 = C.hermite_symmetry(tr, D, x, s, w0)
        row.update({"index": tr.index, "poly": tr.P, "A": bstr(Q["A"]), "B": bstr(Q["B"]), "norm2": bstr(Q["norm2"]),
                    "edge": bstr(Q["edge"]), "Tc": bstr(Tc), "T_L": bstr(TL), "zero_sum": bstr(Z), "bound": bstr(bound),
                    "bound_upper": fmt_up(bound), "root_free_certified": Q["signed_integrals"],
                    "sym_w0": w0, "sym_diff": bstr(sym), "sym_contains_0": bool(sym.contains(0)),
                    "sym_rel": fmt_up(abs(sym) / abs(e1), 4) if not (e1 == 0) else None})
    else:
        c = 2 * arb.pi() * x / q
        lam = (x / q).sqrt()
        cands = []
        for db in t["dbeta"]:
            assert db >= 0
            tr = C.KBTrial(t["m"], c - db, lam, kappa, s, pole)
            Q = C.kb_quantities(tr, D, x, s)
            Z, Tc, TL = C.zero_sum_bound(Q["A"], Q["B"], D)
            bound = Z / C.lo(Q["norm2"])
            w0 = float(x.mid()) ** 0.25
            sym, e1 = C.kb_symmetry(tr, D, x, s, w0)
            cand = {"dbeta": db, "beta": bstr(tr.beta), "m": t["m"], "a_t2": bstr(tr.a), "A": bstr(Q["A"]), "B": bstr(Q["B"]),
                    "norm2": bstr(Q["norm2"]), "edge": bstr(Q["edge"]), "J0": bstr(Q["J0"]), "J1": bstr(Q["J1"]),
                    "pieces": Q["pieces"], "Tc": bstr(Tc), "T_L": bstr(TL), "zero_sum": bstr(Z), "bound": bstr(bound),
                    "bound_upper": fmt_up(bound), "sym_w0": w0, "sym_diff": bstr(sym),
                    "sym_contains_0": bool(sym.contains(0)), "sym_rel": fmt_up(abs(sym) / abs(e1), 4)}
            if pole:
                p0 = tr.psi0()
                cand.update({"psi0": bstr(p0), "psi0_contains_0": bool(p0.contains(0))})
            cands.append((bound, cand))
        best = min(cands, key=lambda bc: float(C.up(bc[0]).mid()))
        row.update({"index": kappa + 2 * s + 4 * pole, "candidates": [cd for _, cd in cands], "best_dbeta": best[1]["dbeta"]})
        row.update({k: best[1][k] for k in ("A", "B", "norm2", "Tc", "T_L", "zero_sum", "bound", "bound_upper",
                                            "sym_contains_0", "sym_rel")})
        if pole:
            row["psi0_contains_0"] = all(cd["psi0_contains_0"] for _, cd in cands)
    if xq_label is not None:
        ref = r2_reference(cons, D, xq_label, parity)
        if ref:
            row["r2_bound_float"] = ref["bound"]
            row["measured_ritz"] = ref.get("measured")
    row["seconds"] = round(time.time() - t0, 1)
    return row


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--construction", choices=("hermite", "kb"), required=True)
    ap.add_argument("--d", default="1")
    ap.add_argument("--xq", default="2")
    ap.add_argument("--zeta-supports", action="store_true", help="zeta at x = e^{1.6}, e^{2.38}, e^{2.6}, e^{2.99}")
    ap.add_argument("--dbeta", default="0,2")
    ap.add_argument("--m", type=int, default=2)
    ap.add_argument("--prec", type=int, default=C.PREC)
    ap.add_argument("--workers", type=int, default=1)
    ap.add_argument("--json", required=True)
    args = ap.parse_args()
    assert args.workers <= 4
    tasks = []
    base = {"construction": args.construction, "prec": args.prec, "m": args.m, "dbeta": [int(v) for v in args.dbeta.split(",")]}
    if args.zeta_supports:
        for sup in SUPPORTS:
            for s in (0, 1):
                tasks.append(dict(base, D=1, s=s, support=sup))
    else:
        for D in [int(v) for v in args.d.split(",")]:
            for xq in args.xq.split(","):
                for s in (0, 1):
                    tasks.append(dict(base, D=D, s=s, xq=xq))
    job = Job(f"r2 certify {args.construction}", total=len(tasks), args=vars(args))
    rows = []
    with mproc.Pool(args.workers) as pool:
        for row in pool.imap_unordered(task_point, tasks):
            rows.append(row)
            rows.sort(key=lambda r: (r["D"], float(r["xq"]) if r["xq"] else 0, r.get("support_2a") or "", r["s"]))
            job.step(f"D={row['D']} x={row['x'][:20]} {row['parity']}: B_cert <= {row['bound_upper']} "
                     f"(R2 float {row.get('r2_bound_float')}) sym0={row['sym_contains_0']} ({row['seconds']}s)")
            with open(args.json, "w") as fh:
                json.dump({"construction": args.construction, "args": vars(args), "flint": __import__("flint").__version__,
                           "rows": rows}, fh, indent=1)
    job.done()


if __name__ == "__main__":
    main()
