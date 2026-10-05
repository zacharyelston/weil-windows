#!/usr/bin/env python3
"""Validation (a): the trial function's own Rayleigh quotient in our basis must lie below the certified bound.

For each point, the E-map window function f_win of the trial (Hermite, and Kaiser-Bessel with the dbeta that gave the
best certified bound) is projected onto our N-mode basis of the sector (even: 1/sqrt L, sqrt(2/L) cos(2 pi k u/L),
k = 0..N; odd: sqrt(2/L) sin(2 pi k u/L), k = 1..N), and RQ = c^T M c / c^T c is evaluated with our explicit-formula
matrix M = decay_law_mp.zeros_side(D, x, N, parity) (archimedean + primes + pole; no zeros enter). RQ is a Rayleigh
quotient of an explicit vector, hence itself an upper bound for the basis minimum; it is NOT certified [N].
Precision-stable: everything is recomputed at dps and dps + 20 and must agree to 1e-6 (relative).

The check that can fail: RQ <= B_cert (upper endpoint, data/research_r2_rigorous/*_cert*.json). By Theorems 1-2,
Q(f_win)/||f_win||^2 <= B_cert exactly; RQ differs from Q(f_win)/||f_win||^2 only by the basis truncation (reported:
1 - |c|^2/||f_win||^2) and by dropping out-of-band terms n > k + n_out on each band interval (convergence in n_out
is reported at selected points).

Numerics. f_win = g + h on each band interval I_k = [sqrt x/(k+1), sqrt x/k] (in w = e^u): g = in-band terms (n <= k,
smooth, high precision, coarse Gauss-Legendre panels), h = out-of-band terms k < n <= k + n_out (oscillatory but of
relative size e^{-c}: lower precision on fine panels). Hermite trials need no split (Gaussian decay).

Usage: .venv/bin/python scripts/research_r2_rigorous/rq_check.py --cert data/research_r2_rigorous/kb_cert.json \
           --herm data/research_r2_rigorous/hermite_cert.json --nmax 240 --json data/research_r2_rigorous/rq_check.json
"""
import argparse
import glob
import json
import math
import multiprocessing as mproc
import os
import sys
import time

import mpmath as mp

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
import cert_lib as C  # noqa: E402
from flint import arb, arb_mat, ctx  # noqa: E402
from decay_law_mp import zeros_side  # noqa: E402
from progress import Job  # noqa: E402

SUPPORT_N = {"1.6": 60, "2.38": 100, "2.6": 120, "2.99": 180}


def scan_n2(D, xq, parity):
    for f in glob.glob(f"data/connes/decay/scan*_D{D}.json"):
        if f.split("/")[-1].split("_")[0] not in ("scanhi", "scanext", "scan8"):
            continue
        for r in json.load(open(f)):
            if abs(float(r["xq"]) - float(xq)) < 1e-9 and r["parity"] == parity:
                return r["n2"], r["dps"], r["lambda_n2"]
    return None, None, None


_GL = {}


def gl(m):
    key = (m, ctx.prec)
    if key not in _GL:
        _GL[key] = [arb.legendre_p_root(m, i, weight=True) for i in range(m)]
    return _GL[key]


def panels(a, b, hmax):
    n = max(1, int(math.ceil(float(((b - a) / hmax).mid()))))
    h = (b - a) / n
    return [(a + i * h, a + (i + 1) * h) for i in range(n)]


def project(nodes, vals, L, N, parity):
    """2 int_0^{L/2} f b_k = sum over nodes of weight * f * b_k (b_k as in native)."""
    pi2L = 2 * arb.pi() / L
    sL, s2L = 1 / L.sqrt(), (2 / L).sqrt()
    if parity == "even":
        c = [arb(0)] * (N + 1)
        for (u, wt), fv in zip(nodes, vals):
            g = 2 * wt * fv
            c[0] += g * sL
            th = pi2L * u
            c1, s1 = th.cos(), th.sin()
            ck, sk = arb(1), arb(0)
            for k in range(1, N + 1):
                ck, sk = ck * c1 - sk * s1, sk * c1 + ck * s1
                c[k] += g * s2L * ck
        return c
    c = [arb(0)] * N
    for (u, wt), fv in zip(nodes, vals):
        g = 2 * wt * fv
        th = pi2L * u
        c1, s1 = th.cos(), th.sin()
        ck, sk = arb(1), arb(0)
        for k in range(1, N + 1):
            ck, sk = ck * c1 - sk * s1, sk * c1 + ck * s1
            c[k - 1] += g * s2L * sk
    return c


def coefficients(cons, D, x, s, dbeta, N, n_out, m_gl=40):
    """Projection coefficients of f_win and ||f_win||^2 (quadrature) at the current ctx.prec."""
    chi = C.chi_fn(D)
    q = 1 if D == 1 else abs(D)
    sq = arb(q).sqrt()
    kappa = 0 if D > 0 else 1
    pole = D == 1
    parity = ("even", "odd")[s]
    L = x.log()
    rx = x.sqrt()
    hmode = L / (2 * N + 8)                       # <= half an oscillation of the top mode per panel
    if cons == "hermite":
        tr = C.HermiteTrial(kappa, s, pole)
        nodes, vals = [], []
        for a, b in panels(arb(0), L / 2, hmode):
            for t, wt in gl(m_gl):
                u = a + (t + 1) * (b - a) / 2
                w = u.exp()
                nmax = int(math.ceil(12 * math.sqrt(q) / float(w.mid()))) + 2
                S = sum((chi(n) * tr.psi(n * w / sq) for n in range(1, nmax + 1) if chi(n)), arb(0))
                nodes.append((u, wt * (b - a) / 2))
                vals.append((u / 2).exp() * S)
        cvec = project(nodes, vals, L, N, parity)
        norm2 = 2 * sum((wt * v * v for (u, wt), v in zip(nodes, vals)), arb(0))
        return cvec, norm2
    c = 2 * arb.pi() * x / q
    lam = (x / q).sqrt()
    tr = C.KBTrial(2, c - dbeta, lam, kappa, s, pole)
    K = int(math.floor(float(rx.mid())))
    if arb(K) > rx:
        K -= 1
    # band intervals in u: [log(sqrt x/(k+1)), log(sqrt x/k)] cap [0, L/2]
    ivs = []
    for k in range(1, K + 1):
        lo_u = (rx / (k + 1)).log() if (rx / (k + 1) > 1) else arb(0)
        hi_u = (rx / k).log()
        if hi_u > lo_u:
            ivs.append((k, lo_u, hi_u))
    gnodes, gvals = [], []
    hnodes, hvals = [], []
    prec_hi = ctx.prec
    for k, a0, b0 in ivs:
        cin = [(n, chi(n)) for n in range(1, k + 1) if chi(n)]
        cout = [(n, chi(n)) for n in range(k + 1, k + n_out + 1) if chi(n)]
        for a, b in panels(a0, b0, hmode):
            for t, wt in gl(m_gl):
                u = a + (t + 1) * (b - a) / 2
                w = u.exp()
                S = sum((cc * tr.psi_full(n * w / sq) for n, cc in cin), arb(0))
                gnodes.append((u, wt * (b - a) / 2))
                gvals.append((u / 2).exp() * S)
        # h: fastest phase rate on this interval ~ c (k + n_out) w/sqrt x <= c (k + n_out)/k; ~3 rad per panel, 20 nodes
        rate = float(c.mid()) * (k + n_out) / k
        hfine = min(float(hmode.mid()), 3.0 / rate)
        ctx.prec = 96
        nodes20 = gl(20)
        for a, b in panels(a0, b0, arb(hfine)):
            for t, wt in nodes20:
                u = a + (t + 1) * (b - a) / 2
                w = u.exp()
                S = sum((cc * tr.psi_out(n * w / sq) for n, cc in cout), arb(0))
                hnodes.append((u, wt * (b - a) / 2))
                hvals.append((u / 2).exp() * S)
        ctx.prec = prec_hi
    cg = project(gnodes, gvals, L, N, parity)
    ctx.prec = 96
    ch = project(hnodes, hvals, L, N, parity)
    ctx.prec = prec_hi
    cvec = [arb(a.mid()) + arb(b.mid()) for a, b in zip(cg, ch)]
    # ||f_win||^2 = 2 int (g + h)^2 ~ 2 int g^2 + 4 int g h + 2 int h^2: cross term on the fine grid
    n_gg = 2 * sum((wt * v * v for (u, wt), v in zip(gnodes, gvals)), arb(0))
    return cvec, n_gg


def rq_point(task):
    cons, D, s, dbeta, N, dps, n_out = task["cons"], task["D"], task["s"], task.get("dbeta", 0), task["N"], task["dps"], task["n_out"]
    parity = ("even", "odd")[s]
    t0 = time.time()
    out = {"construction": cons, "D": D, "parity": parity, "xq": task.get("xq"), "support_2a": task.get("support"),
           "N": N, "dbeta": dbeta, "n_out": n_out}
    rqs = []
    for d in (dps, dps + 20):
        mp.mp.dps = d
        ctx.prec = int(d * 3.33) + 30
        q = 1 if D == 1 else abs(D)
        x = arb(task["support"]).exp() if task.get("support") else arb(task["xq"]) * q
        xm = mp.mpf(x.mid().str(d + 10, radius=False))
        cvec, nf2 = coefficients(cons, D, x, s, dbeta, N, n_out)
        M = zeros_side(D, xm, N, parity)
        A = arb_mat([[arb(mp.nstr(M[i, j], d + 5)) for j in range(M.cols)] for i in range(M.rows)])
        cm = arb_mat([[v] for v in cvec])
        num = (cm.transpose() * A * cm)[0, 0]
        den = sum((v * v for v in cvec), arb(0))
        rqs.append((num / den, den / nf2))
    r1, r2 = rqs[0][0], rqs[1][0]
    rel = abs(r1 - r2) / abs(r2)
    out.update({"rq": r1.str(10, radius=False), "rq_hi": r2.str(10, radius=False), "rel_diff": float(rel.mid()),
                "stable": bool(float(rel.mid()) < 1e-6), "captured": float(rqs[0][1].mid()), "dps": dps,
                "seconds": round(time.time() - t0, 1)})
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--kb", default="data/research_r2_rigorous/kb_cert.json")
    ap.add_argument("--herm", default="data/research_r2_rigorous/hermite_cert.json")
    ap.add_argument("--kb-supports", default="data/research_r2_rigorous/kb_cert_supports.json")
    ap.add_argument("--herm-supports", default="data/research_r2_rigorous/hermite_cert_supports.json")
    ap.add_argument("--only", default="", help="comma list of D to include (default: all in the KB file)")
    ap.add_argument("--xq", default="", help="comma list of x/q to include (default: all)")
    ap.add_argument("--herm-extra", action="store_true", help="also the Hermite points absent from the KB grid")
    ap.add_argument("--supports", action="store_true")
    ap.add_argument("--kb-only", action="store_true", help="skip the Hermite trials (n_out convergence runs)")
    ap.add_argument("--grid", action="store_true")
    ap.add_argument("--nmax", type=int, default=240)
    ap.add_argument("--n-out", type=int, default=10)
    ap.add_argument("--workers", type=int, default=1)
    ap.add_argument("--json", required=True)
    ap.add_argument("--resume", action="store_true", help="keep the rows already in --json and skip their tasks")
    args = ap.parse_args()
    assert args.workers <= 4
    tasks = []
    only = {int(v) for v in args.only.split(",") if v}
    xqs = {float(v) for v in args.xq.split(",") if v}

    def add(cons, row, bound_upper, dbeta=0):
        D, s = row["D"], row["s"]
        q = 1 if D == 1 else abs(D)
        if row.get("support_2a"):
            x = math.exp(float(row["support_2a"]))
            N = SUPPORT_N[row["support_2a"]]
            sdps = 0
        else:
            x = float(row["xq"]) * q
            n2, sdps, _ = scan_n2(D, row["xq"], row["parity"])
            N = min(n2 or args.nmax, args.nmax)
        c = 2 * math.pi * x / q
        dps = max(int(2 * c / math.log(10)) + 40, 50)
        tasks.append({"cons": cons, "D": D, "s": s, "xq": row.get("xq"), "support": row.get("support_2a"), "N": N,
                      "dps": dps, "n_out": args.n_out, "dbeta": dbeta, "bound_upper": bound_upper})

    keys_kb = set()
    if args.grid:
        for r in json.load(open(args.kb))["rows"]:
            if (only and r["D"] not in only) or (xqs and float(r["xq"]) not in xqs):
                continue
            add("kb", r, r["bound_upper"], r["best_dbeta"])
            keys_kb.add((r["D"], float(r["xq"]), r["s"]))
        for r in json.load(open(args.herm))["rows"]:
            k = (r["D"], float(r["xq"]), r["s"])
            if (only and r["D"] not in only) or (xqs and float(r["xq"]) not in xqs):
                continue
            if (k in keys_kb or args.herm_extra) and not args.kb_only:
                add("hermite", r, r["bound_upper"])
    if args.supports:
        for r in json.load(open(args.kb_supports))["rows"]:
            add("kb", r, r["bound_upper"], r["best_dbeta"])
        for r in json.load(open(args.herm_supports))["rows"]:
            if not args.kb_only:
                add("hermite", r, r["bound_upper"])
    rows = []
    if args.resume and os.path.exists(args.json):
        rows = json.load(open(args.json))["rows"]
        done = {(r["construction"], r["D"], r["parity"], r.get("xq"), r.get("support_2a"), r.get("n_out")) for r in rows}
        tasks = [t for t in tasks if (t["cons"], t["D"], ("even", "odd")[t["s"]], t.get("xq"), t.get("support"), t["n_out"]) not in done]
    job = Job("r2 rq check", total=len(tasks), args=vars(args))
    with mproc.Pool(args.workers) as pool:
        for t, out in zip(tasks, pool.imap(rq_point, tasks)):
            out["bound_upper"] = t["bound_upper"]
            rq = mp.mpf(out["rq"])
            bu = mp.mpf(t["bound_upper"])
            out["rq_le_bound"] = bool(rq <= bu)
            out["bound_over_rq"] = mp.nstr(bu / rq, 6) if rq > 0 else None
            rows.append(out)
            job.step(f"{out['construction']} D={out['D']} {out['xq'] or 'e^'+str(out['support_2a'])} {out['parity']} N={out['N']}: "
                     f"RQ={out['rq']} <= B_cert={t['bound_upper']}: {out['rq_le_bound']} (B/RQ {out['bound_over_rq']}) "
                     f"stable={out['stable']} captured={out['captured']:.12f} ({out['seconds']}s)")
            with open(args.json, "w") as fh:
                json.dump({"args": vars(args), "rows": rows}, fh, indent=1)
    job.done()


if __name__ == "__main__":
    main()
