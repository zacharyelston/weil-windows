"""Issue zacharyelston/rh2#7, step 2: the matched-pair decay-law test.

A = L(s, ρ) for the Rédei field of (5, 29), primitive, degree 2, Maass type (Γ_R(s)², N = 145, ε = +1, m₀ = 0).
B = L(s, χ₅)L(s, χ₂₉), a product with the same N, gamma factor, root number and m₀.
Question (R-c, H_same): does the window minimum depend only on (d, N, μ, ε, m₀), or does it see primitivity?

Subcommands:
  prep   (needs cypari2; run on a Mac) Λ(p^k)/p^{k/2} tables for A and B up to x_max, and zero files from step 1.
  check  zero-sum gate (P-DL7 criteria) with the Maass-type tail density.
  scan   P-DL12 protocol: v = √(x/145) on the grid, both sectors, N_basis = ⌊f·k⌋ + 8 for f = 5, 9 (k = 2vL),
         precision-stable (stable_lams), inverse iteration; mpmath/FLINT only (runs in the Docker image).
  fit    R-a, R-b, R-c as registered in docs/BORROMEAN.md.
"""
import argparse, json, math, os, sys, time
import mpmath as mp
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, ".."))
import gl2_form_mp as g2   # noqa: E402
import decay_law_mp as dl  # noqa: E402

OUT = os.path.join(HERE, "..", "..", "data", "borromean")
N_COND = 145
MUS = [mp.mpf(0), mp.mpf(0)]
OBJS = {"rho": "L(s, ρ_{5,29})", "chi5chi29": "L(s, χ₅)L(s, χ₂₉)"}


# ---------------- prep (PARI) ----------------

def cmd_prep(args):
    sys.path.insert(0, HERE)
    import redei_lib as R
    import step1_gates as S1
    R.set_precision(38)
    R.rho_lfun(5, 29, "Lr")
    xmax = int(N_COND * float(args.vmax) ** 2) + 1
    an = R.an_list("Lr", xmax)
    # exact integer coefficients c(p^k) with Λ(p^k) = c(p^k)·log p (the logs are formed at working precision later;
    # storing Λ as float64 put a 1e-17 noise floor under λ ≈ e^{-150} in the first scan, which is void)
    sieve = bytearray([1]) * (xmax + 1)
    sieve[0:2] = b"\x00\x00"
    for i in range(2, int(xmax ** 0.5) + 1):
        if sieve[i]:
            sieve[i * i::i] = bytearray(len(sieve[i * i::i]))
    coef_rho, coef_B = {}, {}
    for p in range(2, xmax + 1):
        if not sieve[p]:
            continue
        c5, c29 = int(R.pari.kronecker(5, p)), int(R.pari.kronecker(29, p))
        pw, k, cs = p, 1, []
        while pw <= xmax:
            ck = k * an[pw] - sum(cs[j - 1] * an[p ** (k - j)] for j in range(1, k))
            cs.append(ck)
            coef_rho[pw] = [p, ck]
            coef_B[pw] = [p, c5 ** k + c29 ** k]
            pw *= p
            k += 1
    out = {"N": N_COND, "mus": [0, 0], "eps": 1, "m0": 0, "xmax": xmax,
           "coef": {"rho": {str(n): v for n, v in coef_rho.items()}, "chi5chi29": {str(n): v for n, v in coef_B.items()}},
           "note": "Λ(p^k) = coef·log p with integer coef; a_p of ρ in coef[p][1]",
           "a_n_check": {"rho_a_p_first": an[2:60]}}
    with open(os.path.join(OUT, "step2_objects.json"), "w") as fh:
        json.dump(out, fh)
    z = json.load(open(os.path.join(OUT, "step1_zeros.json")))
    s1 = json.load(open(os.path.join(OUT, "step1.json")))
    for obj, key, ok in (("rho", "rho_38", s1["rho"]["pass"]),
                         ("chi5chi29", "chi5_chi29_union", s1["chars"]["5"]["g1"]["pass"] and s1["chars"]["29"]["g1"]["pass"])):
        with open(os.path.join(OUT, f"zeros_step2_{obj}.json"), "w") as fh:
            json.dump({"obj": obj, "tmax": z["T"], "zeros": z[key], "central_multiplicity": 0, "g1_pass": ok,
                       "source": "data/borromean/step1_zeros.json"}, fh)
    print("wrote step2_objects.json (xmax", xmax, ") and zeros_step2_{rho,chi5chi29}.json")


# ---------------- forms ----------------

_TABLES = None


def tables():
    global _TABLES
    if _TABLES is None:
        _TABLES = json.load(open(os.path.join(OUT, "step2_objects.json")))
    return _TABLES


def terms(obj, x):
    """[(log p^k, c·log p / p^{k/2})] for p^k ≤ x, built at the working precision from integer coefficients."""
    co = tables()["coef"][obj]
    xi = int(mp.floor(x))
    out = []
    for n_s, (p, c) in sorted(co.items(), key=lambda kv: int(kv[0])):
        n = int(n_s)
        if n > xi:
            break
        assert isinstance(c, int) and isinstance(p, int)
        if c:
            out.append((mp.log(n), c * mp.log(p) / mp.sqrt(n)))
    return out


def form(obj, x, n, parity):
    if int(mp.floor(x)) > tables()["xmax"]:
        raise ValueError("x beyond the prepared Λ table")
    return g2.form(x, n, parity, MUS, log_cond=mp.log(N_COND), terms=terms(obj, x))


def density(t):
    """dN/dt = (1/2π)[log N + Σ_j (−log π + Re ψ((½ + μ_j + it)/2))]."""
    return (mp.log(N_COND) + sum(-mp.log(mp.pi) + mp.re(mp.digamma((mp.mpf(1) / 2 + mu + 1j * t) / 2)) for mu in MUS)) / (2 * mp.pi)


def tail_estimate(Ffun, T, L):
    with mp.workdps(20):
        pts = mp.linspace(mp.mpf(T), 40 * mp.mpf(T), int(39 * T * L / (2 * mp.pi)) + 2)
        return 2 * mp.quad(lambda t: Ffun(t) ** 2 * density(t), pts)


# ---------------- check ----------------

def cmd_check(args):
    """Q(c) = cᵀMc against 2 Σ_{0<γ≤T} F_c(γ)² + tail for edge-vanishing test functions (as gl2_form_mp check).
    Gate (P-DL7): at x = 13, |Q − Σ| ≤ 2e-4 and (Q − Σ)/tail ∈ [0.8, 1.25] for every test; x = 60 reported."""
    out = []
    for obj in OBJS:
        zd = json.load(open(os.path.join(OUT, f"zeros_step2_{obj}.json")))
        assert zd["g1_pass"], f"step 1 G1 did not pass for {obj}"
        for xs in args.x.split(","):
            mp.mp.dps = 30
            x = mp.mpf(xs)
            L = mp.log(x)
            gam = [mp.mpf(g) for g in zd["zeros"]]
            T = mp.mpf(zd["tmax"])
            Mo, Me = form(obj, x, 6, "odd"), form(obj, x, 6, "even")
            cases = []
            for k in (2, 4):
                c = [0] * 6
                c[k - 1] = 1
                cases.append(("odd", c, Mo))
            cases.append(("odd", [1 if i == 1 else -1 if i == 2 else 0 for i in range(6)], Mo))
            cases += [("even", [0, 1, 1, 0, 0, 0, 0], Me), ("even", [0, 0, 1, 0, -1, 0, 0], Me)]
            for par, c, M in cases:
                if par == "odd":
                    Ff = lambda t, c=c: mp.fsum(ci * dl.F_per_odd(i + 1, L, t) for i, ci in enumerate(c) if ci)
                else:
                    Ff = lambda t, c=c: dl.F_per_even(c, L, t)
                cv = mp.matrix(c)
                q = (cv.T * M * cv)[0]
                fz = mp.fsum(2 * Ff(g) ** 2 for g in gam)
                tail = tail_estimate(Ff, T, L)
                ratio = (q - fz) / tail if tail != 0 else mp.inf
                gate = xs == "13"
                ok = (abs(q - fz) <= mp.mpf("2e-4") and mp.mpf("0.8") <= ratio <= mp.mpf("1.25")) if gate else None
                out.append({"obj": obj, "x": xs, "parity": par, "c": c, "Q": mp.nstr(q, 15), "zero_sum": mp.nstr(fz, 15),
                            "abs_diff": mp.nstr(q - fz, 4), "tail": mp.nstr(tail, 4), "ratio": mp.nstr(ratio, 5),
                            "gating": gate, "pass": ok, "n_zeros": len(gam), "T": mp.nstr(T, 6)})
                print(f"{obj:<10} x={xs:<3} {par:<4} c={c}: Q−Σ = {mp.nstr(q - fz, 4):>10}  tail {mp.nstr(tail, 4):>10}  "
                      f"ratio {mp.nstr(ratio, 4):>7}  {'PASS' if ok else 'FAIL' if ok is False else '(not gating)'}", flush=True)
    verdict = all(r["pass"] for r in out if r["gating"])
    print("check gate:", "PASS" if verdict else "FAIL")
    with open(args.json or os.path.join(OUT, "step2_check.json"), "w") as fh:
        json.dump({"rows": out, "pass": verdict}, fh, indent=1)


# ---------------- scan ----------------

def cmd_scan(args):
    obj = args.obj
    vs = [mp.mpf(s) for s in args.v.split(",")]
    path = args.json or os.path.join(OUT, f"scan_step2_{obj}.json")
    rows = json.load(open(path)) if (args.resume and os.path.exists(path)) else []
    done = {(r["v"], r["parity"]) for r in rows}
    parities = args.parity.split(",")
    for v in vs:
        x = v * v * N_COND
        L = mp.log(x)
        kk = float(2 * v * L)
        for parity in parities:
            if (mp.nstr(v, 10), parity) in done:
                continue
            ns = [max(args.nmin, int(f * kk) + 8) for f in (args.f1, args.f2)]
            ns[1] = max(ns[1], ns[0] + 8)
            res = {}
            for n in ns:
                t0 = time.time()

                def build(n=n):
                    return {"L": form(obj, x, n, parity)}

                dps0 = int(40 + 1.5 * 8 * float(mp.pi) * float(v) / 2.303)
                del g2.FALLBACKS[:]
                lam, at = dl.stable_lams(build, dps0, method=g2.lam_min_inv_checked)
                res[n] = {"n": n, "lambda": lam["L"], "dps": at["L"], "secs": time.time() - t0}
                print(f"{obj} v={mp.nstr(v, 4)} {parity} N_basis={n}: λ = {mp.nstr(lam['L'], 12)} at dps {at['L']} ({time.time() - t0:.0f}s)", flush=True)
            r1, r2 = res[ns[0]], res[ns[1]]
            conv = mp.log(r1["lambda"] / r2["lambda"])
            rows.append({"obj": obj, "N": N_COND, "v": mp.nstr(v, 10), "x": mp.nstr(x, 15), "parity": parity, "k": kk,
                         "n1": ns[0], "lambda_n1": mp.nstr(r1["lambda"], 15), "dps_n1": r1["dps"],
                         "n2": ns[1], "lambda_n2": mp.nstr(r2["lambda"], 15), "dps_n2": r2["dps"],
                         "ln_l1_over_l2": mp.nstr(conv, 6), "secs": round(r1["secs"] + r2["secs"], 1)})
            with open(path, "w") as fh:
                json.dump(rows, fh, indent=1)


# ---------------- fit ----------------

def cmd_fit(args):
    """R-a (rate of A), R-b (sector order, both objects), R-c (H_same) as registered in docs/BORROMEAN.md."""
    mp.mp.dps = 30
    four_pi = 4 * mp.pi
    data = {o: json.load(open(os.path.join(args.dir, f"scan_step2_{o}.json"))) for o in OBJS}
    out = {"R-a": [], "R-b": [], "R-c_index": [], "R-c_points": []}
    nstar = {}
    for o, rows in data.items():
        for par in ("even", "odd"):
            Rr = sorted((r for r in rows if r["parity"] == par), key=lambda r: float(r["v"]))
            v = [mp.mpf(r["v"]) for r in Rr]
            y = [mp.log(mp.mpf(r["lambda_n2"])) for r in Rr]
            c, res = dl._lsq([[-vi, mp.log(vi), 1] for vi in v], y)
            cl_, _ = dl._lsq([[-vi, 1] for vi in v], y)
            a4 = c[0] / four_pi
            ns, _, chk = g2.continuous_index(v, [mp.mpf(r["lambda_n2"]) for r in Rr], extrapolate=True)
            nstar[(o, par)] = ns
            out["R-a"].append({"obj": o, "parity": par, "alpha_over_4pi": mp.nstr(a4, 6), "gamma": mp.nstr(c[1], 6),
                               "max_resid": mp.nstr(res, 4), "linear_alpha_over_4pi": mp.nstr(cl_[0] / four_pi, 6),
                               "n_star": mp.nstr(ns, 6), "S_minus_half": mp.nstr(chk, 6),
                               "conv": [r["ln_l1_over_l2"] for r in Rr]})
        ev = {r["v"]: mp.mpf(r["lambda_n2"]) for r in rows if r["parity"] == "even"}
        od = {r["v"]: mp.mpf(r["lambda_n2"]) for r in rows if r["parity"] == "odd"}
        for vv in sorted(ev, key=float):
            lr = mp.log(od[vv] / ev[vv])
            out["R-b"].append({"obj": o, "v": vv, "ln_odd_over_even": mp.nstr(lr, 6), "odd_above_even": bool(lr > 0),
                               "both_below_1e-6": bool(ev[vv] < mp.mpf("1e-6") and od[vv] < mp.mpf("1e-6"))})
    for par in ("even", "odd"):
        d = abs(nstar[("rho", par)] - nstar[("chi5chi29", par)])
        out["R-c_index"].append({"parity": par, "n_star_A": mp.nstr(nstar[("rho", par)], 6),
                                 "n_star_B": mp.nstr(nstar[("chi5chi29", par)], 6), "abs_diff": mp.nstr(d, 5)})
        A = {r["v"]: mp.mpf(r["lambda_n2"]) for r in data["rho"] if r["parity"] == par}
        B = {r["v"]: mp.mpf(r["lambda_n2"]) for r in data["chi5chi29"] if r["parity"] == par}
        for vv in sorted(set(A) & set(B), key=float):
            out["R-c_points"].append({"parity": par, "v": vv, "dln": mp.nstr(mp.log(A[vv] / B[vv]), 6)})
    print(json.dumps(out, indent=1))
    if args.json:
        with open(args.json, "w") as fh:
            json.dump(out, fh, indent=1)


def cmd_validate(args):
    """Gates on the production path (run before any scan):
    V1: ρ's integer coefficients satisfy the Euler-factor recursion c_k = a_p c_{k−1} − χ₁₄₅(p) c_{k−2} (c_0 = 2) at
        unramified p, and c_k = a_p^k at p = 5, 29, for every p^k ≤ xmax.
    V2: B's production-path form equals dl.zeros_side(5) + dl.zeros_side(29) entrywise to ≤ 1e-100 at 120 digits,
        both sectors, x = 13, 1000, 14500 (n = 12). Positive control: the same comparison with a float64 copy of the
        prime terms must FAIL (difference ≥ 1e-30)."""
    co = tables()["coef"]["rho"]
    bad = []
    by_p = {}
    for n_s, (p, c) in co.items():
        by_p.setdefault(p, {})[int(n_s)] = c
    for p, d in by_p.items():
        ks = sorted(d)
        ap = d[p]
        chi = dl.kronecker(145, p)
        prev2, prev1 = 2, ap
        for k_i, n in enumerate(ks, start=1):
            if p in (5, 29):
                exp = ap ** k_i
            elif k_i == 1:
                exp = ap
            else:
                exp = ap * prev1 - chi * prev2
                prev2, prev1 = prev1, exp
            if d[n] != exp:
                bad.append((n, d[n], exp))
    v1 = not bad
    print("V1 coefficient recursion:", "PASS" if v1 else f"FAIL {bad[:5]}", f"({len(co)} prime powers)", flush=True)
    rows, v2 = [], True
    mp.mp.dps = 120
    for x in (13, 1000, 14500):
        for parity in ("even", "odd"):
            A = form("chi5chi29", mp.mpf(x), 12, parity)
            B = dl.zeros_side(5, x, 12, parity) + dl.zeros_side(29, x, 12, parity)
            d = max(abs(A[i, j] - B[i, j]) for i in range(A.rows) for j in range(A.cols))
            ok = d <= mp.mpf("1e-100")
            v2 &= ok
            rows.append({"x": x, "parity": parity, "max_diff": mp.nstr(d, 3), "pass": bool(ok)})
            print(f"V2 x={x} {parity}: max |production − (χ₅ + χ₂₉)| = {mp.nstr(d, 3)}  {'PASS' if ok else 'FAIL'}", flush=True)
    # positive control: float64 prime terms must be caught
    tf = [(lg, mp.mpf(float(v))) for lg, v in terms("chi5chi29", mp.mpf(14500))]
    Af = g2.form(mp.mpf(14500), 12, "even", MUS, log_cond=mp.log(N_COND), terms=tf)
    B = dl.zeros_side(5, 14500, 12, "even") + dl.zeros_side(29, 14500, 12, "even")
    dctl = max(abs(Af[i, j] - B[i, j]) for i in range(Af.rows) for j in range(Af.cols))
    ctl = dctl >= mp.mpf("1e-30")
    print(f"V2 positive control (float64 terms): max diff = {mp.nstr(dctl, 3)}  {'FIRED' if ctl else 'DID NOT FIRE'}", flush=True)
    out = {"V1": v1, "V1_bad": bad[:20], "V2": v2, "V2_rows": rows, "positive_control_diff": mp.nstr(dctl, 3),
           "positive_control_fired": bool(ctl), "pass": bool(v1 and v2 and ctl)}
    with open(os.path.join(OUT, "step2_validate.json"), "w") as fh:
        json.dump(out, fh, indent=1)
    print("validate:", "PASS" if out["pass"] else "FAIL")


def cmd_merge(args):
    """Merge per-point scan files data/borromean/parts/scan_step2_<obj>_<v>_<parity>.json into scan_step2_<obj>.json."""
    import glob
    for o in OBJS:
        rows = []
        for f in sorted(glob.glob(os.path.join(args.dir, "parts", f"scan_step2_{o}_*.json"))):
            rows += json.load(open(f))
        rows.sort(key=lambda r: (float(r["v"]), r["parity"]))
        with open(os.path.join(args.dir, f"scan_step2_{o}.json"), "w") as fh:
            json.dump(rows, fh, indent=1)
        print(o, len(rows), "rows")


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("prep").add_argument("--vmax", default="10")
    c = sub.add_parser("check"); c.add_argument("--x", default="13,60"); c.add_argument("--json")
    s = sub.add_parser("scan")
    s.add_argument("--obj", required=True, choices=list(OBJS))
    s.add_argument("--v", default="6,7,8,9,10")
    s.add_argument("--parity", default="even,odd")
    s.add_argument("--f1", type=float, default=5); s.add_argument("--f2", type=float, default=9)
    s.add_argument("--nmin", type=int, default=24)
    s.add_argument("--resume", action="store_true"); s.add_argument("--json")
    f = sub.add_parser("fit"); f.add_argument("--dir", default=OUT); f.add_argument("--json")
    mg = sub.add_parser("merge"); mg.add_argument("--dir", default=OUT)
    sub.add_parser("validate")
    a = ap.parse_args()
    {"prep": cmd_prep, "check": cmd_check, "scan": cmd_scan, "fit": cmd_fit, "merge": cmd_merge, "validate": cmd_validate}[a.cmd](a)


if __name__ == "__main__":
    main()
