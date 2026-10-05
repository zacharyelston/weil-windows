"""Issue the private research log (issue 7), step 1: pipeline gates for the Rédei L-function L(s, ρ) of (p1, p2) = (5, 29).

Pre-registered in docs/BORROMEAN.md (the review's step 1). Gates only; there is no scientific kill in this step.
Stages, in the registered order:
  g0      inputs: field, D4, conductor, gamma factor, root number, feq, trace rule p ≤ 10^4, reciprocity, (13, 61, 937)
  proxy   L(χ12)L(χ13) (N = 156, Γ_R(s)²): positive control for G1 (precision 19, one lfunzeros call, default divz)
          must FAIL G1; good settings must pass G1; its G2 residual fixes τ = 3 × max residual
  chars   L(χ5), L(χ29) to T, chunks of 100, divz 32, 38 digits, with G1 (for G3's negative control)
  rho     L(ρ) to T, chunks of 100, divz 32, at 38 and 57 digits; G1 at both; zero lists must agree
  g2      |S(x) − P(x)| ≤ τ at 2000 log-spaced x in [2, 1000]
  g3      |R(p) − a_p| ≤ 0.3 for every unlinked p ≤ 1000; negative control on L(χ5)L(χ29)'s zeros
Definitions:
  S(x) = Σ_{0<γ≤T} (1 − γ/T) cos(γ log x);  R(x) = −S(x) / [(T/4π) log x / √x].
  P(x) = (1/2π) ∫_0^T (1 − t/T) cos(t log x) ρ_dens(t) dt − (T/4π) Σ_{n ≤ M} Λ(n)/√n [F(log n − log x) + F(log n + log x)],
         F(v) = (sin(Tv/2)/(Tv/2))², ρ_dens(t) = log N + Σ_j [−log π + Re ψ((½ + μ_j + it)/2)], M = 4·10^5
         (Weil's explicit formula for h(t) = (1 − |t|/T)₊ cos(t log x); no pole; central order m₀ = 0, checked).
  G1: (a) argument principle at T = 100 and 200: N(T) = [θ(T) + arg L(½ + iT)]/π − m₀/2, arg L continuous from
      σ = 2 (principal there: |Im log L| ≤ 2 log ζ(2) < π), θ(T) = (T/2) log N + Σ_j [−(T/2) log π + Im log Γ((½+μ_j+iT)/2)];
      the found counts must equal N(100) and N(200), each within 1e-6 of an integer.
      (b) Turing's method at every boundary b = 100, 200, ..., 1000: D(b) = mean over t ∈ [b − 20, b] (step 0.02) of
      θ(t)/π − N_found(t) must satisfy |D(b)| < 0.5 (a missed zero below b shifts D by about 1).
      Off-line evaluation of L at height ~1000 takes minutes per point in PARI, so (a) cannot reach T = 1000.
"""
import argparse, json, math, random, sys, time, os
import mpmath as mp
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import redei_lib as R
pari = R.pari

P1, P2 = 5, 29
T = 1000
CHUNK = 100
M_PRIME = 400000


def log(*a):
    print(*a, flush=True)


# ---------------- argument principle ----------------

def theta(Tv, N, vga):
    Tv = mp.mpf(Tv)
    th = Tv / 2 * mp.log(N)
    for mu in vga:
        th += -Tv / 2 * mp.log(mp.pi) + mp.im(mp.loggamma((mp.mpf(1) / 2 + mu + 1j * Tv) / 2))
    return th


def arg_L_line(name, Tv, steps=160):
    """arg L(½ + iT), continuous along σ from 2 down to ½ at height T (principal at σ = 2)."""
    def a(sig):
        return float(pari(f"arg(lfun({name}, {sig} + I*{Tv}))"))
    sig = [2 - 1.5 * k / steps for k in range(steps + 1)]
    vals = [a(s) for s in sig]
    total = vals[0]
    for k in range(steps):
        lo, hi, alo, ahi = sig[k], sig[k + 1], vals[k], vals[k + 1]
        d = ahi - alo
        depth = 0
        if abs((d + math.pi) % (2 * math.pi) - math.pi) > math.pi / 3:
            # refine this step
            sub = 64
            ss = [lo + (hi - lo) * j / sub for j in range(sub + 1)]
            av = [a(s) for s in ss]
            for j in range(sub):
                dd = av[j + 1] - av[j]
                total += (dd + math.pi) % (2 * math.pi) - math.pi
            continue
        total += (d + math.pi) % (2 * math.pi) - math.pi
    absL = float(pari(f"abs(lfun({name}, 1/2 + I*{Tv}))"))
    return total, absL


def N_count(name, Tv, N, vga, m0=0):
    th = theta(Tv, N, vga)
    ar, absL = arg_L_line(name, Tv)
    val = float((th + ar) / mp.pi) - m0 / 2
    return val, absL


def central_order(name):
    v = abs(complex(pari(f"lfun({name}, 1/2)")))
    return (0 if v > 1e-8 else int(pari(f"lfunorderzero({name})"))), v


def counts_by_chunk(zs, Tmax=T, chunk=CHUNK):
    out = [0] * (Tmax // chunk)
    for g in zs:
        if 0 < g <= Tmax:
            out[min(int(math.ceil(g / chunk)) - 1, len(out) - 1)] += 1
    return out


def ap_counts(name, N, vga, m0, digits=38, heights=(100, 200)):
    """Argument-principle N(T) at the given heights (only where off-line evaluation is cheap)."""
    R.set_precision(digits)
    out = []
    for Tv in heights:
        v, absL = N_count(name, Tv, N, vga, m0)
        out.append({"T": Tv, "N": v, "N_int": round(v), "dist_int": abs(v - round(v)), "absL_line": absL})
    return out


def turing_D(zs, b, N, vga, h=20.0, step=0.02):
    """Mean over t in [b − h, b] of θ(t)/π − N_found(t) (Turing's method; ≈ number of missed zeros below b)."""
    import bisect
    ts = [b - h + step * (k + 0.5) for k in range(int(round(h / step)))]
    return sum(float(theta(t, N, vga) / mp.pi) - bisect.bisect_right(zs, t) for t in ts) / len(ts)


def zeros_chunked(name, divz=32, N=None, vga=None, repair=True, log_out=None):
    """Zeros in (0, T] by chunks of CHUNK, all chunks sharing one critical-line lfuninit up to T + 50 (PARI loses zeros
    near the top of an lfuninit domain: the proxy dropped about 20 in [900, 1000] with init to T + 1).
    Repair (registered): a chunk (a, b] whose Turing value jumps, D(b) − D(a) > 0.5, is searched again with divz 128 and
    then 512, and the union is kept. Sign-change searches can step over a close pair of zeros (the proxy lost 2 in
    (500, 600] at divz 32). Every repair is logged; the positive control is never repaired."""
    import bisect
    pari(f"{name}_zi = lfuninit({name}, [{T + 50}]);")
    chunks = [R.zeros(f"{name}_zi", CHUNK * k, CHUNK * (k + 1), divz) for k in range(T // CHUNK)]

    def merged():
        out = []
        for z in sorted(g for c in chunks for g in c if 0 < g <= T):
            if not out or z - out[-1] > 1e-9:
                out.append(z)
        return out

    zs = merged()
    if repair and N is not None:
        for k in range(T // CHUNK):
            a_, b_ = CHUNK * k, CHUNK * (k + 1)
            for dz in (128, 512):
                Da = turing_D(zs, a_, N, vga) if a_ >= 20 else 0.0
                Db = turing_D(zs, b_, N, vga)
                if Db - Da <= 0.5:
                    break
                before = len(chunks[k])
                chunks[k] = sorted(set(chunks[k]) | set(R.zeros(f"{name}_zi", a_, b_, dz)))
                zs = merged()
                if log_out is not None:
                    log_out.append({"chunk": [a_, b_], "divz": dz, "D_jump_before": Db - Da, "added": len(chunks[k]) - before})
    return zs


def g1(ap, zs, N, vga):
    """G1 = (a) argument-principle counts at 100, 200 equal the found counts, integral to 1e-6; (b) |D(b)| < 0.5."""
    import bisect
    a_rows = [{**r, "found": bisect.bisect_right(zs, r["T"])} for r in ap]
    a_ok = all(r["dist_int"] < 1e-6 and r["N_int"] == r["found"] for r in a_rows)
    D = [turing_D(zs, b, N, vga) for b in range(CHUNK, T + 1, CHUNK)]
    b_ok = all(abs(d) < 0.5 for d in D)
    return {"argument_principle": a_rows, "turing_D": D, "max_abs_D": max(abs(d) for d in D),
            "chunks_found": counts_by_chunk(zs), "n_zeros": len(zs), "pass_a": a_ok, "pass_b": b_ok, "pass": a_ok and b_ok}


# ---------------- explicit formula ----------------

def lambda_prime_powers(an, M):
    """{p^k: Λ(p^k)} for p^k ≤ M from Dirichlet coefficients: Λ(p^k) = k log p a_{p^k} − Σ_{j<k} Λ(p^j) a_{p^{k−j}}."""
    lam = {}
    sieve = bytearray([1]) * (M + 1)
    sieve[0:2] = b"\x00\x00"
    for i in range(2, int(M ** 0.5) + 1):
        if sieve[i]:
            sieve[i * i::i] = bytearray(len(sieve[i * i::i]))
    for p in range(2, M + 1):
        if not sieve[p]:
            continue
        lp = math.log(p)
        pw, k, vals = p, 1, []
        while pw <= M:
            v = k * lp * an[pw] - sum(vals[j - 1] * an[p ** (k - j)] for j in range(1, k))
            vals.append(v)
            lam[pw] = v
            pw *= p
            k += 1
    return lam


def density_grid(N, vga, dt=0.05):
    ts = [i * dt for i in range(int(T / dt) + 1)]
    lg = [math.log(N) + sum(-math.log(math.pi) + float(mp.re(mp.digamma((mp.mpf(1) / 2 + mu + 1j * t) / 2))) for mu in vga)
          for t in ts]
    G = [(1 - t / T) * g for t, g in zip(ts, lg)]
    return ts, G


def smooth_term(x, ts, G):
    """(1/2π) ∫_0^T G(t) cos(ωt) dt, ω = log x, G piecewise linear on the grid (exact for the interpolant)."""
    w = math.log(x)
    tot = 0.0
    sn = [math.sin(w * t) for t in ts]
    cs = [math.cos(w * t) for t in ts]
    for i in range(len(ts) - 1):
        a, b = ts[i], ts[i + 1]
        s = (G[i + 1] - G[i]) / (b - a)
        # ∫_a^b (G_i + s(t−a)) cos(wt) dt = [(G_i + s(t−a)) sin(wt)/w + s cos(wt)/w²]_a^b
        tot += (G[i + 1] * sn[i + 1] - G[i] * sn[i]) / w + s * (cs[i + 1] - cs[i]) / (w * w)
    return tot / (2 * math.pi)


def prime_term(x, lam_items):
    lx = math.log(x)
    tot = 0.0
    for ln_, c in lam_items:
        v1, v2 = ln_ - lx, ln_ + lx
        z1, z2 = T * v1 / 2, T * v2 / 2
        F1 = 1.0 if abs(z1) < 1e-12 else (math.sin(z1) / z1) ** 2
        F2 = (math.sin(z2) / z2) ** 2
        tot += c * (F1 + F2)
    return -T / (4 * math.pi) * tot


def S_zeros(x, zs):
    lx = math.log(x)
    return sum((1 - g / T) * math.cos(g * lx) for g in zs if 0 < g <= T)


def g2_residuals(zs, an, N, vga, xs):
    lam = lambda_prime_powers(an, M_PRIME)
    items = [(math.log(n), v / math.sqrt(n)) for n, v in lam.items() if v != 0]
    ts, G = density_grid(N, vga)
    res = []
    for x in xs:
        res.append(S_zeros(x, zs) - (smooth_term(x, ts, G) + prime_term(x, items)))
    return res


def R_read(x, zs):
    return -S_zeros(x, zs) / ((T / (4 * math.pi)) * math.log(x) / math.sqrt(x))


# ---------------- stages ----------------

def stage_g0(out):
    t0 = time.time()
    R.set_precision(38)
    S, info = R.redei_field(P1, P2)
    G = pari.galoisinit(S)
    gid = [int(v) for v in pari.galoisidentify(G)]
    sig = [int(v) for v in pari(f"nfinit({S}).sign")]
    disc = int(pari.nfdisc(S))
    R.rho_lfun(P1, P2, "Lr")
    N, vga, w = R.params("Lr")
    feq = float(pari("lfuncheckfeq(Lr)"))
    an = R.an_list("Lr", 10000)
    primes = [int(p) for p in pari("primes([2, 10000])")]
    bad = []
    for p in primes:
        if p in (P1, P2):
            continue
        rule = 2 * R.redei_symbol(P1, P2, p)
        if an[p] != rule:
            bad.append((p, an[p], rule))
    # reciprocity on 25 random admissible triples below 400 (seed 7)
    rng = random.Random(7)
    pool = [int(p) for p in pari("primes([3, 400])") if int(p) % 4 == 1]
    triples, tries = [], 0
    while len(triples) < 25 and tries < 100000:
        tries += 1
        a, b, c = sorted(rng.sample(pool, 3))
        if R.admissible(a, b) and R.admissible(b, c) and R.admissible(a, c) and (a, b, c) not in triples:
            triples.append((a, b, c))
    recip = []
    for a, b, c in triples:
        s1, s2, s3 = R.redei_symbol(a, b, c), R.redei_symbol(b, c, a), R.redei_symbol(c, a, b)
        recip.append({"triple": [a, b, c], "symbols": [s1, s2, s3], "equal": s1 == s2 == s3})
    pub = [R.redei_symbol(13, 61, 937), R.redei_symbol(61, 937, 13), R.redei_symbol(937, 13, 61)]
    m0, absL_half = central_order("Lr")
    out["g0"] = {
        "field": str(S), "construction": list(info), "galoisidentify": gid, "signature": sig, "disc": disc,
        "conductor": N, "vga": vga, "root_number": [w.real, w.imag], "lfuncheckfeq_bits": feq,
        "a5": an[5], "a29": an[29], "trace_rule_mismatches": bad, "trace_rule_primes": len(primes) - 2,
        "reciprocity": recip, "published_13_61_937": pub, "central_order": m0, "abs_L_half": absL_half,
    }
    ok = (gid == [8, 3] and abs(disc) == 145 ** 4 and sig == [8, 0] and N == 145 and vga == [0, 0]
          and abs(w - 1) < 1e-20 and feq < -100 and not bad and all(r["equal"] for r in recip) and len(recip) == 25
          and pub == [-1, -1, -1] and m0 == 0)
    out["g0"]["pass"] = ok
    out["g0"]["seconds"] = round(time.time() - t0, 1)
    log("G0", "PASS" if ok else "FAIL", {k: out["g0"][k] for k in ("galoisidentify", "conductor", "vga", "lfuncheckfeq_bits", "a5", "a29", "central_order")},
        "mismatches", len(bad), "reciprocity", sum(r["equal"] for r in recip), "/", len(recip), "published", pub)
    return ok


def stage_proxy(out, xs):
    t0 = time.time()
    res = {}
    # argument-principle counts (38 digits)
    R.set_precision(38)
    R.lfun_dirichlet(12, "L12"); R.lfun_dirichlet(13, "L13"); R.lfun_product("L12", "L13", "Lp")
    N, vga, w = R.params("Lp")
    m0, _ = central_order("Lp")
    ap = ap_counts("Lp", N, vga, m0)
    res["params"] = {"conductor": N, "vga": vga, "root_number": [w.real, w.imag], "central_order": m0}
    # positive control: precision 19, one call to T, default divz
    R.set_precision(19)
    R.lfun_dirichlet(12, "L12"); R.lfun_dirichlet(13, "L13"); R.lfun_product("L12", "L13", "Lp")
    zs_default = [float(v) for v in pari(f"lfunzeros(Lp, {T})")]
    res["positive_control"] = g1(ap, sorted(zs_default), N, vga)
    res["positive_control"]["must_fail"] = True
    res["positive_control"]["fired"] = not res["positive_control"]["pass"]
    # good settings
    R.set_precision(38)
    R.lfun_dirichlet(12, "L12"); R.lfun_dirichlet(13, "L13"); R.lfun_product("L12", "L13", "Lp")
    rep = []
    zs = zeros_chunked("Lp", N=N, vga=vga, log_out=rep)
    res["good"] = g1(ap, zs, N, vga)
    res["good"]["repairs"] = rep
    an = R.an_list("Lp", M_PRIME)
    resid = g2_residuals(zs, an, N, vga, xs)
    res["g2_max_residual"] = max(abs(r) for r in resid)
    res["tau"] = 3 * res["g2_max_residual"]
    res["seconds"] = round(time.time() - t0, 1)
    out["proxy"] = res
    log("proxy: positive control fired" if res["positive_control"]["fired"] else "proxy: positive control DID NOT fire",
        "| good G1", res["good"]["pass"], res["good"]["n_zeros"], "max|D|", round(res["good"]["max_abs_D"], 3),
        "| default found", res["positive_control"]["n_zeros"], "max|D|", round(res["positive_control"]["max_abs_D"], 3), "| G2 max residual", res["g2_max_residual"], "tau", res["tau"])
    return res["positive_control"]["fired"] and res["good"]["pass"], zs


def stage_chars(out):
    t0 = time.time()
    res, allz = {}, {}
    for d in (5, 29):
        R.set_precision(38)
        R.lfun_dirichlet(d, f"L{d}")
        N, vga, w = R.params(f"L{d}")
        m0, _ = central_order(f"L{d}")
        ap = ap_counts(f"L{d}", N, vga, m0)
        R.set_precision(38)
        R.lfun_dirichlet(d, f"L{d}")
        rep = []
        zs = zeros_chunked(f"L{d}", N=N, vga=vga, log_out=rep)
        res[str(d)] = {"params": {"conductor": N, "vga": vga, "central_order": m0}, "g1": g1(ap, zs, N, vga), "repairs": rep}
        allz[d] = zs
        log(f"L(chi_{d}): G1", res[str(d)]["g1"]["pass"], res[str(d)]["g1"]["n_zeros"], "max|D|", round(res[str(d)]["g1"]["max_abs_D"], 3))
    res["seconds"] = round(time.time() - t0, 1)
    out["chars"] = res
    return all(res[str(d)]["g1"]["pass"] for d in (5, 29)), sorted(allz[5] + allz[29])


def stage_rho(out):
    t0 = time.time()
    res = {}
    R.set_precision(38)
    R.rho_lfun(P1, P2, "Lr")
    N, vga, w = R.params("Lr")
    m0, _ = central_order("Lr")
    ap = ap_counts("Lr", N, vga, m0)
    zsets = {}
    for digits in (38, 57):
        R.set_precision(digits)
        R.rho_lfun(P1, P2, "Lr")
        rep = []
        zsets[digits] = zeros_chunked("Lr", N=N, vga=vga, log_out=rep)
        res[f"g1_{digits}"] = g1(ap, zsets[digits], N, vga)
        res[f"g1_{digits}"]["repairs"] = rep
        log(f"L(rho) {digits} digits: G1", res[f"g1_{digits}"]["pass"], res[f"g1_{digits}"]["n_zeros"], "max|D|", round(res[f"g1_{digits}"]["max_abs_D"], 3))
    a, b = zsets[38], zsets[57]
    res["lists_agree"] = len(a) == len(b) and max((abs(x - y) for x, y in zip(a, b)), default=0) < 1e-15 * T
    res["max_list_diff"] = max((abs(x - y) for x, y in zip(a, b)), default=None) if len(a) == len(b) else None
    res["seconds"] = round(time.time() - t0, 1)
    res["pass"] = res["g1_38"]["pass"] and res["g1_57"]["pass"] and res["lists_agree"]
    res["_z57"] = b
    out["rho"] = res
    return res["pass"], a, N, vga


def main():
    ap_ = argparse.ArgumentParser()
    ap_.add_argument("--json", default="data/borromean/step1.json")
    ap_.add_argument("--stages", default="g0,proxy,chars,rho,g2,g3")
    args = ap_.parse_args()
    stages = args.stages.split(",")
    xs = [2 * 500 ** (j / 1999) for j in range(2000)]
    out = {"issue": "the private research log (issue 7)", "step": 1, "pair": [P1, P2], "T": T, "chunk": CHUNK, "M_prime": M_PRIME,
           "pari": [int(v) for v in pari("version()")][:3], "stages_run": stages}
    t_all = time.time()
    if "g0" in stages:
        stage_g0(out)
    zs_rho = zs_B = None
    if "proxy" in stages:
        stage_proxy(out, xs)
    if "chars" in stages:
        _, zs_B = stage_chars(out)
    if "rho" in stages:
        _, zs_rho, N, vga = stage_rho(out)
    if "g2" in stages and zs_rho is not None:
        R.set_precision(38); R.rho_lfun(P1, P2, "Lr")
        an = R.an_list("Lr", M_PRIME)
        resid = g2_residuals(zs_rho, an, N, vga, xs)
        tau = out["proxy"]["tau"]
        out["g2"] = {"max_residual": max(abs(r) for r in resid), "tau": tau, "pass": max(abs(r) for r in resid) <= tau,
                     "residual_quantiles": sorted(abs(r) for r in resid)[::400]}
        log("G2", out["g2"])
    if "g3" in stages and zs_rho is not None and zs_B is not None:
        R.set_precision(38); R.rho_lfun(P1, P2, "Lr")
        an = R.an_list("Lr", 1000)
        rows, neg = [], []
        for p in [int(v) for v in pari("primes([2, 1000])")]:
            if p in (P1, P2):
                continue
            r = R_read(p, zs_rho)
            rows.append({"p": p, "a_p": an[p], "R": r, "dev": abs(r - an[p]), "borromean": an[p] == -2 and p % 4 == 1})
            if an[p] == -2:
                rB = R_read(p, zs_B)
                neg.append({"p": p, "R_B": rB, "agrees_with_a_p": abs(rB - an[p]) <= 0.3})
        unl = [r for r in rows if r["a_p"] != 0]
        out["g3"] = {"rows": rows, "negative_control": neg,
                     "max_dev_unlinked": max(r["dev"] for r in unl), "n_unlinked": len(unl),
                     "pass_readout": all(r["dev"] <= 0.3 for r in unl),
                     "neg_agreements": sum(n["agrees_with_a_p"] for n in neg), "n_neg": len(neg),
                     "pass_negative_control": sum(n["agrees_with_a_p"] for n in neg) <= 2}
        log("G3 readout", out["g3"]["pass_readout"], "max dev", out["g3"]["max_dev_unlinked"], "| negative control agreements",
            out["g3"]["neg_agreements"], "/", out["g3"]["n_neg"])
    # zero lists for later steps (output only; no gate uses this file)
    zfile = os.path.splitext(args.json)[0] + "_zeros.json"
    with open(zfile, "w") as fh:
        json.dump({"T": T, "note": "zeros 0 < gamma <= T on the critical line, after the registered repair rule",
                   "rho_38": zs_rho, "chi5_chi29_union": zs_B, "rho_57": out.get("rho", {}).get("_z57")}, fh)
    out.get("rho", {}).pop("_z57", None)
    out["seconds_total"] = round(time.time() - t_all, 1)
    os.makedirs(os.path.dirname(args.json) or ".", exist_ok=True)
    with open(args.json, "w") as fh:
        json.dump(out, fh, indent=1)
    log("wrote", args.json)


if __name__ == "__main__":
    main()
