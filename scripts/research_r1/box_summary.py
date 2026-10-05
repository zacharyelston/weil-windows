"""Summaries for the document: (1) per support, per band: min/median/max of the box lower bound R4_lower (δ = 0 and ½)
and the largest sampled height with R4_lower > 1; (2) at a = 0.8, the matrix R4 (δ = 0, N = 100/160) against the box
value at the same heights (the shaping gain); (3) the envelope height 2π exp(4a cosh²(δa) + A_box + ln gain)."""
import json, glob, math, statistics
S = "data/research_r1"

print("== box bands (R4_lower; δ = 0 | δ = ½): min / median / max per band; max height with R4_lower > 1")
for a_s, f in (("0.8", "box_08"), ("1.19", "box_119"), ("1.3", "box_13"), ("1.495", "box_1495"), ("1.6", "box_16"), ("1.75", "box_175"), ("2.0", "box_20")):
    try:
        d = json.load(open(f"{S}/{f}.json"))
    except Exception:
        print(a_s, "missing"); continue
    rows = d["rows"]
    a = float(a_s)
    bands = sorted({r["band"] for r in rows}, key=float)
    print(f"a = {a_s}: x = {math.exp(2*a):.3f}, T* = {2*math.pi*math.exp(2*a):.1f}, 2πe^{{4a}} = {2*math.pi*math.exp(4*a):.0f}, A_full = {d['A_full']}, A_box = {d['A_box']}")
    best = {"0.0": None, "0.5": None}
    for b in bands:
        rb = [r for r in rows if r["band"] == b]
        gs = [float(r["gamma0"]) for r in rb]
        line = f"   band 4a{float(b):+.1f} (γ₀ ≈ {min(gs):.0f}–{max(gs):.0f}): "
        for dk in ("0.0", "0.5"):
            v = [float(r["R4_lower"][dk]) for r in rb]
            line += f"δ={dk}: {min(v):.3f}/{statistics.median(v):.3f}/{max(v):.3f}   "
            for r in rb:
                if float(r["R4_lower"][dk]) > 1 and (best[dk] is None or float(r["gamma0"]) > best[dk]):
                    best[dk] = float(r["gamma0"])
        rho = [float(r["rho_eff"]) for r in rb]
        line += f"ρ_eff range {min(rho):+.2f}…{max(rho):+.2f}"
        print(line)
    Ab = float(d["A_box"])
    for dk in ("0.0", "0.5"):
        dd = float(dk)
        # box weight w(δa) = (sinh δa / δa)²: R4_box = 4a·w/(log(γ₀/2π) − ρ) ⇒ R4 > 1 needs log(γ₀/2π) < g·4a·w + ρ ≤ g·4a·w + A_box
        w = (math.sinh(dd * a) / (dd * a)) ** 2 if dd else 1.0
        env = 2 * math.pi * math.exp(4 * a * w + Ab)
        env_g = 2 * math.pi * math.exp(1.2 * 4 * a * w + Ab)
        typ_g = 2 * math.pi * math.exp(1.2 * 4 * a * w + 1)
        print(f"   δ={dk}: largest sampled height with R4_lower > 1: {best[dk]};  w(δa) = {w:.4f};  box envelope 2πe^{{4a·w + A_box}} = {env:.3g};  "
              f"with shaping gain 1.2: {env_g:.3g};  typical (ρ = +1, g = 1.2): {typ_g:.3g}")

print("\n== a = 0.8: matrix R4 (δ = 0) against the box at the same heights")
try:
    Z = json.load(open(f"{S}/zeta_all.json"))
    cmp = json.load(open(f"{S}/box_08_cmp.json"))["rows"]
    for z in Z:
        if z["a"] != "0.8":
            continue
        prof = z["profiles_log10R4"]["0"]
        line = f"   {z['sector']} N={z['n']} (ω_N = {z['omega_N']}, T_excl(δ=0) = {z['T_excl']['0']}, δ=0.5: {z['T_excl']['0.5']}): "
        for r in cmp:
            g = float(r["gamma0"])
            near = min(prof, key=lambda p: abs(p[0] - g))
            if near[1] is None:
                continue
            line += f"γ₀={g:.0f}: matrix {10**near[1]:.3f} vs box {float(r['R4_lower']['0.0']):.3f} (gain {10**near[1]/float(r['R4_lower']['0.0']):.2f}); "
        print(line)
except Exception as e:
    print("zeta_all/box cmp missing:", e)

print("\n== comb summary per (a, sector, N)")
try:
    for z in Z:
        comb = z["comb_R2_at_zeros"]
        vals = [(int(k), float(g), float(v)) for k, g, v in comb]
        Ts = float(z["Tstar"])
        below = [v for k, g, v in vals if g < Ts]
        print(f"   a={z['a']} {z['sector']} N={z['n']}: λ_min={z['lam_min']}, T*={Ts:.1f}; zeros below T*: {len(below)}, min 2κ there = {min(below):.12f}; "
              f"2κ at γ₁+1e-20 = {z['comb_R2_shifted_from_gamma1'][2][1]}, at +1e-5 = {z['comb_R2_shifted_from_gamma1'][0][1]}; "
              f"last zero γ_{vals[-1][0]}={vals[-1][1]:.1f}: 2κ={vals[-1][2]:.4f}; midpoint 2κ min = {min(float(v) for _, v in z['comb_R2_midpoints']):.3g}; "
              + (f"real-pair R2 (even): {z['real_pair_even_R2']}" if "real_pair_even_R2" in z else ""))
except Exception as e:
    print("comb missing:", e)
