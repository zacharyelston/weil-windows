# R2 Theorems: a certified, unconditional upper bound for the window minimum

Date: 2026-10-03. Branch `r2-rigorous` (from `feat/weil-gram-instrument` at 8c5eb56). This sheet restates the upper-bound half of `docs/RESEARCH_R2_SAMPLING.md` (R2, commit 4dd51fa) so that an auditor can check every step and every number. New code: `scripts/research_r2_rigorous/`. New data: `data/research_r2_rigorous/`. No existing script or document was modified.

**Labels.**
- [R] read in the source (at most one short quotation per source);
- [A] abstract only;
- [S] second-hand (another paper's or a repository document's account);
- [D] derived here (proof; "sketch" where marked);
- [N] computed in floating point (multiprecision, not certified);
- [C] certified: an Arb ball enclosure (python-flint 0.9.0) whose relevant endpoint is the stated bound.

**What this sheet is not.** An upper bound on the window minimum λ_s(x) is consistent with RH and proves nothing about zeros on the line. Every bound below holds whether or not RH is true, and is equally true for functions that violate RH. Davenport–Heilbronn satisfies Theorem 1's identity and still goes negative (R2 §4.1).

## 0. Summary

**Scope.** Degree 1 only: ζ and real primitive Dirichlet L-functions. Every GL(2) or degree-d statement is out of scope here, except as the marked remark in §10.

**Three kinds of statement are kept apart.**
1. **Proved, unconditional** [D], with every hypothesis stated:
   - Theorem 1, the E-map identity Q(f_win) = σΣ_ρ F_out(γ_ρ)², over all zeros with multiplicity.
   - Theorem 2, the bound λ_s(x) ≤ 2A²∫N⁺t⁻³/‖f_win‖².
   - Corollary 3 (Hermite) and Theorem 4 (Kaiser–Bessel), as explicit formulas for admissible A, B and ‖f_win‖².

   The only external inputs are:
   - the explicit formula [S];
   - two explicit zero-count theorems, read in source [R];
   - two DLMF Bessel inequalities [R].
2. **Finite, certified** [C]. **Theorem C:** for every (𝓛, x, s) listed in Tables S, Hermite and KB (§7), λ_s(x) ≤ B_cert(x). B_cert is the upper endpoint of an Arb ball, with no floating midpoint anywhere in the chain:
   - 286 Hermite and 84 Kaiser–Bessel points of R2's grids;
   - ζ at the four certified supports, both sectors, both constructions.

   This is a statement about those x only. Arb certification proves nothing about other x and nothing asymptotic.
3. **Rates, stated separately.**
   - **Proposition H** (proved, §5.1). For each object and sector, λ_s(x) ≤ K_H c^{d+1}e^{−c} for **all** x ≥ 2q. K_H is a certified constant (26 values in §5.1).
   - **Conjecture R** (the Kaiser–Bessel rate e^{−2c}, §6.1). It is supported by the certified table: ln B_cert + 2c fits C c^p with residuals ≤ 0.47. It is **not proved**. The missing estimates are named exactly: a uniform lower bound for ‖f_win‖² (M1) and uniform polynomial bounds for A, B (M2).

**Corrections to R2.**
1. **The ζ zero count** was cited from Trudgian 2014, whose proof has a published gap. Bennett–Martin–O'Bryant–Rechnitzer, Hasanalizade–Shen–Wong and Bellotti–Wong all point it out. The constants were used inside their stated range, and Bellotti–Wong show they hold anyway; the chain now uses Bellotti–Wong's Theorem 1.1 directly (§2.1).
2. **The Dirichlet count** (BMOR Theorem 1.1) was used correctly (§2.2).
3. **Monotonicity.** R2 asserted that the Kaiser–Bessel envelope is monotone on an 80-point grid. That step is replaced by an identity (Lemma H) and closed-form tails (Lemma T).
4. **Floating quadratures.** The norm and envelope integrals are replaced by `acb.integral` with rigorous remainders.
5. **Not all of R2's tabulated values were proven bounds.** The certified bounds are 1.006–1.103 × R2's floating values (Hermite) and 0.89–1.26 × (Kaiser–Bessel).
6. **R2's "rate e^{−2c}" for Theorem 4** is a conjecture (Conjecture R), not a theorem. R2's rate e^{−c} for Corollary 3 is now a theorem (Proposition H).

**Validation** (§8):
- (a) the trial's own Rayleigh quotient in our basis is ≤ B_cert: **386/386** (84 + 8 Kaiser–Bessel, 286 + 8 Hermite; all precision-stable);
- (b) B_cert ≥ the certified lower bound at the four supports: **16/16**;
- (c) the trial's Rayleigh quotient ≥ the certified lower bound: **16/16**;
- (d) Theorem 1's symmetry ball encloses 0 at all 386 certified points;
- (e) Proposition H's explicit bound ≥ B_cert at all 286 Hermite points.

**The sharp prolate (Step 4, §9).**
- A C¹ edge restores finite A and B (Proposition P, proved).
- A taper on the prolate's own edge scale gives a trial whose RQ is close to our minima.
- Its Theorem-2 bound is **worse** than Kaiser–Bessel's at every tested point: 5–40× at the best taper width ([N], Table P). The taper pushes out-of-band mass to ω ≈ c²/κ and inflates A.
- At large x, Kaiser–Bessel's distance from the truth is mostly the trial's (§8). A sharper certified bound needs a better closed-form C¹ trial, not the prolate.

**Not addressed.** The lower half of the decay law is open, and for ζ it is at least as strong as RH. Every upper bound here is consistent with RH and proves nothing about zeros on the line.

## 1. Setting

- **Window.** x > 1, L = log x, W = [−L/2, L/2]. For real f ∈ L²(W) the sector is s ∈ {0, 1}, with f(−u) = (−1)^s f(u). The transform is F(z) = ∫ f(u) e^{izu} du.
- **L-functions.** 𝓛 is either ζ (q = 1, parity κ = 0, a pole), or L(s, χ) for a real primitive Dirichlet character χ of conductor q > 1 and parity κ (χ(−1) = (−1)^κ). The characters used are those of the fundamental discriminants D ∈ {5, 8, 12, 13, 17, −3, −4, −7, −8, −11, −15, −20}, with χ(n) = (D/n) the Kronecker symbol (`decay_law_mp.kronecker`) and q = |D|.
- **Zeros.** Nontrivial zeros are written ρ = ½ + iγ_ρ, so γ_ρ = −i(ρ − ½) ∈ ℂ and |Im γ_ρ| < ½. They are counted **with multiplicity**, the central zero included (m₀ = its multiplicity, any value).
- **The form.** Q(f) is Weil's quadratic form: the archimedean, prime and pole terms evaluated on the autocorrelation of f, which is exactly what our matrices compute (`decay_law_mp.zeros_side`). λ_s(x) = inf Q(f)/‖f‖² over real f ∈ L²(W) of sector s.
- **Explicit formula** (assumed, standard [S]: Weil 1952; `docs/RESEARCH_R1_EXCLUSION.md` §1). For real f supported in W and of bounded variation (in particular for f with finitely many jumps),

  Q(f) = Σ_ρ F(γ_ρ) F(−γ_ρ),

  with the sum over all nontrivial zeros, with multiplicity, absolutely convergent. Our matrices were validated against this zero sum in `docs/DECAY_LAW.md` (gates), and Zhu's normalisation of Q agrees with our with factor 1 (`docs/AUDIT_ZHU.md` §1) [S]. Apart from textbook facts (Poisson summation, Gauss sums, DLMF Bessel identities), this is the only analytic input taken without proof. The rigorous certificates (`docs/CERTIFICATE_238.md`) assume it too.
- **Decay-law variables.** c = 2πx/q, λ = √(x/q) (so c = 2πλ²).
- **Fourier on the line.** ψ̂(ξ) = ∫ ψ(y) e^{−2πiyξ} dy.

## 2. Zero counts: what the sources prove, and what R2 used

### 2.1 ζ

**R2 used** (`scripts/research_r2/r2lib.py`, `zeta_Nplus`):

  N_ζ(t) ≤ θ(t)/π + 1 + 0.112 log t + 0.278 log log t + 2.510, t ≥ e,

citing Trudgian, J. Number Theory 134 (2014) 280–292 as [A]. Here N_ζ(t) = #{0 < Im ρ ≤ t}, N(t) = 2N_ζ(t), and N_ζ = 0 below 14.134.

**What the sources say.**

- **Trudgian 2014** [R] (read in full: arXiv:1208.5846v2, the version on arXiv). Theorem 1 there reads |S(T)| ≤ 0.111 log T + 0.275 log log T + 2.450 for T ≥ e, and Corollary 1 is the N(T) form with +0.2/T₀. The journal constants 0.112, 0.278, 2.510 that R2 used are the published version's [S: the tables of HSW and Bellotti–Wong list Trudgian 2014 as (0.1120, 0.2780, 3.3850) in the N(T) form, i.e. 2.510 + 7/8]. I did not read the journal version.
- **The proof has a known gap** [R]. Three sources say the same thing:
  - Bennett–Martin–O'Bryant–Rechnitzer (arXiv:2005.02989 §1) explain the error in Trudgian's 2015 Dirichlet paper and add that the same difficulty recurs in the 2014 ζ paper (their ref. [11]). The error is in tracking the parameter constraints that Backlund's trick needs.
  - Hasanalizade–Shen–Wong (arXiv:2107.06506, footnote 2) say the same and state that their paper's purpose is "to fix the error occurring in [17]" (= Trudgian 2014).
  - Bellotti–Wong (arXiv:2412.15470v2, footnote 1) repeat it.
- **Bellotti–Wong also rescue the constants** (footnote 2) [R]. Their second estimate is sharper than Trudgian's for T > 387 899, and Platt's computation |S(T)| ≤ 2.5167 for 0 ≤ T ≤ 3.06·10¹⁰ (their eq. (1.2), [S] for the underlying computation) covers smaller T. Together these "assure" that Trudgian's bound and the results relying on it remain valid.

**Verdict on R2.**
- **The constants were not wrong.** R2 applied them only for T ≥ 14.134, inside the stated range T ≥ e, so no constant was used outside its stated range.
- **The citation was not a valid proof.** The cited proof is incomplete; validity rests on Bellotti–Wong's later argument.
- **Replacement.** The certified chain below uses Bellotti–Wong directly.

**Used here** [R: Bellotti–Wong, arXiv:2412.15470v2, Theorem 1.1, second estimate]. For every T ≥ e,

  |N_ζ(T) − (T/2π) log(T/2πe)| ≤ 0.11200 log T + 0.12567 log log T + 3.77417,

with N_ζ(T) = #{ρ: ζ(ρ) = 0, 0 < β < 1, 0 < γ ≤ T}. The count is the one produced by the argument principle in their proof, so it is with multiplicity [D: from the proof's first step].
- **Status of the source.** arXiv v2 (7 July 2025) was read. A web search summary states journal publication (Math. Comp. 2025) [A, not verified].
- **Its dependency at our heights** [R, §5 of their paper]. For e ≤ T ≤ 3.06·10¹⁰, Bellotti–Wong derive Theorem 1.1 from Platt's computed bound |S(T)| ≤ 2.5167 [S] and |N(T) − (T/2π)log(T/2πe) − 7/8| ≤ |S(T)| + 1/(50T), via their eq. (5.1). All heights used here (T_L ≤ 694) are in that range. I checked the step for the second estimate: 0.112 log T + 0.12567 log log T + 3.77417 ≥ 3.886 > 2.5167 + 1/(50e) + 7/8 ≈ 3.40 for T ≥ e [D].
- **Positions.** At the heights that matter here (T ≈ 14–2000) this bound lies 0.1–0.5 above R2's (`zero_count_check.json`), so the certified bounds are slightly weaker than R2's on this account.
- **Fallback with a refereed source.** Hasanalizade–Shen–Wong, J. Number Theory 235 (2022) [R: arXiv:2107.06506, Corollary 1.2]: 0.1038 log T + 0.2573 log log T + C₃ for T ≥ e, with C₃ = 9.4925. The printed C₃ = 9.3675 is a typo per Bellotti–Wong, footnote 1 [R]. Recomputed with that fallback, every ζ bound (both constructions, all ζ points: T_L from 15.6 to 694) would rise by a factor of 1.006–1.66, the largest at the smallest T_L. It is not used.

### 2.2 L(s, χ)

**Used here, and by R2** [R: Bennett–Martin–O'Bryant–Rechnitzer, arXiv:2005.02989 (Math. Comp. 90 (2021) 1455–1482), Theorem 1.1, read in full].

Restated (not quoted): let χ have conductor q > 1, T ≥ 5/7 and ℓ := log(q(T+2)/2π).
- If ℓ ≤ 1.567, then N(T, χ) = 0.
- If ℓ > 1.567, then |N(T, χ) − (T/π) log(qT/2πe) + χ(−1)/4| ≤ 0.22737ℓ + 2 log(1 + ℓ) − 0.5.

The definitions are in §1 of the paper:
- N(T, χ) counts the zeros with 0 < β < 1 and |γ| ≤ T, with multiplicity.
- Z(χ) excludes zeros on the imaginary axis, and Z(χ) = Z(χ*) for the inducing primitive character χ*.
- §2 assumes χ primitive for the derivation.

All characters here are primitive with q ∈ {3, 4, 5, 7, 8, 11, 12, 13, 15, 17, 20}.

**Verdict on R2.** The use was valid.
- R2 applied the inequality only where ℓ > 1.567 and T ≥ 5/7. Below that point it used the bound's value at t_min = max(5/7, 2πe^{1.5671}/q − 2), which is valid because N is nondecreasing.
- R2 cited the paper from its abstract only [A]. It is now read [R], and the constants and ranges are as R2 stated them.
- Here the case ℓ ≤ 1.567 is used as BMOR state it (N = 0 there), so the integral starts at t_low = max(5/7, 2πe^{1.567}/q − 2).

### 2.3 Lemma Z (closed-form tail integrals) [D]

Let N(t) = #{ρ: |Im ρ| ≤ t}, counted with multiplicity, and N⁺ ≥ N on [t_low, ∞) with N = 0 on [0, t_low). Take A, B > 0 and m(t) := min(B², A²/t²). Then

  Σ_ρ m(|Im ρ|) ≤ 2A² ∫_{T_L}^∞ N⁺(t) t⁻³ dt,  T_L := max(A/B, t_low).

*Proof.* N is right-continuous and nondecreasing, m is continuous, nonincreasing and Lipschitz, and m(t)N(t) → 0. By Stieltjes integration by parts,

  Σ_ρ m(|Im ρ|) = ∫_{[0,∞)} m dN = ∫_0^∞ N(t)(−m′(t)) dt,

and −m′(t) = 2A²t⁻³·1{t > A/B} ≥ 0. Then N = 0 below t_low, N ≤ N⁺ above it, and −m′ ≥ 0. A lower limit below the true A/B only adds a nonnegative amount, which is why `cert_lib.zero_sum_bound` uses the lower endpoint of the ball A/B. ∎

The integrals have closed forms (`cert_lib.zeta_tail_integral`, `dirichlet_tail_integral`), all inequalities pointing up:

- **ζ** (N⁺ = 2·[(t/2π)log(t/2πe) + a log t + b log log t + c₃], with (a, b, c₃) from §2.1, T ≥ t_low = 14):

  ∫_T^∞ N⁺ t⁻³ ≤ 2[ log(T/2π)/(2πT) + a(2 log T + 1)/(4T²) + b(log log T/(2T²) + 1/(4T² log T)) + c₃/(2T²) ].

  The log log term is integrated by parts, and ∫_T^∞ dt/(2t³ log t) ≤ 1/(4T² log T).

- **L(s, χ)** (N⁺ the BMOR right side, T ≥ t_low). Set ℓ₀ = ℓ(T) and I_ℓ := ∫_T^∞ ℓ(t)t⁻³dt = ℓ₀/(2T²) + 1/(4T) − log(1 + 2/T)/8, which is exact (by parts, then partial fractions). Then

  ∫_T^∞ N⁺ t⁻³ ≤ log(qT/2π)/(πT) − (χ(−1)/4 + ½)/(2T²) + 0.22737·I_ℓ + 2[(log(1+ℓ₀) − ℓ₀/(1+ℓ₀))/(2T²) + I_ℓ/(1+ℓ₀)].

  The last bracket uses the tangent-line bound log(1+ℓ) ≤ log(1+ℓ₀) + (ℓ − ℓ₀)/(1+ℓ₀) (concavity).

**Checks that can fail** (`zero_count_check.py` → `zero_count_check.json`) [N]:
- **N⁺(t) ≥ N(t) at every zero and just after it.**
  - ζ: the first 200 ordinates (mpmath `zetazero`, heights up to 396), 0 failures, minimum slack 5.4.
  - χ₅, χ₈, χ₁₂, χ₁₃, χ₁₇: all zeros to height 60 found from sign changes, 0 failures, minimum slack 1.0.
- **Closed forms against quadrature.** Every closed form is ≥ the mpmath quadrature of the same N⁺, by at most 0.24 % (from the concavity and by-parts steps).

## 3. Theorem 1 (the E-map identity) [D]

**Class 𝒜_χ(μ).** ψ: ℝ → ℝ such that:
- (A1) ψ is continuous and locally of bounded variation, with ψ(−y) = (−1)^κ ψ(y);
- (A2) |ψ(y)| ≤ C(1+|y|)^{−1−δ} for some C, δ > 0;
- (A3) ψ̂ = μψ (necessarily μ⁴ = 1 and μ² = (−1)^κ);
- (A4) if q = 1, ψ(0) = 0.

Put E_χψ(w) = w^{1/2} Σ_{n≥1} χ(n) ψ(nw/√q), f_ψ(u) = E_χψ(e^u), σ = (i^κ μ)⁻¹ ∈ {±1}, f_win = f_ψ·1_W, f_out = f_ψ − f_win, and F_ψ, F_out the transforms.

**Theorem 1.**
- (i) f_ψ(−u) = σ f_ψ(u).
- (ii) |f_ψ(u)| ≤ C′e^{−(½+δ)|u|}.
- (iii) On |Im z| < ½ + δ, F_ψ(z) = q^{(½+iz)/2} L(½+iz, χ) M_ψ(z) with M_ψ(z) = ∫₀^∞ ψ(y) y^{−½+iz} dy. Hence F_ψ(γ_ρ) = 0 at every nontrivial zero, on or off the line, of any multiplicity.
- (iv) Q(f_win) = σ Σ_ρ F_out(γ_ρ)², the sum over all nontrivial zeros with multiplicity.

*Proof.* This is R2 §3.2, unchanged; the ingredients are listed here for the auditor.

1. **Twisted Poisson summation** for primitive χ: Σ_{n∈ℤ} χ(n)g(n) = (τ(χ)/q) Σ_m χ̄(m)ĝ(m/q), for g continuous with g, ĝ = O((1+|y|)^{−1−δ}) [S, standard; a three-line derivation from Poisson summation and the Gauss sum Σ_a χ(a)e(am/q) = χ̄(m)τ(χ), valid for all m when χ is primitive].
2. **Gauss sum.** τ(χ) = i^κ√q for real primitive χ [S, Gauss].
3. **Applying it.** With g(y) = ψ(yw/√q) this gives E_χψ(w) = i^κμ·E_χψ(1/w). For q = 1 the extra Poisson term (μψ(0)/w − ψ(0))/2 vanishes by (A4). That is (i).
4. **(ii)** follows from (A2) for w ≥ 1 and from (i) for w ≤ 1.
5. **(iii)** follows from Fubini where the Dirichlet series converges (−½ − δ < Im z < −½) and analytic continuation. For ζ, the pole at z = −i/2 meets M_ψ(−i/2) = ∫₀^∞ψ = ψ̂(0)/2 = μψ(0)/2 = 0.
6. **(iv).** F_win = F_ψ − F_out on the strip, so F_win(γ_ρ) = −F_out(γ_ρ). Since −γ_ρ = γ_{1−ρ} is also a zero, F_out(−z) = σF_out(z), and the explicit formula of §1 applies to f_win. ∎

**Checks** [C], two that can fail, run at every certified point:
- **(i) as a ball.** E(1/w₀) − σE(w₀), with all n-tails bounded, must enclose 0 (w₀ = x^{1/4}). With a wrong σ or a wrong ψ̂ it would not.
- **ψ(0) as a ball** must enclose 0 for ζ.

Results are in §7.

## 4. Theorem 2 (unconditional bound) [D]

**Hypotheses.** ψ ∈ 𝒜_χ(μ), and f := f_ψ is absolutely continuous on [L/2, ∞) with ∫_{L/2}^∞ |f′| e^{u/2} du < ∞. Let A, B be any numbers with

  B ≥ B₀ := 2∫_{L/2}^∞ |f| e^{u/2} du = 2∫_{√x}^∞ |S(w)| dw,
  A ≥ A₀ := 2(|f(L/2)| e^{L/4} + ∫_{L/2}^∞ |f′| e^{u/2} du),  with |f(L/2)|e^{L/4} = √x |S(√x)|,

where S(w) = Σ χ(n)ψ(nw/√q). Let N⁺ and t_low be as in §2.3.

**Theorem 2.** In the sector (−1)^s = σ,

  λ_s(x) ≤ Q(f_win)/‖f_win‖² ≤ 2A² ∫_{max(A/B, t_low)}^∞ N⁺(t) t⁻³ dt / ‖f_win‖²  =: B(x).

*Proof.*
1. **The bound by B.** F_out(z) = ∫_{L/2}^∞ f(u)(e^{izu} + σe^{−izu}) du. For |Im z| < ½ and u > 0, |e^{±izu}| ≤ e^{u/2}; this is the only property of the zeros used, 0 < Re ρ < 1. It gives |F_out(γ_ρ)| ≤ B₀.
2. **The bound by A.** Integrating by parts, with the boundary term at ∞ vanishing by Theorem 1(ii), gives |F_out(z)| ≤ A₀/|z| ≤ A₀/|Re z|, and |Re γ_ρ| = |Im ρ|.
3. **Summing.** By Theorem 1(iv), |Q(f_win)| ≤ Σ_ρ min(B², A²/|Im ρ|²). Lemma Z bounds the sum. Since f_win is in the sector, λ_s(x) ≤ Q(f_win)/‖f_win‖². ∎

**Remarks.**
- **Central zeros, multiplicity.** Nothing assumes simple zeros or m₀ ≤ 1. A central zero (or any real zero) has Im ρ = 0 and contributes m(0) = B² for each unit of multiplicity; that is included because N(t) ≥ N(0) counts it. F_ψ vanishes there to at least the multiplicity of the zero.
- **What is unconditional.** Theorems 1–2 use no RH, no zero-density hypothesis and no Euler product. The inputs are the explicit formula (§1) and the zero counts (§2), which are themselves unconditional theorems.
- **Inequality directions.** An upper bound for A, B and N⁺ and a lower bound for ‖f_win‖² give an upper bound for B(x). The code takes the upper endpoints of the balls for A and B, the lower endpoint for ‖f_win‖² and for A/B, and the upper endpoint of the final ball.

## 5. Corollary 3 (Hermite trials): all constants explicit [D + C]

**Trial.** Let n = κ + 2s and h_k(y) = H_k(√(2π)y)e^{−πy²}, so that ĥ_k = (−i)^k h_k.
- If q > 1, take ψ = h_n.
- For ζ, take ψ = H_{n+4}(0)h_n − H_n(0)h_{n+4}, which has ψ(0) = 0 exactly (integer coefficients). For s = 0 this is Riemann's h, and F_ψ ∝ Ξ.

Write ψ(y) = P(z)e^{−z²/2} with z = √(2π)y, P an integer polynomial of degree d, and set Q_P(z) := z(P′(z) − zP(z)), so that yψ′(y) = Q_P(z)e^{−z²/2}. Then ψ ∈ 𝒜_χ(μ) with μ = (−i)^n (Gaussian decay), and f is smooth, so Theorem 2 applies.

**Lemma G** [D]. For Y > 0 and a polynomial R, set I_R(Y) := ∫_Y^∞ |R(√(2π)y)| e^{−πy²} dy and G_k(Z) := ∫_Z^∞ z^k e^{−z²/2} dz = 2^{(k−1)/2} Γ((k+1)/2, Z²/2), with Z = √(2π)Y.
- If the Taylor coefficients of R(Z + t) in t are all ≥ 0, with the constant one > 0 (or all ≤ 0, with the constant one < 0), then R has no root on [Z, ∞) and I_R(Y) = (2π)^{−1/2} |Σ_k r_k G_k(Z)|.
- In every case I_R(Y) ≤ (2π)^{−1/2} Σ_k |r_k| G_k(Z).

**Bounds** [D]:
- **B.** |S(w)| ≤ Σ_{χ(n)≠0} |ψ(nw/√q)| and the substitution y = nw/√q give

  B₀ ≤ 2√q Σ_{n: χ(n)≠0} I_P(nλ)/n.

- **A.** wS′(w) = Σ χ(n)(nw/√q)ψ′(nw/√q), and ∫|S/2 + wS′| ≤ B₀/4 + Σ∫|…| give

  A₀ ≤ 2(√x|S(√x)| + B/4 + √q Σ_{n: χ(n)≠0} I_{Q_P}(nλ)/n).

- **Norm.** ‖f_win‖² = 2∫_1^{√x} S(w)² dw. This uses Theorem 1(i) and f²du = S²dw.

**Certified evaluation** (`cert_lib.hermite_quantities`):
- **Sums over n.** The sums run to N = ⌈12/λ⌉ + 2, with each I evaluated by Arb's `gamma_upper`. The tail n > N is bounded via |R(√(2π)y)| ≤ C_R y^{deg R} for y ≥ 1, ∫_Y^∞ y^d e^{−πy²} ≤ Y^{d−1}e^{−πY²}/π (2πY² ≥ 2(d−1)), and a geometric series with explicit ratio (`gauss_tail`).
- **Edge.** S(√x) is the exact signed sum, plus a Gaussian tail ball.
- **Norm.** S is truncated at N_s = ⌈9√q⌉ + 1, and ∫_1^X S_N² is evaluated by `acb.integral` (S_N is a finite sum of entire functions; X ≤ √x is an exact dyadic). Then ‖S‖ ≥ ‖S_N‖ − ε_S√(X−1), with ε_S ≥ sup_{w≥1}|S − S_N| from the same tail bound, which uses that y^d e^{−πy²} decreases for y ≥ 9.

Every out-of-band integral at every point used the exact (root-free) branch of Lemma G (§7).

### 5.1 Proposition H (the Hermite rate, for all x) [D + C]

**Proposition H.** Fix the object (q, χ) and the sector s, let ψ be the Hermite trial of §5 and d = deg P. For every x ≥ x₀ := 2q,

  λ_s(x) ≤ K_H c^{d+1} e^{−c},  c = 2πx/q,

with the certified constants K_H of the table below.

*Proof.* Put λ = √(x/q) ≥ λ₀ := √2, Z = √(2π)Y and Z₀² = 2πλ₀² = 4π.
1. **Norm.** S(w) does not depend on x, and the window grows with x, so ‖f_win(x)‖² = ∫_{−L/2}^{L/2} f² is nondecreasing in x. Hence ‖f_win(x)‖² ≥ N₀ := ‖f_win(x₀)‖², whose lower endpoint is certified (`hermite_cert.json`, x/q = 2).
2. **Incomplete Gamma.** For X > s − 1 ≥ 0, Γ(s, X) = X^{s−1}e^{−X}∫₀^∞(1 + u/X)^{s−1}e^{−u}du ≤ X^{s−1}e^{−X}/(1 − (s−1)/X), and for s ≤ 1, Γ(s, X) ≤ X^{s−1}e^{−X}. Hence G_k(Z) ≤ g_k Z^{k−1}e^{−Z²/2} with g_k = (1 − (k−1)/Z₀²)⁻¹ for k ≥ 2 (k ≤ 13) and g_k = 1 for k ≤ 1. So for Y ≥ λ₀ and any polynomial R of degree d_R,

   I_R(Y) ≤ K_R Y^{d_R−1}e^{−πY²},  K_R := (2π)^{−1/2} Σ_k |r_k| g_k (2π)^{(k−1)/2} λ₀^{k−d_R}.

3. **Sums over n.** With Θ_j := Σ_{n≥1} n^j e^{−πλ₀²(n²−1)} (nonincreasing in λ, so it bounds the same sum at λ) and C′_P := Σ_k |p_k|(2π)^{k/2}λ₀^{k−d}, the bounds of §5 give
   - B ≤ 2√q K_P Θ_{d−2} λ^{d−1}e^{−πλ²};
   - √x|S(√x)| ≤ √q C′_P Θ_d λ^{d+1}e^{−πλ²};
   - √qΣ_n I_{Q_P}(nλ)/n ≤ √q K_{Q_P} Θ_d λ^{d+1}e^{−πλ²} (deg Q_P = d + 2).

   So A ≤ K_A λ^{d+1}e^{−πλ²}, with K_A := 2√q(C′_P Θ_d + K_P Θ_{d−2}/(2λ₀²) + K_{Q_P}Θ_d).
4. **Zero sum.** The integral is decreasing in its lower limit and T_L ≥ t_low, so it is ≤ 2A²·I(t_low), with I(t) := ∫_t^∞ N⁺ t⁻³ (Lemma Z).
5. **Conclusion.** Theorem 2 then gives λ_s(x) ≤ 2K_A² I(t_low) λ^{2d+2}e^{−2πλ²}/N₀ = K_H c^{d+1}e^{−c}, with K_H := 2K_A² I(t_low)/(N₀(2π)^{d+1}). Every constant is an Arb ball and the upper endpoint is used (`hermite_rate.py`). ∎

| object | x₀ = 2q | K_H (even sector), d | K_H (odd sector), d | min over grid of (K_H c^{d+1}e^{−c})/B_cert, even / odd |
|---|---|---|---|---|
| ζ | 2 | 2.2799342, d = 4 | 0.72589493, d = 6 | 7.26 / 41.6 |
| χ₋₃ | 6 | 10.547339, d = 1 | 4.6487397, d = 3 | 2.61 / 5.26 |
| χ₋₄ | 8 | 11.603261, d = 1 | 7.2755526, d = 3 | 3.65 / 7.31 |
| χ₅ | 10 | 27.873679, d = 0 | 14.578742, d = 2 | 4.12 / 6.01 |
| χ₋₇ | 14 | 20.992526, d = 1 | 22.026738, d = 3 | 9.58 / 19 |
| χ₈ | 16 | 48.495196, d = 0 | 39.426055, d = 2 | 10.4 / 15 |
| χ₋₈ | 16 | 40.205760, d = 1 | 46.940970, d = 3 | 13.6 / 26.9 |
| χ₋₁₁ | 22 | 278.33060, d = 1 | 182.54886, d = 3 | 51.5 / 102 |
| χ₁₂ | 24 | 171.11597, d = 0 | 197.74493, d = 2 | 40 / 57.8 |
| χ₁₃ | 26 | 266.29279, d = 0 | 299.62875, d = 2 | 41.5 / 60 |
| χ₋₁₅ | 30 | 121.86517, d = 1 | 228.08302, d = 3 | 59.6 / 117 |
| χ₁₇ | 34 | 130.71433, d = 0 | 177.45398, d = 2 | 46 / 66.3 |
| χ₋₂₀ | 40 | 228.38271, d = 1 | 433.24019, d = 3 | 63.6 / 125 |

(`data/research_r2_rigorous/hermite_rate.json`; check 286/286.)

**Check that can fail.** The chain above dominates the certified one term by term, so K_H c^{d+1}e^{−c} ≥ B_cert(x) must hold at every grid point. It does at **286/286**; the smallest ratio per object and sector is 2.6–125. An error in either chain would show as a violation.

**Rate.** That is e^{−c} up to the polynomial c^{d+1}, which is half the conjectured rate. It is the elementary version, now with a proof valid for all x.

## 6. Theorem 4 (Kaiser–Bessel trials): all constants explicit [D + C]

**Trial.** Fix an integer m ≥ 2 (m = 2 throughout) and 0 < β ≤ c (β = c or c − 2; the smaller certified bound is kept). Put ν = m + ½ and

  w(t) = (1 − t²)^{m/2} I_m(β√(1 − t²)) = (β/2)^m (1 − t²)^m ₀F̃₁(; m+1; β²(1−t²)/4),  |t| ≤ 1,

where ₀F̃₁ is the regularized ₀F₁, so w is the restriction of an entire function. Then:
- φ(y) = (1 + a t²) w(t) for κ = 0 and φ(y) = t·w(t) for κ = 1, with t = y/λ and φ = 0 for |y| ≥ λ;
- a = 0 except for ζ, where a is the unique real number with ψ(0) = 0 (formula below);
- ψ = (φ + μ̄φ̂)/2 with μ = (−i)^{κ+2s}.

**Lemma K** (the transform) [D]. For ω ∈ ℂ,

  ∫_{−1}^{1} w(t) e^{−iωt} dt = √(2π) β^m K_ν(ω),  K_ν(ω) := 2^{−ν} ₀F̃₁(; ν+1; (β² − ω²)/4),

which is entire, with K_ν′(ω) = −ωK_{ν+1}(ω).

*Proof.*
1. **Expand w.** w(t) = Σ_k (β/2)^{m+2k}(1−t²)^{m+k}/(k!(m+k)!).
2. **One term.** ∫_{−1}^{1}(1−t²)^{μ}e^{−iωt}dt = √π Γ(μ+1) Σ_j (−ω²/4)^j/(j! Γ(μ+j+3/2)) [the standard Poisson integral for J_{μ+½}].
3. **Collect.** Grouping the double sum by k + j and using the binomial theorem on ((β² − ω²)/4)^{k+j} gives √π(β/2)^m Σ_r ((β²−ω²)/4)^r/(r! Γ(m+r+3/2)) = √(2π)β^m K_ν(ω).
4. **Derivative.** Termwise.

The t² and t moments follow by differentiating under the integral: ∫t²w e^{−iωt} = √(2π)β^m(K_{ν+1} − ω²K_{ν+2}) and ∫t w e^{−iωt} = −iω√(2π)β^m K_{ν+1}. R2 cited this as Lewitt 1990 [S] and checked it by quadrature; here it is derived. ∎

**Consequently.** With pref = λ√(2π)β^m and ω = 2πλy:
- φ̂(y) = pref·Σ_i c_i ω^{p_i}K_{ν+j_i}(ω), times −i for κ = 1;
- the (c, p, j) terms are (1, 0, 0), (a, 0, 1), (−a, 2, 2) for κ = 0 (only the first when a = 0), and (1, 1, 1) for κ = 1;
- ψ = (φ + (−1)^s φ̂_r)/2 with φ̂_r real;
- for ζ, a = −(I_m(β) + (−1)^s pref K_ν(0))/((−1)^s pref K_{ν+1}(0)).

**ψ ∈ 𝒜_χ(μ).**
- (A1): φ(±λ) = 0 because w(±1) = 0 (m ≥ 1). φ is C¹ on ℝ with φ″ of bounded variation.
- (A3): φ̂̂ = φ(−·) and μ² = (−1)^κ.
- (A4): by the choice of a. The certified ball ψ(0) encloses 0 at every ζ point (§7).
- (A2) with δ = 2: from Lemma E below.

**Lemma E** (envelopes) [D from R sources]. Let ω ≥ β, s = √(ω² − β²) and k + ½ = ν + j. Then |K_{k+½}(ω)| = |J_{k+½}(s)/s^{k+½}| ≤ min(cap_k, hank_k(s)), where:
- **cap_k** = 2^{−k−½}/Γ(k+3/2) [R: DLMF 10.14.4, |J_ν(z)| ≤ |z/2|^ν e^{|Im z|}/Γ(ν+1), ν ≥ −½];
- **hank_k(s)** = √(2/π) Σ_{j=0}^{k} a_j(k) s^{−k−1−j}, with a_j(k) = (k+j)!/(2^j j!(k−j)!) [R: DLMF 10.49.1–2, the terminating expansion of j_k, with J_{k+½}(s) = √(2s/π) j_k(s)].

For each term (c, p, j), with k = m + j, p − k − 1 ≤ −3 (−3 for (1, 0, 0), (−a, 2, 2), (1, 1, 1); −4 for (a, 0, 1)), so |ψ(y)| = O(y^{−3}).

**Lemma H** (the n-sums without monotonicity) [D]. For g ≥ 0 on [c, ∞),

  Σ_{n≥1} (1/n) ∫_{nc}^∞ g(ω) dω = ∫_c^∞ g(ω) H_{⌊ω/c⌋} dω ≤ ∫_c^∞ g(ω)(1 + log(ω/c)) dω.

The equality is Tonelli; the inequality is H_k ≤ 1 + log k. **This replaces R2's step "env·y^{m+1} is nonincreasing, checked on an 80-point grid", which was not proved**: no monotonicity is used anywhere now.

**Lemma T** (explicit power tails) [D]. For ω ≥ Ω₀ > β, s ≥ r₀ω with r₀ = √(1 − β²/Ω₀²). Hence ω^p hank_k(s) ≤ C·ω^{p−k−1}, with explicit C, and every tail integral and tail sum beyond Ω₀ has a closed form (`cert_lib.C_E`, `env_G`).

**Bounds** [D]. With E(ω) := Σ_i |c_i| ω^{p_i} env_{k_i}(ω) ≥ |Σ_i c_i ω^{p_i}K| (env = the Lemma E minimum), |ψ(y)| ≤ (pref/2)E(2πλy) for y ≥ λ. Lemma H then gives:

  B₀ ≤ (√q β^m/√(2π)) · J₀,  J₀ := ∫_c^∞ E(ω)(1 + log(ω/c)) dω,
  A₀ ≤ 2( √x|S(√x)| + B/4 + (√q β^m/(2√(2π))) · J₁ ),  J₁ := ∫_c^∞ ω E₁(ω)(1 + log(ω/c)) dω,

where ωE₁(ω) = Σ_i |c_i|(p_i ω^{p_i} env_{k_i} + ω^{p_i+2} env_{k_i+1}) bounds ω·|d/dω Σ c_i ω^{p_i}K|.

**Certified evaluation** (`cert_lib.kb_quantities`).

- **J₀ and J₁.** Each term ∫_c^∞ ω^p env_k (1 + log(ω/c)) is split into three pieces:
  - cap on [c, ω*], in closed form;
  - hank on [ω*, Ω], by `acb.integral` with the analytic flag on √ and log;
  - the closed-form tail on [Ω, ∞) from Lemma T.

  The split point ω* (≥ c and > β) is chosen in floating point as the crossover cap = hank. Any choice is valid, because both bounds hold everywhere.

- **The edge value S(√x) = Σ_n χ(n)ψ(nλ).** The terms are summed exactly as balls for n ≤ 200. Out of band K is evaluated as J_{k+½}(s)/s^{k+½} through Arb's Bessel J; Arb's ₀F₁ series loses all precision near argument −10⁴ (tested). The remainder is bounded by (pref/2)C c⁻³/(2·200²) via Lemma T.

- **The norm.** On each band interval I_k = [√x/(k+1), √x/k] ∩ [1, √x], the in-band set is {n ≤ k}. On it, S = S_main + S_rest, where:
  - S_main = Σ_{n≤k} χψ_full + Σ_{k<n≤k+M} χψ_out is a finite sum of entire functions, integrated by `acb.integral` over exact dyadic subintervals of I_k;
  - sup|S_rest| ≤ ε_k, via Lemma T, using ω_n ≥ c·n·a_k/√x.

  Then ‖S‖²_{I_k} ≥ (‖S_main‖ − ε_k|I_k|^{1/2})². M is the smallest of (0, 3, 8, 20, 50) for which this costs < 10⁻⁶ of the interval norm.

- **Ranges and branches.**
  - All sums over n include χ(n) = 0 terms only where that is harmless (the B and A bounds over-count them).
  - The band breakpoint ball must not straddle w = 1; this is asserted.
  - β = c is handled without comparing the ball β with the ball c.

### 6.1 The Kaiser–Bessel rate: Conjecture R (not proved)

**Conjecture R.** For each degree-1 object (ζ or real primitive χ) and each sector s there are x₀, C and p such that the Kaiser–Bessel bound of Theorem 4 (m = 2, β = c − 2) satisfies

  B(x) ≤ C c^p e^{−2c}  for all x ≥ x₀.

In particular λ_s(x) ≤ C c^p e^{−2c}.

**What supports it** [C, finite]. On R2's grid, x/q = 2..16 (c = 12.6..100.5), ln B_cert + 2c fits α + p ln c with:

| class | p | max residual (in ln) |
|---|---|---|
| χ even, n = 0 | 3.45–3.48 | ≤ 0.05 |
| χ odd, n = 1 | 4.98–5.10 | ≤ 0.08 |
| χ even, n = 2 | 7.16–7.21 | ≤ 0.04 |
| χ odd, n = 3 | 7.73–7.84 | ≤ 0.28 |
| ζ, n = 4 / 6 | 9.26 / 11.5 | 0.13 / 0.47 |

These are the KB table of §7. A table of finitely many certified values does not prove a rate.

**What is missing.** Theorem 2 gives B(x) ≤ 2A²I(t_low)/‖f_win‖² for every x, exactly as in Proposition H. So Conjecture R follows from two estimates, uniform in x ≥ x₀.

- **(M1) A lower bound for the norm** (the hard part).

  ‖f_win‖² ≥ e^{2β} c^{−p₂}/C₂.

  Unlike the Hermite case, ψ depends on x through λ and β, so monotonicity in x is not available. The e^{2β} comes from w ∈ [1, O(1)], where many in-band terms contribute. There ψ(y) ≈ e^β (2πλ)⁻¹·(a Hermite function of y) up to relative O(1/c), and S(w) ≈ e^β(2πλ)⁻¹·(a fixed theta-type sum). A proof needs:
  - uniform asymptotics, with explicit error bounds, of I_m(β√(1−t²)) and of K_ν(ω) = I_ν(z)/z^ν for β → ∞ (e.g. Olver's uniform expansions with the error bounds of DLMF §10.41(iv)). This must go to second order in the sector where (−1)^s = −1, because there the leading Gaussians of φ and φ̂ cancel in ψ = (φ − φ̂)/2;
  - a positive lower bound for ∫₁² (Σ_n χ(n) h(nw/√q))² dw for the limiting Hermite function h. That is a single certifiable number per object and sector.
- **(M2) Polynomial upper bounds for the out-of-band quantities.**

  A, B ≤ C₁ c^{p₁}.

  The envelope integrals J₀ and J₁ are polynomial in c by Lemmas E and T: with β = c − 2, s ≥ ω√(1 − β²/c²) ≥ ω·√(4(c−1))/c on [c, ∞). Writing this out also needs:
  - a uniform bound |a| ≤ C c for the ζ coefficient a, a ratio of Bessel values, via standard bounds on I_{ν+1}/I_ν;
  - the edge sum.

  These are routine but not written out here, so M2 is also listed as missing.

With M1 and M2, Theorem 2 gives Conjecture R with p = 2p₁ + p₂ immediately. The constant is then explicit, as in Proposition H.

## 7. Theorem C: the certified numbers [C]

**Theorem C (finite).** For every (𝓛, x, s) listed in Tables S, Hermite and KB below, λ_s(x) ≤ B_cert(x). This is a statement about those x only.

Every B_cert below is the upper endpoint of an Arb ball (256-bit working precision), rounded up in the last printed digit. The data are `data/research_r2_rigorous/{hermite,kb}_cert{,_supports}.json`, printed by `summarize.py`. "cert" is our certified lower bound at the same x (`docs/CERTIFICATE_238.md`, re-read from `data/connes/*certificate*.json` as the lower endpoint of the stored ball).

**Table S.** ζ at the four certified supports (B_cert = certified upper bound; cert = certified lower bound)

| support 2a | x | sector | cert (lower) | Hermite B_cert | KB B_cert | KB β | KB RQ (a) | Hermite RQ (a) | KB B_cert/cert |
|---|---|---|---|---|---|---|---|---|---|
| 1.6 | 4.9530 | even | 1.0277e-17 | 1.657e-7 | 2.937e-13 | c−2 | 5.37e-16 | 8.93e-9 | 2.86e+4 |
| 1.6 | 4.9530 | odd | 9.1183e-15 | 1.2e-5 | 1.045e-10 | c−2 | 2.75e-13 | 8.39e-7 | 1.15e+4 |
| 2.38 | 10.8049 | even | 6.8131e-48 | 5.381e-22 | 4.467e-42 | c−2 | 1.77e-44 | 2.07e-23 | 6.56e+5 |
| 2.38 | 10.8049 | odd | 3.9197e-44 | 2.163e-19 | 6.152e-39 | c−2 | 2.24e-41 | 1.37e-20 | 1.57e+5 |
| 2.6 | 13.4637 | even | 5.7634e-62 | 7.728e-29 | 9.403e-56 | c−2 | 2.8e-58 | 3.45e-30 | 1.63e+6 |
| 2.6 | 13.4637 | odd | 5.3601e-58 | 4.94e-26 | 1.891e-52 | c−2 | 5.68e-55 | 2.34e-27 | 3.53e+5 |
| 2.99 | 19.8857 | even | 3.5011e-96 | 1.233e-45 | 3.448e-89 | c−2 | 1.57e-91 | 4.67e-47 | 9.85e+6 |
| 2.99 | 19.8857 | odd | 8.2563e-92 | 1.773e-42 | 1.257e-85 | c−2 | 4.67e-88 | 7.65e-44 | 1.52e+6 |

**Table Hermite** (`hermite_cert.json`; ranges over x/q, 11 values each; the p-fit column fits ln B_cert + c = α + p ln c)

| object | sector | n | x/q | B_cert range | B_cert / R2 float | ln B_cert + c range | p of C c^p fit (max resid) | sym ⊇ 0 |
|---|---|---|---|---|---|---|---|---|
| ζ | even | 4 | 2..16 | 3.43e-1 .. 1.95e-35 | 1.01 .. 1.08 | 11.5 .. 20.6 | 4.38 (0.026) | True |
| ζ | odd | 6 | 2..16 | 3.01 .. 1.79e-32 | 1.01 .. 1.1 | 13.7 .. 27.4 | 6.62 (0.077) | True |
| χ₋₃ | even | 1 | 2..16 | 2.22e-3 .. 1.79e-40 | 1.01 .. 1.05 | 6.46 .. 9.01 | 1.23 (0.0065) | True |
| χ₋₃ | odd | 3 | 2..16 | 7.69e-2 .. 4.41e-37 | 1.01 .. 1.06 | 10.0 .. 16.8 | 3.28 (0.028) | True |
| χ₋₄ | even | 1 | 2..16 | 1.75e-3 .. 1.35e-40 | 1.01 .. 1.05 | 6.22 .. 8.72 | 1.21 (0.004) | True |
| χ₋₄ | odd | 3 | 2..16 | 8.66e-2 .. 4.73e-37 | 1.01 .. 1.06 | 10.1 .. 16.9 | 3.25 (0.022) | True |
| χ₅ | even | 0 | 2..16 | 2.97e-4 .. 2.93e-42 | 1.01 .. 1.05 | 4.44 .. 4.9 | 0.217 (0.01) | True |
| χ₅ | odd | 2 | 2..16 | 1.68e-2 .. 1.03e-38 | 1.01 .. 1.05 | 8.48 .. 13.1 | 2.2 (0.0054) | True |
| χ₋₇ | even | 1 | 2..16 | 1.21e-3 .. 8.69e-41 | 1.01 .. 1.05 | 5.85 .. 8.29 | 1.18 (0.0027) | True |
| χ₋₇ | odd | 3 | 2..16 | 1.01e-1 .. 5.12e-37 | 1.01 .. 1.06 | 10.3 .. 17.0 | 3.22 (0.015) | True |
| χ₈ | even | 0 | 2..16 | 2.05e-4 .. 1.93e-42 | 1.01 .. 1.05 | 4.07 .. 4.48 | 0.194 (0.0066) | True |
| χ₈ | odd | 2 | 2..16 | 1.81e-2 .. 1.05e-38 | 1.01 .. 1.05 | 8.56 .. 13.1 | 2.18 (0.0023) | True |
| χ₋₈ | even | 1 | 2..16 | 1.63e-3 .. 1.16e-40 | 1.01 .. 1.05 | 6.15 .. 8.57 | 1.17 (0.0037) | True |
| χ₋₈ | odd | 3 | 2..16 | 1.52e-1 .. 7.6e-37 | 1.01 .. 1.06 | 10.7 .. 17.4 | 3.21 (0.014) | True |
| χ₋₁₁ | even | 1 | 2..16 | 2.98e-3 .. 2.06e-40 | 1.01 .. 1.05 | 6.75 .. 9.15 | 1.16 (0.0055) | True |
| χ₋₁₁ | odd | 3 | 2..16 | 1.56e-1 .. 7.57e-37 | 1.01 .. 1.06 | 10.7 .. 17.4 | 3.2 (0.011) | True |
| χ₁₂ | even | 0 | 2..16 | 1.87e-4 .. 1.7e-42 | 1.01 .. 1.04 | 3.98 .. 4.35 | 0.177 (0.0044) | True |
| χ₁₂ | odd | 2 | 2..16 | 2.37e-2 .. 1.32e-38 | 1.01 .. 1.05 | 8.82 .. 13.3 | 2.16 (0.0026) | True |
| χ₁₃ | even | 0 | 2..16 | 2.81e-4 .. 2.53e-42 | 1.01 .. 1.04 | 4.39 .. 4.75 | 0.174 (0.0042) | True |
| χ₁₃ | odd | 2 | 2..16 | 3.46e-2 .. 1.92e-38 | 1.01 .. 1.05 | 9.2 .. 13.7 | 2.16 (0.0031) | True |
| χ₋₁₅ | even | 1 | 2..16 | 1.13e-3 .. 7.59e-41 | 1.01 .. 1.05 | 5.78 .. 8.15 | 1.14 (0.007) | True |
| χ₋₁₅ | odd | 3 | 2..16 | 1.69e-1 .. 7.98e-37 | 1.01 .. 1.05 | 10.8 .. 17.4 | 3.18 (0.0092) | True |
| χ₁₇ | even | 0 | 2..16 | 1.25e-4 .. 1.1e-42 | 1.01 .. 1.04 | 3.58 .. 3.92 | 0.165 (0.0034) | True |
| χ₁₇ | odd | 2 | 2..16 | 1.85e-2 .. 1.01e-38 | 1.01 .. 1.05 | 8.58 .. 13.0 | 2.15 (0.0045) | True |
| χ₋₂₀ | even | 1 | 2..16 | 1.98e-3 .. 1.31e-40 | 1.01 .. 1.05 | 6.34 .. 8.69 | 1.13 (0.0081) | True |
| χ₋₂₀ | odd | 3 | 2..16 | 3.02e-1 .. 1.39e-36 | 1.01 .. 1.05 | 11.4 .. 18.0 | 3.17 (0.0073) | True |

**Table KB** (`kb_cert.json`; x/q ∈ {2, 4, 6, 8, 10, 16}; the better of β = c and c − 2, which is c − 2 at all but a few small-x points; the p-fit column fits ln B_cert + 2c = α + p ln c)

| object | sector | n | x/q | B_cert range | B_cert / R2 float | ln B_cert + 2c range | p of C c^p fit (max resid) | sym ⊇ 0 |
|---|---|---|---|---|---|---|---|---|
| ζ | even | 4 | 2..16 | 7.84e-1 .. 7.59e-69 | 0.887 .. 1.15 | 24.9 .. 44.2 | 9.26 (0.13) | True |
| ζ | odd | 6 | 2..16 | 1.73e+1 .. 1.96e-65 | 0.953 .. 1.05 | 28.0 .. 52.1 | 11.5 (0.47) | True |
| χ₋₃ | even | 1 | 2..16 | 8.48e-4 .. 1.34e-75 | 0.957 .. 1.08 | 18.1 .. 28.7 | 5.1 (0.068) | True |
| χ₋₃ | odd | 3 | 2..16 | 4.11e-2 .. 2.06e-71 | 0.957 .. 1.14 | 21.9 .. 38.3 | 7.84 (0.28) | True |
| χ₋₄ | even | 1 | 2..16 | 6.79e-4 .. 1.0e-75 | 0.958 .. 1.08 | 17.8 .. 28.4 | 5.07 (0.07) | True |
| χ₋₄ | odd | 3 | 2..16 | 4.81e-2 .. 2.21e-71 | 0.958 .. 1.14 | 22.1 .. 38.4 | 7.8 (0.28) | True |
| χ₅ | even | 0 | 2..16 | 3.41e-5 .. 1.89e-78 | 1.09 .. 1.26 | 14.8 .. 22.1 | 3.48 (0.046) | True |
| χ₅ | odd | 2 | 2..16 | 2.0e-3 .. 2.53e-73 | 1.01 .. 1.13 | 18.9 .. 33.9 | 7.21 (0.029) | True |
| χ₋₇ | even | 1 | 2..16 | 4.73e-4 .. 6.41e-76 | 0.96 .. 1.08 | 17.5 .. 27.9 | 5.03 (0.073) | True |
| χ₋₇ | odd | 3 | 2..16 | 5.72e-2 .. 2.37e-71 | 0.96 .. 1.14 | 22.3 .. 38.4 | 7.75 (0.27) | True |
| χ₈ | even | 0 | 2..16 | 2.39e-5 .. 1.23e-78 | 1.1 .. 1.26 | 14.5 .. 21.7 | 3.45 (0.049) | True |
| χ₈ | odd | 2 | 2..16 | 2.23e-3 .. 2.54e-73 | 1.01 .. 1.13 | 19.0 .. 33.9 | 7.16 (0.04) | True |
| χ₋₂₀ | even | 1 | 2..16 | 7.74e-4 .. 9.47e-76 | 0.962 .. 1.09 | 18.0 .. 28.3 | 4.98 (0.075) | True |
| χ₋₂₀ | odd | 3 | 2..16 | 1.58e-1 .. 6.28e-71 | 0.962 .. 1.14 | 23.3 .. 39.4 | 7.73 (0.27) | True |

**Table V** (`kb_cert_supports_m2_var.json`, `kb_cert_supports_m3.json`): certified KB bounds at the supports over (m, β = c − Δβ) (all [C]; Theorem C's table uses m = 2, Δβ ∈ {0, 2})

| support | sector | m = 2: Δβ = 0 | 1 | 2 | 3 | 4 | m = 3: Δβ = 0 | 2 | 4 | best / table value |
|---|---|---|---|---|---|---|---|---|---|---|
| 1.6 | even | 4.96e-12 | 3.45e-13 | 2.94e-13 | 3.21e-12 | 7.09e-12 | 4.15e-11 | 1.16e-12 | 1.34e-11 | 1.0 |
| 1.6 | odd | 7.83e-10 | 8.51e-11 | 1.04e-10 | 1.25e-9 | 2.42e-9 | 4.36e-9 | 2.48e-10 | 2.55e-9 | 0.815 |
| 2.38 | even | 1.67e-39 | 3.57e-42 | 4.47e-42 | 8.69e-41 | 6.31e-41 | 5.75e-38 | 2.37e-41 | 8.9e-40 | 0.8 |
| 2.38 | odd | 1.38e-36 | 4.78e-39 | 6.15e-39 | 8.75e-38 | 2.24e-37 | 3.11e-35 | 2.9e-38 | 5.48e-37 | 0.778 |
| 2.6 | even | 9.54e-53 | 8.91e-56 | 9.4e-56 | 2.25e-54 | 1.37e-54 | 4.97e-51 | 7.35e-55 | 3.23e-53 | 0.947 |
| 2.6 | odd | 1.24e-49 | 1.87e-52 | 1.89e-52 | 3.06e-51 | 8.43e-51 | 4.25e-48 | 1.42e-51 | 2.54e-50 | 0.99 |
| 2.99 | even | 1.5e-85 | 2.63e-89 | 3.45e-89 | 9.0e-88 | 4.12e-88 | 1.67e-83 | 3.18e-88 | 2.18e-86 | 0.762 |
| 2.99 | odd | 4.36e-82 | 1.22e-85 | 1.26e-85 | 2.11e-84 | 6.22e-84 | 3.17e-80 | 1.36e-84 | 2.86e-83 | 0.974 |

The full balls (A, B, ‖f_win‖², T_c, T_L, the zero integral and B_cert) are in the JSON files. Every entry of Table V is also a certified bound, so Theorem C holds with the minimum over Table V's columns: 0.76–1.0 × the Table S values. R2's parameters (m = 2, β ∈ {c, c − 2}) are close to the best of this family; β = c − 1 is the best tested at 7 of the 8 support cases.

**Reading.**

- **Against R2's floating values.**
  - Hermite: the certified bounds are 1.006–1.103 × R2's. That comes from the termwise bound for ∫|S/2 + wS′| in A (R2 integrated |f′| numerically) and from Bellotti–Wong's slightly larger zero count.
  - KB: 0.89–1.26 × R2's. The certified edge value is the exact signed sum, against R2's sum of envelopes, which is lower. The harmonic-sum Lemma H and the termwise envelopes are higher.
  - So R2's tabulated values were not all proven bounds. Where B_cert/R2 > 1, R2's number is not established by anything here. All of R2's qualitative statements (rates, measured/bound ≤ 1) survive with B_cert in place of R2's numbers.
- **Rates.**
  - Hermite: ln B_cert + c stays in a band of width ≤ 14 over x/q = 2..16 for every object. The rate e^{−c} is reproduced pointwise.
  - KB: ln B_cert + 2c grows by 7–24 over the grid while 2c grows by 176. That is the polynomial factor c^p of R2's Table 3.5, now certified pointwise.

## 8. Validation: the checks that can fail

The rule: an upper bound is never compared with another upper bound. B_cert, the trial's Rayleigh quotient RQ and our Ritz minima are all upper bounds for λ_s(x), so B_cert/Ritz validates nothing. The checks below compare an upper bound with a lower bound, or a quantity with the bound that a theorem says must exceed it.

**(a) RQ(trial) ≤ B_cert.**
- **The theorem behind it.** Theorems 1–2 give Q(f_win)/‖f_win‖² ≤ B_cert for the very trial function f_win. That is a statement that can fail: a sign error in σ, a wrong Gauss sum or an error in the transforms would make Q(f_win) O(‖f_win‖²) instead of e^{−2c}-small.
- **The computation** (`rq_check.py`) [N].
  - f_win is projected on our N-mode basis of the sector, with N = min(scan's N, 240), and N = 60, 100, 120, 180 at the supports.
  - RQ = cᵀMc/cᵀc is computed with M = `decay_law_mp.zeros_side` (archimedean + primes + pole; no zeros enter).
  - Everything is redone at dps and dps + 20, with dps = ⌈2c/ln 10⌉ + 40; agreement < 10⁻⁶ is required ("precision-stable").
  - **No eigensolver is involved** (RQ is a quadratic form of a fixed vector), so the inv+full protocol of the scans is replaced by this two-precision check.
- **The KB window function.** It contains out-of-band terms of relative size e^{−c}, which carry Q's whole value. They are kept for n ≤ k + n_out on each band interval, with n_out = 10; convergence in n_out is reported below. They are projected on a fine grid at 96 bits. RQ differs from Q(f_win)/‖f_win‖² by the basis truncation and by n_out; neither is certified.

**Result: 386/386, all precision-stable** (`rq_check_grid.json`, `rq_check_supports.json`).

| construction | points | RQ ≤ B_cert | B_cert/RQ (min .. max) |
|---|---|---|---|
| Kaiser–Bessel | 84 grid + 8 supports | 92/92 | 59.56 .. 3185.0 |
| Hermite | 286 grid + 8 supports | 294/294 | 13.33 .. 100.9 |

- **Every case passes.** The tightest is ζ Hermite at x = 2, even sector, with B_cert/RQ = 13.3.
- **Convergence in n_out.** RQ was recomputed with n_out = 20 at 20 KB cases: ζ, χ₅ and χ₋₂₀ at x/q = 2 and 16, both sectors, plus the 8 support cases. The relative change is at most 0.0017 (χ₋₂₀, x/q = 2, even) and has median 3e-06. That is far inside every margin.
- **Basis truncation.** The captured fraction |c|²/‖f_win‖² (with ‖f_win‖² approximated by its in-band part) is within 10⁻⁶ of 1 at 353 of the 386 cases. Its extremes are 0.9914 (ζ Hermite, x = 2, odd, N = 24) and 1.0025 (χ₋₇ KB, x/q = 2, odd), all at x/q ≤ 3. Every case passes with B_cert/RQ ≥ 13. The `captured` field in the JSON records it per point.

**(b) B_cert ≥ certified lower bound at the four certified supports.** This is an upper bound against a lower bound for the same λ_s(x).

**Result: 16/16.** Table S, last column, and `summary.json`. The smallest ratio B_cert/cert is 1.15·10⁴ (KB, support 1.6, odd sector); for Hermite it is 1.3·10⁹ or more. No violation.

**(c) RQ(trial) ≥ certified lower bound at the four supports.** An upper bound against a lower bound; it would fail if the projection or the matrix were wrong in the other direction.

**Result: 16/16.** RQ/cert lies in [30.2, 4.49e+04] for KB and in [9.2e+07, 1.33e+49] for Hermite (Table D). No violation.

**Table D** where the gap B_cert/λ comes from (supports; diagnostic only, not a validation)

| support | sector | construction | B_cert/RQ (Thm 2 slack) | RQ/Ritz ≤ RQ/λ ≤ RQ/cert (trial excess) |
|---|---|---|---|---|
| 1.6 | even | kb | 5.47e+2 | 3.29e+1 .. 5.23e+1 |
| 1.6 | even | hermite | 1.86e+1 | 5.46e+8 .. 8.69e+8 |
| 1.6 | odd | kb | 3.79e+2 | 1.78e+1 .. 3.02e+1 |
| 1.6 | odd | hermite | 1.43e+1 | 5.42e+7 .. 9.2e+7 |
| 2.38 | even | kb | 2.53e+2 | 1.74e+3 .. 2.59e+3 |
| 2.38 | even | hermite | 2.59e+1 | 2.05e+24 .. 3.05e+24 |
| 2.38 | odd | kb | 2.74e+2 | 3.98e+2 .. 5.73e+2 |
| 2.38 | odd | hermite | 1.58e+1 | 2.43e+23 .. 3.49e+23 |
| 2.6 | even | kb | 3.36e+2 | 3.07e+3 .. 4.86e+3 |
| 2.6 | even | hermite | 2.24e+1 | 3.78e+31 .. 5.99e+31 |
| 2.6 | odd | kb | 3.33e+2 | 6.74e+2 .. 1.06e+3 |
| 2.6 | odd | hermite | 2.11e+1 | 2.78e+30 .. 4.37e+30 |
| 2.99 | even | kb | 2.19e+2 | 2.51e+4 .. 4.49e+4 |
| 2.99 | even | hermite | 2.64e+1 | 7.45e+48 .. 1.33e+49 |
| 2.99 | odd | kb | 2.7e+2 | 3.15e+3 .. 5.65e+3 |
| 2.99 | odd | hermite | 2.32e+1 | 5.17e+47 .. 9.27e+47 |

**(d) Structural checks** [C/N]:
- Theorem 1(i) as a ball: E(1/w₀) − σE(w₀) encloses 0 at every certified point (286 Hermite, 84 KB, 16 supports).
- ψ(0) encloses 0 for every ζ KB trial (12 grid points + 8 supports); the Hermite ζ trials have ψ(0) = 0 exactly (integer coefficients).
- The root-free certificate of Lemma G held for every Hermite out-of-band integral.
- The zero-count checks of §2.3 hold.

**Where the gap comes from** (diagnostic, not a validation).
- B_cert/λ = (B_cert/RQ)·(RQ/λ). The first factor is Theorem 2's slack on that trial; the second is the trial's own excess.
- At the supports, RQ/Ritz ≤ RQ/λ ≤ RQ/cert brackets the excess (Table D).
- **Kaiser–Bessel:** the slack stays within 59.56–3185.0 over all 92 points, with no trend in c at the supports (2.2·10²–5.5·10²). The excess grows with c: RQ/Ritz runs from ≈ 1 at x/q = 2 to 9.235e+4 on the grid, and from 18–52 at x = 4.95 to 2.5·10⁴–4.5·10⁴ (even) at x = 19.9.
- **Hermite:** the slack is 14–26 at the supports, and the excess is the e^{c} of its rate.
- **So the KB bound's distance from the truth at large x is mostly the trial's** (the c^p factor), not Theorem 2's. That is what makes §9 worth asking.

## 9. The sharp prolate trial (Step 4)

**9.1 Why the sharp trial is outside Theorem 2** [D]. Let φ = ξ_n on [−λ, λ], the time-limited prolate, with φ(±λ) ≠ 0.
1. **Decay of ψ.** Integrating by parts once, φ̂(ξ) = −[φ(λ)e^{−2πiλξ} − φ(−λ)e^{2πiλξ}]/(2πiξ) + O(ξ^{−2}). So ψ ~ C sin(2πλy + ϕ)/y out of band, and (A2) fails (δ = 0).
2. **B₀ diverges.** For w > √x, S(w) ≈ (C√q/w) Σ_n χ(n) sin(nθ + ϕ)/n with θ = 2πλw/√q. The n-sum is a bounded, nonvanishing sawtooth-type function of θ, so |S(w)| ≍ 1/w on a set of positive density, and B₀ = 2∫^∞|S(w)|dw = ∞ (logarithmically).
3. **A₀ diverges too.** f_out has jumps where θ ∈ 2πℤ.

Theorem 2 therefore gives nothing for the sharp trial, as R2 said.

**9.2 Proposition P (a C¹ edge is enough)** [D].
- **Hypotheses.** φ is supported in [−λ, λ], of parity κ, C¹ on ℝ (so φ(±λ) = φ′(±λ) = 0), and φ″ is of bounded variation, with total variation V including the jumps at ±λ.
- **Conclusion 1.** |φ̂(ξ)| ≤ V/(2π|ξ|)³ and |φ̂′(ξ)| ≤ 2πV₁/(2π|ξ|)³, with V₁ the variation of (yφ)″.
- **Conclusion 2.** ψ = (φ + μ̄φ̂)/2 satisfies (A1)–(A3) with δ = 2 (and (A4) for ζ via a two-term combination). f_ψ is C¹ on [L/2, ∞), A₀, B₀ < ∞, and Theorem 2 applies.

*Proof.* Integrate by parts twice; the boundary terms vanish because φ, φ′ vanish at ±λ. Then ∫φ″e^{−2πiyξ}dy = (2πiξ)⁻¹∫e^{−2πiyξ}dφ″(y), which is a Stieltjes integral bounded by V. The same holds for yφ. The decay O(y^{−3}) of ψ and ψ′ gives absolute convergence of S, S′ and of the integrals defining A₀, B₀. ∎

So a smoothed edge does give finite A and B. The Kaiser–Bessel trial (m = 2) is one such φ, with closed-form transform.

**9.3 The rate** [heuristic sketch; not claimed].
- **The taper.** Taper the prolate on its own edge scale: τ = S((1−|t|)/ε) with S the C² smoothstep and ε = κ/c². Near t = 1 the prolate behaves like ξ(1)·I₀(c√(2(1−t))), which follows from the regular singular point of the prolate equation at t = 1, so ξ changes by a factor O(1) over 1 − t ~ 1/c².
- **Size of the out-of-band part.** Its transform is ≈ e₀/ω, with e₀ := |ξ(1)| ≍ √(c(1−λ_n(c))) relative to the interior ([S]: Slepian's identity dλ_n/dc = 2λ_n ψ_n(c,1)²/c, not re-read), up to the taper's own bandwidth ω ≈ c²/κ, and O(ω⁻⁴) beyond.
- **A and B.** Hence B ≈ e₀·log(c/κ) and A ≈ e₀·c²/κ (up to factors of λ and q).
- **The bound.** It is ∝ A·B·log(A/B)/‖f_win‖² ∝ (1 − λ_n(c))·poly(c), i.e. heuristically C c^p e^{−2c}, at a polynomial cost that grows like c²/κ from A. Like Conjecture R, this is not proved: it would need the analogues of M1 and M2.
- **The trade-off in κ.** A wider taper (larger κ) lowers A but lets the taper reach into the region where ξ is larger (ξ(1 − κ/c²) ≈ e₀I₀(√(2κ))), which raises the leakage. §9.4 measures it.

**9.4 Evaluation** [N, floating point; `prolate_taper.py` → `data/research_r2_rigorous/prolate_taper.json`].
- **What is computed.** ζ at supports 1.6 and 2.38 (both sectors, taper widths κ ∈ {2, 8, 32, 128}) and at supports 2.6 and 2.99 (both sectors, κ = 8, the best width at the first two). For each: the trial's RQ in our basis (as in (a)), and a floating estimate of Theorem 2's bound.
- **How the estimate is made.**
  - Φ(ω) is evaluated as the Legendre–Bessel series of the prolate minus the edge-taper transform, at dps ≈ c/ln 10 + 40.
  - B and A are trapezoid integrals of |Φ| and ω|Φ′| on [c, 20/ε + 10c] (8 samples per 2π), plus an extrapolated ω⁻⁴ tail.
  - Then the edge sum, the norm, and the certified zero integral.

**Table P** (`prolate_taper.json`; [N] except the KB column)

| support | sector | κ (ε = κ/c²) | taper RQ | RQ/Ritz | taper Thm-2 estimate [N] | KB B_cert [C] | estimate/KB | taper slack (estimate/RQ) | KB slack (B_cert/RQ) | T_c taper |
|---|---|---|---|---|---|---|---|---|---|---|
| 1.6 | even | 2 | 2.52e-17 | 1.54 | 1.25e-11 | 2.94e-13 | 42.7 | 4.98e+5 | 5.47e+2 | 1280.8 |
| 1.6 | even | 8 | 3.11e-17 | 1.9 | 4.58e-12 | 2.94e-13 | 15.6 | 1.47e+5 | 5.47e+2 | 525.97 |
| 1.6 | even | 32 | 9.4e-16 | 57.5 | 2.26e-11 | 2.94e-13 | 77.0 | 2.41e+4 | 5.47e+2 | 224.21 |
| 1.6 | even | 128 | 8.63e-12 | 5.28e+5 | 1.94e-8 | 2.94e-13 | 6.61e+4 | 2.25e+3 | 5.47e+2 | 93.278 |
| 1.6 | odd | 2 | 2.15e-14 | 1.39 | 1.15e-8 | 1.04e-10 | 110.0 | 5.34e+5 | 3.79e+2 | 1250.6 |
| 1.6 | odd | 8 | 2.26e-14 | 1.46 | 3.45e-9 | 1.04e-10 | 33.0 | 1.53e+5 | 3.79e+2 | 493.33 |
| 1.6 | odd | 32 | 4.15e-13 | 26.8 | 9.68e-9 | 1.04e-10 | 92.6 | 2.33e+4 | 3.79e+2 | 210.65 |
| 1.6 | odd | 128 | 1.72e-9 | 1.11e+5 | 2.42e-6 | 1.04e-10 | 2.32e+4 | 1.41e+3 | 3.79e+2 | 87.139 |
| 2.38 | even | 2 | 1.64e-47 | 1.61 | 8.82e-41 | 4.47e-42 | 19.7 | 5.39e+6 | 2.53e+2 | 5650.8 |
| 2.38 | even | 8 | 1.66e-47 | 1.64 | 4.61e-41 | 4.47e-42 | 10.3 | 2.77e+6 | 2.53e+2 | 2483.9 |
| 2.38 | even | 32 | 3.78e-46 | 37.3 | 5.84e-40 | 4.47e-42 | 131.0 | 1.55e+6 | 2.53e+2 | 1143.4 |
| 2.38 | even | 128 | 7.67e-41 | 7.56e+6 | 6.48e-36 | 4.47e-42 | 1.45e+6 | 8.45e+4 | 2.53e+2 | 524.45 |
| 2.38 | odd | 2 | 9.54e-44 | 1.69 | 5.11e-37 | 6.15e-39 | 83.0 | 5.35e+6 | 2.74e+2 | 5599.0 |
| 2.38 | odd | 8 | 9.72e-44 | 1.72 | 2.45e-37 | 6.15e-39 | 39.8 | 2.52e+6 | 2.74e+2 | 2419.3 |
| 2.38 | odd | 32 | 1.13e-42 | 20.1 | 2.43e-36 | 6.15e-39 | 395.0 | 2.15e+6 | 2.74e+2 | 1112.4 |
| 2.38 | odd | 128 | 2.13e-37 | 3.77e+6 | 1.64e-32 | 6.15e-39 | 2.67e+6 | 7.74e+4 | 2.74e+2 | 506.06 |
| 2.6 | even | 8 | 1.59e-61 | 1.74 | 8.56e-55 | 9.4e-56 | 9.1 | 5.39e+6 | 3.36e+2 | 3861.6 |
| 2.6 | odd | 8 | 1.46e-57 | 1.73 | 7.5e-51 | 1.89e-52 | 39.6 | 5.13e+6 | 3.33e+2 | 3786.5 |
| 2.99 | even | 8 | 2.37e-95 | 3.78 | 1.77e-88 | 3.45e-89 | 5.14 | 7.48e+6 | 2.19e+2 | 8297.9 |
| 2.99 | odd | 8 | 2.07e-91 | 1.4 | 3.7e-84 | 1.26e-85 | 29.5 | 1.79e+7 | 2.7e+2 | 8153.4 |

**Reading.**
1. **Smoothing works for the trial.** The tapered trial's RQ is 1.4–1.9 × our Ritz minimum at κ ≤ 8 (3.8 at 2.99 even). It is 10–7·10³ × better than the Kaiser–Bessel trial's RQ (6.6·10³ at 2.99 even): Table D gives KB's RQ/Ritz as 18–33 at x = 4.95 and up to 2.5·10⁴ at x = 19.9. The C² taper on the edge scale 1/c² costs almost nothing in the trial's quality, which confirms R2's near-optimality of the prolate.
2. **It does not work for Theorem 2.** The taper's Theorem-2 slack (estimate/RQ) is 10⁵–10⁷, against Kaiser–Bessel's 2–5·10². The reason is A: the out-of-band mass reaches ω ≈ c²/κ, so A/B = T_c ≈ 500–8000. The zero sum then charges B² to every zero below T_c, with the worst-case strip weight e^{u/2}.
3. **Net effect: the tapered prolate bound is worse at every tested point.**
   - Best κ = 8 in every case.
   - Even sector: 15.6×, 10.3×, 9.1×, 5.1× at supports 1.6, 2.38, 2.6, 2.99.
   - Odd sector: 33×, 40×, 40×, 30×.
   - κ ≥ 32 makes the trial itself much worse (the taper reaches the region where ξ grows like I₀(√(2κ))).
4. **The even-sector trend.** The gap falls in the even sector, because KB's trial excess grows like c^p (Table D) while the taper's slack grows more slowly. A crossover at larger x is plausible in the even sector but **not established**: four floating points, no trend in the odd sector.

**Answer to Step 4.**
- **Yes:** a smoothed edge gives a finite B, and the rate e^{−2c} is retained in the sense of the heuristic §9.3 (not proved, like Conjecture R).
- **No:** it does not sharpen the certified bound at the tested supports. Theorem 2 rewards a compact out-of-band spectrum (small A/B) more than a small L² leakage, and Kaiser–Bessel's analytic edge has exactly that.
- **What the data point to instead.** Table D shows that at large x the Kaiser–Bessel bound is limited by its trial (excess 10³–4·10⁴), not by Theorem 2 (slack ≈ 250). The sharpening worth pursuing is either:
  - a closed-form C¹ trial closer to the prolate in L² that keeps the KB-type compact out-of-band spectrum; or
  - a version of Theorem 2 that charges zeros on the line (RH verified below a height) without the e^{u/2} weight (§9.5, "using verified zeros").

**9.5 What a certified tapered-prolate bound would need** (specified, not run).
- **The obstacle.** Rigorous enclosures of the out-of-band transform Φ on [c, ≈ c²/κ]. Φ is an e^{−c}-small difference of O(1/ω) terms (Legendre–Bessel series of degree ≈ 2c), so ball evaluation on ω-intervals loses everything to dependency.
- **One route.** Exact-point Arb evaluations at a grid with spacing h, plus a Taylor model of order J in ω at each grid point (Φ^{(j)} = transform of (−it)^jφ, again a Legendre–Bessel series). The remainder is bounded trivially by (h/2)^J‖t^Jφ‖₁/J!. J ≈ 45 (c = 125) makes the remainder < e^{−c}; the cost is ≈ J·2c Bessel terms per grid point over ≈ c²/(κh) points. That is 10⁸–10⁹ ball operations per (support, sector), hours on one core.
- **Should it be run?** Only if §9.4 showed a gain over the certified KB bound. It does not at any tested support. If the even-sector trend continues past x ≈ 20, it could become worthwhile at larger x.

**Other routes.**
- **Using verified zeros (remark, not evaluated)** [D]. A different and cheaper sharpening applies to every trial. For zeros below a height H at which RH is verified, |F_out(γ)| ≤ B_line := 2∫_{L/2}^∞|f_out| du, with no e^{u/2} weight, and A_line likewise. Only zeros above H need Theorem 2's strip bound, and those contribute ≤ A²·O(log H/H). The available heights are:
  - for ζ: H = 3·10¹² [S: Platt–Trudgian 2021];
  - for q < 935: H = 2(e⁶π − q)/q [R: BMOR Lemma 6.1(a)];
  - larger heights from Platt's GRH verification [S: Platt, Math. Comp. 85 (2016), cited by BMOR; the height figures were not checked].

  The weight e^{u/2} ≈ √x at the window edge suggests a gain of order x in B·A. That brings in the zero-verification computations as trust assumptions, which is why it is not in the certified chain here.
- **Better closed-form C¹ trials.** Two candidates:
  - KB with other (m, β);
  - the "corrected Kaiser" φ = I₀(z) − 1 − z²/4, z = β√(1−t²), which has a closed-form transform but a larger edge coefficient (β⁴ against β²).

  Neither was found to beat m = 2 in this pass; that is not established either way.

## 10. What is proved, what is assumed, what is a sketch

| item | status |
|---|---|
| Theorem 1 (E-map identity), Theorem 2 (zero-sum bound), Lemmas Z, G, K, E, H, T, Proposition P | [D], proofs above; unconditional |
| explicit formula Q(f) = Σ_ρ F(γ_ρ)F(−γ_ρ) for compactly supported BV f, and our matrices computing its arithmetic side | assumed [S] (as for the certificates) |
| zero counts: Bellotti–Wong Thm 1.1 (ζ; for T ≤ 3.06·10¹⁰ it rests on Platt's computed bound 2.5167 for the absolute value of S(T) [S]), BMOR Thm 1.1 (L(s, χ)) | [R], literature theorems |
| DLMF 10.14.4, 10.49.1–2 | [R] |
| B_cert at the listed x (Theorem C: 286 + 84 grid points, 16 support cases, and Table V) | [C], finite statements only; given Arb's correctness and this code's (`cert_lib.py`, ≈ 600 lines, unaudited) |
| Hermite rate: λ_s(x) ≤ K_H c^{d+1}e^{−c} for all x ≥ 2q (Proposition H) | [D + C]: proved, with certified constants |
| Kaiser–Bessel rate e^{−2c} for all x (Conjecture R) | **conjecture**, supported by Theorem C's table; missing estimates M1 (uniform norm lower bound) and M2 (uniform polynomial bounds on A, B) |
| RQ values, validation (a), prolate-taper numbers | [N] |
| multiplicity / central zeros | counted with multiplicity throughout; nothing assumes simple zeros or m₀ ≤ 1 |
| GL(2), degree d | **out of scope**: see the remark below the table |
| what an upper bound says about RH | nothing: every bound here holds for functions with off-line zeros too |

**Remark (out of scope: GL(2) and higher degree).** Nothing in this sheet is claimed beyond degree 1.
- R2 §6.1's Voronoi/Hankel analogue of Theorem 1 is a sketch.
- No explicit zero count for GL(2) L-functions was sourced, and no GL(2) bound is evaluated.
- P-R389's central-zero rule concerns the true minimum, which an upper bound cannot address. An E-map trial would vanish at the central zero to at least m₀ automatically, as it does at every zero.

## 11. Reproduction

```bash
P=.venv/bin/python; S=scripts/research_r2_rigorous; O=data/research_r2_rigorous
$P $S/zero_count_check.py --json $O/zero_count_check.json
$P $S/certify_bounds.py --construction hermite --d 1,5,8,-3,-4,-7,-20,-8,-11,12,13,-15,17 --xq 2,3,4,5,6,7,8,10,12,14,16 --json $O/hermite_cert.json   # 6 s
$P $S/certify_bounds.py --construction hermite --zeta-supports --json $O/hermite_cert_supports.json
$P $S/certify_bounds.py --construction kb --zeta-supports --workers 3 --json $O/kb_cert_supports.json            # 1 min
$P $S/certify_bounds.py --construction kb --d 1,5,8,-3,-4,-7,-20 --xq 2,4,6,8,10,16 --workers 3 --json $O/kb_cert.json   # 12 min
$P $S/rq_check.py --supports --workers 3 --json $O/rq_check_supports.json
$P $S/rq_check.py --grid --workers 3 --json $O/rq_check_grid.json                                                 # ~1 h (resumable: --resume)
$P $S/rq_check.py --grid --herm-extra --resume --workers 3 --json $O/rq_check_grid.json                          # the other 202 Hermite points
$P $S/rq_check.py --grid --kb-only --only 1,5,-20 --xq 2,16 --n-out 20 --workers 3 --json $O/rq_check_nout20_grid.json
$P $S/rq_check.py --supports --kb-only --n-out 20 --workers 3 --json $O/rq_check_nout20_supports.json
$P $S/hermite_rate.py --json $O/hermite_rate.json                                                                # Proposition H constants
$P $S/certify_bounds.py --construction kb --zeta-supports --dbeta 0,1,2,3,4 --m 2 --workers 3 --json $O/kb_cert_supports_m2_var.json
$P $S/certify_bounds.py --construction kb --zeta-supports --dbeta 0,2,4 --m 3 --workers 3 --json $O/kb_cert_supports_m3.json
$P $S/prolate_taper.py --supports 1.6,2.38 --kappa 2,8,32,128 --json $O/prolate_taper.json                       # [N], ~1 h
$P $S/prolate_taper.py --supports 2.6,2.99 --kappa 8 --resume --json $O/prolate_taper.json
$P $S/summarize.py
```

## 12. Where an auditor should look hardest

1. **The KB norm lower bound** (`cert_lib.kb_quantities`, §6 "the norm").
   - It is the only place where a non-analytic object (the infinite out-of-band sum, a Fourier-type series in w) meets `acb.integral`.
   - The split S = S_main + S_rest per band interval must put each n on the right side of the band edge. The exact dyadic subintervals must lie inside I_k (`up(√x/(k+1))`, `lo(√x/k)`, the w = 1 assertion).
   - ε_k must bound the remainder for every w in the subinterval (Lemma T with ω_n ≥ c n a_k/√x). An error here would inflate ‖f_win‖² and make B_cert too small. (a) catches that only once it exceeds the margin B_cert/RQ (≥ 13 everywhere, 60–3·10³ for KB), and (b) only beyond B_cert/cert ≥ 1.1·10⁴, so smaller errors would go unnoticed.
2. **The envelope chain for A and B** (`env_G`, `C_E`, Lemmas E, H, T).
   - Check that cap and hank are each valid on all of [c, ∞), so that the floating split point is harmless.
   - Check that every (p, k) term has p − k − 1 − j ≤ −2 where a closed-form tail is used (asserted).
   - Check the constant bookkeeping B = (√q β^m/√(2π))J₀ and A = 2(√x|S(√x)| + B/4 + (√q β^m/(2√(2π)))J₁) against the substitution y = nw/√q, ω = 2πλy.
3. **The zero counts and their use.**
   - Bellotti–Wong's Theorem 1.1, second estimate, is from arXiv v2. At the heights used (T_L = 14–700) it rests on Platt's rigorous |S(T)| ≤ 2.5167 [S], through their eq. (5.1).
   - The BMOR inequality must be applied only above t_low.
   - Lemma Z's lower limit must be the lower endpoint of A/B, with N = 0 below t_low.
   - The fallback with a refereed source (HSW 2022, with C₃ corrected to 9.4925) would raise ζ bounds by ≤ 1.66×; it does not change any verdict in §8.
