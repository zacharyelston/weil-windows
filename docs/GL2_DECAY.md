# Weil's Form for a Degree-2 L-Function: P-DL4

This note executes P-DL4, pre-registered at the end of `docs/DECAY_LAW.md`. It extends rh2's zeros-side Weil form from degree 1 to L(E, s) for an elliptic curve E/ℚ, and scans the window minimum λ(x).

Labels used throughout:
- *validated*: checked against an independent computation.
- *computed*: the output of the validated code, with no independent check beyond precision and convergence.
- *could not verify*: stated, but not checked here.

Code: `scripts/gl2_form_mp.py`. It only imports `connes_letter_mp.py` and `decay_law_mp.py`; neither was modified. Data: `data/connes/gl2/`.

## 1. The archimedean block for Γ_ℝ(s + μ)

Notation follows `connes_letter_mp.build_form`:
- the window is [−L/2, L/2], with L = log x;
- G(w) = ∫F(v + w)F(v)dv is supported in [−L, L], and h = |F̂|².

For one factor Γ_ℝ(s + μ), with μ ≥ 0 real, put σ = ¼ + μ/2. In the representation ψ(z) = −γ + ∫₀^∞(e^{−t} − e^{−zt})/(1 − e^{−t})dt, substitute t = 2w. Weil's archimedean term then becomes

  (1/2π)∫h(t)[Re ψ(σ + it/2) − log π]dt = [ψ(σ) − log π]G(0) + ∫₀^L K_σ(w)(G(0) − G(w))dw + G(0)∫_L^∞ K_σ,

where K_σ(w) = 2Σ_m e^{−a_m w} and a_m = 2m + ½ + μ.

This is the μ = 0 kernel 2K(w) = 2e^{w/2}/(e^w − e^{−w}), with every exponent shifted by μ. rh2 subtracts 2L e^{−w/2} inside J_k, and that subtraction is kept. Using ∫₀^∞ K_σ(1 − e^{−w/2})dw = ψ(σ + ¼) − ψ(σ), the ω-independent part collects into

  diag_const(μ) = log π − ψ((1 + μ)/2) − Σ_m 2e^{−(2m+1+μ)L}/(2m + 1 + μ).

- **μ = 0:** this is log 4π + γ + log tanh(L/2), rh2's constant. The tail sums to 2 artanh(e^{−L}) = −log tanh(L/2).
- **μ = 1:** it is γ + log π + log(1 − e^{−2L}).

The closed forms carry over with z_μ = (½ + μ)/2 + iω/2:
- I_k = ½ Im ψ(z_μ) − Σ_m ω e^{−aL}/(a² + ω²);
- J_k = L(ψ((1 + μ)/2) − Re ψ(z_μ)) − ½ Re ψ′(z_μ) + Σ_m [exact − E-free].

The bracketed tail decays like e^{−2mL}, as in rh2.

**The L(E, s) form.** Use the analytic normalisation L(E, s) = Σ a_n n^{−1/2} n^{−s}. Then Λ(s) = N^{s/2}Γ_ℂ(s + ½)L(E, s) = εΛ(1 − s). By duplication, Γ_ℂ(s + ½) = 2(2π)^{−s−½}Γ(s + ½) = Γ_ℝ(s + ½)Γ_ℝ(s + 3/2) exactly. So:
- the zeros-side form is block(½) + block(3/2) + log N·I − (prime terms);
- the prime terms are (log p^k, (α_p^k + β_p^k) log p/p^{k/2}) for p^k ≤ x, in rh2's convention;
- at good p, α_p + β_p = a_p/√p and α_pβ_p = 1;
- at p | N, α_p = a_p/√p and β_p = 0;
- there is no pole.

The conductor enters as log N·G(0), the same convention as log q for characters. The check is the zero sum in §4: an error in this constant would show up there as an O(1) discrepancy.

### Hard validations 1 and 2: pass (*validated*)

`gl2_form_mp.py validate` writes `validate_arch.json`. It compares maximum entrywise differences at dps 100, with N_basis = 30, x ∈ {13, 22, 132, 444}, in both sectors:

| check | max \|Δ\| | max entry |
|---|---|---|
| 1. block(0) against `build_form(x, n, "zeta", terms=[])` | 7e-102 to 6e-101 | 2–4 |
| 2. block(1) + log 5·I against `build_form(x, n, "dh", terms=[])` | 6e-101 to 1.1e-100 | 2–4 |
| block(0) with Λ(n)/√n against `build_form(zeta)` | 1.1e-100 to 3e-100 | |
| block(0) + log 5 + χ₅ primes against `zeros_side(5)` | 8e-101 to 4.3e-100 | |
| block(1) + log 4 + χ₋₄ primes against `zeros_side(−4)` | 1.7e-100 to 3.7e-100 | |

At dps 50 the same differences are 1e-50 to 2e-50. They track the working precision, so the agreement is exact up to rounding. The last three rows also check the matrix assembly (`assemble`), which re-implements `build_form`'s q(j, k) and is used for every form here.

## 2. Curves (recorded before any scan)

**E1 = 11a1** and **E2 = 37b1**. The registration names 37b1 if it is rank 0 with root number +1, and it is (below), so no choice was left.

| | 11a1 | 37b1 |
|---|---|---|
| model [a1, a2, a3, a4, a6] | [0, −1, 1, −10, −20] | [0, 1, 1, −23, −50] |
| Δ | −161051 = −11⁵ | 50653 = 37³ |
| c4 (prime to N) | 496 | 1120 |
| j | −122023936/161051 = −757.6726 | 1120³/37³ = 27736.32 |
| reduction at N | split multiplicative (a₁₁ = 1) | split multiplicative (a₃₇ = 1) |
| conductor | 11 | 37 |
| root number ε (AFE) | +1 | +1 |
| L(E, 1) (arithmetic) = L(E, ½) (analytic) | 0.25384186085591068434 | 0.72568106193615278234 |
| a_p, p = 2, 3, 5, 7, 11, 13 | −2, −1, 1, −2, 1, 4 | 0, 1, 0, −1, 3, −4 |

**Conductor.** Δ = ±N^k with k < 12, so the model is minimal. c4 is prime to N, so the reduction at N is multiplicative, with conductor exponent 1. Δ is prime to every other p, so the reduction is good there. Hence the conductor is N. (*validated* by the Tate criteria, and again by the functional-equation check below, which fails for any other N.)

**a_p by point counting.** a_p = p + 1 − #Ẽ(F_p), counting every projective point of the reduced cubic; the formula holds at the bad prime too. The a_n follow from the Hecke recursion at good p and a_{p^k} = a_p^k at p | N. Checks (*validated*):
- **Known values.** 11a1 gives a₂ = −2, a₃ = −1, a₅ = 1, a₇ = −2, as expected.
- **Eta product.** For 11a1, every a_n with n ≤ 2000 equals the q-expansion coefficient of η(τ)²η(11τ)², which is an independent construction of the weight-2 newform: 0 mismatches.
- **Torsion and Hasse.** #Ẽ(F_p) ≡ 0 (mod 5) for 11a1 and (mod 3) for 37b1, at every good p ≤ 2000. These are the torsion orders Z/5 and Z/3. Hasse's bound |a_p| ≤ 2√p holds throughout.
- **Prime-power coefficients two ways.** The Euler-product coefficients (α_p^k + β_p^k) log p, including the bad prime, agree with the coefficients of −L′/L computed from the Dirichlet series Σ (a_n/√n) n^{−s} by the log-derivative recursion. The agreement is 2e-31 at n ≤ 500, and the recursion's coefficients vanish to 5e-30 off prime powers.

**Root number and rank** (`gl2_form_mp.py afe`, writes `afe_checks.json`). Λ_a(w) = (√N/2π)^w Γ(w)L(E, w) is computed by the smoothed approximate functional equation Λ_a(w) = Σ a_n[c_n^{−w}Γ(w, c_n t) + εc_n^{w−2}Γ(2 − w, c_n/t)], with c_n = 2πn/√N and mpmath's `gammainc`. At five random w per curve, with Re w ∈ [−1.5, 3.5] and |Im w| ≤ 25:
- **Splitting parameter.** t = 1, 1.3 and 0.75 agree to 1e-46 to 6e-41 (working precision).
- **Functional equation.** Λ_a(w) = +Λ_a(2 − w) holds at t = 1.3 to 1e-46 to 1e-40. With ε = −1 the same test fails by 0.4 to 7, so **ε = +1 for both curves** (*validated*).
- **Independent normalisation check.** At Re w = 5 and 4.5, the AFE agrees with the Dirichlet series summed directly to n = 30000, to 1e-20 and 3e-18, which is the series' truncation level. This is independent of the functional equation.
- **Rank.** L(E, 1) ≠ 0 for both curves, so the analytic rank is 0, and hence the algebraic rank is 0 (Kolyvagin). The 11a1 value 0.2538418608559… agrees with the known L(11a1, 1); the 37b1 value is *computed* only.

**Disclosure.** Before this section was committed, a timing and precision pilot computed λ_min at x/N = 2 and 4 with N_basis = 5k*, both sectors, for both curves (8 values). Its purpose was to choose the scan's precision schedule. It could not affect the curve choice, which the registration fixes. Those 8 values are not used anywhere; the scan recomputes everything.

## 3. Zeros of L(E, s), found independently of the form

`gl2_form_mp.py zeros` writes `zeros_11a1.json` and `zeros_37b1.json`, with logs in `runs/`. The method:
- Z_E(τ) = Λ_a(1 + iτ)/((√N/2π)|Γ(1 + iτ)|) is real for ε = +1, and |Z_E| = |L(E, ½ + iτ)|.
- It is evaluated by the AFE at t = 1, with mpmath's `gammainc`.
- The working precision is 30 + 0.7τ + 10 digits, because the incomplete-gamma sum cancels about 0.68τ digits.
- Zeros are bracketed by sign changes on a grid of step 0.05 up to T = 120 and refined to 30 digits with `findroot`.

**Completeness: the count is exact** (*validated*). By the argument principle, N(T) = (θ_E(T) + arg L_a(1 + iT))/π, with θ_E(T) = T log(√N/2π) + Im log Γ(1 + iT). The arg is tracked continuously along w = σ + iT, from σ = 4.5 (where |L_a − 1| < 1) down to σ = 1. At six heights per curve, each chosen midway between found zeros, N(T) comes out an integer to working precision, and it equals the number of sign changes found:

| curve | zeros in (0, 120] | N(T) = found, at T ≈ 20 / 40 / 60 / 80 / 100 / 120 | first zeros |
|---|---|---|---|
| 11a1 | 121 | 9 / 26 / 47 / 70 / 94 / 121 | 6.362613895, 8.603539619, 10.0355091, 11.45125861 |
| 37b1 | 144 | 13 / 34 / 59 / 85 / 114 / 144 | 3.509102943, 5.449734162, 7.599111771, 9.032345852 |

So every zero with 0 < γ ≤ 120 is on the critical line and simple, and none was missed. The largest arg increment per σ-step was 0.22 rad. The first zero of 11a1 is near 6.36, as expected.

**Procedural note.** The first zero run crashed at τ ≈ 68. `findroot`'s residual test was applied at the inflated working precision, which is tighter than Z's real 30-digit accuracy. The refinement was moved to 30 digits, with a bisection fallback (never triggered), and both curves were rerun from scratch.

## 4. Hard validation 3: the form against the zero sum (pass, *validated*)

`gl2_form_mp.py check` writes `check_zero_sum.json`. It compares Q(c) = cᵀQc against 2Σ_{0<γ≤120}F_c(γ)², using the same edge-vanishing test functions as `decay_law_mp.py check`:
- odd sector: single sine modes k = 2, 4, and b₂ − b₃;
- even sector: b₁ + b₂ and b₂ − b₄.

These have F ~ 1/t², so the truncated sum misses about 2∫_T^∞F²dN. The column "tail" evaluates that integral with the smooth density dN/dt = θ_E′(t)/π. It is the accuracy expected from truncating at T = 120.

**x = 13** (the registered check; prime powers up to 13, including 11a1's bad prime 11):

| curve | test | Q | Q − Σ | rel. error | expected tail | rel. error after the tail |
|---|---|---|---|---|---|---|
| 11a1 | odd k = 2 | 0.586273 | 2.21e-5 | 3.8e-5 | 2.06e-5 | 2.5e-6 |
| 11a1 | odd k = 4 | 4.32232 | 8.90e-5 | 2.1e-5 | 8.30e-5 | 1.4e-6 |
| 11a1 | odd b₂ − b₃ | 1.90938 | 1.38e-4 | 7.3e-5 | 1.29e-4 | 4.9e-6 |
| 11a1 | even b₁ + b₂ | 0.245513 | 1.20e-8 | 4.9e-8 | 1.13e-8 | 2.9e-9 |
| 11a1 | even b₂ − b₄ | 5.89958 | 1.93e-7 | 3.3e-8 | 1.82e-7 | 1.9e-9 |
| 37b1 | odd k = 2 | 3.04848 | 2.13e-5 | 7.0e-6 | 2.34e-5 | 6.8e-7 |
| 37b1 | odd k = 4 | 4.95004 | 8.59e-5 | 1.7e-5 | 9.42e-5 | 1.7e-6 |
| 37b1 | odd b₂ − b₃ | 6.34579 | 1.34e-4 | 2.1e-5 | 1.47e-4 | 2.0e-6 |
| 37b1 | even b₁ + b₂ | 5.76610 | 1.16e-8 | 2.0e-9 | 1.28e-8 | 2.1e-10 |
| 37b1 | even b₂ − b₄ | 8.62686 | 1.87e-7 | 2.2e-8 | 2.07e-7 | 2.3e-9 |

**x = 40** (an extra check, so that 37b1's bad prime 37 enters). Relative errors before → after the tail correction:
- 11a1: 8e-9 to 1.9e-4 → 5e-10 to 8e-6. The largest is for odd k = 2, where Q = 0.038 is small.
- 37b1: 6e-10 to 1.6e-5 → 5e-12 to 3e-8.

**Reading.**
- The relative errors are 2e-9 to 7e-5 at x = 13. They are comparable to the degree-1 checks in `docs/DECAY_LAW.md`, which gave 1.4e-9 to 5e-4 at the same T.
- In all 20 cases, Q − Σ has the sign of the predicted truncation tail and is 0.90–1.07 times its size. Subtracting it lowers the error by 9–600×.
- So the residual is truncation, not a defect of the form.
- **This also checks the conductor convention** (log N·I). A shift δ in the diagonal constant would move the even b₁ + b₂ value by 2δ. The tail-corrected residual bounds δ by about 4e-10 for 11a1 (x = 13) and about 1e-11 for 37b1 (x = 40).

## 5. P-DL4, executed as registered

**Protocol.**
- **Grid.** x/N ∈ {2, 4, 6, 8, 12}, both sectors. N_basis = ⌊f·k*⌋ + 8 with f = 5, 9 and k* = xL/N. The +8 is the convention of `decay_law_mp.py` scan and product. N_basis reaches 535 for 11a1 and 666 for 37b1.
- **Precision.** Every λ goes through `decay_law_mp.stable_lams`: the form is rebuilt at dps + 40 until two consecutive minima agree to 1e-8, and the lower-precision value is kept. The starting precision is 40 + 1.5·4π(x/N)/ln 10, as in `decay_law_mp.py product`. Every value was accepted at the first pair, with dps from 56 (x/N = 2) to 138 (x/N = 12).
- **Eigensolver.** λ is found by inverse iteration with FLINT's ball solver, iterated until the Rayleigh quotient is stable to 1e-14. Non-convergence is flagged, never silently accepted; it did not occur.
- **Data.** `scan_11a1.json` and `scan_37b1.json`.
- **Timing.** The scans were launched at 01:08:06. That is after the curve record (e5fc0f2, 01:01:11) and after a preliminary zero-sum check on zeros up to T ≈ 61–68 (relative errors ≤ 3.6e-4, ≤ 3e-5 after the tail). It is before the final T = 120 check in §4, which finished at 01:32 and passed.

**Precision, convergence and PSD** (*computed*):
- **Precision.** All 20 saved 9k* values were recomputed at dps + 80 (`precision_*.json`). They agree to ≤ 4.1e-15 relative, which is all 15 saved digits.
- **Convergence.** |ln(λ_5k*/λ_9k*)| ≤ 0.069, against ≤ 0.17 for the characters.
- **PSD.** FLINT's full eigensolver agrees with inverse iteration to 1e-24 to 1e-34 relative at the 12 values with N_basis ≤ 120. That checks that the minimum is not a negative eigenvalue missed by inverse iteration. The form is a sum over zeros, so it is PSD whenever GRH holds for L(E, s), but that is checked only up to T = 120 (§3).

| curve | x/N | x | N_basis (5k*, 9k*) | λ_even (9k*) | λ_odd (9k*) | ln(λ₅/λ₉) even / odd |
|---|---|---|---|---|---|---|
| 11a1 | 2 | 22 | 38, 63 | 2.168e-10 | 1.035e-7 | 0.061 / 0.028 |
| 11a1 | 4 | 44 | 83, 144 | 2.774e-16 | 2.747e-13 | 0.060 / 0.026 |
| 11a1 | 6 | 66 | 133, 234 | 3.255e-21 | 3.918e-18 | 0.045 / 0.011 |
| 11a1 | 8 | 88 | 187, 330 | 4.520e-25 | 1.158e-21 | 0.052 / 0.012 |
| 11a1 | 12 | 132 | 300, 535 | 8.307e-32 | 3.608e-28 | 0.047 / 0.013 |
| 37b1 | 2 | 74 | 51, 85 | 5.626e-11 | 6.778e-8 | 0.067 / 0.060 |
| 37b1 | 4 | 148 | 107, 187 | 5.817e-17 | 1.459e-13 | 0.069 / 0.022 |
| 37b1 | 6 | 222 | 170, 299 | 9.114e-22 | 3.560e-18 | 0.048 / 0.020 |
| 37b1 | 8 | 296 | 235, 417 | 8.419e-26 | 4.301e-22 | 0.050 / 0.015 |
| 37b1 | 12 | 444 | 373, 666 | 2.243e-32 | 1.510e-28 | 0.039 / 0.010 |

**Fits** (`gl2_form_mp.py fit`, writes `pdl4_fit.json`; the registered rule is applied to the 9k* values). The model is ln λ = −αv + γ ln v + β, fitted over the five points:

| curve | sector | hypothesis (v) | α/4π | γ | β | max resid | rule [0.8, 1.2] |
|---|---|---|---|---|---|---|---|
| 11a1 | even | A′ (x/N) | 0.2105 | −12.92 | −7.85 | 0.52 | killed |
| 11a1 | even | C′ (x/√N) | 0.0635 | −12.92 | 7.64 | 0.52 | killed |
| 11a1 | even | B′ (√(x/N)) | 2.006 | 2.51 | 12.58 | 0.27 | killed |
| 11a1 | odd | A′ | 0.2010 | −12.34 | −2.31 | 0.60 | killed |
| 11a1 | odd | C′ | 0.0606 | −12.34 | 12.48 | 0.60 | killed |
| 11a1 | odd | B′ | 1.916 | 2.40 | 17.20 | 0.47 | killed |
| 37b1 | even | A′ | 0.2048 | −13.32 | −9.07 | 0.44 | killed |
| 37b1 | even | C′ | 0.0337 | −13.32 | 14.98 | 0.44 | killed |
| 37b1 | even | B′ | 1.953 | 0.99 | 10.83 | 0.19 | killed |
| 37b1 | odd | A′ | 0.2089 | −12.06 | −2.75 | 0.42 | killed |
| 37b1 | odd | C′ | 0.0343 | −12.06 | 19.02 | 0.42 | killed |
| 37b1 | odd | B′ | 1.991 | 4.04 | 17.53 | 0.16 | killed |

The 5k* values give the same twelve verdicts, with α/4π within 0.005 of the 9k* fits.

### Verdict on P-DL4: all three hypotheses are killed, for both curves, in both sectors

No registered scaling survives in any of the four (curve, sector) cases, so the outcome is "neither". The verdict does not depend on the 3-parameter fit. Take the local slopes −Δln λ/(4πΔv) between consecutive grid points (`posthoc.json`):
- In x/N they fall monotonically, from 0.51–0.55 to 0.30. A′ is wrong in form, not only in its coefficient: the fitted γ ≈ −13 absorbs the curvature.
- In √(x/N) they stay between 1.71 and 2.01. So B′ is outside its window at every segment, not only in the fit.
- C′ is excluded by the comparison of the two curves. At equal x/N, ln(λ_11a1/λ_37b1) lies between 0.10 and 1.68, while ln λ itself falls by about 50 across the grid. At equal x/√N, which is what C′ requires, the two curves would sit at x/N values a factor √(37/11) = 1.83 apart.

**Exploratory (registered as exploratory): the prolate index.** With c = 2πv and λ/ℓ_n(c) for n = 0–9, there is no surviving scaling to evaluate it under. Under each registered scaling, the flattest index sits at the edge of the range, with no flat ratio:
- n = 0 under B′, where |ln R(12)/R(2)| = 22–24;
- n = 9 under A′ and C′, where |ln R(12)/R(2)| = 59–700.

None of the registered rates is the e^{−2c} prolate rate. So this question has no answer under the registration (*computed*).

## 6. What the data show instead (post hoc, not registered)

These statements were found after the verdict, on the same 20 values. They are *computed*, and they are candidates for a future registration, not results.

1. **The two curves collapse in x/N.** At equal x/N their minima differ by a factor e^{0.1}–e^{1.7}, in a non-monotone way, while λ spans 22 decades. The conductor enters through x/N, as for χ, and not through x/√N.
2. **The rate is about twice B′'s.** With v = √(x/N), the fitted α/4π is 1.92–2.01. With the rate fixed at 8π, ln λ = −8π√(x/N) + γ ln √(x/N) + β fits with γ = 2.3 (even) and 4.3–4.8 (odd), and max residual 0.16–0.58. In the zero-counting variable T* = 2π√(x/N), this reads λ ≈ e^{−4T*}. Degree 1 had λ ≈ e^{−2T*} with T* = 2πx/q.
3. **A prolate index appears at c = 4π√(x/N) = 2T*.** λ/ℓ_n(c) is flattest at n = 2 (even sector) and n = 4 (odd sector) for both curves, with |ln R(12)/R(2)| = 0.006–0.38. R stays between 0.56 and 5.2 over the grid. This matches n = μ₁ + μ₂ + 2s (μ₁ + μ₂ = ½ + 3/2 = 2), the analogue of the degree-1 rule n = κ + 2s. It is not robust in one case: with the metric taken over x/N = 4 against 12, 11a1's even sector prefers n = 1. All four cases dip at x/N = 6. That dip is a real feature of the minima, since they are precision-stable and basis-converged, not noise.

A test that would make (2) and (3) predictions would use fresh curves and an extended grid, with the windows fixed in advance. For example: α/4π ∈ [1.6, 2.4] in v = √(x/N), and flattest n = 2 and 4 at c = 4π√(x/N), for rank-0 curves with other conductors and reduction types.

## 7. Limits

- **Upper bounds only.** These are Rayleigh–Ritz upper bounds on rh2's periodic bases; nothing is certified. The basis size is tied to A′'s crossover k* = xL/N, which is generous for the observed rate: the 5k* and 9k* values differ by at most e^{0.07}. That does not exclude a basis-dependent factor; for ζ, antiperiodic bases sometimes gave lower values.
- **Not asymptotic.** x/N ≤ 12 means √(x/N) ≤ 3.46 and T* ≤ 22. Five points per (curve, sector) leave 2 degrees of freedom in each fit.
- **Two curves.** Both are rank 0 with ε = +1, and both have split multiplicative reduction at a prime conductor. Nothing here tests additive reduction, higher rank or ε = −1.
- **Large-k coverage.** The zero-sum validation exercises only low modes (k ≤ 6). The large-k entries of the μ = ½ and 3/2 blocks share their closed-form code path with μ = 0 and 1, which are exact against `build_form` up to n = 30 at x up to 444. They are not otherwise independently tested at k in the hundreds.

**Reproduce** (outputs in `data/connes/gl2/`, run logs in `runs/`):

```
.venv/bin/python scripts/gl2_form_mp.py validate --x 13,22,132,444 --n 30 --dps 100 --json data/connes/gl2/validate_arch.json
.venv/bin/python scripts/gl2_form_mp.py curves --json data/connes/gl2/curves.json
.venv/bin/python scripts/gl2_form_mp.py afe --json data/connes/gl2/afe_checks.json
.venv/bin/python scripts/gl2_form_mp.py zeros --curve 11a1 --json data/connes/gl2/zeros_11a1.json     # and 37b1
.venv/bin/python scripts/gl2_form_mp.py check --x 13,40 --json data/connes/gl2/check_zero_sum.json
.venv/bin/python scripts/gl2_form_mp.py scan --curve 11a1 --json data/connes/gl2/scan_11a1.json       # and 37b1
.venv/bin/python scripts/gl2_form_mp.py precision --curve 11a1                                         # and 37b1
.venv/bin/python scripts/gl2_form_mp.py fit --json data/connes/gl2/pdl4_fit.json
.venv/bin/python scripts/gl2_form_mp.py posthoc --json data/connes/gl2/posthoc.json
```

---

# P-DL6: fresh GL(2) objects

This part executes P-DL6, pre-registered at the end of `docs/DECAY_LAW.md` (f5697aa). The tested description is λ ≈ R·ℓ_n(c), with c = 4πv and v = √(x/N). The objects are two fresh elliptic curves and Ramanujan's Δ.

## 8. Code changes for weight k (P-DL4 still reproduces)

`gl2_form_mp.py` now handles a weight-k object with Γ_ℂ(s + (k−1)/2) = Γ_ℝ(s + (k−1)/2)Γ_ℝ(s + (k+1)/2), conductor N and coefficients a_n/n^{(k−1)/2}:
- curves have k = 2, μ = ½ and 3/2;
- Δ has k = 12, N = 1, μ = 11/2 and 13/2.

The changes:
- **AFE.** Λ_a(w) = Σ a_n[c_n^{−w}Γ(w, c_n t) + εc_n^{w−k}Γ(k − w, c_n/t)]. This follows from f(1/y) = εy^k f(y) for f(y) = Σa_n e^{−2πny/√N}.
- **Hardy Z.** Evaluated on Re w = k/2.
- **θ.** θ(T) = T log(√N/2π) + Im log Γ(k/2 + iT).
- **Argument path.** Runs from σ = k/2 + 3.5 down to k/2.
- **Hecke recursion.** Uses p^{k−1}.
- **Euler terms.** α_p + β_p = a_p/p^{(k−1)/2}.
- **Tail density.** θ′(t)/π.

For k = 2 every expression and term count is kept exactly as in P-DL4. A regression run reproduces, after the change:
- the saved λ (11a1, x/N = 2, even, 9k*; 37b1, x/N = 4, odd, 5k*) to all 15 digits;
- a saved zero of 37b1 (γ₃₁ = 37.084727667901671939) to 20 digits;
- Z_E(6.3) for 11a1;
- `curves.json` exactly.

## 9. Objects (recorded before any P-DL6 scan)

**Curves: 15a1 and 19a1, the registration's first choices. Both qualify** (`gl2_form_mp.py objects6`, writes `objects6.json`; `afe`, writes `afe_checks6.json`). Every check below is *validated* unless marked otherwise.

| | 15a1 | 19a1 |
|---|---|---|
| model [a1, a2, a3, a4, a6] | [1, 1, 1, −10, −10] | [0, 1, 1, −9, −15] |
| Δ | 50625 = 3⁴·5⁴ | −6859 = −19³ |
| c4 (prime to N) | 481 | 448 |
| j | 2198.2151 | −13109.111 |
| reduction at p \| N | a₃ = −1 (non-split), a₅ = +1 (split) | a₁₉ = +1 (split) |
| conductor (Tate: k_p < 12, c4 prime to N) | 15 | 19 |
| gcd of #Ẽ(F_p), odd good p ≤ 2000 | 8 (torsion Z/2 × Z/4) | 3 (torsion Z/3) |
| a_p, p = 2, 3, 5, 7, 11, 13 | −1, −1, 1, 0, −4, −2 | 0, −2, 3, −1, 3, −4 |
| independent q-expansion | η(τ)η(3τ)η(5τ)η(15τ): a_n agree for n ≤ 2000 (0 mismatches) | none available |
| ε (AFE) | +1 (ε = −1 fails by 1.3–14) | +1 (ε = −1 fails by 0.5–3.9) |
| L(E, 1) arithmetic = L(½) analytic | 0.3501507605831505058 | 0.45325324449610360358 |

- **AFE checks.** Splitting parameters t = 1, 1.3 and 0.75 agree to 4e-46 to 1.2e-41. Λ(w) = Λ(2 − w) holds to 2e-46 to 4e-42. At Re w = 5 the AFE matches the Dirichlet series to its truncation level, 3e-21 to 3e-20.
- **Euler terms.** They match the Dirichlet log-derivative to 1.4e-31 (15a1) and 1e-31 (19a1), and the log-derivative vanishes off prime powers to 5e-30.
- **Rank 0.** L(E, 1) ≠ 0, so both curves have analytic rank 0, and hence algebraic rank 0 (Kolyvagin).
- **Labels.** That these models carry the Cremona labels 15a1 and 19a1 is from memory, not from a database (*could not verify*). It does not matter for the test: conductors 15 and 19 each have a single isogeny class, and the L-function used is the one computed from this model's a_p.

**Ramanujan's Δ.** τ(n) is the coefficient of qⁿ in q∏(1 − qⁿ)²⁴, computed as q(∏(1 − qⁿ)³)⁸ with Jacobi's identity (FLINT polynomials). Checks:
- **Two routes.** It agrees with the naive product ∏(1 − qⁿ)²⁴ for n ≤ 300 (0 mismatches).
- **Known values.** τ(1..10) = 1, −24, 252, −1472, 4830, −6048, −16744, 84480, −113643, −115920.
- **Ramanujan's congruence.** τ(n) ≡ σ₁₁(n) (mod 691) for every n ≤ 1000.
- **Hecke multiplicativity.** Rebuilding τ(n) from τ(p) alone, by multiplicativity and τ(p^{j+1}) = τ(p)τ(p^j) − p¹¹τ(p^{j−1}), reproduces the q-expansion for every n ≤ 2000. This is the Euler product used in the form.
- **Deligne's bound.** |τ(p)| ≤ 2p^{11/2} for every p ≤ 2000.
- **Euler terms.** They match the log-derivative to 2.5e-31.
- **AFE (weight 12, N = 1).** Splitting parameters agree to 4e-48 to 4e-42, and Λ(w) = +Λ(12 − w) holds to 1.4e-48 to 3e-42, so ε = +1. The ε = −1 variant fails by 1.2–3.7. The Dirichlet series at Re w = 10 and 9.5 agrees to 1.5e-20 and 2.6e-18.
- **Central value.** L(Δ, 6) = 0.79212283864603056936 (arithmetic normalisation), which is L(Δ, ½) in the analytic normalisation.

**Zero-sum gate, fixed before any zero or form value of these objects is computed.** P-DL4 achieved relative errors ≤ 7.3e-5 at x = 13 and T = 120, with (Q − Σ)/tail between 0.90 and 1.07. An object passes, and may be scanned, only if at x = 13 and T = 120 all five edge-vanishing tests (odd: k = 2, k = 4, b₂ − b₃; even: b₁ + b₂, b₂ − b₄) satisfy both:
- relative error |Q − Σ|/|Q| ≤ 1e-4;
- (Q − Σ)/(predicted tail) ∈ [0.8, 1.25].

Extra checks are reported but do not gate: x = 40 for the curves (19a1's bad prime enters), and x = 4 and 16 for Δ (its grid ends).

## 10. Zeros (*validated*)

The method is that of §3, at weight k. The argument-principle count is an exact integer, and it equals the number of sign changes found, at six heights per object:

| object | zeros in (0, 120] | N(T) = found, at T ≈ 20 / 40 / 60 / 80 / 100 / 120 | first zeros |
|---|---|---|---|
| 15a1 | 127 | 10 / 28 / 50 / 74 / 100 / 127 | 5.239203926, 7.664880134, 9.523451676, 10.67892245 |
| 19a1 | 131 | 11 / 29 / 52 / 77 / 103 / 131 | 5.039123554, 6.390843061, 9.180168549, 10.24320446 |
| Δ | 78 | 4 / 13 / 27 / 42 / 59 / 78 | 9.2223793999211025222, 13.90754986, 17.44277698, 19.65651314 |

- **Δ's first zero** is 9.22237939992110252 (*computed*). It agrees with the expected t ≈ 9.22, which was checked, not assumed. Z(0) = L(Δ, 6) = 0.79212284, consistent with the central value in §9.
- **A zero just below the cutoff.** Δ has a zero at 119.99454 = T − 0.005. The count at T = 120 was therefore repeated with a finer σ-step: N(120) = 78 and N(119.54) = 77, both exact.
- **Diagnostic run to T = 240.** Δ was later rerun to T = 240 (205 zeros, counts exact at 6 heights up to 240). It reproduces all 78 zeros below 120 exactly.

## 11. The zero-sum gate: the curves pass, Δ fails as fixed in advance

`gl2_form_mp.py check --gate` writes `check6_*.json`. The gate, fixed in §9, is applied at x = 13 and T = 120:

| object | test (x = 13) | Q | rel. error | (Q − Σ)/tail | gate |
|---|---|---|---|---|---|
| 15a1 | odd k = 2 / k = 4 / b₂ − b₃ | 2.299 / 4.250 / 4.538 | 8.9e-6 / 1.9e-5 / 2.8e-5 | 0.958 | |
| 15a1 | even b₁ + b₂ / b₂ − b₄ | 1.820 / 6.982 | 6.1e-9 / 2.6e-8 | 0.956–0.957 | **pass** |
| 19a1 | odd k = 2 / k = 4 / b₂ − b₃ | 2.974 / 4.715 / 3.015 | 6.7e-6 / 1.7e-5 / 4.1e-5 | 0.903–0.904 | |
| 19a1 | even b₁ + b₂ / b₂ − b₄ | 2.467 / 8.787 | 4.3e-9 / 2.0e-8 | 0.892 | **pass** |
| Δ | odd k = 2 / k = 4 / b₂ − b₃ | **0.0234** / 2.315 / **0.323** | **6.0e-4** / 2.4e-5 / **2.7e-4** | 0.926 | |
| Δ | even b₁ + b₂ / b₂ − b₄ | 0.00366 / 2.805 | 2.0e-6 / 4.2e-8 | 0.905 | **FAIL** |

**Extra, non-gating checks.**
- Curves at x = 40: relative errors up to 1.4e-4 (15a1) and 1.6e-4 (19a1), both on the small-Q odd k = 2 test (Q ≈ 0.05). The tail ratios are 0.91–0.94 and 1.02.
- Δ at x = 4: relative errors up to 1.8e-4, tail ratio 0.83–0.86.
- Δ at x = 16: relative errors up to 2.6e-3 (Q = 0.0046), tail ratio 0.99–1.01.

**Δ fails the gate on its relative-error clause, in 2 of the 5 tests.** It passes the tail-ratio clause in all 5. Reading (*computed*):
- **The failing tests have small Q.** Δ's first zero is at 9.22, above the frequencies of these low modes at x = 13 (ω₂ = 4.9), so F is small on every zero.
- **The absolute residuals are smaller than P-DL4's.** They are 1.4e-5 and 8.8e-5, against up to 1.4e-4 for 11a1 and 37b1 at x = 13.
- **They are truncation.** Each is 0.91–0.93× the predicted truncation tail.
- **The rule was mis-stated, but it was fixed in advance.** A relative-error rule penalises small Q. P-DL4's own 11a1 check at x = 40 had relative error 1.9e-4 on such a test. The gate was nonetheless fixed before the data, so it is applied as written: **Δ was not scanned.**

**Post-failure diagnostic, not a substitute for the gate.** The same check with Δ's zeros to T = 240 (`check6_delta_T240_diagnostic.json`) shows:
- at x = 13 the residuals fall by 6.4× (odd tests) and 25× (even tests), each by the same factor as its own predicted tail;
- each residual is still 0.96× the predicted tail;
- the largest relative error becomes 9.4e-5 (odd k = 2); the other four are ≤ 4.3e-5;
- at x = 4 and 16 the residuals are 0.92–0.93× and 1.03–1.05× the predicted tails.

So the behaviour is that of a correct form truncated at T. Whether to accept Δ's form, with the gate run at T = 240 or restated as the tail-ratio clause, is the registrant's decision; it is not taken here. If accepted, `scan6 --obj delta` and `fit6` complete P-DL6 in a few minutes.

**Ordering.** Objects committed at 01:51:26 (276a496); zero runs from 01:51:43; Δ gate failed at 02:06:32; curve gates passed at 02:14:50; curve scans launched at 02:15:19.

## 12. P-DL6 scans of the two curves (*computed*)

`scan6_15a1.json` and `scan6_19a1.json`. The protocol:
- **Grid.** v = √(x/N) ∈ {1.5, 2, 2.5, 3, 3.5}, both sectors.
- **Basis.** N_basis = ⌊f·k⌋ + 8 with f = 5, 9 and k = 2vL, with nmin 24 and n₂ ≥ n₁ + 8, as in `decay_law_mp.py product-v`.
- **Precision.** `stable_lams` (dps + 40, agreement to 1e-8), starting at dps₀ = 40 + 1.5·8πv/ln 10. Every value was accepted at the first pair, with dps 64–97.
- **Eigensolver.** Inverse iteration to 1e-14.

**Precision, convergence and PSD:**
- **Precision recheck.** All 20 9k values recomputed at dps + 80 (`precision6_*.json`) agree to ≤ 3.4e-15.
- **Convergence.** |ln(λ_5k/λ_9k)| ≤ 0.079.
- **PSD.** FLINT's full eigensolver agrees with inverse iteration to ≤ 2.1e-26 at the 12 values with N_basis ≤ 120.

| curve | v | x | N_basis (5k, 9k) | λ_even (9k) | λ_odd (9k) | ln(λ₅/λ₉) even / odd |
|---|---|---|---|---|---|---|
| 15a1 | 1.5 | 33.75 | 60, 103 | 2.232e-11 | 1.270e-8 | 0.037 / 0.026 |
| 15a1 | 2 | 60 | 89, 155 | 2.011e-16 | 2.866e-13 | 0.067 / 0.023 |
| 15a1 | 2.5 | 93.75 | 121, 212 | 1.339e-21 | 2.625e-18 | 0.076 / 0.012 |
| 15a1 | 3 | 135 | 155, 272 | 4.928e-27 | 1.570e-23 | 0.069 / 0.008 |
| 15a1 | 3.5 | 183.75 | 190, 336 | 2.303e-32 | 1.094e-28 | 0.068 / 0.017 |
| 19a1 | 1.5 | 42.75 | 64, 109 | 1.395e-11 | 7.301e-9 | 0.032 / 0.025 |
| 19a1 | 2 | 76 | 94, 163 | 1.243e-16 | 1.494e-13 | 0.058 / 0.017 |
| 19a1 | 2.5 | 118.75 | 127, 222 | 6.745e-22 | 1.509e-18 | 0.069 / 0.013 |
| 19a1 | 3 | 171 | 162, 285 | 2.593e-27 | 1.237e-23 | 0.079 / 0.013 |
| 19a1 | 3.5 | 232.75 | 198, 351 | 1.145e-32 | 6.463e-29 | 0.068 / 0.016 |

## 13. P-DL6a and P-DL6b on the curve cases

`gl2_form_mp.py fit6 --objects 15a1,19a1` writes `pdl6_fit_curves.json`. The rules are applied as registered, to the four curve cases. The two Δ cases are missing (§11).

**P-DL6a (rate).** The model is ln λ = −αv + γ ln v + β, fitted to the 9k values:

| curve | sector | α/4π | γ | β | max resid | linear α/4π | in [1.7, 2.3] |
|---|---|---|---|---|---|---|---|
| 15a1 | even | 2.138 | 6.33 | 13.23 | 0.13 | 1.927 | yes |
| 15a1 | odd | 2.161 | 9.39 | 18.79 | 0.16 | 1.847 | yes |
| 19a1 | even | 2.134 | 5.90 | 12.86 | 0.11 | 1.937 | yes |
| 19a1 | odd | 2.157 | 9.51 | 18.07 | 0.06 | 1.839 | yes |

The 5k values give 2.140–2.158, so all four curve cases lie inside the predicted window and none is near the kill range. **P-DL6a holds for the four curve cases. Its overall verdict needs the two Δ cases.**

**P-DL6b (prolate index).** c = 4πv, n ∈ 0..20, with the flattest metric |ln R(v₅)/R(v₂)|, v₂ = 2 and v₅ = 3.5:

| curve | sector | predicted n | flattest n | metric (flattest) | metric (predicted) | off by |
|---|---|---|---|---|---|---|
| 15a1 | even | 2 | **1** | 0.154 | 0.405 | −1 |
| 15a1 | odd | 4 | **3** | 0.239 | 0.321 | −1 |
| 19a1 | even | 2 | **1** | 0.064 | 0.623 | −1 |
| 19a1 | odd | 4 | 4 | 0.196 | 0.196 | 0 |

**Status of P-DL6b:**
- 1 of the 4 curve cases is exact, the other 3 are off by −1, and none is off by 2 or more. So the curve kill clause does not fire.
- Success needs at least 4 of 6 exact. With 1 exact among the curves, at most 3 of 6 can be exact, so **P-DL6b cannot succeed whatever Δ gives.**
- Δ decides only between "killed" (a flattest n ≤ 8) and "not killed, prediction partly missed".

**Resolution of the index.** Raising n by 1 multiplies ℓ_n(c) by 8c/(n + 1). So the signed ln R(v₅)/R(v₂) moves by exactly ln(v₅/v₂) = 0.56 per unit of n, and the metric resolves n only to about ±0.5. Interpolating the signed value to zero gives an effective continuous index of 1.28 / 3.43 (15a1) and 0.89 / 3.65 (19a1), against the predicted 2 / 4. This is post hoc, not registered. The curves lean 0.4–1.1 below the rule n = μ₁ + μ₂ + 2s, consistently across both curves.

## 14. Limits and fragile points (P-DL6)

- **Δ is unscanned.** Its form failed the gate as written, even though every diagnostic says the form is right (§11). Both overall verdicts wait on the registrant's decision.
- **The index test is coarse.** The registered metric uses two grid points, 0.56 apart in the signed quantity per unit of n. Bumps of ±0.2 in ln R between v₂ and v₅ shift the answer by one index, and three of the four curve cases are one index low. The other metric used in P-DL4 (first against last point) would give a different answer in some cases.
- **Fit freedom.** P-DL6a's window is wide compared with the trade-off between α and γ in a 3-parameter fit on five points. Here α/4π = 2.13–2.16 with γ = 6–10, while P-DL4's √(x/N) fits gave α/4π = 1.92–2.01 with γ = 1–4, over a slightly different range of v. The linear slopes, 1.84–1.94, are below 2.
- **Upper bounds only,** as in P-DL4: the periodic basis, v ≤ 3.5, and two curves, both with multiplicative reduction.

**Reproduce (P-DL6):**

```
.venv/bin/python scripts/gl2_form_mp.py objects6 --json data/connes/gl2/objects6.json
.venv/bin/python scripts/gl2_form_mp.py afe --curves 15a1,19a1,delta --json data/connes/gl2/afe_checks6.json
.venv/bin/python scripts/gl2_form_mp.py zeros --curve 15a1 --json data/connes/gl2/zeros_15a1.json            # and 19a1, delta
.venv/bin/python scripts/gl2_form_mp.py check --curves 15a1 --x 13,40 --gate --json data/connes/gl2/check6_15a1.json   # and 19a1
.venv/bin/python scripts/gl2_form_mp.py check --curves delta --x 13,4,16 --gate --json data/connes/gl2/check6_delta.json
.venv/bin/python scripts/gl2_form_mp.py zeros --curve delta --tmax 240 --count-at 40,80,120,160,200,240 --json data/connes/gl2/zeros_delta_T240.json
.venv/bin/python scripts/gl2_form_mp.py check --curves delta --x 13,4,16 --zeros-tag _T240 --json data/connes/gl2/check6_delta_T240_diagnostic.json
.venv/bin/python scripts/gl2_form_mp.py scan6 --obj 15a1                                                       # and 19a1
.venv/bin/python scripts/gl2_form_mp.py precision --curve 15a1 --file data/connes/gl2/scan6_15a1.json --out data/connes/gl2/precision6_15a1.json
.venv/bin/python scripts/gl2_form_mp.py fit6 --objects 15a1,19a1 --json data/connes/gl2/pdl6_fit_curves.json
# only if Δ's form is accepted:  scan6 --obj delta;  fit6 --json data/connes/gl2/pdl6_fit.json
```

## 15. Δ scan and the final P-DL6 verdict (main session)

After the gate decision recorded in `docs/DECAY_LAW.md` (10d4b24, made before any Δ scan existed), the main session ran:
- `scan6 --obj delta --f1 5 --f2 9 --nmin 24 --acb-max 120`, the same flags as the curve scans. Output: `scan6_delta.json`, log `runs/scan6_delta.run.txt`, convergence ≤ 0.076.
- `fit6 --objects 15a1,19a1,delta` (`pdl6_fit.json`).

**Verdicts.** P-DL6a holds (α/4π 2.102–2.168). P-DL6b is not killed and the prediction partly missed: Δ's flattest n is 14 (even) and 16 (odd), against the predicted 12 and 14. The table is in `docs/DECAY_LAW.md`, "P-DL6 result".

---

# P-DL7: the prolate index against the weight

This part executes P-DL7, pre-registered at the end of `docs/DECAY_LAW.md` (6fd392d). The objects are the level-1 Hecke eigenforms f_k = Δ·E_{k−12}, for k ∈ {16, 18, 20, 22, 26}.

## 16. Code changes and the ε = −1 case

`gl2_form_mp.py` gains the objects `f16 … f26`:
- weight k, N = 1, μ = (k − 1)/2 and (k + 1)/2, root number i^k;
- a_n from the q-expansion Δ·E₄^a E₆^b with 4a + 6b = k − 12, using FLINT integer polynomials.

Three of the forms, k = 18, 22 and 26, have ε = −1. Four changes handle that case:
- **Z on the line.** With ε = −1, Λ(k/2 + iτ) = 2i Im Σ a_n c_n^{−w}Γ(w, c_n), so Z(τ) := Λ/(i·norm) = 2 Im Σ/norm. It is real and odd in τ.
- **Central multiplicity from the count.** The argument principle on the symmetric rectangle gives (θ + arg L)/π = N(T) + m₀/2, where m₀ is the order of the zero at s = ½. Comparing the count with the number of sign changes therefore measures m₀; it must come out the same at every height.
- **Direct order check.** Z(2e-4)/Z(1e-4) should be 2 for a simple zero (8 for a triple). Separately, L′(½) is computed with two splitting parameters.
- **Zero sum.** It includes m₀·F(0)² for the central zero, counted once. F(0) = 0 for the odd-sector tests.

A regression after these changes reproduces P-DL4 (§8's checks) and P-DL6: Δ's Z(9.2) and a saved Δ λ (v = 2.5, odd, 9k), to all printed digits.

## 17. Forms (recorded before any P-DL7 zero, gate or scan) (*validated*)

`gl2_form_mp.py objects7` writes `objects7.json`; `afe` writes `afe_checks7.json`.

| k | ε = i^k | a₁…a₅ | numerator of B_k/k | congruence a_n ≡ σ_{k−1}(n), n ≤ 1000 | L(½) (analytic) |
|---|---|---|---|---|---|
| 16 | +1 | 1, 216, −3348, 13888, 52110 | 3617 | 0 failures | 1.5205616690847 |
| 18 | −1 | 1, −528, −4284, 147712, −1025850 | 43867 | 0 failures | 0 (L′(½) = 0.64960489635628) |
| 20 | +1 | 1, 456, 50652, −316352, −2377410 | 283 · 617 | 0 failures (both primes) | 1.9817354054335 |
| 22 | −1 | 1, −288, −128844, −2014208, 21640950 | 131 · 593 | 0 failures (both primes) | 0 (L′(½) = 1.02878654277693) |
| 26 | −1 | 1, −48, −195804, −33552128, −741989850 | 657931 | 0 failures | 0 (L′(½) = 1.37624801726741) |

Checks, for every form:
- **q-expansion two ways.** The E₄, E₆ products agree with Δ·E_{k−12}, where E_{k−12} is built directly from its Bernoulli coefficient, 1 − (2m/B_m)Σσ_{m−1}(n)qⁿ. Since M_m(1) is one-dimensional the two must coincide, and they do for every n ≤ 2000: 0 mismatches.
- **a₂.** a₂ = −24 + c₁(E_{k−12}), with c₁ = 240, −504, 480, −264, −24.
  - *Correction:* the a₂ values I recalled from memory had the wrong sign for k = 20 and 22. The computed 456 and −288 are what the formula gives, and the congruences confirm them.
- **Hecke multiplicativity.** Rebuilding a_n from a_p alone, with p^{k−1}, reproduces the q-expansion for n ≤ 2000. This is the Euler product used in the form.
- **Deligne's bound.** |a_p| ≤ 2p^{(k−1)/2} for every p ≤ 2000.
- **Eisenstein congruence.** a_n ≡ σ_{k−1}(n) holds modulo every prime dividing the numerator of B_k/k, computed exactly. This is the analogue of Ramanujan's 691.
- **Euler terms.** They match the log-derivative to ≤ 2.2e-31.

AFE checks with the form's own ε:
- **Splitting parameters.** t = 1, 1.3 and 0.75 agree to 2.6e-52 to 1.4e-41. A wrong ε would spoil this as well.
- **Functional equation.** Λ(k − w) = εΛ(w) holds to 1.3e-52 to 1.4e-42. With the opposite sign it fails by 0.18–4.2.
- **Normalisation.** The Dirichlet series at Re w = k/2 + 4 and k/2 + 3.5 agrees to its truncation level (≤ 3.5e-18).
- **ε = −1 forms.** L(½) is 0 exactly at t = 1, which is forced by symmetry, and about 1e-41 at t = 1.3. L′(½) ≠ 0, with the two splitting parameters agreeing to 15 digits, so the central zero is simple.

## 18. Zeros to T = 240, the central zero, and the gate (*validated*)

**Zeros.** In every case the count is exact at six heights, T ≈ 40, 80, …, 240:

| form | ε | zeros in (0, 240] | first zeros | (θ + arg L)/π − found | central zero |
|---|---|---|---|---|---|
| f₁₆ | +1 | 205 | 5.265020228, 11.82395546 | 0 at all 6 heights | none (L(½) = 1.52) |
| f₁₈ | −1 | 206 | 8.14161047, 11.12333425 | ½ at all 6 heights, so m₀ = 1 | simple: Z(2e-4)/Z(1e-4) = 2.000; Z′(0) = 0.64960 = L′(½) |
| f₂₀ | +1 | 207 | 3.607198237, 8.677106764 | 0 | none (L(½) = 1.98) |
| f₂₂ | −1 | 206 | 5.621476062, 10.0025782 | ½, so m₀ = 1 | simple: ratio 2.000; Z′(0) = 1.02879 |
| f₂₆ | −1 | 207 | 4.435321319, 8.332583171 | ½, so m₀ = 1 | simple: ratio 2.000; Z′(0) = 1.37625 |

**Bug caught before any count was used.** A smoke test found it: `count_zeros_argument` had evaluated the AFE with ε = +1 for every object. That made the result an integer even for ε = −1 forms, where it must be half-integral. The fix passes each form's ε (§16). It left the ε = +1 objects unchanged (Δ's N(10) = 1 before and after), and it predates every count reported here.

**Gate** (`check --gate7`, `check7_*.json`), as registered: at x = 13 with zeros to T = 240, every test needs (Q − Σ)/tail ∈ [0.8, 1.25] and |Q − Σ| ≤ 2e-4.

| form | max \|Q − Σ\| (5 tests) | (Q − Σ)/tail | gate | committed (before scan) |
|---|---|---|---|---|
| f₁₆ | 1.32e-5 | 0.919–0.924 | pass | 2e8394c, 03:40:05 |
| f₁₈ | 1.55e-5 | 1.09 | pass | 2e8394c, 03:40:05 |
| f₂₀ | 1.60e-5 | 1.12 | pass | 2e8394c, 03:40:05 |
| f₂₂ | 1.41e-5 | 0.988–0.990 | pass | 7995279, 03:44:16 |
| f₂₆ | 1.31e-5 | 0.913–0.915 | pass | 3589b5a, 04:05:28 |

**The central zero inside the zero sum** (`central7`, `central7_{a,b}.json`; a reported check, not part of the gate). The five gate tests contain no b₀ component, so F(0) = 0 for each of them and none can detect γ = 0. The edge-vanishing combination c = (1, 1/√2, 0, …) has F(0)² = L = 2.56495:
- **ε = −1 forms.** Q exceeds Σ_{γ>0}2F(γ)² by exactly F(0)² = 2.56495. With the central zero counted once, the residual is 1.6–1.9e-11, equal to the predicted tail (ratios 0.91–1.09). A triple zero would leave −2F(0)².
- **ε = +1 forms.** The residual is the tail with no central term (1.6–1.9e-11).

So the form "knows" the central zero through the explicit formula alone, with multiplicity exactly 1.

## 19. Scans (*computed*)

The protocol is P-DL6's:
- v ∈ {2, 2.5, 3, 3.5, 4} with x = v², both sectors;
- N_basis = ⌊f·k′⌋ + 8 with f = 5, 9 and k′ = 2vL;
- `stable_lams` (dps + 40, agreement to 1e-8), starting at dps₀ = 40 + 1.5·8πv/ln 10;
- the scans themselves are `scan7_f*.json`.

Scan start times:
- f₁₆ and f₁₈: 03:40:12;
- f₂₀: 03:44:16 (failed, below), rerun at 04:00:17;
- f₂₂: 03:44:44 (failed), rerun at 04:00:17;
- f₂₆: 04:05:37.

**Eigensolver fallback.** The first f₂₀ and f₂₂ scans stopped because inverse iteration did not converge in 40 steps at v = 2, at any precision up to 600 digits, so `stable_lams` correctly refused a value. These points are far from asymptotic: λ_min = O(1), and the two lowest eigenvalues are close (f₂₀ odd: 1.1555 against 1.2766; f₂₂ even: 1.090 against 1.379).
- **The fix.** `--method inv+full` keeps inverse iteration but falls back to FLINT's full eigensolver (`decay_law_mp.lam_min_acb`, as used for DH) whenever inverse iteration fails. Fallbacks are recorded per row.
- **Where it fired.** f₂₀ v = 2 odd; f₂₂ v = 2 even; f₂₆ v = 2 both sectors and v = 2.5 even. In each case it fired in both precision rounds and at both basis sizes.
- **No effect elsewhere.** f₁₆ and f₁₈ ran with inverse iteration only and had no non-convergence, so they are what the fallback method would give.
- **The measured quantity is unchanged.** It is λ_min, precision-stable as before. The method is not a registered setting.

**Precision, convergence and PSD:**
- **Precision.** All 50 9k values, recomputed at dps + 80 with the same method (`precision7_*.json`), agree to ≤ 4.2e-15.
- **Convergence.** |ln(λ_5k/λ_9k)| ≤ 0.054.
- **PSD.** FLINT's full eigensolver agrees with inverse iteration to ≤ 1.2e-15 wherever both ran (N_basis ≤ 120).

**The forms are pre-asymptotic at the low end of the grid.** At v = 2 (x = 4, L = 1.39), λ_min is O(10⁻³–1) for k ≥ 16, against 1.3e-5 for Δ. At v = 4 it is 3e-19 (f₁₆ even) to 9.5e-8 (f₂₆ even), falling more slowly as k grows.

## 20. P-DL7 verdicts: all three killed (*computed*, `fit7`, `pdl7_fit.json`)

**P-DL7a (rate): KILLED.** The model is ln λ = −αv + γ ln v + β on the 9k values (the 5k values agree to ±0.011):

| k | sector | α/4π | γ | max resid | linear α/4π | status |
|---|---|---|---|---|---|---|
| 16 | even | 2.125 | 23.98 | 0.27 | 1.47 | in [1.7, 2.3] |
| 16 | odd | 2.549 | 44.16 | 0.09 | 1.34 | outside 1.7–2.3, not killed |
| 18 | even | 2.399 | 42.66 | 0.45 | 1.23 | outside 1.7–2.3, not killed |
| 18 | odd | 2.114 | 26.97 | 0.30 | 1.37 | in [1.7, 2.3] |
| 20 | even | 2.499 | 44.21 | 0.28 | 1.29 | outside 1.7–2.3, not killed |
| 20 | odd | **2.915** | 66.72 | 0.50 | 1.09 | **kill** |
| 22 | even | **3.162** | 80.43 | 0.33 | 0.96 | **kill** |
| 22 | odd | 2.525 | 49.40 | 0.22 | 1.17 | outside 1.7–2.3, not killed |
| 26 | even | **3.115** | 89.23 | 0.39 | 0.67 | **kill** |
| 26 | odd | **2.948** | 73.71 | 0.29 | 0.93 | **kill** |

Four cases lie outside [1.4, 2.6]. As k grows, the log-corrected slope rises and the linear slope falls. γ reaches 24–89: the five-point fit absorbs the pre-asymptotic curvature into γ ln v and overshoots α. The forms have not reached a constant rate on v ≤ 4. In this registration that counts as a kill; whether the rate tends to 8π at larger v is not tested here.

**P-DL7b (index slope): KILLED.** c = 4πv, n ∈ 0..40, metric |ln R(v₅)/R(v₂)|:

| k | 12 (Δ, `pdl6_fit.json`) | 16 | 18 | 20 | 22 | 26 |
|---|---|---|---|---|---|---|
| flattest n_even | 14 | 19 | 26 | 24 | 33 | **40** (edge of range; metric 2.22, not flat) |

- The fit n_even = a·k + b gives a = 1.881, b = −9.75 (max residual 3.9), outside the kill band [0.6, 1.5].
- f₂₆'s even index sits at the edge of the registered range, and no n ≤ 40 is flat there. Its true value is therefore ≥ 40, which would only steepen a.

**P-DL7c (sector offset): KILLED.**

| k | ε | n_even | n_odd | n_odd − n_even |
|---|---|---|---|---|
| 12 (Δ, reference) | +1 | 14 | 16 | 2 |
| 16 | +1 | 19 | 22 | 3 |
| 18 | −1 | 26 | 22 | **−4** |
| 20 | +1 | 24 | 29 | **5** |
| 22 | −1 | 33 | 28 | **−5** |
| 26 | −1 | 40 | 35 | **−5** |

The kill fires for k = 18, 20, 22 and 26. Every negative offset belongs to an ε = −1 form, and no ε = +1 form has one.

## 21. What the data show instead (post hoc, not registered)

1. **The sector order is set by the root number.** At v = 4, ln(λ_odd/λ_even) is +7.35, +7.46 and +7.52 for the ε = +1 forms (Δ, f₁₆, f₂₀). For the ε = −1 forms (f₁₈, f₂₂, f₂₆) it is −6.23, −6.23 and −6.23.
   - The sign flips with ε; its size hardly depends on k.
   - A natural reading: an odd test function has F̂(0) = 0 automatically. With ε = −1 the central zero is therefore free for the odd sector and a constraint for the even one, and with ε = +1 it is the reverse. The rule n = Σμ + 2s has no ε in it, which is why P-DL7c fails for exactly the ε = −1 forms.
2. **Pre-asymptotic onset grows with the weight.** At x = 4, λ_min ≈ 1e-5 for Δ, but O(10⁻³–1) for k ≥ 16. The decay starts later as k rises, because the first zeros move up with the gamma factor. A grid of fixed v therefore samples the higher weights earlier in their decay. That depresses the linear slopes (0.67–1.47) and inflates both the log-corrected slopes and the flattest indices.
3. **The flattest index grows faster than k.** The even sector gives 14, 19, 26, 24, 33, ≥ 40, that is n − k ≈ +2 to +14, rising with k. The registered rule n = k + 2s fits Δ and f₁₆ within 2–3, but not beyond. Given point 2, this may be a pre-asymptotic effect, not a property of the index; testing it would need v ranges that scale with the first zero.

## 22. Limits and fragile points (P-DL7)

- **The grid is pre-asymptotic for k ≥ 18.** It is fixed in v (x ≤ 16), and the high-weight forms are still in their O(1)-to-decay transition at v ≤ 2.5. All three kills rest on this grid. They are kills of the registered predictions on this grid, not a test of the asymptotic rate.
- **The eigensolver fallback.** Five (form, v, sector) points needed FLINT's full eigensolver. These are the v ≤ 2.5 points where λ = O(1). The values are precision-stable and agree with `mp.eigsy` where checked. v = 2 enters only the rate fit, not the index metric (v₂ = 2.5, v₅ = 4), but f₂₆'s v = 2.5 even point enters both.
- **The range cap.** f₂₆'s even-sector index hits the registered cap n = 40 without a flat minimum. The registered metric still returns 40, and a larger true value would only steepen P-DL7b's slope.

**Reproduce (P-DL7):**

```
.venv/bin/python scripts/gl2_form_mp.py objects7 --json data/connes/gl2/objects7.json
.venv/bin/python scripts/gl2_form_mp.py afe --curves f16,f18,f20,f22,f26 --json data/connes/gl2/afe_checks7.json
.venv/bin/python scripts/gl2_form_mp.py zeros --curve f16 --tmax 240 --count-at 40,80,120,160,200,240 --json data/connes/gl2/zeros_f16.json   # and f18 … f26
.venv/bin/python scripts/gl2_form_mp.py check --curves f16 --x 13 --gate7 --json data/connes/gl2/check7_f16.json                           # each form
.venv/bin/python scripts/gl2_form_mp.py central7 --objects f16,f18,f20,f22 --json data/connes/gl2/central7_a.json                          # and f26 → central7_b
.venv/bin/python scripts/gl2_form_mp.py scan6 --obj f16 --json data/connes/gl2/scan7_f16.json                                              # f16, f18: inverse iteration
.venv/bin/python scripts/gl2_form_mp.py scan6 --obj f20 --method inv+full --json data/connes/gl2/scan7_f20.json                            # f20, f22, f26
.venv/bin/python scripts/gl2_form_mp.py precision --curve f16 --file data/connes/gl2/scan7_f16.json --out data/connes/gl2/precision7_f16.json
.venv/bin/python scripts/gl2_form_mp.py fit7 --json data/connes/gl2/pdl7_fit.json
```

---

# P-DL9 (agent's part): the fresh weight-2 object with ε = −1, 37a1

P-DL9 is registered at the end of `docs/DECAY_LAW.md` (40f46ee, 04:12:15). It has two parts:
- **The six level-1 forms** on v ∈ {6, …, 10}. The main session runs these (`scan9_{delta,f16,f18,f20,f22,f26}.json`); they are not touched here.
- **The rank-1 curve 37a1,** y² + y = x³ − x, on v = √(x/37) ∈ {2, 2.5, 3, 3.5, 4}. This part runs it, through the P-DL7 gate with the central zero included, and supplies `fit9` for all seven objects.

## 23. 37a1: verification (recorded before any 37a1 gate or scan) (*validated*)

**Code.** `CURVES["37a1"]` = [0, 0, 1, −1, 0], N = 37, ε = −1. The value of ε is entered as data and then checked by the AFE below. Weight 2 with ε = −1 needed two paths that P-DL4 had hard-wired to ε = +1:
- **Z on the line:** Z = 2 Im Σ/norm, which is odd in τ.
- **The central value in `afe`:** it now uses ε and computes L′(½).

The ε = +1 curves keep their exact P-DL4 expressions. A regression reproduces the saved P-DL4 λ values and zero, Δ's Z(9.2) and f18's Z(0.5).

**Splitting the zero search.** `zeros --tmin/--chunk` runs one τ-range, and `zeros-merge` joins contiguous chunks. The merge removes duplicates at the boundaries (a chunk can overrun its end by one step), then runs the central check and the argument-principle counts exactly as `zeros` does.

**Checks** (`objects6 --objects 37a1`, writes `objects9_37a1.json`; `afe --curves 37a1`, writes `afe_checks9_37a1.json`):
- **Conductor.** Δ = 37 and c4 = 48 (prime to 37), so the reduction at 37 is multiplicative and the conductor is 37 (Tate). The reduction is non-split: a₃₇ = −1.
- **a_p.** Point counting gives a₂, a₃, a₅ = −2, −3, −2, matching the values supplied with the task. Also: a₇ = −1, a₁₁ = −5, a₁₃ = −2, a₁₇ = 0, a₁₉ = 0, a₂₃ = 2.
- **Torsion and Hasse.** The gcd of #Ẽ(F_p) over good p ≤ 2000 is 1, consistent with trivial torsion. Hasse's bound holds.
- **Euler terms.** They match the Dirichlet log-derivative to 2.4e-31, and the log-derivative vanishes off prime powers to 1.4e-29.
- **AFE with ε = −1.** Splitting parameters t = 1, 1.3 and 0.75 agree to 6e-45 to 5.8e-41. Λ(2 − w) = −Λ(w) holds to 1.7e-45 to 2.4e-41. **With ε = +1 the same test fails** by 0.018–1.45. The smallest failure, 0.018, is at the near-real point w = 1.75 − 0.08i, and it is still 10³⁹ times the working accuracy. The Dirichlet series at Re w = 5 and 4.5 agrees to its truncation level (9e-20 to 1.7e-17).
- **Rank 1.** L(E, ½) = 0, as ε = −1 forces: exactly 0 at t = 1, and −1.8e-42 at t = 1.3. L′(E, ½) = 0.305999773834052, with the two splitting parameters agreeing to 15 digits. That matches the known value of L′(37a1, 1); the analytic rank is 1.

## 24. 37a1: zeros, gate, scan and the P-DL9 evaluation

**Zeros to T = 240** (*validated*; `zeros_37a1.json`). The search ran as four τ-chunks, [0, 140], [140, 182], [182, 214] and [214, 240], merged by `zeros-merge`. Before the merge they held 174 + 67 + 54 + 45 = 340 zeros, and none were duplicated at the boundaries.
- **First zeros:** 5.003170014, 6.870391217, 8.014330808, 9.933098354.
- **Count.** (θ + arg L)/π = N(T) + ½ exactly at all six heights, T ≈ 40, 80, 120, 160, 200 and 240. The found counts are 34, 85, 143, 206, 272 and 340, so m₀ = 1. Two of those heights, 160 and 200, fall inside chunks away from their ends, and the cumulative counts there check the boundaries at 140 and 182.
- **Central zero.** Z(2e-4)/Z(1e-4) = 2.000000038, and Z′(0) = 0.30599978 = L′(½).

**Gate** (`check --gate7` → `check9_37a1.json`; `central7` → `central9_37a1.json`; `gate9` → `gate9_37a1.json`). The criteria are P-DL7's, with the central-zero test added as P-DL9 requires: at x = 13 and T = 240, every test needs |Q − Σ| ≤ 2e-4 and (Q − Σ)/tail ∈ [0.8, 1.25].

| test (x = 13) | Q | \|Q − Σ\| | (Q − Σ)/tail |
|---|---|---|---|
| odd k = 2 | 2.677 | 3.2e-6 | 0.971 |
| odd k = 4 | 4.529 | 1.29e-5 | 0.971 |
| odd b₂ − b₃ | 6.628 | 2.02e-5 | 0.972 |
| even b₁ + b₂ | 2.381 | 4.5e-10 | 0.975 |
| even b₂ − b₄ | 7.294 | 7.1e-9 | 0.975 |
| even (1, 1/√2, 0, …), central zero once (F(0)² = 2.56495) | 2.566 | 2.5e-11 | 0.975 |

- **Gate: PASS.** Without the central term the last test is off by 2.56495 = F(0)².
- **x = 40 (reported, not gating).** Here 37a1's bad prime enters. The ratios are 0.96–0.98 and the relative errors at most 3.3e-5.
- **Ordering.** Gate committed in **363919f (05:12:10)**; scan started 05:12:15.

**Scan** (*computed*; `scan9_37a1.json`). The registered command:

```
scan6 --obj 37a1 --v 2,2.5,3,3.5,4 --f1 5 --f2 9 --nmin 24 --acb-max 120 --method inv+full
```

- **Precision.** Every λ is stable (dps 72–105). Recomputed at dps + 80 (`precision9_37a1.json`), the values agree to ≤ 1.8e-15.
- **Convergence.** |ln(λ_5k/λ_9k)| ≤ 0.068.
- **Eigensolver.** There were no full-eigensolver fallbacks; λ is well below 1 throughout.
- **PSD.** FLINT's full eigensolver agrees to ≤ 5.4e-30 at the v = 2 points.

| v | x | N_basis (5k′, 9k′) | λ_even | λ_odd | ln(λ_odd/λ_even) |
|---|---|---|---|---|---|
| 2 | 148 | 107, 187 | 8.148e-14 | 1.805e-16 | −6.11 |
| 2.5 | 231.25 | 144, 252 | 5.955e-19 | 7.351e-22 | −6.70 |
| 3 | 333 | 182, 321 | 4.106e-24 | 4.168e-27 | −6.89 |
| 3.5 | 453.25 | 222, 393 | 3.843e-29 | 1.698e-32 | −7.72 |
| 4 | 592 | 263, 467 | 2.269e-34 | 1.174e-37 | −7.57 |

**P-DL9 for 37a1** (`fit9 --objects 37a1` → `pdl9_fit_37a1.json`):

| | even | odd |
|---|---|---|
| P-DL9a: log-corrected α/4π (9k) | **1.889** (γ = 0.34, resid 0.15) | **1.850** (γ = −3.37, resid 0.14) |
| linear α/4π | 1.880 | 1.943 |
| 5k values | 1.900 / 1.880 | 1.842 / 1.942 |

- **P-DL9a holds for 37a1:** both sectors lie in [1.7, 2.3], far from the kill range. The log terms are small, and unlike the high-weight forms the log-corrected and linear slopes agree.
- **P-DL9b holds for 37a1:** sign(ln(λ_odd/λ_even)) = −1 = ε at all 5 grid points. The ratios run from −6.1 to −7.7, and both λ are below 1e-6 at every point, so a mismatch would have been decisive.

**Same conductor, opposite sign: a controlled comparison** (post hoc). 37b1 (N = 37, ε = +1, rank 0) on P-DL4's grid had ln(λ_odd/λ_even) = +7.09, +7.83, +8.27, +8.54, +8.81 at x/N = 2, 4, 6, 8, 12. 37a1 (N = 37, ε = −1, rank 1) has −6.11 at x/N = 4 and −7.72 at x/N = 12.25. The conductor and gamma factor are the same; the sign follows ε.

**All seven objects.** `fit9` takes the list, `fit9 --objects delta,f16,f18,f20,f22,f26,37a1 --json data/connes/gl2/pdl9_fit.json`. It reads `scan9_<obj>.json` and applies both registered verdicts across all of them. It was tested on synthetic data:
- planted α/4π = 2 and γ = 2.5 are recovered exactly;
- a planted wrong-sign object is killed by 9b.

## 25. Fragile points (P-DL9, 37a1)

- **The split zero run.** It relies on the merge. The counts are exact at heights inside every chunk, so no zero was lost or duplicated. The gate's zero sum uses the merged list.
- **ε is entered as data.** The AFE confirms it: the wrong sign fails by 0.018–1.45. The smallest margin, at a near-real point, is still 10³⁹ times the working accuracy.
- **One object, five points.** P-DL9b on 37a1 is a single new ε = −1 object, but its margin is large: |ln ratio| ≥ 6.1 at every point. P-DL9a's fit has two degrees of freedom and v ≤ 4 (x ≤ 592), and the values are upper bounds on the periodic basis, as everywhere here.

---

# P-DL11 (agent's part): validation of the eta-product newforms g4, g6, g8

P-DL11 is registered at the end of `docs/DECAY_LAW.md` (3778445, 05:24:28). This part does only the validation half. The main session runs the scans; no `scan11_*` file is created here.

## 26. Objects, verification and root numbers (*validated*)

`ETA_NEWFORMS` in `gl2_form_mp.py` holds the three objects:

| object | form | weight k | level N | μ |
|---|---|---|---|---|
| g4 | η(τ)⁴η(5τ)⁴ | 4 | 5 | 3/2, 5/2 |
| g6 | η(τ)⁶η(3τ)⁶ | 6 | 3 | 5/2, 7/2 |
| g8 | η(τ)⁸η(2τ)⁸ | 8 | 2 | 7/2, 9/2 |

- The normalisation is a_n n^{−(k−1)/2}, and the conductor term is log N.
- The q-expansion comes from Euler's pentagonal series raised to the powers with FLINT. The naive product, multiplied out term by term, is the cross-check.
- The objects use the curves' Euler-term machinery. At the bad prime, α_p = a_p/p^{(k−1)/2} and β_p = 0, so a_{p^r} = a_p^r. The resulting entries match the hand values α_p log p/√p: −0.3219 for g4 (p = 5), +0.3662 for g6 (p = 3), −0.3466 for g8 (p = 2).
- `scan6 --obj g4|g6|g8` works with the registered flags. At v = 6 the bases are 319/568, 288/513 and 264/469.

`gl2_form_mp.py objects11` writes `objects11.json`. All three objects pass every check:

| check | g4 | g6 | g8 |
|---|---|---|---|
| a₁…a₈ | 1, −4, 2, 8, −5, −8, 6, 0 | 1, −6, 9, 4, 6, −54, −40, 168 | 1, −8, 12, 64, −210, −96, 1016, −512 |
| fast vs naive product, n ≤ 300 | 0 mismatches | 0 | 0 |
| a_mn = a_m a_n, all 3406 coprime pairs with mn ≤ 2000 | 0 failures | 0 | 0 |
| a_{p²} = a_p² − p^{k−1}, good p with p² ≤ 2000 | 0 failures | 0 | 0 |
| a_{N^r} = a_N^r, N^r ≤ 2000 | yes (4 powers) | yes (6) | yes (10) |
| Hecke rebuild from a_p vs q-expansion, n ≤ 2000 | 0 mismatches | 0 | 0 |
| Deligne \|a_p\| ≤ 2p^{(k−1)/2}, good p ≤ 2000 | yes | yes | yes |
| \|a_N\| = N^{(k−2)/2} | a₅ = −5 | a₃ = +9 | a₂ = −8 |
| Euler terms vs log-derivative | 2.4e-31 | 2.3e-41 | 1.4e-41 |

**Root numbers.** The AFE was run at five random w with **both** signs:

| | g4 | g6 | g8 |
|---|---|---|---|
| ε = +1: split rel / FE rel | ≤ 3.9e-41 / ≤ 2.6e-41 | ≤ 5.9e-42 / ≤ 2.9e-42 | ≤ 1.4e-42 / ≤ 8.7e-43 |
| ε = −1: split rel / FE rel | 0.98–1.8 / 1.2–4.2 | 0.66–1.9 / 0.5–6.1 | 0.74–2.4 / 1.1–3.7 |
| **ε (AFE)** | **+1** | **+1** | **+1** |
| Atkin–Lehner prediction i^k·w_N, with w_N = −a_N/N^{k/2−1} | i⁴·(+1) = +1 | i⁶·(−1) = +1 | i⁸·(+1) = +1 |
| AFE vs Dirichlet series | 1.3e-20, 2.4e-18 | 9.6e-21, 1.7e-18 | 3.6e-20, 6.4e-18 |
| L(½), analytic (t = 1 and 1.3 agree) | 0.411861328386 | 0.560038691049 | 0.681292530966 |

All three have ε = +1 and L(½) ≠ 0, so none has a central zero. The ε = −1 cases on P-DL11's grid therefore come from 37a1 alone. Before the AFE was run the registry held ε = None, and `eps_of` refused to return a value; ε = +1 was entered only after the test above.

## 27. Zeros and gates (*validated*)

**Zeros to T = 240.** g4 ran as two τ-chunks, [0, 182] and [182, 240], joined by `zeros-merge`; g6 and g8 ran whole. In every case (θ + arg L)/π is an integer equal to the number of zeros found, at all six heights. That integer value is what ε = +1 with L(½) ≠ 0 requires (m₀ = 0).

| object | zeros in (0, 240] | N(T) at T ≈ 40 / 80 / 120 / 160 / 200 / 240 | first zeros |
|---|---|---|---|
| g4 | 264 | 22 / 60 / 106 / 155 / 209 / 264 | 7.803685993, 9.41501646 |
| g6 | 245 | 19 / 54 / 97 / 143 / 193 / 245 | 8.20343005, 10.16731596 |
| g8 | 230 | 17 / 50 / 89 / 133 / 181 / 230 | 8.27204092, 11.39598699 |

**Gates** (`check --gate7` writes `check11_*.json`; `central7` writes `central11_*.json`; `gate9` writes `gate11_*.json`). The criteria are P-DL7's: at x = 13 and T = 240, every test needs |Q − Σ| ≤ 2e-4 and (Q − Σ)/tail ∈ [0.8, 1.25]. The b₀ test, c = (1, 1/√2, 0, …) with F(0)² = 2.565, is required only when ε = −1. It was run anyway as a sixth test; with ε = +1 it must leave only the tail, with no central term.

| object | max \|Q − Σ\| (5 tests) | (Q − Σ)/tail | b₀ test: Q − Σ (ratio) | P-DL7 gate | 6-test gate |
|---|---|---|---|---|---|
| g4 | 1.61e-5 | 0.931–0.936 | 1.9e-11 (0.931) | PASS | PASS |
| g6 | 1.76e-5 | 1.09 | 2.1e-11 (1.09) | PASS | PASS |
| g8 | 1.57e-5 | 1.01 | 1.9e-11 (1.01) | PASS | PASS |

The bad prime of each form (5, 3, 2) is below 13, so it enters every gate test.

## 28. `fit11` (written for the main session's scans; tested on synthetic data)

`fit11 --objects g4,g6,g8,11a1,37a1 --json …` reads `scan11_<obj>.json` and applies P-DL11 as registered:
- **11a (rate).** The log-corrected α/4π against [1.7, 2.3], with kill outside [1.4, 2.6]. The linear slope and the 5k values are reported alongside.
- **11b (index).** The continuous index n* is the n at which the least-squares slope of ln(λ/ℓ_n(4πv)) against ln(4πv) crosses zero, interpolated linearly between integers n = 0..80. The comparison is n* − k against δ(ε, s) = 0.7, 3.1, 3.8, 1.3 for (+1, even), (+1, odd), (−1, even), (−1, odd). The result holds if at least 8 of the 10 cases are within ±1.2, and is killed if more than 2 are off by more than 2.
- **11c (ε rule).** sign ln(λ_odd/λ_even) = ε at every point; a mismatch kills only where both λ < 1e-6.
- **The weights and root numbers used.** k = 4, 6, 8, 2, 2. ε = +1 for g4, g6, g8 (AFE, §26) and for 11a1 (§2); ε = −1 for 37a1 (§23).
- **A cross-check, not the method.** ln ℓ_n(c) = const_n + (n + ½) ln c − 2c, so the slope is exactly S − n − ½, where S is the slope of ln λ + 2c against ln c. The interpolated n* must therefore equal S − ½. `fit11` prints both.

**Synthetic tests** (scratchpad only). The planted data are ln λ = const + (m + ½) ln c − 2c, with m = k + δ + offset.
- **Scenario A:** offsets of 0 to ±1.1 in 8 cases and ±1.5 in 2, plus a wobble of 0.05 sin v. n* comes back at m − 0.019 (the wobble's slope) in every case, and equals S − ½. α/4π = 2.02. The verdicts are 11a HOLDS, 11b HOLDS (8/10 within ±1.2, 0 off by more than 2) and 11c HOLDS (25/25).
- **Scenario B:** offsets of 2.5 in 3 even cases, and 37a1 planted with the wrong sign. n* comes back exactly. The verdicts are 11b KILLED and 11c KILLED.
  - The wrong-sign planting also made 37a1's odd index equal its even one, a fourth case off by 2.5, so 4 cases are off by more than 2.
  - The planted even offsets of 3.2 exceed the odd 3.1 for g4 and g8, so those sign mismatches are reported too.
  
  Both effects come from the planting and are correct outputs.

## 29. Fragile points (P-DL11 validation)

- **No ε = −1 form among the new objects.** All three are ε = +1, so P-DL11b's (−1, ·) cells and P-DL11c's negative sign rest on 37a1 alone.
- **The structure of the index metric.** The continuous index is exactly S − ½, a smoothed local log-slope of λe^{2c} in c. It inherits any curvature left on v ∈ [6, 10], and the registered δ values were themselves estimated from P-DL9 data with the same metric.
- **Validation covers the form, not the minimum.** The gates check the form against the zero sum at x = 13. The scans at v = 6–10 (x = 72–500 for g4–g8) use the same code, but nothing here validates them independently. The remaining checks there are precision stability, basis convergence and the PSD check.

---

# P-DL12 (agent's part): validation of four fresh curves

P-DL12 is registered at the end of `docs/DECAY_LAW.md` (b22392c, 07:10:46). This part does only the validation half. The main session runs the scans; no `scan12_*` file is created here.

## 30. Curves, root numbers and ranks (*validated*)

The Weierstrass models are from memory. Each was checked by Tate's criteria (`objects6`, `objects12.json`) and by the AFE with both signs (`afe-both`, `afe_both12.json`). The registry kept ε = None, and `eps_of` refused to return a value, until the AFE had decided it. The registry labels were all kept; no substitute was needed.

| | 14a1 | 17a1 | 43a1 | 53a1 |
|---|---|---|---|---|
| model [a1, a2, a3, a4, a6] | [1, 0, 1, 4, −6] | [1, −1, 1, −1, −14] | [0, 1, 1, 0, 0] | [1, −1, 1, 0, 0] |
| Δ | −21952 = −2⁶·7³ | −83521 = −17⁴ | −43 | −53 |
| c4 (prime to N) | −215 | 33 | 16 | −15 |
| conductor (Tate) | 14 | 17 | 43 | 53 |
| a_p at p \| N | a₂ = −1 (non-split), a₇ = +1 (split) | a₁₇ = +1 (split) | a₄₃ = −1 (non-split) | a₅₃ = −1 (non-split) |
| a₂, a₃, a₅, a₇ | −1, −2, 0, 1 | −1, 0, −2, 4 | −2, −2, −4, 0 | −1, −3, 0, −4 |
| gcd of #Ẽ(F_p), good p ≤ 2000 | 6 (torsion Z/6) | 4 (Z/4) | 1 | 1 |
| independent q-expansion | η(τ)η(2τ)η(7τ)η(14τ): a_n agree for n ≤ 2000 | — | — | — |
| ε = +1: split / FE rel | ≤ 2.7e-41 / ≤ 4.0e-41 | ≤ 1.2e-42 / ≤ 2.1e-42 | 0.069–9.9 / 0.040–3.4 | 2.5–8.6 / 0.98–2.4 |
| ε = −1: split / FE rel | 0.99–2.0 / 1.3–11.5 | 0.41–2.2 / 0.51–3.8 | ≤ 3.7e-42 / ≤ 4.4e-42 | ≤ 1.6e-43 / ≤ 2.2e-43 |
| **ε (AFE)** | **+1** | **+1** | **−1** | **−1** |
| Atkin–Lehner prediction ε = −w_N, w_p = −a_p | +1 | +1 | −1 | −1 |
| L(½), analytic (t = 1 and 1.3) | 0.330223659344 | 0.386769938388 | 0, forced (−1.6e-42) | 0, forced (−1.3e-42) |
| L′(½) (t = 1 and 1.3 agree to 15 digits) | — | — | 0.343523974618 | 0.435863824178 |
| **analytic rank** | **0** | **0** | **1** | **1** |

- **Hasse and Euler terms.** Hasse's bound holds for all four. The Euler terms match the Dirichlet log-derivative to ≤ 3.9e-31.
- **Normalisation.** The AFE agrees with the Dirichlet series to its truncation level, ≤ 1.8e-17.
- **The smallest failure of a wrong sign** is 0.040 (43a1, ε = +1, at a near-real point). That is still 10³⁹ times the agreement of the right sign.
- **Scan objects.** The four curves are `scan6` objects through `CURVES`. With the registered flags, v = 6 gives bases of 381/680, 393/701, 448/801 and 461/823. At v = 10 the 9k′ basis reaches 1311–1551, about as large as 37a1's on the P-DL11 grid.

## 31. Zeros and gates (*validated*)

**Zeros to T = 240.** The τ-chunks ran under a 4-worker scheduler: [0, 182] and [182, 240] for 14a1 and 17a1, and four cost-balanced chunks for 43a1 and 53a1. `zeros-merge` joined them and ran the counts at six heights.

| curve | ε | zeros in (0, 240] | first zeros | (θ + arg L)/π − found | central zero |
|---|---|---|---|---|---|
| 14a1 | +1 | 303 | 5.579286817, 7.575711000 | 0 at all 6 heights (found 28 / 73 / 125 / 181 / 241 / 303) | none (L(½) ≠ 0) |
| 17a1 | +1 | 310 | 4.741993155, 7.819103955 | 0 at all 6 heights (found 29 / 76 / 129 / 186 / 247 / 310) | none (L(½) ≠ 0) |
| 43a1 | -1 | 345 | 4.494720273, 6.828717445 | ½, so m₀ = 1 at all 6 heights (found 35 / 87 / 146 / 209 / 276 / 345) | simple: Z(2e-4)/Z(1e-4) = 2.000; Z′(0) = 0.34352 = L′(½) |
| 53a1 | -1 | 353 | 4.508628350, 6.043478899 | ½, so m₀ = 1 at all 6 heights (found 36 / 89 / 151 / 215 / 283 / 353) | simple: Z(2e-4)/Z(1e-4) = 2.000; Z′(0) = 0.43586 = L′(½) |

The merges of 17a1 and 43a1 each dropped one zero that both neighbouring chunks had found at a boundary (a chunk can overrun its end by one step); the exact counts at heights inside every chunk confirm nothing was lost or double-counted.

**Gates** (`check --x 13,60 --gate7` writes `check12_*.json`; `central7` writes `central12_*.json`; `gate9` writes `gate12_*.json`). The criteria are P-DL7's: at x = 13 and T = 240, every test needs |Q − Σ| ≤ 2e-4 and (Q − Σ)/tail ∈ [0.8, 1.25]. The central-zero test (b₀ test, F(0)² = 2.565, central zero counted once) is required for ε = −1. It was run for all four; for ε = +1 it must leave only the tail.

| curve | ε | max \|Q − Σ\| (5 tests, x = 13) | (Q − Σ)/tail | b₀ test: Q − Σ (ratio) | x = 60 ratios (not gating) | 6-test gate |
|---|---|---|---|---|---|---|
| 14a1 | +1 | 1.76e-05 | 0.920–0.923 | 2.1e-11 (0.920) | 1.012–1.018 | PASS |
| 17a1 | +1 | 1.87e-05 | 0.968–0.970 | 2.3e-11 (0.970) | 0.994–1.014 | PASS |
| 43a1 | -1 | 1.92e-05 | 0.907–0.913 | 2.3e-11 (0.907) | 0.992–0.997 | PASS |
| 53a1 | -1 | 2.04e-05 | 0.950–0.951 | 2.5e-11 (0.950) | 1.006–1.019 | PASS |

The x = 60 runs are reported but do not gate; there the bad primes 17, 43 and 53 enter as well.

## 32. `fit12` (written for the main session's scans; tested on synthetic data)

`fit12 --objects 14a1,17a1,43a1,53a1 --json …` reads `scan12_<obj>.json` and uses the P-DL11 metric (`continuous_index`, c = 4πv):
- **12a (index).** n* − k is compared with the hypothesis: −0.64 (even) and +1.81 (odd) for ε = +1; +1.32 (even) and −0.67 (odd) for ε = −1. It holds if at least 6 of 8 cases are within ±0.5, and is killed if more than 2 are off by more than 1.0. Otherwise the verdict is "not killed, partly missed".
- **12b (rate).** The log-corrected α/4π against [1.7, 2.3], with kill outside [1.4, 2.6]. The linear slope and the 5k values are reported alongside.
- **12c (ε rule).** sign ln(λ_odd/λ_even) = ε at every point; a mismatch kills only where both λ < 1e-6.
- **ε values used:** +1 for 14a1 and 17a1, −1 for 43a1 and 53a1 (§30).
- **Extrapolation.** If the slope's zero crossing lies outside n = 0..80, `fit12` extrapolates it linearly (the slope is exactly linear in n) and flags it. `fit11` is unchanged.

**Synthetic tests** (scratchpad only). The planted data are ln λ = −10 + (m + ½) ln c − 2c, with m = k + hypothesis + offset.
- **A:** offsets of 0 to ±0.45 in 6 cases and ±0.8 in 2, plus a wobble of 0.03 cos v. Every n* comes back shifted by the wobble's slope (−0.124), and equals S − ½. The verdicts are 12a HOLDS (6/8 within ±0.5, 0 off by more than 1.0), 12b HOLDS (α/4π = 2.000) and 12c HOLDS (20/20).
- **B:** offsets of ±1.5 in 3 cases, one of which puts n* at −0.14, below 0. n* comes back exactly; the negative one is extrapolated and flagged, and equals S − ½. The verdict is 12a KILLED (5/8 within, 3 off by more than 1.0).
- **C:** A with 53a1 given the wrong sector order. The verdicts are 12c KILLED (15/20) and 12a NOT KILLED, partly missed, because the flipped sign also moved 53a1's odd index. That exercises the third verdict branch.

## 33. Fragile points (P-DL12 validation)

- **The models come from memory.** Each is verified by its conductor, ε and rank, which is what the test needs. Only 14a1 has an independent q-expansion (its eta product).
- **The rank-1 curves** have their central zero verified three ways:
  - the half-integral count, (θ + arg L)/π = N(T) + ½ at all six heights;
  - the Z ratio at the centre;
  - the b₀ test.

  But L′(½) ≠ 0 is shown only numerically, at about 0.34–0.44.
- **The large bases at v = 10** (up to 1551) are new territory for the eigensolver, in the main session's scans.
