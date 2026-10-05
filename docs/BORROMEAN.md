# Borromean primes and the Rédei L-function (issue the private research log (issue 7))

**Status:** step 1 pre-registered (c595803) and run: **all gates pass**. Step 2 pre-registered (e94635f), first scan void (float64 inputs, fixed in 8422f34), rerun: **R-a, R-b and R-c all hold**.

## Background (from the issue and its review)

- **The objects.** For primes p₁, p₂ ≡ 1 mod 4 with (p₁/p₂) = 1, let R be the Rédei field: the unique D₄ octic containing Q(√p₁, √p₂) and ramified only at p₁ and p₂. Let ρ be the 2-dimensional irreducible representation of Gal(R/Q).
- **The trace rule** [D]. a_p = tr ρ(Frob_p) is 0 when p is linked to p₁ or p₂ (a Legendre symbol is −1). When p is pairwise unlinked, a_p = 2·[p₁, p₂, p], where [p₁, p₂, p] = +1 iff p splits completely in R.
- **Arithmetic topology** [L]. [p₁, p₂, p₃] = (−1)^{μ₂(123)}, the mod-2 Milnor triple linking number (Morishita). A Borromean triple is pairwise unlinked with symbol −1.
- **The L-function** [D, N]. ρ is even (det ρ = χ_{p₁p₂}), so L(s, ρ) is of Maass type: conductor p₁p₂ and gamma factor Γ_R(s)² or Γ_R(s+1)². It is **not** a weight-1 form.
- **What step 1 can and cannot show.** By the explicit formula, the Landau readout of a_p from the zeros is decided in advance. Step 1 is therefore a pipeline check with gates only, and it has no scientific kill.
- **Do not use.** The repository's older Rédei code and `data/borromean_triples.db` are wrong, and are not used here.

## Step 1 pre-registration (gates only)

**Code.** `scripts/borromean/redei_lib.py` and `scripts/borromean/step1_gates.py` (`--stages g0,proxy,chars,rho,g2,g3`). Output: `data/borromean/step1.json`.

**Environment.** Python 3.14 `.venv`, with cypari2 2.2.4 bundling PARI 2.17.2. Homebrew `gp` 2.19.0 is installed but not used by these scripts.

**Objects.**
- **(p₁, p₂) = (5, 29).** R = polredabs of the D₄ octic built from the solution (11, 1, 2) of x² = 5y² + 29z², with twist 2. ρ comes from `lfunartin` on the unique character of degree 2.
- **Controls.** L(χ₅) and L(χ₂₉), whose product has the same conductor 145 and gamma factor Γ_R(s)². The proxy is L(χ₁₂)L(χ₁₃): N = 156, Γ_R(s)².

**G0 (inputs; must all hold).**
- Gal(R/Q) is [8, 3] = D₄, |disc R| = 145⁴, and the signature is [8, 0].
- L(ρ) has N = 145, Vga = [0, 0], root number +1, and `lfuncheckfeq` ≤ −100 bits.
- The order at s = ½ is m₀ = 0.
- a_p = 2·[5, 29, p] for every unramified p ≤ 10⁴, with a_p = 0 for linked p.
- Rédei reciprocity holds on 25 random admissible triples below 400 (seed 7).
- [13, 61, 937] = −1 in all three orderings.

**Zeros.**
- L(ρ) to T = 1000 in chunks of 100 with divz 32, at 38 and at 57 digits. All chunks share one critical-line `lfuninit` up to T + 50.
- L(χ₅) and L(χ₂₉) to T = 1000 in the same chunks, at 38 digits.

**G1 (completeness).** This deviates from the review on purpose.
- **Why.** The review asked for an argument-principle count in every chunk up to T = 1000. Each off-line evaluation of L at height about 1000 takes more than five minutes in PARI 2.17.2: one evaluation of L(1.5 + 1000i) for the proxy did not finish in over five minutes. So G1 has two parts.
- **(a) Argument principle at T = 100 and 200,** where it is fast.
  - N(T) = [θ(T) + arg L(½+iT)]/π − m₀/2, with arg L continuous along σ from 2 down to ½ at height T. At σ = 2 the principal value is correct, because |Im log L| ≤ 2 log ζ(2) < π.
  - θ(T) = (T/2) log N + Σ_j [−(T/2) log π + Im log Γ((½+μ_j+iT)/2)].
  - Each N(T) must lie within 1e-6 of an integer and equal the number of zeros found below T.
- **(b) Turing's method at every boundary b = 100, 200, …, 1000.**
  - D(b) is the mean over t ∈ [b − 20, b], on a grid of step 0.02, of θ(t)/π − N_found(t). It must satisfy |D(b)| < 0.5.
  - D(b) ≈ (number of missed zeros below b) − mean S(t), and S has mean near 0, so one missed zero shifts D by about 1.
- **Repair (registered).** A chunk (a, b] whose Turing value jumps, D(b) − D(a) > 0.5, is searched again with divz 128 and then 512, and the union is kept. Every repair is logged in the output. Sign-change searches step over close pairs of zeros. G1 is judged after repair. The positive control is never repaired.
- **Agreement.** The 38- and 57-digit zero lists of L(ρ) must agree.
- **Positive control (must fire).** The proxy at 19 digits, with one `lfunzeros(L, 1000)` call at the default divz, must fail G1. The proxy at the good settings (38 digits, divz 32) must pass.

**Statistics.** S(x) = Σ_{0<γ≤1000} (1 − γ/1000) cos(γ log x), and R(x) = −S(x) / [(1000/4π) log x/√x].

**G2 (explicit formula).** |S(x) − P(x)| ≤ τ at 2000 log-spaced x in [2, 1000].
- P(x) = (1/2π)∫₀^T (1 − t/T) cos(t log x) ρ_dens(t) dt − (T/4π) Σ_{n ≤ 4·10⁵} Λ(n)/√n [F(log n − log x) + F(log n + log x)].
- F(v) = (sin(Tv/2)/(Tv/2))².
- ρ_dens(t) = log N + Σ_j [−log π + Re ψ((½ + μ_j + it)/2)].
- Λ comes from the a_n by the prime-power recursion.
- τ = 3 × the proxy's maximum residual, computed earlier in the same run, before L(ρ)'s G2.

**G3 (readout).**
- |R(p) − a_p| ≤ 0.3 for every unlinked p ≤ 1000 other than 5 and 29. The review forecasts at most 0.08, from the a_n alone.
- **Negative control.** On the zeros of L(χ₅)L(χ₂₉) (the union of the two lists), R(p) lies within 0.3 of a_p at no more than 2 of the 18 primes with a_p = −2. There the product's coefficient is +2.

**Any gate failure** means the pipeline must be fixed before step 2. There is no scientific kill in step 1.

**Pre-registration development.** These calculations were run before this commit, and only on non-outcome objects:
- G0 on (5, 29), which passed (inputs only);
- the counting, chunking and explicit-formula code on the proxy alone, at T = 200;
- a Turing-window test of G1(b) on the proxy to T = 1000, at both good and default settings. The default settings missed 73 zeros, and D(b) rose from 0.04 to 73.9, so the check fires. The good settings gave |D| ≤ 0.04 at b = 100–900 but D(1000) = 4.7: PARI loses zeros near the top of an `lfuninit` domain. The zero searches therefore share one critical-line `lfuninit` up to T + 50 and keep only zeros ≤ T. This was re-tested on the proxy against the separately computed zeros of L(χ₁₂) and L(χ₁₃), which number 1043 + 1056 = 2099 to T = 1000.
  - **Margin.** With the margin the product search found 2097, and D jumped by exactly 2 between 500 and 600.
  - **The missed pair.** The two missed zeros are a close pair from different factors, 504.498119 and 504.502780, 0.0047 apart.
  - **Repair.** divz 128 recovered both (223 of 223 in that chunk, 145 s), and divz 512 confirmed it. This is the repair rule above;
- timing attempts for off-line evaluation at T = 1000, which did not finish;
- the positive control at T = 200, which fired: 316 of 318 zeros found.

No zero of L(s, ρ), L(χ₅) or L(χ₂₉) had been computed.

## Step 1 results: all gates pass

**Run.** `data/borromean/step1.json`, `step1_zeros.json` and `step1.run.txt`, on M3 in 44 minutes, from commit 0517f62. That commit only adds the zero-list output.

| gate | result |
|---|---|
| G0 (inputs) | pass: [8, 3] = D₄, \|disc\| = 145⁴, signature [8, 0], N = 145, Vga [0, 0], W = +1, feq −131 bits, m₀ = 0; trace rule 0 mismatches on 1227 primes; reciprocity 25/25; [13, 61, 937] = −1 in all orderings |
| G1, positive control | fired: the proxy at 19 digits with one call found 2005 zeros, and max\|D\| = 73.9 |
| G1, proxy (good settings) | pass: 2099 zeros, which equals the separate count for L(χ₁₂) and L(χ₁₃); one logged repair, chunk (500, 600], divz 128, +2 zeros (the close pair at 504.50); max\|D\| = 0.035 |
| G1, L(χ₅) and L(χ₂₉) | pass: 904 and 1183 zeros; argument principle exact at 100 and 200; max\|D\| ≤ 0.028; no repairs |
| G1, L(ρ) | pass at 38 and at 57 digits: 2087 zeros each, with identical lists; the argument principle gives 135 and 315 at T = 100 and 200, integral to 0; max\|D\| = 0.011; no repairs |
| G2 (explicit formula) | pass: max \|S − P\| = 4.6e-4 ≤ τ = 1.4e-3 (the proxy's maximum was 4.7e-4) |
| G3 (readout) | pass: max \|R(p) − a_p\| = 0.082 over the 34 unlinked p ≤ 1000 (forecast ≤ 0.08). Negative control: 0 of 18 agreements |

**The Borromean partners of (5, 29) below 1000**, read from the zeros of L(ρ):

| p₃ | 181 | 241 | 349 | 401 | 661 | 701 | 761 |
|---|---|---|---|---|---|---|---|
| R(p₃) | −1.969 | −2.082 | −1.995 | −1.997 | −1.997 | −2.003 | −1.998 |

**What this shows.** These 7 primes are pairwise unlinked from 5 and 29, but Borromean with them. On the torus of L(ρ)'s zeros they sit at a = −2: the zeros gather at θ_{p₃} = 0 instead of avoiding it. On the zeros of L(χ₅)L(χ₂₉), the same primes read +2, which is ordinary avoidance. The largest leakage at a linked prime is 1.22, at 823, from its unlinked neighbour 821, as forecast.

**Status of the result.** This is the pipeline check registered as step 1 [N]. It adds no new mathematics: the explicit formula decides it in advance. It certifies the zeros (2087 to T = 1000, complete by Turing's method) and the coefficient tables that step 2 uses.

## Step 2 pre-registration (the matched-pair test)

**Question (H_same).** Does the window minimum depend only on the invariants (d, N, μ, ε, m₀)? Or does it distinguish a primitive L-function from a product with the same invariants?

**Objects.** Both have degree 2, N = 145, μ = (0, 0), so Γ_R(s)², ε = +1, m₀ = 0, and no pole.
- **A = L(s, ρ).** The primitive Maass-type L-function of the Rédei field of (5, 29).
- **B = L(s, χ₅)L(s, χ₂₉).** The product with the same invariants.

**Code.** `scripts/borromean/step2_matched.py`. The `prep` stage uses PARI and writes the Λ(p^k) tables and the zero files from step 1. The `check`, `scan` and `fit` stages use only mpmath and FLINT.

**Validation already done.** These are not outcomes; no eigenvalue was computed.
- B's degree-2 form (`gl2_form_mp.form` with μ = (0, 0), log 145 and the product's prime terms) equals the sum of our established degree-1 forms for χ₅ and χ₂₉ to 1e-40, at working precision, at x = 13, 60 and 300 in both sectors.
- A's form differs from B's only in its prime terms. Those come from PARI's a_n by the prime-power recursion, which step 1's G2 tested.

**Gate (before any scan).** `check` compares Q(c) = cᵀMc with 2Σ_{0<γ≤1000} F_c(γ)² plus the Maass-type tail, for the five edge-vanishing test functions of `gl2_form_mp check`.
- The criteria are P-DL7's, at x = 13: |Q − Σ| ≤ 2e-4 and (Q − Σ)/tail ∈ [0.8, 1.25] for every test.
- x = 60 is reported but does not gate.
- The zeros come from step 1, after its G1.
- A failed gate stops step 2.

**Scan (P-DL12 protocol).**
- v = √(x/145) ∈ {6, 7, 8, 9, 10}, in both sectors, for both objects.
- N_basis = ⌊f·k⌋ + 8 with f = 5 and 9, k = 2vL.
- Every λ must be precision-stable (`stable_lams`, inverse iteration).
- Convergence is reported as ln(λ_5k/λ_9k). The run is inconclusive, not killed, if |ln(λ_5k/λ_9k)| > 0.3 at v = 10.

**Calibration of R-c, from existing data (no new outcome).** The P-DL12 curves come in same-class pairs: 14a1 with 17a1 (ε = +1), and 43a1 with 53a1 (ε = −1). Each pair shares μ, ε and m₀ but not N.
- At equal v, their |Δ ln λ| is at most 0.69 over 20 point comparisons, with a median of about 0.2.
- Their |Δn*| is at most 0.26 over 4 comparisons.
- A and B share N as well, so H_same predicts agreement at least this close.

**Predictions.**

| prediction | holds | killed |
|---|---|---|
| **R-a** (rate of A, both sectors; log-corrected fit ln λ = −αv + γ ln v + β) | α/4π ∈ [1.7, 2.3] | α/4π outside [1.4, 2.6] in either sector |
| **R-b** (sector order: m₀ = 0, so the central-multiplicity conjecture predicts λ_odd > λ_even at every point, for A and for B) | no mismatch | any mismatch with both λ < 1e-6 |
| **R-c** (H_same) | \|n*_A − n*_B\| ≤ 0.5 in both sectors, and \|ln λ_A − ln λ_B\| ≤ 1.5 at every (v, sector) | \|n*_A − n*_B\| > 1.0 in either sector, or \|ln λ_A − ln λ_B\| > 3.0 at 3 or more of the 10 points |

- **Outcomes outside both columns** are reported as "partly missed".
- **The R-c bands.** The review proposed 3 and 6 nats. Here they are 1.5 and 3, about 2× and 4× the calibration maximum.
- **Index.** n* is P-DL12's continuous index (`gl2_form_mp.continuous_index`, c = 4πv), on the 9k basis.
- **No index prediction for μ = (0, 0).** P-DL12's weight rule covers holomorphic weights only, so n* itself is reported without one.

**Hosts.** ml, in the pinned Docker image. The code and `data/borromean/step2_objects.json` are copied in with rsync, and results are copied back with scp. There are 20 points; the 9k bases reach about 1,730 at v = 10.

**Gate results (before any scan): PASS, 10 of 10 at x = 13** (`data/borromean/step2_check.json`).
- For both objects and all five test functions, |Q − Σ| ≤ 3.9e-7 against the limit of 2e-4.
- (Q − Σ)/tail lies in 0.9988–1.002, so Q equals the zero sum up to the predicted tail.
- The first launch was stopped by the session's 30-minute background limit after L(ρ)'s five tests had all passed (`step2_check_partial_killed.run.txt`). The gate was then run again in full.
- The x = 60 checks, which do not gate, are running separately and will be reported in `step2_check_x60.json`.

**Incident: the first scan was void, and the fix was added before the rerun.**
- **What happened.** The first scan (from e94635f) was void: `prep` stored Λ(p^k) as float64, which put a noise floor of about 1e-17 under eigenvalues near e^{−150}. Every value came out around 1e-18 with random signs, some negative.
- **Why the gates missed it.** The pre-scan validation built B's prime terms inline at full precision, not through the stored tables. The zero-sum gate at x = 13 works at a scale of 1e-7, so it could not see an error of 1e-17.
- **The fix.** Integer coefficients, with Λ(p^k) = c·log p formed at working precision.
- **New gates on the production path** (`data/borromean/step2_validate.json`), all passed:
  - V1, ρ's coefficients satisfy the Euler-factor recursion: passed on all 1755 prime powers;
  - V2, B's production-path form equals χ₅ + χ₂₉ to ≤ 1e-100 at 120 digits: it agreed to 3e-120 at x = 13, 1000 and 14500;
  - V2's positive control, a float64 copy of the terms, must fail: it was off by 4.6e-16, so it fired.
- **Unchanged.** The predictions and statistics stay as registered. The void files are kept in `data/borromean/void_float_lambda/`.

## Step 2 results: R-a, R-b and R-c all hold

**Run.** 20 points on ml, in the pinned Docker image, from code at 8422f34, with checksums matching. The per-point files and logs are in `data/borromean/parts/`, merged into `scan_step2_{rho,chi5chi29}.json`. The fit is `data/borromean/step2_fit.json`.
- There were no errors, no negative eigenvalues and no inverse-iteration failures.
- Convergence: |ln(λ_5k/λ_9k)| ≤ 0.10 at every point, within the 0.3 limit.

| prediction | result | verdict |
|---|---|---|
| R-a: rate of A (log-corrected) | α/4π = 2.000 (even, γ = 0.01, max resid 0.05) and 1.996 (odd, γ = 1.83, max resid 0.08); linear 2.000 and 1.978. B, reported: 2.009 and 2.069 | **holds** ([1.7, 2.3]) |
| R-b: sector order | λ_odd > λ_even at all 10 points for A and B; ln(λ_odd/λ_even) = 10.4–11.5 (A) and 9.9–10.8 (B); every λ < 1e-6 | **holds** (no mismatch) |
| R-c: H_same | \|n*_A − n*_B\| = 0.41 (even: −0.50 against −0.09) and 0.15 (odd: 1.70 against 1.55); \|ln λ_A − ln λ_B\| ≤ 0.47 at all 10 points | **holds** (≤ 0.5 and ≤ 1.5) |

**Per-point ln(λ_A/λ_B):**

| sector | v = 6 | 7 | 8 | 9 | 10 |
|---|---|---|---|---|---|
| even | −0.15 | −0.33 | −0.40 | −0.19 | −0.47 |
| odd | +0.29 | −0.04 | +0.01 | +0.32 | +0.23 |

**Reading [E].** On these windows (c = 4πv = 75–126), the window minimum of the primitive Maass-type L(s, ρ) is indistinguishable from that of the product L(χ₅)L(χ₂₉), which has the same invariants (d = 2, N = 145, Γ_R(s)², ε = +1, m₀ = 0).
- **Size of the gap.** The largest gap is 0.47 nats. That is within the spread of same-class pairs in P-DL12 (0.69), even though those pairs also differ in conductor.
- **Rate and sector order.** The rate is the decay law's e^{−2dT*} (d = 2), and the sector order is the one m₀ = 0 predicts.
- **What it says about the instrument.** At this resolution the instrument sees the invariants, not primitivity. The Borromean data a_{p₃} = ±2 that separate A from B, invisible to the window minimum, are visible on the torus of zeros (step 1).
- **Caveats.** These are finite-basis upper bounds, with v ≤ 10 and a single matched pair.

## Step 2 results (detail)

(see `data/borromean/step2_fit.json`)
