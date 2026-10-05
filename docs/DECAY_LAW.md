# The Decay Law of the Window Minimum

## Status at a glance (2026-10-02)

All values are finite-basis Rayleigh–Ritz upper bounds, in multiprecision floating point, over finite ranges. None is certified. The claims manifest (`scripts/decay_claims_check.py`) regenerates every quoted number from the committed data.

| test | object | outcome |
|---|---|---|
| P1 | antiperiodic ("720°") basis converges faster | **failed**: about 1.2 modes of reach, same limit |
| P-DL1 | rate e^{−4πx/q} for 6 real characters | rate confirmed with a log term (α/4π 0.984–1.004); the strict linear version **failed** in the odd sectors |
| P-DL2 | prolate-index R bounded out to x/q = 16 | **holds** (R/R(10) ∈ [0.79, 1.20]; index 12/12) |
| P-DL3 | degree-2 products: worst-factor floor against product zero counting | worst-factor **killed**; zero counting survives |
| P-DL4 | primitive GL(2) (11a1, 37b1): three scalings | **"neither"**: the zero-counting variable is right but the rate is doubled, e^{−4T*} |
| P-DL5 | degree-2 products: e^{−2d·T*} | **holds** (α/4π 1.79–2.29) |
| P-DL6 | fresh 15a1, 19a1, Δ: rate; index n = Σμ + 2s | rate **holds** (2.10–2.17); index **partly missed** (1 of 6 exact) |
| P-DL7 | level-1 forms k = 16–26 on v ≤ 4: rate, index slope, offset | all **killed**: the grid was pre-asymptotic for high weight |
| P-DL8 | fresh characters −8, −11, 12, 13, −15, 17 | **holds**: index 11/12, rate α/4π 0.984–1.026 |
| P-DL9 | level-1 k = 12–26 on v = 6–10, plus fresh 37a1 (ε = −1): rate; sector order = ε | **holds**: α/4π 1.85–2.07; ε rule 35/35 |
| P-DL10 | degree-3 products: e^{−6T*} | **holds** (α/4π 2.74–3.20) |
| P-DL11 | fresh g4, g6, g8 (+11a1, 37a1) on v = 6–10: rate, index n = k + δ(ε, s), ε rule | rate and ε rule **hold**; index **missed** (0/10 within ±1.2; a near-constant level offset) |
| P-DL12 | four fresh curves (14a1, 17a1 ε = +1; 43a1, 53a1 ε = −1): index offsets within ±0.5, rate, ε rule | **holds**: 8/8, α/4π 1.88–2.06, 20/20 |
| R1 | exclusion regions from certified windows (`docs/RESEARCH_R1_EXCLUSION.md`) | exclusion is a witness (consistency) problem; the certified constant plays no role [D]. Under the hypothesis that all other zeros lie on the line, the computed reach of supports up to 4.0 is far below the verified height: a **forecast** [N], not a theorem (the computed statistic is a lower bound; g ≈ 1.2 measured at one support) |
| P-R389 | rank-2 curve 389a1 (ε = +1, m₀ = 2): does the sector order follow ε or the central zero? | **H_m holds** as registered: odd below even at 5/5 points, ln(λ_odd/λ_even) −7.4 to −8.8; the root-number-only prediction (R2's reading) is contradicted in the tested bases. Consistent with central multiplicity controlling sector order (m₀ > 0 against m₀ = 0; the central weight saturates, so the size of m₀ does not enter); a conjecture, not a universal rule. Rate α/4π 1.95 / 1.94 |
| R2 | the decay law as a sampling problem (`docs/RESEARCH_R2_SAMPLING.md`) | **degree-1 upper bound** (`docs/R2_THEOREMS.md`, not yet externally reviewed): Theorems 1–2 [D]; **Theorem C** certified finite table [C] (Hermite 294, Kaiser–Bessel 92 cases incl. the four certified supports); **Proposition H** proves the Hermite rate e^{−c}·c^{d+1} for all x ≥ 2q; the Kaiser–Bessel rate e^{−2c} is **Conjecture R** (missing: uniform norm and envelope estimates M1, M2); no GL(2) bound. Validations that can fail all pass: trial quotient ≤ B_cert 386/386, B_cert ≥ certified value 16/16. Validated by the checks that can fail: every bound lies above the certified lower bound at all four certified supports (16/16), and the Hermite trials' own Rayleigh quotients lie below their bounds (8/8). **Lower half open** and RH-hard; the tested products keep the fitted rate, and the role of zero spacing in the sharp lower bound is open |
| P-ST | square torus x² + y² (Euler, class number 1, conductor 4) vs Z₁ torus x² + 5y² (non-Euler, class number 2, conductor 20), same reflection symmetry | **holds**: Z_sq's tested matrices positive everywhere and on the degree-2 law (α/4π 2.06 / 2.13); Z₁ negative from x ≈ 19.9 (even) and 24.5 (odd). Not an isolation of the Euler product (conductors differ); that is ζ_K vs Z₁ |

**Empirical laws, measured past the transition:**
1. **Rate:** λ_min ≈ e^{−2d·T*}, with T* = 2π(x/Q)^{1/d}, the height where the zero density reaches the window's L/2π. It holds for degree 1, degree 2 (products, elliptic curves, level-1 eigenforms of weight 12–26) and degree 3 (products).
2. **Degree-1 prolate index:** λ ≈ R·ℓ_n(2πx/q), with n = κ + 2·(sector) + 4·(pole) and R = O(1–40), bounded to x/q = 16. In degree 2 the index depends on the gamma factor, but not as n = Σμ + 2s.
3. **Root-number rule (GL(2)):** sign ln(λ_odd/λ_even) = ε.
4. **Euler-blindness:** DH, which has no Euler product, sits inside the odd-character band until it fails abruptly at x ∈ (30, 31.25).

The sections below are in chronological order.


**Observation (ζ).** Let λ(x) be the smallest eigenvalue of ζ's zeros-side Weil form on test functions supported in a window of length L = log x (our normalisation, x = e^{2a}). The best committed upper bound in each sector, at five points (x = 4.95 to 19.9, λ from 1.6e-17 down to 6.3e-96), fits

  ln λ_even(x) = −α x + γ ln x + β,  α = 12.5556 (α/4π = 0.99915), γ = 4.99, β = 15.55, max residual 0.035;
  ln λ_odd(x):  α = 12.549 (α/4π = 0.99859), γ = 7.20, β = 18.85, max residual 0.040.

The ratio λ_odd/λ_even grows like x^{2.29}. These are found after the fact: this observation came before any pre-registration, and it is what the conductor test below was designed to probe.

**Scope.** This is a fit of one model over x ≤ 20; it does not determine the asymptotic rate.
- Zhu's conjecture (arXiv:2608.24827v2, Conj. 12.1) is −ln λ* = C·N(T*)/ln N(T*)·(1 + o(1)), with C ≈ 2π². Its local slope over this range is about 11–12.
- At x = e^{2.99}, our N = 180 upper bound differs from Zhu's literal leading term by 0.59% in −ln λ.
- That is a disagreement with the leading-term prediction at one finite point, not an exclusion of an asymptotic statement. Equally, the fit here does not show that the rate stays 4π as x → ∞.

The rate e^{−4πx} is Slepian's prolate concentration rate e^{−2c} at c = 2πx = 2πλ², with x = λ² in Connes' notation. That is the time–band product of Connes' Sonin picture, where the support is [−λ, λ] in the additive variable and the Fourier transform is also confined to [−λ, λ]. A zero-counting heuristic gives the same rate. F has exponential type L/2, so it can vanish on the zeros only up to the height T* where their density (1/2π)log(T/2π) reaches L/2π, which gives T* = 2πx. The minimum then costs e^{−2T*}.

**Data and reproduction:** `decay_law_mp.py fit-zeta`, which writes `data/connes/decay/fit_zeta.json`. The rule is the best committed upper bound per point:

| a | even | odd |
|---|---|---|
| 0.8 | cosine N = 800, `zhu_window_ext.json` | Legendre 200 modes, `zhu_legendre_upper_odd.json` |
| 1.0, 1.19, 1.3 | min over periodic and antiperiodic N = 140, `antiperiodic/sweep_a*.json` | same |
| 1.495 | N = 180, `zhu_window_a1495_n180.run.txt` | same |

- **These are Rayleigh–Ritz upper bounds, not certified values.** The certified lower bounds at 1.6, 2.38 and 2.6 (`docs/CERTIFICATE_238.md`) lie 1.6–1.8× below them in the even sector.
- **Convergence is not uniform.** At a = 1.495, N = 140 → 180 still lowers λ by 24% (even) and 17% (odd). A uniform factor of 2 error would move ln λ by 0.7, enough to shift γ by about 0.5. It would leave α nearly unchanged, since α is fixed by the 180-unit span of ln λ.

## Pre-registration: the conductor test (written before the scan)

For a real primitive Dirichlet character χ_D of conductor q = |D|:
- the zero density at height T is (1/2π)log(qT/2π), so the crossover moves to T* = 2πx/q;
- in the self-dual variable n/√q, the prolate parameter becomes c = 2πx/q.

Both pictures predict a **scaling collapse**:

**Prediction P-DL1.** For each D ∈ {−3, −4, 5, −7, 8, −20}, the slope of ln λ_χ against x/q is −4π to within 3%, fitted over x/q ∈ [3, 8], in both sectors.

**Kill.** A slope in x/q more than 10% away from −4π for any D. For example, a rate of −4πx/√q would give a slope of −4π√q in x/q.

**Not pre-committed.** Whether the prefactor (γ, β) depends on q, or on the character's parity κ (Γ_R(s) versus Γ_R(s+1)).

**Implementation** (`scripts/decay_law_mp.py`):
- The zeros-side form of L(s, χ_D) is assembled from connes_letter_mp blocks:
  - even χ uses ζ's archimedean block plus log q on the diagonal;
  - odd χ uses the DH block (Γ_R(s+1), built for conductor 5), shifted by log q − log 5;
  - prime terms are Λ(n)χ_D(n)/√n;
  - there is no pole.
- **Validation** (`check`; outputs in `data/connes/decay/checks/`). Zeros are found from sign changes of Λ(½ + it, χ) up to T = 120. They are tested against the form at x = 13, on edge-vanishing test functions only, for D = 5, 8, −3, −4, −7, −20:
  - **odd sector** (single sine modes and b₂ − b₃): relative error 3.9e-6 to 5.0e-4;
  - **even sector** (b₁ + b₂ and b₂ − b₄): relative error 1.4e-9 to 4.3e-7.

  The absolute error is ≲ 1e-4 throughout, consistent with truncating the zero sum at T = 120. The largest relative error, D = −4's b₂ − b₃, is large only because that form value is small (0.155).

## Results

**Data.** `data/connes/decay/scanhi_D*.json`, from the M4. Grid x/q ∈ {2, 3, 4, 5, 6, 7, 8, 10}, both sectors, at N = 5k* and 9k*, where k* = xL/q is the mode index of the crossover frequency T*. λ is found by inverse iteration with FLINT's solver. That finds the eigenvalue of smallest |λ|, which is the minimum when the form is positive semidefinite. The forms here are zeros-side sums Σ_ρ |F(γ)|² over zeros that are all on the line at these heights. The method was spot-checked against mp.eigsy at D = −4, x/q = 3, where every printed digit agrees. Convergence between the two basis sizes is |ln(λ_5k*/λ_9k*)| ≤ 0.17. The first pass at 1.5k* and 2.5k* (`scan_D*.json`) was off by up to e^{1.6}. Fits: `decay_law_mp.py fit` (`fit_hi_3to8.json` for the registered range, `fit_hi.json` for x/q ∈ [3, 10]).

**Verdict on P-DL1 (x/q ∈ [3, 8]).**

| D | q | χ parity | even: linear slope/4π | even: slope/4π with γ ln(x/q) | γ | odd: linear | odd: with log | γ |
|---|---|---|---|---|---|---|---|---|
| 1 (ζ) | 1 | even + pole | 0.920 | 1.002 | 5.32 | 0.885 | 1.017 | 8.57 |
| 5 | 5 | even | 0.994 | 0.994 | −0.02 | 0.962 | 0.987 | 1.66 |
| 8 | 8 | even | 0.989 | 0.999 | 0.61 | 0.958 | 0.994 | 2.32 |
| −3 | 3 | odd | 0.974 | 1.004 | 1.93 | 0.940 | 1.000 | 3.90 |
| −4 | 4 | odd | 0.975 | 1.003 | 1.82 | 0.943 | 0.998 | 3.55 |
| −7 | 7 | odd | 0.968 | 0.987 | 1.26 | 0.940 | 0.984 | 2.83 |
| −20 | 20 | odd | 0.975 | 0.996 | 1.36 | 0.941 | 0.996 | 3.57 |

- **The rate is confirmed.** In the form used for ζ's own fit, ln λ = −α(x/q) + γ ln(x/q) + β, the twelve character cases give α/4π between 0.984 and 1.004.
- **The strict linear version of P-DL1 fails in the odd sectors.** Their linear slopes are 4–6% short of 4π, because of the γ ln term. Five of the six even sectors pass at 3%; χ₋₇ is at 3.2%.
- **The kill did not fire.** No slope is more than 10% from 4π, and a √q law is excluded by a wide margin, since it would give slopes of 4π√q.

**The same rate for all seven functions tested.** Over the sampled windows, λ_min ≈ e^{−4πx/q}, with c = 2πx/q as the prolate parameter. The subleading term differs. The seven functions fall into three approximate groups, which coincide with the archimedean factor and the pole. The table gives ln λ + 4π·x/q:

| class | members | even sector, x/q = 2 / 5 / 10 | odd sector, x/q = 2 / 5 / 10 |
|---|---|---|---|
| even χ, Γ_R(s) | D = 5, 8 | 5.3–5.7 / 5.8–6.6 / 6.3–6.9 | 12.2–12.4 / 14.7–14.9 / 16.8–16.9 |
| odd χ, Γ_R(s+1) | D = −3, −4, −7, −20 | 8.3–9.3 / 9.8–11.2 / 11.0–12.2 | 15.1–15.9 / 18.9–19.3 / 21.5–22.1 |
| ζ, Γ_R(s) + pole | D = 1 | 18.5 / 23.7 / 27.1 | 22.6 / 30.7 / 35.6 |

Within a class the conductors collapse to within about ±0.7 in ln λ, over a range in which ln λ itself falls by about 100. Across classes the offsets are large. At equal x/q, ζ's window minimum exceeds an even character's by a factor of e^{13} (x/q = 2) to e^{21} (x/q = 10), so the pole is worth far more than the conductor. The odd sector sits above the even sector in ln λ: by about 7 at x/q = 2, rising to about 10 at x/q = 10, for the characters; and by 4 rising to 8.5 for ζ.

**Reading** (what these scans support, and what they don't):
- **The rate.** For the six characters tested, it depends on q through x/q over x/q ∈ [2, 10]. Six characters and one control are not enough to establish that the rate is independent of the coefficients in general.
- **The prefactor groups approximately, but is not a function of the class alone.** The power γ and the band ln λ + 4π x/q cluster by gamma factor and pole. But the constant R varies within a class: about 3.5× among the odd characters and about 2× between the two even ones (see "Prolate index").
- **What is not established.** The scans do not show that the law carries no arithmetic beyond q. Nor do they show that DH's minimum equals that of an Euler product with its functional equation. They show only that DH's values fall inside the band spanned by four odd characters at the sampled points (below).

### Prolate index: n = κ + 2·(sector) + 4·(pole) (empirical hypothesis)

The literature review (`docs/LIT_TANGENTS.md`, decay section) found that Connes (arXiv:2602.04022 §6.4, citing Fuchs 1964) gives 1 − χ₂ ∼ (2¹⁴/3)√2 π⁵ e^{−4πx + (9/2) ln x} for the prolate h₄. It also found that our certified ζ minima track the Fuchs–Slepian leakage

  ℓ_n(c) = ½·4√π 8ⁿ c^{n+½} e^{−2c}/n!,  c = 2πx,

with n = 4 in the even sector and n = 6 in the odd sector. The pole conditions f(0) = f̂(0) = 0 remove the lowest prolate from each Fourier class. That suggests a rule for functions without a pole:
- an even character (weight ½) should start at h₀ (even sector) and h₂ (odd sector);
- an odd character (Γ_ℝ(s+1), weight 3/2) should start at h₁ and h₃.

The rule was written into the comparison code before the ratios were computed. It was not registered in a file. With c = 2πx/q, the ratio λ/ℓ_n(c) on the high-N data is:

| D | assigned n (even, odd) | even: λ/ℓ_n at x/q = 2…10 | odd: λ/ℓ_n | flattest n (even, odd) |
|---|---|---|---|---|
| 5 | 0, 2 | 23 → 35 | 3.8 → 6.2 | 0, 2 |
| 8 | 0, 2 | 15 → 20 | 3.0 → 5.4 | 0, 2 |
| −3 | 1, 3 | 8.6 → 14 | 2.7 → 5.6 | 1, 3 |
| −4 | 1, 3 | 6.7 → 9.8 | 2.7 → 4.6 | 1, 3 |
| −7 | 1, 3 | 3.1 → 4.6 | 1.7 → 3.7 | 1, 3 |
| −20 | 1, 3 | 3.3 → 4.1 | 3.9 → 6.4 | 1, 3 |
| DH | 1, 3 | 4.5 → 8.0 (to x/q = 6) | 2.4 → 3.9, then 0.51 at 6.25 | — |
| 1 (ζ) | 4, 6 | 2.1 → 7.9 | 0.37 → 4.6 | 5, 7 (still drifting) |

"Flattest" is the index n ∈ 0..7 that minimises |ln(ratio at x/q = 10 / ratio at x/q = 3)|.

- **The rule holds for every character case.** The assigned index is the flattest in all 12, with ratios O(1–35) and nearly constant over x/q = 3–10.
- **An empirical hypothesis.** Over the sampled range, each tested character's window minimum stays within a nearly constant factor R of a single prolate leakage: λ_min ≈ R·ℓ_n(2πx/q), with n = κ + 2s. Here κ is the character's parity, s = 0 (even sector) or 1 (odd sector), and R runs from 1.55 to 38.5 depending on the character, sector and x/q. Whether R stays bounded as x/q → ∞ has not been tested.
- **ζ fits less cleanly.** It belongs at n = 4/6 according to the pole rule, but its ratio still grows over this range: 2.1 → 7.9 in the even sector here. In `docs/LIT_TANGENTS.md` over x = 5–13.5, it grows 5.6 → 11.0 for the upper bounds and 3.5 → 6.1 for the certified lower bounds. So it sits between h₄ and h₅. A plausible reason, not checked: the pole vectors cosh(u/2) and sinh(u/2) are not exactly the removed prolates, so ζ's constraint removes them only approximately.
- **The constant R is not universal within a class.** Among odd characters, R falls from about 14 (q = 3) to about 4 (q = 20).

### DH: no Euler product, same functional equation as odd χ mod 5

The earlier control scans (`data/controls/upper_dh_*.json`, N = 64) give, at x/q = 2.16 (x = 10.80), ln λ + 4π x/q = 9.01 (even) and 15.96 (odd). The odd-χ class, linearly interpolated between x/q = 2 and 3, gives 8.87 and 15.83 there.

**DH scan** (`data/connes/decay/scanhi_dh_*.json`; FLINT's full eigensolver, so negative eigenvalues are seen; N = 4k* and 7k*). Values are at N = 7k*:

| x/q | x | even λ_min | ln λ + 4π x/q | odd-χ class (even) | odd λ_min | ln λ + 4π x/q | odd-χ class (odd) |
|---|---|---|---|---|---|---|---|
| 2 | 10 | 8.58e-8 | 8.86 | 8.29–9.29 | 6.25e-5 | 15.45 | 15.12–15.94 |
| 3 | 15 | 4.42e-13 | 9.25 | 8.91–10.33 | 1.23e-9 | 17.18 | 17.05–17.40 |
| 4 | 20 | 4.21e-18 | 10.26 | 8.94–10.67 | 1.35e-14 | 18.33 | 17.44–18.64 |
| 5 | 25 | 1.42e-23 | 10.22 | 9.82–11.21 | 1.02e-19 | 19.10 | 18.92–19.34 |
| 6 | 30 | 9.39e-29 | 10.86 | 10.24–11.46 | 6.21e-25 | 19.66 | 19.62–20.03 |
| 6.25 | 31.25 | **−8.33e-30** | — | | 4.56e-27 | 17.89 | |
| 6.5 | 32.5 | **−1.90e-23** | — | | **−2.98e-26** | — | |

**DH stays inside the Euler class, in both sectors, at every point up to x/q = 6.** Then it fails abruptly:
- **Even sector.** The finite-basis minimum changes sign in x ∈ (30, 31.25), which fits the earlier witness at 30.745.
  - A negative λ_N is a witness: a test function with Q < 0, up to floating-point accuracy.
  - A positive λ_N is only an upper bound. So "positive up to x = 30" is a numerical indication, not a proof.
  - The rigorous bracket for the first failure stays [10.80, 30.745); the numerics suggest it lies in (30, 30.745).
- **Odd sector.** It changes sign in (31.25, 32.5). At 31.25 it has already dropped about 1.9 below the class trend, and N = 94 against 158 differ by a factor of 3, a sign that it is about to cross.
- **After the crossing.** The negative eigenvalue grows by six decades between x/q = 6.25 and 6.5, far faster than the e^{−4πx/q} floor shrinks.

**Reading.** At the sampled points, DH's window minimum gives no early warning of its off-line zeros:
- up to x/q = 6, its values fall inside the band spanned by the four odd characters;
- the sign change follows within Δ(x/q) = 0.25 (even), or after one point visibly below the band (odd).

This is consistent with the semilocal null results (`docs/SEMILOCAL_COUPLING.md`). It does not show that DH's minimum *equals* an Euler function's: the band is about ±0.7 wide in ln λ, and DH was sampled at 9 points per sector before its crossover.

## Pre-registration: does R stay bounded? (P-DL2, written before the run)

The prolate-index hypothesis (H3) was found on x/q ∈ [2, 10]. Extend the scan to x/q ∈ {12, 14, 16} for D = 5, 8, −3, −4, −7, −20, in both sectors, at N = 5k* and 9k*. Every saved value must be precision-stable.

**Prediction P-DL2.**
- For each (D, sector), R(x/q) = λ/ℓ_n(2πx/q), with the assigned n = κ + 2s, stays within [0.5, 2]× its value at x/q = 10.
- With the flattest-index metric extended to x/q = 3 against 16, the assigned n is still the flattest in all 12 cases.

**Kill.**
- Any R outside [0.25, 4]× its x/q = 10 value.
- A different n is the flattest in more than 2 of the 12 cases.

Either one means the single-prolate description is a mid-range coincidence.

### P-DL2 result: holds

**Data.** `data/connes/decay/scanext_D*.json`, from the M4, at N = 5k* and 9k*, up to N = 838 (D = −20, x/q = 16).
- Convergence: |ln(λ_5k*/λ_9k*)| ≤ 0.144.
- Precision: all 36 points were recomputed at +40 digits and agree to all 15 saved digits, a relative change ≤ 3.8e-15 (`precision/precext_*.json`).
- Evaluation: `decay_law_mp.py pdl2`, which writes `pdl2.json`. It applies exactly the registered criteria.
- Timing: P-DL2 was committed at 23:53:44 (69f09f7) and the runs started at 23:53:46.

| D | sector | n | R(10) | R/R(10) at x/q = 12 / 14 / 16 | flattest n, x/q = 3 against 16 |
|---|---|---|---|---|---|
| 5 | even | 0 | 35.4 | 0.93 / 0.98 / 0.96 | 0 |
| 5 | odd | 2 | 6.23 | 0.94 / 1.03 / 0.90 | 2 |
| 8 | even | 0 | 20.0 | 0.89 / 0.79 / 1.01 | 0 |
| 8 | odd | 2 | 5.35 | 0.93 / 0.81 / 0.99 | 2 |
| −3 | even | 1 | 14.2 | 1.03 / 1.09 / 1.13 | 1 |
| −3 | odd | 3 | 5.55 | 1.02 / 1.00 / 1.08 | 3 |
| −4 | even | 1 | 9.84 | 1.09 / 1.09 / 1.15 | 1 |
| −4 | odd | 3 | 4.58 | 1.20 / 1.11 / 1.18 | 3 |
| −7 | even | 1 | 4.62 | 1.00 / 0.99 / 0.98 | 1 |
| −7 | odd | 3 | 3.71 | 0.99 / 1.00 / 0.97 | 3 |
| −20 | even | 1 | 4.08 | 1.08 / 1.07 / 1.10 | 1 |
| −20 | odd | 3 | 6.35 | 1.04 / 1.01 / 1.08 | 3 |

- **Both parts of P-DL2 hold.** Every R/R(10) lies in [0.79, 1.20], well inside the predicted [0.5, 2]. The assigned index is the flattest in all 12 cases.
- **This is the first pre-registered, out-of-sample test of the prolate-index hypothesis.** Over x/q ∈ [3, 16], λ falls from 1e-14–1.5e-9 to 3e-85–1e-77, and each character's minimum stays within about ±20% of a fixed multiple of one Fuchs–Slepian leakage.
- **Limits.** These are still finite-basis upper bounds; nothing here is certified. x/q = 16 is not asymptotic. The odd characters' R drifts upward by up to 20% over x/q = 10 → 16, which is within the band but is a trend worth watching.

## Pre-registration: products of L-functions (P-DL3, written before the run)

For a single character, zero counting (T* = 2πx/q) and the prolate picture (c = 2πx/q) predict the same rate. For a product they differ, and the zeros-side forms add, so the test is cheap. Let ρ(x) = λ(product)/λ(largest-conductor factor); by superadditivity, ρ ≥ 1.

| hypothesis | statement | consequence |
|---|---|---|
| **A, worst factor** | The product's window minimum stays at its largest-conductor factor's prolate floor. | ρ stays bounded |
| **B, zero counting of the product** | The product's zero density, degree d and conductor Q = Π q_i, sets T* = 2π(x/Q)^{1/d}, so λ ~ e^{−2T*}. | ρ grows like e^{4π x/q_max − 4π(x/Q)^{1/d}}: about e^{60}–e^{80} at x/q_max = 8 |

**Products** (degree 2): Π1 = ζ·L(χ₋₂₀) = ζ_K (q_max = 20, Q = 20), Π2 = L(χ₋₃)L(χ₋₄) (q_max = 4, Q = 12), and Π3 = L(χ₅)L(χ₈) (q_max = 8, Q = 40). Both sectors, at x/q_max ∈ {2, 4, 6, 8}, with N = 5k* and 9k* at the q_max scale. Values must be precision-stable.

**No directional prediction.** Each hypothesis has its own kill, judged per product:
- **A is killed** for a product if ρ > 30 at any point.
- **B is killed** for a product if ρ < 10⁶ at x/q_max = 8, in either sector.
- If both are killed, the outcome is reported as "neither".

## Pre-registration: a primitive degree-2 L-function (P-DL4, written before any code exists)

Take L(E, s) for a rank-0 elliptic curve over ℚ, in the analytic normalisation (critical line Re s = ½). The gamma factor is Γ_ℂ(s + ½) = Γ_ℝ(s + ½)Γ_ℝ(s + 3/2), the conductor is N, and there is no pole.

**Curves.** E1 = 11a1 (N = 11). E2 is a rank-0 curve with root number +1 and 26 ≤ N ≤ 60, chosen and recorded before any scan (37b1 if it qualifies).

**Grid.** x/N ∈ {2, 4, 6, 8, 12}, both sectors, at N_basis = 5k* and 9k* with k* = xL/N. Values must be precision-stable.

**Three hypotheses, one variable v each.** Fit ln λ = −α·v + γ ln v + β per curve and sector.

| hypothesis | variable v | reasoning |
|---|---|---|
| **A′, character-like** | x/N | conductor scaling as for χ |
| **C′, self-dual** | x/√N | the window rescaled by the self-dual point |
| **B′, zero counting** | √(x/N) | T* = 2π√(x/N), λ ~ e^{−2T*} |

**Rule.** A hypothesis survives for a curve and sector if α/4π ∈ [0.8, 1.2]; otherwise it is killed there. Since √N ≠ 1, at most one of A′ and C′ can survive on the same data. No directional prediction is made.

**Exploratory, not registered:** the flattest prolate index for GL(2).

### P-DL3 result

**Data.** `data/connes/decay/product_*.json`, from the M4. Convergence: |ln(λ_5k*/λ_9k*)| ≤ 0.09. Every value is precision-stable. The pre-registration was committed at 00:38:43; the runs started at 00:40:23. A smoke test at x/q_max = 2, run after the commit, already showed ρ ≈ 1.5e3 for L(χ₋₃)L(χ₋₄).

| product | sector | ρ at x/q_max = 2 | 4 | 6 | 8 |
|---|---|---|---|---|---|
| ζ·L(χ₋₂₀) | even | 2.83 | 4.9e5 | 8.5e11 | 9.1e18 |
| ζ·L(χ₋₂₀) | odd | 1.54 | 1.7e5 | 1.3e11 | 1.0e18 |
| L(χ₋₃)L(χ₋₄) | even | 1.6e3 | 1.5e10 | 9.5e17 | 9.4e25 |
| L(χ₋₃)L(χ₋₄) | odd | 285 | 1.9e9 | 8.9e16 | 6.9e24 |
| L(χ₅)L(χ₈) | even | 1.7e4 | 6.5e11 | 5.8e20 | 4.1e29 |
| L(χ₅)L(χ₈) | odd | 2.2e3 | 6.3e10 | 2.3e19 | 2.0e28 |

- **A is killed for all three products,** in both sectors (ρ > 30).
- **B survives for all three** (ρ ≥ 10¹⁸ at x/q_max = 8).
- **So a product's window minimum is not its worst factor's floor.** The product's denser zero set raises it by many orders of magnitude. Earlier, at x/20 ≤ 1, ζ_K looked as if it sat on its L(χ₋₂₀) floor (the tower test). That was the pre-asymptotic regime.
- **Exploratory, not registered.** A linear fit of ln λ(product) against √(x/Q), on these four points per product, gives α/4π ≈ 1.6–2.0. That is steeper than B's literal e^{−2T*}. Four points cannot fix the form.

## Pre-registration: the degree-2 rate (P-DL5, written before the run)

**Hypothesis.** λ ≈ e^{−2d·T*}, with T* = 2π(x/Q)^{1/d}. For degree 1 this is e^{−4πx/q}, as established above. For a degree-2 product it is e^{−8π√(x/Q)}.

**Run.** The same three products, on a grid v = √(x/Q) ∈ {1.5, 2, 2.5, 3, 3.5}, both sectors. N_basis = 5k and 9k with k = vL, the mode index of T*. Only the product's own λ is computed, and it must be precision-stable.

**Prediction P-DL5.**
- The log-corrected slope α (in ln λ = −α v + γ ln v + β) satisfies α/4π ∈ [1.7, 2.3], i.e. d = 2 within 15%, for every product and sector.
- The linear slope is reported alongside.

**Kill.**
- Any log-corrected α/4π outside [1.4, 2.6].
- Convergence worse than |Δ ln λ| = 0.3 at the largest v, in which case the run is inconclusive rather than killed.

### P-DL5 result: holds

**Data.** `data/connes/decay/productv_*.json`, from the M4.
- Convergence: |ln(λ_5k/λ_9k)| ≤ 0.20, inside the 0.3 limit.
- Every value is precision-stable.
- Timing: the pre-registration was committed at 00:50:20 (496131a) and the runs started at 00:50:38.

| product | sector | log-corrected α/4π | γ | max resid | linear α/4π |
|---|---|---|---|---|---|
| ζ·L(χ₋₂₀) | even | 2.080 | 8.84 | 0.22 | 1.784 |
| ζ·L(χ₋₂₀) | odd | 2.288 | 17.93 | 0.45 | 1.688 |
| L(χ₋₃)L(χ₋₄) | even | 2.182 | 8.78 | 0.37 | 1.889 |
| L(χ₋₃)L(χ₋₄) | odd | 1.788 | −0.77 | 0.23 | 1.813 |
| L(χ₅)L(χ₈) | even | 2.040 | −0.72 | 0.50 | 2.064 |
| L(χ₅)L(χ₈) | odd | 2.052 | 2.39 | 0.74 | 1.972 |

- **All six log-corrected slopes lie in [1.79, 2.29],** inside the predicted [1.7, 2.3]. None is outside the kill range [1.4, 2.6]. One is close to the edge (ζ_K odd, 2.288).
- **So degree-2 products decay like e^{−8π√(x/Q)}.** Together with degree 1 that gives one empirical rate: λ ≈ e^{−2d·T*}, T* = 2π(x/Q)^{1/d}, where T* is the height at which the zero density of the product reaches the window's own density L/2π.
- **This is less sharp than in degree 1.** The residuals (0.22–0.74 in ln λ) are 3–10× larger than for single characters, and the linear slopes spread over 1.69–2.06. The product curves are bumpier, presumably because individual low zeros of the factors matter more. Five points per curve, v ≤ 3.5.
- **The P-DL4 test on a primitive degree-2 L-function** (an elliptic curve) is being built and run separately. It checks whether the same rate holds without any factorisation.

### P-DL4 result: "neither", as registered

The GL(2) agent's work is in `docs/GL2_DECAY.md`, `scripts/gl2_form_mp.py` and `data/connes/gl2/`, merged from its worktree. The main session reviewed it:
- the arch-block validations reproduce;
- the P-DL4 fit reproduces exactly from the scan JSON;
- at x/N = 6 (even sector) both curves' scan eigenvalues match FLINT's full eigensolver to 1e-16.

**Validation.**
- **Archimedean blocks.** The general Γ_ℝ(s+μ) block reproduces `build_form`'s ζ block (μ = 0) and DH block (μ = 1, plus log 5) to about 1e-100. Those two blocks were built in different ways.
- **Zero sums.** The full L(E, s) forms for 11a1 and 37b1 match zero sums over 121 and 144 zeros (T = 120; the argument-principle counts are exact) to 2e-9 to 7e-5. The residual is the predicted truncation tail.
- **Ordering.** The curve choice (11a1, 37b1) was committed at 01:01:11. The scans started at 01:08:06.

**Verdict.** Every hypothesis is killed for both curves in both sectors, so the registered outcome is "neither":
- A′ (x/N): α/4π = 0.20–0.21;
- C′ (x/√N): α/4π = 0.03–0.06;
- B′ (√(x/N), registered at rate e^{−2T*}): α/4π = 1.92–2.01, outside [0.8, 1.2].

**What the data show** (post hoc for P-DL4):
- **Same variable, doubled rate.** B′'s variable is right but its rate is doubled: λ ≈ e^{−8π√(x/N)} = e^{−4T*}. That is the rate λ ≈ e^{−2d·T*} written into the P-DL5 registration (496131a, 00:50:20) for products. It was written before any GL(2) scan existed, and without the agent knowing it, since the agent's worktree predates it. The registered P-DL5 test covered products only, so for elliptic curves this confirms a stated hypothesis rather than a registered test.
- **Collapse in x/N.** At equal x/N the two curves agree to within ln-ratio 0.10–1.68, while λ spans 22 decades.
- **Prolate index (post hoc).** At c = 4π√(x/N) = d·T*, the flattest index is n = 2 (even) and n = 4 (odd) for both curves, with R between 0.56 and 5.2. That is n = μ₁ + μ₂ + 2s with μ = ½ and 3/2, the analogue of the degree-1 rule n = κ + 2s. It is not fully robust: with a different pair of grid points, 11a1 even prefers n = 1.

**Unified empirical description (to be tested).** λ ≈ R·ℓ_n(c), with:
- c = d·T* = 2πd(x/Q)^{1/d};
- n = Σ_j μ_j + 2s + 4·(pole), where Γ = Π_j Γ_ℝ(s + μ_j) and s is the sector;
- R = O(1–40).

## Pre-registration: fresh GL(2) objects (P-DL6, written before any of their data exist)

**Objects.**
- Two rank-0 elliptic curves with root number +1, not 11a1 or 37b1: 15a1 and 19a1 if they qualify, otherwise the next qualifying conductors in [14, 60]. The choice must be committed before any scan.
- Ramanujan's Δ: weight 12, level 1, L(Δ, s) = Σ τ(n) n^{−11/2} n^{−s}. Gamma factor Γ_ℂ(s + 11/2) = Γ_ℝ(s + 11/2)Γ_ℝ(s + 13/2), N = 1, ε = +1, so Σμ = 12.
- Each form must pass the zero-sum validation, at the accuracy achieved for 11a1 and 37b1, before its scan.

**Grid.** v = √(x/N) ∈ {1.5, 2, 2.5, 3, 3.5} for the curves and {2, 2.5, 3, 3.5, 4} for Δ, both sectors. N_basis = 5k and 9k with k = 2vL. Every value must be precision-stable.

**P-DL6a (rate).**
- **Prediction:** ln λ = −αv + γ ln v + β has α/4π ∈ [1.7, 2.3] for every object and sector.
- **Kill:** any α/4π outside [1.4, 2.6].

**P-DL6b (prolate index).** Use c = 4πv, n ∈ {0, …, 20}, and the flattest metric |ln R(v₅)/R(v₂)| between the second and last grid points.
- **Prediction:**
  - curves: n = 2 (even) and n = 4 (odd);
  - Δ: n = 12 (even) and n = 14 (odd).
  - It succeeds if at least 4 of the 6 (object, sector) cases match exactly and all 6 are within ±1.
- **Kill:** for Δ, a flattest n ≤ 8 in either sector, which would mean no μ-dependence; or for a curve, any case off by 2 or more.

### P-DL6: gate decision for Δ (main session, recorded before any Δ scan exists)

The agent turned P-DL6's "zero-sum validation at the accuracy achieved for 11a1 and 37b1" into a numeric gate: relative error ≤ 1e-4, and (Q − Σ)/predicted tail ∈ [0.8, 1.25], at x = 13, T = 120. It committed that gate in 276a496, before any scan. Δ fails the relative-error clause in two of five tests: 6.0e-4 and 2.7e-4.
- **Those two tests have small form values** (Q = 0.023 and 0.32), because Δ's first zero, 9.22, sits above the modes tested. Their absolute errors (1.4e-5 and 8.8e-5) are no larger than those 11a1 and 37b1 passed with.
- **All five residuals equal the predicted truncation tail,** with ratios 0.905–0.926.
- **With Δ's zeros extended to T = 240,** the residuals fall by exactly their tails (6.4× odd, 25× even), and the worst relative error is 9.4e-5. That passes the agent's numeric gate as written.

**Decision.** Δ is accepted as validated, on the T = 240 zero set. This was decided with no Δ scan data in existence; `scan6_delta.json` did not exist at this commit. The curves' gate passes were unaffected.

**State of P-DL6 before Δ.**
- **P-DL6a holds on all four curve cases:** α/4π = 2.13–2.16 (15a1 and 19a1, both sectors).
- **P-DL6b can no longer succeed.** The flattest indices are 1/3 (15a1) and 1/4 (19a1), against the predicted 2/4: one exact, three off by −1, none off by ≥ 2. Δ now decides only between "killed" (flattest n ≤ 8 in either sector) and "not killed, prediction partly missed".

### P-DL6 result

**Data.**
- Curves: `data/connes/gl2/scan6_{15a1,19a1}.json`, from the agent.
- Δ: `scan6_delta.json`, run by the main session after the gate decision above (10d4b24, 02:45:27).
- Fit: `pdl6_fit.json` (`gl2_form_mp.py fit6`).
- Convergence: |ln(λ_5k/λ_9k)| ≤ 0.079 for the curves and ≤ 0.076 for Δ. Every value is precision-stable.

| object | Σμ | sector | α/4π (9k) | γ | linear α/4π | flattest n | predicted n |
|---|---|---|---|---|---|---|---|
| 15a1 | 2 | even | 2.138 | 6.33 | 1.93 | 1 | 2 |
| 15a1 | 2 | odd | 2.161 | 9.39 | 1.85 | 3 | 4 |
| 19a1 | 2 | even | 2.134 | 5.90 | 1.94 | 1 | 2 |
| 19a1 | 2 | odd | 2.157 | 9.51 | 1.84 | 4 | 4 |
| Δ | 12 | even | 2.102 | 17.93 | 1.61 | 14 | 12 |
| Δ | 12 | odd | 2.169 | 23.48 | 1.52 | 16 | 14 |

**P-DL6a holds.** All six log-corrected slopes lie in [2.10, 2.17], inside [1.7, 2.3]; the 5k values agree.
- Together with P-DL4's 11a1 and 37b1 (post hoc, 1.92–2.01), every primitive degree-2 L-function tested decays at λ ≈ e^{−8π√(x/N)} = e^{−4T*}. These are four elliptic curves and one weight-12 cusp form.
- **Caveat:** Δ's linear slopes are only 1.52–1.61, and its log terms are large (γ = 18–23). The registered model credits much of the curvature to the prefactor.

**P-DL6b: not killed, prediction partly missed.**
- Only 1 of 6 cases is exact. The curves sit 1 below n = Σμ + 2s (3 of 4 cases); Δ sits 2 above (both sectors).
- Δ's kill (flattest n ≤ 8) is far from triggered. The index clearly depends on the gamma factor: it moves from 1–4 (Σμ = 2) to 14–16 (Σμ = 12).
- The data reject the exact rule n = Σμ + 2s. Between Σμ = 2 and 12 the index rises by about 12.5, against the rule's 10.
- The metric is coarse: one index step corresponds to about 0.56 in ln R. Two values of Σμ cannot fix the true dependence.

**Standing of the decay-law hypotheses after P-DL1–P-DL6.**
- **Rate:** λ ≈ e^{−2d·T*}, T* = 2π(x/Q)^{1/d}. It held in every registered test:
  - degree 1: P-DL1 (log-corrected), P-DL2;
  - degree-2 products: P-DL5;
  - primitive degree 2: P-DL6a.
  
  P-DL4 failed only because its zero-counting hypothesis was registered at half this rate.
- **Prolate index:** holds exactly for degree 1, where n = κ + 2s matched 12/12 in sample and 12/12 out to x/q = 16. It is only approximate for degree 2: same direction, wrong slope in Σμ.

## Pre-registration: the prolate index against the weight (P-DL7, written before any code or data)

**Objects.** The level-1 normalised Hecke eigenforms of weight k ∈ {16, 18, 20, 22, 26}: f_k = Δ·E_{k−12}, where E_14 = E_4²E_6, E_10 = E_4E_6 and E_8 = E_4². Each space S_k(1) is one-dimensional.
- L(f_k, s) = Σ a_n n^{−(k−1)/2} n^{−s}.
- Gamma factor Γ_ℂ(s + (k−1)/2) = Γ_ℝ(s + (k−1)/2)Γ_ℝ(s + (k+1)/2), N = 1, ε = +1 (k ≡ 0 mod 4) or −1 (k ≡ 2 mod 4).
- **Correction on registration:** the root number for level 1 is i^k, which is −1 for k = 18, 22, 26. Then L(f_k, ½) = 0. Those three forms have a central zero and are kept, with γ = 0 included in the zero sum. Σμ = k.
- Anchors from P-DL6: Δ (k = 12) and the curves (k = 2).

**Validation gate,** applied to each form before its scan, at x = 13, with zeros to T = 240:
- the residual must equal the predicted truncation tail, (Q − Σ)/tail ∈ [0.8, 1.25];
- the absolute error must be ≤ 2e-4.

The gate is stated in absolute terms because the relative criterion misfired on small form values in P-DL6.

**Protocol (P-DL6's).** v = √x ∈ {2, 2.5, 3, 3.5, 4}, both sectors, N_basis = 5k′ and 9k′ with k′ = 2vL. Values must be precision-stable. c = 4πv, n ∈ {0, …, 40}, and the flattest metric |ln R(v₅)/R(v₂)|.

**Predictions.**
- **P-DL7a (rate):** the log-corrected α/4π ∈ [1.7, 2.3] for every form and sector. Kill: outside [1.4, 2.6].
- **P-DL7b (index slope):** fit the flattest even-sector index n_even(k) = a·k + b over k ∈ {12, 16, 18, 20, 22, 26}.
  - The rule n = Σμ + 2s predicts a ∈ [0.8, 1.2]. The two P-DL6 points suggest a ≈ 1.25, which would count as a miss.
  - Kill: a outside [0.6, 1.5].
- **P-DL7c (sector offset):** n_odd − n_even ∈ {1, 2, 3} for every k.
  - The rule predicts 2.
  - Kill: an offset ≤ 0 or ≥ 5 for any k.

## Pre-registration: fresh characters (P-DL8, written before any code or data)

**Objects.** Six real primitive characters not used before: D ∈ {−8, −11, 12, 13, −15, 17}. Each D is a fundamental discriminant; D > 0 is even (κ = 0) and D < 0 is odd (κ = 1).

**Gate.** Before its scan, each form must pass `decay_law_mp.py check` (x = 13, zeros to T = 120, both sectors) with absolute error ≤ 2e-4 in every test.

**Protocol (scanhi's).** x/q ∈ {2, 3, 4, 5, 6, 7, 8, 10}, both sectors, N = 5k* and 9k* with k* = xL/q, inverse iteration. Then a precision recheck at +40 digits of the x/q = 10 points.

**Predictions.**
- **P-DL8a (index):** with c = 2πx/q and n ∈ 0..7, the flattest index under |ln R(10)/R(3)| is n = κ + 2s in at least 11 of the 12 (D, sector) cases. Kill: 9 or fewer.
- **P-DL8b (rate):** the log-corrected slope of ln λ against x/q, over x/q ∈ [3, 8], satisfies α/4π ∈ [0.95, 1.05] in all 12 cases. Kill: any outside [0.9, 1.1].

### P-DL8 result: both parts hold

**Data.** `data/connes/decay/scan8_D*.json`, from the M4. Gate outputs are in `checks/check8_D*.run.txt`; the evaluation (`decay_law_mp.py index-test`) is in `pdl8.json`.
- **Gates:** every form passed, with max absolute error 6.4e-5 to 9.9e-5 against the limit of 2e-4.
- **Ordering:** each pipeline ran check → gate → scan, so a scan could not start unless its gate exited 0. The registration was committed at 02:50:20 (08c8249), the checks started at 02:51:09, and the scans at 02:51:21.
- **Convergence and precision:** |ln(λ_5k*/λ_9k*)| ≤ 0.156. The x/q = 10 points were rechecked at +40 digits, with relative change ≤ 3.3e-15.

| D | κ | even: flattest (pred.) | R(3) → R(10) | α/4π | odd: flattest (pred.) | R(3) → R(10) | α/4π |
|---|---|---|---|---|---|---|---|
| −8 | 1 | 1 (1) | 5.11 → 5.64 | 1.007 | 3 (3) | 3.78 → 5.42 | 1.012 |
| −11 | 1 | 1 (1) | 6.52 → 6.64 | 1.015 | **4 (3)** | 2.14 → 4.20 | 1.026 |
| 12 | 0 | 0 (0) | 13.5 → 14.8 | 1.003 | 2 (2) | 4.33 → 5.88 | 1.002 |
| 13 | 0 | 0 (0) | 17.0 → 14.9 | 0.984 | 2 (2) | 4.10 → 6.01 | 1.012 |
| −15 | 1 | 1 (1) | 1.98 → 2.75 | 1.021 | 3 (3) | 3.12 → 3.63 | 1.007 |
| 17 | 0 | 0 (0) | 4.64 → 7.25 | 1.001 | 2 (2) | 2.24 → 3.24 | 0.997 |

- **P-DL8a holds:** 11 of 12 cases are flattest at the predicted n = κ + 2s. χ₋₁₁ odd is one high, with R growing 2.1 → 4.2.
- **P-DL8b holds:** every log-corrected α/4π lies in [0.984, 1.026], inside [0.95, 1.05].
- **Degree-1 tally.** The index rule now has 12/12 in sample (P-DL1 characters), 12/12 out to x/q = 16 (P-DL2) and 11/12 on fresh conductors (P-DL8). The fresh set includes the composite conductors 12 and 15 and the even conductor 8.

### P-DL7 result: all three predictions killed, as registered

The agent's work is in `docs/GL2_DECAY.md` §16–22 and `data/connes/gl2/scan7_f*.json`. The main session reviewed it:
- refitting from the raw scans reproduces every α/4π and every flattest index;
- the ordering is clean: forms recorded in 0c6e9b6 (02:56:57), after the registration (02:49:34); each gate committed before its scan (f16/f18/f20 at 03:40:05, scans from 03:40:12; f22 at 03:44:16, scan 03:44:44; f26 at 04:05:28, scan 04:05:37).

**Validation.**
- **Coefficients:** f_k = Δ·E_{k−12}, checked by two q-expansion routes, Hecke multiplicativity, Deligne's bound, and the Eisenstein congruence modulo the primes in the numerator of B_k/k.
- **Functional equation:** the AFE confirms ε = i^k. For ε = −1 there is a simple central zero.
- **Zeros:** to T = 240, with exact argument-principle counts.
- **Gates:** |Q − Σ| ≤ 1.6e-5 for every form. A separate test confirms that the ε = −1 forms need exactly one central zero.

**Results.**

| prediction | outcome | detail |
|---|---|---|
| P-DL7a (rate) | **killed** | α/4π for (even, odd) at k = 16: 2.13, 2.55; 18: 2.40, 2.11; 20: 2.50, **2.92**; 22: **3.16**, 2.53; 26: **3.12**, **2.95** (bold: outside the kill limits). γ = 24–89. Linear slopes only 0.67–1.47. |
| P-DL7b (index slope) | **killed** | n_even = 14, 19, 26, 24, 33, ≥ 40 at k = 12, 16, 18, 20, 22, 26; a = 1.88 |
| P-DL7c (sector offset) | **killed** | n_odd − n_even = 3, −4, 5, −5, −5 at k = 16, 18, 20, 22, 26 |

**Why: the grid was pre-asymptotic for high weight.** That was a design error in the registration.
- At v = 2 the minima of f20, f22 and f26 are O(0.1–1.4), and at v = 4 f26's even minimum is still 9.5e-8.
- The leakage asymptotics need c ≫ n, and with n ≈ 20–40 that means v well above 3. Five points with large curvature push the 3-parameter fits to large γ and inflated α.
- So P-DL7 is a valid kill of the predictions *on that grid*. It is not a test of the asymptotic rate.
- Δ (k = 12) on the same kind of grid gave α/4π 2.10–2.17 (P-DL6a), with smaller curvature (γ = 18–23). That result is also less secure than the elliptic-curve ones.

**Post-hoc finding: the root number sets the order of the sectors.**
- At v = 4, ln(λ_odd/λ_even) = +7.35 to +7.52 for every ε = +1 form (Δ, f16, f20), and −6.23 for every ε = −1 form (f18, f22, f26).
- A plausible mechanism: with ε = −1 the central zero γ = 0 constrains F(0). An even test function must make F(0) small, which is costly, while an odd one has F(0) = 0 for free.
- That is also why P-DL7c's offsets flip sign exactly for ε = −1.

## Pre-registration: asymptotic grid and the root-number rule (P-DL9, written before any data)

**Objects.**
- The six level-1 eigenforms k ∈ {12, 16, 18, 20, 22, 26} (Δ and f16–f26, already validated).
- One fresh weight-2 object with ε = −1: **37a1**, the rank-1 curve y² + y = x³ − x (N = 37). Its form must pass the P-DL7 gate, with the central zero included, before its scan.

**Grid.**
- Level-1 forms: v = √x ∈ {6, 7, 8, 9, 10}, so c = 4πv ∈ [75, 126]. That is well past the transition for every weight here.
- 37a1: v = √(x/37) ∈ {2, 2.5, 3, 3.5, 4}.
- Both sectors. N_basis = 5k′ and 9k′ with k′ = 2vL. Precision-stable, with the inv+full eigen-method.

**Predictions.**
- **P-DL9a (rate):** the log-corrected α/4π ∈ [1.7, 2.3] for every object and sector.
  - Kill: any value outside [1.4, 2.6].
  - The linear slope is reported alongside.
- **P-DL9b (root-number rule):** sign(ln(λ_odd/λ_even)) = ε at every grid point, for all seven objects. So 37a1, with ε = −1, must have λ_odd < λ_even.
  - Kill: any mismatch at a point where both λ < 1e-6.

**Not registered.** The prolate index stays exploratory; the metric is coarse and was capped at n = 40.

### P-DL9 result, level-1 part (37a1 pending)

**Data.** `data/connes/gl2/scan9_{delta,f16,f18,f20,f22,f26}.json`, from the M4.
- The scans started at 04:12:42, after the registration (04:12:15).
- v ∈ {6, …, 10}, N_basis up to 836, inv+full eigen-method, every value precision-stable.
- Convergence: |ln(λ_5k/λ_9k)| ≤ 0.065.

| form | ε | even α/4π | γ | resid | odd α/4π | γ | resid | linear (even, odd) | sign(ln λ_odd/λ_even) at v = 6…10 | ln(λ_odd/λ_even) at v = 10 |
|---|---|---|---|---|---|---|---|---|---|---|
| Δ (k = 12) | +1 | 2.012 | 14.2 | 0.09 | 2.021 | 17.3 | 0.05 | 1.87, 1.85 | + + + + + | +9.45 |
| f16 | +1 | 2.058 | 22.8 | 0.17 | 2.074 | 26.8 | 0.17 | 1.83, 1.80 | + + + + + | +9.95 |
| f18 | −1 | 2.049 | 26.7 | 0.08 | 2.003 | 19.5 | 0.07 | 1.78, 1.81 | − − − − − | −8.51 |
| f20 | +1 | 1.980 | 19.6 | 0.11 | 1.952 | 19.3 | 0.08 | 1.78, 1.76 | + + + + + | +9.82 |
| f22 | −1 | 2.020 | 28.3 | 0.03 | 2.020 | 26.1 | 0.09 | 1.73, 1.76 | − − − − − | −8.78 |
| f26 | −1 | 2.046 | 35.4 | 0.13 | 2.059 | 34.1 | 0.06 | 1.69, 1.71 | − − − − − | −8.83 |

- **P-DL9a holds for all twelve level-1 cases:** α/4π = 1.95–2.07, inside [1.7, 2.3] and tightly around 2.
  - On the asymptotic grid the scatter shrinks from P-DL7's 2.1–3.2 to ±0.06, and the residuals from up to 0.6 to ≤ 0.17.
  - So P-DL7's kill came from the grid: past the transition, every weight from 12 to 26 decays at e^{−4T*} = e^{−8πv}.
  - The linear slopes (1.69–1.87) fall slowly with k. That is the c^{n+½} prefactor, with n growing with the weight.
- **P-DL9b holds for all six forms:** 0 mismatches at 30 points. The gap is about +9.5 for ε = +1 and −8.7 for ε = −1, nearly independent of k.
- **Final verdict:** waits for 37a1 (ε = −1, weight 2).

## Pre-registration: degree 3 (P-DL10, written before any data or code change)

**Hypothesis.** The rate law λ ≈ e^{−2d·T*}, with T* = 2π(x/Q)^{1/d}, extends to d = 3. Take v = (x/Q)^{1/3}. Then ln λ = −α v + γ ln v + β should have α = 2d·2π = 12π, i.e. α/4π = 3.

**Objects** (products of three degree-1 factors; the forms add):
- Π_A = L(χ₋₃)L(χ₋₄)L(χ₅), Q = 60;
- Π_B = L(χ₋₄)L(χ₅)L(χ₈), Q = 160;
- Π_C = ζ·L(χ₋₃)L(χ₋₄), Q = 12, with ζ's pole.

**Grid.** v ∈ {1.5, 2, 2.5, 3, 3.5}, both sectors, with x = Q·v³ and N_basis = 5k and 9k, k = vL. Only the product's own λ is computed, and it must be precision-stable.

**Prediction P-DL10.**
- The log-corrected α/4π ∈ [2.55, 3.45] (d = 3, within 15%) for every product and sector.
- Kill: any α/4π outside [2.1, 3.9].
- Inconclusive: convergence worse than |Δ ln λ| = 0.3 at the largest v.

### P-DL10 result: holds (degree 3)

**Data.** `data/connes/decay/productv3_*.json`, from the M4.
- The runs started at 04:45:41, after the registration (59d9deb, 04:45:24).
- Convergence: |ln(λ_5k/λ_9k)| ≤ 0.23. Every value is precision-stable.

| product | Q | sector | log-corrected α/4π | γ | max resid | linear α/4π |
|---|---|---|---|---|---|---|
| L(χ₋₃)L(χ₋₄)L(χ₅) | 60 | even | 2.743 | −5.74 | 0.68 | 2.935 |
| | | odd | 2.782 | −3.01 | 0.86 | 2.883 |
| L(χ₋₄)L(χ₅)L(χ₈) | 160 | even | 2.763 | −8.47 | 0.16 | 3.046 |
| | | odd | 2.794 | −5.46 | 0.82 | 2.977 |
| ζ·L(χ₋₃)L(χ₋₄) | 12 | even | 3.203 | 12.99 | 0.35 | 2.769 |
| | | odd | 3.035 | 10.05 | 0.13 | 2.699 |

- **All six slopes lie in [2.74, 3.20],** inside the predicted [2.55, 3.45].
- **The linear slopes, 2.70–3.05,** are close to 3 as well.
- **Residuals are large, as for degree-2 products** (up to 0.86): the product curves are bumpy.

**The rate law after P-DL1 to P-DL10.** In every registered test measured past the transition, λ ≈ e^{−2d·T*}, T* = 2π(x/Q)^{1/d}. That is, the window minimum falls like e^{−2d} per unit of T*, where T* is the height at which the zero density reaches the window's density L/2π. The tests cover:

| degree | objects | tests |
|---|---|---|
| d = 1 | twelve characters and ζ | P-DL1, P-DL2, P-DL8 |
| d = 2 | three products | P-DL5 |
| d = 2 | four elliptic curves and six level-1 eigenforms of weight 12–26 | P-DL6a, P-DL9a |
| d = 3 | three triple products | P-DL10 |

These are finite-basis upper bounds over finite ranges, not proofs.

### P-DL9 final result: both parts hold (all seven objects)

**37a1** (rank 1, ε = −1, N = 37) was validated by the agent (`docs/GL2_DECAY.md` §23–25). The main session reviewed it: `fit9` over all seven objects reproduces, and the v = 3 eigenvalues match FLINT's full eigensolver to ≤ 7e-16.

**Validation.**
- **Coefficients:** a_p by point counting; Euler terms match the Dirichlet log-derivative to 2e-31.
- **Root number:** ε = −1 from the AFE. The opposite sign fails by ≥ 0.018, about 10³⁹ times the working accuracy.
- **Rank 1:** L(½) = 0, and L′(½) = 0.3059997738, the known value.
- **Zeros:** 340 to T = 240, with exact argument-principle counts N(T) + ½ (a simple central zero).
- **Gate:** all 6 tests pass, including the central-zero test, with |Q − Σ| ≤ 2e-5 and tail ratios 0.97–0.98.
- **Ordering:** verification committed at 04:17:49, gate at 05:12:10, scan started at 05:12:15.

| v | λ_even | λ_odd | ln(λ_odd/λ_even) |
|---|---|---|---|
| 2 | 8.15e-14 | 1.81e-16 | −6.11 |
| 2.5 | 5.96e-19 | 7.35e-22 | −6.70 |
| 3 | 4.11e-24 | 4.17e-27 | −6.89 |
| 3.5 | 3.84e-29 | 1.70e-32 | −7.72 |
| 4 | 2.27e-34 | 1.17e-37 | −7.57 |

**P-DL9a holds for all seven objects:** α/4π = 1.85–2.07. For 37a1 the values are 1.889 (even) and 1.850 (odd), with small log terms (γ = 0.3 and −3.4) and linear slopes of 1.88 and 1.94.

**P-DL9b holds:** sign(ln λ_odd/λ_even) = ε at all 35 grid points, and every point had both λ < 1e-6.
- The weight-2 pair 37a1 (ε = −1) and 37b1 (ε = +1, P-DL4) share the same conductor and gamma factor. Their ln(λ_odd/λ_even) are about −7 and +8. The sector order follows ε alone.
- **Caveat:** the ε = −1 evidence at weight 2 is one curve.

## Interpretation: the zero deficit (heuristic, not derived)

Let D(t) = (1/2π) log(Q (t/2π)^d) be the smooth zero density of a degree-d L-function of conductor Q. A test function on a window of length L = log x has exponential type L/2, so its own zero density is L/2π. T* is where the two meet. The deficit of zeros below T* is

  N_def = ∫₀^{T*} (L/2π − D(t)) dt = (d/2π) ∫₀^{T*} log(T*/t) dt = d·T*/2π.

So the rate law reads **λ_min ≈ e^{−4π·N_def}**: a factor of e^{−4π} for each zero the window could accommodate but the L-function does not supply. The prolate parameter is c = 2π·N_def in every case: 2πx/q for characters, and 4π√(x/N) = d·T* for degree 2.

**Not derived.** It comes from matching the fitted rates. N_def uses the smooth density; individual zeros near T* do not produce jumps of 4π in ln λ, which are not seen.

## Exploratory: the GL(2) prolate index on the asymptotic grid (post hoc, P-DL9 data)

**Metric.** A continuous index: the n at which the least-squares slope of ln(λ/ℓ_n(c)) against ln c, over the five v points, crosses zero. c = 4πv.

| form | k | ε | n_even | n_odd | n_even − k | n_odd − k |
|---|---|---|---|---|---|---|
| Δ | 12 | +1 | 12.49 | 14.77 | 0.49 | 2.77 |
| f16 | 16 | +1 | 16.59 | 19.03 | 0.59 | 3.03 |
| f20 | 20 | +1 | 21.06 | 23.57 | 1.06 | 3.57 |
| f18 | 18 | −1 | 21.33 | 18.66 | 3.33 | 0.66 |
| f22 | 22 | −1 | 25.77 | 23.57 | 3.77 | 1.57 |
| f26 | 26 | −1 | 30.38 | 27.74 | 4.38 | 1.74 |

- **Within each ε class the index grows about one-for-one with k:** slopes of about 1.07 (ε = +1) and 1.14 (ε = −1).
- **ε = +1:** n ≈ k + 0.7 (even) and k + 3.1 (odd), close to the rule n = Σμ + 2s = k + 2s, offset by about +0.7 to +1.1.
- **ε = −1:** the sectors swap, n ≈ k + 3.8 (even) and k + 1.3 (odd). The central zero raises the even index by about 3 and lowers the odd one by about 2.
- **The earlier pooled slopes were artefacts:** P-DL7b's 1.88 (short grid, two-point metric) and 1.29 here mix the two ε classes.

## Pre-registration: fresh GL(2) objects for the index rule (P-DL11, written before any data)

**Objects** (newforms with a single eta-product, dimension-1 spaces):
- g4 = η(τ)⁴η(5τ)⁴ ∈ S_4(Γ₀(5));
- g6 = η(τ)⁶η(3τ)⁶ ∈ S_6(Γ₀(3));
- g8 = η(τ)⁸η(2τ)⁸ ∈ S_8(Γ₀(2)).

For each, the root number ε is determined from the AFE (both signs tested) and committed before its scan, together with the zero-sum gate (P-DL7 gate, with the central zero included if ε = −1). Also on the new grid: 11a1 (ε = +1) and 37a1 (ε = −1). Their forms are validated already; their data on this grid do not exist.

**Grid.** v = √(x/N) ∈ {6, 7, 8, 9, 10}, both sectors, N_basis = 5k′ and 9k′ with k′ = 2vL, precision-stable, inv+full.

**Predictions.**
- **P-DL11a (rate):** α/4π ∈ [1.7, 2.3]. Kill: outside [1.4, 2.6].
- **P-DL11b (index):** using the continuous metric above, n − k lies within ±1.2 of δ(ε, s), where:
  - δ(+1, even) = 0.7 and δ(+1, odd) = 3.1;
  - δ(−1, even) = 3.8 and δ(−1, odd) = 1.3.

  It holds if at least 8 of the 10 (object, sector) cases are within ±1.2. Kill: more than 2 cases off by more than 2.
- **P-DL11c (ε rule):** sign ln(λ_odd/λ_even) = ε at every point. Kill: any mismatch where both λ < 1e-6.

### P-DL11 result

**Data.** `data/connes/gl2/scan11_{g4,g6,g8,11a1,37a1}.json`, from the M4. The evaluation (`gl2_form_mp.py fit11`) is in `pdl11_fit.json`.
- **Ordering:** the g4, g6 and g8 scans started at 06:22:53, after the validation and gates were committed (b916318, 06:22:02). The 11a1 and 37a1 scans started at 05:24:51, after the registration (3778445, 05:24:28).
- **Convergence:** |ln(λ_5k/λ_9k)| ≤ 0.069. Every value is precision-stable.
- **Validation (agent):** g4, g6 and g8 are Hecke eigenforms (checked to n = 2000), satisfy Deligne's bound, and have ε = +1 by the AFE with both signs tested. Their zeros to T = 240 have exact counts, and all three pass their gates.

| object | k | N | ε | α/4π (even, odd) | n* − k even (pred. δ) | n* − k odd (pred. δ) |
|---|---|---|---|---|---|---|
| g4 | 4 | 5 | +1 | 1.95, 2.00 | −0.59 (0.7) | +1.88 (3.1) |
| g6 | 6 | 3 | +1 | 2.00, 2.03 | −0.66 (0.7) | +1.87 (3.1) |
| g8 | 8 | 2 | +1 | 2.00, 1.96 | −0.59 (0.7) | +1.87 (3.1) |
| 11a1 | 2 | 11 | +1 | 2.06, 2.05 | −0.71 (0.7) | +1.62 (3.1) |
| 37a1 | 2 | 37 | −1 | 2.03, 2.02 | +1.32 (3.8) | −0.67 (1.3) |


- **P-DL11a holds:** α/4π = 1.95–2.06.
- **P-DL11c holds:** 25 of 25 points.
- **P-DL11b: not killed, prediction missed.** None of the 10 cases is within ±1.2; one is off by more than 2 (37a1 even, 2.48), and the kill needed more than two.

**Post hoc.** The miss is a nearly constant shift.
- **ε = +1, level N > 1:** the four objects g4, g6, g8 and 11a1 span weights 2–8 and levels 2–11. They give n − k = −0.64 ± 0.06 (even) and +1.81 ± 0.13 (odd).
- **Level 1:** the forms sit about +1.3 higher in both sectors (+0.5 to +1.1 even, +2.8 to +3.6 odd).
- **ε = −1:** 37a1 gives +1.32 (even) and −0.67 (odd). Those are the ε = +1 level-N offsets with the sectors swapped, the odd value lowered by about 0.5.
- **Within a class the index tracks k one-for-one; between level 1 and level N > 1 there is an offset of about 1.3.** Its cause is not known.
  - One candidate is that the window variable L = log(Nv²) enters the prefactor besides c = 4πv. For level 1, ln L moves with ln c; for N > 1 it moves less.
  - This has not been tested.

## Pre-registration: the level-N GL(2) index offsets on fresh curves (P-DL12, written before any data)

**Hypothesis** (from the P-DL11 post-hoc pattern). For weight-k newforms of level N > 1, with the continuous index n* (the P-DL11 metric, c = 4πv, v = √(x/N) ∈ {6, …, 10}):
- ε = +1: n* − k = −0.64 (even) and +1.81 (odd);
- ε = −1: n* − k = +1.32 (even) and −0.67 (odd), from 37a1.

**Objects.** Four fresh elliptic curves, never scanned:
- rank 0, ε = +1: **14a1** and **17a1**;
- rank 1, ε = −1: **43a1** and **53a1**.

If a label does not qualify, the next qualifying conductor in the same class is used, recorded before any scan. Each curve must be verified (conductor, a_p, ε by the AFE with both signs, rank by L(½) and L′(½)) and pass the P-DL7 gate, with the central-zero test if ε = −1. That is committed before its scan.

**Protocol.** P-DL11's: v ∈ {6, …, 10}, both sectors, N_basis = 5k′ and 9k′, precision-stable, inv+full.

**Predictions.**
- **P-DL12a (index):** n* − k lies within ±0.5 of the hypothesis value in at least 6 of the 8 (curve, sector) cases. Kill: more than 2 cases off by more than 1.0.
- **P-DL12b (rate):** α/4π ∈ [1.7, 2.3]. Kill: outside [1.4, 2.6].
- **P-DL12c (ε rule):** sign ln(λ_odd/λ_even) = ε at every point. Kill: any mismatch with both λ < 1e-6.

### P-DL12 result: holds (all three parts)

**Data.** `data/connes/gl2/scan12_{14a1,17a1,43a1,53a1}.json` from ml (Docker, pinned image); gates committed at 20:05:57 (de84949), scans started 20:07:35 local. 40 points, N_basis up to 1551, |ln(λ_5k/λ_9k)| ≤ 0.067, every value precision-stable. Evaluation `gl2_form_mp.py fit12` → `pdl12_fit.json`.

| curve | ε | sector | n* − k | predicted | |off| |
|---|---|---|---|---|---|
| 14a1 | +1 | even | -0.54 | -0.64 | 0.10 |
| 14a1 | +1 | odd | +1.50 | +1.81 | 0.31 |
| 17a1 | +1 | even | -0.78 | -0.64 | 0.14 |
| 17a1 | +1 | odd | +1.61 | +1.81 | 0.20 |
| 43a1 | -1 | even | +1.69 | +1.32 | 0.37 |
| 43a1 | -1 | odd | -0.52 | -0.67 | 0.15 |
| 53a1 | -1 | even | +1.43 | +1.32 | 0.11 |
| 53a1 | -1 | odd | -0.57 | -0.67 | 0.10 |

- **P-DL12a holds:** 8 of 8 within ±0.5 of the offsets found post hoc in P-DL11 (−0.64/+1.81 for ε = +1; +1.32/−0.67 for ε = −1), now on curves never scanned.
- **P-DL12b holds:** α/4π = 1.88–2.06 (linear 1.96–1.98).
- **P-DL12c holds:** sign ln(λ_odd/λ_even) = ε at 20 of 20 points; the ε = −1 evidence now rests on three curves (37a1, 43a1, 53a1).

## Direction (2026-10-03)

With R1 answered (exclusion is a witness computation, and under the stated hypothesis the forecast reach of larger certificates is far below the verified height), the program proceeds with **cheap instruments**: the decay law and its sampling-theory reading (R2), witnesses and controls, and small pre-registered tests. No further certificate runs are planned; support 3.2 stays unrun. New agents start from `docs/PRIMER.md`.

## Pre-registration: the square torus (P-ST, written before any data)

**Question (user's).** Are the zeros on the critical line the zeros a torus produces? The spectral zeta of a flat torus is the Epstein zeta of its lattice form. Two tori with the same seam (functional equation) but different arithmetic:
- **Z_sq**, the torus of x² + y² (class number 1): Z_sq(s) = 4ζ(s)L(s, χ₋₄), an Euler product; zeros-side window form = ζ's form + L(χ₋₄)'s form (forms add over factors; the constant 4 is irrelevant). Degree 2, conductor 4.
- **Z₁**, the torus of x² + 5y² (class number 2): ½[ζ_K + L(χ₋₄)L(χ₅)], no Euler product, off-line zeros (Davenport–Heilbronn 1936; our census in `docs/EPSTEIN.md`). Degree 2, conductor 20. Measured first positivity failure: x ∈ [15.800, 19.844), i.e. 2a ∈ [2.760, 2.988) [variable corrected after review; the registration originally wrote 2a for x].

**Run.** λ_min of the zeros-side form, both sectors, at x = e^{2a} for 2a ∈ {1.6, 2.0, 2.38, 2.6, 2.99, 3.2, 3.5}, at N = 64 and N = 96 (full eigensolver with fallback; precision-stable). `scripts/square_torus_mp.py` → `data/connes/decay/square_torus.json`.

**Predictions.**
- **P-ST1 (sign).** Z_sq's λ_min > 0 at all seven supports in both sectors. Z₁'s even-sector λ_min < 0 at x ∈ {19.89, 24.5, 33.1} and > 0 at x ≤ 13.46. Kill: any negative Z_sq value (a witness against RH for ζ·L(χ₋₄), or a bug), or Z₁ even positive at x ≥ 19.89 with N = 96.
- **P-ST2 (rate).** Z_sq follows the degree-2 law at its own conductor: ln λ = −αv + γ ln v + β with v = √(x/4) has α/4π ∈ [1.7, 2.3] in both sectors (fitted over the seven points). Kill: outside [1.4, 2.6].
- **P-ST3 (sectors).** Z_sq has λ_odd > λ_even at every support (the pole is present). Kill: any reversal.

**What the outcome means.** If P-ST1–3 hold, "zeros on the seam" is decided by the Euler product (class number), not by torus geometry or the symmetry: two tori with the same seam differ exactly as ζ_K and Z₁ do.

### P-ST result: holds (all three parts)

**Data.** `data/connes/decay/square_torus.json` (M3; registered 07:25:33 in 77ce486, run started 07:25:45). N = 64 and 96, precision-stable; Z_sq convergence |ln(λ₆₄/λ₉₆)| ≤ 0.074.

| 2a | x | Z_sq even | Z_sq odd | Z₁ even | Z₁ odd |
|---|---|---|---|---|---|
| 1.6 | 4.95 | 6.89e-4 | 1.02e-1 | 1.85 | 1.49 |
| 2.0 | 7.39 | 5.58e-6 | 1.05e-3 | 1.05 | 1.43 |
| 2.38 | 10.80 | 1.33e-8 | 9.01e-6 | 0.353 | 1.20 |
| 2.6 | 13.46 | 1.34e-10 | 6.83e-8 | 0.121 | 0.780 |
| 2.99 | 19.89 | 3.48e-14 | 3.31e-11 | **−5.75e-5** | 0.180 |
| 3.2 | 24.53 | 1.29e-16 | 1.72e-13 | **−0.124** | **−0.0220** |
| 3.5 | 33.12 | 8.29e-21 | 1.38e-17 | **−0.500** | **−0.475** |

- **P-ST1 holds.** Z_sq is positive at all 14 points. Z₁'s even sector changes sign between x = 13.46 and 19.89, inside its measured bracket [15.800, 19.844); its odd sector between 19.89 and 24.53, where R1's σ₁ analysis put it ((20, 22)).
- **P-ST2 holds.** With v = √(x/4): α/4π = 2.061 (even, γ = 7.2) and 2.133 (odd, γ = 11.6); linear slopes 1.76 and 1.64 (the small-x points are pre-asymptotic at this conductor).
- **P-ST3 holds.** λ_odd > λ_even at every support.

**Reading (corrected after review).** Both tori share the reflection symmetry s ↔ 1 − s; they differ in lattice shape, in conductor (4 against 20, so their completed functional equations differ), and in class number, i.e. in whether the spectral zeta is an Euler product. The positive values for Z_sq are Rayleigh–Ritz upper bounds in the tested matrices: they are consistent with RH for ζ·L(χ₋₄), which remains conjectural, and they do not establish that its zeros lie on the line. Z₁'s negative values are witnesses: it has zeros off the line (Davenport–Heilbronn), as its census shows. So the result is: two tori with the same symmetry behave differently, and the symmetry alone does not keep zeros on the seam. It does **not** isolate the Euler product as the sole cause, because the conductors differ. The equal-conductor isolation in this work is ζ_K against Z₁ (`docs/EPSTEIN.md`): same degree, gamma factor and conductor 20, and ζ_K = Z₁ + Z₂ is the spectral zeta of the union of the two discriminant-−20 tori; it is certified positive at x = e^{2.99} where Z₁ is certified negative. A single-torus version of that comparison would need Z₂ (2x² + 2xy + 3y²) in the window instrument, which `docs/EPSTEIN.md` notes is blocked by a₁ = 0.

**External review of 0a83c93 (three findings, all accepted).** (1) Positive finite matrices are upper bounds and cannot establish critical-line zeros: the reading now says "consistent with RH", not "keeps its zeros on the line". (2) The comparison does not isolate the Euler product, since the conductors are 4 and 20: the reading now states this and points to ζ_K vs Z₁ for the isolation. (3) The failure bracket was labelled 2a but is in x: corrected in the registration with a marker. The reviewer independently recomputed twelve eigenvalues at 120 digits (agreement 5e-15) and checked the square-lattice coefficient identity exactly.

### R2: external review of the brief, and the corrected validation (2026-10-03)

An auditor reviewed the R2 *brief* and raised five points; how each fared against what R2 produced:

1. **"measured/bound ≤ 1" is not a validation.** Correct. The Ritz minimum U_N and the bound B are both upper bounds on λ*, so neither orders the other and the check cannot fail. The R2 doc and the merge note quoted "334/334 ratios ≤ 1" as validation; that was wrong. The checks that can fail, now done: (a) the trial function's own Rayleigh quotient ≤ B, 8/8 Hermite points (ratio 0.036–0.081; the Kaiser–Bessel files do not record their trial quotients, so (a) is open for them); (b) B ≥ the certified lower bound at the same x, **16/16** at x = 4.953, 10.805, 13.464, 19.886, both sectors, both constructions (`data/research_r2/bounds_at_certified_supports.json`; Kaiser–Bessel B/cert 1.1e4–8.5e6, Hermite 1.3e9–3.5e50).
2. **The ordinary prolate on the log window has parameter c_PW = L·T*/2, not d·T*.** Correct, and borne out: R2's conditional lower-bound sketch is exactly 1 − λ₀(L·T*/2), rate e^{−2πx log x} for ζ, while c = d·T* comes only from the E-map construction in the additive variable. Connecting the two is the open step.
3. **Landau density is not a quantitative sampling inequality; the zeros are not uniformly separated.** Correct. R2 did not infer from "super-Nyquist": its conditional bound needs an explicit separated-subsequence hypothesis that is not known even under RH.
4. **Central zeros: the zero sum needs m₀|F(0)|², and the root number fixes only the parity of m₀.** Correct and not yet handled. Every object tested has m₀ ∈ {0, 1} (ranks 0 and 1, verified by exact counts), so the root-number rule as tested is really a statement about m₀ ∈ {0, 1}. An ε = +1 curve of rank 2 (m₀ = 2, e.g. 389a1) is the sharp test of whether the rule follows ε or the central multiplicity.
5. **Separate the leading-rate conjecture from the strong leakage conjecture; non-Euler controls go negative.** R2 did this (D1 rate, D2 index, D3 bounded ratio; DH satisfies the upper bound and still fails).

The upper-bound theorems stand; the validation claim is corrected. The brief-writing lessons (gates must be able to fail; no `reset --hard` in agent briefs) are recorded for future briefs.

## Pre-registration: rank 2 — does the sector order follow ε or the central multiplicity? (P-R389, written before any code or data)

**Why.** Every object tested so far has m₀ ∈ {0, 1}, and ε = −1 always came with m₀ = 1 (rank 1). So the root-number rule (80/80) cannot tell whether the sector order follows ε or the central zero. The R2 review (point 4) noted this. In the Weil form the central zero enters as m₀·|F(0)|², and F(0) = 0 automatically in the odd sector.

**Object.** **389a1**, y² + y = x³ + x² − 2x, conductor 389, the smallest-conductor rank-2 curve: ε = +1 and m₀ = 2. Before any scan it must be verified and the verification committed:
- the conductor from the discriminant and c4;
- a_p by point counting;
- ε by the AFE with both signs;
- L(½) = 0 and L′(½) = 0 (vanishing forced by ε = +1 and L(½) = 0) and L″(½) ≠ 0;
- zeros to T = 120 with the argument-principle count N(T) + m₀/2 = N(T) + 1 at every height;
- the zero-sum gate, |Q − Σ| ≤ 2e-4 and tail ratio in [0.8, 1.25] in every test, with the central term m₀|F(0)|² = 2|F(0)|²;
- and the central-zero test, which must need exactly m₀ = 2 (with m₀ = 1 or 0 it must fail by |F(0)|² or 2|F(0)|²).

**Grid.** P-DL9's 37a1 protocol: v = √(x/389) ∈ {2, 2.5, 3, 3.5, 4}, both sectors, N_basis = 5k′ and 9k′, k′ = 2vL, inv+full, precision-stable.

**Two hypotheses, decided by sign(ln λ_odd/λ_even) at the five grid points.**

| hypothesis | statement | prediction for 389a1 |
|---|---|---|
| **H_ε** (R2's reading: the ground state sits in the sector (−1)^s = ε) | the order follows the root number | sign = +1 (odd above even) |
| **H_m** (the central zero penalises even test functions through F(0)) | the order follows whether a central zero is present | sign = −1 (odd below even), as for the rank-1 curves |

**Decision.**
- H_ε holds if sign = +1 at ≥ 4 of 5 points; H_m holds if sign = −1 at ≥ 4 of 5 points; otherwise inconclusive.
- No directional prediction is made by the main session.
- Recorded alongside, not decisive: |ln(λ_odd/λ_even)| against 37a1's 6.1–7.7 and 37b1's 7.1–8.8.

**Also registered (rate).** The log-corrected α/4π ∈ [1.7, 2.3] in both sectors (kill outside [1.4, 2.6]), as for every other GL(2) object.

### P-R389 result: the sector order follows the central zero, not the root number

**Ordering.** Registered at 11:16:20 (417ed29). Verification committed at 11:21:07 (cafd9b2): conductor 389; ε = +1, with the wrong sign failing by 0.2–4.7; analytic rank 2, with L(½) ≈ 1e-41, L′(½) ≈ 4e-25 and L″(½) = 1.51863300058, twice the BSD leading coefficient 0.7593. Gate committed at 11:58:12 (d15b7c9): 188 zeros to T = 120, argument-principle count N(T) + 1 at all six heights (m₀ = 2.0), Z(2e-4)/Z(1e-4) = 4.00000005; zero-sum tests 6/6, tail ratio 0.953–0.955; the central test needs exactly m₀ = 2 (off by F(0)² = 2.565 with m₀ = 1, by 5.13 with m₀ = 0). Scan started at 11:58:14 on the M4. Data: `data/connes/gl2/scanR389.json`.

| v | x | λ_even | λ_odd | ln(λ_odd/λ_even) |
|---|---|---|---|---|
| 2 | 1556 | 6.76e-17 | 4.21e-20 | −7.38 |
| 2.5 | 2431 | 4.21e-22 | 1.29e-25 | −8.09 |
| 3 | 3501 | 1.42e-27 | 4.67e-31 | −8.02 |
| 3.5 | 4765 | 6.46e-33 | 1.19e-36 | −8.60 |
| 4 | 6224 | 3.71e-38 | 5.85e-42 | −8.76 |

Convergence |ln(λ_5k/λ_9k)| ≤ 0.097; every value precision-stable.

- **Decision: H_m holds** (sign −1 at 5/5 points). **H_ε is refuted.** An ε = +1 curve with a double central zero orders its sectors like the rank-1 (ε = −1) curves, with a somewhat larger gap (−7.4 to −8.8 against 37a1's −6.1 to −7.7).
- **Corrected 2026-10-05: the central weight does not explain the larger gap.** This entry first said the larger gap was "consistent with the weight m₀|F(0)|² = 2|F(0)|²". That is withdrawn.
  - **Why.** The central term is a rank-one update of the even matrix, and it saturates. Once λ is small, any weight acts as the hard constraint F(0) = 0, up to a relative deficit of order λ/L (from the secular equation).
  - **Evidence on ζ.** At x = 5, 8 and 11, weight 1 and the hard constraint differ by 2e−12, 8e−27 and 6e−42; weight 2 halves that.
  - **Consequence for 389a1.** At its points, λ_even runs from 7e−17 down to 4e−38 with L = log x ≈ 7.3–8.7, so the relative deficit is at most about 1e−17. One central zero and two give the same even minimum. The larger gap must come from the other terms: the conductor (389 against 37), the other zeros, the primes, and the different x at equal v.
  - **Not checked** on 389a1 itself.
- **Rate holds:** α/4π = 1.953 (even) and 1.937 (odd).
- **Consequence (worded after review).** The "root-number rule" (80/80 on ranks 0 and 1) is better described as *consistent with central multiplicity controlling sector order*: the central term m₀|F(0)|² charges only even test functions. It saturates: to relative order λ/L, any m₀ ≥ 1 acts as the constraint F(0) = 0. So the conjecture can only concern m₀ = 0 against m₀ > 0, not the size of m₀. This is a conjecture supported by the tested bases, not a universal rule: with the other terms fixed the central term raises the even minimum, but that alone does not guarantee it exceeds the odd one. For ε = −1, m₀ is odd, so the two readings always coincided before this test. R2's GL(2) statement that "the ground state sits in the sector (−1)^s = ε" (`docs/RESEARCH_R2_SAMPLING.md` §6.2) holds only when m₀ ≤ 1; the index assignment there (m = 0 if (−1)^s = ε) should read m₀ > 0 for "ε = −1". This came from the auditor's point 4 on the R2 brief.
