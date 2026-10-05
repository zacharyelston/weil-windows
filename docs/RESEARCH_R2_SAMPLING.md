# R2: The Decay Law as a Sampling Problem

Date: 2026-10-03. Branch `feat/weil-gram-instrument` (worktree), base `fbac9a8`. Brief: `docs/BRIEF_RESEARCH_HANDOFF.md` §5, question R2. New code in `scripts/research_r2/`, data in `data/research_r2/`; no existing script or document was modified. Every table below is printed by `scripts/research_r2/summarize.py` from the committed JSON.

**Labels.** [R] read in the source; [A] abstract or search summary only; [S] second-hand (another paper's or an rh2 document's account); [D] derived here (a proof, or a sketch where marked); [N] computed here (multiprecision floating point: finite-basis Rayleigh quotients are upper bounds, nothing is certified); [E] rh2's registered empirical laws; [U] unverified.

## 0. Six lines

1. **One construction explains the upper half of the law.** Truncate the Poisson (theta) image E(ψ)(w) = w^{1/2} Σ χ(n) ψ(nw/√q) of a Fourier eigenfunction ψ to the window. Then the Weil form of the truncation is **exactly** σ Σ_ρ F_out(γ_ρ)², the out-of-window tail's transform summed over *all* nontrivial zeros, on or off the line (Theorem 1 [D]). This is Connes' "near-radical" k_λ (arXiv:2602.04022 §6.4 [R]) made into an identity. It uses no Euler product.
2. **Unconditional theorems.** Bounding that sum with only 0 < Re ρ < 1 and an explicit zero count gives λ_min ≤ (2A²/‖f‖²) ∫_{A/B}^∞ N(t) t⁻³ dt (Theorem 2 [D]). With Hermite ψ (for ζ's even sector this is Riemann's own k, whose transform is Ξ) the bound is fully explicit and decays like e^{−c} (Corollary 3). With Kaiser–Bessel ψ it decays at the **conjectured rate e^{−2c}**, c = 2πx/q, with a polynomial loss c^{3.5–12} (Theorem 4). Measured/bound ≤ 1 at all 244 (Hermite) and 82 (Kaiser–Bessel) degree-1 grid points. The KB ratio lies in [6.3e-8, 4.8e-3].
3. **The sharp construction is near-optimal.** With time-limited prolates of index n = κ + 2s + 4·pole, the E-map trial function's Rayleigh quotient (rh2's explicit-formula matrices, so an unconditional upper bound) lies within a factor 1.05–2.48 of the measured minimum at all 106 tested points: median 1.28, 92 of 106 below 1.5 [N]. The degree-1 index rule is therefore the index of the optimal E-map trial: parity from the Γ-factor, Fourier class from the sector, one class member removed by the pole (§5).
4. **GL(2) is the same construction with Voronoi in place of Poisson.** The self-reciprocal transform is the Hankel transform of order k − 1, so the leakage is that of Slepian's generalised prolates at c = 4π√(x/N) = d·T*. The exact generalised-prolate leakage reproduces the measured continuous index of all 15 GL(2) objects within ±0.6 (|shift| ≤ 0.59), and the root-number rule: the ground state sits in the sector (−1)^s = ε (15/15 signs, gap within 2 in ln λ) [N]. Two thirds of the "level-1 offset" (~1.3) is a finite-c weight effect; ~0.4 remains [N].
5. **The lower bound is the whole open content.** rh2's numbers are upper bounds, so they can refute but never confirm it. It implies RH (Yoshida), and it cannot follow from the functional equation plus the smooth zero density: Davenport–Heilbronn satisfies the same upper-bound theorem and goes negative. Under RH it is a lower sampling inequality for the zeros in PW_{L/2}. RH plus a local sampling hypothesis above T* would give only the weaker rate e^{−L·T*} [D, sketch]. The sharp rate needs a "zero-constrained uncertainty principle", equivalent to the near-optimality in item 3: open (Connes' §6.6 "remaining step"). Pair correlation does not enter the rate (products superpose independent zero sets and obey it).
6. **ζ's drift is mostly the Fuchs approximation.** Against the exact prolate eigenvalue 1 − λ_n(c), ζ's index excess drops from +0.60/+0.89 to +0.32/+0.26 (even/odd). That is inside the characters' range [−0.01, +0.32] [N].

## 1. Setting and notation

- **Window.** L = log x, W = [−L/2, L/2]. Real f ∈ L²(W) of sector s ∈ {0, 1} means f(−u) = (−1)^s f(u). F(z) = ∫ f(u) e^{izu} du.
- **Weil's form.** 𝓛 = L(s, χ) for a real primitive character χ mod q of parity κ (χ(−1) = (−1)^κ), or ζ (q = 1, κ = 0, a pole r = 1). Write the nontrivial zeros as ρ = ½ + iγ_ρ, γ_ρ ∈ ℂ, |Im γ_ρ| < ½. For f compactly supported and of bounded variation, the explicit formula gives, with no hypothesis on the zeros, Q(f) = Σ_ρ F(γ_ρ) F(−γ_ρ) [S: Weil 1952; `docs/RESEARCH_R1_EXCLUSION.md` §1]. rh2's matrices E ± 2vvᵀ evaluate the other side of the explicit formula (archimedean + primes + pole) and were validated against zero sums (`docs/DECAY_LAW.md`). λ_s(x) is the infimum of Q(f)/‖f‖² in sector s.
- **Decay-law variables.** T* = 2π(x/Q)^{1/d}; c = d·T* (c = 2πx/q in degree 1). Fuchs–Slepian leakage ℓ_n(c) = ½·4√π 8ⁿ c^{n+½} e^{−2c}/n! (rh2's convention, `decay_law_mp._ell`). Exact leakage ℓ^ex_n(c) = ½(1 − λ_n(c)), λ_n(c) the n-th eigenvalue of the time- and band-limiting operator.
- **Fourier on the line.** ψ̂(y) = ∫ ψ(t) e^{−2πity} dt. Hermite functions h_k(y) = H_k(√(2π) y) e^{−πy²} satisfy ĥ_k = (−i)^k h_k. A time-limited prolate ξ_n on [−λ, λ] (c = 2πλ²) satisfies ξ̂_n = μ_n ξ_n on [−λ, λ] with μ_n = (−i)^n √λ_n(c) [S: Connes 2602.04022 §6.3 [R] for the even case].
- **Generalised prolates.** The finite Hankel transform of order ν, (Hφ)(r) = ∫₀¹ J_ν(crρ)√(crρ) φ(ρ) dρ, has eigenvalues γ_{ν,m}, with concentrations λ_{ν,m}(c) = c γ²_{ν,m} [S: Slepian 1964, PSWF IV, not read: U]. ν = −½ and +½ are the even and odd 1D prolates, with n = 2m and 2m + 1. `scripts/research_r2/hankel_prolate.py` computes them from the commuting differential operator, which is tridiagonal in a Jacobi basis. Its values equal Bouwkamp's 1D values to all printed digits (`--check`) [N].

## 2. (a) The conjecture

### 2.1 Statement

Let 𝓛 be a degree-1 L-function above, or a holomorphic newform of weight k and level N (normalised, root number ε), with conductor Q and degree d. Define the **order** ν and the **level** m:

| 𝓛 | Γ-factor | ν | m in sector s |
|---|---|---|---|
| L(s, χ), χ real of parity κ (ε = 1) | Γ_ℝ(s + κ) | κ − ½ | s + 2r (r = 1 for ζ's pole, else 0) |
| weight-k newform | Γ_ℂ(s + (k−1)/2) | k − 1 | 0 if (−1)^s = ε, else 1 |

The equivalent 1D index is n = 2m + ν + ½. That gives n = κ + 2s + 4r in degree 1 (rh2's rule) and n = k − ½ + 2m for GL(2). The leakage is ℓ^ex_{ν,m}(c) = ½(1 − λ_{ν,m}(c)), with c = d·T*.

**Conjecture D (window decay law).** For x/Q beyond the transition (§2.3), in each sector:

- **(D1, rate)** ln λ_s(x) = −2c (1 + o(1)).
- **(D2, index)** ln λ_s(x) = ln ℓ^ex_{ν,m}(c) + o(ln c). Equivalently, the continuous index n* (P-DL11's metric) converges to 2m + ν + ½ as c → ∞.
- **(D3, strong)** There are constants 0 < r₁ ≤ r₂ < ∞, depending on 𝓛 and the sector (the data show R varies with the conductor inside a Γ-class, by about 3.5× among odd characters), with r₁ ≤ R_s(x) := λ_s(x)/ℓ^ex_{ν,m}(c) ≤ r₂.

D3 ⇒ D2 ⇒ D1. **Status:**

- **The upper half of D1:**
  - Theorem 2 with the Kaiser–Bessel trial gives an explicit unconditional bound for every x.
  - Its value is evaluated [N] at 82 grid points (Theorem 4).
  - Its shape C c^p e^{−2c} for all x is [D, sketch].
  - Corollary 3 (Hermite) is the fully elementary version, at rate e^{−c}.
- **The upper halves of D2–D3** are supported numerically, with R ≤ R_trial (§3.6). They are not proved.
- **The lower halves** are open (§4). For ζ they imply RH.

**What "≈" means in rh2's tables.** rh2's measured R is λ_N/ℓ_n, with λ_N a finite-basis upper bound. So a bounded R in the data is consistent with D3 but supports only its upper half. The lower half has direct evidence only at the four certified supports, where R_ex(cert) ≥ 4.59 (even) and ≥ 2.67 (odd) for ζ (§5.3).

### 2.2 Sampling-theory reformulation (under RH)

Under RH, γ_ρ ∈ ℝ, Q(f) = Σ_γ |F(γ)|², and ‖f‖² = (2π)⁻¹‖F‖²_{L²(ℝ)}. Then

  λ_s(x) = 2π · inf { Σ_γ |F(γ)|² / ‖F‖² : F ∈ PW_{L/2}, F(−t) = (−1)^s F(t) },

the **lower sampling constant of the zero set** for the Paley–Wiener space of type L/2. In this language D1 says:

> Under RH the zeros of ζ satisfy a lower sampling inequality in every Paley–Wiener space PW_τ, and its constant is A(τ) = exp(−4π e^{2τ}(1 + o(1))).

This sharpens the qualitative positivity. That is strict for each τ under RH: an F of exponential type that vanishes on a set of density growing like log t is 0 (Jensen) [D], and the minimum is attained (Bombieri, Thm. 3 [S: `docs/LIT_TANGENTS.md`]). The zero set is not a frame (its density grows like log t, so there is no upper Bessel bound), and Landau's density theorem is satisfied with room to spare (D⁻ = ∞). The decay law says how fast the lower constant degenerates as the band grows: doubly exponentially in τ. The cause is that the zeros are sub-Nyquist below T* = 2πe^{2τ}.

### 2.3 Range of validity (the transition)

The law needs the relevant prolate to be in the plateau of its spectrum, where 1 − λ_{ν,m}(c) ≪ 1. By Landau–Widom, about 2c/π eigenvalues are near 1 [S: via Bonami–Jaming–Karoui eq. (2.8) [R]]. So the transition is roughly c ≳ (π/2)(2m + ν + 1).

- **Degree 1** (n ≤ 6): c ≳ 10, i.e. x/q ≳ 1.6. That matches the observed onset at x/q ≈ 2–3.
- **Weight k = 26** (ν = 25): c ≳ 40, i.e. v ≳ 3.2. P-DL7's grid v ∈ [2, 4] straddled this, which explains that kill.

Using the *exact* generalised-prolate leakage, rather than its large-c asymptotic, removes the need for c ≫ ν² (§6.3).

### 2.4 What would falsify it

- **D1 (upper half):** proved in the sense of §2.1, so it can fail only through an error in Theorems 1–2 or in their evaluation. Two checks guard against that: measured ≤ bound at all 326 points, and F_win = −F_out at zeros to 40 digits.
- **D1 (lower half):** a sustained local slope −d ln λ/dc > 2, for example a drift toward Zhu's π (Zhu's Conj. 12.1 gives −ln λ ≈ 2π² x = π c asymptotically [S: `docs/LIT_TANGENTS.md`]). Detectable by upper bounds alone: a converged finite-basis λ_N whose ratio to ℓ^ex falls exponentially in c refutes it.
- **D2:** the continuous index n* moving away from 2m + ν + ½ as c grows. Concretely, for level-1 forms on v ∈ [14, 18] it must move toward the exact-Hankel values of §6.3 (registered in §7).
- **D3:** R_ex leaving every fixed band, e.g. a registered extension with R_ex(x/q = 32)/R_ex(16) outside [¼, 4]. Or any trial function with Q/‖f‖² < r₁·ℓ^ex, with r₁ fixed in advance from the certified points.
- **Any of D1–D3 for a function with an Euler product and RH-verified zeros** (all rh2 objects): a negative λ_N. That would be a witness against RH.

## 3. (b) The upper bound

### 3.1 The E-map (Connes–Consani; Connes 2026 §6.4)

- **The map.** Connes writes E(f)(u) = u^{1/2} Σ_{n≥1} f(nu) [R: 2602.04022 eq. (12)]. Riemann's k = E(h), with h the combination of h₀ and h₄ (Hermite) that has vanishing integral, and Ξ is its Fourier transform [R: Fact 6.2, eq. (13)].
- **The radical.** The range of E lies in the radical of the global Weil form; Poisson with f(0) = f̂(0) = 0 gives E(f̂)(x) = E(f)(1/x) [R: eq. (18)].
- **The near-radical vector.** k_λ := E(h_λ) restricted to [λ⁻¹, λ], with h_λ built from prolates. QW_λ takes "non-zero, but extremely small values" on it [R: §6.4]. Its convergence to Ξ is Fact 6.4 [R], and the comparison of k_λ with the minimiser is listed as a remaining step [R: §6.6].
- **What is new here.** Connes gives no inequality. The following turns the construction into an exact identity and an explicit bound, for twisted E-maps as well.

### 3.2 Theorem 1 (the E-map identity) [D]

**Hypotheses (class 𝒜_χ(μ)).** ψ: ℝ → ℝ satisfies:

- **(A1)** ψ is continuous and locally of bounded variation, with ψ(−y) = (−1)^κ ψ(y);
- **(A2)** |ψ(y)| ≤ C(1 + |y|)^{−1−δ} for some C, δ > 0;
- **(A3)** ψ̂ = μψ (necessarily μ² = (−1)^κ);
- **(A4)** if q = 1, ψ(0) = 0.

Put E_χψ(w) = w^{1/2} Σ_{n≥1} χ(n) ψ(nw/√q), f_ψ(u) = E_χψ(e^u), and σ = (i^κ μ)⁻¹ ∈ {±1}.

**Theorem 1.**

- **(i)** f_ψ(−u) = σ f_ψ(u): f_ψ lies in the sector with (−1)^s = σ.
- **(ii)** |f_ψ(u)| ≤ C′ e^{−(½+δ)|u|}, so F_ψ is analytic on |Im z| < ½ + δ.
- **(iii)** There, F_ψ(z) = q^{(½+iz)/2} L(½ + iz, χ) M_ψ(z), with M_ψ(z) = ∫₀^∞ ψ(y) y^{−½+iz} dy. Hence **F_ψ(γ_ρ) = 0 for every nontrivial zero**, on or off the line, of any multiplicity.
- **(iv)** For every L > 0, with f_win = f_ψ·1_W, f_out = f_ψ − f_win and F_out the transform of f_out:

  **Q(f_win) = σ Σ_ρ F_out(γ_ρ)².**

**Proof.**

*(i)* For primitive χ mod q and continuous g with g, ĝ = O((1 + |y|)^{−1−δ}), Σ_{n∈ℤ} χ(n) g(n) = (τ(χ)/q) Σ_{m∈ℤ} χ̄(m) ĝ(m/q) [S: standard twisted Poisson summation]. Take g(y) = ψ(yw/√q), so ĝ(ξ) = (√q/w) μ ψ(ξ√q/w). Use χ real, τ(χ) = i^κ√q [S: Gauss sums], and the parity of χ(n)ψ(n·) (each side is twice its sum over n ≥ 1). This gives Σ_{n≥1} χ(n)ψ(nw/√q) = (i^κ μ/w) Σ_{m≥1} χ(m) ψ(m/(w√q)), i.e. E_χψ(w) = i^κ μ · E_χψ(1/w). For q = 1, Poisson leaves the extra term (μψ(0)/w − ψ(0))/2, which (A4) removes.

*(ii)* For w ≥ 1, |E(w)| ≤ C q^{(1+δ)/2} ζ(1+δ) w^{−½−δ}. For w ≤ 1, use (i).

*(iii)* On −½ − δ < Im z < −½ the Dirichlet series converges absolutely and ∫|ψ(y)| y^{−½−Im z} dy < ∞, so Fubini gives the product. M_ψ is analytic on −½ − δ < Im z < ½, since ψ is bounded near 0. L(½ + iz, χ) is entire for q > 1. Both sides are analytic on the strip, hence equal. For q = 1 the pole of ζ at z = −i/2 meets M_ψ(−i/2) = ∫₀^∞ ψ = ψ̂(0)/2 = μψ(0)/2 = 0. Nontrivial zeros have |Im γ_ρ| < ½.

*(iv)* f_win has compact support and bounded variation, so Q(f_win) = Σ_ρ F_win(γ_ρ) F_win(−γ_ρ), absolutely convergent [S: explicit formula]. F_out converges on |Im z| < ½ + δ, and F_win = F_ψ − F_out there. By (iii), F_win(γ_ρ) = −F_out(γ_ρ), and f_out's symmetry gives F_out(−z) = σ F_out(z). ∎

**Remarks.**

- **Euler-blind.** No multiplicativity is used. (i)–(iv) hold verbatim for any real-coefficient Dirichlet series with the same Poisson relation. Davenport–Heilbronn is one: its E-map vector stays near-radical, and positive, past DH's crossover, exactly as `docs/CONNES_LETTER.md` step 5 observed (QW(k_λ) = +3.7e-39 at x = 40, where DH's ε₀ = −0.011) [S].
- **Check.** For the Hermite trials of §3.4 at x = 4, F_win(γ) + F_out(γ) vanishes at the first three zeros of ζ (both sectors), L(χ₋₄) and L(χ₅). The largest residual is 1.1e-40 at 40 digits, while F_win itself is 10⁻³ to 0.5 [N] (`check_identity.py` → `data/research_r2/check_identity.json`).

### 3.3 Theorem 2 (an unconditional bound) [D]

Let ψ ∈ 𝒜_χ(μ). Assume f_ψ is absolutely continuous on [L/2, ∞) with ∫ |f_ψ′| e^{u/2} < ∞. Put

  B = 2∫_{L/2}^∞ |f_ψ(u)| e^{u/2} du,  A = 2(|f_ψ(L/2)| e^{L/4} + ∫_{L/2}^∞ |f_ψ′(u)| e^{u/2} du).

Then |F_out(γ_ρ)| ≤ min(B, A/|Im ρ|) for every nontrivial zero, and

  **λ_s(x) ≤ Q(f_win)/‖f_win‖² ≤ (2A²/‖f_win‖²) ∫_{A/B}^∞ N(t) t⁻³ dt,**

with N(t) = #{ρ: |Im ρ| ≤ t}, counted with multiplicity. Any upper bound for N may be inserted.

*Proof.* F_out(z) = ∫_{L/2}^∞ f(u)(e^{izu} + σe^{−izu}) du. For |Im z| ≤ ½ and u > 0, |e^{±izu}| ≤ e^{u/2}; this is the only property of the zeros used, namely 0 < Re ρ < 1. That gives |F_out| ≤ B. Integrating by parts, with the boundary term at ∞ vanishing by (ii), gives |F_out(z)| ≤ A/|z|. Also |γ_ρ| ≥ |Re γ_ρ| = |Im ρ|. So |Q(f_win)| ≤ Σ_ρ m(|Im ρ|) with m = min(B², A²/t²). The Stieltjes sum ∫ m dN equals 2A² ∫_{A/B}^∞ N t⁻³ dt, since m is constant below A/B. ∎

The right side increases in both A and B, so upper bounds for A and B may replace them.

**Explicit zero counts used** (both [A]):

- ζ: N(t) = 2N_ζ(t), with N_ζ(t) ≤ θ(t)/π + 1 + 0.112 log t + 0.278 log log t + 2.510 (Trudgian 2014, J. Number Theory 134, for t ≥ e). N_ζ = 0 below 14.134.
- L(s, χ): |N(t, χ) − (t/π) log(qt/2πe) + χ(−1)/4| ≤ 0.22737ℓ + 2 log(1 + ℓ) − 0.5, with ℓ = log(q(t+2)/2π) > 1.567 and t ≥ 5/7 (Bennett–Martin–O'Bryant–Rechnitzer, arXiv:2005.02989). N is nondecreasing, which covers smaller t.

### 3.4 Corollary 3: Hermite trial functions (Riemann's k), fully explicit, rate e^{−c} [D + N]

Take ψ = h_n with n = κ + 2s. For ζ, take ψ = H_{n+4}(0) h_n − H_n(0) h_{n+4} with n = 2s, which has ψ(0) = 0.

- For ζ's even sector this is Riemann's h ∝ y²(2πy² − 3)e^{−πy²}. f_ψ is Riemann's k(e^u), and F_ψ ∝ Ξ.
- (A1)–(A4) and Theorem 2's hypotheses hold, because of the Gaussian decay. A, B and ‖f_win‖² are one-dimensional integrals of explicit functions (`hermite_bound.py`).
- The tails decay like e^{−πy²}, so the bound decays like c^p e^{−c} with c = 2πx/q, which is half the conjectured rate.

**Table 3.4** (`data/research_r2/hermite_bound.json`; measured values are rh2's best finite-basis minima, from scanhi, scanext and scan8):

| object | sector | n | x/q | measured/bound (min .. max) | p in bound ≈ C c^p e^{−c} |
|---|---|---|---|---|---|
| ζ | even | 4 | 2..10 | 2.6e-24 .. 4.2e-3 | 4.45 |
| ζ | odd | 6 | 2..10 | 3.8e-23 .. 2.9e-2 | 6.75 |
| χ₅, χ₈, χ₁₂, χ₁₃, χ₁₇ | even | 0 | 2..16 (2..10) | ≥ 1.8e-43 .. ≤ 1.2e-5 | 0.19–0.23 |
| same | odd | 2 | same | ≥ 2.8e-42 .. ≤ 1.9e-4 | 2.17–2.22 |
| χ₋₃, χ₋₄, χ₋₇, χ₋₂₀, χ₋₈, χ₋₁₁, χ₋₁₅ | even | 1 | 2..16 (2..10) | ≥ 4.7e-43 .. ≤ 6.2e-5 | 1.15–1.25 |
| same | odd | 3 | same | ≥ 7.3e-42 .. ≤ 9.8e-4 | 3.19–3.30 |

All 244 points: measured/bound ∈ [1.8e-43, 2.9e-2], every ratio ≤ 1. The ratio falls like e^{−c}, because the bound's rate is half the measured one.

**Theorem 1 + 2 check.** At eight points (ζ and χ₅, x/q = 2, 4, both sectors) the trial functions' Rayleigh quotients in rh2's basis lie below the bound by factors 12–28 (`hermite_bound.py --matrix` → `data/research_r2/hermite_matrix_check.json`) [N].

For every x this is an explicit, unconditional bound. Its shape is λ_s(x) ≤ C c^p e^{−2πx/q}, with fitted p ≈ n + ¼ for the characters (4.45 and 6.75 for ζ) [N, and an elementary Gaussian-tail sketch]. It is much stronger than Zhu's RH-conditional bound λ* ≤ exp(−L e^L), where Zhu's L is the half-width a = ½ log x, so the bound is exp(−a√x) [S: `docs/LIT_TANGENTS.md`]: e^{−2πx} against e^{−(√x log x)/2}, and it needs no hypothesis. The e^{−c} barrier is the same one Bonami–Jaming–Karoui meet for 1 − λ_n(c) itself. Their explicit non-asymptotic bound λ_n(c) ≥ 1 − (7/√c)(2c)ⁿe^{−c}/n! (n ≤ c/2.7, Thm 3.2) also comes from Hermite trial functions [R].

### 3.5 Theorem 4: Kaiser–Bessel trial functions, rate e^{−2c} [D + N]

**The trial functions** (`kb_bound.py`). Let w(t) = (1 − t²)^{m/2} I_m(β√(1 − t²)) on [−1, 1], with m = 2 and β ∈ {c, c − 2}, the better kept. Set λ = √(x/q) and

- φ(y) = p(y/λ) w(y/λ) for κ = 0, with p(t) = 1 + a t²; a = 0 except for ζ, where a is fixed by ψ(0) = 0;
- φ(y) = (y/λ) w(y/λ) for κ = 1;
- ψ = (φ + μ̄ φ̂)/2, with μ = (−i)^{κ+2s}.

**Closed forms.** ∫₋₁¹ w e^{−iωt} dt = √(2π) β^m K_{m+½}(ω), with K_ν(ω) = I_ν(z)/z^ν for z = √(β² − ω²), or J_ν(s)/s^ν for s = √(ω² − β²) [S: Lewitt 1990, JOSA A 7; verified here by quadrature to 12 digits, including the t² and t moments [N]]. The derivatives follow from d/dω K_ν = −ωK_{ν+1}.

**Envelopes.** Beyond the band, |ψ| and |ψ′| are bounded by envelopes that use two facts:

- |J_ν(s)/s^ν| ≤ 1/(2^ν Γ(ν+1)) [S: DLMF 10.14.4];
- Hankel's terminating expansion for half-integer order, |J_{k+½}(s)| ≤ √(2/(πs)) Σ_{j≤k} (k+j)!/(j!(k−j)! 2^j) s^{−j} [D from the finite expansion].

**Bounds on A and B.** Write env₀ ≥ |ψ| and env₁ ≥ |ψ′| beyond the band. Then B ≤ 2√q Σ_n (1/n) ∫_{nλ}^∞ env₀ and A ≤ 2[√x Σ_n env₀(nλ) + B/4 + √q Σ_n (1/n) ∫_{nλ}^∞ y·env₁]. The n-tails beyond n = 40 are bounded using the envelopes' y^{−m−1} decay (`bounds_AB`).

**Rigour.** Three ingredients are floating point [N]: ‖f_win‖² (composite Gauss–Legendre at 25 digits), the envelope integrals (mpmath quadrature), and the monotonicity of env·y^{m+1} used for the n-tails. That monotonicity is asserted in the code on an 80-point geometric grid spanning y/(40λ) ∈ [1, 2^{20}], not proved. Every other step is an inequality. A ball-arithmetic pass would make each table entry a theorem.

**Table 3.5** (`data/research_r2/kb_bound.json`, `kb_bound_b.json`; m = 2, the better of β = c and c − 2):

| object | sector | n | x/q | measured/bound (min .. max) | p in bound ≈ C c^p e^{−2c} |
|---|---|---|---|---|---|
| ζ | even | 4 | 2..10 | 2.8e-6 .. 2.0e-3 | 9.4 |
| ζ | odd | 6 | 2..10 | 1.2e-5 .. 4.8e-3 | 11.8 |
| χ₅ / χ₈ | even | 0 | 2..16 | 3.2e-7 .. 1.3e-4 | 3.5 |
| χ₅ / χ₈ | odd | 2 | 2..16 | 1.2e-7 .. 1.7e-3 | 7.2–7.3 |
| χ₋₃ / χ₋₄ / χ₋₇ / χ₋₂₀ | even | 1 | 2..16 | 6.3e-8 .. 1.6e-4 | 5.0–5.1 |
| χ₋₃ / χ₋₄ / χ₋₇ / χ₋₂₀ | odd | 3 | 2..16 | 1.6e-7 .. 2.0e-3 | 7.8–7.9 |

All 82 points: measured/bound ∈ [6.3e-8, 4.8e-3], every ratio ≤ 1.

**Reading.**

- Over the grid the measured minimum falls by about 50 (ζ, x = 2 → 10) to 75 (characters, x/q = 2 → 16) orders of magnitude, while measured/bound changes by about 3. So the bound has the measured exponential rate, and the gap is the polynomial c^{p − n − ½}.
- **Logic of the comparison.** Both columns are upper bounds on the true λ_s(x). So measured ≤ bound is a consistency requirement, not a logical necessity: a violation would mean the finite basis is unconverged. It holds everywhere.

**Theorem 4 (statement).** For every tested (𝓛, x, s), λ_s(x) is at most the tabulated bound. That holds unconditionally, modulo the quadrature of ‖f_win‖².

**Asymptotic form (sketch [D]).** A and B are O(β^{m+O(1)}): the out-of-band transform is algebraic, with no exponential factor. ‖f_win‖² ≥ const · I_m(β)² c^{−2j}, where j counts the projections the class and pole impose. With β = c − O(1) this gives λ_s(x) ≤ C c^p e^{−2c}. The fitted p is in the table: it runs from 3.5 to 12 against the conjectured n + ½ (0.5–6.5). That polynomial loss is the price of a smooth edge (§3.6).

### 3.6 The sharp construction: time-limited prolates [N, with a sketch [D]]

Replace φ by the time-limited prolate ξ_n on [−λ, λ] (c = 2πλ²), with n = κ + 2s. For ζ, use Connes' combination β₀^{(n+4)} ξ_n − β₀^{(n)} ξ_{n+4}, which has ∫φ = 0. Project the E-map window function onto the sector and evaluate its Rayleigh quotient with rh2's explicit-formula matrix in the scan's own N-mode basis (`prolate_trial.py`). This is an unconditional upper bound for λ_s(x): no zeros enter.

**Table 3.6** (`data/research_r2/prolate_trial_{a,b,c}.json`; R = λ/ℓ_n(c) with Fuchs' ℓ_n):

| object | sector | n | x/q | trial/measured | R_trial | R_measured |
|---|---|---|---|---|---|---|
| ζ | even | 4 | 2..8 | 1.07 .. 1.49 | 2.21 .. 11.0 | 2.06 .. 7.41 |
| ζ | odd | 6 | 2..8 | 1.20 .. 1.32 | 0.44 .. 5.22 | 0.37 .. 3.94 |
| χ₅ | even / odd | 0 / 2 | 2..8 | 1.14 .. 1.66 / 1.10 .. 1.39 | 27.9 .. 53.2 / 4.2 .. 7.6 | 22.8 .. 34.7 / 3.8 .. 5.6 |
| χ₈ | even / odd | 0 / 2 | 2..8 | 1.17 .. 2.07 / 1.15 .. 1.34 | 17.9 .. 40.9 / 3.5 .. 5.9 | 15.4 .. 19.8 / 3.0 .. 4.8 |
| χ₋₃ | even / odd | 1 / 3 | 2..8 | 1.11 .. 1.46 / 1.08 .. 1.27 | 9.5 .. 20.2 / 2.9 .. 6.3 | 8.6 .. 13.9 / 2.7 .. 5.0 |
| χ₋₄ | even / odd | 1 / 3 | 2..8 | 1.05 .. 1.34 / 1.07 .. 1.28 | 7.3 .. 13.1 / 2.9 .. 5.7 | 6.7 .. 9.8 / 2.7 .. 4.6 |
| χ₋₇ | even / odd | 1 / 3 | 2..8 | 1.25 .. 2.37 / 1.16 .. 2.48 | 4.1 .. 7.1 / 2.3 .. 4.7 | 3.1 .. 4.5 / 1.7 .. 3.1 |
| χ₋₂₀ | even / odd | 1 / 3 | 2..8 | 1.06 .. 1.38 / 1.12 .. 1.34 | 3.5 .. 5.7 / 4.4 .. 7.6 | 3.3 .. 4.1 / 3.9 .. 6.1 |
| χ₋₈, χ₋₁₁, χ₋₁₅ (fresh, P-DL8) | even / odd | 1 / 3 | 2..6 | 1.06 .. 1.51 / 1.16 .. 1.51 | | |
| χ₁₂, χ₁₃, χ₁₇ (fresh, P-DL8) | even / odd | 0 / 2 | 2..6 | 1.05 .. 1.76 / 1.10 .. 1.64 | | |

All 106 points: trial/measured ∈ [1.05, 2.48], median 1.28, 92 of 106 below 1.5.

**Reading.**

- The E-map prolate trial of the predicted index is within a factor of 1.05–2.5 of the minimum rh2 measured in the same basis.
- So the near-minimiser is, to that accuracy, E(ξ_n) truncated. This is the twisted, two-sector generalisation of `docs/CONNES_LETTER.md`'s k_λ result for ζ even (1.19–1.48).
- The larger ratios (up to 2.48, χ₋₇ at x/q = 4) are not monotone in x/q. They come from the quadrature of the projection: any coefficient vector is a legitimate trial, so these remain valid upper bounds, merely less sharp.

**Why Theorem 2 does not yet give the sharp constant.** The time-limited prolate jumps at ±λ. So ψ ~ (edge value)·sin(2πλy)/y out of band, E(ψ) carries a sawtooth Σ sin(nα)/n beyond the window, and B = ∞: ∫ |f_out| e^{u/2} diverges logarithmically at the edge of the critical strip.

- **Fix (a):** use RH-verified zeros below T₀ = 3·10¹² [S: Platt–Trudgian 2021, via `docs/RESEARCH_R1_EXCLUSION.md`]. The weights e^{u/2} then disappear for every relevant zero. But the total variation of the sawtooth also diverges, so Theorem 2's integration by parts must be replaced by an L² (Gallagher-type) count.
- **Fix (b):** taper the edge. That loses c^{m+½}, as in Theorem 4.

The sharp upper bound λ_s ≤ C·L^a·ℓ^ex_n(c) is therefore [D, sketch only, constant not established]. Numerically, λ_s ≤ R_trial·ℓ_n(c) on the grid, with Fuchs' ℓ_n and R_trial ≤ 53.2 (χ₅ even), as a finite-basis Rayleigh quotient [N].

## 4. (c) The lower bound: what hypothesis gives it?

### 4.1 Necessary conditions [D]

1. **It implies RH.** If λ_odd(x) > 0 for every x, then Q ≥ 0 on all odd test functions, and Yoshida's Prop. 1 gives RH [S: Suzuki arXiv:2606.09096 §1.1, via `docs/LIT_TANGENTS.md`]. So the lower half of D1 for ζ, for all x, is at least as strong as RH. For a single x it is a finite-window positivity statement; rh2 has certified those at four supports (§5.3).
2. **The functional equation plus the smooth zero density are not enough.** DH has conductor 5, Γ_ℝ(s+1), a functional equation and the smooth zero density of an odd character mod 5.
   - By Theorem 1's remark, its E-map vectors satisfy the *same* upper bound.
   - Yet λ_DH changes sign at x ≈ 30.8 (certified bracket [10.805, 30.745) for the first failure, `docs/CONTROL_CERTIFICATES.md` [S]), through its quadruple at height 85.7 (R1's mirror test σ₁ = 1.0002 at 30.745 [S]).
   - So any hypothesis that yields the lower bound must exclude off-line zeros, at least those within the window's reach (R1 §4: a few times 10³–10⁶ in height for these supports [S]).
3. **Exact quantisation by zero counts is not right.** The zero-deficit heuristic λ ≈ e^{−4πN_def} with the *actual* count N(T*) would jump by e^{4π} each time T* passes a zero. That is not seen (`docs/DECAY_LAW.md`, "Interpretation"). So a correct lower bound depends on the zeros only through averaged density, which also makes fine spacing statistics unlikely to matter (§4.4).

### 4.2 The sampling problem under RH [D]

By §2.2, λ_s(x)/2π is the lower sampling constant of Γ = {γ} for PW_{L/2}. Its two ingredients pull in different directions.

- **Above T*** the zeros are super-Nyquist. Unconditionally #(Γ ∩ [t, t+R]) = (R/2π) log(t/2π) + O(log t) [S: Riemann–von Mangoldt]; under RH the error is O(log t / log log t) [S]. So in long windows their density exceeds L/2π.
- **Landau 1967** (Acta Math. 117) [S: standard; not re-read]. A sampling set for PW_τ needs D⁻ ≥ τ/π. Γ satisfies this trivially (D⁻ = ∞), so Landau's theorem gives positivity-type information only, no constant.
- **Beurling's sufficiency theorem** (separated sets with D⁻ > τ/π are sampling) [S] and the **Ortega-Cerdà–Seip characterisation** of Fourier frames through de Branges spaces (Ann. Math. 155 (2002) 789–806) [A]. Both need uniform separation, which ζ's zeros do not have: their spacing tends to 0 and simplicity is unproved.
- **Below T*** the zeros are sub-Nyquist. A function of type L/2 can nearly vanish on them, with about 2·N_def = 2dT*/2π spare degrees of freedom (Connes–Consani count about 2x small eigenvalues [S: `docs/LIT_TANGENTS.md`]).

### 4.3 A conditional lower bound with a weaker rate [D, sketch]

**Hypothesis LSH(T₁, ε).** Above T₁ = (1+ε)T* the zeros contain a uniformly separated subsequence Γ′ whose lower Beurling density relative to [T₁, ∞) exceeds (1+ε)L/2π. A Beurling–Landau sampling inequality then holds on the half-line, with the boundary leakage into t < T₁ controlled.

**Claim (sketch).** Under RH + LSH, λ_s(x) ≥ A(ε) · (1 − λ₀(L·T₁/2)) · (1 − o(1)) ≈ exp(−(1+ε) L·T*) up to polynomial factors.

*Sketch.* Q(f) ≥ Σ_{γ∈Γ′} |F(γ)|² ≥ A ∫_{|t|>T₁} |F|² − (boundary). The time–frequency concentration of PW_{L/2} on [−T₁, T₁] in the (u, t) plane has parameter c′ = (L/2)T₁, so ∫_{|t|>T₁} |F|² ≥ (1 − λ₀(c′))‖F‖² [S: Slepian–Pollak].

**Status.**

- The rate is e^{−2πx log x} for ζ, weaker than the conjectured e^{−4πx} by the factor log x/2 in the exponent.
- LSH itself is not known even under RH, because separation of zeros is not known. GUE heuristics make it overwhelmingly plausible [U].
- A Turán–Nazarov route, with mesh conditions on [T₁, 2T₁] under RH, would give the same rate class with unspecified constants [U].

### 4.4 What the sharp rate needs: a zero-constrained uncertainty principle [D, conjecture]

**ZUP.** If F ∈ PW_{L/2} satisfies Σ_{γ<T*} |F(γ)|² ≤ η‖F‖², then F carries at least c′ e^{−2c} poly(c)·‖F‖² of energy on the zeros above T*.

Theorem 1 shows that truncated E-images, whose transforms are ζ·M_ψ up to the leakage, *attain* the upper half. ZUP (with Theorems 2–4) is therefore a restatement of the claim that **no function beats the E-map prolate construction by more than a polynomial factor**. In the tested bases nothing beats it by more than 1.05–2.48 (§3.6), but a basis minimum is not the true infimum. That is Connes' remaining step "k_λ is a sufficiently good approximation of θ_x" [R: §6.6], and for ζ, uniformly in x, it implies RH. The natural framework is the de Branges space attached to ξ under RH (Lagarias, arXiv:math/0601653 [A]), whose spectrum is the zeros; or Burnol's Sonine spaces [A, via CCM]. I found no theorem there that gives the constant.

**Pair correlation is irrelevant to the rate** [E + D]:

- The rate e^{−2dT*} holds for products (P-DL5, P-DL10), whose zero sets are superpositions of *independent* sets with no mutual repulsion.
- DH's minimum stays inside the odd-character band until it fails.
- The rate is fixed by c = dT*, a function of the smooth density alone.

Spacing statistics could affect R (the constant), which the data cannot separate from conductor effects.

### 4.5 Literature checked (labels)

| source | status | what it gives here |
|---|---|---|
| Connes, arXiv:2602.04022 §6.2–6.6 | [R] | E-map, Poisson (18), k_λ, near-radical in words, remaining steps |
| Connes–Consani, arXiv:2106.01715 | [A] | small eigenvalues from prolates, numerically |
| CCM, arXiv:2310.18423 Prop. 5.1 | [S: `docs/LIT_TANGENTS.md`] | metaplectic W_λ; weight ½ and 3/2 pieces |
| Fuchs 1964 | [S: BJK eq. (2.9)] | 1 − λ_n(c) ∼ 4√π 8ⁿ c^{n+½} e^{−2c}/n! |
| Bonami–Jaming–Karoui, arXiv:1804.01257 | [R: Thm 3.2, eqs. (2.8)–(2.9)] | explicit λ_n ≥ 1 − (7/√c)(2c)ⁿe^{−c}/n!; Landau–Widom count |
| Landau–Widom 1980 | [S: BJK] | plunge width (2/π²) log((1−ε)/ε) log c |
| Slepian 1964 (PSWF IV) | [U] | generalised prolates (computed here directly, §1) |
| Landau 1967; Beurling | [S] | necessary density; sufficiency for separated sets |
| Ortega-Cerdà–Seip 2002 | [A] | sampling sequences for PW via de Branges spaces |
| Lagarias, math/0601653 | [A] | de Branges spaces from L-functions under RH |
| Trudgian 2014; Bennett–Martin–O'Bryant–Rechnitzer 2021 | [A] | explicit N(T), N(T, χ) |
| Lewitt 1990 | [S, verified [N]] | Kaiser–Bessel transforms |
| Zhu, arXiv:2608.24827 | [S: `docs/LIT_TANGENTS.md`] | Conj. 12.1; Thm 1.3 (RH ⇒ λ* ≤ exp(−L e^L)) |

## 5. (d) The degree-1 index rule from the construction

### 5.1 Derivation [D]

Theorem 1 forces three properties on an admissible ψ:

1. **Parity κ.** The twisted Poisson relation needs ψ(−y) = χ(−1)ψ(y). The Mellin transform of y^κ e^{−πy²} is the Γ_ℝ(s + κ) factor. In CCM's metaplectic picture, even ψ live in the lowest-weight-½ piece L²(ℝ)_ev and odd ψ in the weight-3/2 piece L²(ℝ)_odd [S: CCM Prop. 5.1]. **That is why an odd character starts one index higher:** its lowest admissible Fourier eigenfunction is h₁ (or ξ₁), not h₀.
2. **Fourier class ↔ sector.** E_χψ(1/w) = (i^κμ)⁻¹ E_χψ(w). With μ = (−i)^n this is (−1)^{(n−κ)/2}, so the sector is s with n ≡ κ + 2s (mod 4).
3. **The pole removes one class member.** For ζ the Poisson constant terms (μψ(0)/w − ψ(0))/2 are the pole at s = 1 and the value ζ(0) = −½. They make E(ψ) grow like e^{−u/2} beyond the window unless ψ(0) = 0 (and then ψ̂(0) = 0).
   - Functions of class μ that are concentrated in both time and frequency lie, up to their leakage, in the span of the plateau prolates (about 2c/π of them). On that span, evaluation at 0 is a bounded functional.
   - So ψ(0) = 0 removes one dimension. The least-leakage admissible ψ is dominated by ξ_{n+4}, with leakage ≈ ℓ_{n+4}.
   - (In L² alone the constraint would be free, through a narrow spike at 0. But ψ̂ = μψ turns the spike into a broad tail that is O(1) in the window. This is why the constraint costs exactly one class member, as in Connes' k_λ, which uses h₀ and h₄.)

Hence the lowest admissible index in sector s is **n = κ + 2s + 4r**. Lower indices are either in the wrong sector or violate (A4). Higher ones leak more, by a factor of order c^{n′−n}.

**What this proves and what it does not.** It shows that the rule is the index of the best E-map trial function, an upper-bound statement. It is not a proof that the true minimiser is an E-image; that is the ZUP of §4.4.

### 5.2 The numbers [N]

- **Index check.** The trial of exactly this index reaches 1.05–2.5 times the measured minimum (§3.6).
- **Against the exact leakage.** Measured against ℓ^ex rather than the Fuchs asymptotic, the characters' continuous index is n*_ex − n ∈ [−0.01, +0.32]; with Fuchs it is [0.01, 0.45]. The R_ex drifts over x/q ≥ 3 are 0.86–1.59.

**Table 5.2** (`data/research_r2/exact_leakage.json`; scanhi + scanext + scan8; drift = R(last x/q)/R(x/q = 3); n* fitted over x/q ≥ 3):

| object | sector | n | R_F range | R_F drift | R_ex range | R_ex drift | n*_F | n*_ex |
|---|---|---|---|---|---|---|---|---|
| χ₅ | even / odd | 0 / 2 | 22.8–38.5 / 3.8–6.4 | 0.98 / 1.14 | 23.6–39.1 / 4.3–6.6 | 0.96 / 1.00 | 0.01 / 2.17 | −0.00 / 2.10 |
| χ₈ | even / odd | 0 / 2 | 15.4–20.8 / 3.0–5.4 | 1.22 / 1.35 | 15.9–21.0 / 3.8–5.6 | 1.20 / 1.19 | 0.03 / 2.16 | 0.02 / 2.09 |
| χ₁₂ | even / odd | 0 / 2 | 13.0–15.4 / 4.1–5.9 | 1.10 / 1.36 | 13.5–15.6 / 5.1–6.1 | 1.08 / 1.22 | 0.08 / 2.22 | 0.06 / 2.13 |
| χ₁₃ | even / odd | 0 / 2 | 12.1–17.8 / 3.8–6.2 | 0.88 / 1.46 | 12.4–18.1 / 4.7–6.8 | 0.86 / 1.31 | 0.01 / 2.25 | −0.00 / 2.16 |
| χ₁₇ | even / odd | 0 / 2 | 4.6–7.3 / 1.8–3.2 | 1.56 / 1.45 | 4.8–7.3 / 2.2–3.4 | 1.54 / 1.30 | 0.33 / 2.24 | 0.32 / 2.15 |
| χ₋₃ | even / odd | 1 / 3 | 8.6–16.2 / 2.7–6.0 | 1.22 / 1.80 | 9.5–16.7 / 4.4–6.3 | 1.16 / 1.41 | 1.11 / 3.35 | 1.08 / 3.21 |
| χ₋₄ | even / odd | 1 / 3 | 6.7–11.3 / 2.7–5.5 | 1.34 / 1.37 | 7.4–11.5 / 4.4–5.9 | 1.27 / 1.08 | 1.15 / 3.21 | 1.11 / 3.08 |
| χ₋₇ | even / odd | 1 / 3 | 2.1–4.8 / 1.6–3.7 | 1.41 / 1.25 | 2.2–4.9 / 1.9–4.0 | 1.34 / 0.98 | 1.32 / 3.31 | 1.29 / 3.17 |
| χ₋₈ | even / odd | 1 / 3 | 4.8–5.8 / 3.1–5.4 | 1.10 / 1.44 | 5.0–6.0 / 4.9–5.9 | 1.05 / 1.16 | 1.07 / 3.23 | 1.03 / 3.06 |
| χ₋₁₁ | even / odd | 1 / 3 | 4.9–7.8 / 1.9–4.2 | 1.02 / 1.97 | 5.4–8.0 / 2.9–4.6 | 0.97 / 1.59 | 1.03 / 3.45 | 0.99 / 3.28 |
| χ₋₁₅ | even / odd | 1 / 3 | 2.0–3.1 / 2.5–4.0 | 1.39 / 1.17 | 2.1–3.2 / 3.9–4.5 | 1.33 / 0.95 | 1.22 / 3.17 | 1.18 / 3.00 |
| χ₋₂₀ | even / odd | 1 / 3 | 3.3–4.5 / 3.9–6.9 | 1.19 / 1.67 | 3.7–4.6 / 5.5–7.2 | 1.13 / 1.31 | 1.12 / 3.29 | 1.09 / 3.15 |
| **ζ** | even / odd | 4 / 6 | 2.1–8.1 / 0.37–4.6 | **1.92 / 2.86** | 4.5–9.8 / 2.5–6.2 | **1.35 / 1.32** | **4.60 / 6.89** | **4.32 / 6.26** |

Characters: n*_F − n ∈ [0.01, 0.45], n*_ex − n ∈ [−0.01, 0.32].

### 5.3 ζ's drift [N]

**Table 5.3** (ζ's best upper bounds from `data/connes/decay/fit_zeta.json`; certified lower bounds from `docs/CERTIFICATE_238.md`, via `docs/BRIEF_RESEARCH_HANDOFF.md` §3.1):

| sector | n | x | R_F (upper bound) | R_ex (upper bound) | R_ex (certified lower bound) |
|---|---|---|---|---|---|
| even | 4 | 4.95, 7.39, 10.80, 13.46, 19.89 | 5.55, 6.73, 8.89, 9.63, 12.77 | 7.31, 8.05, 10.03, 10.61, 13.62 | 4.59, −, 6.74, 6.70, 7.61 |
| odd | 6 | same | 2.55, 3.68, 5.03, 5.83, 9.05 | 4.53, 5.34, 6.46, 7.11, 10.34 | 2.67, −, 4.49, 4.52, 5.77 |

- **The drift shrinks against the exact leakage.** ζ's even ratio to ℓ₄ drifts by 2.3× over x = 5 → 20 when measured against the Fuchs asymptotic (5.55 → 12.8). Against the exact 1 − λ₄(c) it drifts 1.86× (7.3 → 13.6, upper bounds), or 1.66× (4.59 → 7.61) using the certified lower bounds. The odd sector goes from 3.6× to 2.3× (cert 2.2×). On the scanhi grid (x = 2–10) the continuous index drops from n*_F = 4.60/6.89 to n*_ex = 4.32/6.26.
- **The cause** is that Fuchs' formula is accurate for n = 4, 6 only at large c: (1 − λ₄)/Fuchs is 0.32 at c = 10 and 0.71 at c = 25 (`hankel_prolate.py --check`).
- **The residual excess** (+0.3) is the same size as the characters', and the E-map trial reproduces it (R_trial/R_meas 1.07–1.49 for ζ). So it is a property of the leakage-to-Q map (the zero sum weighting the edge tail), not a different minimiser.
- **Certified points.** The certified lower bounds give R_ex(cert) ∈ [4.59, 7.61] (even) and [2.67, 5.77] (odd). These are the only *lower*-bound evidence for D3 in rh2. They are consistent with r₁ ≈ 4.5 (even) and 2.6 (odd) for ζ over x ∈ [4.95, 19.9].

## 6. (e) GL(2): Voronoi, Hankel prolates, root number, level offset

### 6.1 Construction [D, sketch]

For a weight-k newform with normalised coefficients λ_f(n), the Voronoi formula replaces Poisson. In the form Σ λ_f(n) g(n) = (2πε/√N) Σ λ_f(n) ∫ g(y) J_{k−1}(4π√(ny/N)) dy [S: Iwaniec–Kowalski; U on constants], the self-reciprocal transform is 𝓗φ(z) = 2π ∫₀^∞ φ(w) J_{k−1}(4π√(zw)) dw.

- In r = √w it is the Hankel transform of order ν = k − 1. It is an involution, and r^ν e^{−2πr²} is a +1 eigenfunction [D: Weber's integral]. That eigenfunction gives Γ_ℂ(s + (k−1)/2).
- The E-map h(v) = Σ λ_f(n) ψ(nv) satisfies h(1/(Nv)) ∝ ε·μ·h(v).
- Time-limiting ψ to [0, Λ] and band-limiting 𝓗ψ to [0, Λ′] needs ΛΛ′N = x. In √-variables this is the finite Hankel transform with **c = 4π√(x/N) = 2·T* = d·T***.

So the GL(2) law is the generalised-prolate leakage of order k − 1 at c = dT*, in the sector picked by ε. All constants in this subsection are sketch-level; the numbers below test the consequences, not the constants.

### 6.2 The root-number rule [D + N]

> **Note added 2026-10-03 (main session, P-R389).** The rank-2 curve 389a1 (ε = +1, m₀ = 2) orders its sectors like the ε = −1 curves (odd below even, 5/5 points). The rule below holds for m₀ ≤ 1; in general the ground state moves to the odd sector whenever m₀ > 0, so "(−1)^s = ε" should read "s = 1 iff m₀ > 0". See `docs/DECAY_LAW.md`, P-R389.

The ground state (m = 0, the 𝓗-eigenvalue +1 member with the least leakage) lands in the sector (−1)^s = ε. The other sector starts at m = 1. Hence λ_{−ε}/λ_ε ≈ (1 − λ_{ν,1})/(1 − λ_{ν,0}) > 1, and **sign ln(λ_odd/λ_even) = ε**. This is rh2's rule (80/80 points), now derived at the level of the near-optimal trial.

**Table 6.2** (`data/research_r2/root_number_gap.json`; v = 10, c = 40π):

| object | N | k | ε | measured ln(λ_odd/λ_even) | ε·ln(ℓ^ex_{ν,1}/ℓ^ex_{ν,0}) |
|---|---|---|---|---|---|
| 11a1 / 14a1 / 17a1 | 11 / 14 / 17 | 2 | +1 | 10.46 / 10.55 / 10.71 | 11.71 |
| 37a1 / 43a1 / 53a1 | 37 / 43 / 53 | 2 | −1 | −9.71 / −9.86 / −9.92 | −11.71 |
| g4 / g6 / g8 | 5 / 3 / 2 | 4 / 6 / 8 | +1 | 10.08 / 9.96 / 9.87 | 10.99 / 10.56 / 10.25 |
| Δ / f16 / f20 | 1 | 12 / 16 / 20 | +1 | 9.45 / 9.95 / 9.82 | 9.79 / 9.45 / 9.17 |
| f18 / f22 / f26 | 1 | 18 / 22 / 26 | −1 | −8.51 / −8.78 / −8.83 | −9.30 / −9.04 / −8.82 |

The sign agrees in 15/15. The magnitudes agree to within 0.01–2.0: the prolate ratio overestimates weight 2 by about 1.2–2.0, i.e. R_{m=1}/R_{m=0} ≈ e^{−1.2…−2.0} there.

The central zero for ε = −1 is not a separate mechanism: the E-image vanishes there automatically, like every other zero.

### 6.3 The index and the "level offset" [N]

The P-DL11 continuous index of the exact generalised-prolate leakage on v = 6..10 was compared with the measured index (`hankel_prolate.py`, `exact_leakage.py`):

**Table 6.3** (`data/research_r2/exact_leakage.json`, `hankel_index.json`; shift = measured − prolate):

| object | N | k | ε | sector (m) | n* − k measured | n* − k Hankel prolate | shift |
|---|---|---|---|---|---|---|---|
| 11a1 | 11 | 2 | +1 | even (0) / odd (1) | −0.71 / +1.62 | −0.48 / +1.57 | −0.23 / +0.05 |
| 14a1 | 14 | 2 | +1 | even (0) / odd (1) | −0.54 / +1.50 | −0.48 / +1.57 | −0.06 / −0.07 |
| 17a1 | 17 | 2 | +1 | even (0) / odd (1) | −0.78 / +1.61 | −0.48 / +1.57 | −0.30 / +0.04 |
| 37a1 | 37 | 2 | −1 | even (1) / odd (0) | +1.32 / −0.67 | +1.57 / −0.48 | −0.24 / −0.19 |
| 43a1 | 43 | 2 | −1 | even (1) / odd (0) | +1.69 / −0.52 | +1.57 / −0.48 | +0.12 / −0.04 |
| 53a1 | 53 | 2 | −1 | even (1) / odd (0) | +1.43 / −0.57 | +1.57 / −0.48 | −0.14 / −0.09 |
| g4 | 5 | 4 | +1 | even (0) / odd (1) | −0.59 / +1.88 | −0.44 / +1.64 | −0.15 / +0.24 |
| g6 | 3 | 6 | +1 | even (0) / odd (1) | −0.66 / +1.87 | −0.38 / +1.74 | −0.28 / +0.13 |
| g8 | 2 | 8 | +1 | even (0) / odd (1) | −0.59 / +1.87 | −0.30 / +1.86 | −0.30 / +0.01 |
| Δ | 1 | 12 | +1 | even (0) / odd (1) | +0.49 / +2.77 | −0.06 / +2.18 | +0.55 / +0.59 |
| f16 | 1 | 16 | +1 | even (0) / odd (1) | +0.59 / +3.03 | +0.29 / +2.62 | +0.30 / +0.41 |
| f18 | 1 | 18 | −1 | even (1) / odd (0) | +3.33 / +0.66 | +2.88 / +0.50 | +0.45 / +0.15 |
| f20 | 1 | 20 | +1 | even (0) / odd (1) | +1.06 / +3.57 | +0.74 / +3.17 | +0.31 / +0.40 |
| f22 | 1 | 22 | −1 | even (1) / odd (0) | +3.77 / +1.57 | +3.50 / +1.02 | +0.26 / +0.55 |
| f26 | 1 | 26 | −1 | even (1) / odd (0) | +4.38 / +1.74 | +4.27 / +1.66 | +0.11 / +0.08 |

- **Level N > 1** (weights 2–8): shifts lie in [−0.30, +0.24], mean −0.08.
- **Level 1** (weights 12–26): shifts lie in [+0.08, +0.59], mean +0.35.
- **The "level-1 offset" is mostly a weight effect.** Before the correction it was ≈ 1.3 (P-DL11). The generalised prolate's finite-c correction grows with ν²/c: n* − k = −0.48 at ν = 1, −0.06 at ν = 11 and +1.66 at ν = 25 (m = 0, v = 6..10), against the asymptotic −½. That accounts for about 0.9 of the 1.3.
- **What remains** (~0.4) is unexplained. The candidate "L = log(Nv²) enters the prefactor" would need R ∝ L^{≈2}; untested.
- **The index rule for GL(2) is n* = n*_prolate(ν = k − 1, m, c) + O(0.6)**, with no free per-class offsets. P-DL11/12's level-N offsets (−0.64/+1.81 for ε = +1, fitted on weights 2–8) agree, within 0.3, with this function's values at those weights: −0.48 to −0.30 and +1.57 to +1.86.

### 6.4 The GL(2) data against a bound

**No rigorous bound was evaluated for GL(2).** Two pieces are missing:

- the Voronoi analogue of Theorem 1 is a sketch (§6.1);
- I did not source an explicit zero count N(t) for GL(2) L-functions [U].

The natural first version is the analogue of Corollary 3. Its trial function is the cusp form itself: f(u) = e^{uk/2} f_k(i e^u/√N), whose Mellin transform is Λ(f, ½ + iz). By Theorem 1's mechanism its window truncation has Q = ε′ Σ_ρ F_out(γ_ρ)², with a rate e^{−c} [D, sketch].

What can be stated now [N]: against the exact generalised-prolate leakage, the measured GL(2) minima (scan9, scan11, scan12, v = 6..10) have R_ex = λ/ℓ^ex_{ν,m} in

| class | R_ex |
|---|---|
| level N > 1, m = 0 | [13.7, 82.5] |
| level N > 1, m = 1 | [4.4, 22.1] |
| level 1, m = 0 | [3.2, 15.5] |
| level 1, m = 1 | [2.3, 10.7] |

Overall R_ex ∈ [2.3, 83] (`exact_leakage.json`). This is the GL(2) counterpart of the degree-1 R_ex ∈ [1.9, 39] of Table 5.2.

## 7. Pre-registration-ready experiment (most sharpening): the weight effect out of sample

**Question.** Is the GL(2) index the exact generalised-prolate index of order k − 1 (§6.3), or a constant per-class offset (P-DL11/12)?

**Objects.** The six validated level-1 forms Δ, f16, f18, f20, f22, f26, plus 11a1 as a weight-2 control.

**Grid.** v ∈ {14, 15, 16, 17, 18} (c = 176–226), both sectors, N_basis = 5k′ and 9k′, precision-stable, inv+full. These are P-DL11's protocol and metric. The registration must be committed before any run; these predictions already exist in `data/research_r2/hankel_index_v14_18.json`.

- **Cost.** For level 1, N_basis reaches 9·2vL = 9·2·18·log(324) ≈ 1870. For 11a1, with x = 11v², it reaches about 2650. That is multi-hour per form on ml (P-DL12 ran to 1551), so it must be run by the main session.

**Prediction P-R2a.** For each (form, sector), with Δn* := n*(v = 14..18) − n*(v = 6..10):

- Δn*_measured = Δn*_prolate ± 0.3, where Δn*_prolate is computed from the exact Hankel leakage (`hankel_index_v14_18.json` minus `hankel_index.json`). That is, each form keeps its v = 6..10 shift.

| form | ν | ε | sector (m) | Δn*_prolate |
|---|---|---|---|---|
| Δ | 11 | +1 | even (0) / odd (1) | −0.24 / −0.37 |
| f16 | 15 | +1 | even (0) / odd (1) | −0.43 / −0.61 |
| f18 | 17 | −1 | odd (0) / even (1) | −0.55 / −0.76 |
| f20 | 19 | +1 | even (0) / odd (1) | −0.68 / −0.93 |
| f22 | 21 | −1 | odd (0) / even (1) | −0.84 / −1.13 |
| f26 | 25 | −1 | odd (0) / even (1) | −1.21 / −1.58 |
| 11a1 (control) | 1 | +1 | even (0) / odd (1) | −0.01 / −0.04 |

**Kill.** Either of the following kills the exact-prolate index rule for GL(2):

- more than 2 of the 12 level-1 cases outside ±0.3 of Δn*_prolate;
- the level-1 indices staying put: |Δn*_measured| < 0.3 in at least 5 of the 8 cases where |Δn*_prolate| > 0.5.

In either case the per-class constant-offset description is preferred. 11a1 must stay within ±0.3 of 0, or the instrument is suspect.

**Secondary (D3 lower half).** R_ex = λ/ℓ^ex_{ν,m} must stay within [¼, 4]× its v = 10 value. Kill: outside.

**Why this experiment.** It is out of sample in c and uses no new objects or code. It separates the two GL(2) descriptions by up to 1.4 in n*. It also extends the asymptotic grid that tests D1's lower half (via upper bounds: a value far below R_ex·ℓ^ex would refute D3).

## 8. Reproduction

```bash
.venv/bin/python scripts/research_r2/hankel_prolate.py --check        # generalised prolates = 1D Bouwkamp values
.venv/bin/python scripts/research_r2/check_identity.py --json data/research_r2/check_identity.json   # Theorem 1 at zeros
.venv/bin/python scripts/research_r2/hermite_bound.py --d 1,5,8,-3,-4,-7,-20,-8,-11,12,13,-15,17 --xq 2,3,4,5,6,7,8,10,12,14,16 --json data/research_r2/hermite_bound.json
.venv/bin/python scripts/research_r2/hermite_bound.py --d 5,1 --xq 2,4 --matrix --json data/research_r2/hermite_matrix_check.json
.venv/bin/python scripts/research_r2/kb_bound.py --d 1,5,8 --xq 2,4,6,8,10,16 --dbeta 0,2 --json data/research_r2/kb_bound.json
.venv/bin/python scripts/research_r2/kb_bound.py --d -3,-4,-7,-20 --xq 2,4,6,8,10,16 --dbeta 0,2 --json data/research_r2/kb_bound_b.json
.venv/bin/python scripts/research_r2/prolate_trial.py --d 1,5,8,-3 --xq 2,3,4,6,8 --json data/research_r2/prolate_trial_a.json
.venv/bin/python scripts/research_r2/prolate_trial.py --d -4,-7,-20 --xq 2,3,4,6,8 --json data/research_r2/prolate_trial_b.json
.venv/bin/python scripts/research_r2/prolate_trial.py --d -8,-11,12,13,-15,17 --xq 2,4,6 --json data/research_r2/prolate_trial_c.json
.venv/bin/python scripts/research_r2/hankel_prolate.py --nu 1,3,5,7,11,15,17,19,21,25 --json data/research_r2/hankel_index.json
.venv/bin/python scripts/research_r2/hankel_prolate.py --nu 1,11,15,17,19,21,25 --v 14,15,16,17,18 --json data/research_r2/hankel_index_v14_18.json
.venv/bin/python scripts/research_r2/exact_leakage.py --json data/research_r2/exact_leakage.json
.venv/bin/python scripts/research_r2/root_number_gap.py --json data/research_r2/root_number_gap.json
.venv/bin/python scripts/research_r2/summarize.py
```

## 9. What I trust least

1. **Theorem 4's numbers are not certified.**
   - ‖f_win‖² is a 25-digit Gauss–Legendre value. The envelope integrals use mpmath quadrature.
   - The tails of the n-sums rely on env(y)·y^{m+1} being nonincreasing beyond y = 40λ. That is asserted on a grid, not proved.
   - The zero-count constants (Trudgian 2014; BMOR 2021) were taken from abstracts [A]. BMOR note an error in Trudgian's 2015 Dirichlet constant, which is why theirs is used.
   - The *inequality* of Theorem 2 is proved. Each tabulated *value* is [N]. A ball-arithmetic pass would make Table 3.5 a theorem.
2. **The GL(2) Voronoi/Hankel construction is a sketch** (§6.1): constants, the ε convention and the √-variable normalisation were not checked against a source. The evidence for it is that the exact generalised-prolate leakage reproduces the measured c-slopes (the continuous index) and the sign and size of the sector gap. That index metric tests only c-dependence, not R. Two further points are weak:
   - the residual level-1 shift (+0.35 mean) is unexplained;
   - the Hankel computation was validated against 1D prolates (ν = ±½) and by a basis-size check, but not against an independent generalised-prolate source.
3. **"Near-optimal" is measured against finite-basis minima.** Trial and measured values are both upper bounds in the *same* N-mode basis. `docs/DECAY_LAW.md` reports |ln(λ_{5k*}/λ_{9k*})| ≤ 0.17, so the true minima may sit somewhat lower. A trial/measured ratio near 1 then means "near the basis optimum", not "near the true infimum". The largest ratios (2.1–2.5) come from projection quadrature and are not physics. And ZUP (§4.4), the lower half of the law, is untouched by all of it.
