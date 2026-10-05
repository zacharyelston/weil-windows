# The Prime Grid Through a Window: a Sharper Envelope for Zhu's Reduction

Step 2 of the plan after `docs/AUDIT_ZHU.md`. Scripts: `scripts/grid_norm.py` (the window norm) and `scripts/grid_certificate_flint.py` (the reduction, an Arb/FLINT prototype).

**Status: numerical prototype, not a certificate.** Every number here is computed at high precision but not enclosed.

## The question

The audit found that Zhu's reduction needs Ψ = H − P ≥ β on [T♯, ∞), so T♯ ≳ 2πe^A.
- **H:** the archimedean term.
- **P:** the prime comb, Σ (2c_n/√n) cos(t log n).
- **A = Σ 2|c_n|/√n:** the comb's supremum.

The supremum is approached because the line t·(log 2, log 3, log 5, …) winds around the prime torus and returns arbitrarily close to the origin (Kronecker). That is the spiral on the torus. Its returns force the doubly exponential cost of Zhu's Theorem 1.4.

A test function supported in a window cannot use those returns fully. In u-space, P acts as Σ (c_n/√n)(τ_{log n} + τ_{−log n}), and a shift by log n pushes part of the function out of the window. What such a function can see is

  μ(b) = sup { ⟨f, P f⟩ / ‖f‖² : supp f ⊂ [−b, b] },

the top of the spectrum of the prime grid's weighted shift graph, restricted to the window. Pre-registered questions:
1. **Is μ well below A?**
2. **Does μ/A separate Euler products from the controls?** Kill: the same ratio for ζ and DH.

## Result 1: μ(b) ≈ A/3

Computed by piecewise-constant Galerkin (`scripts/grid_norm.py`): a lower bound that converges as the cells shrink, with K = 1000–4000 cells and agreement to about 3 digits.

| b | 0.8 | 0.9 | 1.0 | 1.19 | 1.25 | 1.5 | 2.0 | 2.5 |
|---|---|---|---|---|---|---|---|---|
| A (ζ) | 2.942 | 4.381 | 5.852 | 7.075 | 8.521 | 13.02 | 24.38 | 42.55 |
| μ (ζ) | 1.219 | 1.646 | 1.943 | 2.668 | 2.959 | 4.292 | 7.974 | 13.64 |
| μ/A | 0.414 | 0.376 | 0.332 | 0.377 | 0.347 | 0.330 | 0.327 | 0.321 |
| 2πe^A | 119 | 502 | 2,187 | 7,427 | 31,536 | 2.8e6 | 2.4e11 | 1.9e19 |
| 2πe^μ | 21 | 33 | 44 | 90 | 121 | 459 | 18,241 | 5.3e6 |

**Question 1: yes.** For ζ the ratio settles near 1/3 out to b = 2.5. The threshold therefore drops from e^A to roughly e^{A/3}. That is still doubly exponential in the support, but with an exponent about three times smaller.

**Question 2: no; the kill criterion is met.** The controls give the same ratio:

| b | F_{t*} | DH | Z₁ | ζ |
|---|---|---|---|---|
| 0.8 | 0.402 | 0.369 | 0.500 | 0.414 |
| 1.0 | 0.340 | 0.323 | 0.381 | 0.332 |
| 1.19 | 0.333 | 0.318 | 0.338 | 0.377 |

The window effect is generic geometry: shifts leave the window. It is not something the Euler product supplies. Like every earlier instrument, it does not separate ζ from the controls.

### The geometry predicts μ/A → 1/4

A function confined to [−b, b] does best by putting its mass at both ends of the window. The longest prime shifts, with log n near 2b, then link the two ends. By the prime number theorem the comb's weight Σ Λ(n)/√n has density about e^{ℓ/2} in ℓ = log n. An end-concentrated trial function therefore gives μ ≈ e^b = √x, while the torus supremum is A ≈ 2∫₀^{2b} e^{ℓ/2} dℓ ≈ 4√x.

| b | 3.0 | 3.5 | 4.0 |
|---|---|---|---|
| μ (Galerkin, K = 1000) | 22.50 | 36.54 | 59.08 |
| e^b | 20.09 | 33.12 | 54.60 |
| μ/e^b | 1.12 | 1.10 | 1.08 |
| μ/A | 0.299 | 0.287 | 0.277 |

So μ/A → 1/4: the window lets a test function use about a quarter of the torus supremum. Zhu's threshold 2πe^{A} ≈ 2πe^{4√x} becomes about 2πe^{√x}. The exponent drops by a factor of 4, but the cost is still doubly exponential in the support (x = e^{2a}).

## Result 2: using μ in the reduction

**Derivation.** Let K be the Fourier transform of a kernel κ supported in [−δ, δ], and set m = 1 − K and w = 1 − m². For f supported in [−a, a], the function h = f − f∗κ has ĥ = mF and is supported in [−a−δ, a+δ]. Then

  Q(f) = pole + (1/π)∫ Ψ(w + m²)|F|²
       = pole + (1/π)∫ Ψ w|F|² + (1/π)∫ H m²|F|² − ⟨h, P h⟩
       ≥ pole + (1/π)∫ Φ |F|²,  where  Φ = Ψ w + (H − μ(a+δ)) m².

For t ≥ T♯, where m ≈ 1, this gives Φ ≥ H − μ − |w|(A + μ). So Zhu's reduction applies with Φ in place of Ψ, and the envelope only needs 2πe^{μ(a+δ)}, not 2πe^A.

**Version 1 failed.** A polynomial bump κ = (1 − u²/δ²)⁸ at a = 0.8 gave λ₁ = −0.077. Its m switches on near t ≈ 10, where H − μ ≈ −1.5. The leakage there swamps the form's 1e-17 margin.

**The design constraint.** m must stay below √(margin) until H exceeds μ, at t ≈ 2πe^μ, and then rise to 1. That calls for a sharp low-pass κ: a windowed sinc with a Kaiser window of parameter β_K, so the ripple is about 1/I₀(β_K) and the transition width is 2β_K/δ.

**Version 2: windowed-sinc kernel, Kaiser window** (`--mode kaiser`). The code reproduces the audited a = 0.8 block (`--mode zhu`) exactly, λ₁ = 1.02768956005196e-17. Results:

| a | δ | μ(a+δ) used | β_K | T♯ | modes | λ₁(R′) | rh2 upper bound (N = 100) | Zhu's T♯ would need |
|---|---|---|---|---|---|---|---|---|
| 1.0 | 0.25 | 2.98 (Galerkin 2.959) | 40 | 465 | 280 | **4.574e-30** | 6.498e-30 | ≈ 2,300 |
| 1.19 | 0.30 | 4.26 (Galerkin 4.230) | 62 | 900 | 600 | **6.805e-48** (λ₂ = 1.234e-40) | 1.139e-47 (ε₁ = 2.1e-40) | ≈ 7,800 |

Both blocks are positive and lie below rh2's upper bounds, as lower bounds must: the ratios are 0.70 at a = 1.0 and 0.60 at a = 1.19. The cutoffs are 5× smaller (a = 1.0) and 9× smaller (a = 1.19) than Zhu's reduction needs. At a = 1.19, Zhu's retracted support of 2.38, the run took 7m48s on the M3.

**Open check: mode convergence.** At 600 modes the top degree is 1198, only about 130 above a·T♯ = 1071. A run with 760 modes and leading-block λ₁ values is in progress.

**The μ values are rigorous upper bounds.** A Schur test on a piecewise-constant weight, with the weight tuned by Collatz–Wielandt iteration of the Schur majorant and the cell indices certified in Arb (`grid_norm.py --schur`), gives:

| b | Galerkin (lower) | Schur, rigorous upper (K = 1000 / 2000 / 4000) | value used |
|---|---|---|---|
| 1.25 | 2.96024 | 2.99222 / 2.97861 / **2.97271** | 2.98 |
| 1.49 | 4.23080 | 4.26904 / 4.25200 / **4.24517** | 4.26 |

Both values used in the runs above exceed the certified upper bound, so item 1 below is closed for ζ.

**Restricting μ to the comb.** The comb only contains n < e^{2a}, but μ(b) above was computed with every prime power below e^{2b}. For nonnegative shift weights (the operator is positivity-preserving, though not positive semidefinite), dropping shifts can only lower the norm, and the Schur bound bounds the norm, so the restricted operator gives a smaller rigorous bound. At a = 1.19, b = 1.49, using the shifts 2, 3, 4, 5, 7, 8, 9: Galerkin 3.4603, Schur upper bound **μ̄ ≤ 3.46920**. That lowers the cutoff from 900 to 700.

**Correction to the a = 1.19 prototype row.** Its eigenvalue came from inverse iteration, which finds the eigenvalue smallest in absolute value and could miss a negative one. In a deliberately under-resolved test (8-point quadrature), a full eigensolve did find a spurious negative eigenvalue that inverse iteration had missed. It disappeared with 32-point quadrature. The rigorous run decides the question by Cholesky: success at shift μ_s proves A − μ_s I positive definite.

## The rigorous certificate (`scripts/grid_certificate_rigorous.py`)

Each numerical shortcut in the prototype is replaced by an enclosure:

| Ingredient | Rigorous treatment |
|---|---|
| μ̄ | Schur test with a Collatz–Wielandt-tuned piecewise-constant weight; cell indices certified in Arb |
| Kaiser m(t), passband | \|m\| ≤ 2ε, ε = 2/(π β_K I₀(β_K)), from the second mean value theorem on each sidelobe (after r = √(δ²s² − β²): ∫ sin r/√(r²+β²)) |
| Kaiser m(t), band | m = (1 − main)/2 + (1/2π)∫_{−β/δ}^{t−Ω} Ŵ + right tail. The main lobe uses acb.integral of the entire 0F1 form. The right tail uses an integration-by-parts expansion with remainder (K−1)!/r₀^K, from \|g^{(K)}\| ≤ K!/r^{K+1} with g = ∫J₀(βs)e^{−rs}ds. |
| Kaiser beyond T♯ | \|K\| ≤ ε, so \|w\| ≤ 2ε + ε², entering β′ |
| Ψ | acb.digamma, plus the comb in Arb |
| j_n | Upward recurrence below the turning point at precision 400 + 1.25x bits; ratio enclosure r_n ∈ (0, 1) above it |
| Quadrature | Exact Gauss–Legendre balls; Bernstein-ellipse truncation bound with a uniform \|Φ − β′\| bound on the strip (\|m\| ≤ 1 + e^{δb}(2/π)(1 + log Ωδ)) |
| Coupling | ε_B and ε_D as in the audit, with G_real ≥ sup\|Φ − β′\| on [0, T♯] |
| Eigenvalue | Floating Cholesky at μ_s; exact integer residual (fmpz_mat); λ_min ≥ μ_s − ‖R‖_∞ − N·(max entry radius) − c_err |

Plan for a = 1.19 (`--plan`):
- **Parameters:** δ = 0.3, β_K = 64, T♯ = 700, 680 modes, 56-point Gauss.
- **Budgets:** c_err ≤ 1.8e-58, ε_B ≤ 1.2e-103, ε_D ≤ 7.4e-217, all far below the expected margin of about 5e-48.

## What a certificate still needs

1. ~~A rigorous upper bound on μ(b).~~ Done for ζ by Schur's test (above). For signed coefficients (the controls), apply the test to |P|.
2. **Enclosures for the kernel.** G(τ) = ∫Ŵ is integrated with 4-point Gauss here. That needs `acb.integral`, plus a rigorous ripple bound in the passband and beyond T♯, which this prototype sets to exactly 0 and 1.
3. **Bessel enclosures.** The Miller recurrence runs in plain MPFR, because ball arithmetic overestimates its error growth. Use the audit's ratio enclosure instead.
4. **The remaining certificate steps:** the coupling bounds (ε_B, ε_D) with Φ in place of Ψ, an ellipse bound for Φ (m is entire), and a ball Cholesky.

## Prior art (from a first search; not yet read in full)

- **Yoshida (1992)** and **Bombieri (2000):** positivity for small support.
- **Suzuki** (arXiv:2606.09096): Weil's form on (−a, a) as an operator. Its numerical realisation, arXiv:2607.24830, is by Kim et al. Its decomposition reportedly includes "finitely many prime-power partial translations", which are the compressed shifts used here.
- **V. Liu, "Certified Weil Positivity Beyond the Unit Window"** (15 Sep 2026, alphaxiv): certifies a = 1 (bound 2^-151) and a = 17/16 (bound 2^-49162). It uses "weighted Schur estimates" for the prime terms and Legendre orders up to 2815. This is the closest prior work. The method here may overlap with it, and that must be checked before any novelty is claimed. What this prototype adds is sharpness: 4.6e-30 at a = 1 with 280 modes, and the attempt at a = 1.19.
  - **Settled by the audit** (`docs/AUDIT_CERTIFICATE_238.md`): Liu's Theorem B and Appendix B.4 explicitly use compressed prime translations with a weighted Schur estimate. **The window-compressed prime bound is prior art.** The Kaiser frequency split inside Zhu's reduction, and the support 2.38 it reaches, are what this repository adds.
