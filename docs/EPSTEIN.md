# Discriminant −20: Euler Product On and Off, Same Functional Equation

A degree-2 control pair with an identical completed functional equation:

| Function | Definition | Euler product |
|---|---|---|
| ζ_K | ζ(s)L(s, χ₋₂₀), the Dedekind zeta of ℚ(√−5) | yes |
| Z₁ | ½[ζ_K(s) + L(s, χ₋₄)L(s, χ₅)], the Epstein zeta of x² + 5y² (class number 2) | no |

Both satisfy Λ(s) = 20^{s/2}(2π)^{−s}Γ(s)F(s) = Λ(1 − s), and both have a simple pole at s = 1.

Scripts:
- `scripts/epstein_mp.py`, steps 1–2: functions, checks, zero census.
- `scripts/epstein_connes_mp.py`, step 3: Connes' form in degree 2.

## Pre-registration

- **Prediction:** Z₁ has off-line zeros, and its Connes form turns negative at some support while ζ_K's stays positive.
- **Kill:** no off-line zeros for Z₁ up to height 300.

## Step 1: the functions

- Λ(s) = Λ(1 − s) holds to 1e-30 for ζ_K, L(χ₋₄)L(χ₅), Z₁ and Z₂ (the form 2x² + 2xy + 3y²).
- The Hardy functions are real on the critical line.
- Z₁'s coefficients equal r(n)/2, half the number of representations by x² + 5y², for every n ≤ 200.

## Step 2: zero census, γ ∈ [1, 60]

| Function | Argument principle on [−2, 3] | On the line | Off-line pairs |
|---|---|---|---|
| ζ_K | 54 | 54 | 0 |
| Z₁ | 53 | 39 | 7 |

No zeros were found with γ < 1, and none on the real segment (0, 1), where all of these functions are negative.

Z₁'s off-line zeros (β > ½; the mirror points 1 − β + iγ are implied):

| β | γ |
|---|---|
| 0.9329697 | 15.66824953 |
| 0.9376669 | 29.98339524 |
| 0.6969271 | 36.37406369 |
| 0.8231873 | 44.00011318 |
| 0.6349076 | 46.75840866 |
| 0.5974373 | 52.74326641 |
| 0.8111721 | 58.50773733 |

These come from an argument-principle census with Newton refinement, not interval-certified. Literature on Epstein off-line zeros (Potter–Titchmarsh, Bateman–Grosswald, Stark) is listed in `docs/LITERATURE_PASS.md`. These specific values have not been checked against a published table.

## Step 3: Connes' form in degree 2

**Construction.**
- Archimedean part: Γ_ℝ(s)Γ_ℝ(s+1) at conductor 20 gives Q_arch = Q_ζ,arch + Q_DH,arch + log 4 · I, since the DH archimedean part already contains log 5 · I.
- Prime side: c_n from −F′/F, for n ≤ x.
- Zeros-side form: QW = Q + 2vvᵀ, adding the pole correction as for ζ.

**Explicit-formula check.** The test function is φ = (1 + cos 2πu/L)², which vanishes to fourth order at the window edges. The form is compared with the sum over Z₁'s census zeros, with off-line pairs at complex arguments, minus the pole correction:
- x = 5: relative difference 8e-11.
- x = 13: relative difference 3.5e-13.

Both are at the scale set by stopping the zero list at height 60.

**Smallest eigenvalue ε₀ of QW_λ** (N = 48; 50 digits for x ≤ 10, 30 digits beyond):

| x | 2 | 3 | 5 | 7 | 10 | 15 | 20 | 30 | 40 | 60 |
|---|---|---|---|---|---|---|---|---|---|---|
| ζ_K | 2.08 | 1.70 | 1.08 | 0.80 | 0.37 | 0.051 | 4.6e-3 | 2.2e-5 | 1.5e-7 | 4.1e-10 |
| Z₁ | 2.08 | 1.91 | 1.85 | 1.20 | 0.52 | 0.074 | **−1.3e-4** | −0.23 | −0.73 (n₋ = 2) | −1.66 (n₋ = 2) |

## Reading

1. **The prediction holds.**
   - Z₁'s form has a negative finite Rayleigh quotient at x = 20. Finite trigonometric minima decrease to the full-space infimum (Connes–Consani, Cor. 2.4), so this is a genuine negative direction of the full form. The evidence is numerical: the explicit formula is validated, but the arithmetic is not interval-certified.
   - Z₁'s first full-space failure is therefore at x ≤ 20. Its positive finite values at x ≤ 15 do not certify a lower bound.
   - A second negative direction appears by x = 40.
2. **ζ_K stays positive in the tested finite basis up to x = 60.** These positive values are upper bounds on its full-space minimum, not proofs of positivity.
3. **The pattern is the same as DH's.** Z₁ has off-line zeros, the lowest at 0.933 + 15.67i with δ = 0.43, and its form has a negative direction once the support is large enough. The Euler-product partner with the same functional equation shows none in the tested bases. That is consistent with the off-line zeros causing the failure, but no single zero has been isolated as the cause. Z₁ is a natural arithmetic function, not a constructed one like DH, and it shows the effect at small x.
4. **Degree 2 sits far from saturation.**
   - At small support the forms have O(1) minima. ζ_K's ε₀ falls only like e^{−0.5x} between x = 20 and 40, against roughly e^{−4πx} for ζ.
   - This is consistent with the Γ(s) kernel, an exponential rather than a Gaussian, concentrating less well in a finite window. The degree-2 analogue of the prolate floor has not been worked out.

## Crossover check with larger bases

| x | N = 48 | N = 64 | N = 80 |
|---|---|---|---|
| 15 | +0.07450 | +0.07435 | +0.07419 |
| 17.5 | | +0.01556 | +0.01550 |
| 20 | −1.33e-4 | −1.58e-4 | −1.70e-4 |

ζ_K at the same points (N = 64; N = 80): +0.0505, +0.0168, +0.00457; +0.0504, +0.0168, +0.00456.

Z₁'s finite-basis crossover lies in (17.5, 20] and moves only slightly with N. The negative value at x = 20 grows in magnitude as N increases, as it should variationally. The full-space first failure is at x ≤ 20.

## Caveats and next steps

- **Convergence.** N = 64 and 80 confirm the finite-basis transition in (17.5, 20]. The full-space first failure is at x ≤ 20, with no lower bound.
- **Odd test functions.** The even cosine basis detects complex off-line zeros, as Z₁'s are. Real off-line zeros need odd test functions, because an even φ gives a pair term +2φ̂(iδ)² ≥ 0 (`docs/LITERATURE_PASS.md`, recommended next shot).
- **Z₂** (2x² + 2xy + 3y²) has a₁ = 0, so it needs a different normalisation before it can be run through the same instrument.

## Rigorous bracket (`docs/CONTROL_CERTIFICATES.md`)

Z₁'s first full-space failure lies in **[15.79984, 19.844)**.
- **Positive end:** both sectors are certified, with even c = 0.049235, at T♯ = 1000.
- **Negative end:** a certified witness at x = 19.844, with Q/‖f‖² = −1.534·10⁻⁷. The x = 20 witness (−1.3261·10⁻⁴) is certified too.

## Rigorous same-support separation: ζ_K positive where Z₁ fails

Z₁ and its Euler partner ζ_K = ζ·L(s, χ₋₂₀) share the completed functional equation (conductor 20, Γ_ℝ(s)Γ_ℝ(s+1)) and the pole at s = 1. They differ only in the Euler product. At support 2a = 2.99, that is x = e^{2.99} ≈ 19.886:

| Function | Euler product | Certified statement on [−1.495, 1.495] |
|---|---|---|
| ζ_K | yes | **Q(f) ≥ 4.74884875180·10⁻³ ‖f‖²** (even) and ≥ 0.598720158932 ‖f‖² (odd), for every f. |
| Z₁ | no | **Q(f) < 0 for some C_c^∞ f.** The certified witness lies at x = 19.8439 < 19.886 (`docs/CONTROL_CERTIFICATES.md`). |

**ζ_K certificate.** Zhu's reduction with T♯ = 1000 and 1162 modes, in Arb, using `scripts/control_certificate.py --function zetak`.
- **Comb:** exactly Λ(n)(1 + χ₋₂₀(n)) on n ∈ {2, 3, 4, 5, 7, 8, 9, 16}. The primes 11, 13, 17 and 19 have χ₋₂₀ = −1 and drop out.
- **Normalisation:** agrees with our ζ_K matrices to 1e-39 in both sectors.
- **Runs:** the even sector ran on the M3 (930 s) and the odd sector on the M4 (516 s). The exact Cholesky residual is ≤ 5.2·10⁻⁷¹ in each.
- **Consistency:** Our finite-basis upper bounds at the same x are 4.866e-3, 4.845e-3 and 4.834e-3 at N = 48, 64 and 80. So **4.749·10⁻³ ≤ λ*_even(ζ_K) ≤ 4.834·10⁻³**.

This is the pre-registered prediction of this document (Z₁'s form turns negative while ζ_K's stays positive), now as a pair of certificates at the same support. It is the cleanest Euler-versus-non-Euler separation in the repository: identical archimedean data, opposite outcomes. The separation is exhibited, not explained; why the Euler product protects ζ_K is the open question.
