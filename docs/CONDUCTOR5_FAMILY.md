# Conductor-5 Family: Real Off-Line Zeros and the Odd Sector

This follows the recommended next shot in `docs/LITERATURE_PASS.md` (Q5). The instrument is `scripts/conductor5_family_mp.py`; the odd sector lives in `scripts/connes_letter_mp.py` (`parity="odd"`).

## The family

F_t(s) = (1 − t)L(s, χ₅) + t(1 + √5·5^{−s})ζ(s), with χ₅ = (·/5). For every t,

Λ_t(s) = (5/π)^{s/2} Γ(s/2) F_t(s) = Λ_t(1 − s),

checked to 1e-30 at t = 0, t*, 0.5 and 1.

| t | Euler product | Pole at s = 1 |
|---|---|---|
| 0 | yes, L(s, χ₅) | no |
| 0 < t < 1 | no (a₂a₃ − a₆ = −0.21 at t*) | yes |
| 1 | yes, (1 + √5·5^{−s})ζ(s) | yes |

The special value t* = L(¾, χ₅)/(L(¾, χ₅) − G(¾)) = 0.0552921974… gives F_{t*}(¾) = F_{t*}(¼) = 0. That is a real off-line pair at γ = ∓i/4.

## Why odd test functions are needed

For a real test function φ, a real pair at γ = ±iδ contributes:
- **+2φ̂(iδ)² ≥ 0** if φ is even;
- **−2φ̂(iδ)² ≤ 0** if φ is odd.

So only the odd sector can detect a real off-line pair, and every earlier our instrument was even-only. The odd sector uses the basis √(2/L)·sin(2πku/L) on [−L/2, L/2]. Its zeros-side form is E − 2wwᵀ, where w_k = ∫ b_k(u) e^{u/2} du; the pole correction has the opposite sign from the even sector.

## Checks

- **Odd-sector explicit formula for ζ** (x = 13, test function sin θ + ½ sin 2θ): E = 2Σ_γ |φ̂(γ)|² + 2ĝ(i/2)² to 12 digits, using 300 zeros.
- **Census of F_{t*}**, γ ∈ [0.5, 40]: 16 zeros by the argument principle, all 16 on the critical line, and no complex off-line zeros. On (0, 1) the only real zeros are ¼ and ¾.
- **Explicit formula for F_{t*}, both sectors, against the census zeros** (real pair included):

  | x | Sector | Zeros-side form | Sum over census zeros | Real pair's contribution |
  |---|---|---|---|---|
  | 5 | even | 12.70969697 | 12.70969697 | **+11.69** |
  | 5 | odd | 0.7320131 | 0.7320130 | **−0.012** |
  | 7 | even | 17.50934026 | 17.50934026 | **+17.12** |
  | 7 | odd | 0.38168008 | 0.38168006 | **−0.026** |

  The remaining difference comes from stopping the zero list at height 40.

## Results: smallest eigenvalue of the zeros-side form (N = 32, 40 digits)

| t | Sector | x = 2 | 3 | 5 | 7 | 10 | 15 |
|---|---|---|---|---|---|---|---|
| 0 | even | 0.32 | 0.044 | 6.9e-4 | 5.0e-6 | 3.8e-9 | 2.9e-14 |
| 0 | odd | 1.70 | 0.76 | 0.060 | 1.4e-3 | 3.1e-6 | 3.9e-11 |
| **t*** | even | 1.61 | 1.37 | 0.57 | 0.082 | 6.3e-4 | 2.0e-8 |
| **t*** | **odd** | 1.69 | 0.77 | 0.029 | **−0.045** | **−0.081** | **−0.169** |
| 0.5 | even | 1.61 | 1.54 | 1.02 | 0.15 | 9.4e-4 | 3.5e-8 |
| 0.5 | odd | 1.69 | 1.19 | 0.83 | 0.75 | 0.070 | 7.9e-6 |
| 1 | even | 1.61 | 1.61 | 1.61 | 0.027 | 2.5e-4 | 1.3e-8 |
| 1 | odd | 1.69 | 1.61 | 1.61 | 0.14 | 0.017 | 2.1e-6 |

**Witness** (`data/conductor5/witness_tstar_x7_odd.json`): the N = 32 odd minimiser at x = 7 is an explicit odd trigonometric polynomial on [−½ log 7, ½ log 7]. Its Rayleigh quotient is −0.045279565835127186296, identical at 40, 60 and 80 digits. It is computed from the archimedean closed forms and the prime-side terms for n ≤ 7 alone, with no zeros involved.

## Reading

1. **The literature pass's prediction holds.**
   - The odd sector of F_{t*}'s form is negative at x = 7, 10 and 15 (N = 32), with one negative direction.
   - The even sector is positive in the tested finite basis (N = 32) at every tested x. These values are upper bounds on its full-space minimum, not proofs of positivity. The real pair enters even test functions with a + sign, so it cannot by itself make the even sector negative. Zeros above height 40, which have not been censused, could still do so.
   - A negative finite Rayleigh quotient proves a negative direction of the full form (Connes–Consani, Cor. 2.4). So F_{t*}'s first full-space failure is at x ≤ 7.
2. **The real pair is the likely cause.** The census covers the real segment (0, 1) and γ ∈ [0.5, 40]. In that range the only off-line zeros are the real pair. The explicit-formula checks above show that the pair contributes −0.026 to the odd form at x = 7, and the census zeros reproduce the form to 7 digits. That is consistent with the pair causing the failure, provided the zeros above height 40 are on the line. Without a census to greater height, or a bound on the remainder of the zero sum, it is not a proof of attribution.
3. **The witness is a numerical certificate, not an interval-arithmetic one.** It is an explicit test function at an explicit support, x = 7 (support length log 7), with a negative form value stable to 20 digits across 40–80-digit precision. The archimedean series tails are summed to below 10^(−dps). Making it rigorous would need interval evaluation of the digamma closed forms and the prime terms.

   The test function is a trigonometric polynomial cut off at the window edges. Its sine terms vanish there, so it is continuous, but it is not smooth; the form extends to it by density, using trigonometric polynomials as a core. A smoothed version would be needed to fit the strict C_c^∞ setting.
4. **The other family members stay positive in the tested finite basis:** t = 0 and t = 1 (Euler products), and t = 0.5, whose zeros have not been censused. These positive values are only upper bounds.

## Caveats

- **The crossover has an upper end only.** Only N = 32 has been run. Its finite-basis transition lies in (5, 7], and the negative witness at x = 7 gives numerical evidence of full-space failure by x = 7. The positive N = 32 value at x = 5 does not show that the full form is positive there, so the full-space first failure is at x ≤ 7 with no lower bound. Larger N (48, 64) can only move the finite transition down.
- The witness is not interval-certified.

## Rigorous bracket (`docs/CONTROL_CERTIFICATES.md`)

F_{t*}'s first full-space failure lies in **[5.31217, 5.3625)**, a bracket 1% wide. Both ends are certified in Arb:
- **Positive end:** Q ≥ 1.1666·10⁻³‖f‖² for every f supported in [−0.835, 0.835], in both sectors.
- **Negative end:** a witness with Q/‖f‖² = −9.889·10⁻⁷ at x = 5.3625.

Only the odd sector fails, as the real off-line pair predicts: at a = 0.9 the even sector is still certified above 0.2786. This supersedes the N = 32 statement "(5, 7]" above.
