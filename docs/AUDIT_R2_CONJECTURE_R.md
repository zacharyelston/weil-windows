# Audit of GPT's R2 review and Conjecture R attempt (branch gpt/r2-conjecture-r)

**Auditor:** Claude Opus 5.5, for Zac Elston. **Audited:** commits ac89d4d, fcab696 and ed4fd8b on `gpt/r2-conjecture-r`, based on 2c936d6. **Date:** 2026-10-03.

**Verdict.** The review's three corrections are right, and the M2 proof and the conditional reduction hold. The finite gates point in the right directions, but their registered targets are so loose that they could not have failed. The 92/92 passes carry no evidence and should not be quoted as support. M1 is open as GPT reports. The certified norm data independently favour it and match GPT's cancellation sketch [N]. Recommendation: merge after the four changes listed at the end.

## 1. Reproduction and process

- All three scripts were rerun with the pinned Python. Every output is byte-identical to the committed files except one field, which records the commit current at run time.
  - That field is named `registration_commit`, so a rerun overwrites it with a later commit. Rename it to `run_commit` and store the registration commit as a constant.
- The registration was committed before the gates ran: fcab696 at 17:58, then the gate file in ed4fd8b at 18:11. The gate file records fcab696 as the commit current when the gates ran. The registered constants did not change afterwards.
- The manifest passes 265/265. Every pre-existing claim group is still run, and `cert_lib.py` is unchanged.

## 2. The review (`docs/REVIEW_R2_THEOREMS.md`)

| Finding | Audit | Action |
|---|---|---|
| Theorem 1(iii): the ordinary Mellin integral is claimed on the strip of half-width ½ + δ | **Confirmed.** The sheet states the strip with the integral (§3, item iii). The integrand behaves like y^{−½−Im z} near 0, so the integral diverges for Im z ≥ ½ whenever ψ(0) ≠ 0. For characters ψ(0) ≠ 0 in general. The narrower strip −½ − δ < Im z < ½ still contains every nontrivial zero, so no conclusion changes. | Edit the sheet's §3 |
| Theorem 2 and Lemma Z: the BMOR count used below its validity threshold for large q | **Confirmed for the general statement.** For q ≥ 12 the zero-free height in the code falls to 5/7, where BMOR gives no zero-free interval. The integration by parts needs only a count bound valid above the cutoff A/B. The code asserts that the cutoff lies in range, and GPT found 0 failures in the stored cases. | Edit the statement; no number changes |
| §8(a) is not a check of the theorem's inequality | **Confirmed.** §8(a) projects f_win onto our N-mode basis. Projection does not preserve the Rayleigh quotient, and GPT's two-by-two example shows the effect. A captured fraction near 1 does not control the form error, because the form is unbounded. §8(c) remains a valid comparison, since the projected vector is itself a trial. Both stay [N]. | Relabel §8(a) as an [N] diagnostic |
| ζ zero-count source: keep Bellotti–Wong v2, with Platt's S(T) bound taken as [S] | Agreed. Use HSW with C₃ = 9.4925 only if the paper adopts a refereed-only policy. | Author's decision |

The remaining verdicts (Lemmas G, K, E, H, T hold; Propositions P and H hold with stated additions) were read and not re-derived.

## 3. M2 and the conditional reduction (`docs/R2_CONJECTURE_R.md`)

I re-derived each step.

- **|a| ≤ 10c.**
  - The Gaussian bound on I₂ follows from 1 − cos θ = 2 sin²(θ/2) ≥ 2θ²/π² on [0, π].
  - The half-integer formulas have the stated signs. The cubic factor in I_{7/2} is at least 0.535 at b = 10 and increases with b.
  - The bound obtained is (5/2)(√(2π) + 1)c ≈ 8.8c.
  - The stored ζ coefficients are about 2.02c in magnitude in the even sector and below 0.35 in the odd sector, so the bound holds with room.
- **Envelope integrals.**
  - √(ω² − β²) ≥ ω√(4(c − 1))/c ≥ ω/√c on [c, ∞) holds for c ≥ 4/3.
  - The weighted Hankel integral has the closed form stated. The coefficient sums are at most 2431 for k = 5.
  - The exponent condition p − (k − 1)/2 ≤ 2 holds for every term actually used; the largest case is the ω⁴ env₅ derivative term.
  - The cap piece is bounded under an explicit assumption that the split point is at most 100c. That assumption is stated, and it is satisfied, since the true crossover is below s = 100.
- **Edge sum.** ψ(nλ) is out of band, with ω = nc. The terms decay like n⁻³ with a factor c^{−½}, and √x·pref/2 ≤ √q c³/2. The bound 5·10⁵ √q c^{7/2} follows.
- **The reduction.** With A = B = H the cutoff is 1. The integration by parts is exact. The count majorant 5(1 + log 3q)·t(1 + log t) follows from BMOR for t ≥ 1, with room in the constant 5. The zero sum is at most 4K_qH², which gives the stated C and p = 24.

**Comment.** The proof is correct and deliberately coarse. The cap split alone costs a factor of about c³ in A and B. Below the true crossover the cap is the smaller envelope, so bounding the cap piece by the Hankel integral gives p₁ = 5 and p = 10 + p₂. That is an improvement to make, not an error.

## 4. The gates: correct directions, no power to fail

Each gate compares the right endpoints in the right direction. But the registered targets sit far from the data:

| Gate | Smallest margin over 92 rows (log₁₀) |
|---|---|
| G1: C₁c⁸ against the upper ends of A and B | 20.4 |
| G2: the norm target against the new upper enclosure | 16.9 |
| G3: C c²⁴ e^{−2c} against `bound_upper` | 64.8 |

A gate that only fails when the data move by 17 orders of magnitude cannot fail on anything this code produces. The brief asked for gates that can fail. These can in principle, but not in practice. **The 92/92 should be reported as "passed, with margins of at least 10¹⁷; not evidence".**

- The slack comes from the registration: C₂ = 10¹² e^{4πq} is e^{251} times 10¹² at q = 20, and C ≈ 2·10⁴³ q(1 + log 3q)e^{4+4πq}.
- For a future registration, set exploratory targets within one or two orders of magnitude of the data, so that a pass means something. Constants that come from a proof are whatever the proof gives; in that case state the margin next to the pass.

The new upper norm enclosure is sound and tight. Its relative excess over the certified lower bound has median 3·10⁻¹² and maximum 1.5·10⁻⁴. The two enclosures overlap at every row, and none contradicts the other. That pins ‖f_win‖² to better than 10⁻⁴ relative error, which makes the next section possible.

## 5. Independent diagnostic: what the certified norms say about M1 [N]

Fit ln‖f_win‖² − 2β = α − p₂ ln c per object and sector, over the six grid points (c = 12.6 to 100.5, candidate β = c − 2). Add a linear term κc to test for exponential loss.

| Object | Sector | p₂ (fit) | Largest residual (ln) | κ |
|---|---|---|---|---|
| χ₅, χ₈ | s = 0 | 0.93–0.94 | ≤ 0.03 | −0.002 |
| χ₋₃, χ₋₄, χ₋₇, χ₋₂₀ | s = 0 | 1.88–1.91 | ≤ 0.05 | −0.003 |
| χ₅, χ₈ | s = 1 | 3.70–3.74 | ≤ 0.18 | +0.013 |
| χ₋₃, χ₋₄, χ₋₇, χ₋₂₀ | s = 1 | 4.35–4.41 | ≤ 0.11 | +0.008 |
| ζ | s = 0 | 3.37 | 0.11 | +0.008 |
| ζ | s = 1 | 8.37 | 0.44 | +0.033 |

- **No exponential loss is visible.** Every κ is small, and where it is not negligible it is positive, so the decay slows as c grows. The quantity M1 bounds behaves like a power of c on this range. An exponential loss is what would make M1 false.
- **The sector pattern matches GPT's cancellation sketch.**
  - The χ exponents sit near 1 + κ_parity + 2k, where k is the order of the leading cancellation: k = 0 for s = 0 and k = 1 for s = 1. Here κ_parity = 1 for odd χ, because the trial carries an extra factor t. That predicts 1, 2, 3 and 4. For s = 1 the local slopes at c = 100 are 2.9 (even) and 3.9 (odd), closer to 3 and 4 than the whole-range fits.
  - For ζ with s = 1, GPT's sketch says the whole first-order term cancels once a → 1/8. That predicts k = 2 and an asymptotic p₂ near 5. The fitted local slope falls from about 9.3 at c = 12.6 to about 6.4 at c = 100, which is consistent with that prediction.
  - The stored ζ coefficient in that sector is 0.349, 0.217, 0.183, 0.167, 0.158 and 0.145 on the grid, and 0.141 at support 2.99. That is close to 1/8 + 2/c, so it confirms the formal limit.
  - ζ with s = 0 is not explained by this reading. Its local slope is near 2.9 at c = 100.
- This is a reading of six points per class, not a law. It suggests the target for M1: an asymptotic p₂ per class, derived from the expansion and then proved with an explicit remainder.

## 6. GPT's three least-sure points

1. **The remainder after cancellation in the ζ minus sector.** This is the real open step. The data agree that the first-order term cancels.
2. **An explicit norm lower bound for all parameters.** Open; the data show no obstruction.
3. **How conservative the cap-split exponent is.** It is valid but costs about c³; see §3.

## Changes requested before merge

1. Rename the gate file's `registration_commit` field to `run_commit`, and record fcab696 as the registration.
2. In `docs/R2_CONJECTURE_R.md`, state the margins from §4 next to the 92/92, and say that the passes are not evidence.
3. Apply the three sheet corrections in `docs/R2_THEOREMS.md`: the Theorem 1 strip, the count range in Theorem 2 and Lemma Z, and §8(a) as [N]. Then update the matching lines in `docs/PRIMER.md` and `docs/PAPER_OUTLINE.md`.
4. Optional, and the natural next step for M1:
   - Pre-register an asymptotic p₂ for each class, derived from the expansion.
   - Then compute the norm at x/q = 32 and 64. The norm alone is cheap; no certificate is needed.
   - Only then attempt the explicit remainder.
