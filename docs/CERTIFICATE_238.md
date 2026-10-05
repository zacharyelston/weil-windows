# Weil Positivity on Support 2.38: a Certificate

**Claim.** For every complex f ∈ L²(ℝ) supported in [−1.19, 1.19], with Q Weil's quadratic form for ζ in Zhu's normalisation (`docs/AUDIT_ZHU.md`, §1), the following holds. On all of L², Q is understood as the extended form given by its full-line frequency integral, with values in (−∞, +∞]; not every compactly supported L² function has finite Q.

  Q(f) ≥ c ‖f‖²,  c = 6.81311646960 × 10⁻⁴⁸

(from the corrected code's rerun; the earlier published constant, 6.81311646951 × 10⁻⁴⁸, remains valid).

That is support 2a = 2.38, x = e^{2.38} ≈ 10.80, with the prime powers 2, 3, 4, 5, 7, 8 and 9 in the comb. Zhu attempted this support, then withdrew the claim (arXiv:2608.24827 v2, §7): its envelope used a lower bound on the prime comb where an upper bound was needed. The largest support certified in the literature we found is Liu's 2.125 (17/16 on each side, 15 Sep 2026), with constant 2^-49162.

**Status:** computer-assisted, completed 2026-10-01. Every number that feeds the bound is an Arb ball. The trust assumptions are listed at the end.

**Audit: verified with corrections** (`docs/AUDIT_CERTIFICATE_238.md`, 64758f0, independent auditor). Two implementation defects were found and are now fixed at the source:
- **Schur maximum.** It was selected through binary64, understating the true maximum upper endpoint by 1.4e-16. It is now chosen by exact comparison of Arb endpoints. The audit's independent bound, μ ≤ 3.46615908 on a larger window, already validated the value used.
- **Integer residual.** It omitted the rounding of A_mid − μ before flooring. A_mid, μ′ ≤ μ and L̃ are now converted to integers exactly.

The audit also:
- supplied a sharper quadrature debit;
- found that the odd decimal had been rounded up;
- corrected the L² and complex-f wording;
- settled the prior art (below).

## Method

The method is Zhu's finite reduction with the prime comb's torus supremum A replaced by μ̄: the top of the comb's shift operator on the slightly widened window [−b, b], b = a + δ (`docs/GRID_NORM.md`).

1. **Split each f by frequency.** Let κ = (sin Ωu/πu)·W(u), where W is the Kaiser window I₀(β_K√(1 − u²/δ²))/I₀(β_K) on [−δ, δ]. Put m = 1 − κ̂ and h = f − f∗κ, so ĥ = mF and supp h ⊂ [−b, b].
2. **Bound the high-frequency prime terms.** ⟨h, P h⟩ ≤ μ̄‖h‖² gives Q(f) ≥ pole + (1/π)∫₀^∞ Φ|F|², with Φ = Ψ + m²(P − μ̄).
3. **Envelope.** For t ≥ T♯ = Ω + β_K/δ, Φ ≥ log(t/2π) − 1/t − μ̄ − (2ε + ε²)(A + μ̄) ≥ β′. This uses Re ψ(¼ + it/2) ≥ log(t/2) − 1/t for t ≥ ¾ (`docs/AUDIT_ZHU.md`) and the sidelobe bound ε = 2/(π β_K I₀(β_K)).
4. **Reduce to a finite block.** Zhu's reduction then gives Q ≥ min(λ_min(A_N), β′ − ε_D) − ε_B, where A_N is the Legendre block of R′(f) = pole + (1/π)∫₀^{T♯}(Φ − β′)|F|² + β′‖f‖².
5. **Complex f.** For complex f, Q(f) = Q(Re f) + Q(Im f), and each part splits into even and odd sectors with no cross term. So the constant is the minimum over the two sectors.

## Parameters and certified components

| Quantity | Value | How it is obtained |
|---|---|---|
| a, δ, b | 1.19, 0.3, 1.49 | — |
| comb | n ∈ {2, 3, 4, 5, 7, 8, 9}, A = 7.07500562853 | Arb |
| μ̄ on [−1.490000001, 1.490000001] | ≤ 3.46920322110 | Schur test, K = 4000 cells, Collatz–Wielandt weight, cell indices certified (`scripts/grid_norm.py`) |
| β_K, Ω, T♯ | 64, 486.67, 700 | — |
| sidelobe bound ε | 3.19e-29 | Second mean value theorem |
| β′ | 1.24257147610376 | exact dyadic below the Arb bound |
| m(t) in the band | Arb; cross-checked to 22 digits against direct u-space integration at 5 points | acb.integral (main lobe), integration-by-parts tail expansion with remainder (K−1)!/r₀^K |
| nodes | 156,800 (2,800 panels × 56-point Gauss) | exact Gauss–Legendre balls |
| quadrature c_err | ≤ 1.77e-58 | Bernstein ellipse, ρ = 2 + √5, \|Φ − β′\| ≤ 316.05 on the strip |
| modes | 680 per sector (even: degrees 0–1358; odd: 1–1359) | — |
| coupling ε_B, ε_D | ≤ 1.2e-103, ≤ 7.4e-217 | as in the audit, with G_real = 26.91 |

## Results

| Sector | λ₁(A_mid) | Cholesky shift μ_s | exact residual | max entry radius | certified c |
|---|---|---|---|---|---|
| even | 6.81313469666e-48 | 6.81312788353e-48 | 2.96e-118 | 1.68e-56 | **6.81311646960e-48** (corrected-code rerun: 6.8131164696007…; first run 6.81311646951e-48) |
| odd | 3.91967253723e-44 | 3.91966861756e-44 | 2.20e-118 | 2.15e-56 | **3.91966861609e-44** (corrected-code rerun: 3.9196686160954…; previously misprinted …610, rounded up) |

The constant for complex f is the minimum over the two sectors: **c = 6.81311646960e-48**. The audit-corrected code produced these numbers: exact Schur maximum, exact integer residual, and the Cauchy–Schwarz quadrature debit (c_err ≤ 8.86e-59). They no longer depend on the audit's independent bounds. The certified values sit just below the floating λ₁ of the block, and the debits are:
- the residual, plus N × (max entry radius) ≈ 1.1e-53 (even) and 1.5e-53 (odd);
- c_err ≤ 1.8e-58;
- ε_B ≤ 1.2e-103.

The tail branch β′ − ε_D = 1.2426 is not active.

**Run history.** The first even-sector run gave the same Cholesky success at 6.8131e-48, but a max entry radius of 3.5e+97. The cause was the Bessel recurrence amplifying the 400-bit radius of its input x = L·t by about e^{0.47x}. The recurrence now runs at the exact dyadic midpoint, and the node radius is added afterwards using |j_n′| ≤ ½ (DLMF 10.54.2). The Schur window was also rounded up from the float 1.19 + 0.3, which lies below 1.49, to 1.490000001. μ̄ was unchanged to 12 digits.

**Our upper bounds at a = 1.19** (`scripts/zhu_window_mp.py`, N = 100, 160 digits): even 1.139e-47, odd 6.37e-44. Each certified lower bound must lie below the corresponding upper bound.

## Trust assumptions and limits

- **Arb's ball arithmetic:** python-flint 0.9.0, including acb.integral, digamma, 0F1, the Bessel functions and Gauss–Legendre roots. It is assumed to enclose correctly.
- **The new scripts are correct.** That covers `scripts/grid_certificate_rigorous.py` and the Schur test in `scripts/grid_norm.py`. Their checks so far:
  - the prototype reproduces the audited a = 0.8 block (λ₁ = 1.02768956005196e-17);
  - the odd-sector conventions reproduce Zhu's 9.11833845e-15;
  - the Bessel enclosures agree with mpmath to 70 digits;
  - the Kaiser filter agrees with direct u-space integration to 22 digits.
- **Analytic steps, derived in this repository and not externally reviewed:**
  - the frequency-split inequality;
  - the Bonnet sidelobe bound;
  - the remainder of the tail expansion;
  - the derivative bound |j_n′| ≤ ½;
  - the Re ψ bounds used for G_real.
- **The explicit formula itself** is taken as standard, as in the audit.
- **Prior art (from the audit's reading of the papers).**
  - **Liu** (Theorem B, Appendix B.4) already bounds the prime terms by compressed prime translations with a weighted Schur estimate, and certifies widths 2 and 2.125. **The window-compressed prime bound is therefore not new.**
  - **Suzuki** (arXiv:2606.09096, §2.3–2.4) writes the compressed translations explicitly, but certifies only small windows. arXiv:2607.24830 is a numerical realisation of Suzuki's operator by Kim, Hong, Kim, Choi, Jang and Kim.
  - **What this certificate adds:** support 2.38, which none of these works certifies, reached by a different construction. That construction inserts a Kaiser frequency split into Zhu's reduction, so the compressed norm replaces the torus supremum only at high frequency.

## Reproduce

```bash
.venv/bin/python scripts/grid_certificate_rigorous.py --a 1.19 --delta 0.3 --beta-k 64 --tsharp 700 --sector even --nmodes 680 --gl 56 --workers 4 --json data/connes/grid_certificate_a119_even.json
.venv/bin/python scripts/grid_certificate_rigorous.py --a 1.19 --delta 0.3 --beta-k 64 --tsharp 700 --sector odd  --nmodes 680 --gl 56 --workers 4 --json data/connes/grid_certificate_a119_odd.json
```

Each sector takes about 11 minutes on an M3 with 4 workers.

## Extension: support 2.6 (a = 1.3, x = e^{2.6} ≈ 13.46)

The same method was run in the pinned Docker environment (`docker/`) on the 32-thread host. **The code is the audited code with the audit's fixes (e864488), plus two later changes the audit did not cover:**
- 2ce5975: `--prec`, file-based Gram collection, a parallel m-table, and a blocked Cholesky;
- 686e1f1: Bessel ratio iteration at base precision.

See "Post-audit code changes" below.

**Parameters:**
- δ = 0.3 and b = 1.6;
- comb n ∈ {2, 3, 4, 5, 7, 8, 9, 11, 13};
- Schur bound μ̄ ≤ 4.4106016523106 (K = 4000);
- Kaiser β_K = 80, with sidelobe bound ε = 3.2·10⁻³⁶;
- T♯ = 1240, Ω = 973.33, β′ = 0.87358;
- 1250 modes, 64-point Gauss (317,440 nodes);
- c_err ≤ 9.2·10⁻⁶⁸, ε_B ≤ 4.2·10⁻¹³⁸.

| Sector | Certified lower bound | Our upper bound (cosine basis, N = 100) |
|---|---|---|
| even | **5.76344479222·10⁻⁶²** | 1.04·10⁻⁶¹ |
| odd | 5.36011304048·10⁻⁵⁸ | 9.74·10⁻⁵⁸ |

So **Q(f) ≥ 5.76344479222·10⁻⁶² ‖f‖² for every complex f supported in [−1.3, 1.3]**.
- **Checks:** exact residuals ≤ 7.7·10⁻¹¹⁸ and max entry radius ≤ 2.7·10⁻⁷⁰.
- **Run time:** 16 workers per sector; Gram 1100 s, Cholesky 168 s.
- **Data:** `data/connes/grid_certificate_a130_{even,odd}.json`.

The audit (`docs/AUDIT_CERTIFICATE_238.md`) covered the 2.38 run. This run uses the same audited code at different parameters, and has not been independently re-audited.


## Extension: support 2.99 (a = 1.495, x = e^{2.99} ≈ 19.89)

Same method, same post-audit code (686e1f1) and same Docker environment, at 600-bit working precision with 28 workers.

**Parameters** (`data/connes/grid_certificate_a1495_even.json`):
- δ = 0.2 and b = 1.695;
- comb n ∈ {2, 3, 4, 5, 7, 8, 9, 11, 13, 16, 17, 19};
- Schur bound μ̄ ≤ 5.27794743005 (K = 8000);
- Kaiser β_K = 120, with sidelobe bound 1.1·10⁻⁵³;
- T♯ = 2870, Ω = 2270, β′ = 0.845894;
- 3200 modes (max degree 6398), 96-point Gauss (1,102,080 nodes);
- c_err ≤ 1.84·10⁻¹⁰⁶, ε_B ≤ 1.29·10⁻²⁵², ε_D ≤ 4.7·10⁻⁵¹⁷.

| Sector | Certified lower bound | Floating λ₁ | Our upper bound |
|---|---|---|---|
| even | **3.50114217641·10⁻⁹⁶** | 3.52079517811·10⁻⁹⁶ | 6.27·10⁻⁹⁶ (cosine N = 180) |
| odd | 8.25626494001·10⁻⁹² | 8.25627515668·10⁻⁹² | 1.48·10⁻⁹¹ (N = 180) |

So **Q(f) ≥ 3.50114217641·10⁻⁹⁶ ‖f‖² for every complex f supported in [−1.495, 1.495]**: the even sector is the minimum, the odd sector gives 8.25626494001·10⁻⁹².
- **Checks:** exact residual ≤ 3.5·10⁻¹³⁸ at μ_s = 3.5207917·10⁻⁹⁶, and max entry radius 6.1·10⁻¹⁰².
- **Run time:** Gram 18047 s (5.0 h, 28 workers), Cholesky 189 s.
- **Data:** run log `data/connes/grid_certificate_a1495_even.run.txt`.
- **Prediction check:** the value lies inside the bracket [2.9, 5.9]·10⁻⁹⁶ that `docs/LIT_TANGENTS.md` (Test 2) predicted from the prolate h₄ leakage before the run finished. So 3.501·10⁻⁹⁶ ≤ λ*(1.495) ≤ 6.27·10⁻⁹⁶.
- **Odd sector:** exact residual ≤ 4.9·10⁻¹⁴⁰ at μ_s = 8.2562669·10⁻⁹², max entry radius 6.1·10⁻¹⁰²; Gram 17943 s, Cholesky + residual 234 s, total 6h27m. Also inside its predicted bracket [5.7, 13]·10⁻⁹² (`docs/LIT_TANGENTS.md`, Test 2). Data: `data/connes/grid_certificate_a1495_odd.json`, log `data/connes/grid_certificate_a1495_odd.run.txt`.

## Post-audit code changes

The 2.38 audit covered 6a1bccd, and its fixes were applied in e864488. The 2.6 and 2.99 runs use 686e1f1, which adds 144 lines to `scripts/grid_certificate_rigorous.py` (`git diff e864488 686e1f1`):
- `sph_j_enclosure`: the ratio iteration now runs at base precision, and only the upward recurrence uses extra bits.
- **Gram collection:** worker partials are written to marshal'd files of (mid mantissa, mid exponent, radius mantissa, radius exponent) and rebuilt with `_from_pack`.
- `_mtab_segment`: the m-table is computed in parallel segments.
- `cholesky_blocked`: a blocked floating Cholesky.

Soundness does not depend on how the floating factor is computed, because the exact integer residual bounds any factor. The serialisation, m-table and Bessel changes do affect soundness if they are wrong. **They have not been independently audited.** A regression check reruns the audited a = 1.19 certificate with both versions on the same machine.

**Regression result** (M4, `data/connes/regress/`). The audited parameters were rerun with the post-fix audited code e864488 and the current code 686e1f1:

| run | even certified bound | odd certified bound |
|---|---|---|
| audited JSON (M3, e864488, 4 workers) | [6.8131164696007119625e-48 ± 2.94e-68] | [3.9196686160954221947e-44 ± 1.36e-64] |
| M4, e864488, 4 workers | [6.8131164696007119625e-48 ± 2.94e-68] | — |
| M4, 686e1f1, 4 workers | [6.8131164696007119625e-48 ± 2.94e-68] | — |
| M4, e864488, 2 workers | [6.8131164695990460601e-48 ± 3.47e-68] | [3.9196686160954218615e-44 ± 3.31e-64] |
| M4, 686e1f1, 2 workers | [6.8131164695990460601e-48 ± 3.47e-68] | [3.9196686160954218615e-44 ± 3.31e-64] |

- **The post-audit changes do not alter the certificate.** At equal worker count the two code versions agree bit for bit, in both sectors.
- **The result is reproducible across machines and Python versions:** M3 with Python 3.14.3, and M4 with 3.14.4.
- **It depends on the worker count, deterministically.** The node range is split into one contiguous block per worker. The partial Gram matrices are summed in ball arithmetic, so a different partition gives slightly different midpoints and radii, and so a slightly different μ_s. Each run is a rigorous certificate; they differ by a relative 2.4e-13.
- **To reproduce a quoted bound digit for digit, use its recorded `--workers`.**
- **The check's reach is limited.** It exercises the post-audit code at the 2.38 parameters (400 bits, 680 modes). It does not replace an audit of the 144-line diff, which is listed for the next certificate audit.

## Code changes after the 2.99 run (performance; no certified value uses them yet)

The cross-host benchmark (`docs/PERF_BENCHMARK.md`) found that the Gram phase is 87–89% FLINT matrix product and that the two Python-dominated costs are the static partition (a 74-minute straggler at 2.99) and the per-entry collection (36 minutes). Three changes to `scripts/grid_certificate_rigorous.py` (and `smallest_two` in `scripts/grid_certificate_flint.py`):

1. **`--partition strided`** (new default). Chunks stay contiguous in t, so each chunk's cost is unchanged, and are dealt round-robin to the workers. `--partition contiguous` is the original behaviour. The partition is recorded in every run's `args`.
2. **Byte-record collection.** Workers write fixed-width two's-complement records (mantissa width derived from the precision; `int.to_bytes` raises on overflow rather than truncating), and the main process does one exact rebuild and one matrix add per worker file instead of a per-entry Python loop.
3. **`--lambda2`** opt-in. λ₂ was a diagnostic for the gap; it is skipped by default.

**Soundness.** The serialisation is exact and the accumulation order per worker file is unchanged, so with the same partition the result must be bit-identical. The strided partition changes the summation order, as a different worker count already did.

**Regression** (M4, audited a = 1.19 parameters, 4 workers; `data/connes/regress/new_*`):

| run | certified bound |
|---|---|
| old code, contiguous (reference) | [6.8131164696007119625e-48 ± 2.94e-68] |
| new code, contiguous, λ₂ skipped | [6.8131164696007119625e-48 ± 2.94e-68] **bit-identical** |
| new code, strided, `--lambda2` | [6.8131164696001308337e-48 ± 2.47e-68] (relative 8.5e-14) |

λ₁ agrees to all printed digits in every run; the strided Gram phase took 296 s against 312 s. Projected at support 2.99 on the P620: 6.9 h → about 3.8–4.25 h per sector. These changes are listed for the certificate audit's round 2.

## Planning support 3.2 (a = 1.6, x = e^{3.2} ≈ 24.5): not run

`--plan` runs on ml (`data/connes/plans/plan_a160*.run.txt`), with the Schur bound μ̄ ≤ 5.90239716658 at b = 1.8 (K = 8000) and the comb n ∈ {2, …, 23}:

| variant | T♯ | Ω | β′ | nodes | modes | gl | bits | c_err | ε_B |
|---|---|---|---|---|---|---|---|---|---|
| naive scaling | 4300 | 3650 | 0.626 | 1,651,200 | 5200 | 96 | 720 | 1.0e-105 | 4.9e-474 |
| A | 2500 | 1850 | 0.0834 | 1,200,000 | 3000 | 120 | 720 | **1.4e-136** | 4.4e-251 |
| B | 2500 | 1850 | 0.0834 | 1,120,000 | 3000 | 112 | 700 | 1.5e-126 | 4.4e-251 |

The decay law (`docs/DECAY_LAW.md`) predicts the even margin near 6e-121 at this support, so the naive plan's quadrature budget does not fit, while A and B do.

**β′ was over-provisioned.** The envelope condition is β′ = log(T♯/2π) − 1/T♯ − μ̄ − w(A + μ̄) > 0, and the certificate takes min(λ_A, β′ − ε_D) − ε_B with λ_A ≈ 1e-96 or smaller, so any β′ > λ_A serves. The 2.6 and 2.99 runs carried β′ ≈ 0.85 from Zhu's setup; β′ ≈ 0.08 needs T♯ ≈ 2π e^{μ̄ + 0.08}, which at a = 1.6 is 2500 instead of 4300. Nodes and modes scale with T♯ and the matrix product with modes², so the saving is about 3× here and would have been about 2× at 2.99 (T♯ 2870 → ~1340). Nothing certified changes; the earlier runs were merely more expensive than necessary.

**Decision.** A support-3.2 run is about the size of the 2.99 run (about 5 h per sector on ml with the strided code; the odd sector can run concurrently on ml2). The exclusion-region analysis (`docs/RESEARCH_R1_EXCLUSION.md`) answered that question as far as it can be answered: exclusion is a witness computation in any case, and under the hypothesis that all other zeros lie on the line the forecast reach of supports up to 4.0 is far below the verified height (a forecast, not a theorem). **Support 3.2 is not run.** Further certificates would be justified only as tests of the decay law.

