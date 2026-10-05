# Audit: Zhu's certified Weil-positivity bound on [−0.8, 0.8]

Date: 2026-10-01. Worktree `<review-worktree>`, branch `audit/zhu-certificate`,
base `9bd1805` (`feat/weil-gram-instrument`). Nothing in the main checkout was touched.

**Claim under audit.** Xuefeng Zhu, *Weil positivity in compact windows: a finite reduction, certified
two-sided bounds, and a Landau–Widom decay law*, arXiv:2608.24827 **v2** (2 September 2026):
Theorem 1.1 (reduction), Theorem 1.2 (real even f), Corollary 6.3 (complex f):
Q(f) ≥ 8.9·10⁻¹⁸ ‖f‖₂² for supp f ⊂ [−0.8, 0.8], prime powers 2, 3, 4 included; §7 and Remark 3.3
(retraction of a support-2.38 claim).

**How the paper was read.** The arXiv abstract page and the v2 HTML rendering were read with WebFetch,
which returns a model-mediated extraction rather than the raw page. The PDF and TeX source were not
downloaded. Every item marked *taken from the paper* rests on that extraction. The decisive inputs
(A_L, β*, T♯, basis, quadrature grid, the odd-sector eigenvalue and both Cholesky shifts) were
reproduced independently below, which corroborates the extraction where it matters.

**Rigour labels.** **[R]** re-derived rigorously (interval arithmetic with `mpmath.iv`, exact integer
arithmetic, or a written analytic argument); **[N]** reproduced numerically (high precision, no
enclosure); **[P]** taken from the paper; **[U]** could not verify. [R] means rigorous *given* that
`mpmath.iv` rounds outward correctly and that the new scripts below are correct. They have not been
reviewed by a second person.

## Verdict

**Verified with corrections.** Zhu's reduction was re-executed at his parameters (T♯ = 200,
200 even Legendre modes, 800 panels × 32-point Gauss), with every integrand value enclosed. It gives
λ_min(A) ≥ 1.02768955902·10⁻¹⁷ for the even block. Cholesky also succeeds at Zhu's shift 9·10⁻¹⁸. The
odd block reproduces Zhu's λ_min = 9.1183·10⁻¹⁵ to all five printed digits, and certifies Q ≥ 9.118·10⁻¹⁵
on odd f. So **Q(f) ≥ 1.0276·10⁻¹⁷ ‖f‖₂² ≥ 8.9·10⁻¹⁸ ‖f‖₂² for every complex f supported in [−0.8, 0.8]**
[R]. rh2's upper bounds, now extended to N = 800, stay above it: 1.6356·10⁻¹⁷ at N = 800 [N].

The corrections concern the paper's rigour claims and one internal inconsistency. None changes the
constant:

1. The paper bounds the quadrature truncation, but evaluates the integrand values (j_n by Miller
   recurrence, ψ, nodes and weights) in 50-digit floating point without enclosure. That is a formal gap,
   more than 25 orders of magnitude below the margin, and this audit closes it (§2).
2. The constant |Ψ_L − β*| ≤ 20 on Zhu's ρ = 6.55 ellipses is stated, not derived. Numerically the
   maximum is 14.9 [N]. The audit uses a smaller ellipse with a proved bound of 19.41 [R].
3. Theorem 6.2 quotes the upper bound λ₁^even ≤ 2.523·10⁻¹⁶. §11 and the abstract give the sharper
   2.27·10⁻¹⁷. Both are upper bounds, so they are not contradictory, but Theorem 6.2's "clearance factor
   ≥ 32" uses the weaker one.
4. The explicit-formula identity is stated for continuous f, while Theorems 1.2 and 6.3 are stated for
   L² f. The bound extends to L² if Q is *defined* on L² by the frequency integral (value in (−∞, +∞]).
   The proof of Theorem 1.1 is pointwise in frequency, so it then applies verbatim. This is a wording
   point, not a gap.

The decisive reason for "verified" is the independent interval re-execution together with the exact
odd-sector match. "With corrections" refers only to items 1–4.

## 1. Normalisation map

| | Zhu (v2, §§1–2) [P] | rh2 (`docs/CONNES_LETTER.md`, `scripts/connes_letter_mp.py`) |
|---|---|---|
| Variable | additive u; f on [−L, L], **L = a = 0.8** | additive u = log v of Connes' multiplicative variable; window [−L/2, L/2], **L = log x** (rh2's L is the full width) |
| Window ↔ x | a = 0.8 | x = e^{2a} = e^{1.6} = 4.953032424395115 |
| Prime powers | log n < 2a: n = 2, 3, 4 | n ≤ x: 2, 3, 4 (e^{1.6} < 5 [R]) |
| Fourier | F(r) = ∫ f e^{iru} du | φ̂(z) = ∫ φ e^{−izu} du (t ↦ −t; irrelevant, since \|F(t)\|² is even for real f) |
| Norm | ‖f‖₂² = ∫\|f\|² du | Euclidean norm of coefficients in an L²(du)-orthonormal basis: the same norm |
| Form | Q(f) = ĝ(i/2) + ĝ(−i/2) + (1/2π)∫\|F\|²Ψ_L, g = f⋆f̃, Ψ_L = Re ψ(¼ + it/2) − log π − Σ_{n ≤ 4} (2Λ(n)/√n) cos(t log n) | QW (zeros side, "nopole") = geometric form E **plus** the pole correction |
| Pole term | 2F(i/2)F(−i/2) = **+2F(i/2)²** (real even), **−2(∫f sinh(u/2))²** (real odd) | even: E **+ 2vvᵀ**, v_k = b̂_k(i/2); odd: E **− 2wwᵀ**, w_k = ∫b_k e^{u/2} |
| Which side | zeros side (Σ_ρ \|F(γ)\|² under RH) | zeros side |
| Complex f | Q(f) = Q(Re f) + Q(Im f), and each splits into even + odd | two sectors computed separately |

**The map is the identity.** rh2's QW at x = e^{2a} is Zhu's Q on [−a, a], with the same norm. The
minimum over complex f is the minimum over the two real sectors.

*Why complex f reduces to the two real sectors [R, algebraic].* The kernel of Q is real and symmetric,
so Q(a + ib) = Q(a) + Q(b) for real a and b. For real f = e + o (even plus odd), F_e is real and F_o is
imaginary, so \|F\|² = F_e² + \|F_o\|² has no cross term. The pole term is 2(E² − O²) with
E = ∫e·cosh(u/2) and O = ∫o·sinh(u/2), so it has none either. Hence
inf_complex Q/‖f‖² = min(inf_even, inf_odd).

**Numerical confirmation** (`scripts/zhu_normalisation_check.py`, `data/connes/zhu_normalisation.json`,
40 digits, 300 zeta zeros) [N]. The test functions are f = (1 − u²/a²)⁴ (even) and
g = (u/a)(1 − u²/a²)⁴ (odd). Q was computed three ways: on the zeros side (Σ over zeros, which are on
the line), with Zhu's frequency formula, and with rh2's matrix applied to the test function's cosine or
sine coefficients:

| | Zeros side, 300 zeros (γ ≤ 541.85, tail ≤ 1.2·10⁻¹⁸) | Zhu's frequency side | rh2 matrix, N = 60 | rh2 matrix, N = 120 |
|---|---|---|---|---|
| Q(f), even | 3.58899889612999·10⁻⁶ | 3.58899889613007·10⁻⁶ | 3.58899890315083·10⁻⁶ | 3.58899889618617·10⁻⁶ |
| Q(g), odd | 2.50814981521881·10⁻⁵ | 2.50814981521882·10⁻⁵ | 2.50814935005184·10⁻⁵ | 2.50814980038527·10⁻⁵ |
| ‖f‖², ‖g‖² (exact) | 0.479261392202568673, 0.0252242838001351933 | Plancherel: same to 20 digits | Σc_k²: same to 20 digits | same |
| Q(T₀), T₀ = 1/√(2a) | independent u-space quadrature: 0.0802740549593 | 0.08027403 (neglected tail ≤ 3.5·10⁻⁷) | (E + 2vvᵀ)₀₀ = 0.0802740549593 | — |

Zhu's formula and the zeros sum agree to 2·10⁻¹⁴ relative in the even sector, which is within the
zeros-tail bound, and to 2·10⁻¹⁸ in the odd sector. rh2's matrix agrees to within its truncation error:
it converges from N = 60 to N = 120 towards the same value, and Σc_k² equals ‖f‖² to 20 digits. For the constant function T₀ = 1/√(2a), which is both rh2's b₀ and Zhu's first
Legendre mode, an independent u-space evaluation gives Q(T₀) = 0.0802740549593. That equals rh2's
(E + 2vvᵀ)₀₀ to all printed digits, and also matches Zhu's frequency formula (row T₀ in the table). A
first run of that frequency-side check gave 0.07669 because `mpmath.quadosc` mis-summed the 1/t² tail;
the corrected script integrates the tail directly.

**Conversion of rh2's table: multiply by 1.** The table in the brief is already in Zhu's normalisation.

## 2. Re-execution of the certificate (Theorem 1.1 at a = 0.8)

Script: `scripts/zhu_certificate_iv.py`. Outputs: `data/connes/zhu_certificate_even.json` and
`data/connes/zhu_certificate_odd.json`. Run parameters: 110-digit intervals, 200-bit fixed point,
80-digit floating Cholesky, 8 processes, 76 s (even) and 56 s (odd).

**The reduction.**
- Plancherel gives ‖f‖² = (1/π)∫₀^∞\|F\|².
- If Ψ_L ≥ β on [T♯, ∞), then Q(f) ≥ R(f) = pole + (1/π)∫₀^{T♯}(Ψ_L − β)\|F\|² + β‖f‖².
- On the Legendre basis T_n(u) = P̄_n(u/L)/√L, the transform is T̂_n(t) = iⁿ S_n(t), with
  S_n(t) = √(2L(2n+1)) j_n(Lt).
- In sign-adjusted coordinates, R's block is A = s·2ppᵀ + Σ_q c_q S(t_q)S(t_q)ᵀ + βI. Here
  c_q = w_q(Ψ_L(t_q) − β)/π, p_n = ±√(2L(2n+1)) i_n(L/2), and s = +1 (even) or −1 (odd).
- Splitting f into the block part and its complement, B the coupling block and D' = D − βI the
  tail-block deviation, gives Q ≥ (min(λ_min(A), β − ε_D) − ε_B)‖f‖², where ε_B ≥ ‖B‖ and ε_D ≥ ‖D'‖
  are the bounds computed below.

### Even sector (Theorem 1.2)

| Component | Zhu v2 [P] | This audit | Level |
|---|---|---|---|
| A_L = Σ 2Λ(n)/√n, n = 2, 3, 4 | 2.9419735… | 2.941973525223620455503909 | [R] |
| T₁ = 2π e^{A_L} | — | 119.0865561876757 | [R] |
| T♯ | 200 | 200 | — |
| β* = log(T♯/2π) − 1/T♯ − A_L | 0.5134667… | 0.5134667749150707383886 (used: a dyadic β̃ just below) | [R] |
| Envelope (Lemma 3.1): Re ψ(¼ + it/2) − log π ≥ log(t/2π) − 1/t | for t ≥ 15/4 (Binet) | for t ≥ 3/4 (Binet; proof below) | [R] |
| sup P_L = A_L (Lemma 3.2) | stated | at t = 0, trivially | [R] |
| Ψ_L − β* ≥ 0 on [200, ∞) | follows | follows; also min 0.109 on a 0.01-grid of [200, 4000] | [R] / [N] |
| Basis | 200 even Legendre modes, degree 0…398 | same | — |
| Quadrature | 800 panels of width ¼, 32-point Gauss | same grid | — |
| Ellipse | ρ = 6.55, semi-minor ≈ 0.400 | ρ = 2 + √5, semi-minor ¼ | — |
| \|Ψ_L − β\| on the ellipses | ≤ 20 (stated) | ≤ 19.41 (proof below); actual max on Zhu's ellipses 14.90 | [R] / [N] |
| M (integrand bound) | ≤ 4.9·10⁴ | ≤ 1.175·10⁴ | [R] |
| Per-entry quadrature error ε_Q | ≤ 1.03·10⁻⁴⁴ | ≤ 9.2·10⁻³⁴ (includes the node and weight defects, Σ\|def_k\| ≤ 4.3·10⁻⁸⁹) | [R] |
| c_err = N ε_Q | ≤ 2.06·10⁻⁴² | ≤ 1.84·10⁻³¹ | [R] |
| Integrand values j_n(Lt_q), Ψ_L(t_q) | Miller recurrence, 50-digit floating point, not enclosed | enclosed (ratio continued fraction plus recurrence; shifted Stirling with Binet remainder) | [R] |
| Gram sum from enclosures | — | exact integer sum of midpoints; entry radius ≤ 5.4·10⁻⁵⁷, spectral ≤ 1.1·10⁻⁵⁴ | [R] |
| λ_min(Ã) (floating) | not stated for even | **1.027689560051963·10⁻¹⁷** (inverse iteration, residual 1.5·10⁻⁴¹) | [N] |
| λ₂(Ã) | Theorem 6.2: λ₂^even ≥ 2.085·10⁻¹² | 5.140282089·10⁻¹² (not certified) | [N] |
| Cholesky residual at μ = 9·10⁻¹⁸ | r = 1.06·10⁻⁵⁰, slack s = 3.6·10⁻⁴³ | ‖A − μI − L̃L̃ᵀ‖_∞ = 1.16·10⁻⁸⁰, computed **exactly** (no slack needed) | [R] |
| λ_min(A) ≥ | 9·10⁻¹⁸ − 4·10⁻⁴³ | 9·10⁻¹⁸ − 1.9·10⁻³¹ at Zhu's shift; **1.02768955902·10⁻¹⁷ − 1.9·10⁻³¹** at μ = λ₁(1 − 10⁻⁹) | [R] |
| ε_B (block ↔ degrees ≥ 400) | ≤ 10⁻¹⁰⁰ | ≤ 2.2·10⁻¹⁰² | [R] |
| ε_D (tail block − β I) | ≤ 10⁻¹⁰⁰ | ≤ 2.3·10⁻²¹² | [R] |
| **Certified even bound** | **Q ≥ 8.9·10⁻¹⁸ ‖f‖²** | **Q ≥ 1.0276895590·10⁻¹⁷ ‖f‖²** | [R] |

**Convergence of the reduced block in the Legendre degree** [N]: λ₁ of the leading k × k block of A.

| k (max degree) | 20 (38) | 40 (78) | 60 (118) | 80 (158) | 100 (198) | 120–200 (238–398) |
|---|---|---|---|---|---|---|
| even, T♯ = 200 | 1.24830·10⁻¹⁷ | 1.05028·10⁻¹⁷ | 1.03527·10⁻¹⁷ | 1.027789·10⁻¹⁷ | 1.02768956005·10⁻¹⁷ | 1.02768956005·10⁻¹⁷ |
| odd, T♯ = 150 | 1.09626·10⁻¹⁴ | 9.37894·10⁻¹⁵ | 9.11923·10⁻¹⁵ | 9.1183384545·10⁻¹⁵ | 9.1183384545·10⁻¹⁵ | 9.1183384545·10⁻¹⁵ |

The block stabilises to 12 digits once the degree exceeds L·T♯ (160 and 120). That is expected: S_n
is negligible on [0, T♯] for n ≫ LT♯. Zhu's 200 modes leave a wide margin, and that margin is what makes
ε_B ≤ 10⁻¹⁰⁰.

### Odd sector (Section 6.2) and complex f (Corollary 6.3)

| Component | Zhu v2 [P] | This audit | Level |
|---|---|---|---|
| T♯, β* | 150; — | 150; 0.2241180357966231443 | [R] |
| Basis | 200 odd modes, degree 1…399 | same | — |
| Pole sign | −2(∫f sinh(u/2))² | −2ppᵀ | — |
| λ_min(M_odd) | **9.1183·10⁻¹⁵** (inverse iteration, residual 2.3·10⁻³⁰) | **9.118338454501·10⁻¹⁵** (residual 1.3·10⁻⁴¹) | [N] |
| λ₂ | — | 1.7492098·10⁻⁹ | [N] |
| Cholesky at 8.2065·10⁻¹⁵ | succeeds, r = 7.9·10⁻⁵¹ | succeeds, exact residual 8.3·10⁻⁸¹ | [R] |
| ε_Q, c_err, ε_B, ε_D | — | 6.7·10⁻³⁴, 1.3·10⁻³¹, 2.5·10⁻¹⁵³, 4.2·10⁻³¹⁴ | [R] |
| **Certified odd bound** | Q ≥ 8.2065·10⁻¹⁵ ‖f‖² | Q ≥ 9.1183384·10⁻¹⁵ ‖f‖² (and ≥ 8.2065·10⁻¹⁵ at Zhu's shift) | [R] |
| **Complex f** | Q ≥ 8.9·10⁻¹⁸ ‖f‖² | Q ≥ min(even, odd) = **1.0276·10⁻¹⁷ ‖f‖²** | [R] |

The odd value matches Zhu's to five digits. That confirms independently that the matrix assembly
(basis, Ψ_L, the pole sign, the β shift and the quadrature) is the same as Zhu's. For the even sector,
Zhu's 9·10⁻¹⁸ is a chosen Cholesky shift: 0.876 λ_min(A), just as the odd shift 8.2065·10⁻¹⁵ is
0.900 λ_min. **The only number that differs from the paper is the even constant: 1.0277·10⁻¹⁷ here,
against Zhu's shift-limited 8.9·10⁻¹⁸.** The difference lies in the choice of λ₀, not in any component
of the matrix.

### The rigorous steps, briefly

- **Envelope [R].** Binet's second formula gives ψ(z) = log z − 1/(2z) − 2∫₀^∞ s ds/((s² + z²)(e^{2πs} − 1))
  for Re z > 0.
  - Take z = ¼ + iy with y = t/2. Then \|s² + z²\| ≥ \|Im z²\| = y/2, and ∫ s/(e^{2πs} − 1) = 1/24, so the
    remainder is at most 1/(3t).
  - Also Re(1/2z) ≤ 1/(2t²) and log\|z\| ≥ log(t/2).
  - Together, Re ψ ≥ log(t/2) − 1/(2t²) − 1/(3t) ≥ log(t/2) − 1/t for t ≥ ¾.
  - With P_L ≤ A_L, this gives Ψ_L(t) ≥ log(t/2π) − 1/t − A_L, which increases in t. So Ψ_L ≥ β* on [T♯, ∞).
- **Quadrature [R].** On each panel, the integrand G(x) = (Ψ_L − β)S_mS_n/π, mapped to [−1, 1], is
  analytic in the Bernstein ellipse E_ρ and bounded there by M. Its Chebyshev coefficients satisfy
  \|a_k\| ≤ 2Mρ^{−k}.
  - For the rule actually used (stored nodes and weights), the error is at most
    2M[Σ_{k<64}\|def_k\| + (2 + Σw)ρ^{−64}/(1 − 1/ρ)], scaled by h/2 and summed over the panels.
  - def_k = Σw_iT_k(x_i) − ∫T_k are interval enclosures.
  - The bound M follows from |j_n(w)| ≤ e^{\|Im w\|}, ψ(w) = ψ(1 + w) − 1/w with Re w ≥ ⅛, Binet for
    ψ(1 + w), and \|cos(z log n)\| ≤ cosh(¼ log n).
- **Bessel enclosures [R].**
  - Above the turning point (n ≥ ⌊x⌋ + 1, so 2n + 1 > 2x), the ratios r_n = j_n/j_{n−1} lie in (0, 1)
    and satisfy r_n = 1/((2n+1)/x − r_{n+1}). This holds by downward induction from large n, where
    r_n ≈ x/(2n+1).
  - Iterating this map from r ∈ [0, 1] at n = max(2x, 399) + 60 therefore encloses the true ratios.
  - Below the turning point, upward recurrence from j₀ and j₁ is used. Spot checks against
    `mpmath.besselj` confirm containment at x from 10⁻⁴ to 160.
- **Re ψ(¼ + iy) [R].** Shift by 60, then 24 Stirling terms. The remainder is
  \|R_J\| ≤ \|B_{2J+2}\|/(2(J+1)\|w\|^{2J}(Re w)²), from Binet with \|s² + w²\| ≥ (Re w)².
- **Cholesky [R].** L̃ is a floating Cholesky factor of Ã − μI at 80 digits. Its binary entries are
  converted exactly to integers, so L̃L̃ᵀ and the residual are exact rationals. Then
  λ_min(Ã) ≥ μ − ‖residual‖_∞ (Zhu's Lemma 5.2 with s = 0).
- **Coupling [R].** For real t ∈ [0, T♯], with G = 19.41 the bound on \|Ψ_L − β\|:
  - Bounds: \|S_n\| ≤ d_n = √(2L(2n+1))(LT♯)ⁿ/(2n+1)!! and \|p_n\| ≤ √(2L(2n+1))(L/2)ⁿ/(2n+1)!!·e^{(L/2)²/(2(2n+3))}.
  - Entries: \|B_nm\| ≤ 2\|p_n p_m\| + (GT♯/π)d_n√(2L(2m+1)), and similarly for D'.
  - Combine with Schur's test ‖B‖ ≤ (‖B‖₁‖B‖_∞)^{1/2}; the tails sum geometrically (ratio ≤ 0.040).

### Where the paper, as written, falls short of a proof

Each gap below is closed by this audit's re-execution, so none affects the theorem.

- **Lemma 5.1 (ellipse bound).**
  - The bounds \|j_n(z)\| ≤ e^{\|Im z\|} and \|T̂_nT̂_m\| ≤ 4L√(ν_nν_m)e^{2L·0.4} are standard [R].
  - The bound \|Ψ_L − β*\| ≤ 20 on the ρ = 6.55 ellipses is only asserted [P]. The first panel's ellipse
    comes within 0.114 of the pole of ψ(¼ + iz/2) at z = i/2. Crude estimates (1/\|w\| plus Binet) give
    only about 22–32 there, so a sharper local estimate is needed and is not supplied. The true maximum on
    the sampled boundaries is 14.9 [N], so the stated constant is correct.
  - The per-panel formula 4Mρ/((ρ − 1)(ρ⁶⁴ − 1)) agrees, up to a factor below 2, with the Chebyshev-coefficient
    argument used here, for exact Gauss nodes.
- **Assembly (§§5.4, 6.1).**
  - j_n is computed by Miller's backward recurrence, and ψ, the nodes and the weights in 50-digit floating
    point [P]. No enclosure or a-priori error bound is given for these values. Lemma 5.1 bounds only the
    truncation error, and Lemma 5.2 only the Cholesky rounding.
  - What is missing is an evaluation-error term in c_err. At 50 digits it would be of order 10⁻⁴⁵ or smaller, so
    the gap is formal.
- **Theorems 1.2 and 6.3 for f ∈ L².** The explicit formula is invoked for continuous f [P]. For general
  L² f, Q has to be defined by the frequency integral; the proof of Theorem 1.1 then goes through
  unchanged.
- **Theorem 6.2** quotes λ₁^even ≤ 2.523·10⁻¹⁶, against 2.27·10⁻¹⁷ in §11 and the abstract [P].
  Zhu's certified upper bounds were not re-executed [U]. rh2's uncertified 1.6356·10⁻¹⁷ is below both.

### Steps that remain [P] or [U]

- The explicit formula itself: Q(f) = pole + (1/π)∫₀^∞Ψ_L\|F\|² equals Weil's functional. This is
  standard [P], and confirmed numerically on smooth test functions against 300 zeros (§1) [N].
- `mpmath.iv` (not Arb) is trusted for outward rounding of +, ×, ÷, exp, log, sin, cos and √.
- The new scripts were checked against Zhu's odd eigenvalue, against `mpmath.besselj` and
  `mpmath.digamma`, and against rh2's T₀ entry, but have not been reviewed.

## 3. Cross-check against rh2's upper bounds

Each value is an upper bound on λ*(0.8) in its sector, if the arithmetic is accurate [N]. The values
are in Zhu's normalisation, with conversion factor 1.

| N | even, E + 2vvᵀ | odd, E − 2wwᵀ | even minimiser f(L/2)/‖f‖ |
|---|---|---|---|
| 32 | 1.960658·10⁻¹⁷ | 1.919500·10⁻¹⁴ | 1.92·10⁻⁸ |
| 48 | 1.912353·10⁻¹⁷ | 1.703366·10⁻¹⁴ | 1.90·10⁻⁸ |
| 64 | 1.847615·10⁻¹⁷ | 1.676498·10⁻¹⁴ | 1.96·10⁻⁸ |
| 80 | 1.787868·10⁻¹⁷ | 1.671697·10⁻¹⁴ | 1.96·10⁻⁸ |
| 100 | 1.735738·10⁻¹⁷ | 1.666067·10⁻¹⁴ | 1.91·10⁻⁸ |
| 128 | 1.696034·10⁻¹⁷ | 1.656409·10⁻¹⁴ | 1.83·10⁻⁸ |
| 160 | 1.672598·10⁻¹⁷ | 1.644316·10⁻¹⁴ | 1.75·10⁻⁸ |
| 200 | 1.658019·10⁻¹⁷ | 1.631741·10⁻¹⁴ | 1.66·10⁻⁸ |
| 256 | 1.648165·10⁻¹⁷ | 1.617754·10⁻¹⁴ | 1.58·10⁻⁸ |
| 320 | 1.642939·10⁻¹⁷ | 1.606546·10⁻¹⁴ | 1.51·10⁻⁸ |
| 400 | 1.639712·10⁻¹⁷ | 1.596405·10⁻¹⁴ | 1.45·10⁻⁸ |
| 512 | 1.637550·10⁻¹⁷ | 1.586813·10⁻¹⁴ | 1.39·10⁻⁸ |
| 640 | 1.636351·10⁻¹⁷ | 1.579562·10⁻¹⁴ | 1.34·10⁻⁸ |
| 800 | **1.635586·10⁻¹⁷** | 1.573490·10⁻¹⁴ | 1.30·10⁻⁸ |

`scripts/zhu_window_ext.py` uses the same form as `scripts/zhu_window_mp.py` and finds the smallest
eigenvalue by LU and inverse iteration, at 50 digits on 7 processes. N = 800 took 27 minutes. Output:
`data/connes/zhu_window_ext.json`.

- Every value is above 8.9·10⁻¹⁸ and above the certified 1.0277·10⁻¹⁷ (even) and 9.118·10⁻¹⁵ (odd),
  and both sequences decrease monotonically, as finite minima must.
- At 50 digits the N ≤ 100 values reproduce `data/connes/zhu_window.json` (130 digits) to all 20
  printed digits. At N = 256, 70 digits gives the same 20 digits.
- **Convergence.** The brief's "3% per step" at N = 64–100 does not persist.
  - Even-sector differences shrink like N^{−2.0}. A three-point fit λ∞ + CN^{−p} gives p = 2.07 and
    λ∞ = 1.63431·10⁻¹⁷ through N = 400, 512, 640, and p = 2.01 and λ∞ = 1.63423·10⁻¹⁷ through
    N = 512, 640, 800.
  - The even minimiser almost vanishes at the window edge: f(±L/2)/‖f‖ ≈ 1.3·10⁻⁸.
  - The odd sector converges more slowly (p = 0.80, extrapolation 1.542·10⁻¹⁴), plausibly because the
    sine basis forces f(±L/2) = 0.

**An edge-adapted basis** (`scripts/zhu_legendre_upper.py`,
`data/connes/zhu_legendre_upper_{even,odd}.json`; 45 digits, 32 737 frequency nodes on [0, 4000]) [N].
This computes Q itself, not R, on the Legendre basis:

- The log(t/2π) part of the archimedean term is computed exactly from the derivative of the
  Weber–Schafheitlin integral. Spot checks against direct quadrature agree to 10⁻⁷, the accuracy of
  the quadrature.
- The prime comb is computed exactly in u-space, from cross-correlations of Legendre polynomials.
- The remainder ε(t) = Re ψ(¼ + it/2) − log(t/2) is integrated numerically up to T_E = 4000. Beyond T_E
  it is negative and is dropped, so the computed form U satisfies U ≥ Q, and each λ_min(U_N) is again an
  upper bound on λ*. For the minimiser, the overshoot is at most sup_{t≥T_E}\|ε\| × (tail mass) ≤ 5·10⁻²⁹.
- Check: U₀₀ = 0.080274054959356, against rh2's (E + 2vvᵀ)₀₀ = 0.080274054959270. The difference,
  8.6·10⁻¹⁴, is within the dropped-tail bound of 2.6·10⁻¹³.

| Legendre modes (max degree even/odd) | even λ_min(U) | odd λ_min(U) |
|---|---|---|
| 20 (38/39) | 2.04656·10⁻¹⁷ | 1.81692·10⁻¹⁴ |
| 40 (78/79) | 1.67870·10⁻¹⁷ | 1.59346·10⁻¹⁴ |
| 60 (118/119) | 1.65044·10⁻¹⁷ | 1.56627·10⁻¹⁴ |
| 80 (158/159) | 1.64500·10⁻¹⁷ | 1.55773·10⁻¹⁴ |
| 100 (198/199) | 1.64243·10⁻¹⁷ | 1.55410·10⁻¹⁴ |
| 140 (278/279) | 1.63931·10⁻¹⁷ | 1.55105·10⁻¹⁴ |
| 200 (398/399) | **1.637201·10⁻¹⁷** | **1.549161·10⁻¹⁴** |

- **Even sector.** Legendre degree 398 does about as well as cosine N ≈ 550, and both bases decrease
  towards the same limit near 1.634·10⁻¹⁷.
- **Odd sector.** The Legendre basis is clearly better: 1.5492·10⁻¹⁴, against 1.5735·10⁻¹⁴ for the sine
  basis at N = 800, and close to the sine extrapolation of 1.542·10⁻¹⁴.
- No value in either basis comes near 8.9·10⁻¹⁸.

**Resulting bracket for ζ at a = 0.8** (even sector, which is the ground state):

  1.0276·10⁻¹⁷ [R, this audit] ≤ λ*(0.8) ≤ 1.635586·10⁻¹⁷ [N, rh2 cosine N = 800]

The extrapolated value is about 1.634·10⁻¹⁷ [N]. This bracket lies inside Zhu's
[8.9·10⁻¹⁸, 2.27·10⁻¹⁷]. In the odd sector the bracket is 9.1183·10⁻¹⁵ ≤ λ*_odd ≤ 1.549161·10⁻¹⁴
(Legendre degree 399). The x = 5 value in `docs/CONNES_LETTER.md` (1.005·10⁻¹⁷) is not an upper bound
at x = e^{1.6}, and it is not needed here.

## 4. The retracted support-2.38 claim (§7, Remark 3.3)

**What failed [P].** An earlier draft certified support 2.38 (a = 1.19, n ≤ 10). It sharpened the
envelope constant prime by prime: the powers of one prime share a phase, so they cannot all be −1
together. That gives A_eff = −inf_t P_L(t), a **lower** bound for the comb. But Q ≥ R needs
Ψ_L = H − P_L ≥ β on [T♯, ∞), which needs an **upper** bound for P_L. Its supremum is A_L, attained at
t = 0 and approached again at large t, because log 2, log 3, log 5 and log 7 are linearly independent
(Kronecker). §7 itself says that the sharpening bounds the comb in the wrong direction. It keeps a
950 × 950, 70-digit positive-definite matrix only as evidence, not as a certificate.

**Numbers** (`scripts/zhu_retraction_check.py`, `data/connes/zhu_retraction_check.json`; Ψ_L scanned in
double precision on a 0.01-grid of [T♯, 20T♯]) [N]:

| a | constant | A | T₁ = 2πe^A | T♯ | β | min (Ψ_L − β) on the scan |
|---|---|---|---|---|---|---|
| 1.19 | A_eff (retracted) | 4.694796 | 687.2 | 722 | 0.0480 | **−1.79** at t = 897.7: envelope fails |
| 1.19 | A_L (correct) | 7.075006 | 7427 | 7799 | 0.0487 | +0.40 |
| 0.8 | A_L (v2, even) | 2.941974 | 119.1 | 200 | 0.5135 | +0.109 |
| 0.8 | A_L (v2, odd) | 2.941974 | 119.1 | 150 | 0.2241 | +0.076 |
| 0.8 | A_eff (counterfactual) | 2.135002 | 53.1 | 200 | 1.3204 | **−0.70** at t = 217.5: would fail |

The A_eff = 4.69 and A_L = 7.08 at a = 1.19 match §7 [P]. With the correct constant, a = 1.19 needs
T♯ > 7427, so Legendre degree about L·T♯ ≈ 9000, which is the doubly exponential barrier of Theorem 1.4.

**Does the failure mode affect a = 0.8? No [R].** The v2 certificate at a = 0.8 uses A_L, not A_eff:
the β* = 0.5134667… that Zhu prints equals log(200/2π) − 1/200 − A_L with A_L = 2.9419735 exactly, and
this audit's re-execution uses the same constant. With A_L the envelope is proved (Lemma 3.1 plus
sup P_L = A_L), so Q ≥ R holds on [200, ∞). The table shows the per-prime sharpening *would* break a = 0.8
as well. It is absent from the v2 a = 0.8 certificate, and no other step of that certificate uses a
lower bound where an upper bound is needed.

## Reproduce

```bash
PY=<repo>/.venv/bin/python   # mpmath 1.4.1; no numpy, no flint
$PY scripts/zhu_certificate_iv.py --sector even --tsharp 200 --nmodes 200 --workers 8 --mu 9e-18,auto --blocks 20,40,60,80,100,120,160 --json data/connes/zhu_certificate_even.json
$PY scripts/zhu_certificate_iv.py --sector odd  --tsharp 150 --nmodes 200 --workers 8 --mu 8.2065e-15,auto --blocks 20,40,60,80,100,120,160 --json data/connes/zhu_certificate_odd.json
$PY scripts/zhu_normalisation_check.py --zeros 300 --n 120 --dps 40 --json data/connes/zhu_normalisation.json
$PY scripts/zhu_window_ext.py --n 32,48,64,80,100,128,160,200,256,320,400,512,640,800 --dps 50 --workers 7 --json data/connes/zhu_window_ext.json
$PY scripts/zhu_legendre_upper.py --sector even --nmodes 200 --te 4000 --dps 45 --workers 4 --json data/connes/zhu_legendre_upper_even.json
$PY scripts/zhu_legendre_upper.py --sector odd  --nmodes 200 --te 4000 --dps 45 --workers 5 --json data/connes/zhu_legendre_upper_odd.json
$PY scripts/zhu_retraction_check.py --json data/connes/zhu_retraction_check.json
```

Logs go to `logs/` through `scripts/progress.py`.

## Five-line summary

1. Normalisation: rh2's zeros-side QW at x = e^{1.6} is exactly Zhu's Q on [−0.8, 0.8] (same L² norm, pole included, both parities), so the conversion factor is 1. Confirmed three ways to ≤ 2·10⁻¹⁴ [N].
2. Certificate: an interval re-execution of Theorem 1.1 at Zhu's parameters certifies λ_min ≥ 1.0277·10⁻¹⁷ (even) and 9.118·10⁻¹⁵ (odd, matching Zhu's 9.1183·10⁻¹⁵), so Q ≥ 1.0276·10⁻¹⁷ ‖f‖² ≥ 8.9·10⁻¹⁸ ‖f‖² for complex f [R; mpmath.iv, unreviewed code].
3. Upper bounds: rh2's cosine minima, extended to N = 800, and a new Legendre basis keep falling (≈ N⁻², limit ≈ 1.634·10⁻¹⁷) but stay ≥ 1.6356·10⁻¹⁷, which brackets λ*(0.8) in [1.0276, 1.6356]·10⁻¹⁷ [R / N].
4. Retraction: the support-2.38 draft used a lower bound for the prime comb (A_eff = 4.69) where an upper bound (A_L = 7.08) is needed, and the envelope fails at t ≈ 898 [N]. The a = 0.8 certificate uses A_L, so it is unaffected [R].
5. Verdict: verified with corrections. The paper leaves its integrand values unenclosed and asserts the ellipse constant 20 without proof, and its Theorem 6.2 upper bound is inconsistent with §11. None of this changes 8.9·10⁻¹⁸.

## Review of this audit (main session, 2026-10-01)

Every check below passed. The audit is accepted into `feat/weil-gram-instrument`.

| Check | Method | Result |
|---|---|---|
| **Paper numbers** | Second, independent WebFetch read of the v2 HTML (also model-mediated; the PDF was not downloaded) | Matches the audit on 8.9·10⁻¹⁸, λ₀ = 9·10⁻¹⁸ (r = 1.06·10⁻⁵⁰), λ_min(M_odd) = 9.1183·10⁻¹⁵, the odd shift 8.2065·10⁻¹⁵, T♯ = 200/150, ρ = 6.55, \|Ψ_L − β*\| ≤ 20, β* = 0.5134667…, A_L = 2.9419735…, A_eff = 4.6948 against 7.0750, and 2.27·10⁻¹⁷. Theorem 6.2's bound reads 2.5223·10⁻¹⁶ (the audit rounds to 2.523). The paper assembles the odd sector with 160-point Gauss; the audit uses 32-point Gauss with its own enclosed error, and the eigenvalue still matches. |
| **Normalisation, factor 1** | Own code, not the audit's: Q(T₀) for T₀ = 1/√(1.6) on [−0.8, 0.8], from Zhu's frequency formula, integrated to t = 20000 plus an averaged tail | 0.0802740550566, against rh2's (E + 2vvᵀ)₀₀ = 0.0802740549593. The difference, 1·10⁻¹⁰, is within the neglected oscillatory tail. |
| **Certificate reproduces** | Reran `zhu_certificate_iv.py` for both sectors with 4 workers | Even: λ₁ = 1.02768956005196·10⁻¹⁷; Cholesky at 1.02768955902·10⁻¹⁷ succeeds, residual 1.15·10⁻⁸⁰. Odd: 9.11833845450086·10⁻¹⁵; both shifts succeed. Identical to the committed JSON. |
| **Bound assembly** | Read the code | λ_A ≥ μ − ‖residual‖_∞ − N·(entry radius) − N·ε_Q, then min(λ_A, β̃ − ε_D) − ε_B. That is the correct Schur-type split, and ‖R‖₂ ≤ ‖R‖_∞ holds for the symmetric residual. |
| **Ellipse constant G = 19.41** | Re-derived each step | It is a uniform analytic bound over the whole ellipse, not a sample: ψ(w) = ψ(1 + w) − 1/w with Re w ≥ ⅛; \|log v\| ≤ log\|v\| + π/2; Binet remainder ≤ 1/(12(Re v)²) from \|s² + v²\| ≥ (Re v)²; \|j_n(z)\| ≤ e^{\|Im z\|} from j_n = ½(−i)ⁿ∫e^{izx}P_n(x)dx. |
| **Envelope (t ≥ ¾)** | Re-derived | Correct: remainder ≤ 1/(3t), Re(1/2z) ≤ 1/(2t²), and 1/(2t²) + 1/(3t) ≤ 1/t exactly when t ≥ ¾. |
| **rh2 extension** | `scripts/zhu_window_mp.py` (full `mp.eigsy`, not the audit's inverse iteration) at N = 200, 60 digits | Even 1.65801850225·10⁻¹⁷, odd 1.63174134563·10⁻¹⁴; identical to `zhu_window_ext.json`. |
| **Legendre upper bound** | Reran the committed script at 20 modes | 2.04655644·10⁻¹⁷, as tabulated. The first even run in the audit's logs (15:20) was broken: λ_min ≈ −10⁴⁴ and a negative tail mass. It was fixed before the committed run at 15:24, but the doc does not mention it. |

**Remaining trust assumptions.**
- `mpmath.iv` rounds outward correctly.
- The paper was read only through model-mediated extractions.
- No one has reviewed `zhu_certificate_iv.py` line by line beyond the parts above: the Bessel ratio enclosure, the fixed-point conversion and the defect sums were not re-derived here.

**Consequence for rh2.** Zhu's window is the first support at which rh2 has a two-sided bracket of the Weil minimum: 1.0277·10⁻¹⁷ ≤ λ*(0.8) ≤ 1.6356·10⁻¹⁷.
- rh2's N = 100 value of 1.736·10⁻¹⁷ is only about 6% above the extrapolated limit of 1.634·10⁻¹⁷.
- That is a calibration point for how far rh2's finite-basis upper bounds sit above the truth.
- It also corrects the 9bd1805 brief's remark that the trigonometric basis "converges slowly" here: past N ≈ 400 the even values converge like N⁻².
