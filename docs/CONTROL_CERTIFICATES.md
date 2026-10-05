# Control Certificates: Two-Sided Rigorous Statements for F_{t*}, DH and Z₁

Date: 2026-10-01. Worktree `<review-worktree>`, branch `cert/controls`, base `201d64a`
(`feat/weil-gram-instrument`). Brief: `docs/BRIEF_CONTROL_CERTIFICATES.md`. Nothing outside the worktree
was written. Every run used at most 2 worker processes.

**Rigour labels.** **[R]** re-derived rigorously: Arb ball arithmetic (python-flint 0.9.0), exact integer or
rational arithmetic, or a written argument given here. **[N]** reproduced numerically (high precision, no
enclosure). **[P]** taken from a source. **[U]** could not verify. [R] means rigorous *given* that Arb and
python-flint are correct and that the new scripts are correct; the scripts have not been reviewed by a
second person.

**What "positive at a" means.** Q(f) ≥ c‖f‖₂² with c > 0 for every complex f ∈ L² supported in [−a, a],
i.e. Our window at x = e^{2a}. This needs both parity sectors (§1.1), so every function was certified in
both, including DH and Z₁, for which the brief lists only the even sector. Positivity at a implies
positivity, with the same c, at every smaller support.

## Result

| Function | Certified positive [R] (both sectors) | Certified negative [R] (explicit witness) | Our finite-basis crossover [N] |
|---|---|---|---|
| F_{t*} | x = e^{1.67} = 5.31217, Q ≥ 1.1666·10⁻³‖f‖² | x = 5.36250 (Q/‖f‖² ≤ −9.889·10⁻⁷); x = 7 (≤ −0.045280) | N = 32: (5.366, 5.474); N = 64: (5.36248, 5.36250] |
| DH | x = e^{2.38} = 10.80490, Q ≥ 1.0441·10⁻⁸‖f‖² | x = 30.74497 (≤ −2.292·10⁻³²); x = 32 (≤ −1.3999·10⁻²⁹) | N = 64: ≈ 30.83 (this work); N = 80: (30.74404, 30.74497] |
| Z₁ | x = e^{2.76} = 15.79984, Q ≥ 0.049235‖f‖² (T♯ = 1000); at T♯ = 400, x = e^{2.62} = 13.7357, Q ≥ 0.11169‖f‖² | x = 19.84390 (≤ −1.5342·10⁻⁷); x = 20 (≤ −1.3261·10⁻⁴) | N = 48–80: (17.5, 20] (this work); N = 80: (19.84375, 19.84390] |

So the first full-space failure x_c of each function (the infimum of the x at which Q is not ≥ 0) lies,
rigorously, in
- F_{t*}: [5.31217, 5.36250), a window of 1%;
- DH: [10.805, 30.745);
- Z₁: [15.800, 19.844).

The given witnesses (x = 7, 32, 20) are certified, and so are three new ones at the finite-basis
crossovers. For F_{t*} the positive side reaches to within 1% of the crossover. For DH and Z₁ it is stopped
by T♯ (the size of the comb A), as Zhu's Theorem 1.4 predicts, not by the margin.

## 1. The generalised reduction

### 1.1 Setting and the form

Each control F has real Dirichlet coefficients a_n with a_1 = 1, and −F′/F = Σ c_n n^{−s}, with c_n from
a_n log n = Σ_{d|n} c_d a_{n/d}. For f supported in [−a, a] (our window [−L/2, L/2], L = 2a = log x),
F̂(t) = ∫ f e^{itu} du and g(y) = ∫ f(v) f(v + y) dv. For real f:

  Q(f) = pole(f) + (1/π) ∫_0^∞ Ψ(t) |F̂(t)|² dt,
  Ψ(t) = Σ_κ Re ψ(κ + it/2) + K − Σ_{log n < 2a} (2c_n/√n) cos(t log n),
  pole(f) = 2F̂(i/2)F̂(−i/2): +2(∫ f cosh(u/2))² for even f, −2(∫ f sinh(u/2))² for odd f; 0 if F has no pole.

| Function | Λ(s) | κ | K | Pole at s = 1 | a_n |
|---|---|---|---|---|---|
| ζ (regression) | π^{−s/2}Γ(s/2)ζ(s) | ¼ | −log π | yes | 1 |
| F_{t*} | (5/π)^{s/2}Γ(s/2)F_{t*}(s) | ¼ | log(5/π) | yes | (1 − t*)χ₅(n) + t*(1 + √5·[5 ∣ n]) |
| DH | (5/π)^{s/2}Γ((s+1)/2)f(s) | ¾ | log(5/π) | no | Re((1 − iκ)χ(n)), χ(2) = i mod 5 |
| Z₁ | 20^{s/2}(2π)^{−s}Γ(s)Z₁(s) | ¼ and ¾ | log 20 − 2 log π | yes | r(n)/2, r(n) = #{x² + 5y² = n} |

t* = L(¾, χ₅)/(L(¾, χ₅) − G(¾)) is enclosed in Arb through Hurwitz zeta values: [0.05529219744234932650223828 ± 10⁻³⁷] [R].
κ_DH = (√(10 − 2√5) − 2)/(√5 − 1) is exact algebraic. Coefficients of Z₁ are exact rationals. The c_n are balls;
exact zeros (DH: n ≡ 0 mod 5; Z₁: most n) drop out of the comb.

That this Q is Weil's functional, equal to the zeros sum Σ_ρ ĝ(γ_ρ) with real and off-line zeros
included, is the explicit formula for Λ with these gamma factors [P: standard; for ζ it is the formula Zhu's reduction starts from, `docs/AUDIT_ZHU.md` §1].
Every certificate below is a statement about this Q, and §2 shows that it is the form this work computes.

Two facts hold for every function [R, algebraic, as in `docs/AUDIT_ZHU.md` §1]. The kernel is real and
symmetric, so Q(a + ib) = Q(a) + Q(b) for real a, b. For real f = e + o (even plus odd) there is no cross
term: F̂_e is real, F̂_o is imaginary, and the pole term is 2(E² − O²) with E = ∫ e cosh(u/2),
O = ∫ o sinh(u/2). Hence inf Q(f)/‖f‖² over all f supported in [−a, a] is min(inf_even, inf_odd).

### 1.2 Positive side: Zhu's Theorem 1.1, unchanged

If Ψ ≥ β on [T♯, ∞), then Ψ|F̂|² ≥ β|F̂|² pointwise there, and Plancherel ‖f‖² = (1/π)∫_0^∞|F̂|² gives
Q(f) ≥ R(f) = pole + (1/π)∫_0^{T♯}(Ψ − β)|F̂|² + β‖f‖² [R]. The proof uses nothing about ζ. The Legendre
basis T_n(u) = P̄_n(u/a)/√a, S_n(t) = √(2a(2n+1)) j_n(at), the pole vector p_n = ±√(2a(2n+1)) i_n(a/2), the
sign-adjusted block A = s·2ppᵀ + Σ_q c_q S(t_q)S(t_q)ᵀ + βI and the Schur split
Q ≥ (min(λ_min(A), β − ε_D) − ε_B)‖f‖² are the audit's (§2), with a and Ψ replaced. For DH, p = 0.

**Comb bound** [R]. Σ (2c_n/√n) cos(t log n) ≤ A := Σ_{log n < 2a} 2|c_n|/√n for all real t, because
|cos| ≤ 1. This is the upper bound the reduction needs. (Audit §4: a lower bound in this place was the
retracted error.)

**Binet envelope for general κ ∈ (0, 1)** [R]. Binet's second formula,
ψ(z) = log z − 1/(2z) − 2∫_0^∞ s ds/((s² + z²)(e^{2πs} − 1)) for Re z > 0, at z = κ + iy with y = t/2 > 0:
- Re log z = log|z| ≥ log y;
- Re 1/(2z) = κ/(2|z|²) ≤ κ/(2y²);
- |s² + z²| ≥ |Im z²| = 2κy, and ∫_0^∞ s/(e^{2πs} − 1) ds = 1/24, so the remainder is at most 1/(24κy).

So Re ψ(κ + it/2) ≥ log(t/2) − 2κ/t² − 1/(12κt). The right side is ≥ log(t/2) − 1/t exactly when
t ≥ 24κ²/(12κ − 1):
- κ = ¼ (ζ, F_{t*}): log(t/2) − 1/(2t²) − 1/(3t) ≥ log(t/2) − 1/t for t ≥ ¾, the audit's envelope;
- κ = ¾ (DH): log(t/2) − 3/(2t²) − 1/(9t) ≥ log(t/2) − 1/t for t ≥ 27/16;
- degree 2 (Z₁): the sum, Re ψ(¼ + it/2) + Re ψ(¾ + it/2) ≥ 2 log(t/2) − 2/t for t ≥ 27/16. (Duplication,
  as 2 Re ψ(½ + it) − 2 log 2, gives the sharper 2 log(t/2) − 1/(2t²) − 1/(6t); it is not used.)

Each bound increases in t. With n_κ the number of gamma factors,

  Ψ(t) ≥ β* := n_κ [log(T♯/2) − 1/T♯] + K − A   on [T♯, ∞),

and the certificate uses a dyadic β̃ ≤ β*. For ζ this is the audit's β* = log(T♯/2π) − 1/T♯ − A_L.
β* > 0 needs T♯ > T₁ ≈ 2e^{(A − K)/n_κ}; this is what limits the support (§3).

**Ellipse bound, generalised** [R]. On each panel's Bernstein ellipse (z = t + is, |s| ≤ b = ¼):
- w = κ ± iz/2 has Re w ≥ κ − b/2 > 0;
- ψ(w) = ψ(1 + w) − 1/w, and |ψ(v)| ≤ log|v| + π/2 + 1/(2 Re v) + 1/(12(Re v)²) for Re v ≥ 1;
- |cos(z log n)| ≤ cosh(b log n), and |j_n(az)| ≤ e^{ab}.

So |Ψ − β| ≤ G := |K| + |β| + Σ_κ [1/(κ − b/2) + log|v|_max + π/2 + 1/(2 Re v) + 1/(12 (Re v)²)]
+ Σ (2|c_n|/√n) cosh(b log n), and the panel integrands are bounded by M = G·2a(2n_max + 1)e^{2ab}/π.

**Quadrature of the exact Gauss rule** [R]. Nodes and weights are Arb enclosures of the exact 32-point
Gauss–Legendre rule (`arb.legendre_p_root`), so the computed sums enclose the exact-rule sums. For f
analytic in E_ρ with |f| ≤ M, |a_k| ≤ 2Mρ^{−k}. The exact rule is exact for k < 64 and kills odd k by
symmetry, |∫T_k| ≤ 2/(k² − 1) and |Σ w_i T_k(x_i)| ≤ 2. Summed over the panels,

  ε_Q = T♯·M·(2 + 2/(4·32² − 1))·ρ^{−64}/(1 − ρ^{−2}),   c_err = N·ε_Q.

The audit instead bounded its stored floating rule through the exactness defects.

**Spherical Bessel enclosures** [R, given Pincherle's theorem].
- For n > x, the ratio r_n = j_n(x)/j_{n−1}(x) is the value of the continued fraction
  1/((2n+1)/x − 1/((2n+3)/x − …)) (Pincherle; j_n is the minimal solution of the recurrence).
- Every partial denominator is ≥ 2, so every approximant lies in (0, 1].
- The exact identity r_n = 1/((2n+1)/x − r_{n+1}), started from the ball [0, 1] at n = max(n_max, 2x) + 60,
  therefore encloses every r_n down to n = ⌊x⌋ + 1. Below that, the upward recurrence from j_0, j_1 is exact.
- Ball arithmetic overestimates the upward recurrence by about 0.7 bits per unit of x, so nodes and
  Bessel values are computed with 0.9·a·T♯ + 64 extra bits and rounded afterwards.

**λ_min(A)** [R]. A floating Cholesky factor L̃ of mid(A) − μI (gmpy2, 256 bits) is converted exactly to Arb.
The ball residual R = A − μI − L̃L̃ᵀ gives λ_min(A′) ≥ μ − max_i Σ_j |R_ij| for every symmetric A′ in the
ball matrix A. The ball A already holds every entry's enclosure radius. μ = λ₁(mid A)(1 − 10⁻⁹) by
inverse iteration; where an eigenvalue cluster at β̃ stalls the iteration (Z₁ odd at a = 1.1 and 1.2), μ
is stepped down until the factor exists.

### 1.3 Negative side: a u-space formula for Q

**Identity** [R, written argument]. For Re κ > 0, the series ψ(z) = −γ + Σ_{m≥0}[1/(m+1) − 1/(m+z)] gives
Re ψ(κ + it/2) − ψ(κ) = Σ_m [1/(κ+m) − (κ+m)/((κ+m)² + t²/4)], and each term equals
∫_0^∞ 2e^{−2(κ+m)u}(1 − cos tu) du. All terms are non-negative, so by monotone convergence

  Re ψ(κ + it/2) = ψ(κ) + ∫_0^∞ φ_κ(u)(1 − cos tu) du/u,   φ_κ(u) = 2u e^{−2κu}/(1 − e^{−2u}) = u e^{(1−2κ)u}/sinh u.

Multiply by |F̂(t)|²/π and integrate over t > 0. The integrand is non-negative, so Tonelli allows the swap,
and (1/π)∫_0^∞|F̂|²(1 − cos tu) dt = g(0) − g(u) for f ∈ L². Since g = 0 on [L, ∞) and
∫_L^∞ φ_κ(u) du/u = Σ_m e^{−2(κ+m)L}/(κ + m),

  Q(f) = pole(f) + g(0)[K + Σ_κ (ψ(κ) + Σ_{m≥0} e^{−2(κ+m)L}/(κ+m))]
         + ∫_0^L ((g(0) − g(u))/u) Σ_κ φ_κ(u) du − Σ_{log n < L} (2c_n/√n) g(log n),

with both sides in (−∞, +∞] for every f ∈ L²[−L/2, L/2]. If x = e^L is an integer, n = x drops out because
g(L) = 0. For ζ, φ_{1/4}(u)/u = e^{u/2}/sinh u is 2K(u) of `connes_letter_mp.py`; for DH,
φ_{3/4} − φ_{1/4} = −u sech(u/2) is its sech term.

**Closed form of g** [R, algebra; checked against quadrature]. For an even polynomial
f = A_0 + Σ_{k≥1} A_k cos ω_k u on [−L/2, L/2], ω_k = 2πk/L,

  g(u) = (L − u)[A_0² + ½Σ A_k² cos ω_k u] + Σ_{m≥1} σ_m sin ω_m u   (0 ≤ u ≤ L),
  σ_m = 2A_m S_m,   S_m = −A_0(−1)^m/ω_m + ½Σ_{k≥1} A_k(−1)^{k+m}[[k ≠ m]/(ω_k − ω_m) − 1/(ω_k + ω_m)].

For an odd polynomial f = Σ A_k sin ω_k u the first bracket is ½Σ A_k² cos ω_k u, and
σ_m = A_m Σ_{k≥1} A_k(−1)^{k+m}[[k ≠ m]/(ω_k − ω_m) + 1/(ω_k + ω_m)]. Then (g(0) − g(u))/u is entire, and the
integrand is analytic near [0, L] (u/sinh u has poles only at ±iπ, ±2iπ, …), so `acb.integral` encloses
it. The pole term uses ∫ cos(ω_k u) e^{su} du = 2(−1)^k sinh(sL/2)·s/(s² + ω_k²) and
∫ sin(ω_k u) e^{su} du = −2(−1)^k sinh(sL/2)·ω_k/(s² + ω_k²).

## 2. Normalisation checks

`scripts/control_normalisation.py`, `data/controls/normalisation.json`. Three independent evaluations of
the same numbers:
- the u-space formula of §1.3 in Arb (192 bits) [R];
- the frequency side (1/π)∫_0^{2000} Ψ|F̂|² by `acb.integral`, plus pole, plus a rigorous tail bound from
  |F̂(t)| ≤ ‖f⁗‖₁/t⁴ (≤ 10⁻¹⁵ here) [R];
- Our own matrix, `conductor5_family_mp.forms`, `connes_letter_mp.build_form(…, "dh")` and
  `epstein_connes_mp.forms("Z1", …)`, at 40 digits [N].

The test function is (1 + cos θ)² = 3/2 + 2cos θ + ½cos 2θ (even) or sin θ(1 + cos θ)² = 5/4 sin θ + sin 2θ +
¼ sin 3θ (odd), θ = 2πu/L. Both vanish to 4th order at ±L/2 and are exact in our basis, so our value
has no truncation error. The entry is T₀ = 1/√(2a) (even) or b₁ = √(2/L) sin(2πu/L) (odd). For Z₁'s odd
sector this work has no routine; it was assembled from `build_form(…, parity="odd")` pieces exactly as
`epstein_connes_mp.forms` does, with −2wwᵀ.

| Function | Sector | x | Entry, u-space [R] | Entry, this work [N] | Test fn, u-space [R] | Test fn, frequency side [R] | Test fn, this work [N] | \|u-space − this work\| |
|---|---|---|---|---|---|---|---|---|
| ζ | even | e^{1.6} | 0.0802740549592696 | 0.0802740549592696 | 1.15222457863408·10⁻⁴ | [1.152224578634·10⁻⁴ ± 8·10⁻¹⁷] | 1.15222457863408·10⁻⁴ | 8.9·10⁻⁴⁰ |
| ζ | odd | e^{1.6} | 0.00560377065319998 | 0.00560377065319998 | 0.00235976564834408 | [0.00235976564834 ± 5·10⁻¹⁵] | 0.00235976564834408 | 9.4·10⁻⁴¹ |
| F_{t*} | even | 7 | 4.05245867755139 | 4.05245867755139 | 17.5093402555127 | [17.50934025551270 ± 5·10⁻¹⁷] | 17.5093402555127 | 4.3·10⁻³⁹ |
| F_{t*} | odd | 7 | −0.0407491455528991 | −0.0407491455528991 | 1.82471922423064 | [1.824719224230635 ± 5·10⁻¹⁶] | 1.82471922423064 | 9.0·10⁻⁴¹ |
| F_{t*} | even | e^{1.7} | 3.63807311596912 | 3.63807311596912 | 13.8688426265302 | [13.86884262653022 ± 5·10⁻¹⁷] | 13.8688426265302 | 3.7·10⁻³⁹ |
| F_{t*} | odd | e^{1.7} | 0.045497925446554 | 0.045497925446554 | 2.24364892609466 | [2.243648926094663 ± 7·10⁻¹⁶] | 2.24364892609466 | 5.4·10⁻⁴⁰ |
| DH | even | 32 | 0.0899091651355688 | 0.0899091651355688 | 0.00486618438060776 | [0.00486618438060776 ± 5·10⁻¹⁸] | 0.00486618438060776 | 2.8·10⁻⁴⁰ |
| DH | even | e² | 0.239207058035372 | 0.239207058035372 | 1.91726210085129 | [1.917262100851288 ± 5·10⁻¹⁷] | 1.91726210085129 | 2.9·10⁻⁴⁰ |
| DH | odd | e² | 0.270865878491405 | 0.270865878491405 | 3.33406161996123 | [3.33406161996123 ± 3·10⁻¹⁶] | 3.33406161996123 | 2.3·10⁻⁴⁰ |
| Z₁ | even | 20 | 1.14017376166522 | 1.14017376166522 | 27.523336231247 | [27.52333623124700 ± 8·10⁻¹⁸] | 27.523336231247 | 3.5·10⁻³⁹ |
| Z₁ | even | e^{2.4} | 1.72268123877973 | 1.72268123877973 | 20.8911769100259 | [20.89117691002585 ± 5·10⁻¹⁷] | 20.8911769100259 | 2.2·10⁻³⁹ |
| Z₁ | odd | e^{2.4} | 2.00160350285505 | 2.00160350285505 | 5.32937848319338 | [5.329378483193379 ± 3·10⁻¹⁶] | 5.32937848319338 | 2.3·10⁻³⁹ |

**Reading.**
- **No mismatch.** For every function and sector, our zeros-side matrix equals the Q of §1.1, with
  factor 1, to our 40-digit working precision.
- **The rigorous checks overlap.** The frequency-side enclosure (built from the same κ, K, comb and pole
  data as the positive-side code) contains the u-space value in every row. That checks the generalised Ψ
  independently of the u-space derivation, and the derivation independently of this work.
- **ζ.** The T₀ entry 0.0802740549592696 is the audit's (E + 2vvᵀ)₀₀ = 0.0802740549593.
- **The witnesses** (§4) are a third check: the u-space value of each stored witness equals our Rayleigh
  quotient to all 20 printed digits.
- **The positive-side assembly.** `scripts/control_entry_check.py` recomputes one entry of the certificate
  matrix, A₀₀ = R(T₀) or R(T₁), with `acb.integral` instead of the Gauss panels and Bessel recurrences.
  Seven certificates were checked: F_{t*} odd at 0.835 and even at 1.09, DH even and odd at 1.19, Z₁ even
  and odd at 1.31, and Z₁ even at 1.38. In each the two enclosures agree to all 25 printed digits [R]
  (`data/controls/entry_checks.json`). The two ζ JSONs predate the stored A₀₀, and ζ is covered by the
  regression (§3.1).
- **The positive side's Ψ against the u-space Q.** For the even mode, Q(T₀) − A₀₀ =
  (1/π)∫_{T♯}^∞ (Ψ − β̃)|T̂₀|² must be positive and close to (n_κ + A)/(πaT♯). It is 0.0044164 (F_{t*}, a = 1.09),
  0.0044208 (DH, 1.19), 0.0073832 (Z₁, 1.31), 0.0034512 (Z₁, 1.38) and 0.0078749 (ζ, 0.8), against the
  asymptotic values 0.004420, 0.004424, 0.007378, 0.003452 and 0.007842 [N]. That covers every κ-set.

## 3. Positive side

`scripts/control_certificate.py` (one run), `scripts/control_cert_grid.py` (grids),
`data/controls/grid_*.json` and `data/controls/certs/*.json`. Fixed choices: 32-point Gauss on panels of
width ¼ (4T♯ panels), ρ = 2 + √5 (b = ¼), 256-bit balls, Legendre degree ≥ 1.5·a·T♯ + 80. A run counts as
certified when min(λ_min(A), β̃ − ε_D) − ε_B > 0.

### 3.1 Regression: ζ at a = 0.8

| | Audit (`zhu_certificate_iv.py`, mpmath.iv) | This code (Arb) |
|---|---|---|
| even λ₁(mid A), T♯ = 200, 200 modes | 1.027689560051963·10⁻¹⁷ | 1.02768956005196·10⁻¹⁷ |
| even certified, μ = 1.02768955902·10⁻¹⁷ | 1.0276895590199816·10⁻¹⁷ | 1.0276895590199925·10⁻¹⁷ |
| even certified, μ = 9·10⁻¹⁸ | 8.99999999999982·10⁻¹⁸ | 8.99999999999993·10⁻¹⁸ |
| odd λ₁(mid A), T♯ = 150, 200 modes | 9.118338454501·10⁻¹⁵ | 9.11833845450086·10⁻¹⁵ |
| odd certified, μ = λ₁(1 − 10⁻⁹) | 9.11833844538·10⁻¹⁵ | 9.11833844538·10⁻¹⁵ |
| leading blocks k = 20, 40, 60, 80 (even) | 1.24830, 1.05028, 1.03527, 1.027789 (·10⁻¹⁷) | 1.24830, 1.05028, 1.03527, 1.027789 (·10⁻¹⁷) |
| A_L, β* even / odd | 2.9419735252236204555, 0.51346677491507 / 0.22411803579662 | identical |
| G, M (even) | 19.406408, 11751.448 | identical |
| ε_B, ε_D (even) | 2.2154·10⁻¹⁰², 2.3364·10⁻²¹² | identical |
| ε_Q, c_err (even) | 9.2151·10⁻³⁴, 1.843·10⁻³¹ | 3.7285·10⁻³⁴, 7.457·10⁻³² (exact-rule bound, §1.2) |
| ε_Q, c_err (odd) | 6.7244·10⁻³⁴, 1.3449·10⁻³¹ | 2.7207·10⁻³⁴, 5.4415·10⁻³² |
| time | 76 s, 8 processes | 13 s, 2 processes |

The regression passes: the even and odd targets of the brief, 1.02768955902·10⁻¹⁷ and 9.11833845·10⁻¹⁵,
are reproduced, and every matrix ingredient agrees to all printed digits [R]. The certified numbers differ
in the 14th digit only because the quadrature bound is computed differently.

### 3.2 Grids

The columns are: support a and x = e^{2a}; T♯; Legendre modes (degrees); A; β̃; the ellipse bounds G and M;
the error terms ε_Q, c_err = Nε_Q, ε_B, ε_D; λ₁ of the midpoint block [N]; the certified c [R]; and our
finite-basis minimum at the same x, N = 32 / 64 (an upper bound on the true minimum) [N]; and the run time
in seconds (s).

**F_{t*}, odd sector** (`data/controls/grid_ftstar_odd.json`)

| a | x = e^{2a} | T♯ | modes (degrees) | A | β̃ | G | M | ε_Q | c_err | ε_B | ε_D | λ₁(mid A) [N] | certified c [R] | Our N = 32 / 64 [N] | s |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 0.6 | 3.32012 | 400 | 221 (1–441) | 2.00014 | 3.76038 | 21.67 | 9.87e+03 | 6.3e-34 | 1.4e-31 | 3.8e-55 | 3.4e-118 | 0.54648626 | 0.54648626 | 0.562832 / 0.556132 | 39 |
| 0.7 | 4.0552 | 400 | 251 (1–501) | 2.83811 | 2.92241 | 21.72 | 1.38e+04 | 8.7e-34 | 2.2e-31 | 4.3e-57 | 2.8e-122 | 0.17066476 | 0.17066476 | 0.180371 / 0.176419 | 44 |
| 0.8 | 4.95303 | 400 | 281 (1–561) | 2.83811 | 2.92241 | 21.72 | 1.85e+04 | 1.2e-33 | 3.3e-31 | 3.5e-59 | 1.3e-126 | 0.029092637 | 0.029092637 | 0.0337654 / 0.0318202 | 65 |
| 0.81 | 5.05309 | 400 | 284 (1–567) | 3.09569 | 2.66484 | 21.74 | 1.91e+04 | 1.2e-33 | 3.4e-31 | 2.1e-59 | 4.7e-127 | 0.020009132 | 0.020009132 | 0.0242917 / 0.0224413 | 53 |
| 0.82 | 5.15517 | 400 | 287 (1–573) | 3.09569 | 2.66484 | 21.74 | 1.96e+04 | 1.2e-33 | 3.6e-31 | 1.3e-59 | 1.7e-127 | 0.011983979 | 0.011983979 | 0.0158534 / 0.01428 | 82 |
| 0.83 | 5.25931 | 400 | 290 (1–579) | 3.09569 | 2.66484 | 21.74 | 2.02e+04 | 1.3e-33 | 3.7e-31 | 7.8e-60 | 6e-128 | 0.0046263155 | 0.0046263155 | 0.00819929 / 0.00674786 | 74 |
| 0.835 | 5.31217 | 400 | 291 (1–581) | 3.09569 | 2.66484 | 21.74 | 2.04e+04 | 1.3e-33 | 3.8e-31 | 2.1e-59 | 4.4e-127 | 0.0011665772 | 0.0011665772 | 0.004611 / 0.00320548 | 56 |
| 0.84 | 5.36556 | 400 | 293 (1–585) | 3.09569 | 2.66484 | 21.74 | 2.07e+04 | 1.3e-33 | 3.9e-31 | 4.8e-60 | 2.1e-128 | not PD | fails | 0.00117747 / -0.000190821 | 54 |
| 0.84 | 5.36556 | 800 | 545 (1–1089) | 3.09569 | 3.35923 | 23.12 | 4.1e+04 | 5.2e-33 | 2.8e-30 | 8.5e-80 | 9.3e-169 | not PD | fails | 0.00117747 / -0.000190821 | 226 |
| 0.85 | 5.47395 | 400 | 296 (1–591) | 3.09569 | 2.66484 | 21.74 | 2.13e+04 | 1.4e-33 | 4e-31 | 2.9e-60 | 7.5e-129 | not PD | fails | — | 73 |

**F_{t*}, even sector** (`data/controls/grid_ftstar_even.json`)

| a | x = e^{2a} | T♯ | modes (degrees) | A | β̃ | G | M | ε_Q | c_err | ε_B | ε_D | λ₁(mid A) [N] | certified c [R] | Our N = 32 / 64 [N] | s |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 0.8 | 4.95303 | 400 | 281 (0–560) | 2.83811 | 2.92241 | 21.72 | 1.85e+04 | 1.2e-33 | 3.3e-31 | 1.2e-58 | 1.6e-125 | 0.58241833 | 0.58241833 | 0.586135 / 0.585495 | 34 |
| 0.9 | 6.04965 | 400 | 311 (0–620) | 3.40136 | 2.35917 | 21.77 | 2.43e+04 | 1.5e-33 | 4.8e-31 | 7.9e-61 | 4.9e-130 | 0.27863023 | 0.27863023 | 0.281733 / 0.281431 | 54 |
| 1.0 | 7.38906 | 400 | 341 (0–680) | 4.70966 | 1.05086 | 21.93 | 3.13e+04 | 2e-33 | 6.8e-31 | 4.4e-63 | 1.1e-134 | 0.044037988 | 0.044037988 | 0.0459712 / 0.0456485 | 63 |
| 1.05 | 8.16617 | 400 | 356 (0–710) | 5.05451 | 0.706017 | 21.98 | 3.53e+04 | 2.2e-33 | 8e-31 | 3.1e-64 | 5e-137 | 0.011532089 | 0.011532089 | 0.0123507 / 0.0121974 | 63 |
| 1.09 | 8.84631 | 400 | 368 (0–734) | 5.05451 | 0.706017 | 21.98 | 3.86e+04 | 2.5e-33 | 9e-31 | 3.7e-65 | 6.3e-139 | 0.0037952943 | 0.0037952943 | 0.00404756 / 0.00398826 | 54 |

**DH, even sector** (`data/controls/grid_dh_even.json`)

| a | x = e^{2a} | T♯ | modes (degrees) | A | β̃ | G | M | ε_Q | c_err | ε_B | ε_D | λ₁(mid A) [N] | certified c [R] | Our N = 32 / 64 [N] | s |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 0.8 | 4.95303 | 400 | 281 (0–560) | 2.08108 | 3.67945 | 15.15 | 1.29e+04 | 8.2e-34 | 2.3e-31 | 8.5e-59 | 1.1e-125 | 0.0065246949 | 0.0065246949 | 0.00662227 / 0.00659947 | 46 |
| 0.9 | 6.04965 | 400 | 311 (0–620) | 3.6621 | 2.09842 | 15.31 | 1.71e+04 | 1.1e-33 | 3.4e-31 | 5.5e-61 | 3.4e-130 | 0.00064519287 | 0.00064519287 | 0.000661581 / 0.000660192 | 47 |
| 1.0 | 7.38906 | 400 | 341 (0–680) | 4.07998 | 1.68055 | 15.36 | 2.19e+04 | 1.4e-33 | 4.7e-31 | 3.1e-63 | 7.9e-135 | 3.2965874e-05 | 3.2965874e-05 | 3.53133e-05 / 3.45756e-05 | 57 |
| 1.05 | 8.16617 | 400 | 356 (0–710) | 4.09121 | 1.66931 | 15.36 | 2.47e+04 | 1.6e-33 | 5.6e-31 | 2.2e-64 | 3.5e-137 | 4.5034832e-06 | 4.5034832e-06 | 4.9447e-06 / 4.76475e-06 | 80 |
| 1.1 | 9.02501 | 400 | 371 (0–740) | 5.61514 | 0.14539 | 15.6 | 2.8e+04 | 1.8e-33 | 6.6e-31 | 1.5e-65 | 1.5e-139 | 6.0390737e-07 | 6.0390737e-07 | 7.45326e-07 / 7.08324e-07 | 80 |
| 1.15 | 9.97418 | 400 | 386 (0–770) | 5.61514 | 0.14539 | 15.6 | 3.13e+04 | 2e-33 | 7.7e-31 | 1.1e-66 | 6.1e-142 | 7.5917148e-08 | 7.5917148e-08 | 9.32273e-08 / 8.73322e-08 | 112 |
| 1.19 | 10.8049 | 400 | 398 (0–794) | 5.61514 | 0.14539 | 15.6 | 3.4e+04 | 2.2e-33 | 8.6e-31 | 1.2e-67 | 7.3e-144 | 1.0440918e-08 | 1.0440918e-08 | 1.50903e-08 / 1.31096e-08 | 58 |

**DH, odd sector** (`data/controls/grid_dh_odd.json`)

| a | x = e^{2a} | T♯ | modes (degrees) | A | β̃ | G | M | ε_Q | c_err | ε_B | ε_D | λ₁(mid A) [N] | certified c [R] | Our N = 32 / 64 [N] | s |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 0.8 | 4.95303 | 400 | 281 (1–561) | 2.08108 | 3.67945 | 15.15 | 1.29e+04 | 8.2e-34 | 2.3e-31 | 2.4e-59 | 9e-127 | 0.44593877 | 0.44593877 | 0.467712 / 0.458208 | 54 |
| 0.9 | 6.04965 | 400 | 311 (1–621) | 3.6621 | 2.09842 | 15.31 | 1.71e+04 | 1.1e-33 | 3.4e-31 | 1.6e-61 | 2.8e-131 | 0.09539845 | 0.09539845 | 0.107394 / 0.101351 | 64 |
| 1.0 | 7.38906 | 400 | 341 (1–681) | 4.07998 | 1.68055 | 15.36 | 2.2e+04 | 1.4e-33 | 4.8e-31 | 9e-64 | 6.7e-136 | 0.011772846 | 0.011772846 | 0.0136697 / 0.0129366 | 74 |
| 1.05 | 8.16617 | 400 | 356 (1–711) | 4.09121 | 1.66931 | 15.36 | 2.47e+04 | 1.6e-33 | 5.6e-31 | 6.5e-65 | 3e-138 | 0.0025740039 | 0.0025740039 | 0.00318713 / 0.0029447 | 77 |
| 1.1 | 9.02501 | 400 | 371 (1–741) | 5.61514 | 0.14539 | 15.6 | 2.81e+04 | 1.8e-33 | 6.6e-31 | 4.6e-66 | 1.3e-140 | 0.00037322851 | 0.00037322851 | 0.000470825 / 0.00044726 | 83 |
| 1.15 | 9.97418 | 400 | 386 (1–771) | 5.61514 | 0.14539 | 15.6 | 3.13e+04 | 2e-33 | 7.7e-31 | 3.1e-67 | 5.4e-143 | 4.9711159e-05 | 4.9711159e-05 | 6.83667e-05 / 6.29194e-05 | 95 |
| 1.19 | 10.8049 | 400 | 398 (1–795) | 5.61514 | 0.14539 | 15.6 | 3.41e+04 | 2.2e-33 | 8.6e-31 | 3.6e-68 | 6.5e-145 | 1.0931044e-05 | 1.0931044e-05 | 1.46114e-05 / 1.36908e-05 | 61 |

**Z₁, even sector** (`data/controls/grid_z1_even.json`)

| a | x = e^{2a} | T♯ | modes (degrees) | A | β̃ | G | M | ε_Q | c_err | ε_B | ε_D | λ₁(mid A) [N] | certified c [R] | Our N = 32 / 64 [N] | s |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 0.8 | 4.95303 | 400 | 281 (0–560) | 1.38629 | 9.91161 | 36.29 | 3.09e+04 | 2e-33 | 5.5e-31 | 2e-58 | 2.7e-125 | 1.8497693 | 1.8497693 | 1.8502 / 1.85008 | 104 |
| 1.0 | 7.38906 | 400 | 341 (0–680) | 5.75175 | 5.54616 | 36.71 | 5.24e+04 | 3.3e-33 | 1.1e-30 | 7.4e-63 | 1.9e-134 | 1.0486282 | 1.0486282 | 1.05464 / 1.05319 | 219 |
| 1.1 | 9.02501 | 400 | 371 (0–740) | 10.1462 | 1.15171 | 37.39 | 6.72e+04 | 4.3e-33 | 1.6e-30 | 3.7e-65 | 3.6e-139 | 0.73232047 | 0.73232047 | 0.740358 / 0.740044 | 162 |
| 1.2 | 11.0232 | 400 | 401 (0–800) | 10.1462 | 1.15171 | 37.39 | 8.33e+04 | 5.3e-33 | 2.1e-30 | 1.7e-67 | 5.8e-144 | 0.31010063 | 0.31010063 | 0.320234 / 0.318323 | 88 |
| 1.3 | 13.4637 | 400 | 431 (0–860) | 10.1462 | 1.15171 | 37.39 | 1.02e+05 | 6.5e-33 | 2.8e-30 | 7.1e-70 | 8.3e-149 | 0.11907223 | 0.11907223 | 0.121419 / 0.120904 | 78 |
| 1.31 | 13.7357 | 400 | 434 (0–866) | 10.1462 | 1.15171 | 37.39 | 1.04e+05 | 6.6e-33 | 2.9e-30 | 4.1e-70 | 2.7e-149 | 0.11168787 | 0.11168787 | 0.113794 / 0.113321 | 102 |
| 1.38 | 15.7998 | 1000 | 1076 (0–2150) | 12.9675 | 0.166018 | 41.68 | 3.14e+05 | 5e-32 | 5.4e-29 | 2.5e-123 | 5.4e-257 | 0.04923468 | 0.04923468 | 0.0515496 / 0.0507736 | 1689 |

**Z₁, odd sector** (`data/controls/grid_z1_odd.json`)

| a | x = e^{2a} | T♯ | modes (degrees) | A | β̃ | G | M | ε_Q | c_err | ε_B | ε_D | λ₁(mid A) [N] | certified c [R] | Our N = 32 / 64 [N] | s |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 0.8 | 4.95303 | 400 | 281 (1–561) | 1.38629 | 9.91161 | 36.29 | 3.1e+04 | 2e-33 | 5.5e-31 | 5.8e-59 | 2.2e-126 | 1.4874271 | 1.4874271 | 1.49259 / 1.49029 | 74 |
| 1.0 | 7.38906 | 400 | 341 (1–681) | 5.75175 | 5.54616 | 36.71 | 5.25e+04 | 3.3e-33 | 1.1e-30 | 2.2e-63 | 1.6e-135 | 1.4275906 | 1.4275906 | 1.42949 / 1.42872 | 143 |
| 1.1 | 9.02501 | 400 | 371 (1–741) | 10.1462 | 1.15171 | 37.39 | 6.73e+04 | 4.3e-33 | 1.6e-30 | 1.1e-65 | 3.2e-140 | 1.1589354 | fails | 1.31527 / 1.30789 | 233 |
| 1.1 | 9.02501 | 400 | 371 (1–741) | 10.1462 | 1.15171 | 37.39 | 6.73e+04 | 4.3e-33 | 1.6e-30 | 1.1e-65 | 3.2e-140 | 1.1589354 | 1.1473461 | 1.31527 / 1.30789 | 200 |
| 1.2 | 11.0232 | 400 | 401 (1–801) | 10.1462 | 1.15171 | 37.39 | 8.34e+04 | 5.3e-33 | 2.1e-30 | 5e-68 | 5.2e-145 | 1.1683593 | 1.0515234 | 1.1945 / 1.18877 | 120 |
| 1.3 | 13.4637 | 400 | 431 (1–861) | 10.1462 | 1.15171 | 37.39 | 1.02e+05 | 6.5e-33 | 2.8e-30 | 2.2e-70 | 7.5e-150 | 0.74949749 | 0.74949749 | 0.807909 / 0.787444 | 83 |
| 1.31 | 13.7357 | 400 | 434 (1–867) | 10.1462 | 1.15171 | 37.39 | 1.04e+05 | 6.6e-33 | 2.9e-30 | 1.2e-70 | 2.5e-150 | 0.69226774 | 0.69226774 | 0.743555 / 0.725948 | 73 |
| 1.38 | 15.7998 | 1000 | 1076 (1–2151) | 12.9675 | 0.166018 | 41.68 | 3.14e+05 | 5e-32 | 5.4e-29 | 7.9e-124 | 5.5e-258 | 0.16603101 | 0.16601441 | 0.52632 / 0.521179 | 2747 |

In the first a = 1.1 row the inverse iteration had stalled on the eigenvalue cluster at β̃, so the shift
λ₁(1 − 10⁻⁹) lay above λ_min(A) and the Cholesky factor did not exist. The second row is the rerun with
the step-down of §1.2. For a = 1.1, 1.2 and 1.38, λ₁(mid A) is that stalled estimate, not λ_min.

### 3.3 Where each sector stops, and why

| Function, sector | Largest certified a (x) | T♯ | Certified c [R] | Next support tried | What stops it |
|---|---|---|---|---|---|
| F_{t*} odd | 0.835 (5.31217) | 400 | 1.1666·10⁻³ | 0.84 (5.36556): midpoint block not positive definite at T♯ = 400 and 800 | The crossover. The certified witness at x = 5.36250 < e^{1.68} proves Q < 0 at a = 0.84 [R]. |
| F_{t*} even | 1.09 (8.84631) | 400 | 3.795·10⁻³ | not run | T♯: n = 9 enters at x = 9, raising T₁ from 198 to 478. Not needed, since the odd sector fails first. |
| DH even | 1.19 (10.80490) | 400 | 1.0441·10⁻⁸ | not run | T♯: n = 11 enters at x = 11 and T₁ = 1466 (n = 10 adds nothing, c₁₀ = 0). |
| DH odd | 1.19 (10.80490) | 400 | 1.0931·10⁻⁵ | not run | As DH even. |
| Z₁ even | 1.38 (15.79984) | 1000 | 0.049235 | not run | T♯: n = 14 enters at x = 14 (T₁ = 920, hence T♯ = 1000 for a = 1.38), and n = 16 at x = 16 (T₁ = 1300). At T♯ = 400 the largest support is a = 1.31 (x = 13.7357), c = 0.11169. |
| Z₁ odd | 1.38 (15.79984) | 1000 | 0.16601 (β-limited, β̃ = 0.166018) | not run | As Z₁ even. At T♯ = 400: a = 1.31, c = 0.69227. |

**The gap between R and Q.** Zhu's R charges the frequencies above T♯ only β̃, so λ_min(R) ≤ λ*.
- In the even sectors (F_{t*}, DH, Z₁) the certified c is within 0.5–20% of our N = 64 upper bound at
  every support, so the brackets on the minimum itself are tight. For example, DH even at a = 1.0 has
  3.2966·10⁻⁵ ≤ λ* ≤ 3.4576·10⁻⁵. DH odd is within 21% at a = 1.15.
- In F_{t*}'s odd sector the gap is about 0.002–0.003 in absolute terms and does not shrink with the
  minimum. That is why the last certified support, a = 0.835, is still 0.005 below the crossover in a.
- When λ* > β̃, R is β-limited, and c is just below β̃. This happens for Z₁ odd at a = 1.1 and 1.2, where
  our N = 64 minimum (1.31, 1.19) exceeds β̃ = 1.152; there inverse iteration stalls on the eigenvalue
  cluster at β̃, and μ is stepped down (§1.2).


## 4. Negative side: witness certificates

**Route.** The u-space formula of §1.3, evaluated in Arb (`scripts/control_uspace.py witness`). It has no
frequency cut-off and no tail, so the smoothness class of a witness enters no error term. A witness is the
function whose coefficients are exactly the decimal strings in its JSON file; Arb encloses those decimals,
so each certificate is about that function and does not depend on how the function was found.

| Function | x | Witness (source) | Q(f_w)/‖f_w‖², Arb enclosure [R] | Our Rayleigh quotient [N] |
|---|---|---|---|---|
| F_{t*} | 7 | `data/conductor5/witness_tstar_x7_odd.json`, odd, 32 sine modes | −0.0452795658351271862957293947 ± 4·10⁻⁴² | −0.045279565835127186296 |
| DH | 32 | `data/controls/witness_dh_x32_n80.json`, even, N = 80 minimiser of `build_form(32, 80, "dh")` | −1.39985114183258896882769888·10⁻²⁹ ± 9·10⁻⁷⁰ | −1.3998511418325889688·10⁻²⁹ |
| Z₁ | 20 | `data/controls/witness_z1_x20_n48.json`, even, N = 48 minimiser of `forms("Z1", 20, 48)[1]` | −1.32612579886875655879988·10⁻⁴ ± 5·10⁻⁴⁴ | −1.3261257988687565588·10⁻⁴ |
| F_{t*} | 5.3624996576171875 | `data/controls/witness_ftstar_odd_crossover_n64.json`, odd, N = 64 | −9.8886527992908555289068·10⁻⁷ ± 8·10⁻⁴⁸ | −9.8886527992908555289·10⁻⁷ |
| DH | 30.744970703125 | `data/controls/witness_dh_even_crossover_n80.json`, even, N = 80 | −2.29213630574415050060453·10⁻³² ± 2·10⁻⁷² | −2.2921363057441505006·10⁻³² |
| Z₁ | 19.843902587890625 | `data/controls/witness_z1_even_crossover_n80.json`, even, N = 80 | −1.53420045406999951745·10⁻⁷ ± 4·10⁻⁴⁷ | −1.5342004540699995174·10⁻⁷ |

Outputs are in `data/controls/witness_cert_*.json`. Each run takes under 2 s. The DH margin, 10⁻²⁹ and
10⁻³², is tiny against the O(1) pieces (pole 0, constant term −0.614, integral +0.584, primes −0.029), and
the enclosures are still 40 orders narrower than the margin.

**The three new witnesses** come from bisecting our finite-basis minimum in x (`scripts/control_upper.py
crossover`: F_{t*} odd at N = 64 between 5.1552 and 5.4739; Z₁ at N = 80 between 17.5 and 20; DH at N = 80
between 30 and 30.95). Each is the minimiser at the smallest negative bisection point. They tighten the
negative side from x = 7 to 5.3625 (F_{t*}), from 32 to 30.745 (DH) and from 20 to 19.844 (Z₁).

**Smoothness classes** [N, from the stored coefficients]:
- F_{t*} at x = 7: an odd sine polynomial. f(±L/2) = 0 exactly, but f′(±L/2) = −8.03‖f‖, so f is continuous
  with corners at the window edges, and F̂(t) = O(t⁻²).
- DH at x = 32: an even cosine polynomial with f(±L/2) = −2.53·10⁻¹⁴‖f‖. That is a jump at the edges, so
  F̂(t) = O(t⁻¹).
- Z₁ at x = 20: f(±L/2) = 0.0585‖f‖, also a jump, so F̂(t) = O(t⁻¹).

In both classes ∫ log(2 + t)|F̂|² < ∞, so Q(f_w) is finite and both sides of §1.3 apply. A frequency-side
certificate with a tail bound would have needed a tail of about (4J²/π)(n_κ log(T/2) + O(1))/T for a jump J.
For Z₁ that means T ≈ 10⁴. For DH it also means T ≈ 10⁴, but with the head integral accurate to 10⁻³¹
absolute. The u-space route avoids both.

**From a witness to a smooth test function** [R, written argument].
- Let f ∈ L²[−a, a] with ∫ w|F̂|² < ∞, where w(t) = 1 + log(1 + |t|).
- Put f_δ(u) = f(u/(1 − δ)), supported in [−(1 − δ)a, (1 − δ)a], and f_{δ,ε} = f_δ ∗ η_ε for a mollifier
  of width ε < δa. Then f_{δ,ε} ∈ C_c^∞ is supported in [−a, a].
- F̂_{δ,ε}(t) = (1 − δ)F̂((1 − δ)t) η̂(εt) → F̂ in L²(w dt) as δ, ε → 0. Dilation is strongly continuous on
  L²(w dt), because w(t/(1 − δ)) ≤ 2w(t) for δ ≤ ½; then dominated convergence handles η̂(εt) → 1.
- |Ψ| ≤ Cw and the pole term is continuous on L², so Q(f_{δ,ε}) → Q(f).

Since Q(f_w) < 0, some C_c^∞ function supported in the same window has Q < 0.

## 5. Bracket table

x_pos is the largest certified support with both sectors positive; x_neg is the smallest certified witness.
Both are [R]. The first full-space failure x_c lies in [x_pos, x_neg): x_c ≥ x_pos by positivity, and
x_c < x_neg because Q(f_w) < 0 survives a small dilation of f_w (the continuity argument of §4).

| Function | x_pos [R] | c at x_pos [R] | x_neg [R] | Q(f_w)/‖f_w‖² at x_neg [R] | Brief's witness, certified [R] | Our numerical crossover [N] |
|---|---|---|---|---|---|---|
| F_{t*} | e^{1.67} = 5.31217 (odd limits; even certified to e^{2.18} = 8.84631) | 1.1666·10⁻³ | 5.36249966 | −9.889·10⁻⁷ | x = 7: −0.045280 | (5, 7] at N = 32 (this work); (5.36248, 5.36250] at N = 64 |
| DH | e^{2.38} = 10.80490 (even limits; odd c = 1.0931·10⁻⁵) | 1.0441·10⁻⁸ | 30.74497 | −2.292·10⁻³² | x = 32: −1.3999·10⁻²⁹ | ≤ 30.95, x_c ≈ 30.83 at N = 64 (this work); (30.74404, 30.74497] at N = 80 |
| Z₁ | e^{2.76} = 15.79984 (T♯ = 1000; even limits; odd c = 0.16601, β-limited) | 0.049235 | 19.84390 | −1.534·10⁻⁷ | x = 20: −1.3261·10⁻⁴ | (17.5, 20] at N = 48–80 (this work); (19.84375, 19.84390] at N = 80 |
| ζ (reference) | e^{1.6} = 4.95303 (Zhu, audit, regression above) | 1.0277·10⁻¹⁷ | none | — | — | no failure |

**Reading.**
1. **F_{t*}.** The odd sector's first failure is pinned to 1%: positive for every f supported in
   [−0.835, 0.835], negative at x = 5.3625. The bracket on the minimum itself at a = 0.835 is
   1.1666·10⁻³ ≤ λ*_odd ≤ 3.205·10⁻³ (our sine basis, N = 64). The even sector is far from failure where
   the odd one fails: at a = 0.835 it is above 0.2786, the value certified at a = 0.9. This is what the
   real off-line pair predicts (`docs/CONDUCTOR5_FAMILY.md`): only odd test functions see it.
2. **DH.** The positive side stops at x ≈ 11 because A grows, not because the minimum is near zero. At
   a = 1.15 the certified minimum 7.59·10⁻⁸ is within 13% of our N = 64 value 8.73·10⁻⁸. The comb is
   constant on 10 < x < 11 (c₁₀ = 0), which is why a = 1.19 still certifies at T♯ = 400. Beyond x = 11,
   β* > 0 needs T♯ > T₁ = 2e^{A − K}: 1466 for 11 < x < 12, and 3.8·10⁷ at the numerical crossover
   (A = 17.23 at x = 30.74). That is Zhu's doubly exponential barrier (Theorem 1.4). So for DH the rigorous
   bracket stays wide, and the first failure is pinned well only from above.
3. **Z₁.** Degree 2 halves the exponent (T₁ ≈ 2e^{(A − K)/2}), so the positive side reaches x ≈ 13.7 at
   T♯ = 400 and x ≈ 15.8 at T♯ = 1000 (28 and 46 minutes for the two sectors). The negative side is at 19.844.
4. **Our numerical crossovers** all lie inside the rigorous brackets, as they must, since a finite-basis
   minimum is an upper bound on the full one.

## 6. What is not rigorous, and the exact missing pieces

| Item | Status | Exact missing piece |
|---|---|---|
| Q is Weil's functional, Q(f) = Σ_ρ ĝ(γ_ρ) | [P] | The explicit formula for these gamma factors, for F without an Euler product (F_{t*}, DH, Z₁). It needs only the functional equation, finite order and a half-plane where −F′/F converges absolutely; it is standard but was not re-derived here. This work checks it numerically against census zeros (`docs/CONDUCTOR5_FAMILY.md`, `docs/DAVENPORT_HEILBRONN.md`, `docs/EPSTEIN.md`) [N]. Every certificate here is a statement about Q as defined in §1.1. |
| Starting enclosure of the Bessel ratios | [R, given a theorem] | Pincherle's theorem (a minimal solution has a convergent continued fraction) is cited, not proved. `scripts/control_bessel_check.py` compares the recurrence enclosures with Arb's own `bessel_j` at 20 arguments x ∈ [10⁻³, 600], n up to 1.5x + 80: 168 comparisons, all overlapping (`data/controls/bessel_check.json`) [R]. |
| Arb and python-flint 0.9.0 | trusted | Outward rounding, and the correctness of `acb.digamma`, `acb.integral`, `arb.legendre_p_root`, `arb.bessel_i` and `arb_mat` products. Not audited. gmpy2 is used only for the floating Cholesky factor, whose residual is verified in Arb, so it needs no trust. |
| The new scripts | unreviewed | `control_cert_lib.py`, `control_certificate.py` and `control_uspace.py` reproduce the audited ζ numbers exactly, agree with our matrices to 10⁻³⁹ and with an independent frequency-side Arb integral, but have not been read line by line by a second person. |
| L² test functions | wording | As in the audit (item 4), Q is defined on L² by the frequency integral, with values in (−∞, +∞]. Theorem 1.1 is pointwise in frequency and the u-space identity is Tonelli, so both hold on L². §4 passes from the witnesses to C_c^∞. |
| DH positive side beyond x = 11 | not attempted | The next comb step needs T♯ > 1466 (n = 11 enters at x = 11), i.e. about 1300 Legendre modes and 190 000 nodes; estimated well over 30 minutes at 2 processes. |
| Z₁ positive side beyond x = 16 | not attempted | The T♯ = 1000 odd run at a = 1.38 took 46 minutes, over the 30-minute guideline, because another session loaded the machine; it was already running when that became clear, and it was the last extension. n = 16 enters at x = 16, and T₁ = 1300. A run at T♯ ≈ 1500 needs about 1700 Legendre modes and 190 000 nodes. The T♯ = 1000 runs already took 28 and 46 minutes on the loaded machine, so this would exceed the 30-minute limit several times over. |
| F_{t*} even beyond x = 9 | not attempted | n = 9 raises A to 5.94, so T♯ > 478 is needed. The odd sector fails first anyway, so this does not move the bracket. |

## 7. Reproduce

```bash
PY=<repo>/.venv/bin/python   # mpmath 1.4.1 (gmpy2), python-flint 0.9.0; 2 processes
# regression (ζ, a = 0.8)
$PY scripts/control_certificate.py --function zeta --sector even --a 0.8 --tsharp 200 --nmodes 200 --workers 2 --mu 9e-18,auto --blocks 20,40,60,80,100,120,160 --json data/controls/cert_zeta_even_a0.8.json
$PY scripts/control_certificate.py --function zeta --sector odd  --a 0.8 --tsharp 150 --nmodes 200 --workers 2 --mu 8.2065e-15,auto --blocks 20,40,60,80,100 --json data/controls/cert_zeta_odd_a0.8.json
# normalisation
$PY scripts/control_normalisation.py --cases zeta:even:e^1.6,zeta:odd:e^1.6,ftstar:even:7,ftstar:odd:7,ftstar:even:e^1.7,ftstar:odd:e^1.7,dh:even:32,dh:even:e^2,dh:odd:e^2,z1:even:20,z1:even:e^2.4,z1:odd:e^2.4 --tfreq 2000 --json data/controls/normalisation.json
# positive side: one grid per function and sector (rows accumulate in the summary file)
$PY scripts/control_cert_grid.py --function ftstar --sector odd  --a 0.6,0.7,0.8,0.81,0.82,0.83,0.835,0.84 --tsharps 400,800 --stop --workers 2 --summary data/controls/grid_ftstar_odd.json
$PY scripts/control_cert_grid.py --function ftstar --sector even --a 0.8,0.9,1.0,1.05,1.09 --tsharps 400 --stop --workers 2 --summary data/controls/grid_ftstar_even.json
$PY scripts/control_cert_grid.py --function dh --sector even --a 0.8,0.9,1.0,1.05,1.1,1.15,1.19 --tsharps 400 --stop --workers 2 --summary data/controls/grid_dh_even.json
$PY scripts/control_cert_grid.py --function dh --sector odd  --a 0.8,0.9,1.0,1.05,1.1,1.15,1.19 --tsharps 400 --stop --workers 2 --summary data/controls/grid_dh_odd.json
$PY scripts/control_cert_grid.py --function z1 --sector even --a 0.8,1.0,1.1,1.2,1.3,1.31 --tsharps 400 --stop --workers 2 --summary data/controls/grid_z1_even.json
$PY scripts/control_cert_grid.py --function z1 --sector odd  --a 0.8,1.0,1.1,1.2,1.3,1.31 --tsharps 400 --stop --workers 2 --summary data/controls/grid_z1_odd.json
$PY scripts/control_cert_grid.py --function z1 --sector even --a 1.38 --tsharps 1000 --workers 2 --max-minutes 45 --summary data/controls/grid_z1_even.json
$PY scripts/control_cert_grid.py --function z1 --sector odd  --a 1.38 --tsharps 1000 --workers 2 --max-minutes 45 --summary data/controls/grid_z1_odd.json
# witnesses: find with our matrices, certify in u-space
$PY scripts/control_witnesses.py --function dh --x 32 --n 80 --dps 60 --json data/controls/witness_dh_x32_n80.json
$PY scripts/control_witnesses.py --function z1 --x 20 --n 48 --dps 40 --json data/controls/witness_z1_x20_n48.json
$PY scripts/control_upper.py crossover --function ftstar --sector odd --lo 5.1551695 --hi 5.4739474 --n 64 --dps 30 --steps 14 --json data/controls/crossover_ftstar_odd_n64.json --witness data/controls/witness_ftstar_odd_crossover_n64.json
$PY scripts/control_upper.py crossover --function z1 --sector even --lo 17.5 --hi 20 --n 80 --dps 40 --steps 14 --json data/controls/crossover_z1_even_n80.json --witness data/controls/witness_z1_even_crossover_n80.json
$PY scripts/control_upper.py crossover --function dh --sector even --lo 30 --hi 30.95 --n 80 --dps 60 --steps 10 --json data/controls/crossover_dh_even_n80.json --witness data/controls/witness_dh_even_crossover_n80.json
$PY scripts/control_uspace.py witness --function ftstar --json data/conductor5/witness_tstar_x7_odd.json --parity odd --x 7 --prec 256 --out data/controls/witness_cert_ftstar_x7_odd.json
$PY scripts/control_uspace.py witness --function dh --json data/controls/witness_dh_x32_n80.json --parity even --x 32 --prec 320 --out data/controls/witness_cert_dh_x32.json
$PY scripts/control_uspace.py witness --function z1 --json data/controls/witness_z1_x20_n48.json --parity even --x 20 --prec 256 --out data/controls/witness_cert_z1_x20.json
$PY scripts/control_uspace.py witness --function ftstar --json data/controls/witness_ftstar_odd_crossover_n64.json --parity odd --x 5.3624996576171875 --prec 256 --out data/controls/witness_cert_ftstar_crossover.json
$PY scripts/control_uspace.py witness --function z1 --json data/controls/witness_z1_even_crossover_n80.json --parity even --x 19.843902587890625 --prec 256 --out data/controls/witness_cert_z1_crossover.json
$PY scripts/control_uspace.py witness --function dh --json data/controls/witness_dh_even_crossover_n80.json --parity even --x 30.744970703125 --prec 384 --out data/controls/witness_cert_dh_crossover.json
# Our upper bounds at the grid supports, entry and Bessel checks, tables
$PY scripts/control_upper.py points --function ftstar --sector odd --a 0.6,0.7,0.8,0.81,0.82,0.83,0.835,0.84 --n 32,64 --dps 30 --json data/controls/upper_ftstar_odd.json   # likewise for the other five sectors
$PY scripts/control_entry_check.py data/controls/certs/*.json
$PY scripts/control_bessel_check.py --json data/controls/bessel_check.json
python3 scripts/control_report_tables.py
```

Logs go to `logs/` through `scripts/progress.py`.

## 8. Numbers that differ from our existing results or the audited ζ values

No computed value disagrees. Every comparison either matches or is a refinement in the direction the
variational principle requires (a larger basis gives a lower finite minimum and an earlier crossover).

| Quantity | Existing value (source) | This work | Status |
|---|---|---|---|
| ζ even certified bound, μ = 1.02768955902·10⁻¹⁷ | 1.0276895590199816·10⁻¹⁷ (`zhu_certificate_even.json`) | 1.0276895590199925·10⁻¹⁷ | differs from the 14th digit: c_err 7.46·10⁻³² here against 1.84·10⁻³¹, a different bound on the same quadrature (§1.2) |
| ζ ε_Q, c_err | 9.2151·10⁻³⁴, 1.843·10⁻³¹ even; 6.7244·10⁻³⁴, 1.3449·10⁻³¹ odd (audit) | 3.7285·10⁻³⁴, 7.457·10⁻³²; 2.7207·10⁻³⁴, 5.4415·10⁻³² | method (exact-rule bound instead of stored-rule defects) |
| F_{t*} odd crossover | finite-basis transition in (5, 7] at N = 32; full-space failure at x ≤ 7, "no lower bound" (`docs/CONDUCTOR5_FAMILY.md`) | N = 32: (5.366, 5.474); N = 64: (5.36248, 5.36250]; full-space failure certified in [5.31217, 5.36250) | refinement |
| DH crossover | N = 64 fine scan x_c ≈ 30.83, failure at x ≤ 30.95; "for DH one does, at x ≈ 30.5" (`docs/CONNES_LETTER.md`) | N = 80: (30.74404, 30.74497]; failure certified at x ≤ 30.74497 | refinement; the doc's "x ≈ 30.5" is looser than its own fine scan |
| DH minimum at x = 32 | −1.07727·10⁻²⁹ at N = 64 (`data/connes/dh_crossover.txt`) | −1.3998511418·10⁻²⁹ at N = 80, certified | consistent (larger N) |
| Z₁ crossover | (17.5, 20] at N = 48–80 (`docs/EPSTEIN.md`) | N = 80: (19.84375, 19.84390]; failure certified at x ≤ 19.84390 | refinement |
| Z₁ minimum at x = 20 | −1.3261·10⁻⁴ at N = 48; −1.70·10⁻⁴ at N = 80 (`docs/EPSTEIN.md`) | −1.32612579887·10⁻⁴ certified (N = 48); −1.6962·10⁻⁴ at N = 80 (bisection run) | match |
| F_{t*} witness, x = 7 | −0.045279565835127186296 (`witness_tstar_x7_odd.json`) | −0.04527956583512718629573 ± 4·10⁻⁴² | match |
| Brief's ζ targets | 1.02768955902·10⁻¹⁷, 9.11833845·10⁻¹⁵ | 1.02768955902·10⁻¹⁷, 9.11833844538·10⁻¹⁵ (λ₁ = 9.11833845450·10⁻¹⁵) | match |

## Five-line summary

1. **Reduction:** Zhu's Theorem 1.1 carries over unchanged to Ψ = Σ_κ Re ψ(κ + it/2) + K − Σ (2c_n/√n) cos(t log n), with A = Σ 2|c_n|/√n and Binet envelopes re-derived for κ = ¾ and for degree 2. The Arb code reproduces the audited ζ certificates exactly (1.02768955902·10⁻¹⁷ even, 9.11833845·10⁻¹⁵ odd) [R].
2. **Normalisation:** Our zeros-side matrices for F_{t*} (both sectors), DH and Z₁ are this Q with factor 1, to 10⁻³⁹ against a rigorous u-space evaluation, and inside independent frequency-side Arb enclosures. There is no mismatch [R/N].
3. **Positive:** Q ≥ c‖f‖² for every f supported in [−a, a]: F_{t*} at x = 5.31217 (c = 1.1666·10⁻³), DH at x = 10.80490 (c = 1.0441·10⁻⁸), Z₁ at x = 15.79984 (c = 0.049235, T♯ = 1000; x = 13.7357 at T♯ = 400) [R].
4. **Negative:** all three given witnesses are certified negative by an exact u-space formula in Arb (−0.0452796 at x = 7, −1.39985·10⁻²⁹ at x = 32, −1.32613·10⁻⁴ at x = 20), and so are new witnesses at x = 5.36250, 30.74497 and 19.84390 [R].
5. **Brackets:** first full-space failure in [5.312, 5.3625) for F_{t*} (1% wide), [10.80, 30.745) for DH and [15.80, 19.844) for Z₁. DH and Z₁ are stopped on the positive side by T♯ ≈ 2e^{(A−K)/n_κ} (DH would need T♯ ≈ 4·10⁷ at its crossover). The remaining trust points are the explicit formula for non-Euler-product F [P], Pincherle's theorem [P] and Arb.
