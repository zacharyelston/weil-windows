# The Quarter Law: a Window Sees a Quarter of the Torus

**Statement.** Let b > 0 and X = e^{2b}. Define the window-compressed prime-shift operator on L²[−b, b] (functions extended by zero):

  (C_b f)(u) = Σ_{2 ≤ n < X} (Λ(n)/√n) [f(u + log n) + f(u − log n)],  u ∈ [−b, b].

This is the multiplier P(t) = Σ (2Λ(n)/√n) cos(t log n) of Weil's prime comb, seen through the window (`docs/GRID_NORM.md`). Its supremum over t, attained at the torus origin by Kronecker, is A(b) = 2Σ_{n<X} Λ(n)/√n. Then, unconditionally,

  **‖C_b‖ = e^b (1 + o(1)),  A(b) = 4 e^b (1 + o(1)),  hence ‖C_b‖/A(b) → 1/4.**

In terms of x = X: a window-limited test function sees about √x of the comb, while the torus supremum is about 4√x. **Consequence:** the window-norm certificate has threshold about 2π e^{√x}, against 2π e^{4√x} for Zhu's envelope. The exponent drops by a factor of 4 (by e^δ/4 when the window is widened by δ). The cost is still doubly exponential in the support a, since x = e^{2a}.

## Proof

**Notation.** Write ν = Σ_n (Λ(n)/√n) δ_{log n} for the weight measure on (0, 2b), and ν_c = e^{ℓ/2} dℓ for its prime-number-theorem model. By Abel summation from ψ(t) = t + o(t),

  ν((0, ℓ]) = ∫_1^{e^ℓ} t^{−1/2} dψ(t) = 2e^{ℓ/2} + o(e^{ℓ/2}),

so A(b) = 2ν((0, 2b)) = 4e^b(1 + o(1)).

C_b is self-adjoint and positivity-preserving: it has a nonnegative kernel, though it is not positive semidefinite. Hence ‖C_b‖ = sup over f ≥ 0 of ⟨f, C_b f⟩/‖f‖², since |⟨f, Cf⟩| ≤ ⟨|f|, C|f|⟩.

**Lower bound (two ends).** Fix M > 0. Put f = f_L + f_R, with

  f_L(−b + x) = e^{−x/2} and f_R(b − y) = e^{−y/2} for 0 ≤ x, y ≤ M (zero elsewhere).

Then ‖f‖² = 2(1 − e^{−M}). Every term of ⟨f, C_b f⟩ is nonnegative, so it is at least the left–right cross term. Only a shift ℓ = 2b − x − y joins the two pieces, and

  2⟨f_L, C_b f_R⟩ = 2e^{−b} ∫ e^{ℓ/2} L(2b − ℓ) dν(ℓ),

where L(s) is the length of [max(0, s − M), min(M, s)], continuous and bounded by M. The prime number theorem lets us replace dν by ν_c on [2b − 2M, 2b] with error o(e^{2b}), for fixed M. With ν_c the right-hand side is

  2e^{−b} ∫ e^{ℓ} L(2b − ℓ) dℓ = 2e^{b} ∫_0^{2M} e^{−s} L(s) ds = 2e^b (1 − O(Me^{−M})).

Taking M → ∞ slowly gives ‖C_b‖ ≥ e^b(1 − o(1)). By Cauchy–Schwarz, the profile e^{−x/2} is optimal for the two-end coupling: (∫ g e^{−x/2})²/∫g² ≤ ∫e^{−x} = 1.

**Upper bound (Schur test with weight cosh(u/2)).** For φ > 0, ‖C_b‖ ≤ sup_u (C_b φ)(u)/φ(u). Take φ(u) = cosh(u/2). With the model ν_c, the change of variables v = u ± ℓ turns (C_b φ)(u) into the continuous kernel:

  ∫_{−b}^{b} e^{|u−v|/2} cosh(v/2) dv = e^b cosh(u/2) + ½[(u + b)e^{u/2} + (b − u)e^{−u/2}] − cosh(u/2).

So the continuous ratio is e^b + O(b) uniformly in u, and in fact at most e^b + 4b. For ν itself, split the integral at a large fixed ℓ₀ and integrate by parts against ν − ν_c. Each integrand φ(u ± ℓ)·1(cut-off) is piecewise monotone in ℓ, with total variation O(e^{b}cosh(u/2)) relative to its weight. The PNT bound |(ν − ν_c)((0, ℓ])| ≤ ε e^{ℓ/2} + C(ℓ₀) for ℓ ≥ ℓ₀ then gives an error of at most (ε + o(1))·O(e^b cosh(u/2)). Hence ‖C_b‖ ≤ e^b(1 + O(ε) + o(1)), and ε → 0. ∎

**Rigour.** This is a written argument at the level of a careful sketch. The integration-by-parts step in the upper bound is standard but has not been written out with explicit constants. Under RH, the error ψ(t) − t = O(√t log² t) would give ‖C_b‖/A(b) = 1/4 + O(b² e^{−b/2}).

## Numerical confirmation

Galerkin lower bounds (K = 800–1000 cells; `scripts/grid_norm.py`) for the actual primes, and the exact top eigenvalue μ_cont of the continuous kernel e^{|u−v|/2}. μ_cont comes from the even eigenfunction cosh(κu), κ² = ¼ + 1/λ, together with the boundary condition at u = b.

| b | μ(b) | A(b) | μ/A | μ/e^b | A/(4e^b) | μ_cont/e^b |
|---|---|---|---|---|---|---|
| 2.0 | 7.974 | 24.38 | 0.3270 | 1.079 | 0.825 | 1.232 |
| 3.0 | 22.50 | 75.19 | 0.2992 | 1.120 | 0.936 | 1.178 |
| 4.0 | 59.08 | 213.5 | 0.2767 | 1.082 | 0.978 | 1.103 |
| 4.5 | 95.55 | 354.6 | 0.2695 | 1.062 | 0.985 | 1.074 |
| 5.0 | 155.0 | 588.0 | 0.2636 | 1.044 | 0.991 | 1.052 |

The prime ratio μ/e^b tracks the continuous kernel's 1 + O(b e^{−b}) correction from b ≈ 4 on, and A/(4e^b) → 1. So μ/A → 1/4 from above, slowly. This is why the earlier tables (`docs/GRID_NORM.md`) showed ratios between 1/3 and 1/4 at moderate b.

## Reading

The geometry is the one guessed in `docs/GRID_NORM.md`. A function confined to the window does best by sitting at both ends, with e^{−x/2} tails, linked by the longest prime shifts, which carry the heaviest weight e^{ℓ/2}.
- **On the prime torus,** the comb reaches 4√x only where all phases align, at the Kronecker returns.
- **Through the window,** compact support forbids that coherence and leaves √x.

Quantitatively, A = 2ν((0, 2b)) ≈ 4e^b counts every shift at full cosine amplitude, whereas the best window function attains about ν((0, 2b))/2 ≈ e^b.

The law holds for ζ and for every control alike, since it uses only the density of prime powers. That fits the null results on Euler separation (`docs/SEMILOCAL_COUPLING.md`, the window-norm kill criterion).
