# Independent audit of rh2's support-2.38 certificate

Date: 2026-10-01. Claim commit: `6a1bccd`. Audit base: `2d3310d`, whose
certificate code and data are unchanged from that claim. Worktree:
`<review-worktree>`; branch: `audit/certificate-238`.
The main checkout and other worktrees were not modified. No packages were installed.
Python: `<repo>/.venv/bin/python`, python-flint 0.9.0,
mpmath 1.4.1, gmpy2. All FLINT and numerical jobs use one thread; at most two
CPU-intensive jobs run concurrently.

Evidence labels used throughout:

- **[R]** *re-derived*: written mathematical argument or exact arithmetic check.
- **[I]** *recomputed independently*: new code, independent of the authors' evaluation routines;
  numerical checks are explicitly distinguished from rigorous bounds.
- **[A]** *reproduced with the authors' code*: reruns or direct checks of that implementation.
- **[S]** *taken from a source*; paper statements additionally marked *taken from the paper*.
- **[U]** *could not verify*.

## Verdict

**Verified with corrections.** The derivation is valid with
corrections. Two implementation errors were found: a float-based Schur maximum can understate
the upper bound, and the integer residual helper omits rounding before its integer flooring.
The first is closed by an independently constructed, strictly smaller rigorous Schur bound.
The second is paid conservatively without moving the overall constant. The fresh even-sector
rerun, with a sharper quadrature norm estimate, proves that same decimal constant.
Both 450-mode assemblies agree to about 11 significant digits [R, I, A].

The reported even-sector data, with outward allowances and an extra residual debit of
`1e-104`, still give

\[
 Q(f)\ge 6.8131164695121494552\ldots\,10^{-48}\|f\|_2^2
       >6.81311646951\,10^{-48}\|f\|_2^2.
\]

The fresh 680-mode rerun also proves the original overall constant after the sharper
quadrature norm estimate below: its corrected lower bound is
`6.8131164695951201359…e-48`. The unmodified one-worker script's raw lower bound is
`6.8131164695065578584…e-48`, slightly below the brief's last decimal because of a larger
radius after serial accumulation. This is an enclosure difference, not an eigenvalue
disagreement [R, I, A]. The odd-sector decimal in the brief,
`3.91966861610e-44`, was rounded upward from its recorded lower bound. Use
`3.91966861609e-44`. This correction does not change the minimum over sectors [R].

For arbitrary L² functions, interpret Q as the extended form defined by its frequency
integral, permitting `+infinity`, or restrict the finite-valued statement to its logarithmic
form domain. Not every compactly supported L² function has finite Q [R].

## L1–L10 derivation table

| Step | Assessment | Reason and effect on the constant | Evidence |
|---|---|---|---|
| L1 | Verified with a correction | Explicit-formula coefficients and normalization factor 1 are correct. The half-line formula is for real f; complex f requires the full-line integral before L10. Q on all L² must be understood as an extended form. No change to the overall constant. | [R, S] |
| L2 | Verified | Compact convolution widens support by δ. Plancherel gives the compressed shift form with weights Λ(n)/√n. The operator retains the original seven shifts when acting on the wider window. | [R] |
| L3 | Verified | The lower envelope is increasing for t>0. Its exact dyadic β′ is below the bound at 700. The independent Schur bound validates the μ actually used despite the authors' Schur-selection bug. | [R, I] |
| L4 | Verified | Transform, total mass, sidelobe interval estimate, both passband tails, and integration-by-parts remainder have the stated constants. | [R, I] |
| L5 | Verified | Legendre and pole bounds, parity tail sums and Schur bounds yield the displayed finite-to-infinite inequality. The real-axis bound on m follows from the nonnegative main lobe and two signed sidelobe tails. | [R] |
| L6 | Verified | Use the analytic symmetrization of digamma, not the literal real-part operation at complex z. Every ellipse stays in the pole-free strip. The Chebyshev truncation argument gives the stated Gauss error. | [R] |
| L7 | Verified | The ratio interval is valid above the turning point; the upward recurrence is enclosed below it. Exact midpoint evaluation and derivative padding recover the node ball. All ten 160-digit reference values were contained. | [R, I, A] |
| L8 | Verified with corrections | “Nonnegative” means positivity preserving, not positive semidefinite. The Schur row maximum is incorrectly selected through binary64. Independent weight: μ≤3.466159084208 on a larger window, validating the original μ and the brief's target. | [R, I, A] |
| L9 | Verified with a correction | The abstract Gram-residual argument is valid. The implementation's +N flooring units omit preceding 400-bit rounding. An extra 1e-104 covers that error without changing the even decimal constant. | [R, I, A] |
| L10 | Verified | The real symmetric form splits over real/imaginary parts and reflection parity. Both sector bounds therefore cover complex f. | [R] |

## Detailed derivation

### L1, L2 and L10: normalization, support, and complex functions

Use F(t)=∫f(u)e^{itu}du and the unnormalized Lebesgue L² norm. For real f,
Plancherel is `||f||²=(1/pi)∫_0^infinity |F|²`. The comb acts as

\[
 C=\sum_{n\in\{2,3,4,5,7,8,9\}}
       \frac{\Lambda(n)}{\sqrt n}(\tau_{\log n}+\tau_{-\log n}),
\]

with zero extension and compression to the indicated interval. Its multiplier is exactly
P(t)=Σ 2Λ(n)n^{-1/2}cos(t log n). Since log 9<2.38<log 11 and Λ(10)=0,
these are precisely the contributing prime powers [R]. The explicit-formula convention
is Zhu v2 §2; no rescaling of f or the norm occurs, so the conversion factor is 1
[S, R; [Zhu v2](https://arxiv.org/html/2608.24827v2)].

The pole contribution for real f is
`2(∫f cosh(u/2))²−2(∫f sinh(u/2))²`.
For complex f replace the two squares by absolute squares and use the full-line multiplier
integral `(1/(2pi))∫_R Ψ|F|²`. The multiplier and pole kernels are real and symmetric.
Writing f=p+iq therefore gives Q(f)=Q(p)+Q(q); equivalently, the mixed Fourier term is
odd in t and integrates to zero. Reflection commutes with the form, giving the even/odd
split. The positive-half-line identity alone is not valid for arbitrary complex f [R].

W and κ are real and even, hence K and m are real and even on the real frequency axis.
κ is integrable and supported in [−δ,δ], so h=f−f*κ is L², supported in [−a−δ,a+δ],
and has transform mF. For real f,

\[
 \langle h,Ch\rangle=\frac1\pi\int_0^\infty P(t)m(t)^2|F(t)|^2dt
 \le\bar\mu\,\frac1\pi\int_0^\infty m(t)^2|F(t)|^2dt.
\]

Adding this inequality to Q gives precisely Φ=H−P+m²(P−μ). C here contains the original
comb, not all prime powers below exp(2b). No pointwise bound P≤μ is asserted or needed [R].

The archimedean symbol grows logarithmically at infinity, while its negative part is bounded.
Consequently its frequency integral is defined for every supported L² function with value
in `(−infinity,+infinity]`. This supplies the L² extension; an assertion that Q is finite
on all of L² would require correction [R].

### L3 and L4: filter, envelope and rigorous tails

The Kaiser transform follows, for example, by integrating the power series of
`I0(β sqrt(1−u²/δ²))`; it is entire when written as the 0F1 expression. Its real main lobe is

\[
 \widehat W(s)=\frac{2\delta}{I_0(\beta)}
 \frac{\sinh\sqrt{\beta^2-\delta^2s^2}}{\sqrt{\beta^2-\delta^2s^2}}.
\]

Fourier inversion gives `(1/(2pi))∫ W_hat=1`, with the improper oscillatory tails understood
as limits. Multiplication by sin(Ωu)/(πu) convolves with the frequency interval indicator,
so `K(t)=(1/(2pi))∫_{t−Ω}^{t+Ω} W_hat(s) ds` [R].

In a sidelobe, set r=sqrt(δ²s²−β²). The normalized integral becomes

\[
 \frac1{\pi I_0(\beta)}\int \frac{\sin r}{\sqrt{r^2+\beta^2}}dr.
\]

The amplitude is positive and decreasing. The second mean value theorem, or integration
against the primitive of sine, bounds every interval by `2/(pi β I0(β))=ε`.
This also bounds an infinite tail by taking a limit. For t≤Ω−β/δ, both omitted tails
are sidelobes, hence |m|≤2ε. For t≥Ω+β/δ, K is one sidelobe interval, so |K|≤ε [R].

In the band write c=(1−main_total)/2 and
`m=c+∫_{−β/δ}^{t−Ω} W_hat/(2pi)+right_tail`.
The main-lobe integrand is nonnegative, `|c|≤ε`, and the right tail has absolute value ≤ε.
Thus throughout the band `−2ε≤m≤1+2ε`, a bound that is also valid in the passband.
This justifies the real-axis estimate used in L5 [R].

For the right tail, g(r)=(r²+β²)^−1/2 has the Laplace representation
`∫_0^infinity J0(βs)e^{−rs}ds`. Differentiating K times and using |J0|≤1 yields
`|g^(K)(r)|≤K!/r^(K+1)`. Repeated integration by parts gives
`e^{ir0} Σ_{k=0}^{K−1} i^{k+1}g^(k)(r0)` plus a remainder whose absolute value is at most
`∫_{r0}^infinity K!/r^(K+1)dr=(K−1)!/r0^K`.
The code's Taylor-coefficient recurrence follows by equating coefficients in
`p f'=−(p'/2)f`, p=p0+p1x+x². Its factorial and i-power indexing is correct [R].

Now `w=1−m²=2K−K²`, so |w|≤2ε+ε² beyond 700 and
`Φ=H−μ−w(P−μ)`. With |P−μ|≤A+μ this gives the claimed envelope.
For z=1/4+it/2, Binet's remainder is bounded by 1/(3t), since
`|s²+z²|≥2 Re(z) Im(z)=t/4` and `∫s/(e^{2pi s}−1)ds=1/24`.
Also Re(1/(2z))≤1/(2t²) and log|z|≥log(t/2).
For t≥3/4 these losses sum to at most 1/t. The derivative of
`log(t/(2pi))−1/t−μ−constant` is `1/t+1/t²>0` [R].

The exact dyadic selected by the implementation is

\[
 \beta'=\frac{1432587375811093321}{2^{60}}
       =1.242571476103756130705824\ldots.
\]

The floor-minus-one construction and a certain Arb comparison put it below the envelope
at 700. The Schur μ used here is legitimate because the independent rigorous upper bound
in §Recomputations is strictly smaller [R, I, A].

### L5: finite reduction and coupling

For the orthonormal Legendre basis the transform is i^n S_n with
`S_n(t)=sqrt(2L(2n+1)) j_n(Lt)`. Rodrigues' formula and integration by parts give

\[
 j_n(x)=\frac{x^n}{2^{n+1}n!}\int_{-1}^1e^{ixs}(1-s^2)^n ds,
 \qquad |j_n(x)|\le\frac{|x|^n}{(2n+1)!!}\quad(x\in\mathbb R).
\]

For the pole vector, the positive series for the modified spherical Bessel function gives
`i_n(z)≤z^n/(2n+1)!! exp(z²/(2(2n+3)))`. These yield exactly d_n and e_n in lines 224–229.
The ratio d_{n+2}/d_n equals
`sqrt((2n+5)/(2n+1)) y²/((2n+3)(2n+5))`; it decreases in n.
For e the extra exponential ratio is ≤1. Geometric summation with that ratio is therefore
valid. The first omitted degrees are 1360 even and 1361 odd; stepping by two is correct [R].

After the diagonal sign change removing i^n, the pole matrix is `+2pp^T` in the even sector
and `−2pp^T` in the odd sector. For a retained degree m and omitted degree n,

\[
 |B_{mn}|\le 2|p_m|e_n+(G_{\rm real}T/\pi)\sqrt{2L(2m+1)}d_n.
\]

The code's B∞ and B1 bound row and column sums; B1 uses the larger sum over omitted
indices where their maximum would suffice. For the tail block the first omitted e_n,d_n
bound their maxima, producing ε_D. Symmetry then gives
`||B||≤sqrt(B∞B1)` and `||D−β'I||≤ε_D` [R].

The block-diagonal comparison has lower bound min(λ_min(A),β′−ε_D).
The off-diagonal matrix with B,B* has operator norm ||B||, so subtracting ε_B gives the
stated full-form lower bound. The pole term remains intact [R].

For G_real, the digamma series gives `Re ψ(1/4+it/2)≥ψ(1/4)`.
Binet gives the upper bound `log|z|+2+4/3` on the entire real interval.
Together with |P|≤A and |m|≤1+2ε this proves line 416's G_real, approximately 26.908323.
The two arguments of the float-selected H_abs maximum are widely separated, so that
selection is unambiguous here [R, A].

### L6: analyticity and quadrature

The holomorphic extension of H is
`(ψ(1/4+iz/2)+ψ(1/4−iz/2))/2−log pi`. The nearest poles are at z=±i/2;
one must not extend H by taking a real part of a complex argument.
The panel ellipses have |Im z|≤1/4. Both digamma arguments have real part ≥1/8;
m, the comb and S_n are entire. The true integrand is analytic on and inside every ellipse [R].

Writing ψ(w)=ψ(1+w)−1/w, Re v≥9/8 for v=1+w. Binet gives remainder at most
`1/(12(Re v)²)`: indeed `|s²+v²|≥(Re v)²`.
Since |v|>1 and |arg v|<π/2,
`|log v|≤log|v|+π/2`, yielding the code's psi_abs bound [R].

On the real u interval 0≤W≤1. Splitting the sinc integral at |u|=1/Ω gives
`∫|κ|≤(2/pi)(1+log(Ωδ))` when Ωδ≥1.
Hence on the strip `|m(z)|≤1+exp(δb_e)(2/pi)(1+log(Ωδ))`.
The Legendre integral representation gives `|S_n(z)|≤sqrt(2L(2n+1)) exp(Lb_e)`.
These establish G_strip and Mb in lines 403–408 [R].

With ρ=2+sqrt(5), the ellipse's imaginary semiaxis is
`h(ρ−1/ρ)/4=h=1/4`. Truncate the Chebyshev series of each mapped integrand at degree
2n−1. Its sup error is at most `2Mρ^(−2n)/(1−1/ρ)`.
Gauss n-point quadrature integrates that polynomial exactly. Positivity of the weights,
whose sum is h, bounds the difference of integral and quadrature by twice h times this
sup error. Summing the panels gives `ε_Q=4TMρ^(−2n)/(1−1/ρ)`.
A symmetric matrix with entry error ≤ε_Q has norm ≤Nε_Q [R].

For the fresh rerun the audit uses a sharper version of this same estimate.
Let u_i=sqrt(2L(2n_i+1)). Before replacing every u_i by its maximum, the entry errors obey
`|E_ij|≤C u_i u_j`, with C independent of i,j. Cauchy–Schwarz gives
`||E||₂≤C Σ_i u_i²`; it does not require E to be positive semidefinite.
Thus the even-sector quadrature debit can be multiplied by
`Σ(2n_i+1)/(N(2nmax+1))=(2N−1)/(4N−3)=1359/2717`.
For odd degrees the factor is `(2N+1)/(4N−1)=1361/2719`.
This reduces the even debit from about 1.77190e-58 to at most 8.86277e-59.
That saved allowance is sufficient to preserve the original decimal constant in the
one-worker rerun, after also paying the residual-rounding error [R, I].

### L7: Bessel enclosures

For n>x, the real ratio r_n=j_n/j_{n−1} lies in (0,1).
It is the recessive continued-fraction solution
`r_n=1/((2n+1)/x−r_{n+1})`; the denominator is >1 when r_{n+1}∈[0,1].
Starting at a sufficiently high n with [0,1] therefore encloses the true tail ratio;
downward interval iteration preserves containment. The code correctly stops this recurrence
above the turning point and multiplies the enclosed ratios onto the upward value [R].

The upward recurrence is an exact identity evaluated with Arb, so its accuracy depends on
precision but containment does not. The 400+1.25x+96 precision rule is an accuracy heuristic,
not an unenclosed truncation. Evaluating at the exact dyadic midpoint avoids amplification
of the node radius [R].

From `j_n(x)=((-i)^n/2)∫_{−1}^1 e^{ixs}P_n(s)ds`, differentiation and |P_n|≤1 give
`|j_n'(x)|≤(1/2)∫_{−1}^1|s|ds=1/2` on the real axis. Thus adding half the input radius
is sufficient. This representation is [DLMF 10.54.2](https://dlmf.nist.gov/10.54.E2)
[S, R]. The ten 160-digit checks below include both recurrence regimes and high degrees [I, A].

### L8 and L9: Schur and exact Gram arithmetic

C is self-adjoint and positivity preserving. It is not positive semidefinite: take two
small disjoint bumps separated by log 9, with opposite signs and supports small enough that
none of the other shifts meets them. Their quadratic form is negative. The weighted Schur
test needs nonnegative translation weights, not nonnegative spectrum [R].

For a positive step weight φ, the row ratio Cφ/φ bounds ||C|| by weighted Cauchy–Schwarz
(or weighted Young). The certified strict bounds q<log n/h<q+1 ensure a shifted cell meets
only i+q,i+q+1, and the negative shift meets i−q−1,i−q. Omitting off-window cells is correct
up to measure-zero endpoints. Float tuning is harmless because any positive exact-dyadic
weight is admissible. The final maximum, however, must be selected without dropping larger
Arb upper endpoints; the authors' implementation does not do that [R, A].

For any trial factor V, `A_mid−μI=VV^T+R`, so λ_min(A_mid)≥μ−||R||₂.
R is symmetric, hence ||R||₂≤||R||∞. The interval entry radii and quadrature errors cost
N max_rad and Nε_Q. A Cholesky proposal is only a convenient way of obtaining V;
its floating computations need not themselves be rigorous [R].

Floors of exact scaled matrix entries cost <one integer unit per entry. Lines 201–203
first round the subtraction to a 400-bit midpoint, however, and floor that rounded value.
The +N units pay only the final floor. A precise counterexample and a conservative repair
are given next [R, A].

## Code findings with line numbers

Line numbers refer to the audited base, not the new audit scripts [A].

1. **Actual Schur maximum understated** — `scripts/grid_norm.py:157–159`.
   The routine compares `float(ratio.upper())`, discarding distinctions between larger
   endpoints mapping to the same float. Recomputing every row with the same author's
   weight gives

   - returned: `3.46920322110173162911152813316…`;
   - actual maximum upper endpoint: `3.46920322110173176899711048500…`;
   - understatement: `1.3988558235183424…e-16`.

   This invalidates that helper's general rigorous-upper-bound claim [A, R]. The audit's
   independent bound `3.46615908420798001…` is smaller even on a larger window, and therefore
   validates the μ used throughout the present certificate [I]. Repair: compare Arb upper
   endpoints directly, or compute an outward max, without binary64 conversion.

2. **Missing rounding debit in “exact residual”** —
   `scripts/grid_certificate_rigorous.py:201–211`. The integer multiplication is exact,
   but the scaled input has already been rounded by `.mid()` after subtraction.
   At precision 400 set `ell=1−2^−400`, A=[ell], μ=2^−402, trial factor [ell], P=600.
   The helper reports approximately `1.49969681390e-241`, while the exact residual is
   approximately `2.90444393614e-121` [R, A]. Both ell and μ are exact dyadics.
   Repair: convert A_mid and a dyadic lower shift separately to exact scaled integers,
   then subtract the integers. Alternatively pay the pre-floor midpoint error.

   For the actual blocks, |A_mid_ij|<2^30 from the independent analytic entry bound
   `G_real T/pi * 2L(2nmax+1)+4L(2nmax+1)e^L+β′`.
   A per-entry allowance 2^−360 more than covers 400-bit subtraction rounding and the
   very small radius of the decimal μ_s. Its row debit at N=680 is less than
   `2.8955e-106`. Replace the recorded residual by the larger `1e-104` allowance.
   The independently recomputed final lower bounds remain
   `6.8131164695121494552…e-48` even and `3.9196686160954133410…e-44` odd [R, I].
   This correction does not move the claimed overall decimal constant.

3. **Odd decimal rounded upward** — `docs/CERTIFICATE_238.md`, results table, and the brief.
   `3.91966861610e-44` exceeds the stored certified bound
   `3.9196686160954133411…e-44`. Print `3.91966861609e-44` as a lower bound.
   The overall even-sector constant is already rounded down [R].

4. **No radius-loss finding at max_rad** — `grid_certificate_rigorous.py:461–466`.
   Initially suspicious float conversion is exact here: Arb radius magnitudes have at most
   31 significand bits, and the maximum radii are ordinary non-subnormal binary64 values.
   Thus conversion to a 53-bit binary64 and back preserves the radius exactly [R, A].
   FLINT's [magnitude type](https://github.com/flintlib/flint/blob/main/src/mag.h#L109)
   uses `MAG_BITS=30`; the actual JSON radius mantissas are also checked [S, I].
   This observation does not justify float conversion of the 400-bit Schur endpoints.

5. **Other reviewed sites passed at these parameters** [R, A]:
   comb construction (48–58, 310, 315) and wider-window shifts use precisely the seven
   original prime powers; outward Schur window (306–309) is asserted with exact Fractions;
   certain Arb band comparisons (356–361) cover all nodes; band ordering through floats is
   safe because distinct Gauss nodes are well separated at these magnitudes; cumulative
   integration endpoints telescope, so slack need not be paid per intermediate segment.
   At 400 bits their first and last endpoint errors, including both main-lobe endpoints,
   are far below 1e-100. The global transform peak bounds all endpoint slivers (351, 363–380).
   Φ (126) and pole signs (385–395, 464) are correct. β′ (328–330) is dyadic and checked
   with certain Arb comparisons. Parity degrees and n0 (216, 385–386) are correct.
   The temporary floating Cholesky (171–186) and inverse iteration (475) are proposals,
   not proof steps; their correctness is checked through the residual argument.
   The floating selections in the fixed-parameter coupling maxima pick separated values.

## Independent recomputations and reproduction

### Rigorous independent μ

`audit238_components.py` starts from the constant weight and uses a new lazy nonlinear
iteration, without calling the authors' Schur-weight or Galerkin routines for this bound. It uses
8000 cells, 2200 iterations, and the *larger exact-decimal window* b=1.490000002.
All shift indices are certified with 400-bit Arb. The resulting positive binary64
weights are treated as exact dyadics, and every row ratio is evaluated in Arb;
the maximum is selected by comparing exact outward upper endpoints [I].

\[
 \|C\|\le3.4661590842079800092756837108161355993\ldots
            <3.46920322110.
\]

Monotonicity under enlarging the window makes this valid for the actual 1.49 window and
the authors' outward 1.490000001 window. The data retain all 8000 weights and certified
shift indices so the final inequality can be checked without trusting the iteration [I, R].

### Direct u-space m and Bessel containment

The filter reference is computed by direct integration of κ(u)cos(tu) over 24 subintervals
of [0,δ], with mpmath at 160 decimal digits. It does not use W_hat, the band decomposition,
or the authors' integration-by-parts expansion. The comparison enclosure is evaluated
with the authors' main-lobe and tail methods at 400 bits. Ten midpoint frequencies partition
the band [273.333…,700]; every 160-digit reference ball is contained [I, A].

| t (rounded for display) | Direct u-space m(t), about 12 digits | Contained in the 400-bit filter ball? |
|---|---|---|
| 294.6666667 | 1.14792548794e-17 | Yes |
| 337.3333333 | 7.89430987520e-10 | Yes |
| 380 | 1.79284902707e-5 | Yes |
| 422.6666667 | 0.00769029035395 | Yes |
| 465.3333333 | 0.212021409573 | Yes |
| 508 | 0.787978590427 | Yes |
| 550.6666667 | 0.992309709646 | Yes |
| 593.3333333 | 0.999982071510 | Yes |
| 636 | 0.999999999211 | Yes |
| 678.6666667 | 1.00000000000 | Yes |

Bessel references use the half-integer ordinary Bessel function in mpmath at 160 digits.
The x reference is the exact dyadic midpoint used by the recurrence, reconstructed from
its integer mantissa and exponent; the comparison includes the original node-radius padding.
References are for S_n, including the sqrt normalization, not just j_n [I, A].

| n | t | Contained? |
|---|---|---|
| 0 | 0.01 | Yes |
| 1 | 0.1 | Yes |
| 10 | 5 | Yes |
| 50 | 40 | Yes |
| 100 | 90 | Yes |
| 100 | 100 | Yes |
| 101 | 100 | Yes |
| 400 | 500 | Yes |
| 850 | 700 | Yes |
| 898 | 700 | Yes |

The full 155-digit values, 145-digit printed enclosures and containment outcomes are
in `data/connes/audit238_components.json` [I, A]. Neither these finite samples nor their
agreement replace the analytic containment arguments [R].

### Independently assembled 450-mode even block

Completed in 9m50s [I]. `audit238_assembly.py` imports neither certificate nor prototype.
It uses 360-bit MPFR, its own Newton-generated Gauss nodes and weights, its own downward
Miller Bessel recurrence normalized against elementary j0/j1, direct u-space Kaiser quadrature
at every frequency node, and 1400 panels of width 1/2 with 90-point Gauss.
The resulting matrix arithmetic and digamma use FLINT. This shares the underlying arithmetic
library, but not the authors' Bessel, filter, quadrature or assembly implementation [I].

The u-space rule uses 256 Gauss nodes. Doubling it to 512 at ten band frequencies changes
m by less than 3.7e-106. This is a numerical convergence check, not an interval certificate.
An independent numerical LDL inertia check and new inverse iteration assess the matrix.
The μ scalar and exact dyadic β′ are held equal to the certificate's inputs for comparison [I].

The independently assembled midpoint gives
`lambda_1=6.81313469650521178076295538072464…e-48`.
All 450 numerical LDL pivots are positive, and the inverse iteration stabilizes by its
fifth step. Against the published 680-mode midpoint `6.81313469666455e-48`, the relative
difference is about `2.34e-11`, or agreement to roughly 11 significant digits [I, A].
This difference is smaller than the enclosure debits; interval midpoints at different
truncations are not themselves exact matrices for applying eigenvalue interlacing [R].
The authors' 450-mode run gives `6.81313469666455e-48`, identical to its 680-mode
midpoint at every printed digit. The relative difference between the independent and
author 450-mode assemblies is about `2.34e-11`, using the printed author midpoint. This completes the requested
same-truncation assembly comparison [I, A].

### Author reruns

The 680-mode even sector completed unchanged with `--workers 1` in 1658 seconds (27m37s
in the progress log), preserving all original parameters. It reproduces both midpoint
eigenvalues to every printed digit: λ₁=6.81313469666455e-48 and
λ₂=1.23364078837902e-40. Cholesky succeeds at the same shift; its recorded integer residual
is about 2.933e-118 [A].

The max entry radius is `1.6785067801613226e-56`, versus the original four-worker
`1.6785059578451558e-56`. Their difference costs about 5.59175e-60 after multiplying
by 680, enough to move the raw lower bound's final printed digit. Using the rigorously
sharpened quadrature bound and the extra residual payment, this fresh output gives
`6.8131164695951201359…e-48 > 6.81311646951e-48` [R, I, A].

The 450-mode author run completed in 1041 seconds (17m20s in the progress log) [A].
It additionally records a new exact-integer residual
computed by separate conversion of A_mid and the dyadic lower shift. The author's calculations
and returned residual remain unchanged for reproducibility [A, I]. Its independently
computed exact residual is `2.950974911912226…e-118`, compared with the helper's
`2.963797570563459…e-118`. The helper happens to overestimate for this matrix; the
counterexample proves that this is not guaranteed in general [R, I, A].

A 450-mode run tests the finite block only: its elementary complement bounds are too large
to certify the full form at that truncation. The published whole-form certificate uses
680 modes [R]. No fresh odd-sector rerun was requested; its recorded components and conventions
are audited, and the final debit calculation is repeated independently [A, R, I].

## Prior art

These findings are **taken from the paper [S]**, not an independent audit of those authors'
certificates. The primary HTML/PDF texts were read directly; no AI overview was used as evidence.

- **Liu**, *Certified Weil Positivity Beyond the Unit Window*, Theorem B and Appendix B.4:
  **yes**, explicitly uses compressed prime translations and a weighted Schur estimate,
  with a quadratic positive weight on the support interval. Thus a support-window prime
  bound is prior art. Theorems A/B certify full widths 2 and 17/8=2.125, not ≥2.38.
  Its construction differs from this Kaiser frequency split; these readings do not establish
  priority for every variant. [Paper](https://www.alphaxiv.org/abs/2609.weil-positivity-riemann-zeta-bounds.pdf)
- **Suzuki**, arXiv:2606.09096, §2.3 equations (2.5), §2.4, and Theorems 1.1/1.4:
  explicitly writes the compressed finite-window prime translations. I found no weighted
  compressed-prime norm substituted for a torus supremum in a quantitative certificate.
  Its positivity theorem concerns sufficiently small windows; it does not certify width 2.38.
  [Paper](https://arxiv.org/html/2606.09096)
- **arXiv:2607.24830 is by Kim, Hong, Kim, Choi, Jang and Kim, not Suzuki**.
  It numerically realizes Suzuki's operator using finite elements. It does not present the
  requested rigorous window-prime envelope or a support-2.38 positivity certificate.
  Its §4 explains its numerical and conditional limitations.
  [Paper](https://arxiv.org/html/2607.24830)

The broad claim that compressed prime translations or weighted support-window Schur bounds
are new would therefore be incorrect [S].

## Reproduce and retained evidence

All paths below are relative to this audit worktree. Logs use `scripts/progress.py`.
Run at most two one-thread jobs concurrently. The listed long jobs were backgrounded [A].

```sh
PY=<repo>/.venv/bin/python
$PY scripts/audit238_components.py
$PY scripts/audit238_bookkeeping.py
$PY scripts/audit238_assembly.py --nmodes 450 --gl 90 --bits 360
$PY scripts/grid_certificate_rigorous.py --a 1.19 --delta 0.3 --beta-k 64 --tsharp 700 --sector even --nmodes 680 --gl 56 --workers 1 --json data/connes/audit238_author_even680.json
$PY scripts/audit238_author450.py --a 1.19 --delta 0.3 --beta-k 64 --tsharp 700 --sector even --nmodes 450 --gl 56 --workers 1 --json data/connes/audit238_author_even450.json
$PY scripts/audit238_verify_debits.py
$PY scripts/audit238_compare.py
```

The new scripts and resulting JSON data are committed with this report. The original certificate
scripts and original JSON files are preserved. Trust assumptions remain: FLINT/Arb correct
outward arithmetic and special functions, the explicitly re-derived analytic inequalities,
and correctness of this audit's new scripts. Numerical reference and assembly agreement are
cross-checks, not full replacements for the ball certificate [R, I, A].

`data/connes/audit238_comparison.json` retains the acceptance checks, relative differences,
run times, source SHA-256 hashes and limitations. All comparison assertions and all six
audit-script syntax checks passed. The two author jobs and independent component/assembly
jobs exited successfully [I, A].

## Five-line summary

1. Derivation: L1–L10 survive with the stated complex-form, domain and Schur terminology corrections [R].
2. Independent bounds: μ≤3.466159084208; ten filter and ten 160-digit Bessel values are contained [I, A].
3. Computation: the two 450-mode assemblies agree to about 11 digits; the 680-mode rerun matches both published eigenvalues; Schur selection and residual-rounding defects are closed by independent bounds and paid errors [R, I, A].
4. Prior art: Liu already uses a compressed-prime Schur bound and certifies width 2.125; neither compared work certifies 2.38 [S].
5. Verdict: verified with corrections. Fresh even output plus the sharper quadrature norm proves the original overall constant 6.81311646951e-48; the odd printed lower bound must be rounded down [R, I, A].

Commit: recorded in the final delivery message and the audit branch history (a commit cannot
contain its own hash).

## Three most fragile places, ranked

1. Schur certification: exact positive weight values are harmless, but choosing the final
   maximum through floats loses rigor. This audit found an actual understatement [A, R].
2. Gram residual conversion: exact integer multiplication does not undo earlier midpoint
   rounding. Every input-to-integer conversion must have its own error payment [A, R].
3. Kaiser transition and tiny eigenvalue: cancellation near m=0 and the approximately
   1e-48 spectral margin demand both enclosed filter values and adequate quadrature;
   sample agreement alone would be insufficient [R, I].
