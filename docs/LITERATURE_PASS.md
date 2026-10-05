# Literature pass: Weil forms, prolates, and the DH control

Date: 2026-10-01. Base: `7f3956e` on `feat/weil-gram-instrument`.
Worktree: `<review-worktree>`; branch: `litpass/weil-prolate`.
This is a reading and analytic review. No numerical experiment, including the
optional quick check, was run. Repository values below are reported evidence,
not independently re-certified here. No changes were made in the main checkout.

The main finding is that the E-map reconstruction is already documented,
principally as numerical evidence. A leakage-weighted cross-term estimate was
not found in the inspected sources. There are also material qualifications to
the crossover, pole terminology, and coupling-split verdict. The assertion
about all degree-1 conductor-5 families is too broad: an explicit even family
exists, although the odd, real, periodic DH equation is rigid.

## Notation and evidence standard

Write `x = λ²`. Repository test functions have additive support
`[-(log x)/2, (log x)/2]`; their autocorrelation has support
`[-log x, log x]`. In papers using `[-a,a]`, therefore `a = (log x)/2`.
This distinction matters for all support comparisons.

References [A]–[O] are specified at the end. **Proved** means a theorem with
a proof in the inspected source, not that this pass has independently audited
every step. **Numerical** includes high precision without interval certification.
For the recent preprints [H] and [I], the distinction between a theorem as
presented by the author and an independently verified certificate is explicit.
“Not found” is confined to the sources and sections inspected here.

## Classification of C1–C8

| Claim | Status | Source and evidence | Consequence for rh2 |
|---|---|---|---|
| C1: form, table reproduction, entry 48 | **Partly known**; literal unconditional modulus-square wording needs correction | [B], (1.1), (2.8), Proposition 2.1: **proved** explicit-formula identity. [A], §5, p.25: **numerical** table, described as upper bounds. Reproduction and the smaller entry 48 are rh2 numerics. | The zeros form is defined on the full test-function space, but becomes a sum of absolute squares only under RH. “Pole-free” names the zeros side, not omission of the geometric pole correction. No published correction to entry 48 was found. |
| C2: own gamma/conductor leakage and factors 1.5–15 | **Partly known** | [A], §6.4, Figure 1: **numerical** comparison for ζ. [B], §3: near-radical construction. The DH conductor rescaling and quantitative factors were **not found** in [A]–[D], [G]–[I]. | These finite-range ratios remain numerical observations, not universal constants or an asymptotic comparison theorem. |
| C3: whole low spectrum in a Fourier class | **Partly known** | [B], §3, (3.4), Definition 3.1 and Figures 26–36: **numerical** reconstruction of several even and odd eigenfunctions by orthogonalized E-images. [C], §3.6, (29): **proved** construction of the Fourier-class vectors. | Multi-vector reconstruction is already documented. Neither the stated DH sequence nor a theorem identifying the entire eigenvalue sequence with leakages was found there. |
| C4: k, angles, Rayleigh ratios, zero-error equality | **Partly known** | [A], §6.4, (17): construction; §6.5, Fact 6.4: stated convergence result. [B], §3: **numerical** reconstruction. The DH angles and the three-way error comparison were **not found** in those sources. | A convergent candidate is known; its quantitative identification with the ground state is not established. The proposed equalities of error scales are empirical and root-dependent. |
| C5: DH positivity up to 30 and continuum crossover | **Partly known**, with unproved upgrades | [E], Theorem 6.1: **proved, conditional** real-zero theorem; it does not require positivity of the minimum. [J], pp.747–748 and [K], §§1–2: **numerical** DH zeros. rh2 supplies the support scan. | Finite-dimensional positivity cannot certify the full-space lower endpoint 30. Simplicity and evenness must be checked before applying the continuum real-zero theorem to DH. |
| C6: leakage-relative margins decay comparably | **Not found** | No corresponding metric comparison in inspected [A]–[D], [F]–[I]. Evidence is the numerical one-shot and review in this repository. | Raw near-radical eigenvalues and generalized margins are different quantities. No inference that the relative margin tends to zero follows just from the raw floor tending to zero. |
| C7: positive diagonal blocks, negative coupling Schur complement | **Not found** as this specific numerical claim | rh2 `SCHUR_BLOCKS.md`: **numerical**. Schur positivity equivalence is an exact algebraic fact, not an RH-specific result of the cited papers. | The diagnosis is valid for the chosen finite decomposition if its reported blocks are accurate. It does not isolate an invariant arithmetic cause; see the lower-bound correction below. |
| C8: grid attribution is ill-posed | **Contradicted as stated**; conditioning claim supported numerically | rh2 explicitly verifies the additive split. [B], (2.11)–(2.12), defines individual local contributions. No impossibility theorem for stable attribution was found. | A specified finite split is well-defined. Recovering its small residual is ill-conditioned, and interpreting an artificial mask causally needs a genuine functional-equation control. |

## Q1. Positive-operator inequalities with primes

[D], Theorem 1, gives, for smooth multiplicative `g` supported in
`[2^(-1/2),2^(1/2)]` with `ĝ(i/2)=ĝ(0)=0`,

\[
 W_\infty(g*g^*)\geq\operatorname{Tr}(\vartheta(g)S\vartheta(g)^*).
\]

Here `S` projects onto even Sonin functions vanishing, together with their
Fourier transforms, on `[-1,1]`. Theorem 6.11 permits `ĝ(-i/2)=0` alone at
the price of `-c|ĝ(0)|²`, where `c=4γ/log 2`, `γ≈2.94355` from Lemma 6.10;
this γ is **not** Euler's constant. Remark 6.12 brackets the optimal c by 13
and 17. These are proved archimedean estimates. The prime terms vanish in
the interior of this autocorrelation support; they do not handle x=13.

[C], Theorem 4.6 in §4.7, proves a Hilbert-space isomorphism between classical
and semilocal Sonin spaces for every finite S containing infinity and every
λ>0. An isomorphism alone gives no ordering between QW and a positive trace.

A relevant recent result is [I], v2, Theorems 1.1–1.2 and Corollary 6.3.
It claims the unconditional bound

\[
 Q(f)\geq8.9\,10^{-18}\|f\|_2^2,
 \qquad\operatorname{supp}f\subset[-0.8,0.8],
\]

for complex f, with prime powers 2,3,4 included. Its x is `exp(1.6)≈4.953`.
Theorem 1.1 reduces a frequency envelope to a finite matrix with bounded
tail and coupling errors. This is an unrefereed certificate claim, not a
Sonin domination theorem; normalization was checked, certificate arithmetic
was not independently executed. Its support-2.38 claim is retracted in
§7 and Remark 3.3. No corresponding prime-inclusive Sonin
inequality was found in inspected [A]–[D], [G], [H].

**Takeaway:** there is a concrete prime-inclusive finite-window certificate
claim to audit; a general semilocal Sonin lower bound remains unlocated.

## Q2. What is proved about the smallest eigenvalue?

[F], Theorem 3, proves attainment on a finite union of compact multiplicative
intervals. Its Theorem 12, §12, gives a small additive-support estimate: for
support interval length `b<log 2`,

\[
 Q(F)\geq
 \left(\log(1/b)-\log\log(1/b)-O(1)\right)\|F\|_2^2.
\]

This yields positivity for sufficiently small b; the O(1) notation does not
itself certify every b up to log 2. Yoshida's Lemmas 2–3 and Theorem 2 are
discussed explicitly in [H], §1.1; the original [O1] was inaccessible here.
Their numbering and scope are consequently attributed to that account,
not presented as an independent reading of Yoshida.

[G], Theorem 3.6 and Corollary 3.7, establish compact resolvent, lower bounded
discrete spectrum, and a ground state for each finite λ. [H], Theorems 1.1–1.4,
provides a Friedrichs realization through a continuous screw kernel,
continuity of the minimum in a, and positive, simple, even ground state for
sufficiently small a. Its §4.5 discusses the gap in the older continuity
argument. These are theorems as proved in the September 2026 preprint; they
are not a DH crossover certificate.

The elementary variational upper bound is always

\[
 \epsilon(\lambda)\leq
 \frac{QW_\lambda(k_\lambda)}{\|k_\lambda\|^2}.
\]

To replace its right-hand side by `C(1-χ₂)` one still needs a quantitative
estimate for the projected E-image, including its norm. Such a theorem was
not found in [A]–[D], [G], [H]. The numerical comparison in [A], Figure 1,
does not supply it. Fact 6.4 states uniform substrip convergence of k̂λ
with bound `c λ^(-1/2-α)/(1-2α)` on `Im z=α`, `|α|<1/2`;
this is not the repository's exponentially small square-root-leakage bound.

For [I], distinguish the finite-window certificate just described, its
unconditional variational upper bounds (§11), its RH-conditional asymptotic
upper bound (Theorem 1.3), and its empirical decay law (Conjecture 12.1).
None is a proved uniform comparison with `1-χ₂`.

Finally [B], Corollary 2.4, proves that the finite trigonometric minima
converge downward to the full-space infimum. Therefore a negative finite
Rayleigh quotient proves a negative continuum direction, whereas a positive
finite minimum is an upper bound and leaves the continuum sign unresolved.

**Takeaway:** existence, small-window positivity, continuity, and variational
upper bounds are available; the desired uniform leakage comparison was not found.

## Q3. What can control the weighted cross terms?

For an orthonormal comparison span G with leakage weights ℓj and complement
U, the relevant object is an operator bound, not independent row estimates:

\[
 B_{ju}=\frac{QW(g_j,u)}{\sqrt{\ell_j}},\qquad
 S=A-BC^{-1}B^*,\qquad C>0.
\]

One needs `BC⁻¹B*≤A` or a useful quantitative sufficient condition obtained
without presupposing positivity of QW. Summing pointwise bounds over zeros
is circular if their positivity uses RH; rowwise bounds also need control of
the number and mutual correlations of the near-radical vectors.

**Trace formula.** [A], §7.4, (22), is

\[
 -\sum_{v\in S}W_v(f)=\log(TW)f(1)+
 \operatorname{Tr}\!\left(\vartheta(f)(1-P_T^S-\widehat P_W^S)\right).
\]

Polarizing with `f=g_j*u*` gives a trace expression for a cross term. It does
not bound it: `1-P-Q` can be negative on an intersection or near intersection
of projection ranges. To use the identity one must control the compressed
operator, trace-class domains, and the cancellation with the scalar term.
No ℓj-normalized estimate is supplied by that identity.

**Variational and continuous-kernel methods.** [F], Lemma 1 and §§4–7,
supplies Euler–Lagrange/resolvent tools. [H], (2.9)–(2.10) and §8, replaces
the distribution kernel by

\[
 QW(v,w)=\int_{-a}^{a}\!\int_{-a}^{a}
 g(x-y)v'(y)\overline{w'(x)}\,dy\,dx
\]

on compactly supported smooth functions. This is an analytic way to retain
local cancellations before taking norms. It still needs boundary/domain
control to cover the repository's projected E-images and all complement
vectors. No square-root-leakage residual estimate was found in these sections.

**Toeplitz and Carathéodory–Fejér.** [D], §§6.1–6.7, (134)–(140), bounds a
compact error operator by a finite approximation and controls a rank-one
repair with Lemma 6.9. This is an actual norm-estimate technique. [E],
Corollary 1.1, Theorems 5.6 and 6.1, localizes zeros after subtracting the
lowest eigenvalue. That subtraction makes the form positive semidefinite
even if the original minimum is negative. Thus the real-zero theorem cannot
supply the missing sign or cross bound.

**Hankel estimates.** [N], Theorem 3.7, characterizes bounded positive Hankel
forms on the upper-half-plane Hardy space by a positive measure μ with
`dρ=dμ/(1+t²)`, `ρ((0,r))=O(r)` and `ρ((r⁻¹,∞))=O(r)` as r→0.
Theorem 4.1 constructs a bounded symbol. These results require identifying
a positive Hankel form first. No identification of the rh2 B block, or its
ℓj-weighted norm, with such a form was found; replacing that missing step
with assumed positivity would not advance the argument.

An exact finite-dimensional diagnostic clarifies what a useful estimate
should target. If C>0, put `X=C⁻¹B*`. Completing the square gives

\[
 Q(y,z)=\langle Sy,y\rangle+
 \langle C(z+Xy),z+Xy\rangle.
\]

The two terms are stable mathematical objects for a specified decomposition.
They do not give an invariant division into “Euler failure” and “other failure.”
The user-requested independent numerical reimplementation has been set aside
for this literature brief; this identity does not re-certify the saved blocks.

**Takeaway:** compact-kernel residual and finite-error norm methods give
specific routes to try; no unconditional leakage-weighted cross bound was found.

## Q4. Controls without Euler products

DH's failure of RH is established independently of the repository. [J]
computes four off-line zeros below height 200. [K], §2, computes additional
ones through deformation of periodic Dirichlet series. Its Theorem 1 proves
local persistence of zeros under a continuous coefficient deformation,
not their complete enumeration or a support threshold for Weil negativity.
Neither paper computes the rh2 window crossover.

[O2], Bombieri–Ghosh, is a directly relevant 50-page survey, especially its
sections on zeros and the coefficients of `1/f`. Its full text could not be
retrieved successfully in this pass (MathNet returned blocks/non-PDF pages).
Consequently **its coverage of support-dependent Weil positivity remains an
unclosed reading gap**, rather than evidence of absence. The same limitation
applies to original [O3] Potter–Titchmarsh and [O4] Bateman–Grosswald.
Their bibliography was verified. Stark's indexed introduction, p.47 [O5],
attributes infinitely many critical-line Epstein zeros to Potter–Titchmarsh
and a real off-line zero for a sufficiently elongated binary form to
Bateman–Grosswald. Do not attribute all Epstein off-line-zero results to
Potter–Titchmarsh on that basis.

[M], Introduction, (1.1)–(1.2), supplies a reachable primary account of
Epstein analytic continuation and the dual-lattice functional equation.
Arithmetic Epstein functions associated to a *single class* must be
distinguished from a Dedekind zeta function. A whole-field Dedekind zeta
retains its Euler product regardless of class number. Class sums are linear
combinations of Hecke L-functions; cancellation or genus factorizations
must be checked for the particular combination. Class number by itself is
not a controlled switch with every other datum fixed.

For any proposed no-Euler form, the explicit formula needs all zeros and
poles of the chosen completion and correct logarithmic-derivative
coefficients. An off-line zero entails eventual failure of the full Weil
positivity criterion under its admissibility hypotheses; it does not specify
the first compact support where failure is visible. Positive truncations do
not give such a threshold. Real off-line zeros additionally expose a parity
issue: an even-only test space need not detect them (see Q5).

**Takeaway:** zero-based controls are well established, but no published rh2-style
support threshold was found in the accessible sources; Bombieri–Ghosh remains unread.

## Q5. Genuine functional-equation-preserving families

### The conductor-5 statement needs a parity qualification

[K], Theorem 2, §4, gives dimensions of periodic solution spaces.
At period 5 the positive even equation has dimension 2, while the positive
odd DH equation has dimension 1. Thus the normalized real periodic DH
equation admits no continuous family in that space. “No continuous degree-1,
conductor-5 family” without that qualification is false.

An explicit normalized real family is

\[
 F_t(s)=(1-t)L(s,\chi_5)+t(1+\sqrt5\,5^{-s})\zeta(s),
 \qquad 0\leq t\leq1,
\]

where χ5 is the quadratic character modulo 5. [K], §3, (6)–(7), supplies
the common even functional equation. Equivalently

\[
 \Lambda_t(s)=(5/\pi)^{s/2}\Gamma(s/2)F_t(s),
 \qquad\Lambda_t(s)=\Lambda_t(1-s).
\]

Here is an elementary check of the switch, independent of numerics.
At both endpoints the Dirichlet coefficients are multiplicative: the second
endpoint has ζ's Euler product with a modified factor at 5. For `0<t<1`,
`a_t(2)=a_t(3)=2t-1` and `a_t(6)=1`. Hence
`a_t(6)≠a_t(2)a_t(3)`, so the normalized Dirichlet series is not an Euler
product. The first coefficient is 1 throughout.

Two caveats are essential. First, this is the **even** gamma factor, not DH's
odd one. Second, pole data varies: the L endpoint is entire, and the other
endpoint has a pole at 1; its residue is multiplied by t. This family keeps
the functional equation, conductor, and root number fixed, but not the pole
divisor and residue. Its modified Euler factor also need not satisfy the
strict Selberg local-coefficient axiom.

[L1], Theorem 2, classifies degree 1 in the extended Selberg class as
Dirichlet-polynomial combinations of shifted primitive character L-functions.
Its Theorem 3 imposes the Euler axiom and reduces the ordinary Selberg class
to individual primitive L-functions (or ζ). [L2]'s theorem on p.2 separates
periodicity without the Euler axiom from the character conclusion with it.
The classification does not exclude no-Euler combinations.

### Epstein shapes preserve more geometric data, but change the series grid

For determinant-one positive binary Q, define

\[
 Z_Q(s)=\sum_{(m,n)\ne(0,0)}Q(m,n)^{-s},\qquad
 \Lambda_Q(s)=\pi^{-s}\Gamma(s)Z_Q(s).
\]

The dual-lattice equation [M], (1.2), together with
`Q⁻¹=JᵀQJ`, `J=[[0,-1],[1,0]]`, yields
`Λ_Q(s)=Λ_Q(1-s)`. Thus `Q_r=diag(e^r,e^-r)` gives a continuous
self-dual family with fixed gamma factor and fixed residue π at s=1.
At r=0 the square-lattice identity is `Z_Q=4ζ(s)L(s,χ₋₄)`.
These last binary identities follow directly from change of lattice variables
and the sum-of-two-squares representation formula.

This is a generalized Dirichlet series on the real norm grid, not a family
of integer-index Dirichlet series of fixed arithmetic conductor. Moving
between arithmetic discriminants changes the usual conductor normalization;
one cannot feed these functions unchanged into the current integer-grid
explicit-formula instrument. If literal fixed-pole and fixed-gamma control is
required, this geometric family is preferable, but “Euler product off” must
be formulated and verified at the particular arithmetic or generalized-series
points being compared.

**Takeaway:** the odd DH equation is rigid in the stated periodic real space;
an even conductor-5 switch exists, and determinant-one Epstein shapes offer
a different continuous control with a different grid.

## Q6. Fact checks

### Exact real-zero hypotheses

[E], Theorem 6.1, requires a real distribution D on `[0,L]`, the form defined
on trigonometric polynomials by its equation (6), a lower-bounded
**essentially selfadjoint** operator, and a simple isolated spectral minimum
with an even eigenfunction. The trigonometric polynomials are an operator
core in the proof. [A]'s Theorem 6.1 states the corresponding convolution
kernel formulation. Neither requires the minimum to be nonnegative or an
Euler product. Neither proves its own simplicity/evenness hypotheses for
every DH window. Finite approximants have a separate theorem [E], 5.6;
finite even-sector simplicity alone does not certify the full continuum
ground state. [G], §8, explicitly leaves the ζ simplicity/evenness problem
among the remaining steps.

### DH values and the precision actually published

The following is a **published numerical list**, not an argument-principle
certificate of completeness. Reflection gives the left-half-strip partners;
complex conjugation gives negative heights.

| Right-half-strip value in rh2 | Published value | Source |
|---|---|---|
| 0.808517 + 85.699348i | same six decimals | [J], p.747; repeated [K], p.2045 |
| 0.650830 + 114.163343i | same six decimals | [J], p.747 |
| 0.574356 + 166.479306i | same six decimals | [J], p.747 |
| 0.724258 + 176.702461i | same six decimals | [J], p.747 |
| 0.869531 + 240.404672i | 0.86953 + 240.4046i | [K], p.2046 |
| 0.819550 + 320.876490i | 0.81955 + 320.8764i | [K], p.2046 |

The last two agree at the precision printed, apparently truncated in height;
their additional rh2 digits are not verified by this publication. The next
listed height in [K] is 331.0502, beyond the requested cutoff 330.

### Entry 48

The arXiv record for [A] still lists only v1, dated 2026-02-03, and its §5
table still prints `0.0209081`. No erratum addressing that entry was found
in the arXiv record or targeted title/identifier searches. rh2 reports
`0.00209081` for its N=100 reconstruction. Because the source describes the
entries as **upper bounds**, a smaller reproduced discrepancy does not
logically falsify the printed value. A decimal-place typo is plausible,
not an established correction or author-confirmed erratum.

**Takeaway:** use the complete operator hypotheses, cite the DH values at
their published precision, and label entry 48 a suspected typo.

## Three analytic leads for the cross-term problem

1. **Continuous-kernel residuals, rather than norms of separate places.**
   [H], (2.9)–(2.10), provides a continuous-kernel representation; [F],
   Lemma 1, supplies its variational predecessor. For a projected E-image,
   attempt an independently derived bound for its complement residual in
   the form dual norm. In operator notation the target is
   `‖C⁻¹/² P_U Aλ g_j‖≤K√ℓ_j`, with verified domain or weak-form meaning.
   One must extend from smooth compactly supported vectors to the actual
   E-images, include boundary defects, and control joint columns rather
   than just each j. Such a lemma would be useful even on a restricted
   support interval; no all-λ theorem is asserted here.

2. **Certified compact error plus a finite exceptional space.**
   [D], Lemma 6.9 and (134)–(140), is a specific model for turning an
   operator-norm approximation into a form inequality. Extend the compact
   error analysis beyond the prime-free support interval, including a first
   prime shift, and track the exceptional subspace against E(G), not just
   an unconstrained norm. [I], Theorem 1.1, gives a recent alternative
   finite-tail reduction to audit. Neither technique by itself yields the
   needed leakage-relative estimate.

3. **Use semilocal coordinates before estimating.**
   [C], Theorem 4.6, supplies explicit Sonin transport, while [P], §§3–5,
   gives single-prime moments and Jacobi data through q-series. The missing
   extension is a relation between that scaling-operator representation
   and the *Weil-form* E/complement cross block, with estimates on the
   transport and its support defect. Positive moment matrices for the
   scaling model alone do not prove QW positivity.

## Recommended next shot: an analytic, fixed-equation control

Use the even conductor-5 family in Q5, before attempting another attribution
experiment. Fix a real `σ*=3/4` and put

\[
 G(s)=(1+\sqrt5\,5^{-s})\zeta(s),\qquad
 t_*={L(\sigma_*,\chi_5)\over
 L(\sigma_*,\chi_5)-G(\sigma_*)}.
\]

This specifies a control without a fitted support threshold. In fact
`0<t*<1`: ζ(σ)<0 for 0<σ<1, and the integral representation of the
quadratic character has positive numerator
`z-z²-z³+z⁴=z(1-z)(1-z²)` for `z=e^-u∈(0,1)`, divided by `1-z⁵`.
Thus `L(σ,χ5)>0`. By construction `F_t*(3/4)=0`; the fixed functional
equation also gives a zero at 1/4. This is a no-Euler control with exact
off-line zeros, not a numerical conjecture.

**The falsifiable analytic prediction to attempt next:** construct an explicit
compactly supported smooth **odd** test function and a finite support bound
for which the full Weil form of this F_t* is negative, using a controlled
remainder estimate for the other zeros/local terms. Success requires an
actual negative inequality with an explicit support, not existence alone.
Failure to establish such a remainder estimate is a clean failure of this
specific method; it is not evidence for positivity.

There is a necessary control the current cosine instrument would miss.
Writing `γ=(ρ-1/2)/i`, the manufactured real pair has `γ=±i/4`.
For a real even test function its pair contribution is
`2 f̂(i/4)²≥0`. For a real odd test function it is
`-2 f̂(i/4)²≤0`. This follows directly from reflection parity and the
Hermitian zeros form. Consequently use both parity sectors and include the
varying pole residues. Merely running the existing even-only solver would
not test the claimed prediction. This proposal is grounded in [K], §3,
but the explicit compact-support inequality is not supplied there.

## Corrections and review notes for the repository owner

These are review findings only; the original documentation and code were not
edited on this branch.

- **Pole naming and the zeros sum.** In `CONNES_LETTER.md`, Steps 1–2,
  distinguish the geometric form `arch-primes` from the zeros form obtained
  by **adding** the pole correction. [B], (2.8), and [G], (3.2), include it.
  `connes_letter_mp.py`'s `nopole` variant does add `2vvᵀ`, consistently with
  the zeros form; its name and “drop the pole” prose can mislead. Without RH,
  the Hermitian sum is `Σ φ̂(γ̄) overline(φ̂(γ))`, not unconditional absolute
  squares. This is a terminology/logical correction, not a finding that the
  table-reproducing implementation computes the wrong matrix.

- **Finite positivity and crossover.** Replace continuum positivity “for
  x≤30” by positivity in the tested finite spaces. Negative finite evidence
  near 30.83 supplies an upper bound on a possible first continuum failure;
  it does not supply a certified lower bound 30. Continuity needs an
  operator argument. [H]'s ζ proof is a possible model to extend to DH,
  not a theorem already checked for DH here.
  Nor does the saved angle table support immediate orthogonality at the
  first negative value: at x=32 it reports negative ε with
  `sin²(angle(k,θ))≈4.2e-7`; orthogonality is reported by x=34. Keep the
  first sign change distinct from the later large change of eigenvector.

- **Already documented reconstruction.** Credit [B], §3, and [C], §3.6,
  for multi-vector E/prolate reconstruction. Keep the DH extension and the
  numerical leakage constants separate from those precedents.

- **Uniform convergence and relative margin.** The measured zero error is
  not a theorem of uniform square-root leakage. Nor does raw ε→0 imply
  `μ_min(Q,M)→0`; Q=M gives the elementary counterexample μ=1 while M's
  smallest eigenvalue may tend to zero. Close k/ground-state angle and
  positivity are not literally equivalent: Q can be positive while its
  minimum occurs far from k.

- **Schur lower bound has a sign condition.** `SCHUR_BLOCKS.md` calls
  `min(c,s)/(1+b)²` guaranteed, and `schur_blocks_mp.py:91` records it
  whenever c>0. It is a valid lower bound when `m=min(c,s)≥0`, but not when
  m<0. An exact counterexample is
  `H=[[1,2],[2,1]]`: c=1, s=-3, b=2, and λmin(H)=-1, below the claimed
  bound -1/3. No numerical experiment is needed for this counterexample.
  From `H=T*diag(S,C)T`, `T(y,z)=(y,z+Xy)`, one gets the conservative bounds
  `m/(1+b)²` for m≥0 and `m(1+b)²` for m<0, since both `‖T‖` and
  `‖T⁻¹‖` are at most 1+b. The displayed table leaves negative bounds
  blank, but the code's saved negative “bound” still needs qualification.

- **Ill-conditioned attribution.** The exact split remains well-defined;
  its cancellations do not prove nonexistence of a better representation.
  The rejected axis mask rejects that numerical mask prediction. Because
  it need not preserve a functional equation, it does not settle the broad
  causal role of an Euler product.

- **Family and fact-check corrections.** Qualify conductor-5 rigidity by
  parity, real coefficients, periodicity, normalization, and fixed root
  number. Label the two extended DH decimal values as repository values,
  not digits published in [K]. Label entry 48 as suspected, not corrected.

## References and access record

The sections above use these exact editions. No source was quoted verbatim.

- **[A] Alain Connes**, *The Riemann Hypothesis: Past, Present and a Letter
  Through Time* (2026), [arXiv:2602.04022v1](https://arxiv.org/abs/2602.04022v1).
  Inspected §§4.1, 5–7; table p.25; Theorems 6.1, 7.1; Figure 1; (17), (22).
- **[B] Alain Connes and Caterina Consani**, *Spectral Triples and Zeta-Cycles*
  (2021 preprint; 2023 journal), [arXiv:2106.01715](https://arxiv.org/abs/2106.01715),
  [DOI:10.4171/LEM/1049](https://doi.org/10.4171/LEM/1049).
  Inspected §§2–3; (1.1), (2.8), (2.11), Corollary 2.4; (3.4), Figures 26–36.
- **[C] Alain Connes, Caterina Consani and Henri Moscovici**, *Zeta zeros and
  prolate wave operators: semilocal adelic operators* (2023 preprint; 2024
  journal), [arXiv:2310.18423v2](https://arxiv.org/abs/2310.18423v2),
  [DOI:10.1007/s43034-024-00388-z](https://doi.org/10.1007/s43034-024-00388-z).
  Inspected §§3.6, 4.5–4.8; (29), Theorem 4.6.
- **[D] Alain Connes and Caterina Consani**, *Weil positivity and Trace formula,
  the archimedean place* (2020 preprint; 2021 journal),
  [arXiv:2006.13771v1](https://arxiv.org/abs/2006.13771v1),
  [DOI:10.1007/s00029-021-00689-4](https://doi.org/10.1007/s00029-021-00689-4).
  Inspected Theorem 1, Theorem 4.7, §§6.1–6.7, Lemmas 6.9–6.10,
  Theorem 6.11 and Remark 6.12, pp.47–49 in the preprint.
- **[E] Alain Connes and Walter D. van Suijlekom**, *Quadratic Forms, Real
  Zeros and Echoes of the Spectral Action* (2025),
  [arXiv:2511.23257v1](https://arxiv.org/abs/2511.23257v1),
  [DOI:10.1007/s00220-025-05493-1](https://doi.org/10.1007/s00220-025-05493-1).
  Inspected §§1–3, 5–6; Corollary 1.1, Theorems 5.6, 6.1 (p.13).
- **[F] Enrico Bombieri**, *Remarks on Weil's quadratic functional in the
  theory of prime numbers, I*, Rend. Lincei Mat. Appl. 11 (2000), 183–233
  (published 2001), [primary archive](https://www.bdim.eu/item?id=RLIN_2000_9_11_3_183_0).
  Full PDF retrieved; inspected §§4–7, 12, especially Lemma 1 and
  Theorems 3, 5, 12. The follow-up *A variational approach to the explicit
  formula* (2003), [DOI:10.1002/cpa.10089](https://doi.org/10.1002/cpa.10089),
  was not independently read; its theorem numbering was cross-referenced
  only through [H].
- **[G] Alain Connes, Caterina Consani and Henri Moscovici**, *Zeta Spectral
  Triples* (2025), [arXiv:2511.22755](https://arxiv.org/abs/2511.22755).
  Inspected §3.2, Theorem 3.6/Corollary 3.7, Theorem 5.10, §8.
- **[H] Masatoshi Suzuki**, *Weil's quadratic form via the screw function*
  (2026), [arXiv:2606.09096v3](https://arxiv.org/abs/2606.09096v3),
  version dated September 24. Inspected §§1–5, 8; Theorems 1.1–1.4,
  (2.9)–(2.10). Recent preprint; no independent proof audit.
- **[I] Xuefeng Zhu**, *Weil positivity in compact windows: a finite
  reduction, certified two-sided bounds, and a Landau–Widom decay law*
  (2026), [arXiv:2608.24827v2](https://arxiv.org/abs/2608.24827v2).
  Inspected §§1–7, 11–14; Theorems 1.1–1.4, 6.2, Corollary 6.3,
  Remark 3.3 and Conjecture 12.1. Recent preprint; certificates not executed.
- **[J] Robert Spira**, *Some zeros of the Titchmarsh counterexample* (1994),
  Math. Comp. 63, 747–748,
  [DOI:10.1090/S0025-5718-1994-1254148-8](https://doi.org/10.1090/S0025-5718-1994-1254148-8).
  Full AMS PDF retrieved and read.
- **[K] Eugenio P. Balanzario and Jorge Sánchez-Ortiz**, *Zeros of the
  Davenport–Heilbronn counterexample* (2007), Math. Comp. 76, 2045–2049,
  [DOI:10.1090/S0025-5718-07-01999-0](https://doi.org/10.1090/S0025-5718-07-01999-0).
  Full AMS PDF retrieved; inspected all sections, Theorems 1–2 and (1)–(8).
- **[L1] Jerzy Kaczorowski and Alberto Perelli**, *On the structure of the
  Selberg class, I: 0≤d≤1* (1999), Acta Math. 182, 207–241,
  [DOI:10.1007/BF02392574](https://doi.org/10.1007/BF02392574).
  Full journal PDF via Tsinghua archive; inspected Theorems 1–3 and §8.
  **[L2] K. Soundararajan**, *Degree 1 elements of the Selberg class*
  (2003 preprint; 2005 journal), [arXiv:math/0306300](https://arxiv.org/abs/math/0306300),
  [DOI:10.1016/j.exmath.2005.01.013](https://doi.org/10.1016/j.exmath.2005.01.013).
  Inspected axioms and the theorem on p.2.
- **[M] Andreas Strömbergsson and Anders Södergren**, *On the location of the
  zero-free half-plane of a random Epstein zeta function* (2013 preprint),
  [arXiv:1305.1333](https://arxiv.org/abs/1305.1333); inspected author version
  September 2, 2017, [author PDF](https://www2.math.uu.se/~ast10761/papers/epstein_zf.pdf).
  Inspected Introduction, (1.1)–(1.3); used for the dual-lattice equation,
  not its high-dimensional probability results.
- **[N] Maria Stella Adamo, Karl-Hermann Neeb and Jonas Schober**,
  *Reflection positivity and Hankel operators — the multiplicity free case*
  (2021 preprint; 2022 journal), [arXiv:2105.08522v1](https://arxiv.org/abs/2105.08522v1),
  [DOI:10.1016/j.jfa.2022.109493](https://doi.org/10.1016/j.jfa.2022.109493).
  Inspected §§3.1–3.3, Theorem 3.7 and statement of Theorem 4.1.
- **[O1] Hiroyuki Yoshida**, *On Hermitian forms attached to zeta functions*
  (1992), Adv. Stud. Pure Math. 21, 281–325,
  [DOI:10.2969/aspm/02110281](https://doi.org/10.2969/aspm/02110281).
  Original full text inaccessible; account in [H], §1.1, explicitly marked above.
  **[O2] Enrico Bombieri and Amit Ghosh**, *Around the Davenport–Heilbronn
  function* (2011), Russian Math. Surveys 66(2), 221–270,
  [DOI:10.1070/RM2011v066n02ABEH004740](https://doi.org/10.1070/RM2011v066n02ABEH004740).
  Full text inaccessible in this pass; coverage remains unresolved.
  **[O3] H. S. A. Potter and E. C. Titchmarsh**, *The zeros of Epstein's
  zeta-functions* (1935), Proc. London Math. Soc. s2-39, 372–384,
  [DOI:10.1112/plms/s2-39.1.372](https://doi.org/10.1112/plms/s2-39.1.372).
  **[O4] Paul T. Bateman and E. Grosswald**, *On Epstein's zeta function*
  (1964), Acta Arith. 9(4), 365–373,
  [DOI:10.4064/aa-9-4-365-373](https://doi.org/10.4064/aa-9-4-365-373).
  Both originals inaccessible; metadata verified at their publisher/archive.
  **[O5] H. M. Stark**, *On the zeros of Epstein's zeta function* (1967),
  Mathematika 14, 47–55,
  [DOI:10.1112/S0025579300008007](https://doi.org/10.1112/S0025579300008007),
  [publisher record](https://www.cambridge.org/core/journals/mathematika/article/abs/on-the-zeros-of-epsteins-zeta-function/B5CD35507D2F5544FCB58382D0815822).
  Only indexed p.47 introduction inspected; used solely for the explicitly
  attributed historical summary, not as a reading of [O3]–[O4].
- **[P] Alain Connes, Caterina Consani and Henri Moscovici**, *On q-series
  and moment problem associated to local factors* (2024),
  [arXiv:2403.01247v1](https://arxiv.org/abs/2403.01247v1).
  Inspected Introduction and §§3–5, in particular Propositions 3.1–3.2,
  4.1 and the Hankel determinant construction in §5.1: single-prime
  moment/Jacobi construction.

## Five-line summary

Known: E/prolate reconstruction of several low eigenvectors is already documented, principally numerically.
Known bounds: small-window results and a recent prime-inclusive certificate claim exist; no universal leakage comparison was located.
Not found: an unconditional leakage-weighted E/complement cross bound; inaccessible named sources remain explicit reading gaps.
Corrections: finite positivity, DH theorem hypotheses, negative Schur bounds, “ill-posed,” and conductor-5 rigidity all need qualifications.
Best next shot: use the explicit fixed-even-equation conductor-5 family to seek a certified odd-test negative form at finite support.
