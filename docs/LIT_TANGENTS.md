# Literature tangents: double covers, squares, decay of the window minimum

Date: 2026-10-01. Branch `feat/weil-gram-instrument`. Research only: this file is the only one written, and it is not committed.

**Status labels.**
- **[R]** I read the relevant passage in the source (full text, or the stated section).
- **[A]** I saw only the abstract or metadata.
- **[S]** Taken from another paper's account of the source; that paper is named.
- **[D]** Derived or computed here, not stated in any source. Each is checked as described.
- **[U]** Unverified.

Quotes are at most one per source and under 15 words. Everything else is paraphrase.

---

## Decay rate of the window minimum (priority section)

**Question.** Has anyone stated how fast the smallest eigenvalue of Weil's form on a finite window decays with the support? Specifically:
- an e^{−4πλ²} = e^{−2c} law, c = 2πλ² = 2πx;
- a conductor scaling x/q;
- a power-law prefactor c^5 (even sector) or c^7 (odd sector).

**Answer in brief.**
- **No source proves a decay law for the minimum itself.**
- **The closest statement is Connes 2026, §6.4.** It gives the e^{−4πx} law, with exponent 9/2 on x, for the prolate quantity 1 − χ₂. It compares the Weil minimum ε(λ) with 1 − χ₂ only through a figure.
- **Zhu conjectures a different law.** His Landau–Widom form is −ln λ* ≈ 2π² N(T*)/ln N(T*). Asymptotically it is steeper than 4πx, by a factor of π/2.
- **No source states the x/q collapse for Weil's form.**
- **The Fuchs–Slepian prefactor is c^{n+½}.** That gives 9/2 for h₄ and 13/2 for h₆, and no integer index gives exactly 5 or 7. Our certified minima track h₄ in the even sector and h₆ in the odd sector, each with a slowly growing ratio (table below). That accounts for both fitted exponents, γ ≈ 9/2 + ½ and γ ≈ 13/2 + ½.

### What each source states

| Source | States a rate for the window minimum? | What it does state |
|---|---|---|
| **A. Connes**, *The Riemann Hypothesis: Past, Present and a Letter Through Time*, [arXiv:2602.04022v1](https://arxiv.org/abs/2602.04022), §6.4, Fig. 1, footnote 19 [R] | **For 1 − χ₂ only; ε(λ) is compared by a figure** | (1) Numerics in [Connes–Consani 2023] show that ε(λ) → 0 exponentially in µ = λ². (2) Fig. 1 shows a "striking similarity" between log ε and log(1 − χ₂). (3) The asymptotic 1 − χ₂ ∼ (2¹⁴/3)√2 π⁵ exp(−4πe^L + (9/2)L), L = 2 log λ, so e^L = λ² = x, cited to Fuchs 1964, Theorem 1. (4) χ₂ belongs to h₄ = PS₄,₀(2πλ², x/λ), and the near-radical vector k_λ uses h₀ and h₄. No conductor version. |
| **Connes–Consani**, *Spectral triples and ζ-cycles*, Enseign. Math. 69 (2023) 93–148, [arXiv:2106.01715](https://arxiv.org/abs/2106.01715), Intro and §2.5, Figs. 18–21 [R] | **Qualitative only** | (1) log of the smallest even and odd eigenvalue "decays exponentially" as a function of µ = exp L = λ²; no constant is given. (2) The number of small eigenvalues grows roughly like µ, and the odd matrix has one fewer. (3) Intro: the Slepian–Pollak angle operator has 1 + ν(λ²) ∼ 2λ² extremely small eigenvalues, with prolates PS₂ₘ,₀(2πλ², x/λ). (4) The conditions f(0) = f̂(0) = 0 must be imposed. |
| **Connes–Consani**, *Weil positivity and trace formula, the archimedean place*, Selecta Math. 27 (2021) no. 4, paper 77, [arXiv:2006.13771](https://arxiv.org/abs/2006.13771), §4, eq. (68)–(71) [R] | **No** | (1) Single archimedean place, bandwidth fixed at c = 2π. (2) Eq. (71) bounds the decay of |λ(n)| in the index n at fixed c, not in λ. (3) No e^{−4πλ²} statement and no conductor. |
| **Connes–Consani–Moscovici**, *Zeta zeros and prolate wave operators: semilocal adelic operators*, Ann. Funct. Anal. 15 (2024) no. 4, paper 87, [arXiv:2310.18423](https://arxiv.org/abs/2310.18423), §5.1, Prop. 5.1 [R] | **No** | (1) Prop. 5.1, through the metaplectic representation of the double cover of SL(2,ℝ): **W_λ = h² + 4πλ² k − ¼** in U(sl₂). (2) k is the Hermite operator, with spectrum n + ½. (3) L²(ℝ)_ev and L²(ℝ)_odd are irreducible, of lowest weight ½ and 3/2. (4) No eigenvalue-decay statement and no Dirichlet characters. |
| Connes–Consani–Moscovici, *Zeta Spectral Triples*, [arXiv:2511.22755](https://arxiv.org/abs/2511.22755) [A] | No | Rank-one perturbations whose spectra match the low zeros. No decay law in the abstract. |
| **W. H. J. Fuchs**, J. Math. Anal. Appl. 9 (1964) 317–330, [DOI](https://doi.org/10.1016/0022-247X(64)90017-4); **D. Slepian**, J. Math. Phys. 44 (1965) 99–140, [Wiley](https://onlinelibrary.wiley.com/doi/10.1002/sapm196544199) [U: originals not accessible] | No (prolates only) | The standard asymptotic 1 − λ_n(c) ∼ 4√π 8ⁿ c^{n+½} e^{−2c}/n!. Not read in the originals. [D] With n = 4 and 1 − χ = (1 − λ)/(1 + √λ) ≈ (1 − λ)/2, it reproduces Connes' constant (2¹⁴/3)√2 π⁵ ≈ 2.36354·10⁶ to all printed digits, so the form is consistent with his citation. |
| **H. J. Landau, H. Widom**, *Eigenvalue distribution of time and frequency limiting*, J. Math. Anal. Appl. 77 (1980) 469–481 [S: cited by Osipov, arXiv:1206.4541; and by Zhu] | No | Width of the "plunge region" of the prolate spectrum. Zhu says his constant 2π² matches its plunge rate [S: Zhu]. |
| **X. Zhu**, [arXiv:2608.24827v2](https://arxiv.org/abs/2608.24827), Conj. 12.1, Thm. 1.3, Table 3 [R] | **Conjecture plus one conditional bound** | (1) Conj. 12.1: −ln λ*(L) = C·N(T*)/ln N(T*)·(1 + o(1)), C ≈ 2π², T* = 2πe^{2L}, L = half-width. Fitted on L ≤ 1.2 and checked on L = 1.4–2.0, with residuals < 0.7%. (2) Thm. 1.3, assuming RH: λ*(L) ≤ exp(−L e^L) for large L. (3) Table 3 upper bounds: λ* ≤ 3.2·10⁻²⁸³ at L = 2.0. No conductor discussion. |
| **T. Kim et al.**, Suzuki's operator by finite elements, [arXiv:2607.24830](https://arxiv.org/abs/2607.24830), R1, R2, R7, §3.8, §5.4 [R] | **Qualitative only; listed as open** | (1) R7: λ₁(a) > 0 decays "superexponentially". (2) They say it is consistent with e^{−ca²} and extrapolate λ₁ ∼ 10⁻²⁵ at a = 2. (3) §5.4 lists "the exact form of the decay law" as open. (4) R1/R2 concern only the prime-free regime 2a < log 2: λ_k(a) = log(1/a) + log(k − ½) + B₀ + O(a), B₀ = log q − 2 log 2. |
| **M. Suzuki**, *Weil's quadratic form via the screw function*, [arXiv:2606.09096v3](https://arxiv.org/abs/2606.09096), Thms. 1.3–1.4, §1.1 [R] | No (small a only) | (1) Thm. 1.3: λ_a is continuous in a. (2) Thm. 1.4: for small a, λ_a = log(1/a) + µ₁ − log 2π + ψ(2) − 1 + O(a), with a positive, simple, even ground state. (3) §1.1's account of Yoshida: see the Yoshida row. |
| **H. Yoshida**, *On Hermitian forms attached to zeta functions*, Adv. Stud. Pure Math. 21 (1992) 281–325 [S: Suzuki §1.1] | No | (1) Positivity for small a (Lemma 2). (2) RH ⇔ non-degeneracy on a completed space (Thm. 2). (3) Prop. 1: Q ≥ 0 on all odd v implies RH; on all even v, RH except possibly real zeros. (4) Works on K(a), restrictions of 2a-periodic functions. |
| **E. Bombieri**, Rend. Lincei 11 (2000) 183–233 [S: `docs/LITERATURE_PASS.md`] | No (small support only) | (1) Thm. 12: for support length b < log 2, Q ≥ (log(1/b) − log log(1/b) − O(1))‖F‖². (2) Thm. 3: the minimum is attained. |
| **J.-F. Burnol**: C. R. 333 (2001) 201–206, [arXiv:math/0105120](https://arxiv.org/abs/math/0105120); C. R. 335 (2002) 689–692; JTNB 16 (2004) 65–94 [A; as cited by CCM 2024] | No | Sonine and de Branges spaces from the Fourier transform, and vectors attached to the zeros. Nothing on the minimum versus support in the abstracts. |
| **V. Liu**, *Certified Weil Positivity Beyond the Unit Window*, [alphaxiv](https://www.alphaxiv.org/abs/2609.weil-positivity-riemann-zeta-bounds) [A] | No | (1) Certified constants 2⁻¹⁵¹ at half-width 1 and 2⁻⁴⁹¹⁶² at half-width 17/16. (2) A −1/8 "threshold" beyond log 8/2 belongs to one localisation strategy. [D] 2⁻¹⁵¹ ≈ 3.5·10⁻⁴⁶ is about 16 orders below our 4.6·10⁻³⁰ at a = 1 (`docs/GRID_NORM.md`), so these constants say nothing about the true decay. |

**Conductor scaling.** No paper found states λ_χ(x) ≈ F(x/q) for Weil's window form.
- **Literature mechanism consistent with it.** For primitive χ mod q, the twisted theta and Poisson formula behind the functional equation, with factor (q/π)^{s/2}, puts the self-dual point at 1/√q. That is standard (e.g. Davenport, *Multiplicative Number Theory*, Ch. 9; [U]: not re-read here).
- [D] So the prolate bandwidth becomes c = 2πλ²/q.
- This work already found the DH floor at c = 2πx/5 (`docs/CONNES_LETTER.md`).
- Kim et al.'s B₀ = log q − 2 log 2 is a different regime: small a, additive in log q.
- [D] Consistency check. `docs/EPSTEIN.md` reports ζ_K's ε₀ ∼ e^{−0.5x} for x ∈ [20, 40]. The conductor-20 rate is 4π/20 = 0.628, and a prefactor x^{γ} lowers the local slope to 0.628 − γ/x. That is ≈ 0.51 at x = 30 for γ ≈ 3.5. So ζ_K's minimum decays at its conductor-20 factor's floor rate (see Test 1).

### Our certified minima against the Fuchs–Slepian leakages [D]

The table uses the leakage ℓ_n(x) := ½·4√π 8ⁿ c^{n+½} e^{−2c}/n!, with c = 2πx. Certified lower bounds ("cert") and finite-basis upper bounds ("ub") come from `docs/AUDIT_ZHU.md`, `docs/CERTIFICATE_238.md` and commit 108d678.

| Support 2a | x | Even cert | Even ub | Even cert/ℓ₄ | Even ub/ℓ₄ | Odd cert | Odd ub | Odd cert/ℓ₆ | Odd ub/ℓ₆ | Odd cert/ℓ₂ |
|---|---|---|---|---|---|---|---|---|---|---|
| 1.6 | 4.953 | 1.0277e-17 | 1.6356e-17 | 3.49 | 5.55 | 9.118e-15 | — | 1.50 | — | 1.6e7 |
| 2.38 | 10.805 | 6.813e-48 | 1.139e-47 | 5.98 | 9.99 | 3.920e-44 | 6.37e-44 | 3.50 | 5.68 | 8.5e8 |
| 2.6 | 13.464 | 5.763e-62 | 1.04e-61 | 6.09 | 10.98 | 5.360e-58 | 9.74e-58 | 3.71 | 6.74 | 2.2e9 |

**Reading.**
- **Even sector tracks ℓ₄, odd sector tracks ℓ₆.** The even minimum stays within a factor of 3.5–11 of ℓ₄ (Connes' h₄). The odd minimum stays within a factor of 1.5–6.7 of ℓ₆, the second Fourier-(−1) prolate after h₂. It is 10⁷–10⁹ above ℓ₂, so it is not h₂.
- **This is a hypothesis.** One possible reason (not checked): in each Fourier class, the condition f(0) = 0 removes the lowest prolate, as h₀ is removed in Connes' k_λ.
- **The sector ratio follows.** odd/even ≈ (0.4–0.6)·ℓ₆/ℓ₄, and ℓ₆/ℓ₄ = (64/30)c² ≈ 2.13c². This explains why the odd margins are 10³–10⁴ larger.
- **On question (2).** The fitted γ_even ≈ 5 and γ_odd ≈ 7 are 9/2 + ½ and 13/2 + ½. The common extra ½ is the slow growth of the ratio ε/ℓ_n (3.5 → 6.1 even, 1.5 → 3.7 odd, over x = 5 → 13.5). I found no prolate variant with an integer exponent. Slepian IV (Bell Syst. Tech. J. 43 (1964) 3009–3057, generalized prolates) could not be retrieved (archive.org offline) [U].

**Zhu's law against the prolate law [D].**
- **Asymptotically.** N(2πx)/ln N(2πx) → x, so Zhu's −ln λ* ∼ 2π²x = 19.74x, against 4πx = 12.57x.
- **Over the computed range they mimic each other.** The local slope of Zhu's form is 11.0 at x = 12, 12.0 at x = 30 and 12.8 at x = 100. The prolate form, 4π − γ/x with γ = 5, gives 12.15 at x = 12 and 12.40 at x = 30.
- **At the certified points.** The values exp(−2π²N/ln N) are 1.8e-25 at x = 4.95 (far off), 3.4e-50 at x = 10.8 (below the certified lower bound 6.8e-48) and 9.7e-62 at x = 13.46 (inside the bracket [5.8e-62, 1.04e-61]).
- **At x = e^{2.99} the two predictions separate by a factor of 4–8** (Test 2).

**Kim et al.'s extrapolation is inconsistent [D].** Their a and Zhu's L are both the half-width: Kim's first prime threshold a = log 2/2 is the same convention. Their 10⁻²⁵ at a = 2 lies 258 orders of magnitude above Zhu's certified upper bound 3.2·10⁻²⁸³ at L = 2.0. Their measured slopes d log λ₁/da ≈ −20 → −56 on a ∈ [0.2, 0.5] are close to d/da[−4πe^{2a} + 9a] = −28.5 → −59.3 at the upper end. So the double-exponential law fits their own data better than e^{−ca²}.

**Metaplectic link [D, from CCM Prop. 5.1 and eq. (64)].**
- **The rate is the coefficient of the rotation generator.** e^{−2c} = e^{−4πλ²}, and 4πλ² is the coefficient of the compact (rotation) generator k in the prolate operator.
- **k has half-integer spectrum n + ½.** So the metaplectic lift of a 360° rotation acts as −1, and only the 720° lift is the identity.
- **Our sector split is the n mod 4 class of the Fourier transform.** In Connes–Consani, n ≡ 0 for the even sector and n ≡ 2 for the odd sector, both inside the weight-½ piece L²(ℝ)_ev. Odd Dirichlet characters (Γ_ℝ(s+1)) use the weight-3/2 piece L²(ℝ)_odd.
- **Connection to the user's idea.** This is the one place in the literature where a genuine spin-½ (metaplectic) structure governs our window minimum.

---

## A. Metaplectic double cover and positivity of central values

**Sources.**
- **Waldspurger**, *Sur les coefficients de Fourier des formes modulaires de poids demi-entier*, J. Math. Pures Appl. 60 (1981) 375–484 [A: citation via [Wikipedia](https://en.wikipedia.org/wiki/Waldspurger%27s_theorem)]. Fourier coefficients of half-integral weight forms are identified with central values of twists. Kohnen–Zagier (below) record that Waldspurger's proof uses the representation theory of the metaplectic group and fixes only ratios of coefficients.
- **Kohnen–Zagier**, *Values of L-series of modular forms at the center of the critical strip*, Invent. Math. 64 (1981) 175–198, [scan](https://people.mpim-bonn.mpg.de/zagier/files/doi/10.1007/BF01389166/fulltext.pdf), pp. 175–179 [R].
  - **Thm. 1.** For a level-1 eigenform f of weight 2k and its half-integral weight partner g ∈ S⁺_{k+½}(Γ₀(4)): c(|D|)²/⟨g,g⟩ = ((k−1)!/π^k)|D|^{k−½} L(f,D,k)/⟨f,f⟩, for fundamental D with (−1)^k D > 0.
  - **Cor. 1.** L(f,D,k) ≥ 0. They note one expects ≥ 0 because "the contrary would contradict the Riemann hypothesis for L(f, D, s)".
  - They add that even Ramanujan–Petersson gives no direct proof, since s = k lies half a unit left of absolute convergence. The sign of c(n) is "still utterly mysterious".
- **Shimura**, *On modular forms of half integral weight*, Ann. of Math. 97 (1973) 440–481 [A].
- **A. Weil**, *Sur certains groupes d'opérateurs unitaires*, Acta Math. 111 (1964) 143–211 [A]. The Weil (oscillator) representation, a genuine representation of the metaplectic double cover Mp.
- **Gan–Qiu–Takeda**, *The regularized Siegel–Weil formula (the second term identity) and the Rallis inner product formula*, Invent. Math. 198 (2014) 739–831, [Springer](https://link.springer.com/article/10.1007/s00222-014-0509-0) [A]. The Rallis inner product formula for global theta lifts of any dual pair, and its non-vanishing consequences. The norm of a theta lift is expressed through an L-value, which forces a sign.
- **Lapid–Rallis**, *On the nonnegativity of L(½, π) for SO₂ₙ₊₁*, Ann. of Math. 157 (2003) 891–917, [arXiv:math/0402371](https://arxiv.org/abs/math/0402371), Intro [R].
  - L(½, π) ≥ 0 for cuspidal generic π on SO(2n+1).
  - Method: a residue of an intertwining operator of Eisenstein series is positive semi-definite. This is not a theta or square formula.
  - They note GRH gives L(s, π) > 0 on (½, 1], and that central nonnegativity is unknown even for quadratic Dirichlet characters.
  - Applications cited: subconvexity, and the Gross–Prasad and BSD contexts. Lapid also has *Positivity of L(½, π) for symplectic representations*, C. R. 2002, [ScienceDirect](https://www.sciencedirect.com/science/article/pii/S1631073X02022173) [A].
- **Iwaniec 1987**, *Fourier coefficients of modular forms of half-integral weight*, Invent. Math. 87, 385–401 [A: keywords include Kloosterman and Salié sums]. **Duke 1988**, *Hyperbolic distribution problems and half-integral weight Maass forms*, Invent. Math. 92, 73–90 [A]. Standard attribution, not re-read here: a nontrivial bound for half-integral weight coefficients, combined with Waldspurger, gives subconvexity for quadratic twists, and Duke uses it for equidistribution [U].
- **Conrey–Iwaniec**, *The cubic moment of central values of automorphic L-functions*, Ann. of Math. 151 (2000) 1175–1216 [A]. Uses nonnegativity of the central values in a moment-plus-amplifier argument for subconvexity.
- **Iwaniec–Sarnak**, *The non-vanishing of central values of automorphic L-functions and Landau–Siegel zeros*, Israel J. Math. 120 (2000) 155–177, [Springer](https://doi.org/10.1007/s11856-000-1275-9) [A]. At least 50% of L(½, f) are positive, and any improvement on 50% is tied to Landau–Siegel zeros.
- **Conrey–Soundararajan**, *Real zeros of quadratic Dirichlet L-functions*, Invent. Math. 150 (2002) 1–44 [R, pp. 1–3].
  - Thm. 1: for ≥ 20% of odd squarefree d, L(σ, χ₋₈d) > 0 on [0, 1].
  - They stress that before this, nothing excluded every large-conductor L(s, χ) having a real zero.
- **Shankar–Södergren–Templier**, *Central values of zeta functions of non-Galois cubic fields*, [arXiv:2107.10900](https://arxiv.org/abs/2107.10900) (Invent. Math., 2025) [A]. Infinitely many non-Galois cubic K have ζ_K(½) < 0. [D] Since ζ(½) ≈ −1.4604 < 0, those fields have L(½, ρ_K) = ζ_K(½)/ζ(½) > 0. A negative Dedekind value at the centre is the GRH-expected sign when there is a pole.
- **Connes–Consani–Moscovici 2024**, Prop. 5.1 (see the priority section) [R]. The prolate operator, whose leakage sets our floor, is an element of U(sl₂) acting through the metaplectic representation. Their abstract announces a forthcoming adelic Weil-representation candidate for the semilocal prolate operator.

**Bearing on this work.**
- **Squares give unconditional nonnegativity at one point.** The double-cover "square" mechanism (Waldspurger, Kohnen–Zagier, Rallis) does this at a single point per L-function, the centre of each twist, and Kohnen–Zagier say it is what GRH predicts.
- **Its known payoffs are not explicit-formula positivity.** They are subconvexity (Iwaniec, Conrey–Iwaniec), nonvanishing proportions, and the 50% barrier tied to Siegel zeros (Iwaniec–Sarnak).
- **No route found from a square to Weil positivity.** I found no source that turns a squared coefficient into positivity of Weil's explicit-formula functional, which sums over all zeros, or into anything beyond GRH's consequences.
- **The GL(1) shadow is weak.** Central positivity for quadratic L(½, χ) is open (Lapid–Rallis, Conrey–Soundararajan). So the "spinor" half of ζ_K = ζ·L(χ) has no known square formula.
- **Where the metaplectic double cover does reach this work.** It enters through the archimedean prolate operator (CCM Prop. 5.1). That fixes the e^{−4πx} floor, but it is a statement about the archimedean place, not about primes.

## B. Positivity from squares and tensor powers; the ceiling of "twice around"

**Sources.**
- **Deligne**, *La conjecture de Weil. I*, Publ. Math. IHÉS 43 (1974) 273–307, [numdam](http://www.numdam.org/item/PMIHES_1974__43__273_0/) [A].
- **T. Feng**, *Notes on Deligne's "La conjecture de Weil. I"*, [PDF](https://math.berkeley.edu/~fengt/Weil_I.pdf), §5 [R]. The account of the key step:
  - ⊗^{2k} has traces Tr(F, E)^{2k}. These are rational and so ≥ 0, which makes t d/dt log det(1 − F_x t, ⊗^{2k}E) a series with positive coefficients.
  - A product of positive series converges no further than any factor.
  - Symplectic invariant theory puts the only pole of Z(U, ⊗^{2k}E) at q^{−kβ−1}. Hence |α| ≤ q^{β/2 + 1/(2k)}, and k → ∞ gives the bound.
  - The argument needs the symplectic form, monodromy open in Sp, and rational traces.
- **Milne**, *The Riemann Hypothesis over Finite Fields: From Weil to the Present Day*, [arXiv:1509.00797](https://arxiv.org/abs/1509.00797), §5 [R].
  - Rankin (1939) obtained |c(n)| = O(n^{k−1/5}) by Landau's theorem.
  - Langlands' remark: if every L(s, σ, π) is analytic for Re s > 1, then Rankin's idea with σ = ρ⊗ρ̄ and tensor powers gives |λ(t_p)| = 1, the Ramanujan bound.
  - Deligne (Milne's transcript) took the method from Rankin and Landau: information on poles controls local factors.
- **Sarnak**, *Nonvanishing of L-functions on ℜ(s) = 1* (for Shalika's 60th birthday), [IAS PDF](https://publications.ias.edu/sites/default/files/ShalikaBday2002.pdf), §1 [R].
  - de la Vallée Poussin's auxiliary D(s) = ζ³(s)ζ²(s+it₀)ζ²(s−it₀)ζ(s+2it₀)ζ(s−2it₀) = **L(s, Π×Π̃) with Π = 1 ⊞ α^{−it₀} ⊞ α^{it₀}**.
  - The double-frequency factors are the cross terms of Π×Π̃.
  - For quadratic χ, D = ζ·L(χ), or its square, gives L(1, χ) ≠ 0 by Landau's lemma, but no conductor-uniform zero-free region. That is the Landau–Siegel problem.
  - The method covers L(s, π), L(s, π×π′) and their functorial products: "more or less all that can be handled by Poussin's method".
  - Jacquet–Shalika, Invent. Math. 38 (1976) 1–16, and Hoffstein–Ramakrishnan, IMRN 1995 no. 6, 279–308, are cited there.
- **Hoffstein–Lockhart**, *Coefficients of Maass forms and the Siegel zero*, Ann. of Math. 140 (1994) 161–181, with the appendix by Goldfeld, Hoffstein and Lieman, "An effective zero free region" [A]. The symmetric-square case: positivity of an auxiliary Rankin–Selberg product excludes Siegel zeros.
- **Mossinghoff–Trudgian**, *Nonnegative trigonometric polynomials and a zero-free region for the Riemann zeta-function*, [arXiv:1410.3926](https://arxiv.org/abs/1410.3926), §1 and §6 [R].
  - R₀ = 5.573412.
  - History table: de la Vallée Poussin 30.4679 with (1 + cos φ)² = 3/2 + 2cos φ + ½cos 2φ; Kadiri 5.69693.
  - Landau's extremal quantity V = lim inf (f(0) − a₀)/(√a₁ − √a₀)²: Arestov–Kondrat'ev (Math. Notes 47 (1990) 10–20) proved 34.468305 < V < 34.5035864, and they improve the upper bound to 34.4889920009.
  - Landau gives R₀ = V/2 + ε. Stechkin gives R₀ = (V/2)(1 − 1/√5) + ε.
  - Their Table 4: the near-optimal degree-40 polynomials do not beat the degree-16 one inside Kadiri's method.
- **Mossinghoff–Trudgian–Yang**, Research in Number Theory 10 (2024), [arXiv:2212.06867](https://arxiv.org/abs/2212.06867) [A]: R₀ = 5.558691. **Bellotti–Trudgian–Yang**, [arXiv:2603.21490](https://arxiv.org/abs/2603.21490) [A]: 4.896, building on Heath-Brown's Linnik-constant work. **H. S. Tan**, [arXiv:2411.01385](https://arxiv.org/abs/2411.01385) [A]: optimal cosine polynomials at degrees 7–8, extending Arestov's degree ≤ 6. **Kadiri**, Acta Arith. 117 (2005) 303–339 [S: MT].
- **Lu–Zaman–Zhao**, *Numerical Computations Concerning Landau–Siegel Zeros*, [arXiv:2602.03626](https://arxiv.org/abs/2602.03626) [A]. No real zero with σ ≥ 1 − 1/(5 log q) for every quadratic χ mod q ≤ 10¹⁰.

**Bearing on this work.**
- **The squares are Rankin–Selberg squares.** [D, elementary] 3 + 4cos θ + cos 2θ = 2(1 + cos θ)² = ½|1 + e^{iθ}|⁴. Sarnak's 3 + 4cos θ + 2cos 2θ = |1 + 2cos θ|². So "going around twice" is the cross term of a Rankin–Selberg square.
- **The ceiling has two parts.**
  1. **Shape.** Positivity plus Landau's lemma gives only zero-free regions of width ∝ 1/log t at Re s = 1. Better shapes (Korobov–Vinogradov) need exponential-sum bounds (Sarnak §1; MT intro).
  2. **Constant.** The pure trig-polynomial main term is capped by V > 34.468: about V/2 ≈ 17.2 (Landau), or (V/2)(1 − 1/√5) ≈ 9.53 (Stechkin) [D, arithmetic]. Kadiri-type error analysis reaches 4.9–5.6.
- **Neither part touches the critical line.**
- **The method works on the opposite side from Weil's form.** de la Vallée Poussin's positivity is on the prime side, Σ c_n n^{−σ}P(t log n) ≥ 0, and holds for any Dirichlet series with c_n ≥ 0, with no RH content. Weil positivity is on the zeros side.
- **Coherence is used in opposite directions.** On a window, our comb bound uses the opposite of the de la Vallée Poussin trick: |cos| ≤ 1 caps the Kronecker coherence that the trick exploits. The quarter law (`docs/QUARTER_LAW.md`) says a window sees about ¼ of it.
- **Deligne's tensor-power trick has no number-field analogue for zeros.** Over ℚ it becomes Ramanujan via functoriality (Langlands, Milne). Deligne could bound zeros because over 𝔽_q every zero is a Frobenius eigenvalue on a finite H^i.
- **For reading (2).** "Twice around" is a known technique with a known ceiling (log-width regions at Re s = 1). In the quadratic case ζ·L(χ) it is exactly where the method stalls on Siegel zeros.
- **Unrelated coincidence.** Our normalisation test function (1 + cos θ)², θ = 2πu/L (`docs/CONTROL_CERTIFICATES.md` §2), is de la Vallée Poussin's polynomial in a different variable.

## C. Hodge index and coverings in function fields; descent to number fields

**Sources.**
- **Milne 2015** [R], §1 and §3:
  - Weil's proof rests on the Castelnuovo–Severi inequality (D·D) ≤ 2dd′ on C × C′.
  - Mattuck–Tate derived it from Riemann–Roch for surfaces: *On the inequality of Castelnuovo–Severi*, Abh. Math. Sem. Hamburg 22 (1958) 295–299.
  - Grothendieck derived it from the Hodge index theorem: *Sur une note de Mattuck–Tate*, J. reine angew. Math. 200 (1958) 208–215.
  - Weil's 1940 proof realises σ(X ∘ X′) > 0 as positivity of the Rosati involution on the Jacobian.
  - **A double cover in Weil's own proof.** To treat a general divisor, Weil passes to a finite generically Galois covering V → C × C, where a determinant Φ changes only by a sign under Galois, and **replaces Φ by Φ²** to descend.
  - Serre (email quoted by Milne): Artin L-functions being polynomials needs ℓ-adic methods, since "the positivity result alone is not enough".
  - Stepanov (1969) and Bombieri, *Counting points on curves over finite fields (d'après S. A. Stepanov)*, Sém. Bourbaki 430, LNM 383 (1974) 234–241, give the elementary Riemann–Roch-only proof; see also Schmidt, Acta Arith. 24 (1973) 347–367.
- **Mumford**, *Prym varieties I*, in Contributions to Analysis (Academic Press, 1974) 325–350 [A]. Standard facts [U, not read]: for an étale double cover Y → X, Jac Y ∼ Jac X × Prym, and P_Y(u) = P_X(u)·L(X, χ, u), where χ is the cover's quadratic character.
- **Connes–Consani–Marcolli**, *The Weil proof and the geometry of the adeles class space*, [arXiv:math/0703392](https://arxiv.org/abs/math/0703392), Manin Festschrift (Birkhäuser 2008/9) [A]. A dictionary in which Weil's explicit formula is a Lefschetz trace formula and RH is "positivity of the relevant trace pairing".
- **Connes–Consani**:
  - *Geometry of the scaling site*, Selecta Math. 23 (2017) 1803–1850, [arXiv:1603.03191](https://arxiv.org/abs/1603.03191) [A]: periodic orbits per prime, theta functions, and Riemann–Roch with real dimensions.
  - *The Riemann–Roch strategy, complex lift of the scaling site*, [arXiv:1805.10501](https://arxiv.org/abs/1805.10501) [A]: adapt Weil's proof to characteristic 0.
  - *Riemann–Roch for Spec ℤ-bar*, [arXiv:2205.01391](https://arxiv.org/abs/2205.01391) [A]; *Riemann–Roch for the ring ℤ*, [arXiv:2306.00456](https://arxiv.org/abs/2306.00456) [A].
  - *Knots, primes and the adele class space*, [arXiv:2401.08401](https://arxiv.org/abs/2401.08401) [A]: periodic orbits C_p of length log p, and the adele class space as maximal abelian cover.
  - *The Absolute Twistor Line and the Geometry of Spec ℤ-bar*, [arXiv:2609.00299](https://arxiv.org/abs/2609.00299) [A]: real Hodge structures and complex conjugation. The abstract does not mention positivity.
- **Deninger**, *Number theory and dynamical systems on foliated spaces*, [arXiv:math/0204110](https://arxiv.org/abs/math/0204110), Thm. 2.1 and the last section [R].
  - For a Riemannian foliation with a conformal flow, the generator on H¹ is Θ = α/2 + S, with S skew for the Hodge-∗ inner product (h, h′) = tr(h ∪ ∗h′). That is an "RH" statement coming from the Hodge star.
  - It gives a dynamical proof for ordinary elliptic curves through a CM lift in which Frobenius is a split prime π, ππ̄ = p, of an imaginary quadratic field. He calls the construction misleading in general.
  - Also *Analogies between analysis on foliated spaces and arithmetic geometry*, [arXiv:0709.2801](https://arxiv.org/abs/0709.2801) [R: no positivity argument found], and ICM 1998, Doc. Math. Extra Vol. I, 163–186 [A].
- **A. Weil**, *Sur les "formules explicites" de la théorie des nombres premiers*, Comm. Sém. Math. Univ. Lund (1952) 252–265 [S: search summaries]. **A. Connes**, Selecta Math. 5 (1999) 29–106, [DOI](https://doi.org/10.1007/s000290050042) [A]. Explicit formula for Hecke characters of K: RH for all abelian L-series of K ⇔ positivity of one distribution on the idele class group.

**Derived facts [D].**
1. **Function-field "Epstein".**
   - For any elliptic curve E/𝔽_q with point O, let R = 𝔽_q[E ∖ O]; its class group is E(𝔽_q).
   - By Riemann–Roch, the principal-class partial zeta is Z_O(u) = 1 + Σ_{n≥2} q^{n−1}uⁿ = (1 − qu + qu²)/(1 − qu), **independent of E**.
   - Each non-principal class gives u/(1 − qu).
   - Summing over the h = q + 1 − a classes recovers ζ_R = P_E(u)/(1 − qu).
   - Checked directly for y² = x³ + 2x over 𝔽₅ (#E = 2): degree-2 coefficients 5 and 5, matching (counted by explicit 𝔽₂₅ point enumeration).
   - **The numerator of our closed-seam "fake", 1 − 5u + 5u² (`docs/CLOSED_SEAM.md`), is exactly this partial zeta's numerator at q = 5.** The partial zeta is one sheet of the class-group cover. Its numerator is self-reciprocal (qu²·P(1/(qu)) = P(u)), although its denominator has only the pole at u = 1/q. It has S₁ = √q > 2 for every q ≥ 5, and its zeros are real, the function-field analogue of a Siegel pair. On the zeros side, the failing eigenvector of T₂ = ℓ[[2, S₁], [S₁, 2]] is (1, −1), the antisymmetric direction on the two-point orbit.
2. **ℚ(√−5) genus decomposition.** With our Z₁ normalisation a_n = r₁(n)/2:
   - **Z₁ = ½[ζ·L(χ₋₂₀) + L(χ₋₄)·L(χ₅)] = ½[ζ_K + L_K(ψ)]**, verified coefficientwise for n ≤ 5000. ψ is the class-group character, and H = ℚ(i, √5) is the Hilbert class field, since h(−20) = 2.
   - So ζ_H = ζ_K·L_K(ψ) is the unramified double cover of K.
   - **Z₁ is one sheet of it, the twin of fact 1.** All three share the same Γ_ℝ(s)Γ_ℝ(s+1) and conductor 20.
3. **Descent by subtraction is strictly lossy.**
   - Q_{ζ_K} = Q_ζ + Q_{L(χ)} exactly: same test function, additive in log L.
   - The best unconditional upper bound on Q_{L(χ)} replaces its comb by −A_L.
   - So the bound Q_ζ ≥ Q_{ζ_K} − U_L has frequency weight Ψ_ζ − comb_L(t) − A_L ≤ Ψ_ζ pointwise.
   - Certifying the cover and subtracting can therefore never beat certifying ζ directly.
4. **Inert primes as antiperiodic orbits.** (1 − p^{−2s})^{−1} = (1 − p^{−s})^{−1}(1 + p^{−s})^{−1}, which is L²(ℝ/2ℓℤ) = periodic ⊕ antiperiodic on ℝ/ℓℤ, with ℓ = log p. The orbit of an inert prime "closes after two loops", and its doubled-orbit spectrum is the scalar factor plus the χ(p) = −1 factor.

**Bearing on this work.**
- **In function fields, positivity descends because the cover's positivity covers twisted test functions too.** Rosati or Hodge-index positivity holds for Y as a whole. The deck involution splits H¹(Y) orthogonally, so restricting to ±-isotypic test functions (f and f·χ) gives positivity of each factor.
- **The number-field version of "twisted test functions" is Weil's Hecke-character explicit formula** (Weil 1952; Connes 1999). Positivity for ζ_K over all characters of K is GRH for each character's L-function separately.
- **No number-field theorem found that descends positivity from a cover to its base.** Fact 3 shows the subtraction route fails for certificates. The Connes–Consani and Deninger programmes aim at a Hodge/Riemann–Roch positivity on the base object itself, not at descent.
- **Fact 1 is an exact, fully finite model of our Euler versus non-Euler separation.** The covering curve passes, and its partial-zeta sheet fails on the first two-point window, in the antisymmetric direction. That failure is forced by Riemann–Roch alone.

## D. Weil positivity on finite windows: extremal problems and records

**Sources.**
- Yoshida (1992), Bombieri (2000; 2003 not read), Suzuki (2026), Connes–Consani (2021, 2023), Connes–Consani–Moscovici (2024) and Zhu: see the priority section and `docs/LITERATURE_PASS.md`.
- **Connes–van Suijlekom**, Commun. Math. Phys. 406:312 (2025), [arXiv:2511.23257](https://arxiv.org/abs/2511.23257) [S: LITERATURE_PASS]: if the minimum is simple, the minimiser's Fourier zeros are all real.
- **Groskin**, *A finite Guinand–Weil dictionary and archimedean tail order for the truncated Weil quadratic form*, [arXiv:2607.02828](https://arxiv.org/abs/2607.02828) [A]. Band-limited test functions whose zero sums equal truncated form values, and a two-sided certification rule. The abstract states no support record.
- **Desogus**, *The Three Gates*, [arXiv:2609.20367](https://arxiv.org/abs/2609.20367) [A]. Claims that a "restricted odd Weil criterion" yields RH. This is an extraordinary claim, unrefereed, and not audited. It is not a support certificate.
- **Beurling–Selberg / Carneiro school.** All conditional on RH or GRH:
  - Carneiro–Chandee–Milinovich, *A note on the zeros of zeta and L-functions*, Math. Z. 281 (2015) 315–332 [A]: the Guinand–Weil formula with extremal one-sided bandlimited approximations, giving S(t), central order and lowest-zero bounds.
  - Carneiro–Milinovich–Soundararajan, *Fourier optimization and prime gaps*, Comment. Math. Helv. 94 (2019) 533–568, [arXiv:1708.04122](https://arxiv.org/abs/1708.04122) [A].
  - Carneiro–Chirre–Milinovich, *Hilbert spaces and low-lying zeros of L-functions*, Adv. Math. 410 (2022) 108748, [arXiv:2109.10844](https://arxiv.org/abs/2109.10844) [A]: one-level-density extremal problems solved in reproducing kernel Hilbert spaces.
- **Platt**, *Numerical computations concerning the GRH*, [arXiv:1305.3087](https://arxiv.org/abs/1305.3087); Math. Comp. 86 (2017) no. 307 per the AMS URL, volume not reconciled [A]. GRH verified for primitive χ mod q ≤ 400 000 up to height max(10⁸/q, A·10⁷/q + 200).

**Records for ζ.**

| Who | Certified support 2a | Status |
|---|---|---|
| Zhu, arXiv:2608.24827v2 | 1.6 (8.9e-18); the 2.38 claim was retracted | Preprint; audited in `docs/AUDIT_ZHU.md` |
| Liu, alphaxiv 2609 | 2.125, with constant 2⁻⁴⁹¹⁶² | Preprint [A] |
| This work | 2.38 (6.81e-48, audited); 2.6 (5.76e-62, commit 108d678); 2.99 running | Internal |

I found nothing beyond 2.125 in the literature. The decay laws are in the priority section.

**Bearing on this work.**
- **The Carneiro-type extremal problems use our function class but assume RH.** Their test functions have Fourier support in [−Δ, Δ], exactly our class with 2a = 2πΔ. With RH assumed, their zeros side is a sum of |F|²: they optimise bounds and never certify positivity. Their extremal majorants could serve as structured trial vectors for upper bounds.
- **Our certificates are close to the true minimum; the published ones are not.** Our lower bounds sit within a factor of 2 of finite-basis upper bounds. Liu's are about 10¹⁶ (at a = 1) or more below.

## E. Spin structures and theta characteristics in arithmetic

**Sources.**
- **Atiyah**, *Riemann surfaces and spin structures*, Ann. Sci. ÉNS 4 (1971) 47–62, [numdam](http://www.numdam.org/item/ASENS_1971_4_4_1_47_0/) [A]. **Mumford**, *Theta characteristics of an algebraic curve*, Ann. Sci. ÉNS 4 (1971) 181–192 [A].
  - Theta characteristics (L² ≅ K) are spin structures.
  - The parity h⁰(L) mod 2 is deformation-invariant.
  - **Farkas**, *Theta characteristics and their moduli*, [arXiv:1201.2557](https://arxiv.org/abs/1201.2557) [R, §1]: 2^{g−1}(2^g + 1) even and 2^{g−1}(2^g − 1) odd, matching the Arf invariant of quadratic forms on J[2].
- **Hecke's theorem**: the class of the different is a square in the class group, an arithmetic theta characteristic. **Armitage**, *On a theorem of Hecke in number fields and function fields*, Invent. Math. 2 (1967) 238–246 [A].
- **Carmeli–Shusterman–Zehavi**, *Arithmetic Wu formulas and the generalized Hecke theorem*, [arXiv:2606.06008](https://arxiv.org/abs/2606.06008) [A]. Steenrod squares on étale cohomology of arithmetic schemes, a finite-field analogue of Atiyah's theta-characteristic theorem, and Serre's Riemann–Hurwitz for spin bundles. No L-functions or positivity in the abstract.
- **Deligne**, *Les constantes locales de l'équation fonctionnelle de la fonction L d'Artin d'une représentation orthogonale*, Invent. Math. 35 (1976) 299–316, §1 [R].
  - w²(V) is the class of the central extension induced by **the double cover Spin(V, Q) → SO(V, Q)** (1.3).
  - Theorem (1.5): for V real, virtual, of dimension 0 and determinant 1, W(V) = exp(2πi cl(w²(V))). So W = ±1 according to whether w² is trivial.
  - (1.6) recovers **Fröhlich–Queyrut**, Invent. Math. 20 (1973): the global root number of an orthogonal representation is +1.
- **Armitage**, *Zeta functions with a zero at s = ½*, Invent. Math. 15 (1972) 199–207 [A, plus search summary]. There are Dedekind ζ_K with ζ_K(½) = 0: Serre's degree-48 field, and Fröhlich's quaternion fields. The zero is forced by symplectic Artin factors with W = −1.
- **Spin L-functions.**
  - **Andrianov**, *Euler products corresponding to Siegel modular forms of genus 2*, Russian Math. Surveys 29 (1974) 45–116 [A].
  - **Furusawa–Morimoto**, *On the Gross–Prasad conjecture with its refinement for (SO(5), SO(2)) and the generalized Böcherer conjecture*, Compositio Math. 160 (2024) 2115–2202 [A]: weighted averages of Fourier coefficients of degree-2 Siegel eigenforms relate to central values of twisted spinor L-functions through squared-modulus periods, which forces nonnegativity.

**Bearing on this work.**
- **No source found uses a theta characteristic or arithmetic spin structure in a positivity argument for RH or the Weil conjectures.** Weil's positivity (Castelnuovo, Rosati, Hodge index) uses a polarisation, not a square root of K.
- **Where spin does enter arithmetic, it controls signs, not positivity:**
  - Deligne's w² gives root numbers.
  - Fröhlich–Queyrut gives W = +1 for orthogonal representations. The quadratic χ in ζ_K is orthogonal and one-dimensional, so its central order is even.
  - Symplectic pieces can force ζ_K(½) = 0 (Armitage).
- **Spin-½ positivity relevant to this work is metaplectic, at the archimedean place.**
  - CCM's prolate operator has the half-integer Hermite spectrum n + ½.
  - Our sectors are Fourier classes n ≡ 0 and n ≡ 2 (mod 4) inside the weight-½ piece.
  - [D] Imaginary quadratic K uses both metaplectic pieces, since Γ_ℂ(s) = Γ_ℝ(s)Γ_ℝ(s+1): an even and an odd archimedean factor.
- **Yoshida's Prop. 1 [S: Suzuki] is the natural "antisymmetric sector" theorem.** Odd-sector positivity for all a already implies RH. The even sector cannot see real zeros. [D] A real zero pair at ½ ± δ adds +2(∫f cosh δu)² to the even form and −2(∫f sinh δu)² to the odd form, like our pole term with δ in place of ½ and the opposite sign.

---

## Candidate fail-fast tests for this work

Each test has a pre-registered prediction and a kill criterion. "λ" means the window minimum in a stated sector, certified lower bound or finite-basis upper bound as stated.

### Test 1. Spinor-floor dominance on the quadratic tower (reading 1)

**Basis.** [D: C facts 2–3] from Weil 1952 additivity and the Connes–Consani count of about 2x/q near-null prolates.

**Run.** Compute λ_even and λ_odd at x = e^{2.2}, e^{2.6}, e^{2.99} for:
- ζ;
- L(χ₋₂₀), with Γ_ℝ(s+1) and conductor 20;
- ζ_K;
- L(χ₋₄) and L(χ₅);
- L_K(ψ) = L(χ₋₄)L(χ₅);
- ζ_H = ζ·L(χ₋₂₀)·L(χ₋₄)·L(χ₅).

**Prediction.**
- (a) Q_{ζ_K} = Q_ζ + Q_{L(χ₋₂₀)} and Q_{ζ_H} = Q_{ζ_K} + Q_{L_K(ψ)} as matrices, to 10⁻³⁰ relative.
- (b) Superadditivity, a theorem.
- (c) **λ(ζ_K)/λ(L(χ₋₂₀)) ∈ [1, 3] and λ(ζ_H)/λ(ζ_K) ∈ [1, 1.5]** in both sectors. The cover's margin is the largest-conductor factor's floor. At e^{2.99} this means λ_even(L(χ₋₂₀)) ≥ 1.6·10⁻³ (certified ζ_K even = 4.749·10⁻³).

**Kill.**
- (a) or (b) fails: normalisation bug; stop.
- λ(ζ_K) > 10·λ(L(χ₋₂₀)), or λ(ζ_H) > 3λ(ζ_K): the cover shows positivity beyond its factors' floors. That kills this mundane explanation and revives reading (1) as a real mechanism.
- If (c) holds, reading (1) is reduced to archimedean conductor bookkeeping, and the ζ_K-versus-Z₁ separation comes only from Z₁ being one sheet (C fact 2).

**Result (main session, 2026-10-01; `decay_law_mp.py tower`, `data/connes/decay/tower.json`, `tower_control.txt`).** x = e^{2.2}, e^{2.6}, e^{2.99}; N = 48 and 96 agree to about 1%.
- **(a) passes.** ζ + L(χ₋₂₀) matches the independent control code's ζ_K to 1e-60.
- **(b) passes.** Superadditivity holds everywhere.
- **(c), first half, passes.** λ(ζ_K)/λ(L(χ₋₂₀)) = 1.03, 1.09, 1.28 (even) and 1.00, 1.03, 1.07 (odd).
- **(c), second half, fails.** λ(ζ_H)/λ(ζ_K) = 1.19, 2.48, **9.99** (even) and 1.21, 1.54, **3.02** (odd), so the kill (> 3) fires at e^{2.99}.
- **Control run before accepting the revival of reading (1).** Swap the cover's L(χ₋₄)L(χ₅) for products that are not covers. Boost of λ(ζ_K·extra)/λ(ζ_K) at e^{2.99}, even sector:

  | extra factors | boost |
  |---|---|
  | L(χ₋₄)L(χ₅), the cover | 9.9 |
  | L(χ₋₇)L(χ₁₃) | 23.1 |
  | L(χ₋₄)L(χ₋₃) | 9.5 |
  | L(χ₅)L(χ₈) | 6.6 |
  | L(χ₋₃)L(χ₈) | 2.95 |
  | L(χ₋₄) alone | 7.7 |
  | L(χ₅) alone | 5.1 |

  The odd sector behaves the same way.
- **Verdict.**
  - The kill refutes the "largest-conductor floor" bookkeeping. It does not revive the double cover: being an unramified cover contributes nothing that any other added L-function doesn't.
  - Every zeros-side form is Σ_ρ |F(γ)|² over the union of the factors' zeros. Adding a factor adds zeros that the minimiser must also avoid. The boost depends on where the added factor's low zeros fall relative to ζ_K's minimiser, not on Galois structure.
  - **Reading (1) stays dead.**

### Test 2. Prolate-index law at support 2.99 (decay section)

**Prediction at x = e^{2.99} = 19.8857.** ℓ₄ = 4.907·10⁻⁹⁷ and ℓ₆ = 1.634·10⁻⁹², from the Fuchs–Slepian asymptotics.
- **λ_even ∈ [2.9, 5.9]·10⁻⁹⁶**, i.e. ratio ∈ [6, 12] times ℓ₄.
- **λ_odd ∈ [5.7, 13]·10⁻⁹²**, i.e. ratio ∈ [3.5, 8] times ℓ₆.
- odd/even ∈ [10⁴, 4·10⁴].
- Zhu's form with C = 2π² gives 2.26·10⁻⁹⁵ for the even sector.
- Use the exact 1 − χ from `scripts/prolate_ratio_mp.py` as well; the asymptotic is about 5% off at this x.

**Kill.**
- Even lower bound > 6·10⁻⁹⁶ or even upper bound < 2.9·10⁻⁹⁶: the h₄ law is broken. If the bracket contains 2.26·10⁻⁹⁵, Zhu's law is favoured.
- Odd bracket outside [3, 25]·10⁻⁹²: the odd-equals-h₆ assignment is broken.
- Extension: if the fit at x ∈ [20, 40] gives an α drifting above 4π by more than 1%, the prolate law fails asymptotically.

**Reviewer note (main session, 2026-10-01).** Existing data already settles part of this test:
- **Even sector.** The N = 180 upper bound is 6.27·10⁻⁹⁶, and N = 140 gives 8.27·10⁻⁹⁶ (`data/connes/zhu_window_a1495_n180.run.txt`; also `data/connes/antiperiodic/sweep_a1.495.json` at N = 140). The true minimum is therefore below the literal leading term of Zhu's form, 2.26·10⁻⁹⁵. That is a 0.59% difference in −ln λ. Zhu's Conj. 12.1 carries a (1 + o(1)) factor, so this is a disagreement with the leading-term prediction at one finite x, **not** an exclusion of the conjecture (corrected after external review). The certified lower bound, 3.50114217641·10⁻⁹⁶, lies inside the bracket [2.9, 5.9]·10⁻⁹⁶.
- **Odd sector.** The N = 180 upper bound is 1.48·10⁻⁹¹ (N = 140 gives 1.79·10⁻⁹¹; `data/connes/zhu_window_a1495_n180.run.txt`). It is inside the kill window [3, 25]·10⁻⁹² and above the prediction bracket [5.7, 13]·10⁻⁹². It is not converged. **Certified lower bound (2026-10-02): 8.25626494001·10⁻⁹², inside the prediction bracket.** Test 2 holds in both sectors.
- **Test 3** was run independently before this review (`docs/ANTIPERIODIC.md`). The antiperiodic basis is only a basis: same infimum, ≈ 1.2 modes of extra reach.
- **The prolate-index assignment was extended to conductors** (`docs/DECAY_LAW.md`, "Prolate index"). For every one of the 12 real-character cases, λ/ℓ_n(2πx/q) is flattest at n = κ + 2·(sector) + 4·(pole).

### Test 3. Antiperiodic basis is only a basis (reading 3 guardrail)

**Basis.** Yoshida's K(a) uses 2a-periodic restrictions [S: Suzuki]. Both {e^{2πiku/L}} and {e^{i(2k+1)πu/L}} are complete orthonormal bases of L²[−L/2, L/2]. Finite minima decrease to the full-space infimum (Connes–Consani, Cor. 2.4).

**Prediction.**
- The antiperiodic and periodic truncations converge to the same λ* in each sector, and both stay above the certified lower bound.
- [D, Fourier smoothness] Periodic converges faster in the even sector: an even minimiser with v(a) = v(−a) has a continuous periodic extension. Antiperiodic converges faster in the odd sector: v(a) = −v(−a) has a continuous antiperiodic extension.

**Kill.**
- An antiperiodic value below a certified lower bound, or convergence to a different limit: bug.
- If the speed pattern is reversed, the minimiser's boundary behaviour differs from assumed. Diagnostic only.
- Any "spin effect" seen in the antiperiodic basis must survive this test before it is interpreted.

### Test 4. Real zeros are odd-sector only (Yoshida Prop. 1; Stark; Conrey–Soundararajan)

**Run.** Take our matrices for a GRH-verified spinor factor, L(χ₋₄) or L(χ₋₂₀), at a ∈ {0.5, 0.8, 1.0}. Add a planted real zero pair at ½ ± δ, δ ∈ {0.1, 0.3}: the rank-one terms +2ccᵀ (even), c_k = ∫b_k cosh δu, and −2ssᵀ (odd), s_k = ∫b_k sinh δu.

**Prediction.**
- λ_even does not decrease (PSD update).
- λ_odd decreases and turns negative exactly when 2sᵀM_odd⁻¹s ≥ 1 (Sherman–Morrison), first at the smallest such a.

**Kill.** Any decrease of λ_even, or odd failure at a different a than Sherman–Morrison predicts: the parity and pole-sign bookkeeping is wrong, and every odd/even certificate's interpretation must be rechecked.

**Readout.** Landau–Siegel risk lives in the spinor factor (Stark 1974; Conrey–Soundararajan), and the window form sees it only in the odd sector.

### Test 5. Function-field sheet test (C fact 1; exact)

**Run.** Use `scripts/closed_seam_mp.py` on the zeros-side windowed form of the numerator 1 − qu + qu² of the principal-class partial zeta Z_O = (1 − qu + qu²)/(1 − qu), for q ∈ {2, 3, 4, 5, 7, 13}, alongside ζ_E for curves with those q.

**Prediction.**
- On the first two-point window ℓ/2 < a ≤ ℓ, ℓ = log q, Z_O has λ = 2ℓ(1 − √q/2) exactly.
  - Positive for q = 2 and 3, where the zeros lie on |u| = q^{−1/2}, so all windows are PSD.
  - Zero for q = 4.
  - Negative for q ≥ 5, with failing eigenvector (1, −1).
- ζ_E keeps the Hasse margin 2ℓ(1 − |a_q|/(2√q)).

**Kill.** Any deviation means the closed-seam code or the Riemann–Roch derivation is wrong.

**Extension.** A genus-2 curve with #Jac(𝔽_q) = 2 and an étale double cover, whose Prym L-polynomial is nontrivial, gives the full analogue of (ζ_K, L_K(ψ), Z₁, ζ_H) with exact finite Toeplitz positivity.
