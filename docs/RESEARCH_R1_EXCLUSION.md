# R1: What a Certified Finite Window Says About Where ζ's Zeros Cannot Be

Date: 2026-10-02. Worktree branch `worktree-agent-a1f509ee3e7577f26`, base `b4a055d` (`feat/weil-gram-instrument`). Brief: `docs/BRIEF_RESEARCH_HANDOFF.md` §5, question R1. Nothing outside the worktree was written; no existing script or document was modified. New code lives in `scripts/research_r1/`, data in `data/research_r1/`. Every run used at most four processes and logged through `scripts/progress.py`.

**Labels.** [R] read in the source; [A] abstract only; [S] second-hand (from another paper's or an rh2 document's account); [D] derived here; [N] computed here (multiprecision floating point, finite bases: upper bounds and exact finite-dimensional identities, not certificates); [U] unverified.

## 0. Five lines

1. **The term.** A quadruple ρ = ½ ± δ ± iγ₀ enters Weil's form for a real test function f as T(f) = 4(|f̂_c(γ₀)|² − |f̂_s(γ₀)|²), with f_c = f cosh δu and f_s = f sinh δu: a rank-2 form of signature (+, −) in rh2's bases [D], checked against direct quadrature of Σ F(ρ)F(1−ρ) to 10⁻³⁰ in both sectors and against the conductor-5 document's +17.12 [N].
2. **The logic.** "The pair's most-negative value exceeds c, so the pair is excluded" is wrong: the other zeros compensate. DH is certified positive at support 2.38 with c = 1.04·10⁻⁸ while its own pair's most-negative Rayleigh quotient there is −0.20. The only hypothesis-honest reverse statement is a consistency test: under **H(δ, γ₀)** = "that quadruple, and every other zero on the critical line", Q(f) ≥ T(f) must hold for every f supported in [−a, a], i.e. **R₄ := 4cᵀ(Q_a + 4ssᵀ)⁻¹c ≤ 1**. The region excluded at support a is R(a) = {R₄ > 1}. It depends on the form Q_a, not on its certified minimum c, and c cannot bound it in either direction [D].
3. **Calibration.** On F_{t*}, DH and Z₁ the decomposition Q = Σ_pairs T + P with the census zeros gives P ⪰ 0 at every support tested, P's floor equals the Euler partner's (F_{t*} minus its real pair reproduces L(χ₅)'s floor to two digits at every x), and the mirror test σ₁ = 4sᵀ(P + 4ccᵀ)⁻¹s = 1 lands on all three certified crossings: σ₁ = 1.00003 at x = 5.3625, 1.0002 at 30.745, 0.998 at 19.844 [N]. The naive criterion predicts failure at x < 10 for DH and x ≈ 8–9 for Z₁ (killed); extrapolating ln σ₁ from the certified-positive regime predicts the crossings to a few per cent in x when that regime reaches σ₁ ≳ 0.3 (F_{t*}, Z₁) and 30% early when it stops at σ₁ ≈ 0.08 (DH).
4. **The regions for ζ.** R(a) is essentially independent of δ (a window of half-length a resolves the real part of a zero only to ~1/a; cosh δu ∈ [1, 1.3] over the window): it is a statement about zero *multiplicity and local count*, not about the line. Below T* = 2πx the window pins the zeros' positions (2κ(γ_k) = 1 at the true zeros, to 10⁻¹⁰ at a = 0.8 and 10⁻⁴⁷ at a = 1.495); above T* a quadruple is excluded only while the window's own local zero count near γ₀, read off the primes n < x, is below four zeros' worth, which ends at a height of order 2πe^{4a}·e^{ρ} with ρ the prime-phase fluctuation: 435 at support 1.6 (matrix), excursions seen to 10⁴ and an envelope of 5·10⁵ (δ = 0) to 2·10⁶ (δ = ½) at 2.99, against RH verified to 3·10¹² [A]. **No certified support excludes a point that is not already excluded.**
5. **Go/no-go: no-go.** Supports 3.2, 3.5 and 4.0 move the typical reach to 4·10⁴, 8·10⁴, 2.5·10⁵ (δ = 0) and 2·10⁵, 8·10⁵, 10⁷ (δ = ½) in height under H, and even the envelope (every prime phase aligned, shaping gain included) tops out at 7·10⁹ at support 4.0: still below 3·10¹². And whatever such a window excludes is a local zero-count bound that follows from one prime sum over n < x at that height, needs no positivity certificate, assumes RH for every other zero, and says nothing about the line. The certificates' value lies elsewhere (the decay law, R2–R4), not in exclusion.

## 1. The quadruple's term [D]

**Setting.** f real, supported in [−a, a] (rh2's window [−L/2, L/2], L = 2a = log x). F(s) := ∫ f(u) e^{(s−½)u} du, so F(½ + it) = F̂(t) := ∫ f e^{itu} du, and F(1 − s) = ∫ f(−u) e^{(s−½)u} du. The explicit formula gives Weil's form as a sum over all non-trivial zeros, without any hypothesis on their position [S: standard; `docs/CONTROL_CERTIFICATES.md` §1.1 states it for the controls; `docs/AUDIT_ZHU.md` §1 for ζ]:

  Q(f) = Σ_ρ F(ρ) F(1 − ρ).

Each term is the autocorrelation transform ĝ at γ_ρ = −i(ρ − ½), g = f ⋆ f̃, and for ρ = ½ + iγ on the line it is |F̂(γ)|². Pairing γ with −γ gives rh2's 2Σ_{γ>0}|F̂(γ)|², which is therefore the zeros side **only under RH**.

**Proposition 1 (the quadruple).** Let ρ = ½ + δ + iγ₀ with δ ≠ 0 and γ₀ ≠ 0. Its orbit under s ↦ 1 − s and s ↦ s̄ has four points, and their total contribution to Q(f) is

  T(f; δ, γ₀) = 4 Re[ F̂(γ₀ − iδ) F̂(−γ₀ + iδ) ] = 4 ( |f̂_c(γ₀)|² − |f̂_s(γ₀)|² ),  f_c := f cosh δu, f_s := f sinh δu.

*Proof.* F(ρ) = ∫ f e^{(δ+iγ₀)u} = F̂(γ₀ − iδ) and F(1 − ρ) = F̂(−γ₀ + iδ); the ρ̄ term is the complex conjugate of the ρ term and the 1−ρ, 1−ρ̄ terms repeat them. Write e^{δu} = cosh δu + sinh δu: F̂(γ₀ − iδ) = f̂_c(γ₀) + f̂_s(γ₀) and F̂(−γ₀ + iδ) = conj(f̂_c(γ₀)) − conj(f̂_s(γ₀)) because f is real. The product's real part is |f̂_c|² − |f̂_s|². ∎

No parity is used. In the two sectors (b_k rh2's orthonormal bases), T = 4(ccᵀ − ssᵀ) with

| sector | c_k | s_k |
|---|---|---|
| even, b₀ = 1/√L, b_k = √(2/L) cos ω_k u | ∫ b_k cosh δu cos γ₀u du | ∫ b_k sinh δu sin γ₀u du |
| odd, b_k = √(2/L) sin ω_k u | ∫ b_k cosh δu sin γ₀u du | ∫ b_k sinh δu cos γ₀u du |

all in closed form (`scripts/research_r1/r1lib.py`). The matrix has rank 2 and signature (+1, −1): the positive direction is the sector's natural transform at γ₀ with the cosh weight, the negative one the "cross" transform with the sinh weight.

**Degenerate cases** (multiplicity m = 2 instead of 4): a real pair (γ₀ = 0) gives 2(ccᵀ − ssᵀ) with c_k = ∫ b_k cosh δu, s_k = ∫ b_k sinh δu, i.e. **+2(∫ f cosh δu)² for even f and −2(∫ f sinh δu)² for odd f**, the terms already in `docs/CONDUCTOR5_FAMILY.md` and `docs/LIT_TANGENTS.md`; an on-line pair (δ = 0) gives 2|F̂(γ₀)|² ≥ 0, with s = 0.

**Sign structure and size.** Expanding in δ at fixed f,

  T = 4|F̂(γ₀)|² − 4δ²( |F̂′(γ₀)|² + Re F̂(γ₀) conj F̂″(γ₀) ) + O(δ⁴).

So for a test function whose transform vanishes at γ₀, as every near-minimiser's does at the zeros below T*, the quadruple contributes **−4δ²|F̂′(γ₀)|²**: negative, quadratic in δ, and of polynomial, not exponential, size. For a test function resonant with γ₀ (f ≈ cos γ₀u or sin γ₀u on the window) T ≈ +4|f̂_c(γ₀)|² is positive and O(a²). The most negative value of T alone, min_f T/‖f‖² = −2[(|s|²−|c|²) + √((|c|²−|s|²)² + 4(|c|²|s|² − (c·s)²))], is O(δ²a³) and is reached by the π/2-phase-shifted resonance (f ≈ sin(γ₀|u|)·w(u) in the even sector). Table 3.1 lists it for the controls' pairs.

**Checks [N]** (`scripts/research_r1/check_T.py`). For random trigonometric f in both sectors and (δ, γ₀) ∈ {(0.3085, 85.699), (¼, 0), (0, 14.1347), (0.1, 3.5), (0.45, 1.2)}, the closed form agrees with direct quadrature of Σ F(ρ)F(1−ρ) over the orbit to a relative 8·10⁻³¹. For F_{t*}'s real pair and f = (1 + cos θ)² the even-sector term is +17.1204 at x = 7 and +11.6941 at x = 5, the +17.12 and +11.69 of `docs/CONDUCTOR5_FAMILY.md`.

## 2. What "excluded by a certified window" can mean [D]

### 2.1 The naive reading, and why it is wrong

The brief's framing ("the maximum of the negative term relative to ‖f‖², against c, defines R(a, c)") does not survive contact with the controls. The certificate Q(f) ≥ c‖f‖² bounds the **whole** form, including the hypothetical quadruple and every other zero: Q(f) = T(f) + Σ_{other ρ} F(ρ)F(1−ρ). A negative T(f) is compensated by the other zeros' positive terms for that f, and nothing in the certificate limits those. Concretely: DH is certified positive on [−1.19, 1.19] with c = 1.0441·10⁻⁸ (`docs/CONTROL_CERTIFICATES.md`), and DH *has* the quadruple 0.8085 ± 85.699i; its most-negative Rayleigh quotient at that support is −0.197 (Table 3.1). If "min T/‖f‖² < −c" excluded a pair, it would exclude DH's own zero at a support where DH is rigorously positive. The same holds for Z₁ (c = 0.049 at x = 15.8 against −0.74) and, more mildly, F_{t*} (c = 1.17·10⁻³ at 5.312 against −0.0485). The naive criterion is killed three times over.

### 2.2 The consistency theorem

What survives is a statement under an explicit hypothesis on the other zeros.

**Hypothesis H(δ, γ₀).** The non-trivial zeros of ζ consist of the quadruple ½ ± δ ± iγ₀ together with a set S of zeros on the critical line (positions and multiplicities unknown).

**Proposition 2.** Under H(δ, γ₀), for every real f ∈ L² supported in [−a, a],

  Q(f) ≥ T(f; δ, γ₀),

because Σ_{γ∈S} 2|F̂(γ)|² ≥ 0. In a finite orthonormal basis of either sector this is Q_N − 4(ccᵀ − ssᵀ) ⪰ 0, and when Q_N ≻ 0 (which every certified support guarantees) it is equivalent, by Sherman–Morrison applied to the two rank-one pieces, to

  **R₄(a; δ, γ₀) := 4 cᵀ(Q_N + 4ssᵀ)⁻¹c ≤ 1.**

Since Q_N is a compression of Q, R₄ computed in any finite basis is a lower bound for the full-space quantity, so R₄ > 1 in a finite basis already contradicts H. The excluded region at support a is

  R(a) := { (δ, γ₀) : R₄(a; δ, γ₀) > 1 } = { (δ, γ₀) : ∃ f supported in [−a, a] with T(f) > Q(f) }.

*Proof.* Q − T = Σ_S is a sum of positive rank-one forms. The Sherman–Morrison step: Q + 4ssᵀ − 4ccᵀ ⪰ 0 iff 4cᵀ(Q + 4ssᵀ)⁻¹c ≤ 1, and cᵀ(Q + 4ssᵀ)⁻¹c = cᵀQ⁻¹c − 4(cᵀQ⁻¹s)²/(1 + 4sᵀQ⁻¹s). ∎

**Remarks.**
- **The certified constant c plays no role in R(a).** Exclusion needs Q(f) *small* for some f with |f̂_c(γ₀)| large; c says Q(f) is *at least* c, which is the wrong direction, and the only bound it yields, 4‖c‖²/c ≤ 1 ⇒ not excluded, is vacuous (‖c‖² = O(a) against c ≤ 10⁻¹⁷). What exclusion uses is an **upper** bound on Q(f) for an explicit f, i.e. a Rayleigh quotient, which rh2 already certifies in Arb through the u-space formula (`scripts/control_uspace.py witness`). So every point of R(a) can be certified by a witness pair (f, Q(f) ≤ q, T(f) > q) at the cost of one u-space evaluation, without any positivity certificate. The question's "R(a, c)" is therefore R(a): a property of the window's form, not of its minimum.
- **What H buys, exactly.** By Krein's extension theorem [S: standard], a form on [−a, a] that is positive on autocorrelations of functions supported there is the restriction of ∫|F̂(t)|² dμ(t) for some positive measure μ on ℝ. So R₄ ≤ 1 holds iff *some* on-line configuration (a positive measure, not necessarily realisable as the zeros of a function with ζ's functional equation) is consistent with the window and the quadruple. The test is the complete content of the window under H; it cannot be sharpened without more hypotheses.
- **Stronger hypotheses.** If in addition the zeros below some height are *known* (on the line, at known positions), their sum P_known can be subtracted and the test becomes Q − T − P_known ⪰ 0, which excludes more. This is only available below 3·10¹² and is discussed in §4.4.
- **Real zeros.** For a real pair (γ₀ = 0) the odd-sector term is −2ssᵀ ⪯ 0, so Q_odd − T ⪰ Q_odd ⪰ 0 always: **the odd sector can never exclude a real pair by consistency**, although it is the sector that *fails* for one (Yoshida's Prop. 1 [S: Suzuki §1.1]). The even sector's term is +2ccᵀ and the test 2cᵀQ_even⁻¹c ≤ 1 is non-trivial: it is the same mechanism by which the window sees the pole (+2vvᵀ, v the cosh(u/2) vector). For ζ it rejects every real pair massively (§4.3), consistent with ζ having no real zeros in (0, 1) [S: elementary].
- **The forward mirror.** For a function that *has* the quadruple, write Q = P + 4(ccᵀ − ssᵀ) with P ⪰ 0 the rest. Then Q ⪰ 0 iff σ₁ := 4sᵀ(P + 4ccᵀ)⁻¹s ≤ 1, and the first failure is σ₁ = 1. This is the control statement (§3), and it is exact given P.

### 2.3 Why δ hardly enters

Over the window |u| ≤ a ≤ 1.5 and for |δ| ≤ ½, cosh δu ∈ [1, cosh(a/2)] ⊂ [1, 1.30] and |sinh δu| ≤ 0.82. The cosh-weighted vector c(δ, γ₀) therefore differs from the on-line vector c(0, γ₀) by a bounded multiplicative weight, and the sinh-vector s is a bounded perturbation. R₄ at (δ, γ₀) is within a factor cosh²(δa) ≤ 1.7 of R₄ at (0, γ₀), where the test asks whether a **double** on-line zero at γ₀ is consistent with the window; for the box test function of §4 the factor is exactly w(δa) := (sinh δa/δa)², which is 1.05 at a = 0.8 and 1.20 at a = 1.495 for δ = ½. In the computations below the δ-dependence of the boundary of R(a) is at the 10⁻³ level in the pinned regime and a factor ≤ 1.2 in R₄ above it. To resolve the real part of a zero at height γ₀ the window would need e^{δu} to vary appreciably across it, δa ≳ 1, i.e. supports of order 1/δ: a support 2a = 3 cannot tell δ = 0.05 from δ = 0 at any height. **R(a) is a statement about the multiplicity and local density of zeros, not about the critical line.** This is also why the controls fail through the sinh channel (−4δ²|F̂′(γ₀)|², a *derivative* of the transform at the pair's height) and only once the window's near-null space has leaked enough at γ₀ (§3.4).

## 3. Calibration on the controls [N]

Scripts: `scripts/research_r1/controls.py`, data `data/research_r1/ctrl_*.json`. For each control, sector and x, rh2's zeros-side matrix Q (`control_normalisation.rh2_matrix`, factor 1 against the certified Q by `docs/CONTROL_CERTIFICATES.md` §2) is split as Q = Σ_i T_i + P, with T_i the terms of Proposition 1 for the census off-line zeros (DH: the five pairs of `docs/DAVENPORT_HEILBRONN.md`; Z₁: the seven of `docs/EPSTEIN.md`; F_{t*}: the real pair at ¾, ¼). P is "the on-line part plus the pole". Three quantities are tabulated: λ_min(Q) (negative = a witness), λ_min(P) (must be ≥ 0 if the census is complete to the heights the window sees and Proposition 1 is right), and σ₁ = m sᵀ(P₁ + m ccᵀ)⁻¹s for the first pair, with P₁ = Q − T₁ (so the other pairs sit inside P₁). "naive" is min_f T₁(f)/‖f‖², the quantity the brief's framing would compare with c. Bases: F_{t*} N = 64 (odd), 48 (even); DH N = 80; Z₁ N = 64; 40–70 digits.

### 3.1 Tables

**F_{t*}, odd sector** (real pair δ = ¼, m = 2: T = −2ssᵀ; certified bracket for the first failure [5.31217, 5.36250)):

| x | λ_min(Q) | λ_min(P) = λ_min(Q + 2ssᵀ) | σ₁ = 2sᵀP⁻¹s | naive |
|---|---|---|---|---|
| 2 | 1.685 | 1.688 | 0.0019 | −0.0034 |
| 3 | 0.7657 | 0.7777 | 0.0160 | −0.0137 |
| 4 | 0.1950 | 0.2201 | 0.115 | −0.0277 |
| 4.5 | 0.0866 | 0.1165 | 0.258 | −0.0354 |
| 5 | 0.02730 | 0.06173 | 0.559 | −0.0434 |
| 5.2 | 0.01093 | 0.04705 | 0.769 | −0.0466 |
| 5.31217 (x_pos) | 3.205·10⁻³ | 0.04024 | 0.9208 | −0.0485 |
| 5.36250 (x_neg) | −1.010·10⁻⁶ | 0.03743 | **1.00003** | −0.0493 |
| 6 | −0.03031 | 0.01200 | 3.50 | −0.0599 |
| 7 | −0.04543 | 1.466·10⁻³ | 31.4 | −0.0769 |

At x_pos the certified c = 1.1666·10⁻³ lies below λ_64 = 3.205·10⁻³ as it must; at x = 7 the witness −0.045280 is reproduced. **P is positive throughout**, and its floor is the Euler-class floor: 1.47·10⁻³ at x = 7 against L(χ₅)'s odd 1.4·10⁻³ (`docs/CONDUCTOR5_FAMILY.md`, t = 0 row).

**F_{t*}, even sector** (the pair adds +2ccᵀ; the sector never fails): λ_min(Q) = 1.372, 0.5689, 0.08133, 4.006·10⁻³, 6.175·10⁻⁴ at x = 3, 5, 7, 8.846, 10; λ_min(P) = λ_min(Q − 2ccᵀ) = 0.04587, 7.167·10⁻⁴, 5.372·10⁻⁶, 7.863·10⁻⁸, 3.736·10⁻⁹. Compare L(χ₅)'s even minima 0.044, 6.9·10⁻⁴, 5.0·10⁻⁶, —, 3.8·10⁻⁹ at x = 3, 5, 7, 10: **F_{t*} with its real pair removed has L(χ₅)'s floor to two digits at every x.** The real pair is the whole difference between the two functions' window minima, in both sectors.

**DH, even sector** (five quadruples; the first at δ = 0.3085, γ₀ = 85.699; certified bracket [10.805, 30.745)):

| x | λ_min(Q) | λ_min(P_all) | λ_min(P₁) | σ₁ | naive (pair 1) |
|---|---|---|---|---|---|
| 10 | 8.235·10⁻⁸ | 8.112·10⁻⁸ | 8.254·10⁻⁸ | 0.0696 | −0.197 |
| 13 | 8.375·10⁻¹¹ | 7.397·10⁻¹¹ | 7.774·10⁻¹¹ | 0.1077 | −0.276 |
| 16 | 3.795·10⁻¹⁴ | 3.780·10⁻¹⁴ | 3.820·10⁻¹⁴ | 0.1678 | −0.354 |
| 20 | 4.359·10⁻¹⁸ | 3.431·10⁻¹⁸ | 3.839·10⁻¹⁸ | 0.2548 | −0.448 |
| 25 | 1.483·10⁻²³ | 1.702·10⁻²³ | 1.751·10⁻²³ | 0.511 | −0.558 |
| 28 | 7.609·10⁻²⁷ | 1.077·10⁻²⁶ | 1.135·10⁻²⁶ | 0.730 | −0.617 |
| 30 | 1.137·10⁻²⁸ | 3.450·10⁻²⁹ | 4.321·10⁻²⁹ | 0.844 | −0.655 |
| 30.745 (x_neg) | −2.399·10⁻³² | 4.340·10⁻³⁰ | 4.316·10⁻³⁰ | **1.0002** | −0.681 |
| 31.25 | −7.773·10⁻³⁰ | 8.835·10⁻³¹ | 9.105·10⁻³¹ | 1.237 | −0.685 |
| 32 | −1.3999·10⁻²⁹ | 8.004·10⁻³² | 1.111·10⁻³¹ | 1.874 | −0.691 |
| 33 | −1.803·10⁻²³ | 3.396·10⁻³³ | 5.484·10⁻³³ | 4.35 | −0.725 |

The certified witnesses are reproduced (−1.39985·10⁻²⁹ at x = 32 to all printed digits; −2.29·10⁻³² at 30.74497 against −2.40·10⁻³² at 30.745). **P_all stays positive through and beyond the failure**, decaying like the Euler class (3.45·10⁻²⁹ at x = 30, at the lower edge of the odd-character band e^{−4πx/q + 10.2…11.5} of `docs/DECAY_LAW.md`), while Q's negative eigenvalue grows by six decades per unit of x after the crossing.

**DH, odd sector:** σ₁ = 0.0697, 0.108, 0.161, 0.247, 0.448, 0.642, 0.834, 0.868, 0.961, **1.241**, 2.29 at the same x; λ_min(Q) turns negative between 31.25 (+1.56·10⁻²⁶) and 32 (−1.57·10⁻²⁶), where `docs/DECAY_LAW.md` had the sign change in (31.25, 32.5); λ_min(P_all) > 0 throughout (1.49·10⁻²⁷ at 32).

**Z₁, even sector** (seven quadruples, the first at δ = 0.4330, γ₀ = 15.668; certified bracket [15.800, 19.844)):

| x | λ_min(Q) | λ_min(P_all) | λ_min(P₁) | σ₁ | naive (pair 1) |
|---|---|---|---|---|---|
| 5 | 1.846 | 1.591 | 1.655 | 0.0314 | −0.130 |
| 7 | 1.201 | 0.8027 | 0.8063 | 0.0657 | −0.256 |
| 10 | 0.5191 | 0.2859 | 0.2885 | 0.146 | −0.435 |
| 13 | 0.1374 | 0.1240 | 0.1404 | 0.227 | −0.535 |
| 15 | 0.07435 | 0.03254 | 0.03850 | 0.304 | −0.714 |
| 15.8 (x_pos) | 0.05077 | 0.02204 | 0.02367 | 0.376 | −0.737 |
| 17.5 | 0.01556 | 0.01170 | 0.01294 | 0.560 | −0.746 |
| 19 | 2.159·10⁻³ | 6.166·10⁻³ | 7.000·10⁻³ | 0.820 | −0.817 |
| 19.844 (x_neg) | 1.43·10⁻⁵ (N = 64) | 3.807·10⁻³ | 4.136·10⁻³ | **0.998** | −0.889 |
| 20 | −1.582·10⁻⁴ | 3.516·10⁻³ | 3.743·10⁻³ | 1.025 | −0.904 |
| 22 | −5.445·10⁻⁴ | 1.130·10⁻³ | 1.475·10⁻³ | 1.113 | −1.074 |

At x = 20 the N = 64 value −1.58·10⁻⁴ is `docs/EPSTEIN.md`'s; at 19.844 the N = 80 witness is negative while N = 64 is still +1.4·10⁻⁵, as a finite-basis crossover must behave. P_all's floor at x = 20, 3.5·10⁻³, sits beside the Euler partner ζ_K's 4.6·10⁻³ (N = 48). **Z₁, odd sector:** σ₁ = 0.028, 0.047, 0.099, 0.247, 0.279, 0.285, 0.351, 0.460, 0.553, 0.576, **1.124** at the same x; λ_min(Q) = 0.1717 at x = 20 and −0.0148 at 22, so the odd sector fails in (20, 22), after the even one; λ_min(P_all) ≥ 0.060 throughout.

### 3.2 The ladder

**Rung 1, the naive criterion (killed).** Reading "min_f T(f)/‖f‖² more negative than the certified c (or than the on-line floor)" as a prediction of failure: F_{t*} at x = 5.2 (−0.0466 against the floor 0.047), i.e. 5.2–5.3 against the measured 5.31–5.36, by luck of a real, low-frequency pair whose dip direction nearly coincides with the floor's minimiser; **DH at x < 10** (−0.20 against 8·10⁻⁸) against 30.745; **Z₁ at x ≈ 8–9** (−0.44 at x = 10 against 0.29) against 19.84. For the complex pairs the criterion is wrong by a factor 2–3 in x and by twenty orders of magnitude in the quantity compared. The reason is §2.1: the test function that makes T most negative has a large on-line part P(f), and the floor is reached only by functions whose transform vanishes at the on-line zeros below T*, which have small |F̂′(γ₀)|.

**Rung 2, the exact decomposition (passes).** P ⪰ 0 at every one of the 58 (control, sector, x) points, including beyond the failures; P's floor equals the Euler partner's where one exists; and σ₁ = 1 reproduces the three certified crossings to 3·10⁻⁵ (F_{t*}), 2·10⁻⁴ (DH) and 2·10⁻³ (Z₁) in σ₁, which in x is better than the width of each bracket. This is a check of Proposition 1's factor 4 and sign structure in three settings: a real pair in the odd sector, complex pairs in both sectors of a degree-1 function, and complex pairs in degree 2. (The identity λ_min(Q) = 0 ⇔ σ₁ = 1 is algebra; what is tested is that the P so defined is positive and Euler-like, which fails if T is mis-scaled or mis-signed.)

**Rung 3, prediction from the positive regime.** Fit ln σ₁ = αx + β to the grid points with x ≤ x_pos (the certified-positive support) and solve for σ₁ = 1 (`scripts/research_r1/fit_sigma.py`):

| control, sector | points used (all with x ≤ x_pos) | σ₁ at the last point | least-squares line → x_c | last two points → x_c | measured first failure |
|---|---|---|---|---|---|
| F_{t*} odd | 7 (2 ≤ x ≤ 5.312) | 0.92 | 5.292 | 5.364 | [5.31217, 5.36250) |
| DH even | 4 (6 ≤ x ≤ 10.805) | 0.081 | 21.1 | 24.1 | (30.0, 30.745) in this grid; certified < 30.745 |
| DH odd | 4 | 0.080 | 21.6 | 25.6 | (31.25, 32.0) |
| Z₁ even | 5 (5 ≤ x ≤ 15.8) | 0.38 | 19.73 | 23.1 | (19.844, 20.0) at N = 64; certified < 19.844 |
| Z₁ odd | 5 | 0.29 | 19.6 | 35.5 | (20, 22) |

ln σ₁ is not linear in x: concave far from the crossing, convex close to it (DH's slope per unit x goes 0.30 → 0.24 → 0.15 → 0.11 → 0.14 → 0.12 → 0.07 between 6 and 30, then 0.23–0.84 after the crossing). The least-squares line through the certified-positive regime predicts the crossing to 0.4% when that regime reaches σ₁ ≈ 0.9 (F_{t*}), to 1–2% when it reaches σ₁ ≈ 0.3–0.4 (Z₁ even; Z₁ odd's 19.6 against the (20, 22) grid cell), and is **31% early** when it stops at σ₁ ≈ 0.08 (DH, both sectors). Using DH's numerically positive points up to x = 20 (σ₁ = 0.25) the line gives 27.6 (10% early), and the three points 13–20 alone give 31.1 (1% late). **The protocol is good to a few per cent in x when the data reach σ₁ ≳ 0.3, and to ~30% from σ₁ ≈ 0.1**; the planted GL(2) cases of §6 confirm the second number (17–22% from σ₁ ≤ 0.11). The last-two-points rule is always worse.

### 3.3 DH below the certified support

Three more DH points inside the certified-positive range (`data/research_r1/ctrl_dh_*_low.json`, N = 80), so that rung 3 could be run from x ≤ x_pos alone:

| x | even: λ_min(Q), λ_min(P_all), σ₁ | odd: λ_min(Q), λ_min(P_all), σ₁ |
|---|---|---|
| 6 | 7.205·10⁻⁴, 7.073·10⁻⁴, 0.0253 | 0.1107, 0.1091, 0.0262 |
| 8 | 7.236·10⁻⁶, 7.106·10⁻⁶, 0.0463 | 4.108·10⁻³, 3.940·10⁻³, 0.0468 |
| 10.8049 (x_pos) | 1.2906·10⁻⁸, 1.2819·10⁻⁸, 0.0810 | 1.3488·10⁻⁵, 1.2635·10⁻⁵, 0.0800 |

At x_pos the certified constants 1.0441·10⁻⁸ (even) and 1.0931·10⁻⁵ (odd) lie below these N = 80 upper bounds, as they must. The pair's negative term is already 10⁷ (even) and 10⁴ (odd) times larger than the floor there, in the naive sense, and σ₁ is 0.08: the window is twelve units of x short of seeing the pair.

### 3.4 What sets the forward failure

With σ₁ = m sᵀ(P + m ccᵀ)⁻¹s and P's spectrum λ₁ ≪ λ₂ ≪ …, σ₁ ≈ 4Σ_i ⟨s, e_i⟩²/λ_i is dominated by P's near-null space (the functions whose transforms vanish at the on-line zeros below T* and leak above), weighted by ⟨s, e_i⟩ ≈ δ F̂′_{e_i}(γ₀): **the sinh channel reads the slope of the near-null transforms at the pair's height.** Failure comes when the leakage of P's near-null space at γ₀ has grown to the floor: for DH at x_c, T* = 2πx_c/5 = 38.6 and γ₀/T* = 2.22; for Z₁, T* = 2π√(x_c/20) = 6.26 and γ₀/T* = 2.50; DH's odd sector fails at γ₀/T* ≈ 2.16 and Z₁'s at ≈ 2.4. The planted GL(2) cases of §6 show that this ratio is **not** universal (2.1 for γ₀ = 30, δ = ¼; 1.7 for δ = 0.1; above 2.1 and not yet reached at x/N = 20 for γ₀ = 60), so no closed-form forward law is claimed here; R6 remains open, with σ₁(x) as its measurable.

## 4. The regions for ζ

Scripts: `scripts/research_r1/zeta_region.py` (matrix route: rh2's full zeros-side matrices E + 2vvᵀ and E − 2wwᵀ from `control_normalisation.rh2_matrix("zeta")`, eigendecomposed, then R₄ on a γ₀ grid up to the basis frequency ω_N = 2πN/L) and `scripts/research_r1/box_test.py` (test-function route: f = cos(γ₀u) on [−a, a], Q(f) from Zhu's frequency-side formula Q = pole + (1/π)∫Ψ|F̂|², which needs no matrix and works at any height; it gives the lower bound R₄ ≥ 4⟨c,f⟩²/(Q(f) + 4⟨s,f⟩²)). Data: `data/research_r1/zeta_all.json`, `box_*.json`.

### 4.1 Three regimes [N]/[D]

Write κ_a(γ₀) := sup_f |F̂(γ₀)|²/Q(f) over f supported in [−a, a] (so R₄(0, γ₀) = 4κ and the test for a *simple* on-line zero at γ₀ is 2κ ≤ 1).

**(i) The pinned regime, γ₀ ≲ T* = 2πx.** At the true zeros, 2κ_a(γ_k) = 1 to the working precision, and it leaves 1 as soon as γ_k is moved: at a = 0.8, 2κ(γ₁ + 10⁻⁵) = 317, 2κ(γ₁ + 10⁻¹⁰) = 1.0000000000; at a = 1.495, 2κ(γ₁ + 10⁻²⁰) = 1.7·10⁵⁰. Between zeros 2κ is 10⁵–10¹² (a = 0.8) and 10²⁵–10⁸⁹ (a = 1.495). The window therefore knows the positions and the simple multiplicity of every zero below T* (to about √λ_min in height: 10⁻⁹ at a = 0.8, 10⁻⁴⁸ at a = 1.495), and R₄ ≥ 2 at the zeros themselves and ≫ 1 between them: **every quadruple with γ₀ < T* is excluded under H, at any δ, including a quadruple that would replace one of ζ's simple zeros.** The table below gives the comb at each support.

| a (2a) | sector, N | λ_min(Q_N) | T* | zeros below T* | min 2κ(γ_k) below T* | 2κ(γ₁ + 10⁻⁵) | 2κ(γ₁ + 10⁻²⁰) | 2κ at γ₅₀ = 143.1 | smallest midpoint 2κ | largest midpoint 2κ |
|---|---|---|---|---|---|---|---|---|---|---|
| 0.8 (1.6) | even, 100 | 1.74·10⁻¹⁷ | 31.1 | 4 | 1 − 1.0·10⁻¹⁰ | 351 | 1.0000000000 | 0.629 | 0.419 (above 4T*) | 7.4·10¹¹ (between γ₁, γ₂) |
| 0.8 (1.6) | odd, 100 | 1.67·10⁻¹⁴ | 31.1 | 4 | 1 − 7·10⁻¹⁰ | 7.16 | 1.0000000000 | 0.630 | 0.423 | 2.1·10¹⁰ |
| 0.8 (1.6) | even, 160 | 1.67·10⁻¹⁷ | 31.1 | 4 | 1 − 10⁻¹⁰ | — | 1.0000000000 | — | 0.420 | 7.7·10¹¹ |
| 1.19 (2.38) | even, 100 | 1.14·10⁻⁴⁷ | 67.9 | 16 | 1 − 10⁻¹² or better | 2.2·10³² | 216 (1 + 2·10⁻¹⁰ at +10⁻³⁰) | 1 − 9·10⁻¹⁰ at γ₃₀ = 101.3 (1.5 T*) | 1.0·10⁵ (γ₂₉–γ₃₀) | 2.7·10⁴¹ |
| 1.3 (2.6) | even, 100 | 1.04·10⁻⁶¹ | 84.6 | 22 | 1 − 10⁻¹² or better | 2.0·10⁴⁶ | 2.0·10¹⁶ (1.0002 at +10⁻³⁰) | 1.0000000000 at γ₃₀ (1.2 T*) | 1.3·10¹⁵ | 2.3·10⁵⁵ |
| 1.495 (2.99) | even, 120 | 1.04·10⁻⁹⁵ | 124.9 | 30 (all tabulated; 50 in a first run) | 1 − 10⁻¹² or better | 1.7·10⁸⁰ | **1.7·10⁵⁰** (1.7·10³⁰ at +10⁻³⁰) | 1.0000000000 at γ₃₀; 2κ = 1 to 10 digits at γ₅₀ = 143 (first run) | 2.7·10⁴⁴ | 1.7·10⁸⁹ |

(The last three rows come from `data/research_r1/zeta_high.json`, 30 zeros each, even sector; the pinning sharpness scales as √λ_min: 10⁻⁹, 10⁻²³, 10⁻³⁰, 10⁻⁴⁸ in height at the four supports, read off the shift at which 2κ departs from 1.)

Reading the a = 0.8 rows: 2κ = 1 at γ₁…γ₄ (14.1 to 30.4, below T* = 31.1), 0.9999999994 at γ₅ = 32.9, 0.99997 at γ₈ = 43.3 (1.4 T*), 0.993 at γ₁₂ = 56.4 (1.8 T*), 0.93 at γ₁₃, and 0.43–0.73 from γ₂₇ ≈ 95 (3 T*) on; the pinning fades over [T*, 2T*]. At a = 1.495 every tabulated zero up to 143 (1.15 T*) is pinned at the 10⁻⁵⁰ level: moving γ₁ by 10⁻²⁰ changes 2κ from 1 to 1.7·10⁵⁰. The odd sector pins as tightly as the even one (its λ_min is 10³ larger but so is its sensitivity).

**(ii) The local-count regime, T* ≲ γ₀ ≲ T_excl.** Above T* the zeros are denser than the window's resolution π/a, κ falls below ½, and the test becomes a local zero count. A box test function f = cos γ₀u gives f̂_c(γ₀) ≈ sinh(δa)/δ, ‖f‖² ≈ a and Q(f) ≈ [log(γ₀/2π) − ρ_a(γ₀)]·a, where ρ_a(γ₀) = Σ_{n<x}(2Λ(n)/√n)(1 − log n/2a) cos(γ₀ log n) is the prime side's local fluctuation of the zero density at resolution a (its envelope A_box = Σ(2Λ(n)/√n)(1 − log n/2a) is 1.05 at x = 4.95, 2.5 at 10.8, 3.0 at 13.5, 4.1 at 19.9). Hence

  R₄(δ, γ₀) ≳ 4a·w(δa) / (log(γ₀/2π) − ρ_a(γ₀)),  w(z) = (sinh z/z)² ∈ [1, 1.2] for δ ≤ ½, a ≤ 1.5,

and a quadruple is excluded while the window's own count of zeros within ~π/a of γ₀, read from the primes n < x, is below about four zeros' worth. The matrix optimum improves the box by a shaping gain measured at a = 0.8 (§4.2). The boundary T_excl is where the last excursion of R₄ above 1 ends; it is set by 2πe^{4a·w}·e^{ρ} with ρ an almost-periodic function of γ₀, so it is a band of heights, not a number. **Its envelope is 2π·exp(g·4a·w(δa) + A_box)** with g the shaping gain.

**(iii) Above T_excl: nothing.** R₄ < 1 for every f: the window admits a positive measure on the line that reproduces its form together with the quadruple.

### 4.2 Numbers [N]

**Matrix route at a = 0.8** (the only support where the basis reaches the local-count regime: ω_N = 393 at N = 100, 628 at N = 160). R₄(0, γ₀) oscillates about 1 with the prime phases above ~3T*: at γ₀ = 60, 100, 150, 200, 216 it is 1.46, 1.35, 0.97, 1.34, 1.15 (N = 100), against the box values 1.05, 1.27, 0.90, 1.28, 0.98 at the same heights, a shaping gain g = 1.05–1.17 (1.39 at 60, where the pinned regime still helps). The last excursion above 1: **T_excl(δ = 0) = 435 and T_excl(δ = ½) = 481 at N = 160** (14.0 and 15.5 T*, both well inside ω_N = 628; at N = 100 the scan was cut by the basis at 372–390). The δ = ½ value is higher only because the weight w(δa) = 1.05 lifts the next prime-phase peak over 1; it is the same comb of excursions. Both sectors give the same numbers to three digits at N = 100. The excursions follow ρ_a(γ₀): the box's R₄ ≈ 3.2/(log(γ₀/2π) − ρ) exceeds 1 exactly where the three prime powers 2, 3, 4 align to lower the local density (ρ ≈ +0.5 to +1.0).

**Box route at every support** (`box_test.py`; 8 heights per e-fold band around log(γ₀/2π) = 4a + b, b ∈ {−1, −½, 0, ½, 1, 2, 3, 4}; entries are min / median / max of the lower bound R₄ ≥ 4⟨c,f⟩²/(Q(f) + 4⟨s,f⟩²) over the band's samples, δ = 0 first, δ = ½ second):

| support 2a | T* | 2πe^{4a} | A_box | band 4a−1 | 4a | 4a+1 | 4a+2 | 4a+3 | 4a+4 | largest sampled height with R₄ > 1 (δ = 0 / ½) |
|---|---|---|---|---|---|---|---|---|---|---|
| 1.6 | 31 | 154 | 1.05 | 1.07/1.70/2.53 · 1.13/1.79/2.66 | 0.75/0.94/1.33 · 0.79/0.99/1.40 | 0.71/0.77/0.94 · 0.74/0.81/0.99 | 0.50/0.66/0.76 · 0.53/0.69/0.80 | 0.46/0.49/0.57 · 0.48/0.51/0.60 | 0.39/0.44/0.50 · 0.41/0.47/0.52 | 165 / 180 (matrix, N = 160: 435 / 481) |
| 2.38 | 68 | 734 | 2.52 | 0.97/1.18/1.88 · 1.09/1.32/2.11 | 0.80/0.95/1.67 · 0.90/1.07/1.88 | 0.68/0.83/1.09 · 0.76/0.93/1.23 | 0.63/0.71/0.83 · 0.70/0.80/0.93 | 0.53/0.56/0.63 · 0.59/0.63/0.71 | 0.47/0.57/0.66 · 0.52/0.64/0.74 | 1 830 / 2 973 |
| 2.6 | 85 | 1 139 | 3.04 | 1.07/1.28/1.89 · 1.23/1.47/2.17 | 0.84/1.17/1.34 · 0.97/1.34/1.54 | 0.67/0.71/0.99 · 0.77/0.82/1.14 | 0.62/0.71/0.85 · 0.72/0.81/0.98 | 0.58/0.67/0.72 · 0.67/0.77/0.82 | 0.47/0.56/0.69 · 0.55/0.64/0.79 | 1 850 / 2 842 |
| 2.99 | 125 | 2 485 | 4.06 | 1.04/1.46/2.86 · 1.25/1.76/3.44 | 0.83/0.99/1.95 · 0.99/1.19/2.34 | 0.78/0.89/1.18 · 0.93/1.07/1.42 | 0.63/0.76/0.89 · 0.75/0.91/1.07 | 0.59/0.65/0.70 · 0.70/0.78/0.84 | 0.49/0.62/0.70 · 0.59/0.74/0.84 | 10 071 / 27 375 |

Reading across a row: the median falls like 4a/(log(γ₀/2π)) as it should; the spread within a band is the prime-phase term ρ_a (observed range −2.1 to +3.1 at 2.99, against the envelope A_box = 4.06); the δ = ½ column is the δ = 0 column times w(δa) = (sinh δa/δa)² = 1.05 (a = 0.8) to 1.20 (a = 1.495), never more. At support 2.99 the box alone exhibits inconsistent quadruples (R₄ > 1) up to height 10 071 (δ = 0; where ρ = +1.42) and 27 375 (δ = ½), and the shaping gain pushes the last excursions somewhat higher; the envelope 2π·exp(1.2·4a·w(δa) + A_box) is 4.7·10⁵ (δ = 0) and 2.0·10⁶ (δ = ½).

**Figure-ready summary** (heights in the (δ, γ₀) half-plane; "excluded" means R₄ > 1 under H; both sectors behave alike, the odd sector's values at a = 0.8 agreeing with the even's to three digits):

| support 2a | pinned: all γ₀ below | local-count excursions seen up to (δ = 0 / ½) | envelope beyond which nothing is excluded (δ = 0 / ½) | RH verified to |
|---|---|---|---|---|
| 1.6 | 31 (2κ = 1 − 10⁻¹⁰) | 435 / 481 (matrix) | 8·10² / 1.0·10³ | 3·10¹² |
| 2.38 | 68 | 1.8·10³ / 3.0·10³ (box) | 2.4·10⁴ / 4.8·10⁴ | 3·10¹² |
| 2.6 | 85 | 1.9·10³ / 2.8·10³ (box) | 6.7·10⁴ / 1.7·10⁵ | 3·10¹² |
| 2.99 | 125 (2κ = 1 − 10⁻⁴⁷) | 1.0·10⁴ / 2.7·10⁴ (box) | 4.7·10⁵ / 2.0·10⁶ | 3·10¹² |

Every excluded height at every certified support is below 3·10¹² by at least six orders of magnitude.

### 4.3 Real zeros, even sector [N]

For a real pair ½ ± δ the even-sector consistency test is 2cᵀQ_even⁻¹c with c_k = ∫ b_k cosh δu. At a = 0.8 it is 8.1·10¹⁶ for every δ ∈ {0.05, 0.1, 0.25, 0.4, 0.49}; at a = 1.19, 1.3, 1.495 it is 1.3·10⁴⁷, 1.4·10⁶¹, 1.5·10⁹⁵ (≈ 1.4/λ_min, flat in δ to 1%): the near-minimiser e₁ of Q_even has ∫e₁ cosh δu = O(1) while Q(e₁) = λ_min. **The even sector excludes every real pair under H by consistency**, although it cannot by sign (Yoshida's exception [S]). ζ has no real zeros in (0, 1) [S: elementary], so this is a check of the mechanism, not news; for F_{t*} it is the forward mirror of §3 (the real pair raises the even minimum from L(χ₅)'s floor by 10²–10⁵).

### 4.4 Comparison with what is known

- **RH verified to height 3·10¹²** [A: D. Platt and T. Trudgian, "The Riemann hypothesis is true up to 3·10¹²", Bull. LMS 53 (2021) 792–797; the method isolates zeros by sign changes with Turing's count, so it also gives simplicity — [S], not re-read here]. Every point of R(a) at the four certified supports, in both sectors and for every δ, has γ₀ below the envelopes of §4.2, i.e. below 2·10⁶ at support 2.99 and 10³ at 1.6 (observed excursions: 2.7·10⁴ and 435). All of it lies inside the verified strip.
- **Classical zero-free regions near Re s = 1** [A, via `docs/LIT_TANGENTS.md`]: Mossinghoff–Trudgian–Yang (2024), σ ≥ 1 − 1/(5.558691 log|t|) for |t| ≥ 2; Bellotti–Trudgian–Yang (2026) with 4.896; Ford's Vinogradov–Korobov shape σ ≥ 1 − 1/(57.54 (log t)^{2/3}(log log t)^{1/3}) [S]. These exclude δ close to ½ at every height. R(a) is δ-blind (§2.3) and height-limited, so it is disjoint in kind from these: it says nothing at any height they cover beyond 3·10¹².
- **Under the stronger hypothesis** that all zeros below γ₀ are known and on the line (available only for γ₀ ≤ 3·10¹²), subtracting their sum P_known sharpens the test to Q − T − P_known ⪰ 0 and the box estimate to 4a·w(δa) > ½(log(γ₀/2π) − ρ), i.e. heights up to ~2πe^{8a·w(δa)} ≈ 10⁷ at a = 1.495, δ = ½. This is still below 3·10¹², and above 3·10¹² the hypothesis is unavailable. (Such a bound is a short-interval zero-count statement of the kind Goldston–Gonek-type arguments obtain from the explicit formula [S/U: not re-read here]; it is not a positivity statement.)

**Plainly: no certified support excludes a point (β, γ₀) that is not already excluded**, and none of them can, because the window's reach in height is set by 2πe^{4a}·(prime phases), not by how small its certified minimum is.

## 5. Go/no-go for supports 3.2, 3.5, 4.0

**What a larger support would exclude (under H).** The reach in height is the envelope of §4.1, T_env(a, δ) = 2π·exp(g·4a·cosh²(δa) + A_box(x)), with g ≈ 1.2 the shaping gain measured at a = 0.8 and A_box(x) = Σ_{n<x}(2Λ(n)/√n)(1 − log n/2a) the prime side's maximal local-density excursion. The envelope is pessimistic twice over: it assumes the prime phases align perfectly at the height in question (ρ = A_box, an event of probability ~e^{−A_box²/2σ²} per independent height), and it is an upper bound on where R₄ > 1 can occur, not a height that is excluded throughout. The typical reach is 2π·e^{g·4a cosh²(δa)}·e^{ρ} with |ρ| ≲ 1.

| support 2a | x | T* | 2πe^{4a} | A_box | typical reach (ρ = +1), δ = 0 / ½ | envelope (ρ = A_box), δ = 0 / ½ | exponent g·4a·w(δa) + A_box, δ = 0 / ½ (3·10¹² needs 26.9) |
|---|---|---|---|---|---|---|---|
| 1.6 | 4.95 | 31 | 154 | 1.05 | 8·10² / 1.0·10³ | 8·10² / 1.0·10³ | 4.9 / 5.1 |
| 2.38 | 10.8 | 68 | 734 | 2.52 | 5·10³ / 1.1·10⁴ | 2.4·10⁴ / 4.8·10⁴ | 8.2 / 9.0 |
| 2.6 | 13.5 | 85 | 1.1·10³ | 3.04 | 9·10³ / 2.2·10⁴ | 6.7·10⁴ / 1.7·10⁵ | 9.3 / 10.2 |
| 2.99 | 19.9 | 125 | 2.5·10³ | 4.06 | 2·10⁴ / 1.0·10⁵ | 4.7·10⁵ / 2.0·10⁶ | 11.2 / 12.7 |
| 3.2 | 24.5 | 154 | 3.8·10³ | 4.67 | 4·10⁴ / 2.2·10⁵ | 1.5·10⁶ / 8.6·10⁶ | 12.4 / 14.1 |
| 3.5 | 33.1 | 208 | 6.9·10³ | 5.65 | 8·10⁴ / 8.2·10⁵ | 7.9·10⁶ / 8.6·10⁷ | 14.0 / 16.4 |
| 4.0 | 54.6 | 343 | 1.9·10⁴ | 7.54 | 2.5·10⁵ / 1.0·10⁷ | 1.7·10⁸ / 6.8·10⁹ | 17.1 / 20.8 |

(g = 1.2 throughout; w(z) = (sinh z/z)², the box weight, 1.05–1.38 for δ = ½; "typical" uses ρ = +1, about one RMS of the prime-phase term; the envelope uses ρ = A_box, perfect alignment of every prime power below x.) The box scans at these supports (`box_16.json`, `box_175.json`, `box_20.json`, 64 heights each) found R₄ > 1 up to 3.9·10⁴ / 4.2·10⁴ (support 3.2, δ = 0 / ½), 1.8·10⁴ / 1.3·10⁵ (3.5) and 3.9·10⁴ / 9.4·10⁵ (4.0), consistent with the typical column; the band medians cross 1 between 4a and 4a + 1 in log(γ₀/2π), as at the certified supports.

**Reading.** No support up to 4.0 reaches 3·10¹², even on the envelope: at 4.0 with δ = ½, every prime power below 55 aligned and the shaping gain included, the exponent is 20.8 against the 26.9 needed, a factor 440 short in height. To close that gap the shaping gain would have to be g ≥ 1.75 (δ = ½) or g ≥ 2.0 (δ = 0) at a = 2, against the 1.05–1.2 measured at a = 0.8; or the test function would have to beat the box's local-count reading by more than the explicit formula's own fluctuation allows. Suppose nevertheless that some support did reach a height just above 3·10¹² at a few prime-aligned heights. What it would state there is that a quadruple is inconsistent with "all other zeros on the line", and three things qualify that. (i) It follows from evaluating Q(f) for the single test function f = cos γ₀u at that height: a sum over the prime powers n < x plus an archimedean integral, computable in a second and needing no positivity certificate; the cost of the certificate buys nothing here. (ii) It is a short-interval zero-count bound ("the window-smoothed density at γ₀ is far enough below its mean that four extra zeros do not fit"), not a statement about the critical line: the window excludes the δ = 0 double zero at the same heights just as well (§2.3). (iii) It assumes RH for every other zero; dropping that, the other zeros' sinh terms can be negative and nothing follows. None of it uses the certified minimum c.

**Cost** [D, rough; the real model is R4's]. From the quarter law μ̄ ≈ e^{a+δ_K} and Zhu's threshold T♯ ≈ 2πe^{μ̄}, with Gram cost ∝ nodes × modes² ∝ a²T♯³ and precision ∝ x, scaled from the recorded 2.99 run (μ̄ = 5.28, T♯ = 2870, 6.3 h per sector on 28 threads): support 3.2 needs T♯ ≈ 5·10³ and about 8× the time (≈ 2 days per sector); 3.5 needs T♯ ≈ 1.3·10⁴ and ≈ 200× (≈ 2 months); 4.0 needs T♯ ≈ 9·10⁴ and ≈ 2·10⁵× (a century). These are order-of-magnitude figures for a method that is doubly exponential in a; a prolate envelope (R4b) could shave the constant, not the shape.

**Recommendation (one line): no-go.** Pushing the certified support from 2.99 to 3.2–4.0 would raise the exclusion reach from ~10⁵–10⁶ to ~10⁷–10¹⁰ in height, under a hypothesis that already concedes RH elsewhere and through a test that needs no certificate, all of it below the 3·10¹² that is verified unconditionally, at a cost from days to a century; the certificates' scientific value is as anchors for the decay law (R2, R3), and that is what further supports should be justified by, if at all.

## 6. Pre-registered experiment: planted quadruples in validated GL(2) forms

Written before any 37b1 run exists; the 11a1 numbers below were computed here (`scripts/research_r1/plant_gl2.py`, `data/research_r1/plant_11a1_*.json`) and are part of the registration.

**Objects.** The validated GL(2) forms of `docs/GL2_DECAY.md` for 11a1 (N = 11) and 37b1 (N = 37), both rank 0, ε = +1, Γ_ℂ(s + ½), zeros to T = 120 all on the line (`gl2_form_mp.zeros_side_E`). Into each, plant a quadruple by adding Proposition 1's term: Q̃ = Q_E + 4(ccᵀ − ssᵀ), for (δ, γ₀) ∈ {(0.25, 30), (0.10, 30), (0.25, 60)}, both sectors. Q̃ is the zeros-side form of the entire function Λ_E(s)·Π(s − ρ_i) over the orbit (same functional equation, four extra zeros), so it is a control in the sense of `docs/CONTROL_CERTIFICATES.md`, with P = Q_E known exactly.

**Grid and measurement.** x/N ∈ {1, 1.5, 2, 2.5, 3, 4, 5, 6, 8, 10, 13, 16, 20} and bisection of the sign change of λ_min(Q̃) to 10⁻³ in x/N, at N_basis = 5k* and 9k* (k* = T*L/2π, T* = 2π√(x/N)), precision-stable (`stable_lams`), FLINT's full eigensolver (negative eigenvalues must be seen). The crossing x_c is the 9k* value. For (0.25, 60) the grid extends to x/N = 60 if no sign change appears by 20.

**Computed here (11a1, N_basis = 64, 60 digits):**

| case | sector | x_c (x/N) | σ₁ at x/N = 3, 4, 5, 6 |
|---|---|---|---|
| (0.25, 30) | even | 55.661 (5.0601) | 0.111, 0.220, 0.884, 2.86 |
| (0.25, 30) | odd | 59.128 (5.3753) | 0.100, 0.225, 0.488, 2.48 |
| (0.10, 30) | even | 85.418 (7.7652) | 0.0171, 0.0334, 0.132, 0.428 |
| (0.10, 30) | odd | 90.920 (8.2654) | 0.0154, 0.0341, 0.0732, 0.369 |
| (0.25, 60) | even | > 220 (> 20) | 0.0886, 0.115, 0.150, 0.175; 0.407 at x/N = 20 |
| (0.25, 60) | odd | > 143 (> 13) | 0.0846, 0.117, 0.148, 0.170; 0.322 at x/N = 13 |

**Predictions.**
- **P1 (11a1, consistency of the measurement).** The 9k* crossings for the four γ₀ = 30 cases lie within 3% of the x_c above (finite-basis crossovers move down with N; 64 lies between 5k* ≈ 45 and 9k* ≈ 81 at these x).
- **P2 (sector order).** In every case the odd crossing is above the even one (ε = +1 ⇒ λ_odd > λ_even ⇒ the odd near-null space leaks less at γ₀). Kill: any case with the odd crossing below the even one.
- **P3 (δ² law).** At every grid point before the first crossing, σ₁(0.10, 30)/σ₁(0.25, 30) ∈ [0.14, 0.16] in both sectors (leading order (0.1/0.25)² = 0.16, with O(δ²) corrections; computed values here 0.154 → 0.149 over x/N = 3 → 5). Kill: a ratio outside [0.12, 0.18] at any point with σ₁(0.25, 30) ≤ 0.5.
- **P4 (conductor transfer; the one with risk).** For 37b1 with the same planted (0.25, 30) and (0.10, 30), the crossings in x/N lie within ±10% of 11a1's: even 5.06 and 7.77, odd 5.38 and 8.27. This asserts that the forward failure depends on (x/N, γ₀, δ) and the gamma factor only, not on the conductor or the particular low zeros, as the decay law's x/N collapse (`docs/DECAY_LAW.md`, P-DL4: 11a1 and 37b1 agree to ln-ratio 0.1–1.7 in λ) would suggest. Kill: any of the four 37b1 crossings more than 25% from 11a1's in x/N. A kill would mean the local zero configuration near the pair, not the window's density, controls the failure, and R6 must be formulated per object.
- **P5 (the far pair).** For (0.25, 60) the even crossing of 11a1 lies at x/N ∈ [30, 60]; the universal-ratio rule γ₀/T*(x_c) ≈ 2.1–2.5 (which fits DH, Z₁ and the γ₀ = 30 cases) predicts x/N ≈ 20 and is already contradicted by σ₁(20) = 0.41, so it is **not** registered; what is registered is that σ₁ continues to rise monotonically and crosses 1 below x/N = 60. Kill: no crossing by x/N = 60, which would mean the leakage at 2T* saturates and a pair that far out is never seen in this basis size.

**Not predicted.** The extrapolation protocol of §3.2 from the cheap regime x/N ≤ 3 overshoots the γ₀ = 30 crossings by 17–22% here (LS lines through five points with σ₁ ≤ 0.11 give x/N = 6.18, 6.56, 9.15, 9.69 against 5.06, 5.38, 7.77, 8.27), and for (0.25, 60) it predicts x/N ≈ 6.4 against > 20. **The protocol is only trustworthy when the data reach σ₁ ≳ 0.3**, as they did for the three controls of §3; from further away ln σ₁ is too far from linear. This is recorded as a limitation, not registered as a prediction.

## 7. The three things I trust least

1. **The hypothesis H is doing all the work, and it is strong.** "Every zero other than this quadruple is on the line" concedes RH everywhere else; under weaker hypotheses (other off-line zeros allowed) the quadruple's negative partner terms from those zeros make Q − T indefinite in general and nothing is excluded. I have not found a hypothesis-free statement, and I believe there is none at these supports: positivity of a window form is compatible with off-line zeros (DH, Z₁ at their certified-positive supports are the proof). The Krein-extension remark in §2.2 is a sketch at the level of "standard"; the archimedean kernel's log singularity and the measure's infinite mass are not handled.
2. **The envelope T_env and the shaping gain g.** The box test function is a lower bound on R₄; the matrix optimum exceeds it by g = 1.05–1.17 at a = 0.8 in the range 100–216 (1.39 at 60, still in the pinned transition) and g was *not* measured at the other supports, where the matrix route is basis-limited (ω_N ≈ 250 against reaches of 10³–10⁵). I used g = 1.2 at every support; if g grows with a (plausible: more primes, more room to shape), the envelopes in §4.2 and §5 are too small in the exponent by (g − 1.2)·4a·w(δa), i.e. by a factor e^{0.7} ≈ 2 per 0.1 of g at a = 1.495. The conclusion "all below 3·10¹²" has a margin of e^{14} at 2.99 for δ = ½ and survives g up to 3.2 there; at support 4.0 the margin is e^{6} and the conclusion fails for g ≥ 1.75. A matrix measurement of g at a ≥ 1.2 (N ≈ 2a·γ₀L/2π ≈ 1000 modes at γ₀ = 2πx², i.e. hours of mp.eigsy, or a Lanczos solve in FLINT) is the one cheap computation that would firm this up.
3. **The calibration is of the algebra and the census, not of a forward law.** Rung 2 is exact given P, and P is defined as Q minus the pair terms; its positivity and Euler-class floor are the genuine tests, and they pass. But I have no derived expression for σ₁(x; δ, γ₀) — the forward failure law (R6) — and the two empirical rules I tried (γ₀/T*(x_c) ≈ 2.1–2.5; linear ln σ₁ extrapolation) each fail on at least one case (§3.4, §6). The pre-registration's P4 is the first real test of whether such a law exists within a gamma-factor class. The planted controls also use my own T; an independent construction of a function with a chosen complex off-line quadruple (not available in rh2) would close the loop better than DH/Z₁, whose census zeros above 300 (DH) and 60 (Z₁) are unknown and are absorbed into P — their absence from the census is invisible at the supports tested only because the window's reach was below those heights.

Smaller caveats: all matrix numbers are finite-basis (N = 48–160) and not certified; the a = 0.8 T_excl values at N = 100 for δ ≥ 0.25 sit at 0.99 ω_N and are basis-limited (the N = 160 run is the one to quote); the box route's far tails are averaged (sin² → ½, comb → 0), an O(1/(aW²)) relative error with W = 150; Platt–Trudgian and the zero-free-region constants are cited from abstracts and from `docs/LIT_TANGENTS.md`, not re-read.

## 8. Files and reproduction

New files only; no existing script or document was modified.

| path | what |
|---|---|
| `scripts/research_r1/r1lib.py` | pair vectors c, s in rh2's bases (closed forms), the quadruple matrix, direct-quadrature check, Sherman–Morrison tests (reverse R₄, forward σ₁), naive minimum |
| `scripts/research_r1/check_T.py` | Proposition 1 against quadrature; the conductor-5 numbers |
| `scripts/research_r1/controls.py` | Q = Σ T_i + P for F_{t*}, DH, Z₁ on x-grids; λ_min(Q), λ_min(P), σ₁, naive |
| `scripts/research_r1/fit_sigma.py` | rung 3: ln σ₁ extrapolation from x ≤ x_pos |
| `scripts/research_r1/zeta_region.py` | ζ: comb 2κ(γ_k), R₄ profiles and T_excl up to ω_N, even real-pair test |
| `scripts/research_r1/box_test.py` | ζ at any height: Q(cos γ₀u·1_[−a,a]) from Zhu's Ψ, R₄ lower bound, band sampling of the prime phases |
| `scripts/research_r1/box_summary.py` | the tables of §4.2 and §5 from the JSONs |
| `scripts/research_r1/plant_gl2.py`, `plant_fit.py` | §6: planted quadruples in 11a1 (and 37b1 when run), σ₁(x), bisected crossings, cheap-regime extrapolation |
| `data/research_r1/ctrl_*.json` | §3 (controls); `*_low` are the DH points of §3.3 |
| `data/research_r1/zeta_all.json` | §4: a = 0.8, N = 100, both sectors, seven δ values (the N = 100 T_excl values are basis-limited) |
| `data/research_r1/zeta_08_n160_d0.json`, `zeta_08_n160_d05.json` | §4: a = 0.8, N = 160, δ = 0 and ½ (T_excl = 435 and 481) |
| `data/research_r1/zeta_high.json` | §4: a = 1.19, 1.3, 1.495, even sector, 30 zeros, δ ∈ {0, ½} (combs; R₄ > 1 throughout the basis range) |
| `data/research_r1/box_*.json` | §4.2, §5: box bands at a = 0.8 … 2.0; `box_08_cmp.json` the matrix comparison heights |
| `data/research_r1/plant_11a1_*.json` | §6 |
| `data/research_r1/logs/*.out` | the run logs (including the first a = 1.495 N = 120 run with 50 zeros, `zeta_1495.out`) |

Commands (interpreter `<repo>/.venv/bin/python`, run from the repository root; logs in `logs/` via `scripts/progress.py`):

```bash
PY=.venv/bin/python; R=scripts/research_r1; D=data/research_r1
$PY $R/check_T.py
$PY $R/controls.py --function ftstar --sector odd  --x 2,3,4,4.5,5,5.2,5.31217,5.3625,6,7 --n 64 --dps 40 --json $D/ctrl_ftstar_odd.json
$PY $R/controls.py --function ftstar --sector even --x 3,5,7,8.84631,10 --n 48 --dps 40 --json $D/ctrl_ftstar_even.json
$PY $R/controls.py --function dh --sector even --x 6,8,10.8049,10,13,16,20,25,28,30,30.745,31.25,32,33 --n 80 --dps 70 --json $D/ctrl_dh_even.json   # and --sector odd
$PY $R/controls.py --function z1 --sector even --x 5,7,10,13,15,15.8,17.5,19,19.844,20,22 --n 64 --dps 50 --json $D/ctrl_z1_even.json          # and --sector odd
$PY $R/zeta_region.py --jobs 0.8:even:100:50,0.8:odd:100:50 --json $D/zeta_all.json                       # seven δ values, 50 zeros, grid to ω_N
$PY $R/zeta_region.py --jobs 0.8:even:160:50 --deltas 0 --nzeros 10 --json $D/zeta_08_n160_d0.json         # and --deltas 0.5 → zeta_08_n160_d05.json
$PY $R/zeta_region.py --jobs 1.19:even:100:90,1.3:even:100:110,1.495:even:120:150 --deltas 0,0.5 --nzeros 30 --json $D/zeta_high.json
$PY $R/box_test.py --a 1.495 --delta 0,0.5 --json $D/box_1495.json      # and a = 0.8, 1.19, 1.3, 1.6, 1.75, 2.0; --gammas 60,100,... for the a = 0.8 comparison
$PY $R/plant_gl2.py --cases 0.25:30,0.10:30,0.25:60 --xn 1,1.5,2,2.5,3,4,5,6,8,10,13,16,20 --sector even --n 64 --dps 60 --json $D/plant_11a1_even.json   # and --sector odd; --label 37b1 for §6
```

Each control run takes 1–4 minutes; the ζ matrix scan 5–20 minutes per job; each box band set about 6 minutes; the planted 11a1 run 2 minutes per sector.
