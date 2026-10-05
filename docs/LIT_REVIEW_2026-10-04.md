# Literature review for publishability (2026-10-04)

Three independent read-only reviewers covered three clusters: the certificates and R1, the R2 bound, and the decay law. The orchestrator spot-checked their key claims, as noted below. Labels: [R] read full text or the relevant sections, [A] abstract or metadata only, [S] secondhand, [V] verified by the orchestrator against the primary record.

## Verdict

The work is publishable, but not as originally framed. Fast-moving, unrefereed, AI-assisted work posted between 14 September and 4 October 2026 overlaps three of our four clusters.

- **The decay law across L-functions is the strongest and most novel result.** Its conductor collapse also answers a published open problem: Zhu's open problem (iv).
- **The certificates keep the largest supports**, but no longer "the only results beyond 2.125".
- **R2 is preceded for ζ by Mori.** What remains ours there is the extension to Dirichlet characters and to each sector separately.

Post to arXiv soon, because priority risk is high.

## 1. Certificates, controls and R1 (C1–C4)

| Claim | Verdict | Closest prior work | Still new |
|---|---|---|---|
| C1: certified positivity at supports 2.38, 2.6, 2.99 | partly known | zeta-lab PR #258 (merged 28 Sep [V]): 2.38, even sector only, Q ≥ 5.7179e-48. Liu (Math. Comp. submission, repo 14 Sep [V]): 2.0 and 2.125. Zhu v2: 1.6, with 2.38 withdrawn | 2.6 and 2.99 in both sectors; 2.38 for all complex f with the sharper constant 6.813e-48 |
| C2: methods | partly known | Weighted-Schur compressed prime translations: Liu Thm B. The quarter law is essentially zeta-lab's theory Thm 3. Exact-residual Cholesky is standard | the Kaiser frequency split inside Zhu's reduction |
| C3: ζ_K positive where Z₁ is negative | partly known | Arda PR #631 (open, 25 Sep [V]): the same pair, interval-certified but conditional on four paper-level steps. zeta-lab PR #253 (open, 25 Sep [V]): Epstein and DH negative, DH crossing at x ≥ 30.617 | an unconditional (given stated lemmas) Arb certificate in both sectors; the two-sided brackets |
| C4: R1, exclusion as a consistency test | novel formulation | Booker 2006 (explicit-formula zero location); Kim et al. 2026 (injected off-line zeros) | the rank-2 quadruple term, R₄, "the certified constant plays no role", and the reach forecast, which stays heuristic |

**Correction made.** `docs/PRIMER.md` said "prior art stops at 2.125 (Liu)". That is now false.

**Referees will ask for:**
- an independent reimplementation at 2.99, in another interval library or by another person;
- a second audit round, covering the post-audit diff and the Kaiser and Bessel enclosures;
- written proofs of every inequality in an appendix;
- certificate files and code archived with a DOI;
- the quarter law cited as zeta-lab's, or proved with explicit constants;
- the function class stated precisely (L² on a closed interval, against C_c^∞);
- the AI-assistance disclosure.

**Caution.** A record alone is thin. The barrier grows doubly exponentially, so larger fixed windows say nothing about RH.

## 2. R2: the unconditional upper bound (D1–D6)

**Priority for ζ belongs to Ryota Mori** (Zenodo, unrefereed and AI-assisted; dates checked against Zenodo's API [V]).
- 10.5281/zenodo.23059298, 30 Sep (v2), then 23066578, 1 Oct (v3): *Unconditional doubly exponential upper bounds for the bottom of windowed Weil quadratic forms*. Thm A is the Hermite-type rate. Thm B uses Kaiser–Bessel windows to reach e^{−4πμ}, with μ = x. The existential form is checked in Lean.
- Sequels:
  - 23092542 (2 Oct): the CCM prolate vector, μ⁸ e^{−4πμ};
  - 23092454 (2 Oct): the prolate defect;
  - 23119609 (3 Oct): μ^{9/2}(log μ)⁴;
  - 23134085 (4 Oct): a conditional lower bound under RH plus a local-pairs hypothesis, e^{−Cμ log log μ}.
- Our first R2 commit was 3 Oct.

| Claim | Verdict | Note |
|---|---|---|
| D1: E-map identity | partly known | the classical theta and Poisson relation; Mori A Prop. 3.5 for ζ. The twisted, two-sided, sector form for real χ was not found |
| D2: Theorem 2 | partly known | Mori uses the same strip-only argument for ζ, with Σ Re(1/ρ) in place of explicit zero counts |
| D3: Theorem C | partly known, low added value | Zhu's certified upper bounds for ζ are far sharper; the character values are new |
| D4: Proposition H | known for ζ; new for χ and for each sector | |
| D5: Conjecture R | **a theorem for ζ (Mori A Thm B, Mori II Thm E)** | Mori A Lemma 5.3 is our M2, and his Lemmas 5.5–5.6 prove M1 non-effectively. Adapting them should prove Conjecture R for χ and both sectors. The real gap is the second-order cancellation in the minus sector |
| D6: lower bound RH-hard, conditional sketch | known or superseded | Mori IV's conditional rate is stronger than our §4.3 sketch |

**Also outdated:** R2_THEOREMS §9's "the sharp prolate trial gives nothing" (Mori II handles it). The zero-count chain is unnecessary: the Σ Re(1/ρ) trick replaces Bellotti–Wong and Platt.

**Use.** One section of the paper, reframed as the Dirichlet and per-sector extension of Mori's bounds, reached independently. Either prove Conjecture R by adapting Mori's argument or drop it.

## 3. The decay law (E1–E6)

| Claim | Verdict | Closest prior work |
|---|---|---|
| E1: λ ≈ e^{−2d·T*}, collapse in x/q | **novel beyond ζ** | Zhu Conj. 12.1 covers ζ only, and his open problem (iv) asks for Dirichlet L-functions. Connes 2026 §6.4 treats the prolate quantity only. Kim et al. use ζ only |
| E2: index rule | novel as an empirical rule | its ingredients are classical (Fuchs–Slepian; CCM 2024 Prop. 5.1) |
| E3: sector order follows m₀ | partly known mechanism, novel observation | the central term in rank bounds (Mestre 1986, Bober 2013); Bailleul 2021 |
| E4: products lie far above their worst factor | novel, but nearly a corollary of E1 | |
| E5: Borromean readout | known method; a demonstration | Landau, Ford–Zaharescu; Amano et al. 2013 leave the real case open. LMFDB 2.145.4t3.b.a is our ρ |
| E6: primitive against product (step 2) | novel for Weil's form | Murty–Perelli 1999 (zeros detect primitivity, under a pair-correlation conjecture) |

**Headline: the data discriminate against Zhu's law beyond ζ.** Transposed to L(s, χ), Zhu's −ln λ ≈ 2π² N(T*)/ln N(T*) grows with q at fixed x/q.
- **Degree 1.** At x/q = 10 the orchestrator recomputed both from `kb_cert.json`'s finite-basis minima [V]. The measured −ln λ is 113.5–114.7 (even) and 103.6–104.1 (odd) for q = 3, 4, 7, 20, a spread of about 1 nat. Zhu's transposed law predicts 149, 161, 184 and 226, a spread of 77.
- **Degree 2.** The same holds for the matched curve pairs: an observed gap of at most 0.69 against about 7 nats predicted.
- **Caveat.** These are finite-basis upper bounds. The discrimination rests on their observed convergence, which is [E], not certified.
- **On ζ alone the two laws agree to about 2% up to x ≈ 120.** Our ζ fit reproduces Zhu's certified upper bounds out to x = 54.6.

**Referees will ask for:**
- a certified two-sided bracket at a few points for one character and one curve;
- convergence across basis types;
- a head-to-head comparison of models on out-of-sample error;
- one object pushed further, to x/q ≥ 30 or v ≥ 15;
- data with a DOI;
- the full table of registered tests, including the kills;
- "law" stated as "empirical law".

## Recommended paper plan

1. **Main paper (Experimental Mathematics).** "An empirical decay law for Weil's form across L-functions":
   - **Headline:** conductor and degree collapse, and the discrimination against Zhu's law;
   - **Supporting sections:** the certificates (with the corrected related work), the controls, the R2 extension to characters and sectors (citing Mori), and the R1 formulation;
   - **Before submission:** a certified bracket for one character and one curve.
2. **Optional second paper (Math. Comp.).** The certificates at 2.6 and 2.99, after an independent reimplementation and audit round 2.
3. **Now.** An arXiv posting with an honest related-work section citing Mori, zeta-lab, Arda, Liu, Zhu and Kim et al. by date, plus a Zenodo snapshot of the pre-registrations.

## Sources (with status)

- Mori: Zenodo 23059298, 23066578, 23092542, 23092454, 23119609, 23134085 [V dates and titles; R2 reviewer read A, II and IV].
- zeta-lab: github.com/teal-sea/zeta-lab PR #258 [V], PR #253 [V]; theory and numerics RESULTS [R by reviewer].
- Arda: github.com/DrMurphyIsIn/Arda PR #631 [V]; docs [R by reviewer].
- Liu: github.com/luciferyu666/certified-weil-positivity [V, created 14 Sep]; manuscript [R partial by reviewer].
- Zhu, arXiv:2608.24827v2 [R by reviewer; no v3].
- Connes, arXiv:2602.04022 §6.4 [R].
- Connes–Consani, arXiv:2106.01715 and 2006.13771 [R].
- CCM, arXiv:2310.18423 and 2511.22755 [R partial].
- Kim et al., arXiv:2607.24830v3 [A].
- Bailleul, ANT 15 (2021) [R intro].
- Amano et al., Tokyo J. Math. 36 (2013) [R intro].
- Booker, Exp. Math. 15 (2006) [A].
- Mestre 1986; Bober 2013; Murty–Perelli 1999 [A or S].

The reviewers' scratch copies are kept locally and not committed.
