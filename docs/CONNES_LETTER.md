# Connes' "Letter to Riemann": Reproduction and the Davenport–Heilbronn Control

Source: A. Connes, *The Riemann Hypothesis: Past, Present and a Letter Through Time*, [arXiv:2602.04022](https://arxiv.org/abs/2602.04022), §5–7.

**The construction.**
- QW_λ is Weil's quadratic form on even test functions supported in [λ⁻¹, λ]. In the additive variable u = log x this is [−L/2, L/2], with L = 2 log λ.
- Only the n ≤ λ² enter the prime side.
- Let η be the minimiser. Connes–van Suijlekom (Theorem 6.1) show that, when the minimum is simple, every zero of η̂ is real.
- Numerically, those zeros approximate the zeta zeros.

## Corrections from the literature pass (`docs/LITERATURE_PASS.md`)

- **Finite versus full-space positivity.** Positive minima in a finite basis are upper bounds on the full-space minimum; negative ones prove full-space negativity. "Positive for x ≤ …" below always means in the tested finite bases.
- **The zeros-side form** (the scripts' "nopole") adds the pole correction to the geometric form. It is a sum of squares only under RH.
- **Precedents.**
  - Reconstructing several low eigenvectors from E-mapped prolates is already documented: Connes–Consani, arXiv:2106.01715, §3; Connes–Consani–Moscovici, arXiv:2310.18423, §3.6.
  - A "Landau–Widom decay law" for the minimum is conjectured by Zhu, arXiv:2608.24827, Conjecture 12.1.
  - Zhu also states a certified, prime-inclusive positive lower bound at x = e^{1.6} ≈ 4.953 (Corollary 6.3), in an unrefereed preprint. `docs/AUDIT_ZHU.md` re-executes it in interval arithmetic: verified with corrections, and the even constant improves to 1.0277e-17. Our upper bound then brackets the minimum: 1.0277e-17 ≤ λ*(e^{1.6}) ≤ 1.6356e-17.

  Our DH extension and the numerical constants are separate from these.
- **After DH's crossover:** the minimiser is still close to k_λ at x = 32 (sin² ≈ 4.2e-7) and becomes orthogonal only by x ≈ 34.
- **Entry 48** is a suspected typo, not a confirmed one.

## Step 1: double precision (`src/rh/connes_letter.rs`, `rh-connes-letter`)

The basis is cos(2πku/L), k ≤ N. The local terms W_p and W_ℝ (equations 9–10 of the paper) reduce to O(N) one-dimensional integrals and prime sums.

- **Explicit-formula check:** QW(φ) = Σ_ρ |φ̂(γ)|² − 2φ̂(i/2)² to 1e-6 on a C¹ test function.
- **Which form Connes uses:** the full form has one negative eigenvalue, coming from the pole term. His QW_λ is the zeros-side form, on all even φ: the geometric form (archimedean minus primes) **plus** the pole correction 2φ̂(i/2)². It equals Σ_ρ φ̂(γ_ρ)·conj(φ̂(γ̄_ρ)), which becomes the sum of squares Σ|φ̂(γ)|² only if RH holds. The scripts call it `--variant nopole`; the name refers to the zeros side and does not mean the pole is dropped (see Step 2).
- **Zeros reproduced:** at λ² = 7 (n = 2, 3, 4, 5, 7) the first two zeta zeros come out to 4e-15.

## Step 2: high precision (`scripts/connes_letter_mp.py`, mpmath)

The archimedean integrals are computed as exact accelerated series, using K(w) = Σ_m e^{−(2m+½)w}.

At λ² = 13, N = 100 and 130 digits, the computed differences match the letter's table to the six digits printed:

| k | Paper | This work |
|---|---|---|
| 1 | 2.60179e-55 | 2.60179e-55 |
| 2 | 4.80071e-52 | 4.80071e-52 |
| 10 | 4.24869e-37 | 4.24869e-37 |
| 20 | 3.76751e-24 | 3.76751e-24 |
| 30 | 4.47113e-15 | 4.47113e-15 |
| 40 | 2.10931e-7 | 2.10931e-7 |
| 48 | 0.0209081 | 0.00209081 (same digits; suspected typo, not confirmed, since the paper's entries are upper bounds) |
| 50 | 0.00212727 | 0.00212727 |

The smallest eigenvalue is ε₀ = 3.7e-59. Restricting instead to φ̂(i/2) = 0 gives only 7e-50 for the first zero, which identifies the pole-corrected zeros-side form as Connes' QW_λ.

Full output: `data/connes/letter_lambda13_n100.{txt,json}`.

## Step 3: the same construction on Davenport–Heilbronn

DH has the same kind of functional equation as ζ, no Euler product, and zeros off the line. Its form uses:
- the archimedean term for the odd conductor-5 Γ-factor, which is ζ's term plus log 5·F(0) + ½∫₀^∞ sech(u/2)(F(u) + F(−u)) du;
- the prime side c_n/√n from −f′/f.

That prime side includes points off the axes of the prime grid (6, 12, 14, …). There is no pole term. The form matches DH's zeros, with the off-line pairs entering at complex arguments, to 1e-6.

Smallest eigenvalue ε₀ (high precision, N = 64, 60 digits, unless marked):

| λ² | DH ε₀ | ζ ε₀ |
|---|---|---|
| 13 | +8.6e-11 | +3.7e-59 (N = 100) |
| 20 | +4.5e-18 | below 1e-60 |
| 25 | +1.6e-23 | below 1e-60 |
| 28 | +8.4e-27 | below 1e-60 |
| **30** | **+1.3e-28** | \|ε\| ≲ 3e-60 |
| **31** | **−3.8e-30** | |
| 33 | −1.5e-23 | |
| 34 | −1.3e-18 | \|ε\| ≲ 2e-60 |
| 36 | −2.8e-11 | |
| 40 | −0.022 (f64) | |
| 50 | −0.45 (f64, n₋ = 1) | |
| 100 | −1.43, −0.30 (f64, n₋ = 2) | |

Behaviour of the minimiser's zeros (`rh-connes-letter --function dh`):
- **λ² ≤ 30, in the tested finite bases: the construction works for DH.** A positive finite minimum is only an upper bound on the true one, because finite trigonometric minima decrease to the full-space infimum (Connes–Consani, *Spectral Triples and Zeta-Cycles*, Cor. 2.4). A negative finite value does prove a negative direction in the full space. So the first full-space failure lies at λ² ≤ 30.95, with no certified lower bound. The minimiser reconstructs up to 26 of DH's own on-line zeros, to between 1e-15 and 1e-12. Small-support positivity does not need the Euler product.
- **Near the off-line zero 0.8085 + 85.6993i, before the crossover:** Theorem 6.1 keeps η̂'s zeros real. It replaces the off-line pair by one real zero at 85.73 and pushes the neighbouring on-line zeros out by about 0.2.
- **Crossover at λ² ∈ (30, 31), L ≈ 3.42:** ε₀ changes sign. The minimiser then concentrates on the off-line zero (|η̂(ρ)| goes from 1e-9 to 0.46) and the reconstruction collapses.
- **λ² = 100:** two negative directions, one for each off-line pair within reach (85.70 and 114.16).

## What this shows about the strategy

1. **Real zeros of the approximants prove nothing by themselves.** Theorem 6.1 holds for DH at every λ, but DH's zeros are not all real. The weight of the argument rests on the convergence step (§6.6).
2. **For DH, convergence fails exactly when positivity fails.** Below the crossover, DH behaves like ζ. Above it, the minimiser abandons the on-line zeros. A convergence proof for ζ therefore has to exclude a crossover at every λ, and that is Weil positivity: as strong as RH.
3. **ζ and DH differ only once the crossover arrives.** Each function's minimum rides its own prolate floor (Steps 4–5): ζ at c = 2πx, DH at c = 2πx/5 with an odd prolate. ζ's ε₀ is far smaller than DH's at the same λ only because ζ's floor is lower.

## Step 4: the prolate ratio ε(λ) / (1 − χ₂(λ)) (`scripts/prolate_ratio_mp.py`)

The prolate quantity χ₂ is the eigenvalue of the time- and band-limited Fourier transform on [−λ, λ] belonging to h_{4,λ}. It depends only on the archimedean place. In Slepian's notation χ₂² = λ₄(c) with c = 2πx, where x = λ².

- **How it is computed:** Bouwkamp's tridiagonal Legendre method, then μ_n = √2 β₀/ψ_n(0) and λ_n = (c/2π) μ_n².
- **Checks:** λ₀(1) = 0.5725817806 and λ₀(4) = 0.9958854904, matching Slepian's tables. Connes' asymptotic (2¹⁴/3)√2 π⁵ e^{−4πx + 9/2 log x} matches the exact value with a ratio tending to 1 (0.45 at x = 2, 0.905 at 13, 0.958 at 30).
- **Speed-up:** the ζ archimedean terms now use closed forms in ψ and ψ′ plus geometric tails. They agree with the series to 1e-50, and the N = 100 build fell from 226 s to 1.8 s.

| x = λ² | 1 − χ₂ | ζ: ε₀ (N) | ζ: ε/(1 − χ₂) | DH: ε₀ (N = 64) | DH: ε/(1 − χ₂) |
|---|---|---|---|---|---|
| 3 | 8.7e-9 | 5.6e-8 (100) | 6.4 | 0.21 | 2.4e7 |
| 5 | 1.3e-18 | 1.0e-17 (100) | 7.7 | 6.0e-3 | 4.6e15 |
| 7 | 7.8e-29 | 7.7e-28 (100) | 9.9 | 8.6e-5 | 1.1e24 |
| 10 | 1.7e-44 | 1.6e-43 (180) | 9.2 | 8.3e-8 | 4.8e36 |
| 13 | 2.5e-60 | 3.0e-59 (180) | 11.9 | 8.6e-11 | 3.5e49 |
| 16 | 2.7e-76 | | | 4.0e-14 | 1.5e62 |
| 20 | 1.1e-97 | | | 4.1e-18 (N = 100) | 3.6e79 |
| 25 | 1.6e-124 | | | 1.6e-23 | 1.0e101 |
| 28 | 1.1e-140 | | | 8.4e-27 | 7.4e113 |
| 30 | 1.9e-151 | | | 1.3e-28 | 6.8e122 |

ε₀ can only fall as N grows, so each ζ ratio is an upper bound. At x = 13, ε₀ goes 3.72, 3.19, 2.96 (×1e-59) for N = 100, 140, 180, and extrapolating gives about 10. DH's ε₀ changes by about 10% between N = 64 and 100.

**Reading, corrected.** At first I divided DH's ε by ζ's floor, which made DH look far off. That was the wrong comparison. Each function has to be measured against its own archimedean floor:
- **ζ** (even Γ(s/2), conductor 1): 1 − χ₂ for h₄ at c = 2πx. The ratio is 6–12.
- **DH** (odd Γ((s+1)/2), conductor 5): 1 − χ₁ for the first odd prolate at c = 2πx/5 (`scripts/prolate_floor_dh.py`). The ratio is 3.4, 4.9, 6.0, 6.0, 7.7, 4.9, 8.9, 6.4, 5.3, 11.3 at x = 3 … 30, while ε_DH falls by 27 orders of magnitude.

So the minimum rides the floor set by each function's own Γ-factor and conductor, with or without an Euler product, until DH's sign change at x ≈ 30.5.

## Step 5: Connes' approximant k_λ against the minimiser θ_x (`scripts/k_lambda_mp.py`)

**Construction.**
- **ζ:** k_λ = E(h_λ), with E(f)(v) = v^{1/2} Σ f(nv) and h_λ = β₀⁽⁴⁾h_{0,λ} − β₀⁽⁰⁾h_{4,λ} (vanishing integral).
- **DH:** the analogue is k(u) = u^{1/2} Σ a_n h₁(nu), where h₁ is the first odd prolate function on [−λ/√5, λ/√5] (c = 2πx/5). The window is centred at the symmetry point 1/√5.
- **Quadrature** is split at the breakpoints v = λ/n, where a term switches on.

**ζ** (pole-free form, N = 100; at x = 13 the N = 140 values are 4.56e-59, 3.19e-59 and 4.2e-9):

| x | ε₀ | QW(k_λ)/‖k_λ‖² | ratio | sin²∠(k_λ, θ_x) | ‖k_odd‖/‖k‖ | first zero of k̂_λ |
|---|---|---|---|---|---|---|
| 5 | 1.0e-17 | 1.2e-17 | 1.19 | 2.0e-7 | 1.4e-9 | 7.2e-10 |
| 7 | 7.7e-28 | 1.1e-27 | 1.42 | 4.2e-8 | 1.1e-14 | 6.0e-15 |
| 10 | 1.8e-43 | 2.6e-43 | 1.48 | 2.0e-8 | 1.6e-22 | 1.7e-22 |
| 13 | 3.7e-59 | 5.1e-59 | 1.38 | 5.3e-9 | 2.2e-30 | 2.4e-30 |

**DH** (full form, N = 64):

| x | ε₀ | QW(k_λ)/‖k_λ‖² | ratio | sin²∠(k_λ, θ_x) |
|---|---|---|---|---|
| 5 | 6.0e-3 | 6.3e-3 | 1.06 | 1.3e-4 |
| 13 | 8.6e-11 | 9.1e-11 | 1.06 | 4.2e-8 |
| 20 | 4.5e-18 | 5.5e-18 | 1.24 | 2.7e-8 |
| 25 | 1.6e-23 | 2.9e-23 | 1.80 | 4.9e-9 |
| 30 | 1.3e-28 | 1.6e-28 | 1.25 | 3.8e-9 |
| 32 | −1.1e-29 | +1.0e-30 | — | 4.2e-7 |
| 34 | −1.3e-18 | +8.1e-33 | — | 1.0 |
| 40 | −0.011 | +3.7e-39 | — | 1.0 |

**What this shows.**
1. **For both functions, k_λ is close to the minimiser.** QW(k_λ) is within a factor of 1–2.3 of ε₀ and sin²∠ is about 1e-8. So the near-radical construction does not depend on the Euler product.
2. **k̂_λ's zero error equals its departure from evenness**, which is about √(1 − χ₂). k_λ reconstructs the zeros to (1 − χ₂)^{1/2}, while θ_x does much better (2.6e-55 at x = 13).
3. **The remaining step (§6.6) fails for DH exactly at its crossover.** k_λ^DH keeps sliding down its prolate floor and stays positive, because it is near-radical: its Fourier transform nearly vanishes at all of DH's zeros, including the off-line ones. The minimiser, by contrast, moves to a direction orthogonal to k_λ (sin² = 1) that exploits the off-line zero.
4. **So for ζ, "k_λ approximates θ_x for every λ" amounts to saying no vector ever undercuts the prolate floor.** That is Weil positivity. The data agree with it to sin² ≈ 4e-9 at x = 13, but they cannot settle it, since a crossover beyond the computed range cannot be ruled out numerically.

## Step 6: the whole low spectrum is prolate (`scripts/spectrum_prolate_mp.py`)

The lowest eigenvalues ε₀ < ε₁ < … of QW_λ were compared with the prolate leakages 1 − |χ_n| = (1 − λ_n)/(1 + √λ_n), split by residue class n mod 4. A residue class here corresponds to a Fourier eigenvalue: n ≡ 0 gives +1, n ≡ 1 gives −i, and so on.

| | x | ε_j/(1 − \|χ_n\|), matching class | Other classes |
|---|---|---|---|
| ζ (c = 2πx; n = 4, 8, 12, 16, 20) | 7 | 9.9, 4.4, 4.1, 3.0, 2.8 | 30 … 1.5e7 |
| ζ | 13 | 15.0, 7.0, 5.6, 5.1, 3.9 | 100 … 1.6e8 |
| DH (c = 2πx/5; n = 1, 5, 9, 13, 17) | 13 | 7.7, 4.1, 2.3, 1.5, 1.5 | 0.004 … 950 |
| DH | 25 | 6.4, 7.6, 3.5, 4.3, 2.3 | 7e-4 … 1.6e3 |

ζ's eigenvalues here are upper bounds at N = 100; for example ε₀ at x = 13 falls to 3.0e-59 at N = 180.

**Reading.** The j-th eigenvalue tracks the leakage of the j-th prolate in the function's own Fourier class:
- **ζ:** Fourier eigenvalue +1, starting at n = 4, because the near-radical vectors need ∫h = 0, so h₀ is used up by that condition.
- **DH:** odd, eigenvalue −i, starting at n = 1, with no pole and so no such condition.

The constant is between 1.5 and 15 while the eigenvalues span up to 27 decades.

This supports a picture in which the near-radical subspace of QW_λ is spanned by E_a(h_n) for the prolates h_n of that class, each contributing about its own leakage. The low spectrum of Weil's form is then archimedean, fixed by the Γ-factor and the conductor, whether or not there is an Euler product. What the Euler product decides is whether some other direction ever goes below this floor. For DH one does, at x ≈ 30.5.

## Step 7: what tips DH over

**Fine scan of DH's ε₀ between x = 30 and 31 (N = 64, 60 digits):**

| x | 30.0 | 30.25 | 30.5 | 30.75 | 30.95 | 30.999 | 31.0 |
|---|---|---|---|---|---|---|---|
| ε₀ | +1.29e-28 | +5.93e-29 | +1.81e-29 | +2.07e-30 | −2.99e-30 | −3.80e-30 | −3.81e-30 |

The sign change is continuous, at x_c ≈ 30.83, and it happens before n = 31 enters. A new grid point enters where the autocorrelation vanishes (F(L) = 0), so prime terms switch on continuously. What tips the form is the growth of the support, not any single prime.

**Reach in height.** The cosine basis resolves heights up to about 2πN/L. At x = 25, DH's ε₀ is 1.59e-23, 1.42e-23 and 1.31e-23 at N = 64, 120 and 200 (the last reaching height about 390). So DH's higher off-line zeros, 0.8695 + 240.40i and 0.8196 + 320.88i, do not make the form negative earlier. The low zero at 85.70 decides the crossover.

## Prime grid, torus and semilocal space

Kolossváry's prime grid ([arXiv:1711.02903](https://arxiv.org/abs/1711.02903)), extended to negative exponents, is the multiplicative group of the positive rationals, ≅ ⊕_p ℤ.

- **Duality.** Its Pontryagin dual is the prime torus, and the grid point N pairs with the angles θ_p = γ log p to give N^{iγ}.
- **Landau's formula on the grid.** The zeros' distribution on the torus has Fourier weight only on the grid's axes (prime powers), with weight Λ(n)/√n. DH puts weight off the axes as well.
- **Connes' Γ_S** is the grid restricted to a finite set of primes S, with signs. QW_λ uses the grid points n ≤ λ²: for ζ only those on the axes, for DH also those off them.

## Reproduce

```bash
cargo run --release --bin rh-connes-letter -- --function zeta --lambda-sq 5,7 --n 50
cargo run --release --bin rh-connes-letter -- --function dh --lambda-sq 13,23,30,40,50 --n 90 --t-max 125
.venv/bin/python scripts/connes_letter_mp.py --lambda-sq 13 --n 100 --dps 130 --variant nopole
.venv/bin/python scripts/connes_letter_mp.py --function dh --variant full --lambda-sq 31 --n 64 --dps 60 --zeros 0
.venv/bin/python scripts/prolate_ratio_mp.py --function zeta --x 3,5,7,10,13 --n 60,100 --dps 150
.venv/bin/python scripts/prolate_ratio_mp.py --function dh --x 3,5,7,10,13,16,20,25,28,30 --n 64 --dps 60 --prolate-dps 200
.venv/bin/python scripts/prolate_floor_dh.py
.venv/bin/python scripts/k_lambda_mp.py --x 5,7,10,13 --n 100 --dps 80 --zeros 12
.venv/bin/python scripts/k_lambda_mp.py --function dh --x 5,10,13,20,25,28,30,32,34,36,40 --n 64 --dps 60
.venv/bin/python scripts/spectrum_prolate_mp.py --function zeta --x 7,13 --n 100 --dps 130 --eigs 5
.venv/bin/python scripts/spectrum_prolate_mp.py --function dh --x 13,25 --n 64 --dps 60 --prolate-dps 150 --eigs 5
```

The `.venv` needs mpmath: `python3 -m venv .venv && .venv/bin/pip install mpmath`.
