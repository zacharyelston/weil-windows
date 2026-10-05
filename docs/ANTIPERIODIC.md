# Twice Around: the Antiperiodic Basis

**Question.** rh2's finite bases on the window [−L/2, L/2] (L = log x) use periodic frequencies ω_k = 2πk/L, which close after one trip across the window. The half-integer frequencies ω_k = 2π(k + ½)/L close only after two trips (the "720°" lattice). Does the antiperiodic basis approximate Weil's minimiser better?

The two lattices impose opposite boundary behaviour:

| sector | periodic (ω = 2πk/L) | antiperiodic (ω = 2π(k+½)/L) |
|---|---|---|
| even | cos: zero slope at ±L/2 (Neumann) | cos: vanishes at ±L/2 (Dirichlet) |
| odd | sin: vanishes at ±L/2 (Dirichlet) | sin: zero slope at ±L/2 (Neumann) |

## Pre-registration (written before the sweep)

**Measure.** λ_N = λ_min of ζ's zeros-side form on the first N + 1 modes (Rayleigh–Ritz, so λ_N is an upper bound on the window's infimum and decreases towards it). Sweep: a ∈ {0.8, 1.0, 1.19, 1.3, 1.495} (x = e^{2a}), N ∈ {16, 32, 48, 64, 100, 140}, both sectors, both lattices.

**Prediction P1 (even sector).** The antiperiodic basis converges much faster: at each a, antiperiodic λ_N at N = 48 is at or below periodic λ_N at N = 140. Mechanism: if the even minimiser vanishes at the window edges, Dirichlet modes match it and the periodic cosines suffer the slow 1/k² coefficient decay of a kink-free but slope-mismatched expansion.

**Prediction P2 (odd sector).** The boundary matching reverses, so the antiperiodic basis is no better, and probably worse, in the odd sector.

**Kill.** If antiperiodic λ_N ≥ periodic λ_N at every N and every a in the even sector, the "twice around" basis buys nothing and P1 is dead.

**Comparison values** (certified lower bounds and best upper bounds from earlier runs, even sector):

| a | certified Q ≥ | best upper bound |
|---|---|---|
| 0.8 | 1.0277e-17 (Zhu, audited) | ≈ 1.634e-17 (limit) |
| 1.19 | 6.81311646960e-48 | 1.14e-47 |
| 1.3 | 5.76344479222e-62 | 1.04e-61 |
| 1.495 | (running) | 6.27e-96 at N = 180 |

## Validation (`scripts/antiperiodic_mp.py check`; output in `data/connes/antiperiodic/check.run.txt`)

- With shift 0, the generalised builder reproduces `connes_letter_mp.build_form` exactly (max |Δ| = 0 at x = 13, N = 12, both sectors).
- On the half-integer lattice the archimedean truncation terms change sign. The closed forms cut ∫_0^∞ at x = L, where e^{−zL} = e^{−aL} e^{iωL}, and cos(ωL) = −1 there instead of +1. The first version missed this and failed the zero-sum check (it even gave a negative diagonal entry). It has been fixed.
- Against 2Σ_ρ F(γ)² over the first 400 zeros, the agreement for edge-vanishing test functions (zero-sum tail ~1/T²) is:
  - antiperiodic even diagonals: 4.5e-7 to 4.8e-5;
  - antiperiodic odd b₀ + b₁ (off-diagonal): 1.6e-8 and 5.3e-9;
  - antiperiodic even b₁ − 2b₂ + b₃: 2.5e-6 and 1.3e-5;
  - periodic odd controls: 1e-6 to 2e-5.

  Edge-jumping functions converge only like log T/T, and they show the same few-percent shortfall on both lattices.

## Results (`data/connes/antiperiodic/`, M4, all 120 runs in < 10 min)

**P1 failed.** The table gives the ratio λ_anti/λ_per at each N, plus antiperiodic N = 48 against periodic N = 140:

| a | sector | 16 | 32 | 48 | 64 | 100 | 140 | anti 48 / per 140 |
|---|---|---|---|---|---|---|---|---|
| 0.8 | even | 0.91 | 1.03 | 0.94 | 0.96 | 1.01 | 1.04 | 1.06 |
| 0.8 | odd | 1.00 | 1.08 | 1.12 | 1.06 | 0.98 | 0.96 | 1.16 |
| 1.0 | even | 1.34 | 0.98 | 1.03 | 0.96 | 1.00 | 0.95 | 1.20 |
| 1.0 | odd | 1.23 | 0.99 | 1.04 | 1.01 | 1.01 | 1.06 | 1.31 |
| 1.19 | even | 0.34 | 0.66 | 1.10 | 1.03 | 0.95 | 0.98 | 2.28 |
| 1.19 | odd | 0.40 | 0.76 | 1.15 | 1.03 | 0.99 | 1.02 | 2.35 |
| 1.3 | even | 0.32 | 0.70 | 0.62 | 0.79 | 0.99 | 1.00 | 1.6e3 |
| 1.3 | odd | 0.32 | 0.78 | 0.63 | 0.82 | 1.03 | 0.98 | 8.5e2 |
| 1.495 | even | 0.25 | 0.19 | 0.81 | 0.21 | 0.69 | 1.00 | 1.8e21 |
| 1.495 | odd | 0.29 | 0.21 | 0.87 | 0.21 | 0.72 | 1.03 | 4.3e20 |

The kill as literally written did not fire, since the antiperiodic value is lower at many (a, N). The prediction it guarded did fail. While λ_N is still falling exponentially, the antiperiodic basis is ahead by a factor of 2–5. At a = 1.495, about 0.58 decades per mode, that is worth ≈ 1.2 modes. It is explained by reach: the half-integer lattice extends half a frequency step further, and in the odd sector it also has one more mode. Once both bases are converged, they agree to a few percent with no consistent winner, so the boundary condition doesn't matter.

**What this says.** The window minimum is set by the frequency cutoff alone. A basis is good to the extent that it reaches frequency T* = 2πx (k ≈ xL modes), and its boundary behaviour is irrelevant. Two facts make this consistent:
- the minimiser is already tiny at the edges (edge value 1.9e-8 at a = 0.8, `data/connes/zhu_window_ext.json`), so neither boundary condition is violated in any way that matters;
- the even and odd ratio rows above are nearly identical at each a, so both sectors feel the same cutoff.

The "twice around" reading gives nothing here. The pattern it pointed at, a rate set by a cutoff T* = 2πx, is taken up in `docs/DECAY_LAW.md`.
