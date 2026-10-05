# Weil positivity on finite windows across L-functions

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23163863.svg)](https://doi.org/10.5281/zenodo.23163863)

Code, data and documentation for the technical report **"Weil positivity on finite windows across L-functions: certified windows, a pre-registered decay law, and an unconditional bound for Dirichlet L-functions"** (Zac Elston, version 1, October 2026). The report is in [`docs/report/weil_windows_report_v1.pdf`](docs/report/weil_windows_report_v1.pdf), with its LaTeX source alongside.

For an L-function and x > 1, let λ_s(x) be the smallest value of Weil's quadratic form on unit-norm test functions supported in a window of length log x, in the even (s = 0) or odd (s = 1) sector. The Riemann hypothesis for that L-function is equivalent to λ_s(x) ≥ 0 for every x. Nothing here proves or disproves it. Upper bounds on λ also hold for functions with off-line zeros, and certificates exist only for finitely many windows.

## Results and their status

The labels are:
- **[C]** certified in Arb interval arithmetic, for the stated finite parameters;
- **[D]** derived: proved in the linked document, reviewed by AI agents but not yet by an external referee;
- **[E]** empirical: the outcome of a test registered before its data existed;
- **[N]** numerical: finite-basis Rayleigh–Ritz upper bounds.

| result | status | where |
|---|---|---|
| ζ positive on windows of support 2a = 2.38, 2.6, 2.99, both sectors; λ ≥ 3.501·10⁻⁹⁶ at 2.99 | [C] (2.38 audited; 2.6 and 2.99 regression-tested) | [`docs/CERTIFICATE_238.md`](docs/CERTIFICATE_238.md) |
| two-sided brackets for functions without an Euler product; ζ_K positive where its Epstein partner Z₁ is negative at x = e^{2.99} | [C] | [`docs/CONTROL_CERTIFICATES.md`](docs/CONTROL_CERTIFICATES.md), [`docs/EPSTEIN.md`](docs/EPSTEIN.md) |
| decay law λ ≈ e^{−2d·T*}, T* = 2π(x/Q)^{1/d}, across degrees 1–3, primitive and product, holomorphic and Maass type | [E], with every registered test and failure listed | [`docs/DECAY_LAW.md`](docs/DECAY_LAW.md), [`docs/GL2_DECAY.md`](docs/GL2_DECAY.md) |
| conductor collapse at fixed x/q, which discriminates against Zhu's Landau–Widom law transposed to Dirichlet L-functions | [N]/[E] | report §4 |
| matched pair: the Rédei L-function of (5, 29) and L(χ₅)L(χ₂₉), with the same invariants, have window minima within 0.47 nats | [E] | [`docs/BORROMEAN.md`](docs/BORROMEAN.md) |
| unconditional upper bound for real primitive Dirichlet L-functions, with the Hermite rate proved for all x; the Kaiser–Bessel rate e^{−2c} stays a conjecture for characters | [D] + [C] | [`docs/R2_THEOREMS.md`](docs/R2_THEOREMS.md), [`docs/R2_CONJECTURE_R.md`](docs/R2_CONJECTURE_R.md) |
| an off-line zero is excluded by a finite window only through a consistency test; reach forecast | [D] + [N] | [`docs/RESEARCH_R1_EXCLUSION.md`](docs/RESEARCH_R1_EXCLUSION.md) |

**Related and prior work.** All of the following is recent, unrefereed work that overlaps ours:
- Liu: certified supports 2.0 and 2.125.
- zeta-lab: support 2.38, even sector, on 28 September 2026.
- Arda: a conditional ζ_K against Epstein-zeta separation.
- Mori: ζ upper bounds including the rate e^{−4πx}, on Zenodo from 30 September 2026.
- Zhu (arXiv:2608.24827): the finite reduction, and the Landau–Widom conjecture for ζ.

[`docs/LIT_REVIEW_2026-10-04.md`](docs/LIT_REVIEW_2026-10-04.md) records what overlaps and what is new.

## Layout

- `docs/`: the report, and one document per result, each with its registrations, results and audits.
- `docs/viz/zero_torus.html`: an interactive page of ζ's zeros on the prime torus, the window bounds and the Borromean readout. It is live at **https://zacharyelston.github.io/weil-windows/viz/zero_torus.html**.
- `scripts/`: the code, Python with mpmath, python-flint (Arb) and cypari2 (PARI).
- `data/`: every committed result the documents quote.
- `docker/`: the pinned environment.

## Website

GitHub Pages serves `docs/` through `.github/workflows/pages.yml`:
- the landing page at https://zacharyelston.github.io/weil-windows/;
- the interactive visual at https://zacharyelston.github.io/weil-windows/viz/zero_torus.html;
- the report at https://zacharyelston.github.io/weil-windows/report/weil_windows_report_v1.pdf.

## Reproduce

```bash
python -m venv .venv && .venv/bin/pip install -r docker/requirements.txt
.venv/bin/python scripts/decay_claims_check.py        # regenerates every decay-law number quoted in docs/ from data/: expect 161/161
```

The pinned environment is Python 3.14, python-flint 0.9.0, mpmath 1.4.1, gmpy2 2.3.1 and cypari2 2.2.4 (PARI 2.17.2). Inside Docker:

```bash
docker build -t weil-windows-env:py3.14.4 -f docker/Dockerfile docker
docker/run.sh scripts/decay_claims_check.py
```

Each document names the commands that produced its data. Some runs are long. The certificate at support 2.99 takes several hours per sector, and the GL(2) scans take hours per object.

## Methods note

**Pre-registration.** Every empirical test was committed, with its predictions and kill criteria, before its data existed.

**Gates.** They compare quantities enclosed from opposite sides, and carry positive controls that must fire.

**Failures stay in the record:**
- the registered tests that failed or were killed are reported alongside the ones that held;
- one void run is kept with its explanation (`data/borromean/void_float_lambda/`).

**History.** This repository is a curated snapshot of a larger private research repository. Its history starts at the snapshot. References in `docs/` to `the private research repository` issues or commit hashes point into that private repository. The relevant registrations, reviews and results are reproduced in the documents here.

## AI assistance and responsibility

Computation, derivations, code and drafting were carried out with AI assistance: Claude Opus 5.5, with some work by Claude Fable 5.1 (Anthropic). Independent review and audit were done by GPT Astra and GPT Sol 6.1, and by separate Claude agents working read-only. The author is responsible for all content.

## License

Code (`scripts/`, `docker/`) is under the MIT License. Text, report and data (`docs/`, `data/`) are under CC BY 4.0. See [`LICENSE`](LICENSE).

## Citation

Zac Elston, *Weil positivity on finite windows across L-functions: certified windows, a pre-registered decay law, and an unconditional bound for Dirichlet L-functions*, technical report, Zenodo (2026), [doi:10.5281/zenodo.23163863](https://doi.org/10.5281/zenodo.23163863).

That is the concept DOI; it resolves to the latest version. Each release has its own DOI as well; see [`CHANGELOG.md`](CHANGELOG.md) and [`CITATION.cff`](CITATION.cff).
