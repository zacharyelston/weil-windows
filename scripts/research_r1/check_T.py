"""Check the closed-form quadruple term against direct quadrature of Σ F(ρ)F(1−ρ), both sectors,
and reproduce the conductor-5 doc's real-pair contribution (+17.12, even, x = 7, f = (1 + cos θ)²)."""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import mpmath as mp
import r1lib as R

mp.mp.dps = 30
L = mp.log(7)
n = 5
worst = mp.mpf(0)
for parity in ("even", "odd"):
    dim = n + 1 if parity == "even" else n
    coeffs = [mp.mpf(k % 3 - 1) * mp.mpf("0.7") ** k + mp.mpf("0.3") for k in range(dim)]
    cv = mp.matrix(coeffs)
    f = R.basis_fn(coeffs, L, parity)
    for delta, gamma in [("0.3085", "85.699"), ("0.25", "0"), ("0", "14.1347"), ("0.1", "3.5"), ("0.45", "1.2")]:
        delta, gamma = mp.mpf(delta), mp.mpf(gamma)
        Tm, c, s, m = R.pair_matrix(delta, gamma, L, n, parity)
        closed = (cv.T * Tm * cv)[0]
        direct = R.T_direct(f, delta, gamma, L)
        err = abs(closed - direct)
        worst = max(worst, err / max(abs(direct), mp.mpf(1)))
        print(f"{parity:<4} δ={mp.nstr(delta,4):<6} γ₀={mp.nstr(gamma,7):<8} m={m}: closed {mp.nstr(closed, 15):>22}  direct {mp.nstr(direct, 15):>22}  |Δ| {mp.nstr(err, 3)}")
print("max relative error", mp.nstr(worst, 3))

# conductor-5 doc: real pair at δ = 1/4, even f = 3/2 + 2 cos θ + ½ cos 2θ at x = 7 contributes +17.12
A = [mp.mpf(3) / 2, mp.mpf(2), mp.mpf(1) / 2]
cv = mp.matrix([A[0] * mp.sqrt(L), A[1] * mp.sqrt(L / 2), A[2] * mp.sqrt(L / 2)] + [0] * (n - 2))
Tm, c, s, m = R.pair_matrix(mp.mpf(1) / 4, 0, L, n, "even")
print("F_t* real pair, even, x = 7, (1+cos θ)²:", mp.nstr((cv.T * Tm * cv)[0], 8), " (doc: +17.12)")
cv = mp.matrix([A[0] * mp.sqrt(mp.log(5)), A[1] * mp.sqrt(mp.log(5) / 2), A[2] * mp.sqrt(mp.log(5) / 2)] + [0] * (n - 2))
Tm, c, s, m = R.pair_matrix(mp.mpf(1) / 4, 0, mp.log(5), n, "even")
print("F_t* real pair, even, x = 5, (1+cos θ)²:", mp.nstr((cv.T * Tm * cv)[0], 8), " (doc: +11.69)")
