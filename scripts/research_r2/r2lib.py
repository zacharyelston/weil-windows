"""R2 helpers: Poisson (E-map) trial functions for the window form and the unconditional bound of
docs/RESEARCH_R2_SAMPLING.md (Theorems 1-2).

Conventions (rh2's): window u in [-L/2, L/2], L = log x; F(t) = int f(u) e^{itu} du; the zeros-side form is
Q(f) = sum_rho F(gamma_rho) F(-gamma_rho), gamma_rho = -i(rho - 1/2) (= 2 sum_{gamma>0} |F(gamma)|^2 under RH).
Fourier on the line: psi_hat(y) = int psi(x) e^{-2 pi i x y} dx. Hermite functions h_k(y) = H_k(sqrt(2 pi) y) e^{-pi y^2}
satisfy h_k_hat = (-i)^k h_k.

E-map (self-dual variable w = e^u, conductor q, coefficients a_n = chi(n)):
    G(w) = w^{1/2} sum_{n>=1} a_n psi(n w / sqrt q).
Its Mellin transform int G(w) w^{it} dw/w = q^{(1/2+it)/2} L(1/2+it, chi) M_psi(t) vanishes at every nontrivial zero.
"""
import os
import sys

import mpmath as mp

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from decay_law_mp import kronecker  # noqa: E402


def hermite_coeffs(n):
    """Physicists' Hermite polynomial H_n as a coefficient list (ascending powers)."""
    a, b = [mp.mpf(1)], [mp.mpf(0), mp.mpf(2)]
    if n == 0:
        return a
    for k in range(1, n):
        c = [mp.mpf(0)] + [2 * v for v in b]          # H_{k+1} = 2z H_k - 2k H_{k-1}
        for i, v in enumerate(a):
            c[i] -= 2 * k * v
        a, b = b, c
    return b


def poly_add(p, q, s=1):
    out = [mp.mpf(0)] * max(len(p), len(q))
    for i, v in enumerate(p):
        out[i] += v
    for i, v in enumerate(q):
        out[i] += s * v
    return out


def poly_eval(p, z):
    r = mp.mpf(0)
    for v in reversed(p):
        r = r * z + v
    return r


def poly_deriv(p):
    return [i * p[i] for i in range(1, len(p))]


class HermitePsi:
    """psi(y) = P(z) e^{-pi y^2}, z = sqrt(2 pi) y, P a combination of Hermite polynomials of one class mod 4."""

    def __init__(self, combo):
        P = [mp.mpf(0)]
        for a, k in combo:
            P = poly_add(P, [mp.mpf(a) * v for v in hermite_coeffs(k)])
        self.P, self.dP = P, poly_deriv(P)
        assert len({k % 4 for _, k in combo}) == 1, "one Fourier class only"
        self.n = max(k for _, k in combo)
        self.combo = combo

    def __call__(self, y):
        z = mp.sqrt(2 * mp.pi) * y
        return poly_eval(self.P, z) * mp.exp(-mp.pi * y * y)

    def deriv(self, y):
        z = mp.sqrt(2 * mp.pi) * y
        return mp.sqrt(2 * mp.pi) * (poly_eval(self.dP, z) - z * poly_eval(self.P, z)) * mp.exp(-mp.pi * y * y)


def hermite_choice(kappa, s, pole):
    """Hermite trial psi for index n = kappa + 2s + 4*pole: h_n alone, or (pole) the combination of h_n, h_{n+4}
    with psi(0) = 0 (then psi_hat(0) = 0 too, as psi_hat = mu psi)."""
    n = kappa + 2 * s
    if not pole:
        return HermitePsi([(1, n)])
    return HermitePsi([(hermite_coeffs(n + 4)[0], n), (-hermite_coeffs(n)[0], n + 4)])


class Window:
    """f(u) = G(e^u), G(w) = w^{1/2} sum_n a_n psi(n w / sqrt q); psi decays like e^{-pi y^2}."""

    def __init__(self, psi, D, x, ycut=12):
        self.psi, self.D = psi, D
        self.q = 1 if D == 1 else abs(D)
        self.sq = mp.sqrt(self.q)
        self.x = mp.mpf(x)
        self.L = mp.log(self.x)
        self.a = (lambda n: 1) if D == 1 else (lambda n: kronecker(D, n))
        self.ycut = mp.mpf(ycut)

    def S(self, w, deriv=False):
        tot, n = mp.mpf(0), 1
        while True:
            y = n * w / self.sq
            if y > self.ycut:
                return tot
            a = self.a(n)
            if a:
                tot += a * ((n / self.sq) * self.psi.deriv(y) if deriv else self.psi(y))
            n += 1

    def f(self, u):
        w = mp.exp(u)
        return mp.sqrt(w) * self.S(w)

    def fprime(self, u):
        w = mp.exp(u)                     # d/du [w^{1/2} S(w)] = w^{1/2} (S/2 + w S')
        return mp.sqrt(w) * (self.S(w) / 2 + w * self.S(w, deriv=True))


def zeta_Nplus(t):
    """Upper bound for N(t) = #{0 < Im rho <= t}: 0 below the first ordinate 14.1347; else
    theta(t)/pi + 1 + 0.112 log t + 0.278 log log t + 2.510 (Riemann-von Mangoldt, Trudgian 2014 for S(t), t >= e)."""
    t = mp.mpf(t)
    if t < mp.mpf("14.134"):
        return mp.mpf(0)
    return mp.siegeltheta(t) / mp.pi + 1 + mp.mpf("0.112") * mp.log(t) + mp.mpf("0.278") * mp.log(mp.log(t)) + mp.mpf("2.510")


def dirichlet_Nplus(t, q, kappa):
    """Upper bound for N(t, chi) = #{|Im rho| <= t}, chi primitive mod q > 1 (Bennett-Martin-O'Bryant-Rechnitzer 2021:
    t >= 5/7 and l = log(q(t+2)/2pi) > 1.567 give |N - (t/pi) log(qt/2pi e) + chi(-1)/4| <= 0.22737 l + 2 log(1+l) - 0.5).
    N is nondecreasing, so for t below t_min the bound at t_min is used."""
    q = mp.mpf(q)
    tmin = max(mp.mpf(5) / 7, 2 * mp.pi * mp.exp(mp.mpf("1.5671")) / q - 2)
    t = max(mp.mpf(t), tmin)
    ell = mp.log(q * (t + 2) / (2 * mp.pi))
    chim1 = 1 if kappa == 0 else -1
    return (t / mp.pi * mp.log(q * t / (2 * mp.pi * mp.e)) - mp.mpf(chim1) / 4
            + mp.mpf("0.22737") * ell + 2 * mp.log(1 + ell) - mp.mpf("0.5"))


def zero_sum_bound(A, B, Nplus, zeta, t0=0):
    """Upper bound for sum over all nontrivial zeros of min(B^2, A^2/|Im rho|^2)
    = 2 A^2 int_{Tc}^inf N(t) t^{-3} dt with Tc = A/B (Stieltjes integration by parts; N(t) = 0 below t0).
    zeta=True: Nplus counts 0 < Im <= t, so the total is doubled."""
    lo = max(A / B, mp.mpf(t0))
    g = lambda t: Nplus(t) / t ** 3
    val = 2 * A * A * mp.quad(g, [lo, 2 * lo, 10 * lo, 1000 * lo, mp.inf])
    return 2 * val if zeta else val
