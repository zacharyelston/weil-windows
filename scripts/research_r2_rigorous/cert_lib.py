"""Certified (Arb ball) evaluation of R2's unconditional upper bound (docs/R2_THEOREMS.md, Theorem 2).

For a trial psi in class A_chi(mu) (Theorem 1), f(u) = e^{u/2} S(e^u), S(w) = sum_{n>=1} chi(n) psi(n w / sqrt q), and

    lambda_s(x) <= Q(f_win)/||f_win||^2 <= 2 A^2 int_{max(A/B, t_low)}^inf N+(t) t^-3 dt / ||f_win||^2,

where B >= 2 int_{sqrt x}^inf |S(w)| dw, A >= 2 (sqrt x |S(sqrt x)| + int_{sqrt x}^inf |S/2 + w S'| dw),
||f_win||^2 = 2 int_1^{sqrt x} S(w)^2 dw, and N+ >= N(t) = #{rho : |Im rho| <= t} (with multiplicity).
Every quantity below is an Arb ball; the bound is the upper endpoint of the final ball. Nothing in the chain uses a
floating midpoint: floats only choose split points, truncation orders and quadrature endpoints, each of which is then
used as an exact dyadic number (any choice gives a valid bound).

Inequality directions (the whole point of this file):
    A, B, the zero-count integral   -> enclosed, the UPPER endpoint is used;
    ||f_win||^2                     -> enclosed, the LOWER endpoint is used.

Zero counts (docs/R2_THEOREMS.md section 2, read in source):
    zeta   : |N(T) - (T/2pi) log(T/2pi e)| <= 0.11200 log T + 0.12567 log log T + 3.77417, T >= e
             (Bellotti-Wong, arXiv:2412.15470v2, Theorem 1.1, second estimate; N(T) = #{0 < gamma <= T});
             N(t) = 0 for t < 14 (first ordinate 14.1347...).
    L(s,chi): chi of conductor q > 1, T >= 5/7, l = log(q(T+2)/2pi): N(T,chi) = 0 if l <= 1.567, else
             |N(T,chi) - (T/pi) log(qT/2pi e) + chi(-1)/4| <= 0.22737 l + 2 log(1+l) - 0.5
             (Bennett-Martin-O'Bryant-Rechnitzer, arXiv:2005.02989, Theorem 1.1; N(T,chi) = #{|gamma| <= T}).
"""
import math
import os
import sys

from flint import acb, arb, ctx

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from decay_law_mp import kronecker  # noqa: E402

PREC = 256


def setprec(p=PREC):
    ctx.prec = p


def chi_fn(D):
    return (lambda n: 1) if D == 1 else (lambda n: kronecker(D, n))


def exact(v):
    """An exact dyadic arb from a float (used only for split points and endpoints)."""
    return arb(float(v))


def up(b):
    """Upper endpoint of a ball as an exact arb."""
    return arb(b.upper())


def lo(b):
    return arb(b.lower())


def certainly_pos(b):
    return b > 0


# ------------------------------------------------------------------------------------------- zero counts

BW_ZETA = (arb("0.11200"), arb("0.12567"), arb("3.77417"))
BMOR = (arb("0.22737"), arb("1.567"))
ZETA_T0 = arb(14)          # N(t) = 0 below the first ordinate 14.1347...


def zeta_Nplus(t):
    """Upper bound for N(t) = #{rho : |Im rho| <= t} = 2 N_zeta(t), t >= e (Bellotti-Wong Thm 1.1, 2nd estimate)."""
    a, b, c3 = BW_ZETA
    t = arb(t)
    return 2 * (t / (2 * arb.pi()) * (t / (2 * arb.pi() * arb(1).exp())).log() + a * t.log() + b * t.log().log() + c3)


def zeta_tail_integral(T):
    """Upper bound (ball; use its upper endpoint) for int_T^inf N(t) t^-3 dt, N(t) = 2 N_zeta(t), T >= 14 (> e).
    Closed forms:  int_T^inf (t/2pi) log(t/2pi e) t^-3 = log(T/2pi)/(2 pi T);  int a log t t^-3 = a (2 log T + 1)/(4T^2);
    int b loglog t t^-3 <= b (loglog T/(2T^2) + 1/(4 T^2 log T)) (by parts, b >= 0);  int c t^-3 = c/(2T^2)."""
    a, b, c3 = BW_ZETA
    T = arb(T)
    assert T >= ZETA_T0
    pi = arb.pi()
    lT = T.log()
    main = (T / (2 * pi)).log() / (2 * pi * T)
    ta = a * (2 * lT + 1) / (4 * T * T)
    tb = b * (lT.log() / (2 * T * T) + 1 / (4 * T * T * lT))
    tc = c3 / (2 * T * T)
    return 2 * (main + ta + tb + tc)


def dirichlet_tlow(q):
    """t_low = max(5/7, 2 pi e^{1.567}/q - 2): N(t, chi) = 0 for t <= t_low when t_low > 5/7 (BMOR Thm 1.1), and the
    BMOR inequality holds for every t > t_low."""
    t1 = 2 * arb.pi() * BMOR[1].exp() / q - 2
    t57 = arb(5) / 7
    return up(t1) if t1 > t57 else up(t57)          # exact dyadic >= the true t_low


def dirichlet_Nplus(t, q, chim1):
    """BMOR upper bound for N(t, chi) = #{|gamma| <= t}, valid for t > t_low (where l > 1.567)."""
    c, _ = BMOR
    t = arb(t)
    pi = arb.pi()
    ell = (q * (t + 2) / (2 * pi)).log()
    return t / pi * (q * t / (2 * pi * arb(1).exp())).log() - arb(chim1) / 4 + c * ell + 2 * (1 + ell).log() - arb("0.5")


def dirichlet_tail_integral(T, q, chim1):
    """Upper bound for int_T^inf N+(t, chi) t^-3 dt, T >= t_low.  Closed forms:
       int_T^inf (t/pi) log(qt/2pi e) t^-3 = log(qT/2pi)/(pi T);
       int_T^inf l(t) t^-3 = l(T)/(2T^2) + 1/(4T) - log(1+2/T)/8      (l = log(q/2pi) + log(t+2), by parts);
       log(1+l) <= log(1+l0) + (l-l0)/(1+l0), l0 = l(T) (concavity), then the line above."""
    c, _ = BMOR
    T = arb(T)
    assert T >= dirichlet_tlow(q)
    pi = arb.pi()
    l0 = (q * (T + 2) / (2 * pi)).log()
    main = (q * T / (2 * pi)).log() / (pi * T)
    const = (-arb(chim1) / 4 - arb("0.5")) / (2 * T * T)
    Iell = l0 / (2 * T * T) + 1 / (4 * T) - (1 + 2 / T).log() / 8
    tl = c * Iell
    tlog = 2 * (((1 + l0).log() - l0 / (1 + l0)) / (2 * T * T) + Iell / (1 + l0))
    return main + const + tl + tlog


def zero_sum_bound(A, B, D):
    """Upper bound (ball) for sum_rho min(B^2, A^2/|Im rho|^2), from A, B upper bounds (exact upper endpoints taken):
    = 2 A^2 int_{A/B}^inf N(t) t^-3 dt <= 2 A^2 int_{T_L}^inf N+(t) t^-3 dt, T_L = max(A/B, t_low) (N = 0 below t_low)."""
    A, B = up(A), up(B)
    Tc = A / B
    # sum <= 2A^2 int_{max(Tc, t_low)}^inf N t^-3 <= 2A^2 int_{max(lo(Tc), t_low)}^inf N t^-3 (N >= 0), then N <= N+
    if D == 1:
        tl = ZETA_T0
        TL = lo(Tc) if lo(Tc) >= tl else tl
        I = zeta_tail_integral(TL)
    else:
        q = abs(D)
        chim1 = 1 if D > 0 else -1
        tl = dirichlet_tlow(q)
        assert lo(Tc) > arb(5) / 7, "T_c below 5/7: not covered"
        TL = lo(Tc) if lo(Tc) >= tl else tl
        I = dirichlet_tail_integral(TL, q, chim1)
    Z = 2 * A * A * I
    return Z, Tc, TL


# ------------------------------------------------------------------------------------------- tails

def gauss_tail(C, d, alpha, N):
    """Upper bound for sum_{n > N} C n^d e^{-alpha n^2} (C, alpha > 0): geometric with ratio
    r = ((N+2)/(N+1))^d e^{-alpha (2N+3)} >= t_{n+1}/t_n for all n >= N+1."""
    N = int(N)
    r = (arb(N + 2) / (N + 1)) ** d * (-alpha * (2 * N + 3)).exp()
    assert r < 1, "gauss_tail ratio not < 1"
    return up(C * arb(N + 1) ** d * (-alpha * (N + 1) ** 2).exp() / (1 - r))


# ------------------------------------------------------------------------------------------- Hermite trials

def herm_coeffs(n):
    """Physicists' Hermite polynomial H_n, integer coefficients, ascending powers."""
    a, b = [1], [0, 2]
    if n == 0:
        return a
    for k in range(1, n):
        c = [0] + [2 * v for v in b]
        for i, v in enumerate(a):
            c[i] -= 2 * k * v
        a, b = b, c
    return b


def poly_add(p, q, s=1):
    out = [0] * max(len(p), len(q))
    for i, v in enumerate(p):
        out[i] += v
    for i, v in enumerate(q):
        out[i] += s * v
    while len(out) > 1 and out[-1] == 0:
        out.pop()
    return out


def poly_mul_z(p):
    return [0] + list(p)


def poly_deriv(p):
    return [i * p[i] for i in range(1, len(p))] or [0]


def poly_eval(p, z):
    r = 0 * z
    for v in reversed(p):
        r = r * z + v
    return r


class HermiteTrial:
    """psi(y) = P(z) e^{-z^2/2}, z = sqrt(2 pi) y; P = H_n (no pole) or H_{n+4}(0) H_n - H_n(0) H_{n+4} (zeta, psi(0) = 0),
    n = kappa + 2s. psi_hat = mu psi with mu = (-i)^n (h_k_hat = (-i)^k h_k)."""

    def __init__(self, kappa, s, pole):
        n = kappa + 2 * s
        if pole:
            P = poly_add([herm_coeffs(n + 4)[0] * v for v in herm_coeffs(n)], [herm_coeffs(n)[0] * v for v in herm_coeffs(n + 4)], -1)
            assert P[0] == 0, "psi(0) != 0"
            self.index = n + 4
        else:
            P = herm_coeffs(n)
            self.index = n
        self.P = P
        self.d = len(P) - 1
        # y psi'(y) as a polynomial in z times e^{-z^2/2}: d/dy = sqrt(2 pi) d/dz, so y psi' = z (P' - z P) e^{-z^2/2}
        R = poly_add(poly_deriv(P), poly_mul_z(P), -1)
        self.Q = poly_mul_z(R)
        self.CP = sum(abs(v) * (2 * arb.pi()) ** (arb(k) / 2) for k, v in enumerate(P) if v)
        self.CQ = sum(abs(v) * (2 * arb.pi()) ** (arb(k) / 2) for k, v in enumerate(self.Q) if v)
        self.dQ = len(self.Q) - 1

    def psi(self, y):
        s2p = (2 * arb.pi()).sqrt()
        return poly_eval(self.P, s2p * y) * (-arb.pi() * y * y).exp()


def gk(k, Z):
    """int_Z^inf z^k e^{-z^2/2} dz = 2^{(k-1)/2} Gamma((k+1)/2, Z^2/2)."""
    return arb(2) ** (arb(k - 1) / 2) * (Z * Z / 2).gamma_upper(arb(k + 1) / 2)


def taylor_shift(P, Z):
    from math import comb
    return [sum(P[j] * comb(j, k) * Z ** (j - k) for j in range(k, len(P)) if P[j]) for k in range(len(P))]


def int_abs_poly_gauss_y(P, Y):
    """int_Y^inf |P(sqrt(2 pi) y)| e^{-pi y^2} dy, Y > 0 (ball). If P(Z + t) has Taylor coefficients of one strict sign
    (certified), P has no root on [Z, inf) and the integral is |sum p_k G_k(Z)|/sqrt(2 pi); otherwise the triangle bound."""
    s2p = (2 * arb.pi()).sqrt()
    Z = s2p * Y
    d = taylor_shift(P, Z)
    signed = (d[0] > 0 and all(v >= 0 for v in d)) or (d[0] < 0 and all(v <= 0 for v in d))
    if signed:
        return abs(sum((p * gk(k, Z) for k, p in enumerate(P) if p), arb(0))) / s2p, True
    return sum((abs(p) * gk(k, Z) for k, p in enumerate(P) if p), arb(0)) / s2p, False


def hermite_quantities(trial, D, x, s):
    """Certified A, B (upper) and ||f_win||^2 (lower) for the Hermite E-map trial. Returns dict of balls."""
    chi = chi_fn(D)
    q = 1 if D == 1 else abs(D)
    x = arb(x)
    sq = arb(q).sqrt()
    lam = (x / q).sqrt()
    pi = arb.pi()
    d = trial.d
    # ---- B and the A-integral: sums over n <= N of exact integrals, Gaussian tail beyond
    N = int(math.ceil(12 / float(lam.lower()))) + 2
    sB, sA, all_signed = arb(0), arb(0), True
    for n in range(1, N + 1):
        if not chi(n):
            continue
        i0, ok0 = int_abs_poly_gauss_y(trial.P, n * lam)
        i1, ok1 = int_abs_poly_gauss_y(trial.Q, n * lam)
        all_signed &= ok0 and ok1
        sB += i0 / n
        sA += i1 / n
    # I(Y) <= C int_Y^inf y^d e^{-pi y^2} <= (C/pi) Y^{d-1} e^{-pi Y^2} <= (C/pi) Y^d e^{-pi Y^2} for Y >= 1, 2 pi Y^2 >= 2(d-1)
    assert float((N * lam).lower()) >= 12
    tB = gauss_tail(trial.CP * lam ** d / pi, d, pi * lam * lam, N)
    tA = gauss_tail(trial.CQ * lam ** trial.dQ / pi, trial.dQ, pi * lam * lam, N)
    B = 2 * sq * (sB + arb(0, tB))
    # ---- edge value S(sqrt x) = sum chi(n) psi(n lam)
    Se = sum((chi(n) * trial.psi(n * lam) for n in range(1, N + 1) if chi(n)), arb(0))
    te = gauss_tail(trial.CP * lam ** d, d, pi * lam * lam, N)
    Se = Se + arb(0, te)
    edge = x.sqrt() * abs(Se)
    Aint = up(B) / 4 + sq * (sA + arb(0, tA))
    A = 2 * (edge + Aint)
    # ---- norm: 2 int_1^{X} S_N(w)^2 dw, X <= sqrt x exact; tail eps_S = sup_{w >= 1} |S - S_N|
    Ns = int(math.ceil(9 * math.sqrt(q))) + 1
    epsS = gauss_tail(trial.CP / arb(q) ** (arb(d) / 2), d, pi / q, Ns)
    X = lo(x.sqrt())
    coeffs = [(n, chi(n)) for n in range(1, Ns + 1) if chi(n)]
    s2p = (2 * pi).sqrt()

    def SN(w, _an):
        tot = acb(0)
        for n, c in coeffs:
            y = n * w / sq
            tot += c * poly_eval(trial.P, s2p * y) * (-pi * y * y).exp()
        return tot * tot

    I = acb.integral(SN, 1, X, rel_tol=arb(2) ** -80)
    nS2 = I.real
    nS = lo(nS2).sqrt() if nS2 > 0 else None
    assert nS is not None
    nlow = nS - epsS * (X - 1).sqrt()
    assert nlow > 0
    norm2 = 2 * nlow * nlow
    return {"A": A, "B": B, "norm2": norm2, "edge": edge, "Aint": Aint, "S_edge": Se, "int_SN2": nS2,
            "eps_S": arb(epsS), "signed_integrals": bool(all_signed), "n_terms": N}


def hermite_symmetry(trial, D, x, s, w0):
    """Theorem 1(i): E(1/w) = sigma E(w), E(w) = w^{1/2} S(w). Returns the ball E(1/w0) - sigma E(w0) and |E(w0)|."""
    chi = chi_fn(D)
    q = 1 if D == 1 else abs(D)
    sq = arb(q).sqrt()
    sigma = 1 if s == 0 else -1
    pi = arb.pi()

    def E(w):
        # terms with n w / sqrt q <= 12 exactly, Gaussian tail beyond (|psi(y)| <= C_P y^d e^{-pi y^2}, y >= 1)
        N = int(math.ceil(12 * math.sqrt(q) / float(w.lower()))) + 2
        tot = sum((chi(n) * trial.psi(n * w / sq) for n in range(1, N + 1) if chi(n)), arb(0))
        t = gauss_tail(trial.CP * (w / sq) ** trial.d, trial.d, pi * w * w / q, N)
        return w.sqrt() * (tot + arb(0, t))

    w0 = arb(w0)
    e1, e2 = E(w0), E(1 / w0)
    return e2 - sigma * e1, e1


# ------------------------------------------------------------------------------------------- Kaiser-Bessel trials

def a_hank(k, j):
    return arb(math.factorial(k + j)) / (math.factorial(j) * math.factorial(k - j) * 2 ** j)


def cap_K(k):
    """|K_{k+1/2}(omega)| = |J_nu(s)/s^nu| <= 2^-nu / Gamma(nu+1) for real s >= 0 (DLMF 10.14.4)."""
    nu = arb(2 * k + 1) / 2
    return arb(2) ** (-nu) / (nu + 1).gamma()


def hank_K(k, s):
    """|K_{k+1/2}(omega)| <= sqrt(2/pi) sum_{j<=k} a_j(k) s^{-k-1-j}, s = sqrt(omega^2 - beta^2) > 0
    (DLMF 10.49.1-2: terminating expansion of j_k; J_{k+1/2}(s) = sqrt(2s/pi) j_k(s))."""
    return (2 / arb.pi()).sqrt() * sum((a_hank(k, j) * s ** (-(k + 1 + j)) for j in range(k + 1)), arb(0) * s)


class KBTrial:
    """phi on [-lam, lam] from the Kaiser-Bessel window w(t) = (1 - t^2)^{m/2} I_m(beta sqrt(1 - t^2)):
         kappa = 0: phi(y) = (1 + a t^2) w(t), t = y/lam (a fixed by psi(0) = 0 for zeta, else 0);  kappa = 1: phi = t w(t).
       phi_hat(xi) = pref * sum_i c_i omega^{p_i} K_{nu + j_i}(omega) (times -i for kappa = 1), omega = 2 pi lam xi,
       pref = lam sqrt(2 pi) beta^m, nu = m + 1/2, K_mu(omega) = 2^-mu 0F1~(; mu+1; (beta^2 - omega^2)/4)  [derived in
       docs/R2_THEOREMS.md, Lemma K]. psi = (phi + conj(mu) phi_hat)/2 = (phi + cp * phihat_r)/2, cp = (-1)^s."""

    def __init__(self, m, beta, lam, kappa, s, pole):
        self.m, self.beta, self.lam, self.kappa, self.s = m, arb(beta), arb(lam), kappa, s
        self.nu = arb(2 * m + 1) / 2
        self.cp = 1 if s == 0 else -1
        self.pref = self.lam * (2 * arb.pi()).sqrt() * self.beta ** m
        self.k0 = m                      # nu + j = (m + j) + 1/2
        if kappa == 0:
            self.a = arb(0)
            if pole:
                Im = (self.beta / 2) ** m * (self.beta * self.beta / 4).hypgeom_0f1(m + 1, regularized=True)
                K0, K1 = self.K(0, arb(0)), self.K(1, arb(0))
                self.a = -(Im + self.cp * self.pref * K0) / (self.cp * self.pref * K1)
            self.terms = [(arb(1), 0, 0), (self.a, 0, 1), (-self.a, 2, 2)] if pole else [(arb(1), 0, 0)]
        else:
            self.a = arb(0)
            self.terms = [(arb(1), 1, 1)]

    def K(self, j, omega, an=False):
        """K_mu(omega) = 2^-mu 0F1~(; mu+1; (beta^2 - omega^2)/4), entire in omega. Out of band (Re(omega^2 - beta^2) > 4)
        it is evaluated as J_mu(s)/s^mu, s = sqrt(omega^2 - beta^2) (Arb's Bessel J switches to the asymptotic
        expansion; the 0F1 series loses all precision near (beta^2 - omega^2)/4 = -10^4). Both are enclosures of the same
        entire function (principal branches agree for Re s > 0)."""
        mu = self.nu + j
        d2 = omega * omega - self.beta * self.beta
        d2r = d2.real if isinstance(d2, acb) else d2
        if float(d2r.mid()) > 4:
            sq = d2.sqrt(analytic=an) if isinstance(d2, acb) else d2.sqrt()
            return sq.bessel_j(mu) / sq ** mu
        return (-d2 / 4).hypgeom_0f1(mu + 1, regularized=True) / arb(2) ** mu

    def phihat_r(self, y, an=False):
        om = 2 * arb.pi() * self.lam * y
        return self.pref * sum((c * om ** p * self.K(j, om, an) for c, p, j in self.terms if not (c == 0)), 0 * om)

    def phi_ent(self, y):
        t = y / self.lam
        u = 1 - t * t
        w = (self.beta / 2) ** self.m * u ** self.m * (self.beta * self.beta * u / 4).hypgeom_0f1(self.m + 1, regularized=True)
        return (1 + self.a * t * t) * w if self.kappa == 0 else t * w

    def psi_full(self, y, an=False):
        return (self.phi_ent(y) + self.cp * self.phihat_r(y, an)) / 2

    def psi_out(self, y, an=False):
        return self.cp * self.phihat_r(y, an) / 2

    def psi0(self):
        return self.psi_full(arb(0))

    # envelope pieces: E(omega) = sum_i |c_i| omega^{p_i} env_{k0 + j_i}(omega)  >=  |sum_i c_i omega^p_i K_{nu+j_i}(omega)|
    def env_terms(self):
        return [(abs(c), p, self.k0 + j) for c, p, j in self.terms if not (c == 0)]

    def deriv_env_terms(self):
        """omega * E_1(omega) as sum |c| omega^P env_K: d/domega(omega^p K_mu) = p omega^{p-1} K_mu - omega^{p+1} K_{mu+1}."""
        out = []
        for c, p, k in self.env_terms():
            if p:
                out.append((c * p, p, k))
            out.append((c, p + 2, k + 1))
        return out


def env_G(p, k, beta, c, Omega_factor=64):
    """Upper bound (ball) for int_c^inf omega^p env_k(omega) (1 + log(omega/c)) d omega, env_k = min(cap, hank).
    cap on [c, omega*], hank on [omega*, Omega] (acb.integral) and on [Omega, inf) with s >= r omega (closed form)."""
    pi = arb.pi()
    capk = cap_K(k)
    bf, cf = float(beta.mid()), float(c.mid())
    # crossover s*: hank(s*) = cap (floating, only a split point)
    capf = float(capk.mid())

    def hankf(s):
        return math.sqrt(2 / math.pi) * sum(float(a_hank(k, j).mid()) * s ** (-(k + 1 + j)) for j in range(k + 1))
    slo, shi = 1e-12, 1e6
    for _ in range(200):
        sm = math.sqrt(slo * shi)
        if hankf(sm) > capf:
            slo = sm
        else:
            shi = sm
    om_star = exact(math.sqrt(shi * shi + bf * bf) * (1 + 1e-12))
    if not (om_star >= up(c)):
        om_star = exact(float(up(c).mid()) * (1 + 1e-9))
    assert om_star >= up(c) and om_star > beta, "split point must lie above c and beta"

    def F(om):
        return om ** (p + 1) / (p + 1) * (1 + (om / c).log()) - om ** (p + 1) / (p + 1) ** 2
    part_cap = capk * (F(om_star) - F(c))          # exact integral over [c, omega*] (F' >= 0 there), times cap
    Omega = exact(max(Omega_factor * cf, 2 * float(om_star.mid())))
    sq2pi = (2 / pi).sqrt()
    coeffs = [a_hank(k, j) for j in range(k + 1)]

    def integrand(om, an):
        s = (om * om - beta * beta).sqrt(analytic=an)
        h = sum((cf_ * s ** (-(k + 1 + j)) for j, cf_ in enumerate(coeffs)), acb(0))
        return om ** p * sq2pi * h * (1 + (om / c).log(analytic=an))
    mid = acb.integral(integrand, om_star, Omega, rel_tol=arb(2) ** -60).real
    r = (1 - beta * beta / (Omega * Omega)).sqrt()
    tail = arb(0)
    for j, cf_ in enumerate(coeffs):
        e = k + 1 + j
        ep = e - p
        assert ep > 1
        tail += sq2pi * cf_ * r ** (-e) * Omega ** (1 - ep) * ((1 + (Omega / c).log()) / (ep - 1) + arb(1) / (ep - 1) ** 2)
    return part_cap + mid + tail, om_star


def C_E(terms, beta, Omega0):
    """omega^3 E(omega) <= C_E for omega >= Omega0 > beta (s >= r0 omega, r0 = sqrt(1 - beta^2/Omega0^2))."""
    r0 = (1 - beta * beta / (Omega0 * Omega0)).sqrt()
    tot = arb(0)
    for c, p, k in terms:
        for j in range(k + 1):
            ex = p - k - 1 - j + 3
            assert ex <= 0
            tot += c * (2 / arb.pi()).sqrt() * a_hank(k, j) * r0 ** (-(k + 1 + j)) * Omega0 ** ex
    return tot


M_LADDER = (0, 3, 8, 20, 50)


def kb_quantities(trial, D, x, s, N_edge=200):
    """Certified A, B (upper) and ||f_win||^2 (lower) for the Kaiser-Bessel E-map trial."""
    chi = chi_fn(D)
    q = 1 if D == 1 else abs(D)
    x = arb(x)
    sq = arb(q).sqrt()
    lam = trial.lam
    c = 2 * arb.pi() * x / q
    beta, m = trial.beta, trial.m
    # beta = c - dbeta with dbeta >= 0 (checked by the caller): out of band means omega >= c >= beta
    pi = arb.pi()
    terms = trial.env_terms()
    dterms = trial.deriv_env_terms()
    # ---- J0 = int_c^inf E (1 + log(omega/c)), J1 = int_c^inf omega E_1 (1 + log(omega/c))
    J0 = sum((cc * env_G(p, k, beta, c)[0] for cc, p, k in terms), arb(0))
    J1 = sum((cc * env_G(p, k, beta, c)[0] for cc, p, k in dterms), arb(0))
    pre = sq * beta ** m / (2 * pi).sqrt()
    B = pre * J0                                   # = sqrt q pref/(2 pi lam) J0
    # ---- edge S(sqrt x) = sum_n chi(n) psi_out(n lam), explicit to N_edge, tail (pref/2) C_E c^-3 /(2 N^2)
    Se = sum((chi(n) * trial.psi_out(n * lam) for n in range(1, N_edge + 1) if chi(n)), arb(0))
    Om0 = (N_edge + 1) * c
    te = trial.pref / 2 * C_E(terms, beta, Om0) / c ** 3 / (2 * N_edge ** 2)
    Se = Se + arb(0, up(te))
    edge = x.sqrt() * abs(Se)
    Aint = up(B) / 4 + pre / 2 * J1                # = B/4 + sqrt q pref/(4 pi lam) J1
    A = 2 * (edge + Aint)
    # ---- norm: intervals I_k = [sqrt x/(k+1), sqrt x/k] cap [1, sqrt x], S_in = sum_{n<=k} chi(n) psi_full(n w/sqrt q)
    rx = x.sqrt()
    K = int(math.floor(float(rx.mid())))
    if arb(K) > rx:
        K -= 1
    # On [a_k, b_k] (inside I_k) S = S_main + S_rest, S_main = sum_{n<=k} chi psi_full + sum_{k<n<=k+M} chi psi_out
    # (finite sum of entire functions: acb.integral), |S_rest| <= (pref/2) C_E(Om) (c a_k/sqrt x)^-3 /(2 (k+M)^2)
    # with Om = c (k+M+1) a_k/sqrt x (envelope E(omega) <= C_E omega^-3 for omega >= Om).  Per interval
    # ||S|| >= ||S_main|| - sup|S_rest| sqrt|I|, and ||S||^2_total = sum over intervals.
    nlow2 = arb(0)
    pieces = []
    for k in range(1, K + 1):
        lb = rx / (k + 1)
        assert (lb > 1) or (lb <= 1), "band breakpoint straddles w = 1"
        a_ = up(lb) if lb > 1 else arb(1)
        b_ = lo(rx / k)
        if not (b_ > a_):
            continue
        cin = [(n, chi(n)) for n in range(1, k + 1) if chi(n)]
        for Mk in M_LADDER:                    # smallest M whose remainder costs < 1e-6 of the interval norm
            cout = [(n, chi(n)) for n in range(k + 1, k + Mk + 1) if chi(n)]

            def Smain2(w, an, cin=cin, cout=cout):
                tot = acb(0)
                for n, cc in cin:
                    tot += cc * trial.psi_full(n * w / sq, an)
                for n, cc in cout:
                    tot += cc * trial.psi_out(n * w / sq, an)
                return tot * tot
            I = acb.integral(Smain2, a_, b_, rel_tol=arb(2) ** -60).real
            Om = c * (k + Mk + 1) * a_ / rx
            eps = trial.pref / 2 * C_E(terms, beta, Om) * (c * a_ / rx) ** -3 / (2 * arb(max(k + Mk, 1)) ** 2)
            if up(eps) * (b_ - a_).sqrt() <= arb("1e-6") * lo(I).sqrt():
                break
        nk = lo(I).sqrt() - up(eps) * (b_ - a_).sqrt()
        assert nk > 0, ("interval norm not positive", k)
        nlow2 += nk * nk
        pieces.append({"k": k, "M": Mk, "a": a_.str(12), "b": b_.str(12), "int_Smain2": I.str(12), "eps_rest": up(eps).str(4)})
    nlow = lo(nlow2).sqrt()
    assert nlow > 0
    norm2 = 2 * nlow * nlow
    return {"A": A, "B": B, "norm2": norm2, "edge": edge, "Aint": Aint, "S_edge": Se, "J0": J0, "J1": J1,
            "pieces": pieces, "c": c}


def kb_S(trial, D, x, w, N_out=400):
    """S(w) = sum_n chi(n) psi(n w/sqrt q) for real w > 0 (ball): in-band n < sqrt x/w exactly, out-of-band n up to
    N_out exactly, tail (pref/2) C_E sum (c n w/sqrt x)^-3."""
    chi = chi_fn(D)
    q = 1 if D == 1 else abs(D)
    x = arb(x)
    sq = arb(q).sqrt()
    rx = x.sqrt()
    c = 2 * arb.pi() * x / q
    w = arb(w)
    nin = int(math.floor(float((rx / w).mid())))
    tot = arb(0)
    for n in range(1, max(nin, 0) + N_out + 1):
        if not chi(n):
            continue
        y = n * w / sq
        if (y < trial.lam):
            tot += chi(n) * trial.psi_full(y)
        elif (y >= trial.lam):
            tot += chi(n) * trial.psi_out(y)
        else:                                    # ball straddles the band edge: both formulas agree there (phi = 0)
            tot += chi(n) * trial.psi_out(y) + arb(0, up(abs(trial.phi_ent(y))) / 2)
    Nt = max(nin, 0) + N_out
    Om0 = c * (Nt + 1) * w / rx
    assert Om0 > c
    h = c * w / rx
    t = trial.pref / 2 * C_E(trial.env_terms(), trial.beta, Om0) * (h ** -3) / (2 * arb(Nt) ** 2)
    return tot + arb(0, up(t))


def kb_symmetry(trial, D, x, s, w0):
    sigma = 1 if s == 0 else -1
    w0 = arb(w0)
    e1 = w0.sqrt() * kb_S(trial, D, x, w0)
    e2 = (1 / w0).sqrt() * kb_S(trial, D, x, 1 / w0)
    return e2 - sigma * e1, e1
