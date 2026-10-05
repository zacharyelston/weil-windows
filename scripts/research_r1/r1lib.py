"""R1 helpers: the contribution of a zero quadruple ½ ± δ ± iγ₀ to Weil's form on rh2's window
[−L/2, L/2] (L = 2a = log x), in rh2's orthonormal bases, and the Sherman–Morrison tests.

Zeros side: Q(f) = Σ_ρ F(ρ)F(1−ρ), F(s) = ∫ f(u) e^{(s−½)u} du. For real f and the quadruple
{ρ, ρ̄, 1−ρ, 1−ρ̄}, ρ = ½ + δ + iγ₀ (δ ≠ 0, γ₀ ≠ 0):

    T(f) = 4 Re[F̂(γ₀ − iδ) F̂(−γ₀ + iδ)] = 4 (|f̂_c(γ₀)|² − |f̂_s(γ₀)|²),
    f_c = f cosh δu,  f_s = f sinh δu,  F̂(t) = ∫ f e^{itu} du.

Multiplicity m: 4 for a quadruple, 2 for a pair (δ = 0 or γ₀ = 0), 1 for the point s = ½.
In a real orthonormal basis b_k the term is the rank-2 indefinite matrix m (c cᵀ − s sᵀ) with
  even sector (b_0 = 1/√L, b_k = √(2/L) cos ω_k u):   c_k = ∫ b_k cosh δu cos γ₀u,  s_k = ∫ b_k sinh δu sin γ₀u
  odd sector  (b_k = √(2/L) sin ω_k u, k ≥ 1):         c_k = ∫ b_k cosh δu sin γ₀u,  s_k = ∫ b_k sinh δu cos γ₀u
"""
import mpmath as mp


def I_cc(nu, d, L):
    """∫_{−L/2}^{L/2} cos(νu) cosh(δu) du."""
    h = L / 2
    if nu == 0 and d == 0:
        return L
    return 2 * (d * mp.sinh(d * h) * mp.cos(nu * h) + nu * mp.cosh(d * h) * mp.sin(nu * h)) / (d * d + nu * nu)


def I_ss(nu, d, L):
    """∫_{−L/2}^{L/2} sin(νu) sinh(δu) du."""
    h = L / 2
    if nu == 0 or d == 0:
        return mp.mpf(0)
    return 2 * (d * mp.cosh(d * h) * mp.sin(nu * h) - nu * mp.sinh(d * h) * mp.cos(nu * h)) / (d * d + nu * nu)


def pair_vectors(delta, gamma, L, n, parity):
    """(c, s) as mp column vectors in rh2's basis of the given parity."""
    d, g = mp.mpf(delta), mp.mpf(gamma)
    r = mp.sqrt(2 / L)
    om = lambda k: 2 * mp.pi * k / L
    if parity == "even":
        c = [I_cc(g, d, L) / mp.sqrt(L)] + [r * (I_cc(om(k) + g, d, L) + I_cc(om(k) - g, d, L)) / 2 for k in range(1, n + 1)]
        s = [I_ss(g, d, L) / mp.sqrt(L)] + [r * (I_ss(g + om(k), d, L) + I_ss(g - om(k), d, L)) / 2 for k in range(1, n + 1)]
    else:
        c = [r * (I_cc(om(k) - g, d, L) - I_cc(om(k) + g, d, L)) / 2 for k in range(1, n + 1)]
        s = [r * (I_ss(om(k) + g, d, L) + I_ss(om(k) - g, d, L)) / 2 for k in range(1, n + 1)]
    return mp.matrix(c), mp.matrix(s)


def multiplicity(delta, gamma):
    return (2 if delta != 0 else 1) * (2 if gamma != 0 else 1)


def pair_matrix(delta, gamma, L, n, parity):
    c, s = pair_vectors(delta, gamma, L, n, parity)
    m = multiplicity(delta, gamma)
    return m * (c * c.T - s * s.T), c, s, m


# ---- direct check: F(s) by quadrature ------------------------------------------------------------

def basis_fn(coeffs, L, parity):
    r = mp.sqrt(2 / L)
    if parity == "even":
        return lambda u: coeffs[0] / mp.sqrt(L) + mp.fsum(coeffs[k] * r * mp.cos(2 * mp.pi * k * u / L) for k in range(1, len(coeffs)))
    return lambda u: mp.fsum(coeffs[k - 1] * r * mp.sin(2 * mp.pi * k * u / L) for k in range(1, len(coeffs) + 1))


def F_of_s(f, s, L):
    return mp.quad(lambda u: f(u) * mp.exp((s - mp.mpf(1) / 2) * u), [-L / 2, 0, L / 2])


def T_direct(f, delta, gamma, L):
    """Σ over the distinct points of the orbit of ρ under s → 1−s, s → s̄ of F(ρ)F(1−ρ)."""
    rho = mp.mpc(mp.mpf(1) / 2 + delta, gamma)
    pts = {rho, mp.conj(rho), 1 - rho, 1 - mp.conj(rho)}
    # mpc is unhashable-ish across equal values; dedupe numerically
    uniq = []
    for p in pts:
        if not any(abs(p - q) < mp.mpf(10) ** (-mp.mp.dps + 5) for q in uniq):
            uniq.append(p)
    return mp.re(mp.fsum(F_of_s(f, p, L) * F_of_s(f, 1 - p, L) for p in uniq))


# ---- Sherman–Morrison tests ----------------------------------------------------------------------

class Spectral:
    """Eigendecomposition of a symmetric positive matrix, for quadratic forms in its inverse."""

    def __init__(self, Q):
        vals, vecs = mp.eigsy(Q)
        self.vals = vals
        self.vecs = vecs
        self.n = Q.rows
        self.lam_min = min(vals)

    def coords(self, v):
        return self.vecs.T * v

    def inv_form(self, u, v=None):
        """uᵀ Q⁻¹ v (v defaults to u)."""
        a = self.coords(u)
        b = a if v is None else self.coords(v)
        return mp.fsum(a[i] * b[i] / self.vals[i] for i in range(self.n))


def reverse_test(spec, c, s, m):
    """Consistency of Q with a hypothetical term m(ccᵀ − ssᵀ): returns R := m cᵀ(Q + m ssᵀ)⁻¹ c.
    Q − m(ccᵀ − ssᵀ) ⪰ 0  ⇔  R ≤ 1 (for Q ≻ 0). R > 1 means the quadruple is inconsistent with Q."""
    cQc = spec.inv_form(c)
    sQs = spec.inv_form(s)
    cQs = spec.inv_form(c, s)
    return m * (cQc - m * cQs * cQs / (1 + m * sQs)), cQc, sQs


def forward_test(specP, c, s, m):
    """Positivity of P + m(ccᵀ − ssᵀ) for P ≻ 0: returns σ := m sᵀ(P + m ccᵀ)⁻¹ s; the form fails iff σ > 1."""
    sPs = specP.inv_form(s)
    cPc = specP.inv_form(c)
    cPs = specP.inv_form(c, s)
    return m * (sPs - m * cPs * cPs / (1 + m * cPc))


def naive_min(c, s, m):
    """Most negative eigenvalue of the pair term m(ccᵀ − ssᵀ) alone (its minimum Rayleigh quotient)."""
    cc = (c.T * c)[0]
    ss = (s.T * s)[0]
    cs = (c.T * s)[0]
    tr = cc - ss
    disc = mp.sqrt(tr * tr + 4 * (cc * ss - cs * cs))
    return m * (tr - disc) / 2


def lam_min(M):
    return min(mp.eigsy(M, eigvals_only=True))
