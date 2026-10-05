"""Rédei fields, Rédei symbols and the Rédei L-function L(s, ρ), via PARI (cypari2 2.2.4, PARI 2.17.2).

For primes p1, p2 ≡ 1 mod 4 with (p1/p2) = 1, R is the unique D4 octic field that contains Q(√p1, √p2) and is
ramified only at p1 and p2. ρ is its 2-dimensional irreducible representation, so L(s, ρ) is a degree-2 L-function
of conductor p1·p2 with gamma factor Γ_R(s)² or Γ_R(s+1)² (ρ is even: det ρ = χ_{p1 p2}). Facts used (issue #7 review):
  a_p = tr ρ(Frob_p) = 0 if p is linked to p1 or p2 (a Legendre symbol −1),
  a_p = 2·[p1, p2, p] if p is pairwise unlinked, where [p1, p2, p] = +1 iff p splits completely in R.
The repository's older Rédei code (src/cascade_redei.rs, src/redei_pari.rs, data/borromean_triples.db) is wrong or
incomplete and is not used here.
"""
import math
import cypari2

pari = cypari2.Pari()
pari.allocatemem(4 * 10**9)


def legendre_solutions(p1, p2, bound=200):
    """Primitive (x, y, z) with x² = p1 y² + p2 z², z ≥ 1."""
    out = []
    for y in range(bound):
        for z in range(1, bound):
            x2 = p1 * y * y + p2 * z * z
            x = math.isqrt(x2)
            if x * x == x2 and math.gcd(math.gcd(x, y), z) == 1:
                out.append((x, y, z))
    return out


def redei_field(p1, p2):
    """polredabs polynomial of the D4 octic containing Q(√p1, √p2), unramified outside p1·p2 (|disc| = (p1 p2)^4)."""
    for (x, y, z) in legendre_solutions(p1, p2):
        for d in (1, -1, 2, -2):
            f = pari(f"x^4 - 2*({d})*({x})*x^2 + ({d})^2*(({x})^2 - ({p1})*({y})^2)")
            if not pari.polisirreducible(f):
                continue
            S = pari.polredabs(pari.nfsplitting(f))
            if int(pari.poldegree(S)) != 8:
                continue
            if abs(int(pari.nfdisc(S))) == (p1 * p2) ** 4:
                return S, (x, y, z, d)
    raise ValueError(f"no Rédei field found for ({p1}, {p2})")


_nf_cache = {}


def _nf(p1, p2):
    key = (p1, p2)
    if key not in _nf_cache:
        S, _ = redei_field(p1, p2)
        _nf_cache[key] = (S, pari.nfinit(S))
    return _nf_cache[key]


def redei_symbol(p1, p2, r):
    """[p1, p2, r] from splitting in R: +1 split completely, −1 split in Q(√p1, √p2) but not in R, 0 if linked."""
    if pari.kronecker(r, p1) != 1 or pari.kronecker(r, p2) != 1:
        return 0
    _, nf = _nf(p1, p2)
    n = len(pari.idealprimedec(nf, r))
    if n == 8:
        return 1
    if n == 4:
        return -1
    raise ValueError(f"unexpected splitting of {r} in R({p1},{p2}): {n} primes")


def admissible(p1, p2):
    return p1 % 4 == 1 and p2 % 4 == 1 and p1 != p2 and pari.kronecker(p1, p2) == 1


def set_precision(digits):
    pari(f"default(realprecision, {digits})")


def rho_lfun(p1, p2, name="Lr"):
    """Create the PARI L-function of ρ under the gp name `name` (at the current precision); return field data."""
    S, nf = _nf(p1, p2)
    pari(f"nfR_{name} = nfinit({S}); G_{name} = galoisinit(nfR_{name}); T_{name} = galoischartable(G_{name});")
    ncol = int(pari(f"#T_{name}[1]"))
    cols = [j for j in range(1, ncol + 1) if int(pari(f"vecmax(abs(T_{name}[1][,{j}]))")) == 2]
    if len(cols) != 1:
        raise ValueError(f"expected one 2-dimensional character, found columns {cols}")
    pari(f"{name} = lfunartin(nfR_{name}, G_{name}, T_{name}[1][,{cols[0]}], T_{name}[2]);")
    return S


def lfun_dirichlet(d, name):
    """L(s, χ_d) for a fundamental discriminant d (PARI lfuncreate(d))."""
    pari(f"{name} = lfuncreate({d});")


def lfun_product(a, b, name):
    pari(f"{name} = lfunmul({a}, {b});")


def params(name):
    """(conductor, Vga list, root number) of a PARI L-function (lfunparams = [N, k, Vga])."""
    N = int(pari(f"lfunparams({name})[1]"))
    vga = [int(v) for v in pari(f"lfunparams({name})[3]")]
    w = pari(f"lfunrootres({name})[3]")
    return N, vga, complex(w)


def an_list(name, nmax):
    return [0] + [int(v) for v in pari(f"lfunan({name}, {nmax})")]


def zeros(name, a, b, divz=8):
    """Zeros γ in [a, b] on the critical line (PARI lfunzeros), as floats."""
    return [float(v) for v in pari(f"lfunzeros({name}, [{a}, {b}], {divz})")]
