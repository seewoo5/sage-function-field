"""Dirichlet character"""
from sage.all import *

from typing import Tuple, Optional


load(os.path.dirname(os.path.abspath(__file__)) + "/sage/utils.sage")


def order_char_group(m):
    """Order of the Dirichlet character group modulo m."""
    if not m.is_monic():
        raise ValueError("The modulus must be monic.")
    if not is_square_free(m):
        raise NotImplementedError("The modulus must be square-free.")
    res = 1
    for fac, _ in m.factor():
        res = lcm(res, norm_poly(fac) - 1)
    return res


class DirichletCharacterFF:
    """
    Dirichlet character for function fields modulo
    monic square-free polynomials
    """
    def __init__(self, m, exps: Optional[Tuple[int]] = None):
        if not m.is_monic():
            raise ValueError("The modulus must be monic.")
        if not is_square_free(m):
            raise NotImplementedError("The modulus must be square-free.")
        self._modulus = m
        if exps is not None:
            self._exponents = self._init_exponents(exps)
        else:
            self._init_exponents_trivial()
        self._q = m.parent().base_ring().cardinality()
        self._p = m.parent().characteristic()

    def _init_exponents_trivial(self):
        self._exponents = [0] * len(list(self._modulus.factor()))

    def _init_exponents(self, exps):
        return tuple(
            exp % (norm_poly(fac) - 1) for (fac, _), exp in zip(self._modulus.factor(), exps)
        )

    def is_trivial(self):
        return all(exp == 0 for exp in self._exponents)

    def m(self):
        return self._modulus

    def exps(self):
        return self._exponents

    def q(self):
        return self._q

    def p(self):
        return self._p

    def image_conductor(self):
        """For this N, the image of the character will lie in Q(\zeta_N)"""
        N = 1
        for (fac, _), exp in zip(self._modulus.factor(), self._exponents):
            N_ = norm_poly(fac) - 1
            if exp != 0:
                N_ /= gcd(N_, exp)  # Optimal N
            else:
                N_ = 1
            if N_ == 2:  # Q(\zeta_2) = Q
                N_ = 1
            N = lcm(N, N_)
        return N

    def __call__(self, f):
        CF = CyclotomicField(self.image_conductor())
        res = CF(1)
        for (fac, _), exp in zip(self._modulus.factor(), self._exponents):
            if f % fac == 0:
                return CF(0)
            else:
                # If f \equiv t^k mod fac, then the output is
                # \zeta_N^{k * exp} where N = |fac| - 1 = q^deg(fac) - 1
                N_ = norm_poly(fac) - 1
                zeta_N_ = CyclotomicField(N_).gen()
                t_ = fac.parent().gen()
                for k in range(N_):
                    if (f - t_^k) % fac == 0:
                        res *= zeta_N_ ** (k * exp)
                        break
        return res

    def lfunc(self):
        """
        L-function associated to the Dirichlet character.
        Polynomial in u = q^(-s) of degree at most deg(m) - 1.
        """
        N = self.image_conductor()  # Use the optimal N
        R.<u> = PolynomialRing(CyclotomicField(N), 'u')
        Lfunc = R(1)
        for n in range(1, self._modulus.degree()):
            for f in self._modulus.parent().polynomials(n):
                if f.is_monic():
                    Lfunc += self(f) * (u ^ n)
        return Lfunc

    def __mul__(self, other):
        """
        Multiply two Dirichlet characters.
        """
        if not isinstance(other, DirichletCharacterFF):
            raise TypeError("Can only multiply with another DirichletCharacterFF.")
        if self._modulus != other._modulus:
            raise ValueError("Can only multiply characters with the same modulus.")
        exps = tuple(
            (x + y) % (norm_poly(fac) - 1) for (fac, _), (x, y) in zip(self._modulus.factor(), zip(self._exponents, other._exponents))
        )
        return DirichletCharacterFF(self._modulus, exps)

    def __repr__(self):
        return f"DirichletCharacter of modulus {self._modulus} over GF({self._q}) with exponents={self._exponents}"


def quad_char(f, m):
    # Quadratic character modulo monic square-free polynomial m
    def _quad_char_irred(f, P):
        # Quadratic character of f modulo P
        assert P.is_irreducible(), "Polynomial P must be irreducible"
        assert P.is_monic(), "Polynomial P must be monic"
        # Characteristic has to be odd for quadratic character
        assert f.parent().base_ring().characteristic() != 2, "Characteristic must be odd"
        if f % P == 0:
            return 0
        f = f % P
        P_norm = norm_poly(P)
        if (f ^ ((P_norm - 1) // 2)) % P == 1:
            return 1
        elif (f ^ ((P_norm - 1) // 2)) % P == -1:
            return -1
        else:
            raise ValueError("Unexpected value for quadratic character")
    # Quadratic character of f modulo P
    # Not necessarily irreducible but square-free
    assert m.is_monic(), "Polynomial P must be monic"

    # Factor P into irreducible factors
    res = 1
    for factor, e in m.factor():
        if e > 1:
            raise ValueError("Polynomial P must be square-free")
        res *= _quad_char_irred(f, factor)
    return res


class DirichletCharacterFFQuadratic(DirichletCharacterFF):
    """
    Nontrivial Quadratic Dirichlet character associated to a given
    monic square-free modulus m. If m = P_1 P_2 ... P_r, then
    the character is defined as a product of quadratic characters
    \chi_{P_i}(f) = f ^ ((|P_i| - 1) / 2).
    Only over odd characteristic fields.
    """
    def __init__(self, m):
        if m.parent().characteristic() == 2:
            raise ValueError("Quadratic Dirichlet characters are not defined over characteristic 2 fields.")
        exps = self._init_exponents_quad(m)
        super().__init__(m, exps)

    def _init_exponents_quad(self, m):
        return tuple((norm_poly(fac) - 1) // 2 for fac, _ in m.factor())

    def __call__(self, f):
        return quad_char(f, self._modulus)
