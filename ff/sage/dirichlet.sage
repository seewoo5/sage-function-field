"""Dirichlet character"""
from sage.all import *

from functools import lru_cache
from typing import Tuple, Optional


try:
    load("utils.sage")
except Exception:
    try:
        load(os.path.dirname(os.path.abspath(__file__)) + "/sage/utils.sage")
    except Exception:
        pass


def char_group_exponent(m):
    """Exponent of the Dirichlet character group modulo m."""
    if not m.is_monic():
        raise ValueError("The modulus must be monic.")
    if not is_square_free(m):
        raise NotImplementedError("The modulus must be square-free.")
    res = 1
    for fac, _ in m.factor():
        res = lcm(res, norm_poly(fac) - 1)
    return res


def is_generator(g, m):
    """
    Decide if g is a generator of the unit group modulo m.
    """
    if g.gcd(m).degree() > 0:
        return False
    Mp = euler_totient(m)
    for d in divisors(Mp):
        if d < Mp and (g ^ d) % m == 1:
            return False
    return True


def unit_group_generator(m):
    """
    Generator of the unit group modulo m, with smallest possible degree.
    Raise ValueError if no generator is found, i.e. the unit group is not cyclic.
    """
    R = m.parent()
    for g in R.monics(max_degree=m.degree() - 1):
        if is_generator(g, m):
            return g
    raise ValueError("No generator found for the unit group modulo m.")


class DirichletCharacterFF:
    """
    Dirichlet character for function fields modulo monic polynomials
    m = p_1^r_1 ... p_k^r_k, where (F_q[T] / p_i^r_i)^\times are all cyclic.
    """
    def __init__(self, m, exps: Optional[Tuple[int]] = None):
        if not m.is_monic():
            raise ValueError("The modulus must be monic.")
        self._modulus = m
        self._generators = self._init_generators()
        if exps is not None:
            self._exponents = self._init_exponents(exps)
        else:
            self._init_exponents_trivial()
        self._q = m.parent().base_ring().cardinality()
        self._p = m.parent().characteristic()

    def _init_generators(self) -> Tuple:
        return tuple(
            unit_group_generator(fac ^ e) for (fac, e) in self._modulus.factor()
        )

    def _init_exponents_trivial(self):
        self._exponents = [0] * len(list(self._modulus.factor()))

    def _init_exponents(self, exps):
        return tuple(
            exp % euler_totient(fac ^ e) for (fac, e), exp in \
            zip(self._modulus.factor(), exps)
        )

    def is_trivial(self):
        return all(
            exp % euler_totient(fac ^ e) == 0 for (fac, e), exp in \
            zip(self._modulus.factor(), self._exponents)
        )

    def m(self):
        return self._modulus

    def exps(self):
        return self._exponents

    def q(self):
        return self._q

    def p(self):
        return self._p

    def ring(self):
        return self._modulus.parent()

    def image_conductor(self):
        """For this N, the image of the character will lie in Q(\zeta_N)"""
        N = 1
        for (fac, e), exp in zip(self._modulus.factor(), self._exponents):
            N_ = euler_totient(fac ^ e)
            if exp != 0:
                N_ /= gcd(N_, exp)  # Optimal N
            else:
                N_ = 1
            if N_ == 2:  # Q(\zeta_2) = Q
                N_ = 1
            N = lcm(N, N_)
        return N

    def __call__(self, f):
        f = self.ring()(f) % self._modulus
        CF = CyclotomicField(self.image_conductor())
        res = CF(1)
        for (fac, e), g, exp in zip(
            self._modulus.factor(), self._generators, self._exponents
        ):
            if f % fac == 0:
                return CF(0)
            else:
                # If f \equiv g^k mod fac, then the output is
                # \zeta_N^{k * exp} where N = |fac| - 1 = q^deg(fac) - 1
                N_ = euler_totient(fac ^ e)
                zeta_N_ = CyclotomicField(N_).gen()
                for k in range(N_):
                    if (f - g^k) % (fac ^ e) == 0:
                        res *= zeta_N_ ** (k * exp)
                        break
        return res

    def lfunc(self):
        """
        L-function associated to the Dirichlet character.
        When character is nontrivial, it is a polynomial in u = q^(-s)
        of degree at most deg(m) - 1.
        When character is trivial, it is the rational function
        (\prod_{P|m} (1 - u^{\deg P})) / (1 - qu).
        """
        # TODO: more efficient implementation using functional equation
        if self.is_trivial():
            K.<u> = FunctionField(QQ)
            Lfunc = 1 / (1 - self.q() * u)
            for fac, _ in self._modulus.factor():
                Lfunc *= (1 - u ^ fac.degree())
        else:
            N = self.image_conductor()  # Use the optimal N
            R.<u> = PolynomialRing(CyclotomicField(N), 'u')
            Lfunc = R(1)
            for n in range(1, self._modulus.degree()):
                for f in self._modulus.parent().monics(n):
                    Lfunc += self(f) * (u ^ n)
        return Lfunc

    def is_even(self):
        """
        Check if character is even, i.e. trivial on $\mathbb{F}_q^\times$.
        """
        for a in self._modulus.parent().base_ring():
            if a != 0:
                if self(a) != 1:
                    return False
        return True

    def is_odd(self):
        return not self.is_even()

    def is_primitive(self):
        """
        Check if the character is primitive, i.e. not induced from a character
        of a smaller modulus.
        """
        raise NotImplementedError

    def gauss_sum(self):
        """
        Compute the Gauss sum associated to odd Dirichlet character, defined as
        \sum_{a \in \mathbb{F}_q^\times} \chi(a) \zeta_p^{\text{Tr}(a)}.
        For even character, it simply returns 1.
        """
        if self.is_even():
            return CyclotomicField(self._p)(1)
        zp = CyclotomicField(self._p).gen()
        res = 0
        for a in self._modulus.parent().base_ring():
            if a != 0:
                res += self(a) * zp ^ (a.trace())
        return res

    def gauss_sum_sign(self):
        """
        Sign of Gauss sum, defined as Gauss(chi) * q^(-1/2).
        """
        if self.is_even():
            return CyclotomicField(self._p)(1)
        return self.gauss_sum() * self.q() ** (-1 / 2)

    def __mul__(self, other):
        """
        Multiply two Dirichlet characters.
        """
        if not isinstance(other, DirichletCharacterFF):
            raise TypeError("Can only multiply with another DirichletCharacterFF.")
        if self._modulus != other._modulus:
            raise ValueError("Can only multiply characters with the same modulus.")
        exps = tuple(
            (x + y) % euler_totient(fac ^ e) for (fac, e), (x, y) in zip(self._modulus.factor(), zip(self._exponents, other._exponents))
        )
        return DirichletCharacterFF(self._modulus, exps)

    def __pow__(self, other: Integer):
        """
        Compute power of Dirichlet character. chi2 = chi1 ^ other.
        """
        exps = tuple((x * other) % euler_totient(fac ^ e) for (fac, e), x in zip(self._modulus.factor(), self._exponents))
        return DirichletCharacterFF(self._modulus, exps)

    def conjugate(self):
        """
        Compute the complex conjugate Dirichlet character.
        """
        return DirichletCharacterFF(self._modulus, tuple(-x for x in self._exponents))

    def __eq__(self, other):
        if not isinstance(other, DirichletCharacterFF):
            return TypeError("Can only compare with another DirichletCharacterFF.")
        return (self._modulus == other._modulus and
                self._exponents == other._exponents)

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
        if not m.is_squarefree():
            raise ValueError("Quadratic Dirichlet characters are only defined for square-free moduli.")
        exps = self._init_exponents_quad(m)
        super().__init__(m, exps)

    def _init_exponents_quad(self, m):
        return tuple((norm_poly(fac) - 1) // 2 for fac, _ in m.factor())

    def __call__(self, f):
        f = self.ring()(f) % self._modulus
        return quad_char(f, self._modulus)

    def conjugate(self):
        """
        Compute the complex conjugate Dirichlet character.
        For quadratic (real) characters, this is just the character itself.
        """
        return DirichletCharacterFFQuadratic(self._modulus)


def random_dirichlet_character(m):
    """
    Return random Dirichlet character modulo m.
    """
    q = m.parent().base_ring().cardinality()
    exponents = []
    for fac, e in m.factor():
        # exponent randomly sampled in [0, euler_totient(fac ^ e))
        rexp = randint(0, euler_totient(fac ^ e) - 1)
        exponents.append(rexp)
    return DirichletCharacterFF(m, tuple(exponents))
