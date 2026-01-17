"""Prime counting functions"""
from sage.all import *
from functools import lru_cache
import os
import polars as pl


try:
    load("dirichlet.sage")
except Exception:
    try:
        load(os.path.dirname(os.path.abspath(__file__)) + "/sage/dirichlet.sage")
    except Exception:
        pass


def prime_count(n, q):
    """
    Count the number of irreducible monic polynomials of degree n over GF(q).
    Use Gauss formula: pi(n) = (1 / n) * \sum_{d|n} μ(d) * q^{n/d}
    """
    cnt = 0
    for d in divisors(n):
        cnt += moebius(d) * q^(n/d)
    cnt /= n
    return cnt


def prime_count_cong_naive(n, m, c):
    """
    Count the number of irreducible monic polynomials of degree n over GF(q)
    which is c modulo m.
    Check all monic polynomials of degree n.
    """
    R = m.parent()
    c = R(c)
    cnt = 0
    if n < m.degree():
        for f in R.monics(n):
            if f % m == c and f.is_irreducible():
                cnt += 1
    else:
        # vary g in f = gm + c over degree N - M monic polynomials
        # and check irreducibility, where M = deg m
        for g in R.monics(n - m.degree()):
            f = g * m + c
            if f.is_irreducible():
                cnt += 1
    return cnt


def sum_pow_zeros(chi, n):
    """
    Let L(u, chi) = \prod_j (1 - \alpha_j u) be Dirichlet series.
    Compute \sum_j \alpha_j^n, from the generating function 
    -u * d/du log L(u, chi) = -u * L'(u) / L(u).
    """
    L_chi = chi.lfunc()
    L_chi_deriv = L_chi.derivative()
    u = L_chi.parent().gen()
    R.<u> = PowerSeriesRing(L_chi.parent().base_ring(), default_prec=n+1)
    gen_func = - u * L_chi_deriv / R(L_chi)
    return gen_func[n]


def s(m, n):
    """
    Sum of deg P for irreducible P | m with (deg P) | n.
    """
    res = 0
    for factor, _ in m.factor():
        degP = factor.degree()
        if n % degP == 0:
            res += degP
    return res


def U(m, n):
    # for mth root of unity
    # m by m matrix
    # [1 1 ... 1 / 1 z^n z^2n .. z^(m-1)n / 1 z^2n z^4n ... z^(2(m-1))n / ... / 1 z^(m-1)n ... z^((m-1)(m-1))n ] for z = \zeta_m
    z = CyclotomicField(m).gen()
    n = n % m
    mat = []
    for i in range(m):
        row = []
        for j in range(m):
            row.append(z^(i * j * n))
        mat.append(row)
    return Matrix(mat)


@lru_cache(maxsize=None)
def U_inv(m, n):
    """Mobius inverse of U(m, n)"""
    if n == 1:
        return (1 / m) * U(m, -1)
    else:
        s = 0
        for d in divisors(n):
            if d != n:
                s -= U_inv(m, d) * U(m, n / d)
            else:
                break
        return s * U(m, -1) * (1 / m)


def prime_count_cong_vec(n, m):
    """
    Count the number of irreducible monic polynomials of degree n over GF(q) for
    each congruence class.
    Assume (A / m)^\times is cyclic (e.g. m is irreducible).
    The output vector give counts in order of 1, g, g^2, ..., g^{M'-1} where
    g is a generator of the unit group modulo m.
    """
    q = m.parent().base_ring().cardinality()
    M = m.degree()
    Mp = euler_totient(m)
    cnt_vec = vector([0] * Mp)
    chi = DirichletCharacterFF(m, (1,))
    for d in divisors(n):
        vec = [q^(n/d) - s(m, n/d)]
        for l in range(1, Mp):
            vec.append(-sum_pow_zeros(chi^l, n / d))
        cnt_vec += U_inv(Mp, d) * vector(vec)
    cnt_vec = cnt_vec / n
    cnt_vec = vector(ZZ, cnt_vec)
    return cnt_vec


def prime_count_cong(n, m, c):
    if m.gcd(c).degree() > 0:
        return 0
    cnt_vec = prime_count_cong_vec(n, m)
    Mp = len(cnt_vec)
    t_ = m.parent().gen()
    for k in range(Mp):
        if (c - t_^k) % m == 0:
            return ZZ(cnt_vec[k])
    raise ValueError(f"No solution for {c} modulo {m} in prime_count_cong.")


def prime_count_table(N, m):
    """
    Create a table of prime counts for each congruence class mod m
    up to degree N.
    """
    rows = []
    g = unit_group_generator(m)
    for n in range(1, N+1):
        cnt_vec = prime_count_cong_vec(n, m)
        rows.append({"deg" : n} | {f"({g})^{k} = {g^k % m}" : cnt_vec[k] for k in range(len(cnt_vec))})
    df = pl.from_dicts(rows)
    return df
