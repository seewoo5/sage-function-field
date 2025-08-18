"""Basic functions"""
from sage.all import *


def leading_coeff(f):
    # Leading coefficient of a polynomial f
    return f.coefficient(f.degree())


def normalize(f):
    # Divide by leading coefficient
    return f / leading_coeff(f)


def norm_poly(f):
    # Norm of a polynomial f, defined as q^{deg(f)}
    # where q is the cardinality of the base field
    q = f.parent().base_ring().cardinality()
    return q ** f.degree()


def is_square_free(f):
    # Check if polynomial f is square-free
    for _, e in f.factor():
        if e > 1:
            return False
    return True


def euler_totient(f):
    # Euler's totient function for polynomial f
    # defined as |f| * prod_{P | f} (1 - |P|^{-1}) = \prod_{P | f} (|P| - 1) * |P|^{e-1}
    # where |f| is the norm of f
    assert f.is_monic(), "Polynomial f must be monic"

    res = 1
    for factor, e in f.factor():
        nm = norm_poly(factor)
        res *= (nm - 1) * nm ** (e - 1)
    return res


def moebius(f):
    # Mobius function of f
    if not is_square_free(f):
        return 0
    return (-1) ** omega(f)


def omega(f):
    # Number of prime factors of f, counted with multiplicity
    return sum([e for _, e in f.factor()])


def lamda(f):
    # Lambda function of f
    return (-1) ** omega(f)
