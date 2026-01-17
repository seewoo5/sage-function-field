from sage.all import *
from sage.repl.attach import load_attach_path

from tqdm import tqdm
from functools import lru_cache
import os


HERE = os.path.dirname(os.path.abspath(__file__))
SRC_DIR = os.path.normpath(os.path.join(HERE, "..", "ff", "sage"))
load_attach_path(SRC_DIR)
load("prime_count.sage")


def prime_count_cong_vec_naive(n, m):
    """
    Count the number of irreducible monic polynomials of degree n over GF(q).
    Check all monic polynomials of degree n.
    """
    R = m.parent()
    Mp = euler_totient(m)
    cnt_vec = vector([0] * Mp)
    g = unit_group_generator(m)

    if n < m.degree():
        for f in R.monics(n):
            for k in range(Mp):
                c = g^k % m
                if f % m == c and f.is_irreducible():
                    cnt_vec[k] += 1
    else:
        # vary h in f = hm + c over degree N - M monic polynomials
        # and check irreducibility, where M = deg m
        for h in R.monics(n - m.degree()):
            for k in range(Mp):
                f = h * m + (g^k % m)
                if f.is_irreducible():
                    cnt_vec[k] += 1
    cnt_vec = vector(ZZ, cnt_vec)
    return cnt_vec


def test_count(m, N):
    """
    For each congruence class, check if \pi(n, m, c) matches for 1 \le n \le N.
    """
    M = m.degree()
    Mp = euler_totient(m)
    t_ = m.parent().gen()
    for n in tqdm(range(1, N+1)):
        cnt_vec = prime_count_cong_vec(n, m)
        cnt_vec_naive = prime_count_cong_vec_naive(n, m)
        assert cnt_vec == cnt_vec_naive, f"Failed for deg = {n}, m = {m}: {cnt_vec} != {cnt_vec_naive}"
    print(f"Test passed for m = {m} over F_{m.parent().base_ring().cardinality()}, N={N}")


if __name__ == "__main__":
    # Manual test cases
    R2.<t> = GF(2)['t']
    test_ms = [t^2 + t + 1, t^3 + t + 1, t^4 + t + 1]
    N = 20
    for m in test_ms:
        print(f"[Prime Count] Testing for p = 2, modulus {m}")
        test_count(m, N)

    R3.<t> = GF(3)['t']
    test_ms = [t^2 + 1, t^3 - t + 1, t^2]
    N = 20
    for m in test_ms:
        print(f"[Prime Count] Testing for p = 3, modulus {m}")
        test_count(m, N)

    # Random test cases
    N = 10
    for p in [5, 7, 11]:
        while True:
            m = GF(p)['t'].random_element(2)
            if m.is_irreducible():
                break
        m = normalize(m)
        print(f"[Prime Count] Testing for p = {p}, modulus {m}")
        test_count(m, N)