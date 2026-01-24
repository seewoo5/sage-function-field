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

    # Pre-compute all residue classes and build a lookup dictionary
    residues = []
    c = R(1)
    for k in range(Mp):
        residues.append(c)
        c = (c * g) % m
    residue_to_idx = {c: k for k, c in enumerate(residues)}

    if n < m.degree():
        for f in R.monics(n):
            if f.is_irreducible():
                res = f % m
                if res in residue_to_idx:
                    cnt_vec[residue_to_idx[res]] += 1
    else:
        # vary h in f = hm + c over degree N - M monic polynomials
        # and check irreducibility, where M = deg m
        # Iterate over h once and check all congruence classes
        for h in R.monics(n - m.degree()):
            hm = h * m
            for k in range(Mp):
                f = hm + residues[k]
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
    N = 12
    for m in test_ms:
        print(f"[Prime Count] Testing for p = 3, modulus {m}")
        test_count(m, N)

    # Random test cases for larger characteristics with degree 2 irreducible modulus
    N = 5
    for p in [5, 7, 11]:
        while True:
            m = GF(p)['t'].random_element(2)
            if m.is_irreducible():
                break
        m = normalize(m)
        print(f"[Prime Count] Testing for p = {p}, modulus {m}")
        test_count(m, N)