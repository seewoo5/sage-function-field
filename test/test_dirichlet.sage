from sage.all import *
from sage.repl.attach import load_attach_path

from tqdm import tqdm
from functools import lru_cache
import os


HERE = os.path.dirname(os.path.abspath(__file__))
SRC_DIR = os.path.normpath(os.path.join(HERE, "..", "ff", "sage"))
load_attach_path(SRC_DIR)
load("dirichlet.sage")


def test_gauss_sum(chi):
    """Test if Gauss sum has norm q^(1/2)"""
    gs = chi.gauss_sum()
    if chi.is_even():
        assert gs * gs.conjugate() == 1, f"[Gauss sum] Failed for character {chi}, {gs * gs.conjugate()} != {1}"
    else:
        assert gs * gs.conjugate() == chi.q(), f"[Gauss sum] Failed for character {chi}, {gs * gs.conjugate()} != {chi.q()}"


def test_conjugate(chi):
    """Test if chi * chi_bar is trivial."""
    chi_bar = chi.conjugate()
    assert (chi * chi_bar).is_trivial(), f"[Conjugate] Failed for character {chi}"


def test_trivial(m):
    chi0 = DirichletCharacterFF(m)
    assert chi0.is_trivial(), f"[Trivial] Failed for modulus {m}"


def test_quadratic(m):
    chi_m = DirichletCharacterFFQuadratic(m)
    assert chi_m.conjugate() == chi_m, f"[Quadratic] chi_m is not self-conjugate for modulus {m}"
    assert chi_m.image_conductor() == 1, f"[Quadratic] Conductor of chi_m is not 1 for modulus {m}"

    q = chi_m.q()
    f = GF(q)['t'].random_element(10)
    assert chi_m(f) in [-1, 0, 1], f"[Quadratic] chi_m(f) is not in {-1, 0, 1} for modulus {m} and input {f} over GF({q})"

    test_lfunc(chi_m)
    test_gauss_sum(chi_m)


def test_lfunc(chi):
    # TODO: Test functional equation
    q = chi.q()
    RC.<u> = CC['u']
    Lchi = chi.lfunc()
    if chi.is_even():
        assert Lchi(1) == 0, f"[L-function] L(1, chi) != 0 for even character {chi}"

    err = 1e-8
    for root, _ in PolynomialRing(CC, 'u')(Lchi).roots():
        if abs(root.abs() - 1) > err:
            assert abs(root.abs() - q^(-1/2)) < err, f"[L-function] Root {root} norm {root.norm()} is not close to q^(-1/2) for character {chi}"


def test_suite(m):
    chi = random_dirichlet_character(m)
    q = chi.q()

    test_trivial(m)
    test_conjugate(chi)
    test_gauss_sum(chi)
    if not chi.is_trivial():
        test_lfunc(chi)
    if q % 2 == 1 and m.is_squarefree():
        test_quadratic(m)


if __name__ == "__main__":
    # Manual test cases (modulus)
    print("[Dirichlet] Manual Test Cases")
    R2.<t> = GF(2)['t']
    for m in [t^2, t^2 + t + 1, t^3 + t + 1]:
        test_suite(m)
    R3.<t> = GF(3)['t']
    for m in [t^2, t^2 + 1, t^3 - t + 1]:
        test_suite(m)

    # Random test cases
    print("[Dirichlet] Random Test Cases")
    for q in [2, 3, 4, 5, 7, 9]:
        while True:
            m = GF(q)['t'].random_element(3)
            if m.is_irreducible():
                break
        m = normalize(m)
        test_suite(m)
