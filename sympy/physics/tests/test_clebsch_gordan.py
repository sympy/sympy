from sympy import S, sqrt, pi
from sympy.physics.wigner import clebsch_gordan, wigner_9j, wigner_6j, gaunt
from sympy.core.numbers import Rational

# Todo: more tests should be added from:
# http://en.wikipedia.org/wiki/Table_of_Clebsch-Gordan_coefficients


def test_clebsch_gordan_docs():
    assert clebsch_gordan(S(3)/2, S(1)/2, 2, S(3)/2, S(1)/2, 2) == 1
    assert clebsch_gordan(S(3)/2, S(1)/2, 1, S(3)/2, -S(1)/2, 1) == sqrt(3)/2
    assert clebsch_gordan(S(3)/2, S(1)/2, 1, -S(1)/2, S(1)/2, 0) == -sqrt(2)/2


def test_clebsch_gordan1():
    j_1 = S(1)/2
    j_2 = S(1)/2
    m = 1
    j = 1
    m_1 = S(1)/2
    m_2 = S(1)/2
    assert clebsch_gordan(j_1, j_2, j, m_1, m_2, m) == 1

    j_1 = S(1)/2
    j_2 = S(1)/2
    m = 0
    j = 1
    m_1 = S(1)/2
    m_2 = S(1)/2
    assert clebsch_gordan(j_1, j_2, j, m_1, m_2, m) == 0

    j_1 = S(1)/2
    j_2 = S(1)/2
    m = 0
    j = 1
    m_1 = S(1)/2
    m_2 = -S(1)/2
    assert clebsch_gordan(j_1, j_2, j, m_1, m_2, m) == sqrt(2)/2

    j_1 = S(1)/2
    j_2 = S(1)/2
    m = 0
    j = 0
    m_1 = S(1)/2
    m_2 = -S(1)/2
    assert clebsch_gordan(j_1, j_2, j, m_1, m_2, m) == sqrt(2)/2

    j_1 = S(1)/2
    j_2 = S(1)/2
    m = 0
    j = 1
    m_1 = -S(1)/2
    m_2 = S(1)/2
    assert clebsch_gordan(j_1, j_2, j, m_1, m_2, m) == sqrt(2)/2

    j_1 = S(1)/2
    j_2 = S(1)/2
    m = 0
    j = 0
    m_1 = -S(1)/2
    m_2 = S(1)/2
    assert clebsch_gordan(j_1, j_2, j, m_1, m_2, m) == -sqrt(2)/2


def test_clebsch_gordan2():
    j_1 = S(1)
    j_2 = S(1)/2
    m = S(3)/2
    j = S(3)/2
    m_1 = 1
    m_2 = S(1)/2
    assert clebsch_gordan(j_1, j_2, j, m_1, m_2, m) == 1

    j_1 = S(1)
    j_2 = S(1)/2
    m = S(1)/2
    j = S(3)/2
    m_1 = 1
    m_2 = -S(1)/2
    assert clebsch_gordan(j_1, j_2, j, m_1, m_2, m) == 1/sqrt(3)

    j_1 = S(1)
    j_2 = S(1)/2
    m = S(1)/2
    j = S(1)/2
    m_1 = 1
    m_2 = -S(1)/2
    assert clebsch_gordan(j_1, j_2, j, m_1, m_2, m) == sqrt(2)/sqrt(3)

    j_1 = S(1)
    j_2 = S(1)/2
    m = S(1)/2
    j = S(1)/2
    m_1 = 0
    m_2 = S(1)/2
    assert clebsch_gordan(j_1, j_2, j, m_1, m_2, m) == -1/sqrt(3)

    j_1 = S(1)
    j_2 = S(1)/2
    m = S(1)/2
    j = S(3)/2
    m_1 = 0
    m_2 = S(1)/2
    assert clebsch_gordan(j_1, j_2, j, m_1, m_2, m) == sqrt(2)/sqrt(3)


def test_wigner():
    def tn(a, b):
        return abs((a - b).n(64) < S('1e-64'))
    assert tn(wigner_9j(1, 1, 1, 1, 1, 1, 1, 1, 0, prec=64), S(1)/18)
    assert wigner_9j(3, 3, 2, 3, 3, 2, 3, 3, 2) == 3221*sqrt(
        70)/(246960*sqrt(105)) - 365/(3528*sqrt(70)*sqrt(105))
    assert wigner_6j(5, 5, 5, 5, 5, 5) == Rational(1, 52)
    assert tn(wigner_6j(8, 8, 8, 8, 8, 8, prec=64), -S(12219)/965770)


def test_gaunt():
    def tn(a, b):
        return abs((a - b).n(64) < S('1e-64'))
    assert gaunt(1, 0, 1, 1, 0, -1) == -1/(2*sqrt(pi))
    assert tn(gaunt(
        10, 10, 12, 9, 3, -12, prec=64), (-S(98)/62031) * sqrt(6279)/sqrt(pi))
