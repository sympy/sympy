from sympy import S, sqrt
from sympy.physics.wigner import clebsch_gordan

# Todo: more tests should be added from:
# http://en.wikipedia.org/wiki/Table_of_Clebsch-Gordan_coefficients

def test_clebsch_gordan_docs():
    assert clebsch_gordan(S(3)/2,S(1)/2,2, S(3)/2,S(1)/2,2) == 1
    assert clebsch_gordan(S(3)/2,S(1)/2,1, S(3)/2,-S(1)/2,1) == sqrt(3)/2
    assert clebsch_gordan(S(3)/2,S(1)/2,1, -S(1)/2,S(1)/2,0) == -sqrt(2)/2

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
