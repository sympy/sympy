from __future__ import annotations

from sympy.ntheory import qs, qs_factor
from sympy.ntheory.qs import SievePolynomial, _generate_factor_base, \
    _initialize_first_polynomial, _initialize_ith_poly, \
    _gen_sieve_array, _check_smoothness, _trial_division_stage, _gauss_mod_2, \
    _find_factor
from sympy.testing.pytest import slow


@slow
def test_qs_1():
    assert qs(10009202107, 100, 10000) == {100043, 100049}
    assert qs(211107295182713951054568361, 1000, 10000) == \
        {13791315212531, 15307263442931}
    assert qs(980835832582657*990377764891511, 3000, 50000, seed=123) == \
        {980835832582657, 990377764891511}
    assert qs(18640889198609*20991129234731, 1000, 50000) == \
        {18640889198609, 20991129234731}


def test_qs_2() -> None:
    n = 10009202107
    M = 50
    sieve_poly = SievePolynomial(10, 80, n)
    assert sieve_poly.eval_v(10) == sieve_poly.eval_u(10)**2 - n == -10009169707
    assert sieve_poly.eval_v(5) == sieve_poly.eval_u(5)**2 - n == -10009185207

    idx_1000, idx_5000, factor_base = _generate_factor_base(2000, n)
    assert idx_1000 == 82
    assert [factor_base[i].prime for i in range(15)] == \
        [2, 3, 7, 11, 17, 19, 29, 31, 43, 59, 61, 67, 71, 73, 79]
    assert [factor_base[i].tmem_p for i in range(15)] == \
        [1, 1, 3, 5, 3, 6, 6, 14, 1, 16, 24, 22, 18, 22, 15]
    assert [factor_base[i].log_p for i in range(5)] == \
        [710, 1125, 1993, 2455, 2901]

    g, B = _initialize_first_polynomial(
        n, M, factor_base, idx_1000, idx_5000, seed=0)
    assert g.a == 1133107
    assert g.b == 682543
    assert B == [272889, 409654]
    assert [factor_base[i].soln1 for i in range(15)] == \
        [0, 0, 3, 7, 13, 0, 8, 19, 9, 43, 27, 25, 63, 29, 19]
    assert [factor_base[i].soln2 for i in range(15)] == \
        [0, 1, 1, 3, 12, 16, 15, 6, 15, 1, 56, 55, 61, 58, 16]
    assert [factor_base[i].b_ainv for i in range(5)] == \
        [[0, 0], [0, 2], [3, 0], [3, 9], [13, 13]]

    g_1 = _initialize_ith_poly(n, factor_base, 1, g, B)
    assert g_1.a == 1133107
    assert g_1.b == 136765

    sieve_array = _gen_sieve_array(M, factor_base)
    assert sieve_array[0:5] == [8424, 13603, 1835, 5335, 710]

    assert _check_smoothness(9645, factor_base) == (5, False)
    assert _check_smoothness(210313, factor_base)[0][0:15] == \
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 1]
    assert _check_smoothness(210313, factor_base)[1]

    partial_relations: dict[int, tuple[int, int]] = {}
    smooth_relation, partial_relation = _trial_division_stage(
        n, M, factor_base, sieve_array, sieve_poly, partial_relations,
        ERROR_TERM=25*2**10)

    assert partial_relations == {
        8699: (440, -10009008507),
        166741: (490, -10008962007),
        131449: (530, -10008921207),
        6653: (550, -10008899607)
    }
    assert [smooth_relation[i][0] for i in range(5)] == [
        -250, -670615476700, -45211565844500, -231723037747200, -1811665537200]
    assert [smooth_relation[i][1] for i in range(5)] == [
        -10009139607, 1133094251961, 5302606761, 53804049849, 1950723889]
    assert smooth_relation[0][2][0:15] == [
        1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]

    assert _gauss_mod_2(
        [[0, 0, 1], [1, 0, 1], [0, 1, 0], [0, 1, 1], [0, 1, 1]]
    ) == (
        [[[0, 1, 1], 3], [[0, 1, 1], 4]],
        [True, True, True, False, False],
        [[0, 0, 1], [1, 0, 0], [0, 1, 0], [0, 1, 1], [0, 1, 1]]
    )


def test_qs_3():
    N = 1817
    smooth_relations = [
        (2455024, 637, [0, 0, 0, 1]),
        (-27993000, 81536, [0, 1, 0, 1]),
        (11461840, 12544, [0, 0, 0, 0]),
        (149, 20384, [0, 1, 0, 1]),
        (-31138074, 19208, [0, 1, 0, 0])
    ]
    matrix = [s_relation[2] for s_relation in smooth_relations]
    dependent_row, mark, gauss_matrix = _gauss_mod_2(matrix)
    assert dependent_row == [[[0, 0, 0, 0], 2], [[0, 1, 0, 0], 3]]
    assert mark == [True, True, False, False, True]
    assert gauss_matrix == [
        [0, 0, 0, 1],
        [0, 1, 0, 0],
        [0, 0, 0, 0],
        [0, 1, 0, 0],
        [0, 1, 0, 1]
    ]

    factor = _find_factor(
        dependent_row, mark, gauss_matrix, 0, smooth_relations, N)
    assert factor == 23


def test_qs_4():
    N = 10007**2 * 10009 * 10037**3 * 10039
    for factor in qs(N, 1000, 2000):
        assert N % factor == 0
        N //= factor


def test_qs_factor():
    assert qs_factor(1009 * 100003, 2000, 10000) == {1009: 1, 100003: 1}
    assert qs_factor(1009**2 * 2003**2*30011*400009, 2000, 10000) == \
        {1633856814842812561: 1, 30011: 1}


def test_issue_27616():
    #https://github.com/sympy/sympy/issues/27616
    N = 9804659461513846513 + 1
    assert qs(N, 5000, 20000) is not None
