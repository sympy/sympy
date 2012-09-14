"""Tests for sparse distributed polynomials. """

from sympy.polys.distributedpolys import (
    sdp_LC, sdp_LM, sdp_LT, sdp_del_LT,
    sdp_coeffs, sdp_monoms,
    sdp_sort, sdp_strip, sdp_normal,
    sdp_from_dict, sdp_to_dict,
    sdp_indep_p, sdp_one_p, sdp_one, sdp_term_p,
    sdp_abs, sdp_neg,
    sdp_add_term, sdp_sub_term, sdp_mul_term,
    sdp_add, sdp_sub, sdp_mul, sdp_sqr, sdp_pow,
    sdp_monic, sdp_content, sdp_primitive,
    _term_rr_div, _term_ff_div,
    sdp_div, sdp_quo, sdp_rem,
    sdp_lcm, sdp_gcd,
)

from sympy.polys.monomialtools import (
    lex, grlex, grevlex,
)

from sympy.polys.polyerrors import (
    ExactQuotientFailed, DomainError,
)

from sympy.polys.domains import ZZ, QQ

from sympy import S, Symbol, symbols

from sympy.utilities.pytest import raises, skip, XFAIL

def test_sdp_LC():
    assert sdp_LC([], QQ) == QQ(0)
    assert sdp_LC([((1,0), QQ(1,2))], QQ) == QQ(1,2)
    assert sdp_LC([((1,1), QQ(1,4)), ((1,0), QQ(1,2))], QQ) == QQ(1,4)

def test_sdp_LM():
    assert sdp_LM([], 1) == (0, 0)
    assert sdp_LM([((1,0), QQ(1,2))], 1) == (1, 0)
    assert sdp_LM([((1,1), QQ(1,4)), ((1,0), QQ(1,2))], 1) == (1, 1)

def test_sdp_LT():
    assert sdp_LT([], 1, QQ) == ((0, 0), QQ(0))
    assert sdp_LT([((1,0), QQ(1,2))], 1, QQ) == ((1, 0), QQ(1,2))
    assert sdp_LT([((1,1), QQ(1,4)), ((1,0), QQ(1,2))], 1, QQ) == ((1, 1), QQ(1,4))

def test_sdp_del_LT():
    assert sdp_del_LT([]) == []
    assert sdp_del_LT([((1,0), QQ(1,2))]) == []
    assert sdp_del_LT([((1,1), QQ(1,4)), ((1,0), QQ(1,2))]) == [((1,0), QQ(1,2))]

def test_sdp_coeffs():
    assert sdp_coeffs([]) == []
    assert sdp_coeffs([((1,0), QQ(1,2))]) == [QQ(1,2)]
    assert sdp_coeffs([((1,1), QQ(1,4)), ((1,0), QQ(1,2))]) == [QQ(1,4), QQ(1,2)]

def test_sdp_monoms():
    assert sdp_monoms([]) == []
    assert sdp_monoms([((1,0), QQ(1,2))]) == [(1,0)]
    assert sdp_monoms([((1,1), QQ(1,4)), ((1,0), QQ(1,2))]) == [(1,1), (1,0)]

def test_sdp_sort():
    pass

def test_sdp_strip():
    assert sdp_strip([((2,2), 0), ((1,1), 1), ((0,0), 0)]) == [((1,1), 1)]

def test_sdp_normal():
    pass

def test_sdp_from_dict():
    pass

def test_sdp_indep_p():
    pass

def test_sdp_one_p():
    pass

def test_sdp_one():
    pass

def test_sdp_term_p():
    pass

def test_sdp_abs():
    pass

def test_sdp_neg():
    pass

def test_sdp_add_term():
    pass

def test_sdp_sub_term():
    pass

def test_sdp_mul_term():
    pass

def test_sdp_add():
    pass

def test_sdp_sub():
    pass

def test_sdp_mul():
    pass

def test_sdp_sqr():
    pass

def test_sdp_pow():
    f = sdp_from_dict({(1,): 2, (0,): 3}, grlex)

    assert sdp_pow(f, 0, 0, grlex, ZZ) == sdp_one(0, ZZ)
    assert sdp_pow(f, 1, 0, grlex, ZZ) == f

    assert sdp_pow(f, 2, 0, grlex, ZZ) == \
        sdp_from_dict({(2,): 4, (1,): 12, (0,): 9}, grlex)
    assert sdp_pow(f, 3, 0, grlex, ZZ) == \
        sdp_from_dict({(3,): 8, (2,): 36, (1,): 54, (0,): 27}, grlex)
    assert sdp_pow(f, 4, 0, grlex, ZZ) == \
        sdp_from_dict({(4,): 16, (3,): 96, (2,): 216, (1,): 216, (0,): 81}, grlex)
    assert sdp_pow(f, 5, 0, grlex, ZZ) == \
        sdp_from_dict({(5,): 32, (4,): 240, (3,): 720, (2,): 1080, (1,): 810, (0,): 243}, grlex)

    f = sdp_from_dict({(3,1,0): 1, (1,2,0): -2, (0,0,1): -3, (0,0,0): 1}, grlex)
    g = sdp_from_dict({(6,2,0): 1, (4,3,0): -4, (2,4,0): 4, (3,1,1): -6, (3,1,0): 2,
                      (1,2,1): 12, (1,2,0): -4, (0,0,2): 9, (0,0,1): -6, (0,0,0): 1}, grlex)

    assert sdp_pow(f, 2, 2, grlex, ZZ) == g

    raises(ValueError, lambda: sdp_pow(f, -2, 2, grlex, ZZ))

def test_sdp_monic():
    pass

def test_sdp_content():
    pass

def test_sdp_primitive():
    pass

def test_sdp_div():
    f = sdp_from_dict({(2,1): 4, (1,1): -2, (1,0): 4, (0,1): -2, (0,0): 8}, grlex)

    assert sdp_div(f, [sdp_from_dict({(0,0): 2}, grlex)], 1, grlex, ZZ) == \
        ([sdp_from_dict({(2,1): 2, (1,1): -1, (1,0): 2, (0,1): -1, (0,0): 4}, grlex)], [])

    assert sdp_div(f, [sdp_from_dict({(0,1): 2}, grlex)], 1, grlex, ZZ) == \
        ([sdp_from_dict({(2,0): 2, (1,0): -1, (0,0): -1}, grlex)],
          sdp_from_dict({(1,0): 4, (0,0): 8}, grlex))

    f = sdp_from_dict({(1,0): 1, (0,0): -1}, grlex)
    g = sdp_from_dict({(0,1): 1, (0,0): -1}, grlex)

    assert sdp_div(f, [g], 1, grlex, ZZ) == ([[]], f)

    f = sdp_from_dict({(3,): 1, (2,): -12, (0,): -42}, grlex)
    g = sdp_from_dict({(1,): 1, (0,): -3}, grlex)

    q = sdp_from_dict({(2,): 1, (1,): -9, (0,): -27}, grlex)
    r = sdp_from_dict({(0,): -123}, grlex)

    assert sdp_div(f, [g], 0, grlex, ZZ) == ([q], r)

    f = sdp_from_dict({(2,): QQ(1), (1,): QQ(2), (0,): QQ(2)}, grlex)

    g = sdp_from_dict({(0,): QQ(1)}, grlex)
    h = sdp_from_dict({(0,): QQ(2)}, grlex)

    q = sdp_from_dict({(2,): QQ(1,2), (1,): QQ(1), (0,): QQ(1)}, grlex)

    assert sdp_div(f, [g], 0, grlex, QQ) == ([f], [])
    assert sdp_div(f, [h], 0, grlex, QQ) == ([q], [])

    f = sdp_from_dict({(1,2): 1, (0,0): 1}, grlex)
    G = [sdp_from_dict({(1,1): 1, (0,0): 1}, grlex),
         sdp_from_dict({(0,1): 1, (0,0): 1}, grlex)]

    Q = [sdp_from_dict({(0,1): 1}, grlex),
         sdp_from_dict({(0,0): -1}, grlex)]
    r = sdp_from_dict({(0,0): 2}, grlex)

    assert sdp_div(f, G, 1, grlex, ZZ) == (Q, r)

    f = sdp_from_dict({(2,1): 1, (1,2): 1, (0,2): 1}, grlex)

    G = [sdp_from_dict({(1,1): 1, (0,0): -1}, grlex),
         sdp_from_dict({(0,2): 1, (0,0): -1}, grlex)]

    Q = [sdp_from_dict({(1,0): 1, (0,1): 1}, grlex),
         sdp_from_dict({(0,0): 1}, grlex)]
    r = sdp_from_dict({(1,0): 1, (0,1): 1, (0,0): 1}, grlex)

    assert sdp_div(f, G, 1, grlex, ZZ) == (Q, r)

    G = [sdp_from_dict({(0,2): 1, (0,0): -1}, grlex),
         sdp_from_dict({(1,1): 1, (0,0): -1}, grlex)]

    Q = [sdp_from_dict({(1,0): 1, (0,0): 1}, grlex),
         sdp_from_dict({(1,0): 1}, grlex)]
    r = sdp_from_dict({(1,0): 2, (0,0): 1}, grlex)

    assert sdp_div(f, G, 1, grlex, ZZ) == (Q, r)

def test_sdp_rem():
    f = sdp_from_dict({(2,1): 4, (1,1): -2, (1,0): 4, (0,1): -2, (0,0): 8}, grlex)

    assert sdp_rem(f, [sdp_from_dict({(0,0): 2}, grlex)], 1, grlex, ZZ) == []
    assert sdp_rem(f, [sdp_from_dict({(0,1): 2}, grlex)], 1, grlex, ZZ) == \
          sdp_from_dict({(1,0): 4, (0,0): 8}, grlex)

    f = sdp_from_dict({(1,0): 1, (0,0): -1}, grlex)
    g = sdp_from_dict({(0,1): 1, (0,0): -1}, grlex)

    assert sdp_rem(f, [g], 1, grlex, ZZ) == f

    f = sdp_from_dict({(3,): 1, (2,): -12, (0,): -42}, grlex)
    g = sdp_from_dict({(1,): 1, (0,): -3}, grlex)

    r = sdp_from_dict({(0,): -123}, grlex)

    assert sdp_rem(f, [g], 0, grlex, ZZ) == r

    f = sdp_from_dict({(1,2): 1, (0,0): 1}, grlex)
    G = [sdp_from_dict({(1,1): 1, (0,0): 1}, grlex),
         sdp_from_dict({(0,1): 1, (0,0): 1}, grlex)]

    r = sdp_from_dict({(0,0): 2}, grlex)

    assert sdp_rem(f, G, 1, grlex, ZZ) == r

    f = sdp_from_dict({(2,1): 1, (1,2): 1, (0,2): 1}, grlex)

    G = [sdp_from_dict({(1,1): 1, (0,0): -1}, grlex),
         sdp_from_dict({(0,2): 1, (0,0): -1}, grlex)]

    r = sdp_from_dict({(1,0): 1, (0,1): 1, (0,0): 1}, grlex)

    assert sdp_rem(f, G, 1, grlex, ZZ) == r

    G = [sdp_from_dict({(0,2): 1, (0,0): -1}, grlex),
         sdp_from_dict({(1,1): 1, (0,0): -1}, grlex)]

    r = sdp_from_dict({(1,0): 2, (0,0): 1}, grlex)

    assert sdp_rem(f, G, 1, grlex, ZZ) == r

def test_sdp_lcm():
    pass

def test_sdp_gcd():
    pass
