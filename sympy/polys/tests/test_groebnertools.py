"""Tests for sparse distributed polynomials and Groebner bases. """

from sympy.polys.groebnertools import (
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
    sdp_groebner,
)

from sympy.polys.monomialtools import (
    monomial_lex_key as O_lex,
    monomial_grlex_key as O_grlex,
    monomial_grevlex_key as O_grevlex,
)

from sympy.polys.polyerrors import (
    ExactQuotientFailed, DomainError,
)

from sympy.polys.domains import ZZ, QQ

from sympy import S, Symbol, symbols, groebner

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
    f = sdp_from_dict({(1,): 2, (0,): 3}, O_grlex)

    assert sdp_pow(f, 0, 0, O_grlex, ZZ) == sdp_one(0, ZZ)
    assert sdp_pow(f, 1, 0, O_grlex, ZZ) == f

    assert sdp_pow(f, 2, 0, O_grlex, ZZ) == \
        sdp_from_dict({(2,): 4, (1,): 12, (0,): 9}, O_grlex)
    assert sdp_pow(f, 3, 0, O_grlex, ZZ) == \
        sdp_from_dict({(3,): 8, (2,): 36, (1,): 54, (0,): 27}, O_grlex)
    assert sdp_pow(f, 4, 0, O_grlex, ZZ) == \
        sdp_from_dict({(4,): 16, (3,): 96, (2,): 216, (1,): 216, (0,): 81}, O_grlex)
    assert sdp_pow(f, 5, 0, O_grlex, ZZ) == \
        sdp_from_dict({(5,): 32, (4,): 240, (3,): 720, (2,): 1080, (1,): 810, (0,): 243}, O_grlex)

    f = sdp_from_dict({(3,1,0): 1, (1,2,0): -2, (0,0,1): -3, (0,0,0): 1}, O_grlex)
    g = sdp_from_dict({(6,2,0): 1, (4,3,0): -4, (2,4,0): 4, (3,1,1): -6, (3,1,0): 2,
                      (1,2,1): 12, (1,2,0): -4, (0,0,2): 9, (0,0,1): -6, (0,0,0): 1}, O_grlex)

    assert sdp_pow(f, 2, 2, O_grlex, ZZ) == g

    raises(ValueError, "sdp_pow(f, -2, 2, O_grlex, ZZ)")

def test_sdp_monic():
    pass

def test_sdp_content():
    pass

def test_sdp_primitive():
    pass

def test_sdp_div():
    f = sdp_from_dict({(2,1): 4, (1,1): -2, (1,0): 4, (0,1): -2, (0,0): 8}, O_grlex)

    assert sdp_div(f, [sdp_from_dict({(0,0): 2}, O_grlex)], 1, O_grlex, ZZ) == \
        ([sdp_from_dict({(2,1): 2, (1,1): -1, (1,0): 2, (0,1): -1, (0,0): 4}, O_grlex)], [])

    assert sdp_div(f, [sdp_from_dict({(0,1): 2}, O_grlex)], 1, O_grlex, ZZ) == \
        ([sdp_from_dict({(2,0): 2, (1,0): -1, (0,0): -1}, O_grlex)],
          sdp_from_dict({(1,0): 4, (0,0): 8}, O_grlex))

    f = sdp_from_dict({(1,0): 1, (0,0): -1}, O_grlex)
    g = sdp_from_dict({(0,1): 1, (0,0): -1}, O_grlex)

    assert sdp_div(f, [g], 1, O_grlex, ZZ) == ([[]], f)

    f = sdp_from_dict({(3,): 1, (2,): -12, (0,): -42}, O_grlex)
    g = sdp_from_dict({(1,): 1, (0,): -3}, O_grlex)

    q = sdp_from_dict({(2,): 1, (1,): -9, (0,): -27}, O_grlex)
    r = sdp_from_dict({(0,): -123}, O_grlex)

    assert sdp_div(f, [g], 0, O_grlex, ZZ) == ([q], r)

    f = sdp_from_dict({(2,): QQ(1), (1,): QQ(2), (0,): QQ(2)}, O_grlex)

    g = sdp_from_dict({(0,): QQ(1)}, O_grlex)
    h = sdp_from_dict({(0,): QQ(2)}, O_grlex)

    q = sdp_from_dict({(2,): QQ(1,2), (1,): QQ(1), (0,): QQ(1)}, O_grlex)

    assert sdp_div(f, [g], 0, O_grlex, QQ) == ([f], [])
    assert sdp_div(f, [h], 0, O_grlex, QQ) == ([q], [])

    f = sdp_from_dict({(1,2): 1, (0,0): 1}, O_grlex)
    G = [sdp_from_dict({(1,1): 1, (0,0): 1}, O_grlex),
         sdp_from_dict({(0,1): 1, (0,0): 1}, O_grlex)]

    Q = [sdp_from_dict({(0,1): 1}, O_grlex),
         sdp_from_dict({(0,0): -1}, O_grlex)]
    r = sdp_from_dict({(0,0): 2}, O_grlex)

    assert sdp_div(f, G, 1, O_grlex, ZZ) == (Q, r)

    f = sdp_from_dict({(2,1): 1, (1,2): 1, (0,2): 1}, O_grlex)

    G = [sdp_from_dict({(1,1): 1, (0,0): -1}, O_grlex),
         sdp_from_dict({(0,2): 1, (0,0): -1}, O_grlex)]

    Q = [sdp_from_dict({(1,0): 1, (0,1): 1}, O_grlex),
         sdp_from_dict({(0,0): 1}, O_grlex)]
    r = sdp_from_dict({(1,0): 1, (0,1): 1, (0,0): 1}, O_grlex)

    assert sdp_div(f, G, 1, O_grlex, ZZ) == (Q, r)

    G = [sdp_from_dict({(0,2): 1, (0,0): -1}, O_grlex),
         sdp_from_dict({(1,1): 1, (0,0): -1}, O_grlex)]

    Q = [sdp_from_dict({(1,0): 1, (0,0): 1}, O_grlex),
         sdp_from_dict({(1,0): 1}, O_grlex)]
    r = sdp_from_dict({(1,0): 2, (0,0): 1}, O_grlex)

    assert sdp_div(f, G, 1, O_grlex, ZZ) == (Q, r)

def test_sdp_lcm():
    pass

def test_sdp_gcd():
    pass

def test_sdp_groebner():
    f = sdp_from_dict({(1,2): QQ(2,), (2,0): QQ(1)}, O_lex)
    g = sdp_from_dict({(0,3): QQ(2), (1,1): QQ(1), (0,0): QQ(-1)},  O_lex)

    a = sdp_from_dict({(1,0): QQ(1,1)}, O_lex)
    b = sdp_from_dict({(0,3): QQ(1,1), (0,0): QQ(-1,2)}, O_lex)

    assert sdp_groebner((f, g), 1, O_lex, QQ) == [a, b]

    f = sdp_from_dict({(2,1): QQ(2,), (0,2): QQ(1)}, O_lex)
    g = sdp_from_dict({(3,0): QQ(2), (1,1): QQ(1), (0,0): QQ(-1)},  O_lex)

    a = sdp_from_dict({(0,1): QQ(1,1)}, O_lex)
    b = sdp_from_dict({(3,0): QQ(1,1), (0,0): QQ(-1,2)}, O_lex)

    assert sdp_groebner((f, g), 1, O_lex, QQ) == [b, a]

    f = sdp_from_dict({(0,0,2): QQ(-1), (1,0,0): QQ(1)}, O_lex)
    g = sdp_from_dict({(0,0,3): QQ(-1), (0,1,0): QQ(1)}, O_lex)

    assert sdp_groebner((f, g), 1, O_lex, QQ) == [f, g]

    f = sdp_from_dict({(3,0): QQ(1), (1,1): QQ(-2)}, O_grlex)
    g = sdp_from_dict({(2,1): QQ(1), (0,2): QQ(-2), (1,0): QQ(1)}, O_grlex)

    a = sdp_from_dict({(2,0): QQ(1)}, O_grlex)
    b = sdp_from_dict({(1,1): QQ(1)}, O_grlex)
    c = sdp_from_dict({(0,2): QQ(1), (1, 0): QQ(-1,2)}, O_grlex)

    assert sdp_groebner((f, g), 1, O_grlex, QQ) == [a, b, c]

    f = sdp_from_dict({(2,0,0): -QQ(1), (0,1,0): QQ(1)}, O_lex)
    g = sdp_from_dict({(3,0,0): -QQ(1), (0,0,1): QQ(1)}, O_lex)

    assert sdp_groebner((f, g), 2, O_lex, QQ) == [
        sdp_from_dict({(2,0,0): QQ(1), (0,1,0): -QQ(1)}, O_lex),
        sdp_from_dict({(1,1,0): QQ(1), (0,0,1): -QQ(1)}, O_lex),
        sdp_from_dict({(1,0,1): QQ(1), (0,2,0): -QQ(1)}, O_lex),
        sdp_from_dict({(0,3,0): QQ(1), (0,0,2): -QQ(1)}, O_lex),
    ]

    f = sdp_from_dict({(2,0,0): -QQ(1), (0,1,0): QQ(1)}, O_grlex)
    g = sdp_from_dict({(3,0,0): -QQ(1), (0,0,1): QQ(1)}, O_grlex)

    assert sdp_groebner((f, g), 2, O_grlex, QQ) == [
        sdp_from_dict({(0,3,0): QQ(1), (0,0,2): -QQ(1)}, O_grlex),
        sdp_from_dict({(2,0,0): QQ(1), (0,1,0): -QQ(1)}, O_grlex),
        sdp_from_dict({(1,1,0): QQ(1), (0,0,1): -QQ(1)}, O_grlex),
        sdp_from_dict({(1,0,1): QQ(1), (0,2,0): -QQ(1)}, O_grlex),
    ]

    f = sdp_from_dict({(2,0,0): -QQ(1), (0,0,1): QQ(1)}, O_lex)
    g = sdp_from_dict({(3,0,0): -QQ(1), (0,1,0): QQ(1)}, O_lex)

    assert sdp_groebner((f, g), 2, O_lex, QQ) == [
        sdp_from_dict({(2,0,0): QQ(1), (0,0,1): -QQ(1)}, O_lex),
        sdp_from_dict({(1,1,0): QQ(1), (0,0,2): -QQ(1)}, O_lex),
        sdp_from_dict({(1,0,1): QQ(1), (0,1,0): -QQ(1)}, O_lex),
        sdp_from_dict({(0,2,0): QQ(1), (0,0,3): -QQ(1)}, O_lex),
    ]

    f = sdp_from_dict({(2,0,0): -QQ(1), (0,0,1): QQ(1)}, O_grlex)
    g = sdp_from_dict({(3,0,0): -QQ(1), (0,1,0): QQ(1)}, O_grlex)

    assert sdp_groebner((f, g), 2, O_grlex, QQ) == [
        sdp_from_dict({(0,0,3): QQ(1), (0,2,0): -QQ(1)}, O_grlex),
        sdp_from_dict({(2,0,0): QQ(1), (0,0,1): -QQ(1)}, O_grlex),
        sdp_from_dict({(1,1,0): QQ(1), (0,0,2): -QQ(1)}, O_grlex),
        sdp_from_dict({(1,0,1): QQ(1), (0,1,0): -QQ(1)}, O_grlex),
    ]

    f = sdp_from_dict({(0,2,0): -QQ(1), (1,0,0): QQ(1)}, O_lex)
    g = sdp_from_dict({(0,3,0): -QQ(1), (0,0,1): QQ(1)}, O_lex)

    assert sdp_groebner((f, g), 2, O_lex, QQ) == [
        sdp_from_dict({(1,0,0): QQ(1), (0,2,0): -QQ(1)}, O_lex),
        sdp_from_dict({(0,3,0): QQ(1), (0,0,1): -QQ(1)}, O_lex),
    ]

    f = sdp_from_dict({(0,2,0): -QQ(1), (1,0,0): QQ(1)}, O_grlex)
    g = sdp_from_dict({(0,3,0): -QQ(1), (0,0,1): QQ(1)}, O_grlex)

    assert sdp_groebner((f, g), 2, O_grlex, QQ) == [
        sdp_from_dict({(2,0,0): QQ(1), (0,1,1): -QQ(1)}, O_grlex),
        sdp_from_dict({(1,1,0): QQ(1), (0,0,1): -QQ(1)}, O_grlex),
        sdp_from_dict({(0,2,0): QQ(1), (1,0,0): -QQ(1)}, O_grlex),
    ]

    f = sdp_from_dict({(0,0,2): -QQ(1), (1,0,0): QQ(1)}, O_lex)
    g = sdp_from_dict({(0,0,3): -QQ(1), (0,1,0): QQ(1)}, O_lex)

    assert sdp_groebner((f, g), 2, O_lex, QQ) == [
        sdp_from_dict({(1,0,0): QQ(1), (0,0,2): -QQ(1)}, O_lex),
        sdp_from_dict({(0,1,0): QQ(1), (0,0,3): -QQ(1)}, O_lex),
    ]

    f = sdp_from_dict({(0,0,2): -QQ(1), (1,0,0): QQ(1)}, O_grlex)
    g = sdp_from_dict({(0,0,3): -QQ(1), (0,1,0): QQ(1)}, O_grlex)

    assert sdp_groebner((f, g), 2, O_grlex, QQ) == [
        sdp_from_dict({(2,0,0): QQ(1), (0,1,1): -QQ(1)}, O_grlex),
        sdp_from_dict({(1,0,1): QQ(1), (0,1,0): -QQ(1)}, O_grlex),
        sdp_from_dict({(0,0,2): QQ(1), (1,0,0): -QQ(1)}, O_grlex),
    ]

    f = sdp_from_dict({(0,2,0): -QQ(1), (0,0,1): QQ(1)}, O_lex)
    g = sdp_from_dict({(0,3,0): -QQ(1), (1,0,0): QQ(1)}, O_lex)

    assert sdp_groebner((f, g), 2, O_lex, QQ) == [
        sdp_from_dict({(1,0,0): QQ(1), (0,1,1): -QQ(1)}, O_lex),
        sdp_from_dict({(0,2,0): QQ(1), (0,0,1): -QQ(1)}, O_lex),
    ]

    f = sdp_from_dict({(0,2,0): -QQ(1), (0,0,1): QQ(1)}, O_grlex)
    g = sdp_from_dict({(0,3,0): -QQ(1), (1,0,0): QQ(1)}, O_grlex)

    assert sdp_groebner((f, g), 2, O_grlex, QQ) == [
        sdp_from_dict({(0,0,3): QQ(1), (2,0,0): -QQ(1)}, O_grlex),
        sdp_from_dict({(1,1,0): QQ(1), (0,0,2): -QQ(1)}, O_grlex),
        sdp_from_dict({(0,2,0): QQ(1), (0,0,1): -QQ(1)}, O_grlex),
        sdp_from_dict({(0,1,1): QQ(1), (1,0,0): -QQ(1)}, O_grlex),
    ]

    f = sdp_from_dict({(0,0,2): -QQ(1), (0,1,0): QQ(1)}, O_lex)
    g = sdp_from_dict({(0,0,3): -QQ(1), (1,0,0): QQ(1)}, O_lex)

    assert sdp_groebner((f, g), 2, O_lex, QQ) == [
        sdp_from_dict({(1,0,0): QQ(1), (0,0,3): -QQ(1)}, O_lex),
        sdp_from_dict({(0,1,0): QQ(1), (0,0,2): -QQ(1)}, O_lex),
    ]

    f = sdp_from_dict({(0,0,2): -QQ(1), (0,1,0): QQ(1)}, O_grlex)
    g = sdp_from_dict({(0,0,3): -QQ(1), (1,0,0): QQ(1)}, O_grlex)

    assert sdp_groebner((f, g), 2, O_grlex, QQ) == [
        sdp_from_dict({(0,3,0): QQ(1), (2,0,0): -QQ(1)}, O_grlex),
        sdp_from_dict({(1,0,1): QQ(1), (0,2,0): -QQ(1)}, O_grlex),
        sdp_from_dict({(0,1,1): QQ(1), (1,0,0): -QQ(1)}, O_grlex),
        sdp_from_dict({(0,0,2): QQ(1), (0,1,0): -QQ(1)}, O_grlex),
    ]

    f = sdp_from_dict({(2,2): QQ(4), (1,1): QQ(4), (0,0): QQ(1)}, O_lex)
    g = sdp_from_dict({(2,0): QQ(1), (0,2): QQ(1), (0,0):-QQ(1)}, O_lex)

    assert sdp_groebner((f, g), 1, O_lex, QQ) == [
        sdp_from_dict({(1,0): QQ(1,1), (0,7): QQ(-4,1), (0,5): QQ(8,1), (0,3): QQ(-7,1), (0,1): QQ(3,1)}, O_lex),
        sdp_from_dict({(0,8): QQ(1,1), (0,6): QQ(-2,1), (0,4): QQ(3,2), (0,2): QQ(-1,2), (0,0): QQ(1,16)}, O_lex),
    ]

    raises(DomainError, "sdp_groebner([], 1, O_lex, ZZ)")

def test_benchmark_minpoly():
    x, y, z = symbols('x,y,z')

    I = [x**3+x+1, y**2+y+1, (x+y)*z-(x**2+y)]

    assert groebner(I, x, y, z, order='lex') == [
        -975 + 2067*x + 6878*z - 11061*z**2 + 6062*z**3 - 1065*z**4 + 155*z**5,
        -308 + 159*y + 1043*z - 1161*z**2 + 523*z**3 - 91*z**4 + 12*z**5,
        13 - 46*z + 89*z**2 - 82*z**3 + 41*z**4 - 7*z**5 + z**6,
    ]

    assert groebner(I, x, y, z, order='lex', field=True) == [
        -S(25)/53 + x + 6878*z/2067 - 3687*z**2/689 + 6062*z**3/2067 - 355*z**4/689 + 155*z**5/2067,
        -S(308)/159 + y + 1043*z/159 - 387*z**2/53 + 523*z**3/159 - 91*z**4/159 + 4*z**5/53,
        13 - 46*z + 89*z**2 - 82*z**3 + 41*z**4 - 7*z**5 + z**6,
    ]

@XFAIL
def test_benchmark_coloring():
    skip('takes too much time')

    V = range(1, 12+1)
    E = [(1,2),(2,3),(1,4),(1,6),(1,12),(2,5),(2,7),(3,8),(3,10),
         (4,11),(4,9),(5,6),(6,7),(7,8),(8,9),(9,10),(10,11),
         (11,12),(5,12),(5,9),(6,10),(7,11),(8,12),(3,4)]

    V = [Symbol('x' + str(i)) for i in V]
    E = [(V[i-1], V[j-1]) for i, j in E]

    x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12 = V

    I3 = [x**3 - 1 for x in V]
    Ig = [x**2 + x*y + y**2 for x, y in E]

    I = I3 + Ig

    assert groebner(I[:-1], V, order='lex') == [
        x1 + x11 + x12,
        x2 - x11,
        x3 - x12,
        x4 - x12,
        x5 + x11 + x12,
        x6 - x11,
        x7 - x12,
        x8 + x11 + x12,
        x9 - x11,
        x10 + x11 + x12,
        x11**2 + x11*x12 + x12**2,
        x12**3 - 1,
    ]

    assert groebner(I, V, order='lex') == [1]

def test_benchmark_katsura_3():
    x0, x1, x2 = symbols('x:3')

    I = [x0 + 2*x1 + 2*x2 - 1,
         x0**2 + 2*x1**2 + 2*x2**2 - x0,
         2*x0*x1 + 2*x1*x2 - x1]

    assert groebner(I, x0, x1, x2, order='lex') == [
        -7 + 7*x0 + 8*x2 + 158*x2**2 - 420*x2**3,
        7*x1 + 3*x2 - 79*x2**2 + 210*x2**3,
        x2 + x2**2 - 40*x2**3 + 84*x2**4,
    ]

    assert groebner(I, x0, x1, x2, order='grlex') == [
        7*x1 + 3*x2 - 79*x2**2 + 210*x2**3,
        -x1 + x2 - 3*x2**2 + 5*x1**2,
        -x1 - 4*x2 + 10*x1*x2 + 12*x2**2,
        -1 + x0 + 2*x1 + 2*x2,
    ]

def test_benchmark_katsura_4():
    x0, x1, x2, x3 = symbols('x:4')

    I = [x0 + 2*x1 + 2*x2 + 2*x3 - 1,
         x0**2 + 2*x1**2 + 2*x2**2 + 2*x3**2 - x0,
         2*x0*x1 + 2*x1*x2 + 2*x2*x3 - x1,
         x1**2 + 2*x0*x2 + 2*x1*x3 - x2]

    assert groebner(I, x0, x1, x2, x3, order='lex') == [
        5913075*x0 - 159690237696*x3**7 + 31246269696*x3**6 + 27439610544*x3**5 - 6475723368*x3**4 - 838935856*x3**3 + 275119624*x3**2 + 4884038*x3 - 5913075,
        1971025*x1 - 97197721632*x3**7 + 73975630752*x3**6 - 12121915032*x3**5 - 2760941496*x3**4 + 814792828*x3**3 - 1678512*x3**2 - 9158924*x3,
        5913075*x2 + 371438283744*x3**7 - 237550027104*x3**6 + 22645939824*x3**5 + 11520686172*x3**4 - 2024910556*x3**3 - 132524276*x3**2 + 30947828*x3,
        128304*x3**8 - 93312*x3**7 + 15552*x3**6 + 3144*x3**5 - 1120*x3**4 + 36*x3**3 + 15*x3**2 - x3,
    ]

    assert groebner(I, x0, x1, x2, x3, order='grlex') == [
        393*x1 - 4662*x2**2 + 4462*x2*x3 - 59*x2 + 224532*x3**4 - 91224*x3**3 - 678*x3**2 + 2046*x3,
        -x1 + 196*x2**3 - 21*x2**2 + 60*x2*x3 - 18*x2 - 168*x3**3 + 83*x3**2 - 9*x3,
        -6*x1 + 1134*x2**2*x3 - 189*x2**2 - 466*x2*x3 + 32*x2 - 630*x3**3 + 57*x3**2 + 51*x3,
        33*x1 + 63*x2**2 + 2268*x2*x3**2 - 188*x2*x3 + 34*x2 + 2520*x3**3 - 849*x3**2 + 3*x3,
        7*x1**2 - x1 - 7*x2**2 - 24*x2*x3 + 3*x2 - 15*x3**2 + 5*x3,
        14*x1*x2 - x1 + 14*x2**2 + 18*x2*x3 - 4*x2 + 6*x3**2 - 2*x3,
        14*x1*x3 - x1 + 7*x2**2 + 32*x2*x3 - 4*x2 + 27*x3**2 - 9*x3,
        x0 + 2*x1 + 2*x2 + 2*x3 - 1,
    ]

