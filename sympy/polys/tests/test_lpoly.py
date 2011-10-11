"""Tests for distributed polynomials and series using lpoly"""

from sympy import Symbol, Rational, sympify, QQ, sqrt
from sympy.series import series
from sympy.core.function import expand
from sympy.polys.lpoly import LPoly, lgens, LPolySubs
from sympy.polys.monomialtools import lex, grlex
from sympy.functions.elementary.trigonometric import (cos, sin)
from sympy import Symbol
from sympy.functions.elementary.exponential import exp
from sympy.series.order import O
from sympy.polys.lpoly import monomial_from_sequence

def test_str():
    # str of a Poly object gives different output using QQ in python
    # or in gmpy mode, in one case giving n/1, in the other giving n
    # for QQ(n, 1); in this test these quantities do not appear
    lp = LPoly(list('xyz'), QQ, lex)
    x, y, z = lp.gens()
    p = lp('  +z^4 +1/2*z^2 -1/4')
    assert str(p) == ' +z^4 +1/2*z^2 -1/4'
    p = lp('z^4 -1/2')
    assert str(p) == ' +z^4 -1/2'

def test_read_monom():
    lp = LPoly(list('xyz'), QQ, lex)
    s = 'x^2*y^3*z'
    assert lp.read_monom(s) == monomial_from_sequence((2, 3, 1))

def test_mon_eval():
    lp = LPoly(list('xyz'), QQ, lex)
    s = '31*x^2*y^3*z'
    assert lp.mon_eval(s) == ((2, 3, 1), QQ(31))

def test_gens():
    lp = LPoly(list('xy'), QQ, lex)
    x, y = lp.gens()
    assert x == lp('x')
    #assert x**4 == lp('x^4')

def test_zero():
    lp = LPoly('x, y', QQ, lex)
    p = lp('0')
    assert str(p) == '0'
    assert p == lp(0)

def test_from_mon():
    lp = LPoly('x, y', QQ, lex)
    p = lp.from_mon((1, 4, QQ(7, 3)))
    assert p == lp('7/3*y^4')

def test_variables():
    lp = LPoly('z, y, x', QQ, lex)
    z, y, x = lp.gens()
    p = lp('x*y + 3*y^2')
    assert p.variables() == (1, 2)
    lp = LPoly('x, y, z', QQ, lex)
    x, y, z = lp.gens()
    p = lp('x*y + 3*y^2')
    assert p.variables() == (0, 1)

def test_eq():
    lp, x, y = lgens('x, y', QQ, lex)

    p1a = lp('1 + x')
    p2a = lp('x + 1')
    assert p1a == p2a
    p1 = 1 + x
    p2 = 1 + 2*x
    assert p1 == p1a
    p2 = p2 -x
    assert p1 == p2
    p2 -= x
    assert p2 == lp(1)
    assert p2 == 1

def test_coefficient():
    gens=['x%d' % i for i in range(11)]
    lp = LPoly(gens, QQ, lex)
    x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10 = lp.gens()
    p1 = x0**6 + 6*x0**5*x1 + 15*x0**4*x1**2 + 20*x0**3*x1**3 + 15*x0**2*x1**4 + 6*x0*x1**5 + x1**6 + x0**4*x1**2*x2
    m = x0**2
    p2 = p1.coefficient(m)
    assert p2 == 15*x1**4
    m = x0**2*x1**2
    p2 = p1.coefficient(m)
    assert p2 == lp(0)
    m = x0**2*x1**4
    p2 = p1.coefficient(m)
    assert p2 == lp(15)
    m = x0**4*x1**2
    p2 = p1.coefficient(m)
    assert p2 == x2+15

def test_add():
    lp, x, y = lgens('x, y', QQ, lex)
    zero = lp(0)
    one = lp(1)
    assert zero + zero == zero
    assert zero + one == one
    p1 = lp('-2*x')
    p2 = lp('2*x')
    assert p1+p2 == zero
    assert p1+p1 == lp('-4*x')

    p2 = lp('3/7*y^2')
    p3 = p1 + p2
    assert p3 == lp('-2*x+3/7*y^2')
    p4 = p3
    p3 = p3 + 1
    assert p3 == lp('-2*x+3/7*y^2 + 1')
    assert p4 == lp('-2*x+3/7*y^2')

def test_iadd():
    lp, x, y = lgens('x, y', QQ, lex)
    p3 = lp('-2*x+3/7*y^2')
    # check that p3 is mutable
    p4 = p3
    p3 += 1
    assert p3 == lp('-2*x+3/7*y^2 + 1')
    assert p4 == lp('-2*x+3/7*y^2 + 1')
    p3 = p3 - 1
    assert p3 == lp('-2*x+3/7*y^2')
    p3 -= lp('1')
    assert p3 == lp('-2*x+3/7*y^2-1')
    p3 -= 1
    assert p3 == lp('-2*x+3/7*y^2-2')
    assert p3 - lp('3/7*y^2') == lp('-2*x-2')
    p4 = 2 - p3
    assert p4 == lp('2*x-3/7*y^2+4')
    p4 -= 1
    assert p4 == lp('2*x-3/7*y^2+3')
    p4 -= lp('2*x')
    assert p4 == lp('-3/7*y^2+3')

def test_mul():
    lp, x, y = lgens('x, y', QQ, lex)
    p1 = x*x
    assert p1 * 0 == lp(0)
    #assert 0*p1 == lp(0)
    assert p1 * 1 == p1
    #assert 1 * p1 == p1
    p2 = lp('x + x^2')
    assert p1*p2 == lp('x^3 + x^4')
    assert 2*p2 == lp('2*x + 2*x^2')
    p3 = lp('y + 2/3*y^2')
    p4 = lp('2/3*x^2*y^2 + x^2*y + 2/3*x*y^2 + x*y')
    assert p2*p3 == p4

def test_square():
    lp, x, y = lgens('x, y', QQ, lex)
    p1 = x.square()
    assert p1 == lp('x^2')
    p2 = (x + y).square()
    assert p2 == x*x + 2*x*y + y*y
    p3 = lp('3/7*x + x^3*y + 8*y^2')
    assert p3.square() == lp('x^6*y^2 + 16*x^3*y^3 + 6/7*x^4*y + 64*y^4 + 48/7*x*y^2 + 9/49*x^2')


def test_mul_iadd():
    lp, x, y = lgens('x, y', QQ, lex)
    p = lp('x^2 + y^2 - 1')
    p1 = x + y
    p2 = x - y
    p.mul_iadd(p1, p2)
    assert p == 2*x*x - 1

def test_iadd_mon():
    lp, x, y = lgens('x, y', QQ, lex)
    p1 = x*x + y*y
    m = monomial_from_sequence((1, 2))
    p1.iadd_mon((m, QQ(2)))
    assert p1 == lp('x^2 + y^2 + 2*x*y^2')

def test_iadd_m_mul_q():
    lp, x, y = lgens('x, y', QQ, lex)
    p1 = lp('x^2 + y^2')
    p2 = lp('x^3 + y^3')
    m = x*y
    expva = m.keys()[0]
    c = m[expva]
    p1.iadd_m_mul_q(p2, (expva, c))
    assert p1 == lp('x^2 + y^2 + x^4*y + x*y^4')

def test_leading_expv():
    lp, x, y = lgens('x, y', QQ, lex)
    p = lp('2*x^2 + 3*y^3')
    assert p.leading_expv() == monomial_from_sequence((2, 0))
    lp, y, x = lgens('y, x', QQ, lex)
    p = lp('2*x^2 + 3*y^3 + x^2*y^2')
    assert p.leading_expv() == monomial_from_sequence((3, 0))
    lp, y, x = lgens('y, x', QQ, grlex)
    p = lp('2*x^2 + 3*y^3 + x^2*y^2')
    assert p.leading_expv() == monomial_from_sequence((2, 2))

def test_leading_term():
    lp, x, y = lgens('x, y', QQ, lex)
    p = lp('2*x^2 + 3*y^3')
    assert p.leading_term() == lp('2*x^2')

def test_div():
    lp, x, y = lgens('x, y', QQ, lex)
    p = y**2/4
    assert 4*p == y**2
    p1 = lp('x^2 + y^2 - 1')
    p2 = p1/7
    assert p2 == lp('1/7*x^2 + 1/7*y^2 - 1/7')

def test_pow():
    lp, x, y = lgens('x, y', QQ, lex)
    assert 4*(y**2/4) == y**2
    assert 4*(y**3/4) == y**3
    assert 4*((-y)**2/4) == y**2
    assert 4*((-y)**3/4) == -y**3
    p = 1 + x + y + x**2/2 + x*y - y**2/3 + x**3 + y**3 + x**2*y**2
    p2 = p*p
    p4 = p2**2
    assert p4 == p**4

def test_division():
    lp, x, y = lgens('x, y', QQ, lex)
    f = x**3
    f0 = x-y**2
    f1 = x-y
    qv, r = f.division((f0, f1))
    assert (qv[0], qv[1], r) == (x**2 + x*y**2 + y**4, lp(0), y**6)
    expv0 = f0.leading_expv()
    expv1 = f1.leading_expv()
    r1 = f.mod1(  ((expv0, f0), (expv1, f1))  )
    assert r == r1
    lp, y, x = lgens('y, x', QQ, lex)
    f = x**3*y**2 + x*y**4
    g = x**2 + y
    qv, r = f.division((g, ))
    assert r == x**9 +x**7
    f = x**10*y**2 + x**20*y**4
    qv, r = f.division((g, ))
    assert r == x**28 + x**14
    f = x**20*y**2 + 2*x**25*y**4
    qv, r = f.division((g, ))
    assert r == 2*x**33 + x**24
    f = x**61*y**2 + x**60*y**4
    g = x**2 + y
    qv, r = f.division((g, ))
    assert r == x**68 + x**65

def test_mul_trunc():
    lp, x, y, t = lgens('x, y, t', QQ, lex)
    p = 1 + t*x + t*y
    for i in range(2):
        p = p.mul_trunc(p, 't', 3)

    assert p == 6*x**2*t**2 + 12*x*y*t**2 + 6*y**2*t**2 + 4*x*t + 4*y*t + 1

def test_square_trunc():
    lp, x, y, t = lgens('x, y, t', QQ, lex)
    p = (1 + t*x + t*y)*2
    p1 = p.mul_trunc(p, 'x', 3)
    p2 = p.square_trunc('x', 3)
    assert p1 == p2

def test_trunc():
    lp, x, y, t = lgens('x, y, t', QQ, lex)
    p = (y + t*x)**4
    p1 = p.trunc('x', 3)
    assert p1 == y**4 + 4*y**3*t*x + 6*y**2*t**2*x**2

def test_pow_trunc():
    lp, x, y, z = lgens('x, y, z', QQ, lex)
    p0 = y + x*z
    p = p0**16
    for xx in ('x', 'y', 'z'):
        p1 = p.trunc(xx, 8)
        p2 = p0.pow_trunc(16, xx, 8)
        assert p1 == p2

def test_has_constant_term():
    lp, x, y, z = lgens('x, y, z', QQ, lex)
    p = y + x*z
    assert p.has_constant_term('x')
    p = x + x**4
    assert not p.has_constant_term('x')
    p = 1 + x + x**4
    assert p.has_constant_term('x')

def test_inversion():
    lp, x = lgens('x', QQ, lex)
    p = 2 + x + 2*x**2
    n = 5
    p1 = p.series_inversion('x', n)
    assert (p*p1).trunc('x', n) == lp(1)
    lp, x, y = lgens('x, y', QQ, grlex)
    p = 2 + x + 2*x**2 + y*x + x**2*y
    p1 = p.series_inversion('x', n)
    assert (p*p1).trunc('x', n) == lp(1)

def test_derivative():
    gens=['x%d' % i for i in range(11)]
    lp = LPoly(gens, QQ, lex)
    p = lp('288/5*x0^8*x1^6*x4^3*x10^2 + 8*x0^2*x2^3*x4^3 +2*x0^2-2*x1^2')
    p1 = p.derivative(0)
    assert p1 == lp('2304/5*x0^7*x1^6*x4^3*x10^2 + 16*x0*x2^3*x4^3 + 4*x0')
    p1 = p.derivative('x0')
    assert p1 == lp('2304/5*x0^7*x1^6*x4^3*x10^2 + 16*x0*x2^3*x4^3 + 4*x0')

def test_integrate():
    lp, x, y = lgens('x, y', QQ, lex)
    p = x**2 * y**5 + x + y
    p1 = p.integrate('x')
    assert p1 == lp('1/3*x^3*y^5 +1/2*x^2 + x*y')
    p1 = p.integrate('y')
    assert p1 == lp('1/6*x^2*y^6 + x*y + 1/2*y^2')

def test_nth_root():
    lp, x, y = lgens('x, y', QQ, lex)
    h = 4
    p0 = 1 + x**2 * y**5 + x
    p = p0.pow_trunc(4, 'x', h)
    p1 = p.nth_root(4, 'x', h)
    assert p1 == p0
    p = p0.pow_trunc(-4, 'x', h)
    p1 = p.nth_root(-4, 'x', h)
    assert p1 == p0

def test_sqrt():
    lp, x, y = lgens('x, y', QQ, lex)
    h = 4
    p0 = (1 + x**2 * y**5 + x)**3
    p0 = p0.trunc('x', h)
    p1 = p0.square_trunc('x', h)
    p2 = p1.sqrt('x', h)
    assert p2 == p0

def test_log():
    # log of univariate series
    lp, x = lgens('x', QQ, lex)
    p = 1 + x
    p1 = p.log('x', 4)
    assert p1 == x - x**2/2 + x**3/3
    p = lp('1+x+2/3*x^2')
    p1 = p.log('x', 9)
    assert p1 == lp('-17/648*x^8 + 13/189*x^7 - 11/162*x^6 - 1/45*x^5 + 7/36*x^4 - 1/3*x^3 + 1/6*x^2 + x')
    p2 = p.series_inversion('x', 9)
    p3 = p2.log('x', 9)
    assert p3 == -p1
    lp, x, y = lgens('x, y', QQ, lex)
    p = lp('1+x+2*y*x^2')
    p1 = p.log('x', 6)
    assert p1 == lp('4*x^5*y^2 - 2*x^5*y - 2*x^4*y^2 + 1/5*x^5 + 2*x^4*y - 1/4*x^4 - 2*x^3*y + 1/3*x^3 + 2*x^2*y - 1/2*x^2 + x')

def test_atan():
    lp, x, y = lgens('x, y', QQ, lex)
    p = x*y + x**2
    p1 = p.atan('x', 6)
    assert p1 == x**5*y**5/5 -x**5*y -x**4*y**2 -x**3*y**3/3 +x**2 +x*y

def test_atanh():
    lp, x, y = lgens('x, y', QQ, lex)
    p = x*y + x**2
    p1 = p.atanh('x', 6)
    assert p1 == lp('1/5*x^5*y^5 + x^5*y + x^4*y^2 + 1/3*x^3*y^3 + x^2 + x*y')

def test_tanh():
    lp, x, y = lgens('x, y', QQ, lex)
    p = x*y
    p1 = p.tanh('x', 8)
    assert p1 == lp('-17/315*x^7*y^7 + 2/15*x^5*y^5 - 1/3*x^3*y^3 + x*y')
    p = x*y + x

def test_exp():
    lp, x = lgens('x', QQ, lex)
    p = x + x**4
    for h in [10, 30]:
        p1 = p._exp_series0('x', h)
        p2 = p.exp('x', h)
        assert p1 == p2
        q = (1 + p).series_inversion('x', h) - 1
        p1 = q._exp_series0('x', h)
        p2 = q.exp('x', h)
        assert p1 == p2
        q1 = p1.log('x', h)
        assert q1 == q

def test_sin():
    lp, x, y = lgens('x, y', QQ, lex)
    h = 10
    p = 1 + x*y + x**4
    p = p.series_inversion('x', h) - 1
    p1 = p.sin('x', h)
    assert p1 == lp('224179/72576*x^9*y^9 - 1591/720*x^8*y^8 - 15/4*x^9*y^5 + 6931/5040*x^7*y^7 - 1/24*x^8*y^4 - 5/8*x^6*y^6 - 5/2*x^9*y + 2*x^7*y^3 - 1/120*x^5*y^5 + x^8 - 5/2*x^6*y^2 + 1/2*x^4*y^4 + 2*x^5*y - 5/6*x^3*y^3 - x^4 + x^2*y^2 - x*y')

def test_cos():
    lp, x, y = lgens('x, y', QQ, lex)
    h = 10
    p = 1 + x*y + x**4
    p = p.series_inversion('x', h) - 1
    p1 = p.cos('x', h)
    assert p1 == lp('8791/5040*x^9*y^9 - 16699/8064*x^8*y^8 - 1501/120*x^9*y^5 + 87/40*x^7*y^7 + 55/6*x^8*y^4 - 1501/720*x^6*y^6 + 3*x^9*y - 35/6*x^7*y^3 + 11/6*x^5*y^5 - 1/2*x^8 + 3*x^6*y^2 - 35/24*x^4*y^4 - x^5*y + x^3*y^3 - 1/2*x^2*y^2 + 1')

def test_cos_sin():
    lp, x, y = lgens('x, y', QQ, lex)
    h = 10
    p = 1 + x*y + x**4
    p = p.series_inversion('x', h) - 1
    assert p.cos_sin('x', h) == (p.cos('x', h), p.sin('x', h))

def test_sinh():
    lp, x, y = lgens('x, y', QQ, lex)
    p = x*y
    p1 = p.sinh('x', 8)
    assert p1 == lp('1/5040*x^7*y^7 + 1/120*x^5*y^5 + 1/6*x^3*y^3 + x*y')

def test_cosh():
    lp, x, y = lgens('x, y', QQ, lex)
    p = x*y
    p1 = p.cosh('x', 8)
    assert p1 == lp('1/720*x^6*y^6 + 1/24*x^4*y^4 + 1/2*x^2*y^2 + 1')

def test_lambert():
    lp, x = lgens('x', QQ, lex)
    p1 = x.lambert('x', 10)
    assert p1 == lp('531441/4480*x^9 - 16384/315*x^8 + 16807/720*x^7 - 54/5*x^6 + 125/24*x^5 - 8/3*x^4 + 3/2*x^3 - x^2 + x')

def test_asin():
    lp, x = lgens('x', QQ, lex)
    p1 = x.asin('x', 12)
    assert p1 == lp('63/2816*x^11 + 35/1152*x^9 + 5/112*x^7 + 3/40*x^5 + 1/6*x^3 + x')

def test_asinh():
    lp, x = lgens('x', QQ, lex)
    p1 = x.asinh('x', 12)
    assert p1 == lp('-63/2816*x^11 + 35/1152*x^9 - 5/112*x^7 + 3/40*x^5 - 1/6*x^3 + x')

def test_basic():
    lp, _x, _y = lgens('_x, _y', QQ, lex)
    x = Symbol('x')
    y = Symbol('y')
    p = (_x + _y**2/3)**2
    p2 = p.tobasic(x, y)
    assert p2 == expand((x + y**2/3)**2)
    x = Symbol('x')
    a = Symbol('a')
    lp, A, X = lgens('A, X', QQ, lex)
    h = 5
    p1 = (1 + A*X + 2*A*X**2).nth_root(-2, 'X', h) + X
    p2 = p1.sqrt('X', h)
    assert p2 == 1 + X/2 - A*X/4 - X**2/8 - 3*A*X**2/8 + 5*A**2*X**2/32 + X**3/16 - 15*A**3*X**3/128 + 5*A*X**3/32 + 37*A**2*X**3/64 - 5*X**4/128 - 173*A**3*X**4/256 - 7*A*X**4/64 + 115*A**2*X**4/256 + 195*A**4*X**4/2048
    p2 = p2.tobasic(a, x)
    assert p2 == 1 + x/2 - a*x/4 - x**2/8 - 3*a*x**2/8 + 5*a**2*x**2/32 + x**3/16 - 15*a**3*x**3/128 + 5*a*x**3/32 + 37*a**2*x**3/64 - 5*x**4/128 - 173*a**3*x**4/256 - 7*a*x**4/64 + 115*a**2*x**4/256 + 195*a**4*x**4/2048


def test_SR1():
    x = Symbol('x')
    lp, _x  = lgens('_x', sympify, lex)
    h = 10
    p=(_x*sqrt(2) + _x**2*sqrt(3)).cos('_x', h).pow_trunc(-1, '_x', h)
    p1 = p.tobasic(x)
    p2 = series(1/cos(x*sqrt(2) + x**2*sqrt(3)), x, 0, h)
    p3 = p1 - p2
    assert p3 == O(x**h)

def test_SR2():
    x = Symbol('x')
    lp, _x  = lgens('_x', sympify, lex)
    h = 5
    p = (_x + 1).exp('_x', h)
    p1 = p.tobasic(x)
    assert series(exp(x+1), x, 0, 5) == p1 + O(x**5)

def test_subs():
    lp, x, y, z = lgens('x, y, z', QQ, lex)
    p1 = x**2 + 2*y**2 + 3*z**2
    p2 = p1.subs(z=x+y)
    p3 = x**2 + 2*y**2 + 3*(x+y)**2
    assert p2 == p3

    lp, x, y = lgens('x, y', QQ, lex)
    p = (1+x)**4
    p1 = p.subs(x=(1+y)**5)
    assert p1.coefficient(y**20) == lp('1')

    p = 1 + x**2 + 2*y**2
    p1 = p.subs(x=1+y**2, y=1+x)
    assert p1 == 1 + (1+y**2)**2 + 2*(1+x)**2

def test_subs_trunc():
    lp, x, y = lgens('x, y', QQ, lex)
    p1 = (1+x)**4 + 2*x**3 * y**4
    p2 = p1.subs_trunc('x', 5, y=x+y)
    p3 = ((1+x)**4 + 2*x**3 * (x+y)**4).trunc('x', 5)
    assert p2 == p3
    p2 = p1.subs_trunc('x', 7, x=x**2, y=x+y)
    p3 = ((1+x**2)**4 + 2*x**6 * (x+y)**4).trunc('x', 7)
    assert p2 == p3

def test_LPolySubs():
    lp = LPoly(['x%d' % i for i in range(11)], QQ, lex)
    x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10 = lp.gens()
    lp1, c = lgens('c', QQ, lex)
    rules = {'x0':c, 'x1':c+1, 'x2':c**2+1, 'x3':c+2, 'x4':c**4, 'x5':c+1,
        'x6':c*(c-3), 'x7':c*(c-7), 'x8':c+7, 'x9':c+9, 'x10':c+10}
    sb = LPolySubs(lp, lp1, rules)
    p1 = x0**2*x1 + 2*x0**2 + 3*x1+ 7*x3*x4*x5 + x6*x7*x8
    p = sb.subs(p1)
    assert p == 7*c**6 + 22*c**5 + 11*c**4 - 48*c**3 + 150*c**2 + 3*c + 3

    p = sb.subs_trunc(p1, 'c', 3)
    assert p == 150*c**2 +3*c +3

    lp, x, y, z = lgens('x, y, z', QQ, lex)
    p1 = x**2 + 2*y**2 + 3*z**2
    rules = {'z': x+y}
    sb = LPolySubs(lp, lp, rules)
    p2 = sb.subs(p1)
    p3 = x**2 + 2*y**2 + 3*(x+y)**2
    assert p2 == p3

    lp, x, y = lgens('x, y', QQ, lex)
    rules = {'x':(1+y)**5}
    sb = LPolySubs(lp, lp, rules)
    p = (1+x)**4
    p1 = sb.subs(p)
    assert p1.coefficient(y**20) == lp(1)
    assert p1.coefficient(y**18) == lp(190)

    sb = LPolySubs(lp, lp, {'x':1+y**2, 'y':1+x})
    p = 1 + x**2 + 2*y**2
    p1 = sb.subs(p)
    assert p1 == 1 + (1+y**2)**2 + 2*(1+x)**2
