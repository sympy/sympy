"""Tests for distributed polynomials and series using lpoly"""

from sympy import Symbol, S, symbols, Rational, sympify, sqrt
from sympy.polys.domains import QQ, PythonRationalType
from sympy.series import series
from sympy.core.function import expand
from sympy.polys.lpoly import LPoly, MLPoly, lgens, mlgens, nclgens, LPolySubs, monomial_as_expr, TaylorEvalError, PythonRationalType_new
from sympy.polys.monomialtools import lex, grlex
from sympy.functions.elementary.trigonometric import cos, sin
from sympy.functions.elementary.exponential import exp, log
from sympy.series.order import O
from sympy.polys.lpoly import monomial_from_sequence
from sympy.physics.quantum import Operator
from sympy.utilities.pytest import raises, XFAIL
from sympy import I

def test_PythonRationalType_new():
    p = PythonRationalType(4, 3)
    p1 = PythonRationalType_new(p)
    assert p1 == p
    p1 = PythonRationalType_new('4/3')
    assert p1 == p
    p1 = PythonRationalType_new('4')
    assert p1 == 4
    p1 = PythonRationalType_new(4)
    assert p1 == 4

def test_str():
    # str of a LPolyElement object gives different output using QQ in python
    # or in gmpy mode, in one case giving n/1, in the other giving n
    # for QQ(n, 1); in this test these quantities do not appear
    lp = LPoly(list('xyz'), QQ, lex)
    x, y, z = lp.gens
    p = lp('  +z**4 +1/2*z**2 -1/4')
    assert str(p) == 'z**4 + 1/2*z**2 - 1/4'
    p = lp('z**4 -1/2')
    assert str(p) == 'z**4 - 1/2'
    p = lp('-z - 1')
    assert p == -z - 1
    p = lp('-z - 1 -z')
    assert p == -2*z - 1
    p = -x
    assert str(p) == '-x'
    assert +x == x
    lp2 = LPoly('w', QQ, lex)
    def test1(p):
        w = lp2.gens
        p2 = lp2(p)
    raises(NotImplementedError, 'test1(p)')

def test_read_monom():
    lp = LPoly(list('xyz'), QQ, lex)
    s = 'x**2*y**3*z'
    assert lp.read_monom(s) == monomial_from_sequence((2, 3, 1))

def test_mon_eval():
    lp = LPoly(list('xyz'), QQ, lex)
    s = '31*x**2*y**3*z'
    assert lp.mon_eval(s) == ((2, 3, 1), QQ(31))

def test_monomial_as_expr():
    x, y = symbols('x,y')
    assert monomial_as_expr((1,), x) == x
    assert monomial_as_expr((1, 2), x, y) == x*y**2

def test_gens():
    lp = LPoly(list('xy'), QQ, lex)
    x, y = lp.gens
    assert x == lp('x')
    #assert x**4 == lp('x**4')

def test_zero():
    lp = LPoly('x, y', QQ, lex)
    p = lp('0')
    assert str(p) == '0'
    assert p == lp(0)

def test_from_mon():
    lp, x, y = lgens('x, y', QQ, lex)
    p = lp.from_mon((1, 4, QQ(7, 3)))
    assert p == 7*y**4/3
    p1 = lp.from_mon(('y', 4, QQ(7, 3)))
    assert p1 == p

def test_variables():
    lp = LPoly('z, y, x', QQ, lex)
    z, y, x = lp.gens
    p = x*y + 3*y**2
    assert p.variables() == (1, 2)
    lp = LPoly('x, y, z', QQ, lex)
    x, y, z = lp.gens
    p = x*y + 3*y**2
    assert p.variables() == (0, 1)
    lp = LPoly('z,y,x,', QQ, lex)
    z, y, x = lp.gens
    p = x*y + 3*y**2
    assert p.variables() == (1, 2)

def test_mlgens():
    lp,x = mlgens('x', complex, lex)
    p = 1 + 1j + x*(1 - 1j)
    p2 = p*p
    assert p2 == -2j*x**2 + 4*x + 2j

def test_nclgens():
    lp, x = nclgens('x', sympify, lex)
    A = Operator('A')
    B = Operator('B')
    p1 = x*A + B
    assert p1**2 == x**2 * A**2 + B**2 + x*(A*B + B*A)
    p2 = p1.square_trunc('x', 2)
    assert p2 == B**2 + x*(A*B + B*A)
    p3 = p1.square_trunc(0, 2)
    assert p2 == p3
    p1 = x*A
    p2 = p1**2
    assert p2 == x**2 * A**2
    # Baker Campbell Hausdorff formula
    prec = 4
    pa = (x*A).exp('x', prec)
    pb = (x*B).exp('x', prec)
    p = pa.mul_trunc(pb, 'x', prec)
    p1 = p.log('x', prec)
    p1 = p1.expand()
    assert p1 == (-A*B*A/6 + A*B**2/12 + A**2*B/12 - B*A*B/6 + B*A**2/12 +
                  B**2*A/12)*x**3 + (A*B/2 - B*A/2)*x**2 + (A + B)*x

    p1 = 1 + x*A + x**2*B
    p2 = p1.series_inversion('x', 4)
    assert p2 == 1 - x*A - x**2*B + x**2*A**2 + x**3*(A*B + B*A - A**3)
    p1 = 2 + x*A + x**2*B
    p2 = p1.series_inversion('x', 4).expand()
    assert p2 == QQ(1, 2) - x*A/4 - x**2*B/4 + x**2*A**2/8 + x**3*(A*B/8 +
            B*A/8 - A**3/16)
    p1 = x*A + x**2*B
    p2 = p1.atan('x', 4)
    assert p2 == A*x + B*x**2 - A**3*x**3/3
    p2 = p2.tan('x', 4)
    assert p2 == p1
    p2 = p1.atanh('x', 4)
    assert p2 == A*x + B*x**2 + A**3*x**3/3
    p2 = p2.tanh('x', 4)
    assert p2 == p1

    lp1, x = lgens('x', sympify, lex, commuting=False)
    p1 = x*A + B
    assert p1**2 == x**2 * A**2 + B**2 + x*(A*B + B*A)
    lp2, y = lgens('y', lp1, lex)
    p1 = x*A + y*B
    assert p1**2 == A**2*x**2 + (A*B+B*A)*x*y + B**2*y**2

    p1 = x + A*x**2
    def test1(p):
        p2 = p1.series_reversion('x', 4, 'y')
    raises(NotImplementedError, 'test1(p1)')
    def test2(p):
        p2 = p.asinh('x', 4)
    raises(NotImplementedError, 'test2(p1)')
    def test3(p):
        p2 = p.asin('x', 4)
    raises(NotImplementedError, 'test3(p1)')
    def test4(p):
        p2 = p.nth_root(2, 'x', 4)
    raises(NotImplementedError, 'test4(p1)')
    def test5(p):
        p2 = p.series_inversion('x', 4)
    raises(NotImplementedError, 'test5(p1)')
    def test6(p):
        p2 = p.lambert('x', 4)
    raises(NotImplementedError, 'test6(p1)')
    def test7(p):
        p2 = p.log('x', 4)
    raises(NotImplementedError, 'test7(p1)')

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
    x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10 = lp.gens
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
    assert p1.coeff(m) == 15
    p2 = p1.coeff(x0**10)
    assert p2 == 0
    p2 = x0 + 1
    def test1(p1, p2):
        c = p1.coefficient(p2)
    raises(TypeError, 'test1(p1, p2)')

def test_coefficient_t():
    lp, x, y = lgens('x, y', QQ, lex)
    p = (1 + x + y)**3
    assert p.coefficient_t((0,0)) == y**3 + 3*y**2 + 3*y + 1
    assert p.coefficient_t((0,1)) == 3*y**2 + 6*y + 3
    assert p.coefficient_t((0,2)) == 3*y + 3
    assert p.coefficient_t((0,3)) == 1
    assert p.coefficient_t((1,0)) == x**3 + 3*x**2 + 3*x + 1

def test_add():
    lp, x, y = lgens('x, y', QQ, lex)
    zero = lp(0)
    one = lp(1)
    assert zero + zero == zero
    assert zero + one == one
    p1 = -2*x
    p2 = 2*x
    assert p1+p2 == zero
    assert p1+p1 == -4*x
    p = x
    p -= 1
    assert p == x - 1
    p = x
    p += 1
    assert p == x + 1
    p = 2*x
    p -= 1
    assert p == 2*x - 1
    p += 1
    assert p == 2*x
    p = 2*x + 1
    p -= 1
    assert p == 2*x

    p2 = 3*y**2/7
    p3 = p1 + p2
    assert p3 == -2*x+3*y**2/7
    p4 = p3
    p3 = p3 + 1
    assert p3 == -2*x+3*y**2/7 + 1
    assert p4 == -2*x+3*y**2/7
    p3 = -1 + p3
    assert p3 == -2*x+3*y**2/7
    lp1, x = lgens('x', QQ, lex)
    lp2, y = lgens('y', QQ, lex)
    def test1(x, y):
        p = x + y
    raises(ValueError, 'test1(x, y)')

def test_iadd():
    lp, x, y = lgens('x, y', QQ, lex)
    p3 = -2*x + 3*y**2/7
    # check that p3 is mutable
    p4 = p3
    p3 += 1
    assert p3 == -2*x + 3*y**2/7 + 1
    assert p4 == p3
    p3 = p3 - 1
    assert p3 == -2*x + 3*y**2/7
    p3 -= lp('1')
    assert p3 == -2*x + 3*y**2/7 - 1
    p3 -= 1
    assert p3 == -2*x + 3*y**2/7 - 2
    assert p3 - 3*y**2/7 == -2*x-2
    p4 = 2 - p3
    assert p4 == 2*x - 3*y**2/7 + 4
    p4 -= 1
    assert p4 == 2*x - 3*y**2/7 + 3
    p4 -= 2*x
    assert p4 == -3*y**2/7 + 3
    p4 += -3
    assert p4 == -3*y**2/7
    p1 = 1 + x + x*y + x**2
    p2 = x**2 + x**3
    p1 -= p2
    assert p1 == 1 + x + x*y - x**3
    p1 = x
    p1 -= p2
    assert p1 == x - x**2 - x**3
    assert str(x) == 'x'
    p1 = 1 + x + x*y + 2*x**2
    p1 -= p2
    assert p1 == 1 + x + x*y + x**2 - x**3
    lp2, z = lgens('z', QQ, lex)
    def test1(p):
        p -= z
    raises(ValueError, 'test1(p1)')
    def test2(p):
        p += z
    raises(ValueError, 'test2(p1)')
    def test3(p):
        p = p - z
    raises(ValueError, 'test3(p1)')
    def test4(p):
        p = p + z
    raises(ValueError, 'test4(p1)')

def test_mul():
    lp, x, y = lgens('x, y', QQ, lex)
    p1 = x*x
    assert p1 * 0 == lp(0)
    assert p1 * 1 == p1
    p2 = x + x**2
    assert p1*p2 == x**3 + x**4
    assert 2*p2 == 2*x + 2*x**2
    p3 = y + 2*y**2/3
    p4 = 2*x**2*y**2/3 + x**2*y + 2*x*y**2/3 + x*y
    assert p2*p3 == p4
    p = x
    p1 = p.imul_num(0)
    assert p1 == 0
    assert str(p) == 'x'
    p = 2*x
    p1 = p.imul_num(0)
    assert p == 0
    lp1, z = lgens('z', QQ, lex)
    def test1(p):
        p2 = p*z
    raises(ValueError, 'test1(p)')

def test_square():
    lp, x, y = lgens('x, y', QQ, lex)
    p1 = x.square()
    assert p1 == x**2
    p2 = (x + y).square()
    assert p2 == x*x + 2*x*y + y*y
    p3 = 3*x/7 + x**3*y + 8*y**2
    assert p3.square() == x**6*y**2 + 16*x**3*y**3 + 6*x**4*y/7 + 64*y**4 + 48*x*y**2/7 + 9*x**2/49

def test_mul_iadd():
    lp, x, y = lgens('x, y', QQ, lex)
    p = x**2 + y**2 - 1
    p1 = x + y
    p2 = x - y
    p.mul_iadd(p1, p2)
    assert p == 2*x*x - 1
    p = x
    p = p.mul_iadd(p1, p2)
    assert p == x + x**2 - y**2
    def test1(p):
        p2 = p.mul_iadd(p1, 2)
    raises(NotImplementedError, 'test1(p)')
    lp1, z = lgens('z', QQ, lex)
    def test2(p):
        p2 = p.mul_iadd(p1, z)
    raises(ValueError, 'test2(p)')

def test_iadd_mon():
    lp, x, y = lgens('x, y', QQ, lex)
    p1 = x*x + y*y - 2*x*y**2
    m = monomial_from_sequence((1, 2))
    p1 = p1.iadd_mon((m, QQ(2)))
    assert p1 == x**2 + y**2

def test_iadd_m_mul_q():
    lp, x, y = lgens('x, y', QQ, lex)
    p1 = x**2 + y**2
    p2 = x**3 + y**3
    m = x*y
    expva = m.keys()[0]
    c = m[expva]
    p1 = p1.iadd_m_mul_q(p2, (expva, c))
    assert p1 == x**2 + y**2 + x**4*y + x*y**4
    p1 = x
    p1 = p1.iadd_m_mul_q(p2, (expva, c))
    assert p1 == x + x**4*y + x*y**4
    assert str(x) == 'x'

def test_leading_expv():
    lp, x, y = lgens('x, y', QQ, lex)
    p = 2*x**2 + 3*y**3
    assert p.leading_expv() == monomial_from_sequence((2, 0))
    lp, y, x = lgens('y, x', QQ, lex)
    p = 2*x**2 + 3*y**3 + x**2*y**2
    assert p.leading_expv() == monomial_from_sequence((3, 0))
    lp, y, x = lgens('y, x', QQ, grlex)
    p = 2*x**2 + 3*y**3 + x**2*y**2
    assert p.leading_expv() == monomial_from_sequence((2, 2))
    p = lp(0)
    assert p.leading_expv() == None

def test_leading_term():
    lp, x, y = lgens('x, y', QQ, lex)
    p = 2*x**2 + 3*y**3
    assert p.leading_term() == 2*x**2

def test_expand():
    lp, x = lgens('x', sympify, lex)
    p1 = x + cos(2)*(1 + cos(3))
    p2 = x + cos(2) + cos(2)*cos(3)
    assert p1.expand() == p2
    lp, x = lgens('x', QQ, lex)
    p = x + 1
    assert p.expand() == p

def test_div():
    lp, x, y = lgens('x, y', QQ, lex)
    p = y**2/4
    assert 4*p == y**2
    p1 = x**2 + y**2 - 1
    p2 = p1/7
    assert p2 == x**2/7 + y**2/7 - QQ(1, 7)
    p1 = x**2 + x
    assert p1/x == x + 1
    p1 = x*y - x + y - 1
    assert p1/(y - 1) == x + 1
    def test1(p):
        p2 = p/0
    raises(ZeroDivisionError, 'test1(p1)')
    def test2(p):
        p2 =  p/y
    raises(NotImplementedError, 'test2(p1)')
    lp1, z = lgens('z', QQ, lex)
    def test3(p):
        p2 = p/z
    raises(NotImplementedError, 'test3(p1)')

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
    assert p**0 == 1
    assert p**1 == p
    p = lp(0)
    def test1(p):
        p1 = p**0
    raises(ValueError, 'test1(p)')
    p = 1 + x
    def test2(p):
        p1 = p**-1
    raises(ValueError, 'test2(p)')
    lp, x = lgens('x', QQ, lex)
    p = 1 + x + 3*x**2
    p1 = p**25
    assert p1.coeff(x**50) == 847288609443
    def test1(p):
        p1 = p**QQ(1,3)
    raises(ValueError, 'test1(p)')
    def test2():
        p1 = x**QQ(1,3)
    raises(ValueError, 'test2()')

def test_division():
    lp, x, y = lgens('x, y', QQ, lex)
    f = x**3
    f0 = x-y**2
    f1 = x-y
    qv, r = f.division((f0, f1))
    assert (qv[0], qv[1], r) == (x**2 + x*y**2 + y**4, lp(0), y**6)
    expv0 = f0.leading_expv()
    expv1 = f1.leading_expv()
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
    f = lp(0)
    qv, r = f.division((g, ))
    assert not qv and r == 0
    f = x + 1
    g = lp(1)
    qv, r = f.division((g, ))
    assert qv[0] == f and r == 0
    lp1, z = lgens('z', QQ, lex)
    g = z
    def test1(f, g):
        q, r = f.division([g])
    raises(ValueError, 'test1(f, g)')

def test_mul_trunc():
    lp, x, y, t = lgens('x, y, t', QQ, lex)
    p = 1 + t*x + t*y
    for i in range(2):
        p = p.mul_trunc(p, 't', 3)

    assert p == 6*x**2*t**2 + 12*x*y*t**2 + 6*y**2*t**2 + 4*x*t + 4*y*t + 1
    p = 1 + t*x + t*y + t**2*x*y
    p1 = p.mul_trunc(p, 2, 2)
    assert p1 == 1 + 2*t*x + 2*t*y
    lp1, z = lgens('z', QQ, lex)
    def test1(p):
        p2 = p.mul_trunc(z, 'x', 2)
    raises(ValueError, 'test1(p)')

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
    lp, x = lgens('x', QQ, lex)
    h = 7
    p = x.cos('x', h)
    p1 = p.pow_trunc(Rational(3, 2), 'x', h)
    assert p1 == 1 - 3*x**2/4 + 5*x**4/32 - 19*x**6/1920
    p1 = p.pow_trunc(Rational(-3, 2), 'x', h)
    assert p1 == 1 + 3*x**2/4 + 13*x**4/32 + 379*x**6/1920
    p1 = p.pow_trunc(Rational(1, 2), 'x', h)
    assert p1 == 1 - x**2/4 - x**4/96 - 19*x**6/5760
    p1 = p.pow_trunc(Rational(-1, 2), 'x', h)
    assert p1 == 1 + x**2/4 + 7*x**4/96 + 139*x**6/5760
    p = 1 + x
    p1 = p.pow_trunc(Rational(3,1), 'x', 2)
    assert p1 == 1 + 3*x
    assert p.pow_trunc(0, 'x', 2) == 1
    assert p.pow_trunc(-2, 'x', 2) == 1 - 2*x
    p = lp(0)
    def test1(p):
        p1 = p.pow_trunc(0, 'x', 4)
    raises(ValueError, 'test1(p)')

def test_pow_miller():
    lp, x = lgens('x', QQ, lex)
    n = 5
    for c in [0, 1]:
        p1 = (c + 3*x + 2*x**2)**n
        p2 = (c + 3*x + 2*x**2).pow_miller(n)
        assert p1 == p2

def test_has_constant_term():
    lp, x, y, z = lgens('x, y, z', QQ, lex)
    p = y + x*z
    assert p.has_constant_term('x')
    p = x + x**4
    assert not p.has_constant_term('x')
    p = 1 + x + x**4
    assert p.has_constant_term('x')
    p = x + y + x*z
    def test1(p):
        b = p.has_constant_term(['y', 'z'])
    raises(NotImplementedError, 'test1(p)')

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
    p = x + x**2
    def test1(p):
        p1 = p.series_inversion('x', 4)
    raises(NotImplementedError, 'test1(p)')
    lp, x, y = lgens('x, y', QQ, lex)
    p = 1 + x + y
    def test2(p):
        p1 = p.series_inversion('x', 4)
    raises(NotImplementedError, 'test2(p)')
    def test3(p):
        p1 = p._series_inversion1('x', 4)
    raises(ValueError, 'test3(p)')
    p = x + y
    def test4(p):
        p1 = p._series_inversion1('x', 4)
    raises(ValueError, 'test4(p)')

def test_reversion():
    lp, x, y = lgens('x, y', QQ, lex)
    p = x+x**2
    prec = 10
    p1 = p.series_reversion('x', prec, 'y')
    assert p1 == 1430*y**9 - 429*y**8 + 132*y**7 - 42*y**6 + 14*y**5 - 5*y**4 + 2*y**3 - y**2 + y
    assert p1 + p1.pow_trunc(2, 'y', prec) == y
    p = x.asin('x', prec)
    p1 = p.series_reversion('x', prec, 'y')
    assert p1 == y.sin('y', prec)
    prec = 3
    p = x.asin('x', prec)
    p1 = p.series_reversion(0, prec, 1)
    assert p1 == y.sin('y', prec)
    p = x + 1
    def test1(p):
        p1 = p.series_reversion('x', prec, 'y')
    raises(ValueError, 'test1(p)')

    prec = 10
    lp, x, y, z = lgens('x, y, z', QQ, lex)
    p = x + x**2 + x**2*y + x**3*y
    r = p.series_reversion('x', prec, 'z')
    assert p.subs_trunc('z', prec, x=r) == z

def test_series_from_list():
    lp, x, y = lgens('x, y', QQ, lex)
    c = [1, 3, 5, 7]
    p1 = (x+y).series_from_list(c, 'x', 3, concur=0)
    p2 = (1 + 3*(x+y) + 5*(x+y)**2 + 7*(x+y)**3).trunc('x', 3)
    assert p1 == p2
    lp, x = lgens('x', QQ, lex)
    h = 25
    p = x.exp('x', h) - 1
    p1 = p.series_from_list(c, 'x', h)
    p2 = lp(0)
    for i, cx in enumerate(c):
        p2 += cx*p.pow_trunc(i, 'x', h)
    assert p1 == p2
    h = 10
    p = x.exp('x', h)
    p1 = p.fun('_nth_root1', 2, 'x', h)
    p2 = (x/2).exp('x', h)
    assert p1 == p2
    p = x + x**2
    p1 = p.fun('_exp1', 'x', h)
    p2 = p.exp('x', h)
    assert p1 == p2
    h = 4
    def square(x, iv, h):
        return x*x
    p1 = p.fun(square, 'x', h)
    assert p1 == x**2 + 2*x**3

def test_derivative():
    gens=['x%d' % i for i in range(11)]
    lp = LPoly(gens, QQ, lex)
    p = lp('288/5*x0**8*x1**6*x4**3*x10**2 + 8*x0**2*x2**3*x4**3 +2*x0**2-2*x1**2')
    p1 = p.derivative(0)
    assert p1 == lp('2304/5*x0**7*x1**6*x4**3*x10**2 + 16*x0*x2**3*x4**3 + 4*x0')
    p1 = p.derivative('x0')
    assert p1 == lp('2304/5*x0**7*x1**6*x4**3*x10**2 + 16*x0*x2**3*x4**3 + 4*x0')

def test_integrate():
    lp, x, y = lgens('x, y', QQ, lex)
    p = x**2 * y**5 + x + y
    p1 = p.integrate('x')
    assert p1 == x**3*y**5/3 + x**2/2 + x*y
    p1 = p.integrate('y')
    assert p1 == lp('1/6*x**2*y**6 + x*y + 1/2*y**2')
    assert p1 == x**2*y**6/6 + x*y + y**2/2

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
    p = lp(1)
    p1 = p.nth_root(4, 'x', h)
    assert p1 == p
    p = 1 + x
    p1 = p.nth_root(0, 'x', h)
    assert p1 == lp(1)
    p1 = p.nth_root(1, 'x', h)
    assert p1 == p1

    p = 4 + x**2 + 2*x**3 + x**4/3
    p1 = p.nth_root(2, 'x', 5)
    assert p1 == 13*x**4/192 + x**3/2 + x**2/4 + 2

    lp, x = lgens('x', sympify, lex)
    p = 2 + x
    p1 = p.nth_root(2, 'x', 4)
    c = sqrt(2)
    assert p1 == c + c*x/4 - c*x**2/32 + c*x**3/128
    p = -2 + x
    def test1(p):
        p1 = p.nth_root(2, 'x', 4)
    raises(NotImplementedError, 'test1(p)')
    p = x
    def test2(p):
        p1 = p.nth_root(2, 'x', 4)
    raises(NotImplementedError, 'test2(p)')
    lp, x = lgens('x', QQ, lex)
    p = 2 + x
    def test3(p):
        p1 = p.nth_root(2, 'x', 4)
    raises(TaylorEvalError, 'test3(p)')
    p = lp(0)
    def test4(p):
        p1 = p.nth_root(0, 'x', 4)
    raises(ValueError, 'test4(p)')
    n = QQ(2, 3)
    p = 1 + x
    def test5(p, n):
        p1 = p.nth_root(n, 'x', 4)
    raises(ValueError, 'test5(p, n)')

def test_sqrt():
    lp, x, y = lgens('x, y', QQ, lex)
    h = 4
    p0 = (1 + x**2 * y**5 + x)**3
    p0 = p0.trunc('x', h)
    p1 = p0.square_trunc('x', h)
    p2 = p1.sqrt('x', h)
    assert p2 == p0
    p3 = p0.square_trunc(0, h)
    assert p1 == p3

def test_log():
    # log of univariate series
    lp, x = lgens('x', QQ, lex)
    p = 1 + x
    p1 = p.log('x', 4)
    assert p1 == x - x**2/2 + x**3/3
    p = 1 + x +2*x**2/3
    p1 = p.log('x', 9)
    assert p1 == -17*x**8/648 + 13*x**7/189 - 11*x**6/162 - x**5/45 + \
      7*x**4/36 - x**3/3 + x**2/6 + x

    p2 = p.series_inversion('x', 9)
    p3 = p2.log('x', 9)
    assert p3 == -p1
    lp, x, y = lgens('x, y', QQ, lex)
    p = 1 + x + 2*y*x**2
    p1 = p.log('x', 6)
    assert p1 == 4*x**5*y**2 - 2*x**5*y - 2*x**4*y**2 + x**5/5 + 2*x**4*y - x**4/4 - 2*x**3*y + x**3/3 + 2*x**2*y - x**2/2 + x

    lp, x = lgens('x', sympify, lex)
    p = 2 + x
    p1 = p.log('x', 4)
    assert p1 == log(2) + x/2 - x**2/8 + x**3/24
    p = -2 + x
    def test1(p):
        p1 = p.log('x', 4)
    raises(NotImplementedError, 'test1(p)')
    lp, x = lgens('x', QQ, lex)
    p = 2 + x
    def test2(p):
        p1 = p.log('x', 4)
    raises(TaylorEvalError, 'test2(p)')

def test_acot():
    lp, x = lgens('x', QQ, lex)
    p = x.re_acoth('x', 8)
    assert p == x**7/7 + x**5/5 + x**3/3 + x
    p = x.acot1('x', 8)
    assert p == x**7/7 - x**5/5 + x**3/3 - x
    p = x + 1
    def test1(p):
        p1 = p.acot1('x', 4)
    raises(NotImplementedError, 'test1(p)')
    p = x + 1
    def test1(p):
        p1 = p.re_acoth('x', 4)
    raises(NotImplementedError, 'test1(p)')

def test_atan():
    lp, x, y = lgens('x, y', QQ, lex)
    p = x*y + x**2
    p1 = p.atan('x', 6)
    assert p1 == x**5*y**5/5 -x**5*y -x**4*y**2 -x**3*y**3/3 +x**2 +x*y
    p = x + 1
    def test1(p):
        p1 = p.atan('x', 4)
    raises(NotImplementedError, 'test1(p)')

def test_atanh():
    lp, x, y = lgens('x, y', QQ, lex)
    p = x*y + x**2
    p1 = p.atanh('x', 6)
    assert p1 == x**5*y**5/5 + x**5*y + x**4*y**2 + x**3*y**3/3 + x**2 + x*y
    p = x + 1
    def test1(p):
        p1 = p.atanh('x', 4)
    raises(NotImplementedError, 'test1(p)')

def test_tan():
    lp, x = lgens('x', QQ, lex)
    p = x + x**4
    p1 = p.tan('x', 8)
    assert p1 == x + x**3/3 + x**4 + 2*x**5/15 + x**6 + 17*x**7/315
    p = 1 + x
    def test1(p):
        p1 = p.tan('x', 4)
    raises(NotImplementedError, 'test1(p)')

def test_tanh():
    lp, x, y = lgens('x, y', QQ, lex)
    p = x*y
    p1 = p.tanh('x', 8)
    assert p1 == -17*x**7*y**7/315 + 2*x**5*y**5/15 - x**3*y**3/3 + x*y
    p = x*y + x
    p = x + 1
    def test1(p):
        p1 = p.tanh('x', 4)
    raises(NotImplementedError, 'test1(p)')

def test_exp():
    lp, x = lgens('x', QQ, lex)
    p = x + x**4
    for h in [10, 30]:
        q = (1 + p).series_inversion('x', h) - 1
        p1 = q.exp('x', h)
        q1 = p1.log('x', h)
        assert q1 == q
    p1 = p.exp('x', 30)
    assert p1.coeff(x**29) == QQ(74274246775059676726972369, 353670479749588078181744640000)
    def test1():
        p = (x + 1).exp('x', 3)
    raises(TaylorEvalError, 'test1()')
    prec = 21
    p = (1 + x).log('x', prec)
    p1 = p.exp('x', prec)
    assert p1 == x + 1

def test_sin():
    lp, x, y = lgens('x, y', QQ, lex)
    h = 10
    p = 1 + x*y + x**4
    p = p.series_inversion('x', h) - 1
    p1 = p.sin('x', h)
    assert p1 == lp('224179/72576*x**9*y**9 - 1591/720*x**8*y**8 - 15/4*x**9*y**5 + 6931/5040*x**7*y**7 - 1/24*x**8*y**4 - 5/8*x**6*y**6 - 5/2*x**9*y + 2*x**7*y**3 - 1/120*x**5*y**5 + x**8 - 5/2*x**6*y**2 + 1/2*x**4*y**4 + 2*x**5*y - 5/6*x**3*y**3 - x**4 + x**2*y**2 - x*y')
    lp, x = lgens('x', sympify, lex)
    p = 1 + x
    p1 = p.sin('x', 3)
    c1, s1 = cos(1), sin(1)
    assert p1 == s1 + c1*x - s1*x**2/2
    prec = 42
    p = x.asin('x', prec)
    p1 = p.sin('x', prec)
    assert p1 == x
    p = I + x
    def test1(p):
        p1 = p.sin('x', 3)
    raises(TaylorEvalError, 'test1(p)')

    lp, x = lgens('x', QQ, lex)
    p = 1 + x
    def test2(p):
        p1 = p.sin('x', 3)
    raises(TaylorEvalError, 'test2(p)')

def test_cos():
    lp, x, y = lgens('x, y', QQ, lex)
    h = 10
    p = 1 + x*y + x**4
    p = p.series_inversion('x', h) - 1
    p1 = p.cos('x', h)
    assert p1 == lp('8791/5040*x**9*y**9 - 16699/8064*x**8*y**8 - 1501/120*x**9*y**5 + 87/40*x**7*y**7 + 55/6*x**8*y**4 - 1501/720*x**6*y**6 + 3*x**9*y - 35/6*x**7*y**3 + 11/6*x**5*y**5 - 1/2*x**8 + 3*x**6*y**2 - 35/24*x**4*y**4 - x**5*y + x**3*y**3 - 1/2*x**2*y**2 + 1')

    lp, x = lgens('x', sympify, lex)
    p = 1 + x
    p1 = p.cos('x', 3)
    c1, s1 = cos(1), sin(1)
    assert p1 == c1 - x*s1 - x**2*c1/2
    p = I + x
    def test1(p):
        p1 = p.cos('x', 3)
    raises(NotImplementedError, 'test1(p)')
    lp, x = lgens('x', QQ, lex)
    prec = 42
    p = x.asin('x', prec)
    p1 = p.cos('x', prec)
    assert p1.coeff(x**40) == QQ(-883631595, 274877906944)
    p = 1 + x
    def test2(p):
        p1 = p.cos('x', 3)
    raises(TaylorEvalError, 'test2(p)')

def test_cos_sin():
    lp, x, y = lgens('x, y', QQ, lex)
    h = 10
    p = 1 + x*y + x**4
    p = p.series_inversion('x', h) - 1
    assert p.cos_sin('x', h) == (p.cos('x', h), p.sin('x', h))
    assert p.cosh_sinh('x', h) == (p.cosh('x', h), p.sinh('x', h))

def test_sinh():
    lp, x, y = lgens('x, y', QQ, lex)
    p = x*y
    p1 = p.sinh('x', 8)
    assert p1 == x**7*y**7/5040 + x**5*y**5/120 + x**3*y**3/6 + x*y

def test_cosh():
    lp, x, y = lgens('x, y', QQ, lex)
    p = x*y
    p1 = p.cosh('x', 8)
    assert p1 == x**6*y**6/720 + x**4*y**4/24 + x**2*y**2/2 + 1

def test_lambert():
    lp, x = lgens('x', QQ, lex)
    p1 = x.lambert('x', 10)
    assert p1 == lp('531441/4480*x**9 - 16384/315*x**8 + 16807/720*x**7 - 54/5*x**6 + 125/24*x**5 - 8/3*x**4 + 3/2*x**3 - x**2 + x')
    p = x + 1
    def test1(p):
        p1 = p.lambert('x', 4)
    raises(NotImplementedError, 'test1(p)')

def test_asin():
    lp, x = lgens('x', QQ, lex)
    p1 = x.asin('x', 12)
    assert p1 == 63*x**11/2816 + 35*x**9/1152 + 5*x**7/112 + 3*x**5/40 + x**3/6 + x
    p = x + 1
    def test1(p):
        p1 = p.asin('x', 4)
    raises(NotImplementedError, 'test1(p)')

#slow
def test_asin1():
    lp, x = lgens('x', QQ, lex)
    prec = 42
    p = x.sin('x', prec)
    p1 = p.asin('x', prec)
    assert p1 == x

def test_asinh():
    lp, x = lgens('x', QQ, lex)
    p1 = x.asinh('x', 12)
    assert p1 == -63*x**11/2816 + 35*x**9/1152 - 5*x**7/112 + 3*x**5/40 - x**3/6 + x
    p = x + 1
    def test1(p):
        p1 = p.asinh('x', 4)
    raises(NotImplementedError, 'test1(p)')

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
    lp, y = lgens('y', QQ, lex)
    p = y/3 + 2*y**2/4
    lp1, x = lgens('x', sympify, lex)
    p1 = p.toSR(lp1)
    assert p1 == x/3 + 2*x**2/4
    assert p1.coeff(x).is_Rational
    p2 = p1.toSR(lp1)
    assert p1 == p2

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
    assert p1.coefficient(y**20) == lp(1)
    p = x + x**3
    p1 = p.subs(x=(1+y)**2)
    assert p1 == y**6 + 6*y**5 + 15*y**4 + 20*y**3 + 16*y**2 + 8*y + 2

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
    lp = LPoly('x:11', QQ, lex)
    x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10 = lp.gens
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
    lp1, z = lgens('z', sympify, lex)
    def test1(lp, lp1, z):
        sb = LPolySubs(lp, lp1, {'x':z})
    raises(NotImplementedError, 'test1(lp, lp1, z)')

def test_polypoly():
    lp1, x = lgens('x', QQ, lex)
    lp2, y = lgens('y', lp1, lex)
    assert lp1.__class__ == LPoly
    assert lp2.__class__ == MLPoly
    assert str(x) == 'x'
    assert str(y) == '(1)*y'
    assert str(y**2) == '(1)*y**2'
    assert str(-y**2) == '(-1)*y**2'
    p1 = (x+1)*y + 1
    p2 = x*y + y + 1
    assert p1 == p2
    p3 = x
    p2 = p2 - p3
    assert p2 == x*y + y + 1 - x

    one = lp2(1)
    assert one == lp2(1)
    assert y*one == y
    zero = lp2(0)
    assert not zero*y
    assert zero == lp2(0)
    assert lp2(y) == y
    p2 = 2*x + 2*y
    p1 = lp1(2)
    p2 = p2/p1
    assert p2 == x + y
    def test1(lp):
        p = lp('+( +x^2)*y^4')
    raises(NotImplementedError, 'test1(lp2)')

def test_polypoly_add():
    lp1, x = lgens('x', QQ, lex)
    lp2, y = lgens('y', lp1, lex)
    p3 = 1 + y
    assert p3 - 1 == y
    p = x + y
    assert p == y + x
    assert p - x == y

    p1 = (1+x)*y
    p2 = p1 + lp2(0)
    assert p1 == p2
    p2 += y
    assert p2 == x*y + 2*y
    p3 = p2 - x
    assert p3 == x*y + 2*y - x
    p3 = p3 + x
    assert p3 == x*y + 2*y

def test_polypoly_sub():
    lp1, x = lgens('x', QQ, lex)
    lp2, y = lgens('y', lp1, lex)
    p = x*y**2
    p2 = p*p
    p3 = x**2*y**4
    assert p2 == p3
    p3 = y - 1
    assert p3 + 1 == y

def test_polypoly_mul():
    lp1, x = lgens('x', QQ, lex)
    lp2, y = lgens('y', lp1, lex)
    p = x*y**2
    p2 = p*p
    p2a = x**2*y**4
    assert p2 == p2a
    p3 = x**2 * y**4
    assert p2 == p3
    assert p2a == x**2*y**4

    p = (x+1)*y**2 - x*y + 1
    p2 = p*p
    p3 = (x**2 + 2*x + 1)*y**4 + (-2*x**2 - 2*x)*y**3 + (x**2 + 2*x + 2)*y**2 - 2*x*y + 1
    p4 = p2 - p3
    assert p4 == lp2(0)
    assert p2 == p3
    p = x + y
    p1 = 1 + x + y
    p2 = y
    p = p.mul_iadd(p1,p2)
    assert p == x + 2*y + x*y + y**2

def test_polypoly_div():
    lp1, x = lgens('x', QQ, lex)
    lp2, y = lgens('y', lp1, lex)
    p1 = x*y**2
    assert p1/lp1(2) == x*y**2/2
    def test1(p):
        p2 = p/x
    raises(NotImplementedError, 'test1(p1)')

def test_polypoly_pow():
    lp1, x = lgens('x', QQ, lex)
    lp2, y = lgens('y', lp1, lex)
    p = x + y
    p2 = p**2
    p4 = p**4
    assert p2 == x**2 + 2*x*y + y**2
    assert p**3 == x**3 + 3*x**2*y + 3*x*y**2 + y**3
    assert p4 == x**4 + 4*x**3*y + 6*x**2*y**2 + 4*x*y**3 + y**4

    # poly on poly on poly example
    lp3, z = lgens('z', lp2, lex)
    p = (1 + x)*y + (1 + y)*z
    p1 = p**3
    assert p1 == (y**3 + 3*y**2 + 3*y + 1)*z**3 + ((3*x + 3)*y**3 + (6*x + 6)*y**2 + (3*x + 3)*y)*z**2 + ((3*x**2 + 6*x + 3)*y**3 + (3*x**2 + 6*x + 3)*y**2)*z + (x**3 + 3*x**2 + 3*x + 1)*y**3

def test_polypoly_pow_trunc():
    lp1, x = lgens('x', QQ, lex)
    lp2, y = lgens('y', lp1, lex)
    p = 1 + x + x*y
    p1 = p**4
    p1 = p1.trunc('y',3)
    assert p1 == 6*x**4*y**2 + 4*x**4*y + x**4 + 12*x**3*y**2 + 12*x**3*y + 4*x**3 + 6*x**2*y**2 + 12*x**2*y + 6*x**2 + 4*x*y + 4*x + 1
