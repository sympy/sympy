from sympy.core.symbol import Symbol
from sympy.core.numbers import Rational
from sympy.utilities.pytest import raises
_pyround = round
from sympy.functions.elementary.miscellaneous import sqrt, root, Min, Max, real_root, round
from sympy import S, Float, I, cos, sin, oo, pi, Add

def test_Min():
    from sympy.abc import x, y, z
    n = Symbol('n', negative=True)
    n_ = Symbol('n_', negative=True)
    nn = Symbol('nn', nonnegative=True)
    nn_ = Symbol('nn_', nonnegative=True)
    p = Symbol('p', positive=True)
    p_ = Symbol('p_', positive=True)
    np = Symbol('np', nonpositive=True)
    np_ = Symbol('np_', nonpositive=True)

    assert Min(5, 4) == 4
    assert Min(-oo, -oo) == -oo
    assert Min(-oo, n) == -oo
    assert Min(n, -oo) == -oo
    assert Min(-oo, np) == -oo
    assert Min(np, -oo) == -oo
    assert Min(-oo, 0) == -oo
    assert Min(0, -oo) == -oo
    assert Min(-oo, nn) == -oo
    assert Min(nn, -oo) == -oo
    assert Min(-oo, p) == -oo
    assert Min(p, -oo) == -oo
    assert Min(-oo, oo) == -oo
    assert Min(oo, -oo) == -oo
    assert Min(n, n) == n
    assert Min(n, np) == Min(n, np)
    assert Min(np, n) == Min(np, n)
    assert Min(n, 0) == n
    assert Min(0, n) == n
    assert Min(n, nn) == n
    assert Min(nn, n) == n
    assert Min(n, p) == n
    assert Min(p, n) == n
    assert Min(n, oo) == n
    assert Min(oo, n) == n
    assert Min(np, np) == np
    assert Min(np, 0) == np
    assert Min(0, np) == np
    assert Min(np, nn) == np
    assert Min(nn, np) == np
    assert Min(np, p) == np
    assert Min(p, np) == np
    assert Min(np, oo) == np
    assert Min(oo, np) == np
    assert Min(0, 0) == 0
    assert Min(0, nn) == 0
    assert Min(nn, 0) == 0
    assert Min(0, p) == 0
    assert Min(p, 0) == 0
    assert Min(0, oo) == 0
    assert Min(oo, 0) == 0
    assert Min(nn, nn) == nn
    assert Min(nn, p) == Min(nn, p)
    assert Min(p, nn) == Min(p, nn)
    assert Min(nn, oo) == nn
    assert Min(oo, nn) == nn
    assert Min(p, p) == p
    assert Min(p, oo) == p
    assert Min(oo, p) == p
    assert Min(oo, oo) == oo

    assert Min(n, n_).func is Min
    assert Min(nn, nn_).func is Min
    assert Min(np, np_).func is Min
    assert Min(p, p_).func is Min

    # lists
    raises(ValueError, 'Min()')
    assert Min(x, y) == Min(y, x)
    assert Min(x, y, z) == Min(z, y, x)
    assert Min(x, Min(y, z)) == Min(z, y, x)
    assert Min(x, Max(y, -oo)) == Min(x, y)
    assert Min(p, oo, n,  p, p, p_) == n
    assert Min(p_, n_, p) == n_
    assert Min(n, oo, -7, p,  p, 2) == Min(n, -7)
    assert Min(2, x, p, n, oo, n_,  p, 2, -2, -2) == Min(-2, x, n, n_)
    assert Min(0, x, 1, y) == Min(0, x, y)
    assert Min(1000, 100, -100, x, p, n) == Min(n, x, -100)
    assert Min(cos(x), sin(x)) == Min(cos(x), sin(x))
    assert Min(cos(x), sin(x)).subs(x, 1) == cos(1)
    assert Min(cos(x), sin(x)).subs(x, S(1)/2) == sin(S(1)/2)
    raises(ValueError, 'Min(cos(x), sin(x)).subs(x, I)')
    raises(ValueError, 'Min(I)')
    raises(ValueError, 'Min(I, x)')
    raises(ValueError, 'Min(S.ComplexInfinity, x)')


def test_Max():
    from sympy.abc import x, y, z
    n = Symbol('n', negative=True)
    n_ = Symbol('n_', negative=True)
    nn = Symbol('nn', nonnegative=True)
    nn_ = Symbol('nn_', nonnegative=True)
    p = Symbol('p', positive=True)
    p_ = Symbol('p_', positive=True)
    np = Symbol('np', nonpositive=True)
    np_ = Symbol('np_', nonpositive=True)

    assert Max(5, 4) == 5

    # lists

    raises(ValueError, 'Max()')
    assert Max(x, y) == Max(y, x)
    assert Max(x, y, z) == Max(z, y, x)
    assert Max(x, Max(y, z)) == Max(z, y, x)
    assert Max(x, Min(y, oo)) == Max(x, y)
    assert Max(n, -oo, n_,  p, 2) == Max(p, 2)
    assert Max(n, -oo, n_,  p) == p
    assert Max(2, x, p, n, -oo, S.NegativeInfinity, n_,  p, 2) == Max(2, x, p)
    assert Max(0, x, 1, y) == Max(1, x, y)
    assert Max(x, x + 1, x - 1) == 1 + x
    assert Max(1000, 100, -100, x, p, n) == Max(p, x, 1000)
    assert Max(cos(x), sin(x)) == Max(sin(x), cos(x))
    assert Max(cos(x), sin(x)).subs(x, 1) == sin(1)
    assert Max(cos(x), sin(x)).subs(x, S(1)/2) == cos(S(1)/2)
    raises(ValueError, 'Max(cos(x), sin(x)).subs(x, I)')
    raises(ValueError, 'Max(I)')
    raises(ValueError, 'Max(I, x)')
    raises(ValueError, 'Max(S.ComplexInfinity, 1)')
    # interesting:
    # Max(n, -oo, n_,  p, 2) == Max(p, 2)
    # True
    # Max(n, -oo, n_,  p, 1000) == Max(p, 1000)
    # False


def test_root():
    from sympy.abc import x, y, z
    n = Symbol('n', integer=True)

    assert root(2, 2) == sqrt(2)
    assert root(2, 1) == 2
    assert root(2, 3) == 2**Rational(1,3)
    assert root(2, -5) == 2**Rational(4,5)/2

    assert root(-2, 1) == -2

    assert root(-2, 2) == sqrt(2)*I
    assert root(-2, 1) == -2

    assert root(x, 2) == sqrt(x)
    assert root(x, 1) == x
    assert root(x, 3) == x**Rational(1,3)
    assert root(x, -5) == x**Rational(-1,5)

    assert root(x, n) == x**(1/n)
    assert root(x, -n) == x**(-1/n)

def test_nthroot():
    assert real_root(-8, 3) == -2
    assert real_root(-16, 4) == root(-16, 4)
    r = root(-7, 4)
    assert real_root(r) == r
    r1 = root(-1, 3)
    r2 = r1**2
    r3 = root(-1, 4)
    assert real_root(r1 + r2 + r3) == -1 + r2 + r3

def test_round():
    from sympy.abc import x

    assert round(Float('0.1249999'), 2) == 0.12
    d20 = 12345678901234567890
    ans = round(S(d20), 2)
    assert ans.is_Float and ans == d20
    ans = round(S(d20), -2)
    assert ans.is_Float and ans == 12345678901234567900
    assert round(S('1/7'), 4) == 0.1429
    assert round(S('.[12345]'), 4) == 0.1235
    assert round(S('.1349'), 2) == 0.13
    n = S(12345)
    ans = round(n)
    assert ans.is_Float
    assert ans == n
    ans = round(n, 1)
    assert ans.is_Float
    assert ans == n
    ans = round(n, 4)
    assert ans.is_Float
    assert ans == n
    assert round(n, -1) == 12350

    r = round(n, -4)
    assert r == 10000
    # in fact, it should equal many values since __eq__
    # compares at equal precision
    assert all(r == i for i in range(9984, 10049))

    assert round(n, -5) == 0

    assert round(pi + sqrt(2), 2) == 4.56
    assert round(10*(pi + sqrt(2)), -1) == 50
    raises(TypeError, 'round(x + 2, 2)')
    assert round(S(2.3), 1) == 2.3
    e = round(12.345, 2)
    assert e == _pyround(12.345, 2)
    assert type(e) is float

    assert round(Float(.3, 3) + 2*pi) == 7
    assert round(Float(.3, 3) + 2*pi*100) == 629
    assert round(Float(.03, 3) + 2*pi/100, 5) == 0.09283
    assert round(Float(.03, 3) + 2*pi/100, 4) == 0.0928

    assert round(S.Zero) == 0

    a = (Add(1, Float('1.'+'9'*27, ''), evaluate=0))
    assert round(a, 10) == Float('3.0000000000','')
    assert round(a, 25) == Float('3.0000000000000000000000000','')
    assert round(a, 26) == Float('3.00000000000000000000000000','')
    assert round(a, 27) == Float('2.999999999999999999999999999','')
    assert round(a, 30) == Float('2.999999999999999999999999999','')

    raises(TypeError, 'round(x)')
    raises(TypeError, 'round(1 + 3*I)')
