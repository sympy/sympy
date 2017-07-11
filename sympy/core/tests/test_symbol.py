from sympy import (Symbol, Wild, GreaterThan, LessThan, StrictGreaterThan,
    StrictLessThan, pi, I, Rational, sympify, symbols, Dummy
)

from sympy.utilities.pytest import raises


def test_Symbol():
    a = Symbol("a")
    x1 = Symbol("x")
    x2 = Symbol("x")
    xdummy1 = Dummy("x")
    xdummy2 = Dummy("x")

    assert a != x1
    assert a != x2
    assert x1 == x2
    assert x1 != xdummy1
    assert xdummy1 != xdummy2

    assert Symbol("x") == Symbol("x")
    assert Dummy("x") != Dummy("x")
    d = symbols('d', cls=Dummy)
    assert isinstance(d, Dummy)
    c, d = symbols('c,d', cls=Dummy)
    assert isinstance(c, Dummy)
    assert isinstance(d, Dummy)
    raises(TypeError, lambda: Symbol())


def test_Dummy():
    assert Dummy() != Dummy()


def test_Dummy_force_dummy_index():
    raises(AssertionError, lambda: Dummy(dummy_index=1))
    assert Dummy('d', dummy_index=2) == Dummy('d', dummy_index=2)
    assert Dummy('d1', dummy_index=2) != Dummy('d2', dummy_index=2)
    d1 = Dummy('d', dummy_index=3)
    d2 = Dummy('d')
    # might fail if d1 were created with dummy_index >= 10**6
    assert d1 != d2
    d3 = Dummy('d', dummy_index=3)
    assert d1 == d3
    assert Dummy()._count == Dummy('d', dummy_index=3)._count


def test_as_dummy():
    x = Symbol('x')
    x1 = x.as_dummy()
    assert x1 != x
    assert x1 != x.as_dummy()

    x = Symbol('x', commutative=False)
    x1 = x.as_dummy()
    assert x1 != x
    assert x1.is_commutative is False


def test_lt_gt():
    from sympy import sympify as S
    x, y = Symbol('x'), Symbol('y')

    assert (x >= y) == GreaterThan(x, y)
    assert (x >= 0) == GreaterThan(x, 0)
    assert (x <= y) == LessThan(x, y)
    assert (x <= 0) == LessThan(x, 0)

    assert (0 <= x) == GreaterThan(x, 0)
    assert (0 >= x) == LessThan(x, 0)
    assert (S(0) >= x) == GreaterThan(0, x)
    assert (S(0) <= x) == LessThan(0, x)

    assert (x > y) == StrictGreaterThan(x, y)
    assert (x > 0) == StrictGreaterThan(x, 0)
    assert (x < y) == StrictLessThan(x, y)
    assert (x < 0) == StrictLessThan(x, 0)

    assert (0 < x) == StrictGreaterThan(x, 0)
    assert (0 > x) == StrictLessThan(x, 0)
    assert (S(0) > x) == StrictGreaterThan(0, x)
    assert (S(0) < x) == StrictLessThan(0, x)

    e = x**2 + 4*x + 1
    assert (e >= 0) == GreaterThan(e, 0)
    assert (0 <= e) == GreaterThan(e, 0)
    assert (e > 0) == StrictGreaterThan(e, 0)
    assert (0 < e) == StrictGreaterThan(e, 0)

    assert (e <= 0) == LessThan(e, 0)
    assert (0 >= e) == LessThan(e, 0)
    assert (e < 0) == StrictLessThan(e, 0)
    assert (0 > e) == StrictLessThan(e, 0)

    assert (S(0) >= e) == GreaterThan(0, e)
    assert (S(0) <= e) == LessThan(0, e)
    assert (S(0) < e) == StrictLessThan(0, e)
    assert (S(0) > e) == StrictGreaterThan(0, e)


def test_no_len():
    # there should be no len for numbers
    x = Symbol('x')
    raises(TypeError, lambda: len(x))


def test_ineq_unequal():
    S = sympify

    x, y, z = symbols('x,y,z')

    e = (
        S(-1) >= x, S(-1) >= y, S(-1) >= z,
        S(-1) > x, S(-1) > y, S(-1) > z,
        S(-1) <= x, S(-1) <= y, S(-1) <= z,
        S(-1) < x, S(-1) < y, S(-1) < z,
        S(0) >= x, S(0) >= y, S(0) >= z,
        S(0) > x, S(0) > y, S(0) > z,
        S(0) <= x, S(0) <= y, S(0) <= z,
        S(0) < x, S(0) < y, S(0) < z,
        S('3/7') >= x, S('3/7') >= y, S('3/7') >= z,
        S('3/7') > x, S('3/7') > y, S('3/7') > z,
        S('3/7') <= x, S('3/7') <= y, S('3/7') <= z,
        S('3/7') < x, S('3/7') < y, S('3/7') < z,
        S(1.5) >= x, S(1.5) >= y, S(1.5) >= z,
        S(1.5) > x, S(1.5) > y, S(1.5) > z,
        S(1.5) <= x, S(1.5) <= y, S(1.5) <= z,
        S(1.5) < x, S(1.5) < y, S(1.5) < z,
        S(2) >= x, S(2) >= y, S(2) >= z,
        S(2) > x, S(2) > y, S(2) > z,
        S(2) <= x, S(2) <= y, S(2) <= z,
        S(2) < x, S(2) < y, S(2) < z,
        x >= -1, y >= -1, z >= -1,
        x > -1, y > -1, z > -1,
        x <= -1, y <= -1, z <= -1,
        x < -1, y < -1, z < -1,
        x >= 0, y >= 0, z >= 0,
        x > 0, y > 0, z > 0,
        x <= 0, y <= 0, z <= 0,
        x < 0, y < 0, z < 0,
        x >= 1.5, y >= 1.5, z >= 1.5,
        x > 1.5, y > 1.5, z > 1.5,
        x <= 1.5, y <= 1.5, z <= 1.5,
        x < 1.5, y < 1.5, z < 1.5,
        x >= 2, y >= 2, z >= 2,
        x > 2, y > 2, z > 2,
        x <= 2, y <= 2, z <= 2,
        x < 2, y < 2, z < 2,

        x >= y, x >= z, y >= x, y >= z, z >= x, z >= y,
        x > y, x > z, y > x, y > z, z > x, z > y,
        x <= y, x <= z, y <= x, y <= z, z <= x, z <= y,
        x < y, x < z, y < x, y < z, z < x, z < y,

        x - pi >= y + z, y - pi >= x + z, z - pi >= x + y,
        x - pi > y + z, y - pi > x + z, z - pi > x + y,
        x - pi <= y + z, y - pi <= x + z, z - pi <= x + y,
        x - pi < y + z, y - pi < x + z, z - pi < x + y,
        True, False
    )

    left_e = e[:-1]
    for i, e1 in enumerate( left_e ):
        for e2 in e[i + 1:]:
            assert e1 != e2


def test_Wild_properties():
    # these tests only include Atoms
    x = Symbol("x")
    y = Symbol("y")
    p = Symbol("p", positive=True)
    k = Symbol("k", integer=True)
    n = Symbol("n", integer=True, positive=True)

    given_patterns = [ x, y, p, k, -k, n, -n, sympify(-3), sympify(3),
                       pi, Rational(3, 2), I ]

    integerp = lambda k: k.is_integer
    positivep = lambda k: k.is_positive
    symbolp = lambda k: k.is_Symbol
    realp = lambda k: k.is_real

    S = Wild("S", properties=[symbolp])
    R = Wild("R", properties=[realp])
    Y = Wild("Y", exclude=[x, p, k, n])
    P = Wild("P", properties=[positivep])
    K = Wild("K", properties=[integerp])
    N = Wild("N", properties=[positivep, integerp])

    given_wildcards = [ S, R, Y, P, K, N ]

    goodmatch = {
        S: (x, y, p, k, n),
        R: (p, k, -k, n, -n, -3, 3, pi, Rational(3, 2)),
        Y: (y, -3, 3, pi, Rational(3, 2), I ),
        P: (p, n, 3, pi, Rational(3, 2)),
        K: (k, -k, n, -n, -3, 3),
        N: (n, 3)}

    for A in given_wildcards:
        for pat in given_patterns:
            d = pat.match(A)
            if pat in goodmatch[A]:
                assert d[A] in goodmatch[A]
            else:
                assert d is None


def test_symbols():
    x = Symbol('x')
    y = Symbol('y')
    z = Symbol('z')

    assert symbols('x') == x
    assert symbols('x ') == x
    assert symbols(' x ') == x
    assert symbols('x,') == (x,)
    assert symbols('x, ') == (x,)
    assert symbols('x ,') == (x,)

    assert symbols('x , y') == (x, y)

    assert symbols('x,y,z') == (x, y, z)
    assert symbols('x y z') == (x, y, z)

    assert symbols('x,y,z,') == (x, y, z)
    assert symbols('x y z ') == (x, y, z)

    xyz = Symbol('xyz')
    abc = Symbol('abc')

    assert symbols('xyz') == xyz
    assert symbols('xyz,') == (xyz,)
    assert symbols('xyz,abc') == (xyz, abc)

    assert symbols(('xyz',)) == (xyz,)
    assert symbols(('xyz,',)) == ((xyz,),)
    assert symbols(('x,y,z,',)) == ((x, y, z),)
    assert symbols(('xyz', 'abc')) == (xyz, abc)
    assert symbols(('xyz,abc',)) == ((xyz, abc),)
    assert symbols(('xyz,abc', 'x,y,z')) == ((xyz, abc), (x, y, z))

    assert symbols(('x', 'y', 'z')) == (x, y, z)
    assert symbols(['x', 'y', 'z']) == [x, y, z]
    assert symbols(set(['x', 'y', 'z'])) == set([x, y, z])

    raises(ValueError, lambda: symbols(''))
    raises(ValueError, lambda: symbols(','))
    raises(ValueError, lambda: symbols('x,,y,,z'))
    raises(ValueError, lambda: symbols(('x', '', 'y', '', 'z')))

    a, b = symbols('x,y', real=True)
    assert a.is_real and b.is_real

    x0 = Symbol('x0')
    x1 = Symbol('x1')
    x2 = Symbol('x2')

    y0 = Symbol('y0')
    y1 = Symbol('y1')

    assert symbols('x0:0') == ()
    assert symbols('x0:1') == (x0,)
    assert symbols('x0:2') == (x0, x1)
    assert symbols('x0:3') == (x0, x1, x2)

    assert symbols('x:0') == ()
    assert symbols('x:1') == (x0,)
    assert symbols('x:2') == (x0, x1)
    assert symbols('x:3') == (x0, x1, x2)

    assert symbols('x1:1') == ()
    assert symbols('x1:2') == (x1,)
    assert symbols('x1:3') == (x1, x2)

    assert symbols('x1:3,x,y,z') == (x1, x2, x, y, z)

    assert symbols('x:3,y:2') == (x0, x1, x2, y0, y1)
    assert symbols(('x:3', 'y:2')) == ((x0, x1, x2), (y0, y1))

    a = Symbol('a')
    b = Symbol('b')
    c = Symbol('c')
    d = Symbol('d')

    assert symbols('x:z') == (x, y, z)
    assert symbols('a:d,x:z') == (a, b, c, d, x, y, z)
    assert symbols(('a:d', 'x:z')) == ((a, b, c, d), (x, y, z))

    aa = Symbol('aa')
    ab = Symbol('ab')
    ac = Symbol('ac')
    ad = Symbol('ad')

    assert symbols('aa:d') == (aa, ab, ac, ad)
    assert symbols('aa:d,x:z') == (aa, ab, ac, ad, x, y, z)
    assert symbols(('aa:d','x:z')) == ((aa, ab, ac, ad), (x, y, z))


    # issue 6675
    def sym(s):
        return str(symbols(s))
    assert sym('a0:4') == '(a0, a1, a2, a3)'
    assert sym('a2:4,b1:3') == '(a2, a3, b1, b2)'
    assert sym('a1(2:4)') == '(a12, a13)'
    assert sym(('a0:2.0:2')) == '(a0.0, a0.1, a1.0, a1.1)'
    assert sym(('aa:cz')) == '(aaz, abz, acz)'
    assert sym('aa:c0:2') == '(aa0, aa1, ab0, ab1, ac0, ac1)'
    assert sym('aa:ba:b') == '(aaa, aab, aba, abb)'
    assert sym('a:3b') == '(a0b, a1b, a2b)'
    assert sym('a-1:3b') == '(a-1b, a-2b)'
    assert sym(r'a:2\,:2' + chr(0)) == '(a0,0%s, a0,1%s, a1,0%s, a1,1%s)' % (
        (chr(0),)*4)
    assert sym('x(:a:3)') == '(x(a0), x(a1), x(a2))'
    assert sym('x(:c):1') == '(xa0, xb0, xc0)'
    assert sym('x((:a)):3') == '(x(a)0, x(a)1, x(a)2)'
    assert sym('x(:a:3') == '(x(a0, x(a1, x(a2)'
    assert sym(':2') == '(0, 1)'
    assert sym(':b') == '(a, b)'
    assert sym(':b:2') == '(a0, a1, b0, b1)'
    assert sym(':2:2') == '(00, 01, 10, 11)'
    assert sym(':b:b') == '(aa, ab, ba, bb)'

    raises(ValueError, lambda: symbols(':'))
    raises(ValueError, lambda: symbols('a:'))
    raises(ValueError, lambda: symbols('::'))
    raises(ValueError, lambda: symbols('a::'))
    raises(ValueError, lambda: symbols(':a:'))
    raises(ValueError, lambda: symbols('::a'))


def test_call():
    f = Symbol('f')
    assert f(2)
    raises(TypeError, lambda: Wild('x')(1))

def test_unicode():
    xu = Symbol(u'x')
    x = Symbol('x')
    assert x == xu

    raises(TypeError, lambda: Symbol(1))
