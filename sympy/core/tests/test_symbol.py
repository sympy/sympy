from sympy import (Symbol, Wild, Inequality, StrictInequality, pi, I, Rational,
    sympify, symbols, Dummy, S)

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
    c,d = symbols('c,d', cls=Dummy)
    assert isinstance(c, Dummy)
    assert isinstance(d, Dummy)
    raises(TypeError, 'Symbol()')

def test_Dummy():
    assert Dummy() != Dummy()
    Dummy._count = 0
    d1 = Dummy()
    Dummy._count = 0
    assert d1 == Dummy()

def test_as_dummy_nondummy():
    x = Symbol('x')
    x1 = x.as_dummy()
    assert x1 != x
    assert x1 != x.as_dummy()
    # assert x == x1.as_nondummy()

    x = Symbol('x', commutative = False)
    x1 = x.as_dummy()
    assert x1 != x
    assert x1.is_commutative == False
    # assert x == x1.as_nondummy()

def test_lt_gt():
    x, y = Symbol('x'), Symbol('y')

    assert (x <= y) == Inequality(x, y)
    assert (x >= y) == Inequality(y, x)
    assert (x <= 0) == Inequality(x, 0)
    assert (x >= 0) == Inequality(0, x)

    assert (x < y) == StrictInequality(x, y)
    assert (x > y) == StrictInequality(y, x)
    assert (x < 0) == StrictInequality(x, 0)
    assert (x > 0) == StrictInequality(0, x)

    assert (x**2+4*x+1 > 0) == StrictInequality(0, x**2+4*x+1)

def test_no_len():
    # there should be no len for numbers
    x = Symbol('x')
    xxl = Symbol('xxl')
    raises(TypeError, "len(x)")
    raises(TypeError, "len(xxl)")

def test_Wild_properties():
    # these tests only include Atoms
    x   = Symbol("x")
    y   = Symbol("y")
    p   = Symbol("p", positive=True)
    k   = Symbol("k", integer=True)
    r   = Symbol("r", real=True)
    n   = Symbol("n", integer=True, positive=True)

    given_patterns = [ x, y, p, k, -k, n, -n, sympify(-3), sympify(3), pi, Rational(3,2), I ]

    integerp = lambda k : k.is_integer
    positivep = lambda k : k.is_positive
    symbolp = lambda k : k.is_Symbol
    realp = lambda k : k.is_real

    S = Wild("S", properties=[symbolp])
    R = Wild("R", properties=[realp])
    Y = Wild("Y", exclude=[x,p,k,n])
    P = Wild("P", properties=[positivep])
    K = Wild("K", properties=[integerp])
    N = Wild("N", properties=[positivep, integerp])

    given_wildcards = [ S, R, Y, P, K, N ]

    goodmatch = {
        S : (x,y,p,k,n),
        R : (p,k,-k,n,-n,-3,3,pi,Rational(3,2)),
        Y : (y,-3,3,pi,Rational(3,2),I ),
        P : (p, n,3,pi, Rational(3,2)),
        K : (k,-k,n,-n,-3,3),
        N : (n,3)}

    for A in given_wildcards:
        for pat in given_patterns:
            d = pat.match(A)
            if pat in goodmatch[A]:
                assert d[A] in goodmatch[A]
            else:
                assert d == None

def test_Pure():
    assert (S.Pure == S.Pure) == True

    assert (S.Pure == Symbol('x')) == False
    assert (Symbol('x') == S.Pure) == False

    assert (S.Pure == Dummy('x')) == False
    assert (Dummy('x') == S.Pure) == False

    assert (S.Pure == Symbol('x', commutative=False)) == False
    assert (Symbol('x', commutative=False) == S.Pure) == False

    assert (S.Pure == Symbol('pure')) == False
    assert (Symbol('pure') == S.Pure) == False

    assert (S.Pure == 1) == False
    assert (S.Pure == I) == False

    assert (S.Pure != S.Pure) == False

    assert (S.Pure != Symbol('x')) == True
    assert (Symbol('x') != S.Pure) == True

    assert (S.Pure != Dummy('x')) == True
    assert (Dummy('x') != S.Pure) == True

    assert (S.Pure != Symbol('x', commutative=False)) == True
    assert (Symbol('x', commutative=False) != S.Pure) == True

    assert (S.Pure != Symbol('pure')) == True
    assert (Symbol('pure') != S.Pure) == True

    assert (S.Pure != 1) == True
    assert (S.Pure != I) == True

def test_symbols():
    x = Symbol('x')
    y = Symbol('y')
    z = Symbol('z')

    assert symbols('') is None

    assert symbols('x') == x
    assert symbols('x,') == (x,)
    assert symbols('x ') == (x,)

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

    assert symbols('x,,y,,z') == (x, y, z)
    assert symbols(('x', '', 'y', '', 'z')) == (x, y, z)

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

def test_call():
    f = Symbol('f')
    assert f(2)
