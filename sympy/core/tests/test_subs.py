from sympy import Symbol, Wild, sin, cos, exp, sqrt, pi, Function, Derivative,\
        abc, Integer, Eq, symbols, Add, I, Real, log, Rational, Lambda, atan2

def test_subs():
    n3=Rational(3)
    n2=Rational(2)
    n6=Rational(6)
    x=Symbol("x")
    c=Symbol("c")
    e=x
    e=e.subs(x,n3)
    assert e == Rational(3)

    e=2*x
    assert e == 2*x
    e=e.subs(x,n3)
    assert e == Rational(6)

def test_trigonometric():
    x = Symbol('x')
    n3 = Rational(3)
    e=(sin(x)**2).diff(x)
    assert e == 2*sin(x)*cos(x)
    e=e.subs(x,n3)
    assert e == 2*cos(n3)*sin(n3)

    e=(sin(x)**2).diff(x)
    assert e == 2*sin(x)*cos(x)
    e=e.subs(sin(x),cos(x))
    assert e == 2*cos(x)**2

    assert exp(pi).subs(exp, sin) == 0
    assert cos(exp(pi)).subs(exp, sin) == 1

def test_powers():
    x = Symbol('x')
    assert sqrt(1 - sqrt(x)).subs(x, 4) == I
    assert (sqrt(1-x**2)**3).subs(x, 2) == - 3 * I * sqrt(3)
    assert (x ** Rational(1,3)).subs(x, 27) == 3
    assert (x ** Rational(1,3)).subs(x, -27) == 3 * (-1) ** Rational(1,3)
    assert ((-x) ** Rational(1,3)).subs(x, 27) == 3 * (-1) ** Rational(1,3)

def test_logexppow():   # no eval()
    x = Symbol("x")
    w = Symbol("dummy :)")
    e = (3**(1+x)+2**(1+x))/(3**x+2**x)
    assert e.subs(2**x, w) != e
    assert e.subs(exp(x*log(Rational(2))),w) != e

def test_bug():
    x1=Symbol("x1")
    x2=Symbol("x2")
    y=x1*x2
    y.subs(x1,Real(3.0))

def test_subbug1():
    x=Symbol("x")
    e=(x**x).subs(x,1)
    e=(x**x).subs(x,1.0)

def test_subbug2():
    # Ensure this does not cause infinite recursion
    x = Symbol('x')
    assert Real(7.7).epsilon_eq(abs(x).subs(x, -7.7))

def test_dict():
    x = Symbol('x')
    a,b,c = map(Wild, 'abc')

    f = 3*cos(4*x)
    r = f.match(a*cos(b*x))
    assert r == {a: 3, b: 4}
    e =  a/b * sin(b*x)
    assert e._subs_dict(r) == r[a]/r[b] * sin(r[b]*x)
    assert e._subs_dict(r) == 3 * sin(4*x) / 4


def test_dict_ambigous():   # see #467
    x = Symbol('x')
    y = Symbol('y')
    z = Symbol('z')

    f = x*exp(x)
    g = z*exp(z)

    df= {x:y, exp(x): y}
    dg= {z:y, exp(z): y}

    assert f._subs_dict(df) == y**2
    assert g._subs_dict(dg) == y**2

    # and this is how order can affect the result
    assert f .subs(x,y) .subs(exp(x),y)  == y*exp(y)
    assert f .subs(exp(x),y) .subs(x,y)  == y**2


def test_deriv_sub_bug3():
    x = Symbol("x")
    y = Symbol("y")
    f = Function("f")
    pat = Derivative(f(x), x, x)
    assert pat.subs(y, y**2) == Derivative(f(x), x, x)
    assert pat.subs(y, y**2) != Derivative(f(x), x)

def test_equality_subs1():
    f = Function("f")
    x = abc.x
    eq = Eq(f(x)**2, x)
    res = Eq(Integer(16), x)
    assert eq.subs(f(x), 4) == res

def test_equality_subs2():
    f = Function("f")
    x = abc.x
    eq = Eq(f(x)**2, 16)
    assert bool(eq.subs(f(x), 3)) == False
    assert bool(eq.subs(f(x), 4)) == True

def test_issue643():
    x = Symbol('x')
    y = Symbol('y')

    e = sqrt(x)*exp(y)
    assert e.subs(sqrt(x), 1)   == exp(y)

def test_subs_dict1():
    x, y = symbols('xy')
    assert (1+x*y).subs(x, pi) == 1 + pi*y
    assert (1+x*y).subs({x:pi, y:2}) == 1 + 2*pi

def test_subs_dict2():
    x = Symbol('x')
    a,b,c = map(Wild, 'abc')

    f = 3*cos(4*x)
    r = f.match(a*cos(b*x))
    assert r == {a: 3, b: 4}
    e =  a/b * sin(b*x)
    assert e.subs(r) == r[a]/r[b] * sin(r[b]*x)
    assert e.subs(r) == 3 * sin(4*x) / 4

def test_mul():
    x, y, z = map(Symbol, 'xyz')
    assert (x*y*z).subs(z*x,y) == y**2
    assert (2*x*y).subs(5*x*y,z) == 2*z/5

def test_add():
    a, b, c, d, x = abc.a, abc.b, abc.c, abc.d, abc.x
    assert (a**2 - b - c).subs(a**2 - b, d) in [d - c, a**2 - b - c]
    assert (a**2 - c).subs(a**2 - c, d) == d
    assert (a**2 - b - c).subs(a**2 - c, d) in [d - b, a**2 - b - c]
    assert (a**2 - x - c).subs(a**2 - c, d) in [d - x, a**2 - x - c]
    assert (a**2 - b - sqrt(a)).subs(a**2 - sqrt(a), c) == c - b
    assert (a+b+exp(a+b)).subs(a+b,c) == c + exp(c)
    assert (c+b+exp(c+b)).subs(c+b,a) == a + exp(a)

    # this should work everytime:
    e = a**2 - b - c
    assert e.subs(Add(*e.args[:2]), d) == d + e.args[2]
    assert e.subs(a**2 - c, d) == d - b

def test_subs_issue910():
    assert (I*Symbol("a")).subs(1, 2) == I*Symbol("a")

def test_subs_subs_nums():
    x = Symbol("x")
    assert sin(1).subs(1, 2) == sin(2)
    assert sin(2).subs(1, 3) == sin(2)
    assert (2*x).subs(1, 3) == 2*x
    assert (2*x).subs(2, 3) == 3*x
    assert (2*x).subs(x, 3) == 6

def test_functions_subs():
    x, y = map(Symbol, 'xy')
    f, g = map(Function, 'fg')
    l = Lambda(x, y, sin(x) + y)
    assert (g(y, x)+cos(x)).subs(g, l) == sin(y) + x + cos(x)
    assert (f(x)**2).subs(f, sin) == sin(x)**2
    assert (f(x,y)).subs(f,log) == log(x,y)
    assert (f(x,y)).subs(f,sin) == f(x,y)
    assert (sin(x)+atan2(x,y)).subs([[atan2,f],[sin,g]]) == f(x,y) + g(x)
    assert (g(f(x+y, x))).subs([[f, l], [g, exp]]) == exp(x + sin(x + y))
