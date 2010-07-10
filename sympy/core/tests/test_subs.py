from sympy import (Symbol, Wild, sin, cos, exp, sqrt, pi, Function, Derivative,
        abc, Integer, Eq, symbols, Add, I, Real, log, Rational, Lambda, atan2)

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
    x, y = symbols('x,y')
    assert (1+x*y).subs(x, pi) == 1 + pi*y
    assert (1+x*y).subs({x:pi, y:2}) == 1 + 2*pi
    c2,c3,q1p,q2p,c1,s1,s2,s3= symbols('c2,c3,q1p,q2p,c1,s1,s2,s3')
    test=c2**2*q2p*c3 + c1**2*s2**2*q2p*c3 + s1**2*s2**2*q2p*c3 \
        - c1**2*q1p*c2*s3 - s1**2*q1p*c2*s3
    assert test.subs({c1**2 : 1-s1**2, c2**2 : 1-s2**2, c3**3: 1-s3**2}) \
        == c3*q2p*(1 - s2**2) + c3*q2p*s2**2*(1 - s1**2) - c2*q1p*s3*(1 - s1**2) \
        + c3*q2p*s1**2*s2**2 - c2*q1p*s3*s1**2

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
    x, y, z, a = map(Symbol, 'xyza')
    assert (x*y*z).subs(z*x,y) == y**2
    assert (z*x).subs(1/x,z) == z*x
    assert (x*y/z).subs(1/z,a) == a*x*y
    assert (x*y/z).subs(x/z,a) == a*y
    assert (x*y/z).subs(y/z,a) == a*x
    assert (x*y/z).subs(x/z,1/a) == y/a
    assert (x*y/z).subs(x,1/a) == y/(z*a)
    assert (2*x*y).subs(5*x*y,z) != 2*z/5

def test_subs_simple():
    # Define symbols
    a = symbols('a', commutative = True)
    x = symbols('x', commutative = False)

    """ SIMPLE TESTS """
    assert (2*a ).subs(1,3) == 2*a
    assert (2*a ).subs(2,3) == 3*a
    assert (2*a ).subs(a,3) == 6
    assert sin(2).subs(1,3) == sin(2)
    assert sin(2).subs(2,3) == sin(3)
    assert sin(a).subs(a,3) == sin(3)

    assert (2*x ).subs(1,3) == 2*x
    assert (2*x ).subs(2,3) == 3*x
    assert (2*x ).subs(x,3) == 6
    assert sin(x).subs(x,3) == sin(3)

def test_subs_constants():
    # Define symbols
    a,b = symbols('a,b', commutative = True)
    x,y = symbols('x,y', commutative = False)

    """ CONSTANTS TESTS """
    assert (a*b  ).subs(2*a,1) == a*b
    assert (1.5*a*b).subs(a,1) == 1.5*b
    assert (2*a*b).subs(2*a,1) == b
    assert (2*a*b).subs(4*a,1) == 2*a*b

    assert (x*y  ).subs(2*x,1) == x*y
    assert (1.5*x*y).subs(x,1) == 1.5*y
    assert (2*x*y).subs(2*x,1) == y
    assert (2*x*y).subs(4*x,1) == 2*x*y

def test_subs_commutative():
    # Define symbols
    a,b,c,d,K = symbols('a,b,c,d,K', commutative = True)

    """ COMMUTATIVE TESTS """
    assert (a*b    ).subs(a*b,K) == K
    assert (a*b*a*b).subs(a*b,K) == K**2
    assert (a*a*b*b).subs(a*b,K) == K**2
    assert (a*b*c*d).subs(a*b*c,K) == d*K
    assert (a*b**c ).subs(a,K) == K*b**c
    assert (a*b**c ).subs(b,K) == a*K**c
    assert (a*b**c ).subs(c,K) == a*b**K
    assert (a*b*c*b*a  ).subs(a*b,K) == c*K**2
    assert (a**3*b**2*a).subs(a*b,K) == a**2*K**2

def test_subs_noncommutative():
    # Define symbols
    w,x,y,z,L = symbols('w,x,y,z,L', commutative = False)

    """ NONCOMMUTATIVE TESTS """
    assert (x*y    ).subs(x*y,L) == L
    assert (w*y*x  ).subs(x*y,L) == w*y*x
    assert (w*x*y*z).subs(x*y,L) == w*L*z
    assert (x*y*x*y).subs(x*y,L) == L**2
    assert (x*x*y  ).subs(x*y,L) == x*L
    assert (x*x*y*y).subs(x*y,L) == x*L*y
    assert (w*x*y  ).subs(x*y*z,L) == w*x*y
    assert (x*y**z     ).subs(x,L) == L*y**z
    assert (x*y**z     ).subs(y,L) == x*L**z
    assert (x*y**z     ).subs(z,L) == x*y**L
    assert (w*x*y*z*x*y).subs(x*y*z,L) == w*L*x*y
    assert (w*x*y*y*w*x*x*y*x*y*y*x*y).subs(x*y,L) == w*L*y*w*x*L**2*y*L

def test_subs_basic_funcs():
    # Define symbols
    a,b,c,d,K = symbols('a,b,c,d,K', commutative = True)
    w,x,y,z,L = symbols('w,x,y,z,L', commutative = False)

    """ OTHER OPERATION TESTS"""
    assert (x+y  ).subs(x+y,L) == L
    assert (x-y  ).subs(x-y,L) == L
    assert (x/y  ).subs(x,L) == L/y
    assert (x**y ).subs(x,L) == L**y
    assert (x**y ).subs(y,L) == x**L
    assert ( (a-c)/b  ).subs(b,K) == (a-c)/K
    assert (exp(x*y-z)).subs(x*y,L) == exp(L-z)
    assert (a*exp(x*y-w*z)+b*exp(x*y+w*z)).subs(z,0) == a*exp(x*y)+b*exp(x*y)
    assert ((a-b)/(c*d-a*b)).subs(c*d-a*b,K) == (a-b)/K
    assert (w*exp(a*b-c)*x*y/4).subs(x*y,L) == w*exp(a*b-c)*L/4
    #assert (a/(b*c)).subs(b*c,K) == a/K,'Failed'; print '.' #FAILS DIVISION

def test_subs_wild():
    # Define symbols
    R = Wild('R'); S = Wild('S'); T = Wild('T'); U = Wild('U')

    """ WILD TESTS """
    assert (R*S ).subs(R*S,T) == T
    assert (S*R ).subs(R*S,T) == T
    assert (R+S ).subs(R+S,T) == T
    assert (R**S).subs(R,T) == T**S
    assert (R**S).subs(S,T) == R**T
    assert (R*S**T).subs(R,U) == U*S**T
    assert (R*S**T).subs(S,U) == R*U**T
    assert (R*S**T).subs(T,U) == R*S**U

def test_subs_mixed():
    # Define symbols
    a,b,c,d,K = symbols('a,b,c,d,K', commutative = True)
    w,x,y,z,L = symbols('w,x,y,z,L', commutative = False)
    R,S,T,U = Wild('R'), Wild('S'), Wild('T'), Wild('U')

    """ MIXED TESTS """
    assert (    a*x*y     ).subs(x*y,L) == a*L
    assert (  a*b*x*y*x   ).subs(x*y,L) == a*b*L*x
    assert (R*x*y*exp(x*y)).subs(x*y,L) == R*L*exp(L)
    assert (     a*x*y*y*x-x*y*z*exp(a*b)   ).subs(x*y,L) == a*L*y*x-L*z*exp(a*b)
    assert (c*y*x*y*x**(R*S-a*b)-T*(a*R*b*S)).subs(x*y,L).subs(a*b,K).subs(R*S,U) == c*y*L*x**(U-K)-T*(U*K)

def test_division():
    a,b,c = symbols('a,b,c', commutative = True)
    x,y,z = symbols('x,y,z', commutative = True)
    assert (    1/a   ).subs(a,c)  == 1/c
    assert (   1/a**2 ).subs(a,c)  == 1/c**2
    assert (   1/a**2 ).subs(a,-2) == Rational(1,4)
    assert ( -(1/a**2)).subs(a,-2) == -Rational(1,4)

    assert (    1/x   ).subs(x,z)  == 1/z
    assert (   1/x**2 ).subs(x,z)  == 1/z**2
    assert (   1/x**2 ).subs(x,-2) == Rational(1,4)
    assert ( -(1/x**2)).subs(x,-2) == -Rational(1,4)

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

def test_subs_iter():
    x, y = symbols('x,y')
    assert x.subs(reversed([[x, y]])) == y
    it = iter([[x, y]])
    assert x.subs(it) == y
