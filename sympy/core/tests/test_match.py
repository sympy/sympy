from sympy import abc, Function, Symbol, Wild, Derivative, sin, cos, Real, \
        Rational, exp, I, Integer, diff, Mul, var, oo, S, Add, Poly
from sympy.utilities.pytest import XFAIL


def test_symbol():
    x = Symbol('x')
    a,b,c,p,q = map(Wild, 'abcpq')

    e = x
    assert e.match(x) == {}
    assert e.match(a) == {a: x}

    e = Rational(5)
    assert e.match(c) == {c: 5}
    assert e.match(e) == {}
    assert e.match(e+1) == None

def test_add():
    x,y,a,b,c = map(Symbol, 'xyabc')
    p,q,r = map(Wild, 'pqr')

    e = a+b
    assert e.match(p+b) == {p: a}
    assert e.match(p+a) == {p: b}

    e = 1+b
    assert e.match(p+b) == {p: 1}

    e = a+b+c
    assert e.match(a+p+c) == {p: b}
    assert e.match(b+p+c) == {p: a}

    e = a+b+c+x
    assert e.match(a+p+x+c) == {p: b}
    assert e.match(b+p+c+x) == {p: a}
    assert e.match(b) == None
    assert e.match(b+p) == {p: a+c+x}
    assert e.match(a+p+c) == {p: b+x}
    assert e.match(b+p+c) == {p: a+x}

    e = 4*x+5
    assert e.match(4*x+p) == {p: 5}
    assert e.match(3*x+p) == {p: x+5}
    assert e.match(p*x+5) == {p: 4}


def test_power():
    x,y,a,b,c = map(Symbol, 'xyabc')
    p,q,r = map(Wild, 'pqr')

    e = (x+y)**a
    assert e.match(p**q) == {p: x+y, q: a}
    assert e.match(p**p) == None

    e = (x+y)**(x+y)
    assert e.match(p**p) == {p: x+y}
    assert e.match(p**q) == {p: x+y, q: x+y}

    e = (2*x)**2
    assert e.match(p*q**r) == {p: 4, q: x, r: 2}

    e = Integer(1)
    assert e.match(x**p) == {p: 0}

def test_match_exclude():

    x = Symbol('x')
    y = Symbol('y')
    p = Wild("p", exclude=[x, y])
    q = Wild("q", exclude=[x, y])
    r = Wild("r", exclude=[x, y])

    e = 3/(4*x+5)
    assert e.match(3/(p*x+q)) == {p: 4, q: 5}

    e = 3/(4*x+5)
    assert e.match(p/(q*x+r)) == {p: 3, q: 4, r: 5}

    e = 2/(x+1)
    assert e.match(p/(q*x+r)) == {p: 2, q: 1, r: 1}

    e = 1/(x+1)
    assert e.match(p/(q*x+r)) == {p: 1, q: 1, r: 1}

    e = 4*x+5
    assert e.match(p*x+q) == {p: 4, q: 5}

    e = 4*x+5*y+6
    assert e.match(p*x+q*y+r) == {p: 4, q: 5, r: 6}

def test_mul():
    x,y,a,b,c = map(Symbol, 'xyabc')
    p,q = map(Wild, 'pq')

    e = 4*x
    assert e.match(p*x) == {p: 4}
    assert e.match(p*y) == {p: 4*x/y}

    e = a*x*b*c
    assert e.match(p*x) == {p: a*b*c}
    assert e.match(c*p*x) == {p: a*b}

    e = (a+b)*(a+c)
    assert e.match((p+b)*(p+c)) == {p: a}

    e = x
    assert e.match(p*x) == {p: 1}

    e = exp(x)
    assert e.match(x**p*exp(x*q)) == {p: 0, q: 1}

    e = I*Poly(x, x)
    assert e.match(I*p) == {p: Poly(x, x)}

def test_complex():
    a,b,c = map(Symbol, 'abc')
    x,y = map(Wild, 'xy')

    (1+I).match(x+I) == {x : 1}
    (a+I).match(x+I) == {x : a}
    (a+b*I).match(x+y*I) == {x : a, y : b}
    (2*I).match(x*I) == {x : 2}
    (a*I).match(x*I) == {x : a}
    (a*I).match(x*y) == {x : a, y : I}
    (2*I).match(x*y) == {x : 2, y : I}

def test_functions():
    from sympy.core.function import WildFunction
    x = Symbol('x')
    g = WildFunction('g')
    p = Wild('p')
    q = Wild('q')

    f = cos(5*x)
    notf = x
    assert f.match(p*cos(q*x)) == {p: 1, q: 5}
    assert f.match(p*g) == {p: 1, g: cos(5*x)}
    assert notf.match(g) == None

@XFAIL
def test_functions_X1():
    assert f.match(p*g(q*x)) == {p: 1, g: cos, q: 5}

def test_interface():
    x,y = map(Symbol, 'xy')
    p,q = map(Wild, 'pq')

    assert (x+1).match(p+1) == {p: x}
    assert (x*3).match(p*3) == {p: x}
    assert (x**3).match(p**3) == {p: x}
    assert (x*cos(y)).match(p*cos(q)) == {p: x, q: y}

    assert (x*y).match(p*q) in [{p:x, q:y}, {p:y, q:x}]
    assert (x+y).match(p+q) in [{p:x, q:y}, {p:y, q:x}]
    assert (x*y+1).match(p*q) in [{p:1, q:1+x*y}, {p:1+x*y, q:1}]

def test_derivative1():
    x,y = map(Symbol, 'xy')
    p,q = map(Wild, 'pq')

    f = Function('f',nargs=1)
    fd = Derivative(f(x), x)

    assert fd.match(p) == {p: fd}
    assert (fd+1).match(p+1) == {p: fd}
    assert (fd).match(fd) == {}
    assert (3*fd).match(p*fd) != None
    p = Wild("p", exclude=[x])
    q = Wild("q", exclude=[x])
    assert (3*fd-1).match(p*fd + q) == {p: 3, q: -1}

def test_derivative_bug1():
    f = Function("f")
    x = Symbol("x")
    a = Wild("a", exclude=[f, x])
    b = Wild("b", exclude=[f])
    pattern = a * Derivative(f(x), x, x) + b
    expr = Derivative(f(x), x)+x**2
    d1 = {b: x**2}
    d2 = pattern.matches(expr, d1, evaluate=True)
    assert d2 == None

def test_derivative2():
    f = Function("f")
    x = Symbol("x")
    a = Wild("a", exclude=[f, x])
    b = Wild("b", exclude=[f])
    e = Derivative(f(x), x)
    assert e.match(Derivative(f(x), x)) == {}
    assert e.match(Derivative(f(x), x, x)) == None
    e = Derivative(f(x), x, x)
    assert e.match(Derivative(f(x), x)) == None
    assert e.match(Derivative(f(x), x, x)) == {}
    e = Derivative(f(x), x)+x**2
    assert e.match(a*Derivative(f(x), x) + b) == {a: 1, b: x**2}
    assert e.match(a*Derivative(f(x), x, x) + b) == None
    e = Derivative(f(x), x, x)+x**2
    assert e.match(a*Derivative(f(x), x) + b) == None
    assert e.match(a*Derivative(f(x), x, x) + b) == {a: 1, b: x**2}

def test_match_deriv_bug1():
    n = Function('n')
    l = Function('l')

    x = Symbol('x')
    p = Wild('p')

    e = diff(l(x), x)/x - diff(diff(n(x), x), x)/2 - \
        diff(n(x), x)**2/4 + diff(n(x), x)*diff(l(x), x)/4
    e = e.subs(n(x), -l(x))
    t = x*exp(-l(x))
    t2 = t.diff(x, x)/t
    assert e.match( (p*t2).expand() ) == {p: -Rational(1)/2}

def test_match_bug2():
    x,y = map(Symbol, 'xy')
    p,q,r = map(Wild, 'pqr')
    res = (x+y).match(p+q+r)
    assert (p+q+r).subs(res) == x+y

def test_match_bug3():
    x,a,b = map(Symbol, 'xab')
    p = Wild('p')
    assert (b*x*exp(a*x)).match(x*exp(p*x)) == None

def test_match_bug4():
    x = Symbol('x')
    p = Wild('p')
    e = x
    assert e.match(-p*x) == {p: -1}

def test_match_bug5():
    x = Symbol('x')
    p = Wild('p')
    e = -x
    assert e.match(-p*x) == {p: 1}

def test_match_bug6():
    x = Symbol('x')
    p = Wild('p')
    e = x
    assert e.match(3*p*x) == {p: Rational(1)/3}

def test_behavior1():
    x = Symbol('x')
    p = Wild('p')
    e = 3*x**2
    a = Wild('a', exclude = [x])
    assert e.match(a*x) == None
    assert e.match(p*x) == {p: 3*x}

def test_behavior2():
    x = Symbol('x')
    p = Wild('p')

    e = Rational(6)
    assert e.match(2*p) == {p: 3}

    e = 3*x + 3 + 6/x
    a = Wild('a', exclude = [x])
    assert e.expand().match(a*x**2 + a*x + 2*a) == None
    assert e.expand().match(p*x**2 + p*x + 2*p) == {p: 3/x}

def test_match_polynomial():
    x = Symbol('x')
    a = Wild('a', exclude=[x])
    b = Wild('b', exclude=[x])
    c = Wild('c', exclude=[x])
    d = Wild('d', exclude=[x])

    eq = 4*x**3 + 3*x**2 + 2*x + 1
    pattern = a*x**3 + b*x**2 + c*x + d
    assert eq.match(pattern) == {a: 4, b: 3, c: 2, d: 1}
    assert (eq-3*x**2).match(pattern) == {a: 4, b: 0, c: 2, d: 1}

def test_exclude():
    x,y,a = map(Symbol, 'xya')
    p = Wild('p', exclude=[1,x])
    q = Wild('q', exclude=[x])
    r = Wild('r', exclude=[sin,y])

    assert sin(x).match(r) == None
    assert cos(y).match(r) == None

    e = 3*x**2 + y*x + a
    assert e.match(p*x**2 + q*x + r) == {p: 3, q: y, r: a}

    e = x+1
    assert e.match(x+p) == None
    assert e.match(p+1) == None
    assert e.match(x+1+p) == {p: 0}

    e = cos(x) + 5*sin(y)
    assert e.match(r) == None
    assert e.match(cos(y) + r) == None
    assert e.match(r + p*sin(q)) == {r: cos(x), p: 5, q: y}

def test_floats():
    a,b = map(Wild, 'ab')

    e = cos(0.12345)**2
    r = e.match(a*cos(b)**2)
    assert r == {a: 1, b: Real(0.12345)}

def test_Derivative_bug1():
    f = Function("f")
    x = abc.x
    a = Wild("a", exclude=[f(x)])
    b = Wild("b", exclude=[f(x)])
    eq = f(x).diff(x)
    assert eq.match(a*Derivative(f(x), x) + b) == {a: 1, b: 0}


def test_match_wild_wild():
    p = Wild('p')
    q = Wild('q')
    r = Wild('r')

    assert p.match(q+r)  in  [ {q: p, r: 0} , {q: 0, r: p} ]
    assert p.match(q*r)  in  [ {q: p, r: 1} , {q: 1, r: p} ]


    p = Wild('p')
    q = Wild('q', exclude=[p])
    r = Wild('r')

    assert p.match(q+r) == {q: 0, r: p}
    assert p.match(q*r) == {q: 1, r: p}

    p = Wild('p')
    q = Wild('q', exclude=[p])
    r = Wild('r', exclude=[p])

    assert p.match(q+r) == None
    assert p.match(q*r) == None

def test_combine_inverse():
    x, y = var("x y")
    assert Mul._combine_inverse(x*I*y, x*I) == y
    assert Mul._combine_inverse(x*I*y, y*I) == x
    assert Mul._combine_inverse(oo*I*y, y*I) == oo
    assert Mul._combine_inverse(oo*I*y, oo*I) == y
    assert Add._combine_inverse(oo, oo) == S(0)
    assert Add._combine_inverse(oo*I, oo*I) == S(0)
