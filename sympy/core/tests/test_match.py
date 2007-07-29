from sympy import Rational, Symbol, Derivative, Function, cos, exp


def test_symbol():
    x,y,a,b,c,p,q = map(Symbol, 'xyabcpq')
    e = a
    assert e.match(b,[p]) == None
    assert e.match(a,[p]) == {}
    assert e.match(p,[p]) == {p: a}
    assert e.match(p,[p,q]) == {p: a}
    assert e.match(p,[p]) != {p: b}

    e = x
    assert e.match(x,[b,c]) == {}

    e = Rational(5)
    assert e.match(c, [c]) == {c: 5}
    assert e.match(c, [b,c]) == {c: 5}
    assert e.match(e, [b]) == {}
    assert e.match(e, [b,c]) == {}
    assert e.match(e+1, [b,c]) == None

def test_add():
    x,y,a,b,c,p,q = map(Symbol, 'xyabcpq')
    e = a+b
    assert e.match(p+b,[p]) == {p: a}
    assert e.match(p+a,[p]) == {p: b}
    e = 1+b
    assert e.match(p+b, [p,b]) in [{p:1, b:b}, {p:b, b:1}]
    e = a+b+c
    assert e.match(a+p+c,[p]) == {p: b}
    assert e.match(b+p+c,[p]) == {p: a}
    e = a+b+c+x
    assert e.match(a+p+x+c,[p]) == {p: b}
    assert e.match(b+p+c+x,[p]) == {p: a}
    assert e.match(b,[p]) == None
    assert e.match(b+p,[p]) == {p: a+c+x}
    assert e.match(a+p+c,[p]) == {p: b+x}
    assert e.match(b+p+c,[p]) == {p: a+x}
    e = 4*x+5
    assert e.match(4*x+c,[c]) == {c: 5}
    assert e.match(3*x+c,[c]) == {c: x+5}
    assert e.match(b*x+5,[b]) == {b: 4}
    assert e.match(b*x+c,[b,c]) == {b: 4, c: 5}
    e = 4*x+5*y+6
    assert e.match(a*x+b*y+c, [a,b,c]) == {a: 4, b: 5, c: 6}

def _test_power():
    x,y,a,b,c,p,q = map(Symbol, 'xyabcpq')
    e = (x+y)**a
    assert e.match(p**q,[p,q]) == {p: x+y, q: a}
    assert e.match(p**p,[p]) == None
    e = (x+y)**(x+y)
    assert e.match(p**p,[p]) == {p: x+y}
    assert e.match(p**q,[p,q]) == {p: x+y, q: x+y}

    e = 3/(4*x+5)
    assert e.match(3/(a*x+b), [a,b]) == {a: 4, b: 5}

    e = 3/(4*x+5)
    assert e.match(a/(b*x+c),[a,b,c]) == {a: 3, b: 4, c: 5}

    e = 3*x**2+y*x+p
    # XXX Match made currently == {a: y/x, b: 2*x, c: p}
    assert e.match(a*x**2+b*x+c,[a,b,c]) == {a: 3, b: y, c: p}

    e = 2/(x+1)
    assert e.match(a/(b*x+c),[a,b,c]) == {a: 2, b: 1, c: 1}

    e = 1/(x+1)
    assert e.match(a/(b*x+c),[a,b,c]) == {a: 1, b: 1, c: 1}

    assert ((2*x)**2).match(a*q**p, [a, q, p]) == {p: 2, q: x, a: 4}

def _test_mul():
    x,y,a,b,c,p,q = map(Symbol, 'xyabcpq')
    e = 4*x
    assert e.match(b*x,[b]) == {b: 4}
    assert e.match(b*x,[b,c]) == {b: 4}
    # XXX A match is made {b: 4*x/y}, but is wanted behaviour?
    assert e.match(b*y,[b]) == None
    assert e.match(b*y,[a,b,c]) == None
    e = a*x*b*c
    assert e.match(p*x,[p]) == {p: a*b*c}
    assert e.match(c*p*x,[p]) == {p: a*b}
    e = (a+b)*(a+c)
    assert e.match((p+b)*(p+c),[p]) == {p: a}

    e = x
    assert e.match(a*x,[a]) == {a: 1}
    
def test_complex():
    x,y,a,b,c = map(Symbol, 'xyabc')
    from sympy import I
    (1+I).match(x+I, [x]) == {x : 1}
    (a+I).match(x+I, [x]) == {x : a}
    (a+b*I).match(x+y*I) == {x : a, y : b}
    (2*I).match(x*I) == {x : 2}
    (a*I).match(x*I) == {x : a}
    (a*I).match(x*y) == {x : a, y : I}
    (2*I).match(x*y) == {x : 2, y : I}

def _test_functions():
    # XXX The match currently is {a: 1/cos(5*x)*cos(5*x), b: 5}
    #     but this is clearly a=1. Problem with evaluation?
    x,y,a,b,c = map(Symbol, 'xyabc')
    f = cos(5*x)
    assert f.match(a*cos(b*x), [a,b]) == {a: 1, b: 5}

def test_interface():
    x,y,a,b,c,p,q = map(Symbol, 'xyabcpq')
    #assert (x+1).match(a+1) == {a: x}
    #assert (x*3).match(a*3) == {a: x}
    #assert (x**3).match(a**3) == {a: x}
    assert list((a*cos(b)).atoms(type=Symbol)) in [[a,b], [b, a]]
    #assert (x*cos(y)).match(a*cos(b)) == {a: x, b: y}

    assert (x+1).match(a+1, [a]) == {a: x}
    assert (x*3).match(a*3, [a]) == {a: x}
    assert (x**3).match(a**3, [a]) == {a: x}
    assert (x*cos(y)).match(a*cos(b), [a,b]) == {a: x, b: y}

    assert (x*y).match(p*q,[p,q]) in [{p:x, q:y}, {p:y, q:x}]
    assert (x+y).match(p+q,[p,q]) in [{p:x, q:y}, {p:y, q:x}]
    assert (x*y+1).match(p*q,[p,q]) in [{p:1, q:1+x*y}, {p:1+x*y, q:1}]

def _test_derivative():
    x,y,a,b,c,p,q = map(Symbol, 'xyabcpq')
    class f(Function): pass
    x = Symbol("x")
    fd = Derivative(f(x),x)
    e = fd+1
    #assert e.match(a+1) == {a: fd}
    assert e.match(a+1, [a]) == {a: fd}
    e = fd
    #assert e.match(a) == {a: fd}
    assert e.match(a, [a]) == {a: fd}
    e = fd
    assert e.match(fd,[a]) == {} # XXX this line fails
    e = 3*fd
    assert e.match(a*fd, [a]) != None
    e = 3*fd - 1
    assert e.match(a*fd + b, [a,b]) == {a:3, b:-1}

def _test_match_deriv_bug1():
    class l(Function): pass
    class n(Function): pass

    r = Symbol("r")

    e = Derivative(l(r),r)/r-Derivative(Derivative(n(r),r),r)/2- \
        Derivative(n(r),r)**2/4+Derivative(n(r),r)*Derivative(l(r),r)/4

    e = e.subs(n(r), -l(r)) # XXX this line fails

    t = r*exp(-l(r))

    t2 = ( diff(t, r,2)/t ).expand()

    a = Symbol("a")
    tt = (a*t2).expand()
    r = e.match(tt, [a])
    assert r == {a: -Rational(1)/2}

    a = Symbol("a__8")
    tt = (a*t2).expand()
    r = e.match(tt, [a])
    assert r == {a: -Rational(1)/2}

def _test_match_bug4():
    r = Symbol("r")
    a = Symbol("a", is_dummy = True)
    e = r
    tt = -a*r # XXX this line fails

    r = e.match(tt, [a])
    assert r == {a: -1}

def _test_match_bug6():
    r = Symbol("r")
    a = Symbol("a", is_dummy = True)
    e = r
    tt = 3*a*r # XXX this line fails

    r = e.match(tt, [a])
    assert r == {a: Rational(1)/3}

def _test_match_bug5():
    r = Symbol("r")
    a = Symbol("a", is_dummy = True)
    e = -r
    tt = -a*r # XXX this line fails

    r = e.match(tt, [a])
    assert r == {a: 1}

def test_behavior1():
    x = Symbol("x")
    a = Symbol("a")
    e = 3*x**2
    #assert e.match(a*x,[a], exclude = None) == {a: 3*x}
    #assert e.match(a*x,[a], exclude = [x]) == None
    #assert e.match(a*x,[a]) == None
    assert e.match(a*x,[a]) == {a: 3*x}

def test_behavior2():
    x = Symbol("x")
    a = Symbol("a")
    e = Rational(6)
    assert e.match(2*a, [a]) == {a: 3}
    e = 3*x**2+3*x+6
    p = x
    #assert (e/p).expand().match(a*x**2+a*x+2*a, [a]) == None
    #assert (e/p).expand().match(a*x**2+a*x+2*a, [a], exclude=None) == {a: 3/x}
    assert (e/p).expand().match(a*x**2+a*x+2*a, [a]) == {a: 3/x}

def test_bug2():
    x,y,a,b,c = map(Symbol, 'xyabc')
    r = (x+y).match(a+b+c, [a,b,c])
    assert (a+b+c).subs_dict(r) == x+y


def test_bug3():
     a = Symbol('a'); b = Symbol('b'); x = Symbol('x')

     assert (b*x*exp(a*x)).match(x*exp(a*x),[a,b]) == None
