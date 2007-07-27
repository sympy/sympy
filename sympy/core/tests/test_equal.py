from sympy import Symbol, Rational, exp

def testequal():
    b = Symbol("b")
    a = Symbol("a")
    c = Symbol("c")
    e1 = a+b
    e2 = 2*a*b
    e3 = a**3*b**2
    e4 = a*b+b*a
    assert not e1 == e2
    assert not e1 == e2
    assert e1 != e2
    assert e2 == e4
    assert e2 != e3
    assert not e2 == e3

    x = Symbol("x")
    e1 = exp(x+1/x)
    y = Symbol("x")
    e2 = exp(y+1/y)
    assert e1 == e2
    assert not e1 != e2
    y = Symbol("y")
    e2 = exp(y+1/y)
    assert not e1 == e2
    assert e1 != e2

    e5 = Rational(3)+2*x-x-x
    assert e5 == 3
    assert 3 == e5
    assert e5 != 4
    assert 4 != e5
    assert e5 != 3+x
    assert 3+x != e5

def _test_expevalbug(): ### no eval()
    x = Symbol("x")
    e1 = exp(1*x)
    e2 = e1.eval()
    e3 = exp(x)
    assert e1 == e2
    assert e1 == e3
    assert e2 == e3
