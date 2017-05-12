from sympy.pm import rubi_match
from sympy.core.symbol import Symbol, Wild
from sympy.functions import log, sin

def test_symbols():
    x, a = map(Symbol, 'xa')
    p = Wild('p')

    assert rubi_match(a, p, x) == {p: a}
    assert rubi_match(5, p, x) == {p: 5}


def test_mul():
    x, a, b= map(Symbol, 'xab')
    q = Wild('q')

    e = a*x
    p = q*x
    assert rubi_match(e, p, x) == {q: a}

    e = a*b*x
    p = q*x
    assert rubi_match(e, p, x) == None


def test_add():
    x, a, b, c = map(Symbol, 'xabc')
    q, r = map(Wild, 'qr')

    assert rubi_match(a, q, x) == {q: a}

    e = a*x + b
    p = q*x + r
    assert rubi_match(e, p, x) == {q: a, r: b}

    e = a*x + b + c
    p = q*x + r
    assert rubi_match(e, p, x) == None

    e = a*x + b*x
    p = q*x + r*x
    assert rubi_match(e, p, x) == {q: a, r: b}


def test_pow():
    x, a, b, c, d = map(Symbol, 'xabcd')
    q, r, s, t = map(Wild, 'qrst')

    e = x**a
    p = x**q
    assert rubi_match(e, p, x) == {q: a}

    e = (a*x + b)**c
    p = (q*x + r)**s
    assert rubi_match(e, p, x) == {q: a, r: b, s: c}

    e = (a*x + b)**(c*x + d)
    p = (q*x + r)**(s*x + t)
    assert rubi_match(e, p, x) == {q: a, r: b, s: c, t:d}


def test_function():
    x, a, b, c = map(Symbol, 'xabc')
    q, r, s = map(Wild, 'qrs')

    e = log(x*a) + c
    p = log(x*q) + r
    assert rubi_match(e, p, x) == {q: a, r: c}

    e = a*sin(b*x + c)
    p = q*sin(r*x + s)
    assert rubi_match(e, p, x) == {q: a, r: b, s: c}
