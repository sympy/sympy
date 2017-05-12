from sympy.pm import rubi_match
from sympy.core.symbol import Symbol, Wild

def test_symbols():
    x, a = map(Symbol, 'xa')
    p = Wild('p')

    assert rubi_match(a, p, x) == {p: a}
    assert rubi_match(5, p, x) == {p: 5}


def test_add():
    x, a, b, c = map(Symbol, 'xabc')
    q, r = map(Wild, 'qr')

    assert rubi_match(a,q,x) == {q: a}

    e = a*x + b
    p = q*x + r
    assert rubi_match(e, p, x) == {q: a, r: b}

    e = a*x + b + c
    p = q*x + r
    assert rubi_match(e, p, x) == None

def test_mul():
    x, a, b= map(Symbol, 'xab')
    q = Wild('q')

    e = a*x
    p = q*x
    assert rubi_match(e, p, x) == {q: a}

    e = a*b*x
    p = q*x
    assert rubi_match(e, p, x) == None

def test_function():
    x, a, b, c = map(Symbol, 'xabc')
    q, r = map(Wild, 'qr')
    
