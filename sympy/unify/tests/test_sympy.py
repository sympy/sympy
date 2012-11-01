from sympy import Add, Basic, Wild
from sympy.unify.unify import Compound, Variable
from sympy.unify.unify_sympy import (destruct, construct, unify, is_associative,
        is_commutative, iswild, wildify, wildtoken, patternify, outermost)
from sympy.abc import w, x, y, z, n, m

def test_destruct():
    expr     = Basic(1, 2, 3)
    expected = Compound(Basic, (1, 2, 3))
    assert destruct(expr) == expected

def test_construct():
    expr     = Compound(Basic, (1, 2, 3))
    expected = Basic(1, 2, 3)
    assert construct(expr) == expected

def test_nested():
    expr = Basic(1, Basic(2), 3)
    cmpd = Compound(Basic, (1, Compound(Basic, (2,)), 3))
    assert destruct(expr) == cmpd
    assert construct(cmpd) == expr

def test_unify():
    expr = Basic(1, 2, 3)
    a, b, c = map(Wild, 'abc')
    pattern = Basic(a, b, c)
    assert list(unify(expr, pattern, {})) == [{a: 1, b: 2, c: 3}]

def test_unify_commutative():
    expr = Add(1, 2, 3, evaluate=False)
    a, b, c = map(Wild, 'abc')
    pattern = Add(a, b, c, evaluate=False)

    assert setdicteq(unify(expr, pattern, {}), ({a: 1, b: 2, c: 3},
                                                {a: 1, b: 3, c: 2},
                                                {a: 2, b: 1, c: 3},
                                                {a: 2, b: 3, c: 1},
                                                {a: 3, b: 1, c: 2},
                                                {a: 3, b: 2, c: 1}))

def setsetstr(a):
    return set(frozenset(str(item) for item in ael.items()) for ael in a)
def setdicteq(a, b):
    return setsetstr(a) == setsetstr(b)

def test_listdictseteq():
    a = [{1:1, 2:2, 3:3}, {2:2, 3:3}]
    b = [{2:2, 3:3}, {1:1, 2:2, 3:3}]
    assert setdicteq(a, b)

def test_unify_iter():
    expr = Add(1, 2, 3, evaluate=False)
    a, b, c = map(Wild, 'abc')
    pattern = Add(a, c, evaluate=False)
    assert is_associative(destruct(pattern))
    assert is_commutative(destruct(pattern))

    result   = list(unify(expr, pattern, {}))
    expected = [{a: 1, c: Add(2, 3, evaluate=False)},
                {a: 2, c: Add(1, 3, evaluate=False)},
                {a: 3, c: Add(1, 2, evaluate=False)},
                {a: Add(1, 2, evaluate=False), c: 3},
                {a: Add(1, 3, evaluate=False), c: 2},
                {a: Add(2, 3, evaluate=False), c: 1}]
    assert setdicteq(result, expected)

def test_hard_match():
    from sympy import sin, cos
    expr = sin(x) + cos(x)**2
    p, q = map(Wild, 'pq')
    pattern = sin(p) + cos(p)**2
    assert list(unify(expr, pattern, {})) == [{p: x}]

def test_wildify():
    assert iswild(wildify(1))
    assert wildtoken(wildify(1)) is 1

def test_patternify():
    assert destruct(patternify(x + y, x)) in (Compound(Add, (Variable(x), y)),
                                              Compound(Add, (y, Variable(x))))
    pattern = patternify(x**2 + y**2, x)
    assert list(unify(pattern, w**2 + y**2, {})) == [{x: w}]

def test_matrix():
    from sympy import MatrixSymbol
    X = MatrixSymbol('X', n, n)
    Y = MatrixSymbol('Y', 2, 2)
    Z = MatrixSymbol('Z', 2, 3)
    p = patternify(X, 'X', n)
    assert list(unify(p, Y, {})) == [{'X': 'Y', n: 2}]
    assert list(unify(p, Z, {})) == []

def test_outermost():
    assert set(outermost(x)) == {x}
    assert set(outermost(x, y)) == {x, y}
    assert set(outermost(x, y, 2*x)) == {2*x, y}
    assert set(outermost(x, y, 2*x + y)) == {2*x + y}
