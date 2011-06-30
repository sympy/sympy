from sympy import (EmptySet, FiniteSet, S, Symbol, Interval, exp, erf, sqrt,
        symbols, simplify, Eq, cos, And)
from sympy.statistics import (Die, Bernoulli, Coin, P, E, var, covar, skewness,
        Density)

oo = S.Infinity

def test_dice():
    d1,d2,d3 = Die(6), Die(6), Die(6)
    X,Y,Z = d1.value, d2.value, d3.value
    a,b = symbols('a b')

    assert E(X) == 3+S.Half
    assert var(X) == S(35)/12
    assert E(X+Y) == 7
    assert E(X+X) == 7
    assert E(a*X+b) == a*E(X)+b
    assert var(X+Y) == var(X) + var(Y)
    assert var(X+X) == 4 * var(X)
    assert covar(X,Y) == S.Zero
    assert covar(X, X+Y) == var(X)
    assert Density(Eq(cos(X*S.Pi),1))[True] == S.Half

    assert P(X>3) == S.Half
    assert P(2*X > 6) == S.Half
    assert P(X>Y) == S(5)/12
    assert P(Eq(X,Y)) == P(Eq(X,1))

    assert E(X, X>3) == 5
    assert E(X, Y>3) == E(X)
    assert E(X+Y, Eq(X,Y)) == E(2*X)
    assert E(X+Y-Z, 2*X>Y+1) == S(49)/12

    assert P(X>3, X>3) == S.One
    assert P(X>Y, Eq(Y, 6)) == S.Zero
    assert P(Eq(X+Y, 12)) == S.One/36
    assert P(Eq(X+Y, 12), Eq(X, 6)) == S.One/6

    assert Density(X+Y) == Density(Y+Z) != Density(X+X)
    d = Density(2*X+Y**Z)
    assert d[S(22)] == S.One/108 and d[S(4100)]==S.One/216 and S(3130) not in d

def test_bernoulli():
    p, a, b = symbols('p a b')
    B = Bernoulli(p, a, b, symbol='B')
    X = B.value

    assert E(X) == a*p + b*(-p+1)
    assert Density(X)[a] == p
    assert Density(X)[b] == 1-p

    B = Bernoulli(p, 1, 0, symbol='B')
    X = B.value
    assert E(X) == p
    assert var(X) == -p**2 + p
    E(a*X+b) == a*E(X)+b
    var(a*X+b) == a**2 * var(X)
