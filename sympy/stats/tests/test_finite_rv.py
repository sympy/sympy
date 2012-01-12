from sympy import (EmptySet, FiniteSet, S, Symbol, Interval, exp, erf, sqrt,
        symbols, simplify, Eq, cos, And, Tuple, Or, Dict, sympify)
from sympy.stats import (Die, Bernoulli, Coin, P, E, Var, Covar, Sample,
        Density, Given, independent, dependent, Where, FiniteRV, pspace, CDF)
from sympy.utilities.pytest import raises

oo = S.Infinity
def BayesTest(A,B):
    assert P(A, B) == P(And(A, B)) / P(B)
    assert P(A, B) == P(B, A) * P(A) / P(B)

def test_dice():
    X, Y, Z= Die(6), Die(6), Die(6)
    a,b = symbols('a b')

    assert E(X) == 3+S.Half
    assert Var(X) == S(35)/12
    assert E(X+Y) == 7
    assert E(X+X) == 7
    assert E(a*X+b) == a*E(X)+b
    assert Var(X+Y) == Var(X) + Var(Y)
    assert Var(X+X) == 4 * Var(X)
    assert Covar(X,Y) == S.Zero
    assert Covar(X, X+Y) == Var(X)
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

    assert pspace(X).domain.as_boolean() == Or(
            *[Eq(X.symbol, i) for i in [1,2,3,4,5,6]])

    assert Where(X>3).set == FiniteSet(4,5,6)

def test_given():
    X = Die(6)
    Density(X, X>5) == {S(6): S(1)}
    Where(X>2, X>5).as_boolean() == Eq(X.symbol, 6)
    Sample(X, X>5) == 6

def test_domains():
    x, y = symbols('x y')
    X, Y= Die(6, symbol=x), Die(6, symbol=y)
    # Domains
    d = Where(X>Y)
    assert d.condition == (x > y)
    d = Where(And(X>Y, Y>3))
    assert d.as_boolean() == Or(And(Eq(x,5), Eq(y,4)), And(Eq(x,6), Eq(y,5)),
        And(Eq(x,6), Eq(y,4)))
    assert len(d.elements) == 3

    assert len(pspace(X+Y).domain.elements) == 36

    Z = Die(4, symbol=x)

    raises(ValueError, "P(X>Z)") # Two domains with same internal symbol

    pspace(X+Y).domain.set == FiniteSet(1,2,3,4,5,6)**2

    assert Where(X>3).set == FiniteSet(4,5,6)
    assert X.pspace.domain.dict == FiniteSet(
            Dict({X.symbol:i}) for i in range(1,7))

    assert Where(X>Y).dict == FiniteSet(Dict({X.symbol:i, Y.symbol:j})
            for i in range(1,7) for j in range(1,7) if i>j)

def test_dice_bayes():
    X, Y, Z = Die(6), Die(6), Die(6)

    BayesTest(X>3, X+Y<5)
    BayesTest(Eq(X-Y, Z), Z>Y)
    BayesTest(X>3, X>2)

def test_bernoulli():
    p, a, b = symbols('p a b')
    X = Bernoulli(p, a, b, symbol='B')

    assert E(X) == a*p + b*(-p+1)
    assert Density(X)[a] == p
    assert Density(X)[b] == 1-p

    X = Bernoulli(p, 1, 0, symbol='B')

    assert E(X) == p
    assert Var(X) == -p**2 + p
    E(a*X+b) == a*E(X)+b
    Var(a*X+b) == a**2 * Var(X)

def test_dependence():
    X, Y = Die(), Die()
    assert independent(X, 2*Y)
    assert not dependent(X, 2*Y)
    assert dependent(X, Y+X)

    XX, YY = Given(Tuple(X, Y), X+Y>5) # Create a dependency
    assert dependent(XX, YY)

def test_CDF():
    D = Die(6)
    o = S.One

    assert CDF(D) == sympify({1:o/6, 2:o/3, 3:o/2, 4:2*o/3, 5:5*o/6, 6:o})

def test_coins():
    C, D = Coin(), Coin()
    H, T = sorted(Density(C).keys())
    assert P(Eq(C, D)) == S.Half
    assert Density(Tuple(C, D)) == {(H, H): S.One/4, (H, T): S.One/4,
            (T, H): S.One/4, (T, T): S.One/4}
    assert Density(C) == {H: S.Half, T: S.Half}

    E = Coin(S.One/10)
    assert P(Eq(E, H))==S(1)/10

    d = pspace(C).domain

    assert d.as_boolean() == Or(Eq(C.symbol, H), Eq(C.symbol, T))

    raises(ValueError, "P(C>D)") # Can't intelligently compare H to T


def test_FiniteRV():
    F = FiniteRV({1:S.Half, 2:S.One/4, 3:S.One/4})

    assert Density(F) =={S(1):S.Half, S(2):S.One/4, S(3):S.One/4}
    assert P(F>=2)==S.Half

    assert pspace(F).domain.as_boolean() == Or(
            *[Eq(F.symbol, i) for i in [1,2,3]])

