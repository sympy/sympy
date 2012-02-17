from sympy import (EmptySet, FiniteSet, S, Symbol, Interval, exp, erf, sqrt,
        symbols, simplify, Eq, cos, And, Tuple, Or, Dict, sympify, binomial)
from sympy.stats import (DiscreteUniform, Die, Bernoulli, Coin, Binomial,
        Hypergeometric, P, E, Var, Covar, Skewness, Sample, Density, Given,
        independent, dependent, Where, FiniteRV, pspace, CDF)
from sympy.utilities.pytest import raises

oo = S.Infinity
def BayesTest(A,B):
    assert P(A, B) == P(And(A, B)) / P(B)
    assert P(A, B) == P(B, A) * P(A) / P(B)

def test_discreteuniform():
    # Symbolic
    a, b, c = symbols('a b c')
    X = DiscreteUniform([a,b,c])

    assert E(X) == (a+b+c)/3
    assert Var(X) == (a**2+b**2+c**2)/3 - (a/3+b/3+c/3)**2
    assert P(Eq(X, a)) == P(Eq(X, b)) == P(Eq(X, c)) == S('1/3')

    Y = DiscreteUniform(range(-5, 5))

    # Numeric
    assert E(Y) == S('-1/2')
    assert Var(Y) == S('33/4')
    assert Skewness(Y) == 0
    for x in range(-5, 5):
        assert P(Eq(Y, x)) == S('1/10')
        assert P(Y <= x) == S(x+6)/10
        assert P(Y >= x) == S(5-x)/10

    assert Density(Die(6)) == Density(DiscreteUniform(range(1,7)))

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

def test_binomial_numeric():
    nvals = range(5)
    pvals = [0, S(1)/4, S.Half, S(3)/4, 1]

    for n in nvals:
        for p in pvals:
            X = Binomial(n, p)
            assert Eq(E(X), n*p)
            assert Eq(Var(X), n*p*(1-p))
            if n > 0 and 0 < p < 1:
                assert Eq(Skewness(X), (1-2*p)/sqrt(n*p*(1-p)))
            for k in range(n+1):
                assert Eq(P(Eq(X, k)), binomial(n, k)*p**k*(1-p)**(n-k))

def test_binomial_symbolic():
    n = 10 # Because we're using for loops, can't do symbolic n
    p = symbols('p', positive=True)
    X = Binomial(n, p)
    assert Eq(simplify(E(X)), n*p)
    assert Eq(simplify(Var(X)), n*p*(1-p))
# Can't detect the equality
#    assert Eq(simplify(Skewness(X)), (1-2*p)/sqrt(n*p*(1-p)))

    # Test ability to change success/failure winnings
    H, T = symbols('H T')
    Y = Binomial(n, p, succ=H, fail=T)
    assert Eq(simplify(E(Y)), simplify(n*(H*p+T*(1-p))))

def test_hypergeometric_numeric():
    for N in range(1, 5):
        for m in range(0, N+1):
            for n in range(1, N+1):
                X = Hypergeometric(N, m, n)
                N, m, n = map(sympify, (N, m, n))
                assert sum(Density(X).values()) == 1
                assert E(X) == n * m / N
                if N > 1:
                    assert Var(X) == n*(m/N)*(N-m)/N*(N-n)/(N-1)
                # Only test for skewness when defined
                if N > 2 and 0 < m < N and n < N:
                    assert Eq(Skewness(X),simplify((N-2*m)*sqrt(N-1)*(N-2*n)
                        / (sqrt(n*m*(N-m)*(N-n))*(N-2))))


def test_FiniteRV():
    F = FiniteRV({1:S.Half, 2:S.One/4, 3:S.One/4})

    assert Density(F) =={S(1):S.Half, S(2):S.One/4, S(3):S.One/4}
    assert P(F>=2)==S.Half

    assert pspace(F).domain.as_boolean() == Or(
            *[Eq(F.symbol, i) for i in [1,2,3]])

