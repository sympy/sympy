from sympy import (EmptySet, FiniteSet, S, Symbol, Interval, exp, erf, sqrt,
        symbols, simplify, Eq, cos, And, Tuple)
from sympy.stats import (Die, Normal, Exponential , P, E, Var, Covar,
        Skewness, Density, Given, independent, dependent, Where, pspace,
        random_symbols, Sample)
from sympy.stats.rv import ProductPSpace, rs_swap

def test_where():
    X, Y = Die(), Die()
    Z = Normal(0, 1)

    assert Where(Z**2<=1).set == Interval(-1, 1)
    assert Where(Z**2<1).as_boolean() == And(Z.symbol<1, Z.symbol>-1)
    assert len(Where(X<3).set) == 2
    assert 1 in Where(X<3).set

def test_random_symbols():
    X, Y = Normal(0,1), Normal(0,1)

    assert set(random_symbols(2*X+1)) == set((X,))
    assert set(random_symbols(2*X+Y)) == set((X,Y))
    assert set(random_symbols(2*X+Y.symbol)) == set((X,))
    assert set(random_symbols(2)) == set()

def test_pspace():
    X, Y = Normal(0,1), Normal(0,1)

    assert not pspace(5+3)
    assert pspace(X) == X.pspace
    assert pspace(2*X+1) == X.pspace
    #assert pspace(2*X+Y) == ProductPSpace(Y.pspace, X.pspace)

def test_rs_swap():
    x, y = symbols('x y')
    X = Normal(0, 1, symbol=x)
    Y = Exponential(1, symbol=y)

    XX = Normal(0, 2, symbol=x)
    YY = Normal(0, 3, symbol=y)

    expr = 2*X+Y
    assert expr.subs(rs_swap((X,Y), (YY,XX))) == 2*XX+YY

def test_RandomSymbol():

    X = Normal(0, 1, symbol=Symbol('x'))
    Y = Normal(0, 2, symbol=Symbol('x'))
    assert X.symbol == Y.symbol
    assert X!=Y

def test_ProductPSpace():
    X = Normal(0, 1)
    Y = Normal(0, 1)
    px = X.pspace
    py = Y.pspace
    #assert pspace(X+Y) == ProductPSpace(px, py)
    assert pspace(X+Y) == ProductPSpace(py, px)

def test_E():
    assert E(5) == 5

def test_Sample():
    X = Die(6)
    Y = Normal(0,1)

    assert Sample(X) in [1,2,3,4,5,6]
    assert Sample(X+Y).is_Float

    assert P(X+Y>0, Y<0, numsamples=10).is_Rational
    assert E(X+Y, numsamples=10).is_Float
    assert Var(X+Y, numsamples=10).is_Float
