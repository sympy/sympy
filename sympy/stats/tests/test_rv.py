from sympy import (EmptySet, FiniteSet, S, Symbol, Interval, exp, erf, sqrt,
        symbols, simplify, Eq, cos, And, Tuple, integrate, oo)
from sympy.stats import (Die, Normal, Exponential , P, E, Var, Covar,
        Skewness, Density, Given, independent, dependent, Where, pspace,
        random_symbols, Sample)
from sympy.stats.rv import ProductPSpace, rs_swap
from sympy.utilities.pytest import raises, XFAIL

def test_where():
    X, Y = Die(), Die()
    Z = Normal(0, 1)

    assert Where(Z**2<=1).set == Interval(-1, 1)
    assert Where(Z**2<=1).as_boolean() == Interval(-1,1).as_relational(Z.symbol)
    assert Where(And(X>Y, Y>4)).as_boolean() == And(
            Eq(X.symbol, 6), Eq(Y.symbol, 5))

    assert len(Where(X<3).set) == 2
    assert 1 in Where(X<3).set

    X, Y = Normal(0, 1), Normal(0, 1)
    assert Where(And(X**2 <= 1, X >= 0)).set == Interval(0, 1)
    XX = Given(X, And(X**2 <= 1, X >= 0))
    assert XX.pspace.domain.set == Interval(0, 1)
    assert XX.pspace.domain.as_boolean() == And(0 <= X.symbol, X.symbol**2 <= 1)

    raises(TypeError, "XX = Given(X, X+3)")

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
    assert pspace(2*X+Y) == ProductPSpace(Y.pspace, X.pspace)

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

    assert X.name == X.symbol.name

def test_overlap():
    X = Normal(0, 1, symbol=Symbol('x'))
    Y = Normal(0, 2, symbol=Symbol('x'))

    raises(ValueError, "P(X>Y)")

def test_ProductPSpace():
    X = Normal(0, 1)
    Y = Normal(0, 1)
    px = X.pspace
    py = Y.pspace
    assert pspace(X+Y) == ProductPSpace(px, py)
    assert pspace(X+Y) == ProductPSpace(py, px)

def test_E():
    assert E(5) == 5

def test_Sample():
    X = Die(6)
    Y = Normal(0,1)
    z = Symbol('z')

    assert Sample(X) in [1,2,3,4,5,6]
    assert Sample(X+Y).is_Float

    assert P(X+Y>0, Y<0, numsamples=10).is_Rational
    assert E(X+Y, numsamples=10).is_Float
    assert Var(X+Y, numsamples=10).is_Float

    raises(ValueError, "P(Y>z, numsamples=5)")

def test_Given():
    X = Normal(0, 1)
    Y = Normal(0, 1)
    A = Given(X, True)
    B = Given(X, Y>2)

    assert X == A == B

def test_dependence():
    X, Y = Die(), Die()
    assert independent(X, 2*Y)
    assert not dependent(X, 2*Y)

    X, Y = Normal(0,1), Normal(0,1)
    assert independent(X, Y)
    assert dependent(X, 2*X)

    # Create a dependency
    XX, YY = Given(Tuple(X,Y), Eq(X+Y, 3))
    assert dependent(XX, YY)

@XFAIL
def test_dependent_finite():
    X, Y = Die(), Die()
    # Dependence testing requires symbolic conditions which currently break
    # finite random variables
    assert dependent(X, Y+X)

    XX, YY = Given(Tuple(X, Y), X+Y>5) # Create a dependency
    assert dependent(XX, YY)

def test_normality():
    X, Y = Normal(0,1), Normal(0,1)
    z = Symbol('z', real=True)
    x, density = Density(X-Y, Eq(X+Y, z))

    assert integrate(density, (x, -oo, oo)) == 1
