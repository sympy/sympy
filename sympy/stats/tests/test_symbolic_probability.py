from sympy import symbols, Mul, sin
from sympy.stats import Normal, Poisson, variance
from sympy.stats.rv import probability, expectation
from sympy.stats import Covariance, Variance, Probability, Expectation


def test_literal_probability():
    X = Normal('X', 2, 3)
    Y = Normal('Y', 3, 4)
    Z = Poisson('Z', 4)
    W = Poisson('W', 3)
    x, y, w, z = symbols('x, y, w, z')

    assert Probability(X > 0).doit() == probability(X > 0)
    assert Probability(X > x).doit() == probability(X > x)

    assert Expectation(X).doit() == expectation(X)
    assert Expectation(X**2).doit() == expectation(X**2)
    assert Expectation(x*X) == x*Expectation(X)
    assert Expectation(2*X + 3*Y + z*X*Y) == 2*Expectation(X) + 3*Expectation(Y) + z*Expectation(X*Y)
    assert Expectation(2*X + 3*Y + z*X*Y, evaluate=False).args == (2*X + 3*Y + z*X*Y,)
    assert Expectation(sin(X)) == Expectation(sin(X), evaluate=False)
    assert Expectation(2*x*sin(X)*Y + y*X**2 + z*X*Y) == 2*x*Expectation(sin(X)*Y) + y*Expectation(X**2) + z*Expectation(X*Y)

    assert Variance(w) == 0
    assert Variance(X).doit() == variance(X)
    assert Variance(X + z) == Variance(X)
    assert Variance(X*Y).args == (Mul(X, Y),)
    assert type(Variance(X*Y)) == Variance
    assert Variance(z*X) == z**2*Variance(X)
    assert Variance(X + Y) == Variance(X) + Variance(Y) + 2*Covariance(X, Y)
    assert Variance(X + Y + Z + W) == (Variance(X) + Variance(Y) + Variance(Z) + Variance(W) +
                                       2 * Covariance(X, Y) + 2 * Covariance(X, Z) + 2 * Covariance(X, W) +
                                       2 * Covariance(Y, Z) + 2 * Covariance(Y, W) + 2 * Covariance(W, Z))
    assert Variance(X**2).doit() == variance(X**2)
    assert Variance(X**2, evaluate=False) == Variance(X**2)
    assert Variance(x*X**2) == x**2*Variance(X**2)
    assert Variance(sin(X)).args == (sin(X),)
    assert Variance(sin(X), evaluate=False) == Variance(sin(X))
    assert Variance(x*sin(X)) == x**2*Variance(sin(X))

    assert Covariance(w, z) == 0
    assert Covariance(X, w) == 0
    assert Covariance(w, X) == 0
    assert Covariance(X, Y).args == (X, Y)
    assert type(Covariance(X, Y)) == Covariance
    assert Covariance(z*X + 3, Y) == z*Covariance(X, Y)
    assert Covariance(X, X) == Variance(X)
    assert Covariance(z*X + 3, w*Y + 4) == w*z*Covariance(X,Y)
    assert Covariance(X, Y) == Covariance(Y, X)
    assert Covariance(X + Y, Z + W) == Covariance(W, X) + Covariance(W, Y) + Covariance(X, Z) + Covariance(Y, Z)
    assert Covariance(x*X + y*Y, z*Z + w*W) == (x*w*Covariance(W, X) + w*y*Covariance(W, Y) +
                                                x*z*Covariance(X, Z) + y*z*Covariance(Y, Z))
    assert Covariance(x*X**2 + y*sin(Y), z*Y*Z**2 + w*W) == (w*x*Covariance(W, X**2) + w*y*Covariance(sin(Y), W) +
                                                        x*z*Covariance(Y*Z**2, X**2) + y*z*Covariance(Y*Z**2, sin(Y)))
    assert Covariance(X, X**2) == Covariance(X, X**2, evaluate=False)
    assert Covariance(X, sin(X)) == Covariance(sin(X), X, evaluate=False)
    assert Covariance(X**2, sin(X)*Y) == Covariance(sin(X)*Y, X**2, evaluate=False)


def test_probability_rewrite():
    X = Normal('X', 2, 3)
    Y = Normal('Y', 3, 4)
    Z = Poisson('Z', 4)
    W = Poisson('W', 3)
    x, y, w, z = symbols('x, y, w, z')

    assert Variance(w).rewrite(Expectation) == 0
    assert Variance(X).rewrite(Expectation) == Expectation(X ** 2) - Expectation(X) ** 2
    assert Variance(X, condition=Y).rewrite(Expectation) == Expectation(X ** 2, Y) - Expectation(X, Y) ** 2
    assert Variance(X, Y) != Expectation(X**2) - Expectation(X)**2
    assert Variance(X + z).rewrite(Expectation) == Expectation(X ** 2) - Expectation(X) ** 2
    assert Variance(X * Y).rewrite(Expectation) == Expectation(X ** 2 * Y ** 2) - Expectation(X * Y) ** 2

    assert Covariance(w, X).rewrite(Expectation) == 0
    assert Covariance(X, Y).rewrite(Expectation) == Expectation(X*Y) - Expectation(X)*Expectation(Y)
    assert Covariance(X, Y, condition=W).rewrite(Expectation) == Expectation(X * Y, W) - Expectation(X, W) * Expectation(Y, W)
