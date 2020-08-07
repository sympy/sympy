from sympy import (S, symbols, sqrt, Abs, exp, pi, oo, erf, erfc, Integral,
                Interval, Dummy, factorial, binomial, log, beta, lerchphi,
                Piecewise, lowergamma, gamma)
from sympy.stats import P, E, density, sample, cdf
from sympy.stats.crv_types import NormalDistribution, LaplaceDistribution
from sympy.stats.drv_types import (NegativeBinomialDistribution, PoissonDistribution,
                                    GeometricDistribution, LogarithmicDistribution,
                                    YuleSimonDistribution)
from sympy.stats.frv_types import (BinomialDistribution, HypergeometricDistribution,
                                BetaBinomialDistribution, Binomial)
from sympy.stats.mixture_rv import MixtureDistribution, MixturePSpace, Mixture
from sympy.testing.pytest import raises, skip, ignore_warnings
from sympy.external import import_module

y, z = symbols('y z')
def test_continuous_mixture():
    N = NormalDistribution(0, 1)
    M = NormalDistribution(1, 2)
    Z = LaplaceDistribution(3, 1)
    D = Mixture('D', [2, 5, 3], [N, M, Z])
    assert D.pspace.distribution.is_Continuous
    assert isinstance(D.pspace.distribution, MixtureDistribution)
    assert isinstance(D.pspace, MixturePSpace)
    assert D.pspace.set == Interval(-oo, oo)
    assert density(D)(z) == 3*exp(-Abs(z - 3))/20 + sqrt(2)*exp(-(z - 1)**2/8
                        )/(8*sqrt(pi)) + sqrt(2)*exp(-z**2/2)/(10*sqrt(pi))
    assert E(D).simplify() == S(7)/5
    k = Dummy('k')
    cdf_expr = erf(sqrt(2)*(z - 1)/4)/4 - erfc(sqrt(2)*z/2)/10 + 3*Integral(
                exp(-Abs(k - 3)), (k, -oo, z))/20 + S(9)/20
    assert cdf(D)(z).simplify().dummy_eq(cdf_expr)
    prob = Integral(3*exp(-Abs(k - 3))/20 + sqrt(2)*exp(-(k - 1)**2/8)/(
        8*sqrt(pi)) + sqrt(2)*exp(-k**2/2)/(10*sqrt(pi)), (k, 0, oo))

    with ignore_warnings(UserWarning): ### TODO: Restore tests once warnings are removed
        assert P(D > 0, evaluate=False).rewrite(Integral).dummy_eq(prob)
    scipy = import_module('scipy')
    if not scipy:
        skip('Scipy is not installed. Abort tests')
    else:
        with ignore_warnings(UserWarning): ### TODO: Restore tests once warnings are removed
            assert next(sample(D)) in D.pspace.set
            samp = next(sample(D, size=5))
            for sam in samp:
                assert sam in D.pspace.set

def test_discrete_mixture():
    X = NegativeBinomialDistribution(5, S(1)/3)
    Y = PoissonDistribution(2)
    D = Mixture('D', [2, 3], [X, Y])
    assert D.pspace.distribution.is_Discrete
    assert D.pspace.set == S.Naturals0
    assert density(D)(z).simplify() == 3*2**z*exp(-2)/(5*factorial(z)) + S(64)*3**(-z)*binomial(z + 4,
                                z)/1215
    assert E(D) == S(11)/5
    assert cdf(D)(z).simplify() == Piecewise((1 - 3*lowergamma(z + 1, 2)/(5*gamma(z + 1)
            ) - S(4)*3**(-z)*z**4/3645 - S(64)*3**(-z)*z**3/3645 - S(392)*3**(-z)*z**2/3645 -\
            S(1112)*3**(-z)*z/3645 - S(422)*3**(-z)/1215, z >= 0), (0, True))
    assert P(D > 2).simplify() == S(2813)/3645 - 3*exp(-2)
    X = GeometricDistribution(S(2)/5)
    Y = LogarithmicDistribution(S(2)/3)
    Z = YuleSimonDistribution(3)
    B = Binomial('B', 2, S(1)/10)
    wts = list(density(B).dict.values())
    D = Mixture('D', wts, [X, Y, Z])
    assert density(D)(z).simplify() == S(9)*10**z*15**(-z)/(50*z*log(3)
                        ) + S(3)*beta(z, 4)/100 + S(27)*15**(-z)*3**(2*z)/50
    assert cdf(D)(z).simplify() == Piecewise(((-12*2**z*3**(-z)*lerchphi(S(2)/3, 1, z + 1
            ) + (-81*3**z*5**(-z) - z*beta(z, 4) + 82)*log(3) + log(387420489))/(100*log(3)),
            z >= 1), (0, True))
    assert E(D).simplify() == S(9)/(25*log(3)) + S(51)/25
    scipy = import_module('scipy')
    if not scipy:
        skip('Scipy is not installed. Abort tests')
    else:
        with ignore_warnings(UserWarning): ### TODO: Restore tests once warnings are removed
            assert next(sample(D)) in D.pspace.set
            samp = next(sample(D, size=5))
            for sam in samp:
                assert sam in D.pspace.set

def test_finite_mixture():
    X = HypergeometricDistribution(10, 5, 5)
    Y = BetaBinomialDistribution(5, 1, 2)
    Z = BinomialDistribution(5, S(2)/3, 1, 0)
    B = Binomial('B', 2, S(1)/3)
    wts = list(density(B).dict.values())
    D = Mixture('D', wts, [X, Y, Z])
    assert D.pspace.distribution.is_Finite
    assert D.pspace.domain.set == {0, 1, 2, 3, 4, 5}
    assert sum(density(D).dict.values()).evalf() == 1
    assert E(D) == S(40)*beta(2, 6)/9 + S(40)*beta(6, 2)/9 + S(160)*beta(3, 5
                    )/9 + S(160)*beta(5, 3)/9 + S(80)*beta(4, 4)/3 + S(40)/27
    assert P(D < 5) == S(40)*beta(5, 3)/9 + S(80)*beta(4, 4)/9 + S(80)*beta(3, 5
                    )/9 + S(40)*beta(2, 6)/9 + S(10198)/15309
    scipy = import_module('scipy')
    if not scipy:
        skip('Scipy is not installed. Abort tests')
    else:
        with ignore_warnings(UserWarning): ### TODO: Restore tests once warnings are removed
            assert next(sample(D)) in D.pspace.set
            samp = next(sample(D, size=5))
            for sam in samp:
                assert sam in D.pspace.set

def test_mixture_raises():
    X = HypergeometricDistribution(10, 5, 5)
    Y = BetaBinomialDistribution(5, 1, 2)
    Z = BinomialDistribution(5, S(2)/3, 1, 0)
    raises(TypeError, lambda: Mixture('D', [1, 2], [X, y]))
    Z = BinomialDistribution(2, S(2)/3, 1, 0)
    raises(ValueError, lambda: Mixture('D', [1, 2, 3], [X, Y, Z]))
    raises(ValueError, lambda: Mixture('D', [-2, 3], [X, Y]))
    raises(ValueError, lambda: Mixture('D', [2, 2, 1], [X, Y]))
    N = NormalDistribution(0, 1)
    raises(TypeError, lambda: Mixture('D', [1, 3], [N, Y]))
