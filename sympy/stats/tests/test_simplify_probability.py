from sympy.external.importtools import import_module

from sympy import exp, symbols, log
from sympy.stats import Normal, LogNormal
from sympy.stats.crv_types import NormalDistribution, LogNormalDistribution
from sympy.stats.rv import RandomSymbol
from sympy.stats.simplify_probability import probability_simplify

mu, sigma = symbols("mu sigma")
mu1, mu2, mu3, mu4 = symbols("mu1:5")
sigma1, sigma2, sigma3, sigma4 = symbols("sigma1:5")

matchpy = import_module("matchpy")


def test_simplify_random_var_expressions():
    if matchpy is None:
        return

    N = Normal("N", mu, sigma)
    result = probability_simplify(exp(N))
    assert isinstance(result, RandomSymbol)
    assert result.name == 'exp(N)'
    assert isinstance(result.pspace.distribution, LogNormalDistribution)

    LN = LogNormal("LN", mu, sigma)
    result = probability_simplify(log(LN))
    assert isinstance(result, RandomSymbol)
    assert result.name == 'log(LN)'
    assert isinstance(result.pspace.distribution, NormalDistribution)

    N1 = Normal("N1", mu1, sigma1)
    N2 = Normal("N2", mu2, sigma2)
    N3 = Normal("N3", mu3, sigma3)
    N4 = Normal("N4", mu4, sigma4)
    result = probability_simplify(N1 + N2)
    assert isinstance(result, RandomSymbol)
    assert result.name == 'RandomSymbol<N1 + N2>'
    assert isinstance(result.pspace.distribution, NormalDistribution)

    result = probability_simplify(N1 + N2 + N3)
    assert isinstance(result, RandomSymbol)
    assert result.name == 'RandomSymbol<N1 + N2 + N3>'
    assert isinstance(result.pspace.distribution, NormalDistribution)

    result = probability_simplify(N1 + N2 + N3 + N4)
    assert isinstance(result, RandomSymbol)
    assert result.name == 'RandomSymbol<N1 + N2 + N3 + N4>'
    assert isinstance(result.pspace.distribution, NormalDistribution)
