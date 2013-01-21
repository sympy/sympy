from sympy.stats.simplify import (statsimp, rrs, expression_rrs,
        unpack_Density, rv_eqs)
from sympy.stats import Normal, ChiSquared
from sympy.stats.crv_types import (ChiSquaredDistribution, NormalDistribution,
        Normal)
from sympy.stats.rv import Density, pspace
from sympy import Symbol, simplify
from sympy.rules.branch import exhaust, multiplex, chain, yieldify
from sympy.rules import rebuild

from sympy.unify import unify, rewriterule

expr_rrs = map(rewriterule, *zip(*rv_eqs))
exprsimp = exhaust(multiplex(yieldify(rebuild), *expr_rrs))

x, y = map(Symbol, 'xy')

def test_chisquared():
    X = Normal('x', 0, 1)
    assert next(exprsimp(X**2)).pspace.density == ChiSquaredDistribution(1)
    assert next(statsimp(Density(X**2))) == ChiSquaredDistribution(1)

def test_chisquared_two_degrees():
    X = Normal('X', 0, 1)
    Y = Normal('Y', 0, 1)
    assert next(exprsimp(X**2 + Y**2)).pspace.density == ChiSquaredDistribution(2)
    assert next(statsimp(Density(X**2 + Y**2))) == ChiSquaredDistribution(2)

def test_chisquared_three_degrees():
    X = Normal('X', 0, 1)
    Y = Normal('Y', 0, 1)
    Z = Normal('Z', 0, 1)
    assert next(exprsimp(X**2 + Y**2 + Z**2)).pspace.density ==\
            ChiSquaredDistribution(3)
    assert next(statsimp(Density(X**2 + Y**2 + Z**2))) == \
            ChiSquaredDistribution(3)

def test_unpack_Density():
    assert next(unpack_Density(Density(Normal('x', 0, 1)))) == \
            NormalDistribution(0, 1)

def test_unify():
    assert tuple(unify(NormalDistribution(0, 1), NormalDistribution(0, 1)))
    assert tuple(unify(Normal(x, 0, 1), Normal(x, 0, 1)))
    assert tuple(unify(Normal(y, 0, 1), Normal(x, 0, 1), wilds=[x]))
    assert list(unify(Normal(x, 0, 1), Normal(y, 0, 1),
                      {}, wilds=[y])) == [{y: x}]
