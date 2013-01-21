from sympy.stats.simplify import (statsimp, rrs, expression_rrs,
        unpack_Density, rebuild, yieldify)
from sympy.stats import Normal, ChiSquared
from sympy.stats.crv_types import (ChiSquaredDistribution, NormalDistribution,
        Normal)
from sympy.stats.rv import Density, pspace
from sympy import Symbol, simplify
from sympy.rules.branch import exhaust, multiplex, chain

from sympy.unify import unify

exprstatsimp = chain(exhaust(multiplex(*expression_rrs)), yieldify(rebuild))

x, y = map(Symbol, 'xy')

def test_chisquared():
    X = Normal('x', 0, 1)
    dens = list(exprstatsimp(Density(X**2)))[0]
    assert isinstance(pspace(dens).density, ChiSquaredDistribution)

def test_chisquared_two_degrees():
    X = Normal('X', 0, 1)
    Y = Normal('Y', 0, 1)
    dens = list(exprstatsimp(Density(X**2 + Y**2)))[0].expr.pspace.density
    assert dens == ChiSquaredDistribution(2)

def test_unpack_Density():
    assert next(unpack_Density(Density(Normal('x', 0, 1)))) == \
            NormalDistribution(0, 1)

def test_statsimp():
    expr = Normal('X', 0, 1)**2 + Normal('Y', 0, 1)**2
    assert set(statsimp(Density(expr))) == set([ChiSquaredDistribution(2)])

def test_unify():
    assert tuple(unify(NormalDistribution(0, 1), NormalDistribution(0, 1)))
    assert tuple(unify(Normal(x, 0, 1), Normal(x, 0, 1)))
    assert tuple(unify(Normal(y, 0, 1), Normal(x, 0, 1), wilds=[x]))
    assert list(unify(Normal(x, 0, 1), Normal(y, 0, 1),
                      {}, wilds=[y])) == [{y: x}]
