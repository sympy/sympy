from sympy.stats.simplify import statsimp, unpack_Density, exprrule
from sympy.stats import Normal, ChiSquared, LogNormal
from sympy.stats.crv_types import (ChiSquaredDistribution, NormalDistribution,
        Normal)
from sympy.stats.rv import Density, pspace
from sympy import Symbol, simplify, log, Dummy
from sympy.rules.branch import exhaust, multiplex, chain, yieldify
from sympy.rules import rebuild

from sympy.unify import unify


x, y = map(Symbol, 'xy')

def test_chisquared():
    X = Normal('x', 0, 1)
    assert next(exprrule(X**2)).pspace.density == ChiSquaredDistribution(1)
    assert next(statsimp(Density(X**2))) == ChiSquaredDistribution(1)

def test_chisquared_two_degrees():
    X = Normal('X', 0, 1)
    Y = Normal('Y', 0, 1)
    assert rebuild(next(exprrule(X**2 + Y**2))).pspace.density == \
            ChiSquaredDistribution(2)
    assert next(statsimp(Density(X**2 + Y**2))) == ChiSquaredDistribution(2)

def test_chisquared_three_degrees():
    X = Normal('X', 0, 1)
    Y = Normal('Y', 0, 1)
    Z = Normal('Z', 0, 1)
    assert next(statsimp(Density(X**2 + Y**2 + Z**2))) == \
            ChiSquaredDistribution(3)

def test_unpack_Density():
    assert next(unpack_Density(Density(Normal('x', 0, 1)))) == \
            NormalDistribution(0, 1)

def test_unify():
    assert tuple(unify(NormalDistribution(0, 1), NormalDistribution(0, 1)))
    assert tuple(unify(Normal(x, 0, 1), Normal(x, 0, 1)))
    assert tuple(unify(Normal(y, 0, 1), Normal(x, 0, 1), variables=[x]))
    assert list(unify(Normal(x, 0, 1), Normal(y, 0, 1),
                      {}, variables=[y])) == [{y: x}]

def test_lognormal():
    mu, sigma = Symbol('mu', real=True), Symbol('sigma', positive=True)
    assert rebuild(next(exprrule(log(Normal('X', mu, sigma))))) == \
            LogNormal('X', mu, sigma)
    assert rebuild(next(exprrule(log(Normal('X', mu, sigma)) + 1))) == \
            LogNormal('X', mu, sigma) + 1

def test_exprrule():
    assert next(exprrule(Normal(y, 0, 1)**2)) == ChiSquared(y, 1)
    assert next(exprrule(2*Normal(y, 0, 1)**2)) == 2*ChiSquared(y, 1)
    assert next(exprrule(Density(Normal(y, 0, 1)**2))) == Density(ChiSquared(y, 1))

def test_normal():
    assert rebuild(next(exprrule(Normal(y, 2, 5) + 2))) == Normal(y, 2+2, 5)
    assert rebuild(next(exprrule(Normal(y, 2, 5) * 3))) == Normal(y, 2, 5*3)

    X = Normal(x, 2, 3)
    assert next(statsimp(Density(X*2 + 1))) == NormalDistribution(3, 6)
