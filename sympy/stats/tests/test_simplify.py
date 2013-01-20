from sympy.stats.simplify import statsimp
from sympy.stats import Normal, ChiSquared
from sympy.stats.crv_types import ChiSquaredDensity, NormalDensity, Normal
from sympy.stats.rv import Density, pspace
from sympy import Symbol, simplify

from sympy.unify import unify

x, y = map(Symbol, 'xy')

def rebuild(x):
    try:
        return type(x)(*map(rebuild, x.args))
    except:
        return x

def test_chisquared():
    X = Normal('x', 0, 1)
    dens = list(statsimp(Density(X**2)))[0]
    assert isinstance(pspace(dens).density, ChiSquaredDensity)

def test_chisquared_two_degrees():
    X = Normal('X', 0, 1)
    Y = Normal('Y', 0, 1)
    dens = list(statsimp(Density(X**2 + Y**2)))[0].expr.pspace.density

    assert isinstance(dens, ChiSquaredDensity)
    assert rebuild(dens.k) == 2

def _unifies(a, b):
    assert tuple(unify, a, b)

def test_unify():
    assert tuple(unify(NormalDensity(0, 1), NormalDensity(0, 1)))
    assert tuple(unify(Normal(x, 0, 1), Normal(x, 0, 1)))
    assert tuple(unify(Normal(y, 0, 1), Normal(x, 0, 1), wilds=[x]))
    assert list(unify(Normal(x, 0, 1), Normal(y, 0, 1),
                      {}, wilds=[y])) == [{y: x}]
