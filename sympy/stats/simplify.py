from sympy.stats import Normal, ChiSquared
from sympy.stats.crv import ContinuousDistribution
from sympy.stats.rv import Density, RandomSymbol
from sympy.unify import rewriterule
from sympy.rules.branch import (multiplex, exhaust, chain, yieldify, debug,
        condition)
from sympy.rules.branch.traverse import top_down
from sympy.rules import rebuild, flatten
from sympy import Symbol, Dummy, factor

x, y, z, n, m = map(Symbol, 'xyznm')

# Equivalences of random expressions under density. E.g.
# Density(Normal(x, 0, 1)) == Density(StandardNormal(y))
# We add the `Density` in these patterns later
_rv_eqs = (
    (Normal(x, 0, 1)**2, ChiSquared(x, 1), [x]),
    (ChiSquared(x, m) + ChiSquared(y, n), ChiSquared(x, n+m), [x, y, n, m]),
)

z = Dummy('z')
def additive_eq(src, tgt, wilds):
    """ (X -> Y) -> (X + z -> Y + z) """
    return (src + z, tgt + z, tuple(wilds) + (z,))

rv_eqs = _rv_eqs + tuple(map(additive_eq, *zip(*_rv_eqs)))

def unpack_Density(d):
    if (isinstance(d, Density) and
        isinstance(d.expr, RandomSymbol) and
        isinstance(d.expr.pspace.density, ContinuousDistribution)):
        yield d.expr.pspace.density

expression_rrs = [rewriterule(src, tgt, wilds)
                    for src, tgt, wilds in rv_eqs]

rrs = expression_rrs + [unpack_Density]

rules = [yieldify(factor)] + rrs + [yieldify(rebuild)]

statsimp = condition(lambda x: isinstance(x, Density),
                     exhaust(top_down(multiplex(*rules))))
