from sympy.stats import Normal, ChiSquared, LogNormal
from sympy.stats.crv import ContinuousDistribution
from sympy.stats.rv import Density, RandomSymbol
from sympy.unify.rewrites import rewriterules
from sympy.rules.branch import (multiplex, exhaust, chain, yieldify, debug,
        condition)
from sympy.rules.branch.traverse import top_down
from sympy.rules import rebuild, flatten
from sympy import Symbol, Dummy, factor, log, Abs
import functools

x, y, z, n, m = map(Symbol, 'xyznm')
r = Symbol('r', real=True)

mu, sigma = Symbol('mu', real=True), Symbol('sigma', positive=True)

# Equivalences of random expressions under density. E.g.
# Density(Normal(x, 0, 1)) == Density(StandardNormal(y))
# We add the `Density` in these patterns later
_rv_eqs = (
    (Normal(x, 0, 1)**2, ChiSquared(x, 1), [x], None),
    (ChiSquared(x, m) + ChiSquared(y, n), ChiSquared(x, n+m), [x, y, n, m],
        None),
    (log(Normal(x, mu, sigma)), LogNormal(x, mu, sigma), [x, mu, sigma], None),
    (Normal(x, mu, sigma) + y, Normal(x, mu + y, sigma), [x, y, mu, sigma], None),
    (Normal(x, mu, sigma) * r, Normal(x, mu, sigma*Abs(r)), [x, r, mu, sigma], None),
)

z = Dummy('z')
def additive_eq(src, tgt, wilds, cond):
    """ (X -> Y) -> (X + z -> Y + z) """
    return (src + z, tgt + z, tuple(wilds) + (z,), cond)

rv_eqs = _rv_eqs + tuple(map(additive_eq, *zip(*_rv_eqs)))

def unpack_Density(d):
    if (isinstance(d, Density) and
        isinstance(d.expr, RandomSymbol) and
        isinstance(d.expr.pspace.density, ContinuousDistribution)):
        yield d.expr.pspace.density

from sympy.unify.core import new, is_leaf, children, op
top_down_unify = functools.partial(top_down, new=new, is_leaf=is_leaf,
                                             children=children, op=op)

exprrule = rewriterules(*zip(*rv_eqs),
                    strategy=lambda rules: exhaust(top_down_unify(multiplex(*rules))))

statsimp = condition(lambda x: isinstance(x, Density),
                      exhaust(top_down(multiplex(yieldify(factor),
                                                 exprrule,
                                                 yieldify(rebuild),
                                                 unpack_Density))))
