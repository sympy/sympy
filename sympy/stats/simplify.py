from sympy.stats import Normal, ChiSquared
from sympy.stats.crv import ContinuousDistribution
from sympy.stats.rv import Density, RandomSymbol
from sympy.unify import rewriterule
from sympy.rules.branch import multiplex, exhaust, chain, yieldify
from sympy.rules import rebuild
from sympy import Symbol

x, y, z, n, m = map(Symbol, 'xyznm')


# Equivalences of random expressions under density. E.g.
# Density(Normal(x, 0, 1)) == Density(StandardNormal(y))
# We add the `Density` in these patterns later
expression_equivalences = [
    (Normal(x, 0, 1)**2, ChiSquared(x, 1), [x]),
    (Normal(x, 0, 1)**2 + z, ChiSquared(x, 1) + z, [x, z]),
    (ChiSquared(x, m) + ChiSquared(y, n), ChiSquared(x, n+m), [x, y, n, m])
]

def unpack_Density(d):
    if (isinstance(d, Density) and
        isinstance(d.expr, RandomSymbol) and
        isinstance(d.expr.pspace.density, ContinuousDistribution)):
        yield d.expr.pspace.density

expression_rrs = [rewriterule(Density(src), Density(tgt), wilds)
                    for src, tgt, wilds in expression_equivalences]

rrs = expression_rrs + [unpack_Density]

statsimp = chain(exhaust(multiplex(*rrs)), yieldify(rebuild))
