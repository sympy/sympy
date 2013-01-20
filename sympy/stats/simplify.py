from sympy.stats import Normal, ChiSquared
from sympy.stats.rv import Density
from sympy.unify import rewriterule
from sympy.rules.branch import multiplex, exhaust
from sympy import Symbol

x, y, z, n, m = map(Symbol, 'xyznm')
rr_data = [
        (Density(Normal(x, 0, 1)**2), Density(ChiSquared(x, 1)), [x]),
        (Density(Normal(x, 0, 1)**2 + z), Density(ChiSquared(x, 1) + z), [x,z]),
        (Density(ChiSquared(x, m) + ChiSquared(y, n)),
            Density(ChiSquared(x, n+m)), [x, y, n, m])
        ]

rrs = [rewriterule(src, tgt, wilds) for src, tgt, wilds in rr_data]

statsimp = exhaust(multiplex(*rrs))
