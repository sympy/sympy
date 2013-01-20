from sympy.stats import Normal, ChiSquared
from sympy.stats.rv import Density
from sympy.unify import rewriterule
from sympy.rules.branch import multiplex, exhaust
from sympy import Symbol

x, y, z, n, m = map(Symbol, 'xyznm')
density_equivalences = [
    (Normal(x, 0, 1)**2, ChiSquared(x, 1), [x]),
    (Normal(x, 0, 1)**2 + z, ChiSquared(x, 1) + z, [x, z]),
    (ChiSquared(x, m) + ChiSquared(y, n), ChiSquared(x, n+m), [x, y, n, m])
]

density_rrs = [rewriterule(Density(src), Density(tgt), wilds)
                    for src, tgt, wilds in density_equivalences]

rrs = density_rrs + []

statsimp = exhaust(multiplex(*rrs))
