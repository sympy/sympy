"""The module helps converting sympy expressions into shorter forms of them.

for example:
the expression E**(pi*I) will be converted into -1
the expression (x+x)**2 will be converted into 4*x**2
"""
from simplify import (collect, rcollect, separate, radsimp, ratsimp, fraction,
    simplify, trigsimp, powsimp, combsimp, hypersimp, hypersimilar, nsimplify,
    logcombine, separatevars, numer, denom, powdenest, posify, polarify,
    unpolarify, collect_const, signsimp, besselsimp, ratsimpmodprime)

from fu import (TR0, TR1, TR2, TR3, TR4, TR5, TR6, TR7, TR8, TR9, TR10,
    TR10i, TR11, TR12, TR13, CTR1, CTR2, CTR3, CTR4, RL1, RL2, fu)
FU = dict(zip('''TR0, TR1, TR2, TR3, TR4, TR5, TR6, TR7, TR8, TR9, TR10,
    TR10i, TR11, TR12, TR13, CTR1, CTR2, CTR3, CTR4, RL1, RL2'''.split(', '),
    (TR0, TR1, TR2, TR3, TR4, TR5, TR6, TR7, TR8, TR9, TR10,
    TR10i, TR11, TR12, TR13, CTR1, CTR2, CTR3, CTR4, RL1, RL2)))
del (TR0, TR1, TR2, TR3, TR4, TR5, TR6, TR7, TR8, TR9, TR10,
    TR10i, TR11, TR12, TR13, CTR1, CTR2, CTR3, CTR4, RL1, RL2)

from sqrtdenest import sqrtdenest

from cse_main import cse

from traversaltools import use

from epathtools import epath, EPath

from hyperexpand import hyperexpand
