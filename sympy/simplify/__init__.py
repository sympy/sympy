"""The module helps converting SymPy expressions into shorter forms of them.

for example:
the expression E**(pi*I) will be converted into -1
the expression (x+x)**2 will be converted into 4*x**2
"""
from .simplify import (ratsimp,
    simplify, combsimp, hypersimp, hypersimilar,
    logcombine, separatevars, posify, polarify,
    unpolarify, signsimp, ratsimpmodprime, bottom_up, nsimplify)

from .fu import FU, fu

from .sqrtdenest import sqrtdenest

from .cse_main import cse

from .traversaltools import use

from .epathtools import epath, EPath

from .hyperexpand import hyperexpand

from .radsimp import collect, rcollect, radsimp, collect_const, fraction, numer, denom

from .trigsimp import trigsimp, exptrigsimp

from .besselsimp import besselsimp

from .powsimp import powsimp, powdenest
