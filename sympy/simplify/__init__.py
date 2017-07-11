"""The module helps converting SymPy expressions into shorter forms of them.

for example:
the expression E**(pi*I) will be converted into -1
the expression (x+x)**2 will be converted into 4*x**2
"""
from .combsimp import combsimp
from .cse_main import cse
from .epathtools import EPath, epath
from .fu import FU, fu
from .hyperexpand import hyperexpand
from .powsimp import powdenest, powsimp
from .radsimp import collect, collect_const, denom, fraction, numer, radsimp, \
    rcollect
from .ratsimp import ratsimp, ratsimpmodprime
from .simplify import besselsimp, bottom_up, hypersimilar, hypersimp, \
    logcombine, nsimplify, posify, separatevars, signsimp, simplify
from .sqrtdenest import sqrtdenest
from .traversaltools import use
from .trigsimp import exptrigsimp, trigsimp
