"""The module helps converting SymPy expressions into shorter forms of them.

for example:
the expression E**(pi*I) will be converted into -1
the expression (x+x)**2 will be converted into 4*x**2
"""

__all__ = []

from .simplify import (
    simplify, hypersimp, hypersimilar,
    logcombine, separatevars, posify, besselsimp,
    signsimp, bottom_up, nsimplify
)
__all__ += [
    "simplify", "hypersimp", "hypersimilar",
    "logcombine", "separatevars", "posify", "besselsimp",
    "signsimp", "bottom_up", "nsimplify"
]

from .fu import FU, fu
__all__ += ["FU", "fu"]

from .sqrtdenest import sqrtdenest
__all__ += ["sqrtdenest"]

from .cse_main import cse
__all__ += ["cse"]

from .traversaltools import use
__all__ += ["use"]

from .epathtools import epath, EPath
__all__ += ["epath", "EPath"]

from .hyperexpand import hyperexpand
__all__ += ["hyperexpand"]

from .radsimp import (
    collect, rcollect, radsimp, collect_const,
    fraction, numer, denom
)
__all__ += [
    "collect", "rcollect", "radsimp", "collect_const",
    "fraction", "numer", "denom"
]

from .trigsimp import trigsimp, exptrigsimp
__all__ += ["trigsimp", "exptrigsimp"]

from .powsimp import powsimp, powdenest
__all__ += ["powsimp", "powdenest"]

from .combsimp import combsimp
__all__ += ["combsimp"]

from .gammasimp import gammasimp
__all__ += ["gammasimp"]

from .ratsimp import ratsimp, ratsimpmodprime
__all__ += ["ratsimp", "ratsimpmodprime"]
