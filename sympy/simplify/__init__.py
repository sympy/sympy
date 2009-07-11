"""The module helps converting sympy expressions into shroter forms of them.

for example:
the expression E**(pi*I) will be converted into -1
the expression (x+x)**2 will be converted into 4*x**2
"""
from simplify import collect, separate, together, radsimp, ratsimp, fraction, \
    simplify, trigsimp, powsimp, combsimp, hypersimp, hypersimilar, nsimplify, \
    logcombine

from rewrite import cancel, trim, apart

from sqrtdenest import sqrtdenest

from cse_main import cse
