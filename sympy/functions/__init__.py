
from sympy.core.basic import Basic

import combinatorial
import elementary
import special

from special.polynomials import legendre, hermite, chebyshevt, chebyshevu, \
        chebyshevu_root, chebyshevt_root

# see #391
from combinatorial.factorials import factorial, rf, ff, binomial
from combinatorial.factorials import Factorial, RisingFactorial, FallingFactorial, Binomial

from elementary.miscellaneous import sqrt, min_, max_
from elementary.complexes import re, im, sign, abs, conjugate, arg
from elementary.trigonometric import acot, cot, tan, cos, sin, asin, acos, atan
from elementary.exponential import exp, log
from elementary.hyperbolic import sinh, cosh, tanh, coth, asinh, acosh, atanh, acoth

from special.error_functions import erf
from special.gamma_functions import gamma, lowergamma, uppergamma, polygamma, loggamma
from special.zeta_functions import dirichlet_eta, zeta

ln = log

for _n, _cls in Basic.singleton.items():
    exec '%s = _cls()' % (_n)
