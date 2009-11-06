"""A functions module, inculdes all the standart functions.

Combinatorial - factorial, fibonacci, harmonic, bernoulli...
Elementary - hyperbolic, trigonometric, exponential, floor and ceiling, sqrt...
Special - gamma, zeta,spherical harmonics...
"""
from sympy.core.basic import Basic

import combinatorial
import elementary
import special

from special.polynomials import legendre, assoc_legendre, hermite, chebyshevt, chebyshevu, \
        chebyshevu_root, chebyshevt_root

# see #391
from combinatorial.factorials import factorial, rf, ff, binomial
from combinatorial.factorials import Factorial, RisingFactorial, FallingFactorial, Binomial
from combinatorial.numbers import fibonacci, lucas, harmonic, bernoulli, bell

from elementary.miscellaneous import sqrt, min_, max_
from elementary.complexes import re, im, sign, abs, conjugate, arg
from elementary.trigonometric import acot, cot, tan, cos, sin, asin, acos, atan, atan2
from elementary.exponential import exp, log, LambertW
from elementary.hyperbolic import sinh, cosh, tanh, coth, asinh, acosh, atanh, acoth
from elementary.integers import floor, ceiling
from elementary.piecewise import Piecewise, piecewise_fold

from special.error_functions import erf
from special.gamma_functions import gamma, lowergamma, uppergamma, polygamma, loggamma
from special.zeta_functions import dirichlet_eta, zeta
from special.spherical_harmonics import Ylm, Zlm
from special.tensor_functions import Dij, Eijk
from special.delta_functions import DiracDelta, Heaviside

ln = log

for _n, _cls in Basic.singleton.items():
    exec '%s = _cls()' % (_n)
