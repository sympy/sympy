__version__ = '0.8-svn'

import lib

from mptypes import *
from intervals import mpi

from calculus import diff, diffc, secant, polyval, polyroots, quadts, quadosc
from calculus import sumem, chebyfit, sumsh, sumrich, limit
from calculus import odeint
from calculus import pslq, identify, findpoly

from specfun import lambertw
from specfun import gamma, factorial, lower_gamma, upper_gamma, \
    gammaprod, binomial, rf, ff
from specfun import zeta, bernoulli, bernoulli_range, log_range, jv, jn, j0, j1
from specfun import phi, catalan, euler, khinchin, glaisher, apery
from specfun import hyper, hyp0f1, hyp1f1, hyp2f1, erf, ellipk, ellipe, \
    agm, jacobi, legendre, chebyt, chebyu

