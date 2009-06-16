__version__ = '0.12'

from mptypes import (
    mpnumeric, mpf, mpc, mpi, mpmathify,
    make_mpf, make_mpc, make_mpi,
    mp, extraprec, extradps, workprec, workdps,
    eps, j, inf, nan, isnan, isinf, isint,
    nstr, nprint, fraction, almosteq,
    arange, linspace, rand, absmin, absmax,
    fsum, fprod, fdot,
    chop
)

from usertools import monitor, timing

from functions import (
    pi, degree, e, ln2, ln10, phi, euler,
    catalan, khinchin, glaisher, apery, mertens, twinprime,
    sqrt, cbrt, exp, ln, log, log10, power,
    cos, sin, tan, cosh, sinh, tanh,
    acos, asin, atan, asinh, acosh, atanh,
    sec, csc, cot, sech, csch, coth,
    asec, acsc, acot, asech, acsch, acoth,
    cospi, sinpi, sinc, sincpi,
    fabs, re, im, conj, floor, ceil,
    nthroot, hypot, modf,
    ldexp, frexp, sign, arg,
    degrees, radians, atan2,
    fib, fibonacci,
    zeta, altzeta, gamma, factorial, fac, fac2, beta,
    psi, psi0, psi1, psi2, psi3,
    polygamma, digamma, trigamma, tetragamma, pentagamma,
    harmonic, bernoulli, bernfrac, stieltjes,
    bernpoly, polylog,
    gammainc, gammaprod, binomial, rf, ff,
    hyper, hyp0f1, hyp1f1, hyp2f1,
    erf, erfc, erfi, erfinv, npdf, ncdf,
    ei, li, ci, si, chi, shi,
    fresnels, fresnelc, airyai, airybi,
    ellipe, ellipk, agm, jacobi, legendre, chebyt, chebyu,
    besselj, j0, j1, besseli, bessely, besselk, hankel1, hankel2,
    lambertw,
    barnesg, superfac, hyperfac,
    loggamma, siegeltheta, siegelz, grampoint,
    zetazero, riemannr, primepi, primepi2, primezeta,
    bell, polyexp,
    expm1
)

from elliptic import jtheta, djtheta, jsn, jcn, jdn

from calculus import richardson, shanks, nsum, nprod
from calculus import diff, diffun, diffs, taylor, pade
from calculus import polyval, polyroots
from calculus import fourier, fourierval
from calculus import sumem, chebyfit, limit
from calculus import odeint

from odes import odefun

from quadrature import quad, quadgl, quadts, quadosc, TanhSinh, GaussLegendre

from identification import pslq, identify, findpoly

from matrices import matrix, eye, diag, zeros, ones, hilbert, randmatrix, \
    norm, mnorm
from linalg import lu_solve, inverse, residual, qr_solve, cholesky_solve, det, \
    cond, lu

from visualization import plot, cplot

from optimization import findroot, multiplicity

# be careful when changing this name, don't use test*!
def runtests():
    """
    Run all mpmath tests and print output.
    """
    import os.path
    from inspect import getsourcefile
    import tests.runtests as tests
    testdir = os.path.dirname(os.path.abspath(getsourcefile(tests)))
    importdir = os.path.abspath(testdir + '/../..')
    tests.testit(importdir, testdir)
