__version__ = '0.10'

from mptypes import (\
    mpnumeric, mpf, mpc, mpi, convert_lossless,
    make_mpf, make_mpc, make_mpi,
    mp, extraprec, extradps, workprec, workdps,
    eps, j, inf, nan, isnan, isinf, isint,
    nstr, nprint, fraction, almosteq,
    arange, linspace, rand, absmin, absmax
)

from functions import (\
    pi, degree, e, ln2, ln10, phi, euler,
    catalan, khinchin, glaisher, apery,
    sqrt, cbrt, exp, ln, log, log10, power,
    cos, sin, tan, cosh, sinh, tanh,
    acos, asin, atan, asinh, acosh, atanh,
    sec, csc, cot, sech, csch, coth,
    asec, acsc, acot, asech, acsch, acoth,
    cospi, sinpi,
    floor, ceil,
    nthroot, hypot, modf,
    ldexp, frexp, sign, arg,
    degrees, radians, atan2,
    zeta, gamma, factorial, fac,
    psi, psi0, psi1, psi2, psi3,
    polygamma, digamma, trigamma, tetragamma, pentagamma,
    harmonic, bernoulli, bernfrac, stieltjes,
    lower_gamma, upper_gamma, gammaprod, binomial, rf, ff,
    hyper, hyp0f1, hyp1f1, hyp2f1,
    erf, erfc, erfi, npdf, ncdf,
    ei, li, ci, si, chi, shi,
    fresnels, fresnelc, airyai, airybi,
    ellipe, ellipk, agm, jacobi, legendre, chebyt, chebyu,
    jv, jn, j0, j1,
    lambertw
)

from elliptic import jacobi_theta_1, jacobi_theta_2, jacobi_theta_3, \
    jacobi_theta_4, jacobi_elliptic_sn, jacobi_elliptic_cn, jacobi_elliptic_dn

from calculus import diff, diffc, polyval, polyroots, quadosc
from calculus import sumem, chebyfit, sumsh, sumrich, limit
from calculus import odeint
from quadrature import quad, quadgl, quadts, TanhSinh, GaussLegendre

from identification import pslq, identify, findpoly

from matrices import matrix, eye, diag, zeros, ones, randmatrix, mnorm_1, \
    mnorm_oo, mnorm_F, norm_p
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
