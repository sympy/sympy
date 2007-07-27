"""
Routines for arbitrary-precision numerical quadrature
"""

from float_ import Float
from constants import *
from functions import *

import math

class Quadrature:
    # Change variables to standard interval
    def transform(self, f, a, b, sa, sb):
        if a == 0 and str(b) == 'oo':
            def g(x):
                return f(x)+f(Float(1)/x)/x**2
            a, b = Float(0), Float(1)
        else:
            a = Float(a)
            b = Float(b)
            g = f
        sa = Float(sa)
        sb = Float(sb)
        if (a, b) == (sa, sb):
            return f
        else:
            # linear change of variables
            C = (b-a)/2
            D = (b+a)/2
            def h(x):
                return C * f(D + C*x)
            return h

def _ts_weight(k, h):
    # x[k] = tanh(pi/2 sinh(kh))
    # w[k] = pi/2 cosh(kh) / cosh^2(pi/2 sinh(kh))
    t = Float(k) * h
    a = exp(t)
    pi = pi_float()
    ar = Float(1) / a
    sinht = (a - ar) / 2
    cosht = (a + ar) / 2
    b = exp((pi * sinht) / 2)
    br = Float(1) / b
    sinhb = (b - br) / 2
    coshb = (b + br) / 2
    return sinhb/coshb, (pi/2)*cosht/coshb**2


class TanhSinh(Quadrature):
    """
    The tanh-sinh quadrature rule (also known as the doubly-exponential
    quadrature rule) is well suited for extremely high precision levels
    (tens or hundreds of digits). Unlike Gaussian quadrature, tanh-sinh
    quadrature usually works well with integrands that have integrable
    singularities at the endpoints.

    Simple usage example (integrating 3*x**2 from 2 to 3):

        >>> t = TanhSinh()
        >>> print t(lambda x: 3*x**2, 2, 3)
        19.0000000000000

    Note that the target function should be a regular Python function
    that takes a Float as input and returns a Float. You can also
    integrate from 0 to infinity (use the string 'oo' or the SymPy
    object oo for the upper limit of integration).

    As another example, we can try to compute pi with high precision by
    integrating the area of a circle quadrant (the curve is given by
    the function y = sqrt(1-x**2)). Due to the vertical derivative of
    this function as x approaches 1, Gaussian quadrature would handle
    it poorly. However, the tanh-sinh quadrature rule handles it easily.
    The following calculation should finish in less than a second:

        >>> Float.setdps(50)
        >>> t = TanhSinh()
        >>> print t(lambda x: 4*sqrt(1-x**2), 0, 1)
        3.1415926535897932384626433832795028841971693993751
        >>> print pi_float()     # verifying the result
        3.1415926535897932384626433832795028841971693993751

    Calculating weights for tanh-sinh quadrature is a relatively
    expensive operation, sometimes more expensive than all evaluations
    of the integrand when the rule is applied. It is therefore useful
    to able to generate a TanhSinh object t and then store it away for
    repeated use.

    You can create a quadrature rule with custom parameters:

      TanhSinh()        create TanhSinh rule appropriate for the
                        current working precision
      TanhSinh(dps)     specify a custom precision level
      TanhSinh(dps, m)  custom precision and level m (advanced option)

    The level m is a small integer between, say, 5 for low precision
    levels and 10 for extremely high precision. If you don't specify m,
    SymPy attempts to guess an appropriate value based on the dps value.

    Like all algorithms for numerical integration, there are functions
    for which the tanh-sinh rule gives the wrong result without any
    warning. In particular, it does not generally work very well if the
    integrand has a discontinuity between the endpoints or is
    oscillatory; such integrals should be broken up into multiple
    intervals. A good way to check the result of an integration is to
    repeat it with dps and m both set slightly higher and see how many
    digits agree between the results.

    This implementation of the tanh-sinh algorithm is based on the
    description given in Borwein, Bailey & Girgensohn, "Experimentation
    in Mathematics - Computational Paths to Discovery", A K Peters,
    2003, pages 312-313. It also described in various documents
    available online.
    """

    def __init__(self, dps=None, m=None):
        # Compute abscissas and weights
        prec = dps or Float.getdps()
        Float.store()
        Float.setdps(prec + 10)
        self.prec = prec
        # This estimate needs serious checking
        if m == None:
            m = int(5 + max(0, math.log(prec/30.0, 2)))
        self.m = m
        self.h = h = Float(1) / 2**m
        eps = Float(10)**(-prec-1)
        self.x = []
        self.w = []
        for k in xrange(20 * 2**m + 1):
            x, w = _ts_weight(k, h)
            self.x.append(x)
            self.w.append(w)
            diff = abs(self.x[-1] - Float(1))
            if diff <= eps:
                break
            if m > 6 and k % 300 == 299:
                Float.store(); Float.setdps(5)
                try:
                    print "progress:", -log(diff, 10) / prec
                except:
                    pass
                Float.revert()
        Float.revert()

    def __call__(self, integrand, a, b, extraprec=10):
        Float.store()
        Float.setdps(self.prec + extraprec)
        f = self.transform(integrand, a, b, -1, 1)
        S = Float(0)
        h = Float(1)
        m = self.m
        res = []
        for k in xrange(1, m+1):
            h = h / 2
            for i in xrange(0, len(self.x), 2**(m-k)):
                if i % (2**(m-k+1)) != 0 or k == 1:
                    if i == 0:
                        S = S + self.w[0]*f(Float(0))
                    else:
                        S = S + (self.w[i])*(f(-self.x[i]) + f(self.x[i]))
            res.append(h*S)
        Float.revert()
        return res[-1]

