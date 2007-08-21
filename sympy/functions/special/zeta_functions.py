
from sympy.core import *

class Zeta(DefinedFunction):
    """
    Usage
    =====
        zeta(s) -> Riemann zeta function of s

    Notes
    =====
        * The zeta function has a pole at s = 1
        * For even positive integers n, zeta(n) is a rational number
          times a power of pi
        * For nonpositive integers n, zeta(n) is a rational number

    Examples
    ========
        >>> from sympy.specfun.zeta_functions import *
        >>> zeta(4)
        1/90*pi**4
        >>> zeta(5)  # no closed form known
        zeta(5)

    """
    nofargs = 1

    def _eval_apply(self, s):
        if s.is_integer:
            if s == 0:
                return Rational(-1,2)
            if s == 1:
                return oo
            if s > 1 and int(s) % 2 == 0:
                return abs(S.Bernoulli(s)) * 2**(s-1) / S.Factorial(s) * pi**s
            if s < 1:
                return -S.Bernoulli(-s+1)/(-s+1)


class DirichletEta(DefinedFunction):
    """
    Dirichlet eta function
    """
    nofargs = 1

    def _eval_apply(self, s):
        if s == 1:
            return S.Log(2)
        else:
            return (1-2**(1-s)) * S.Zeta(s)




class PolyGamma(DefinedFunction):
    """
    polygamma(m, z) -- mth order polygamma function of z
    """
    nofargs = 2

    #def __repr__(self):
    #    return "polygamma(%r, %r)" % self._args

    def _eval_apply(self, m, z):
        # TODO: rational arguments, reflection formula
        if m.is_integer and m >= 0 and z == 0:
            return S.Infinity
        if m == 0:
            if isinstance(z, Rational):
                #if z < 0:
                #    return polygamma(0, 1-z) - pi/tan(pi*z)
                if z.is_integer and z > 0:
                    return -S.EulerGamma + S.Harmonic(z-1, 1)
        # if m == 1:
        #    if isinstance(z, Rational) and z < 0:
        #        return -polygamma(1, 1-z) + pi**2 / sin(pi*z)**2
        if m.is_integer and m > 0 and z.is_integer and z > 0:
            return (-1)**(m+1)*S.Factorial(m)*(S.Zeta(m+1)-S.Harmonic(z-1, m+1))
        if m.is_integer and z == Rational(1,2):
            return (-1)**(m+1)*S.Factorial(m)*(2**(m+1)-1)*S.Zeta(m+1)

    def fdiff(self, argindex=2):
        if argindex == 2:
            m = Basic.Symbol('m', dummy=True)
            z = Basic.Symbol('z', dummy=True)
            return Lambda(S.PolyGamma(m+1, z), m, z)
        else:
            raise NotImplementedError


def digamma(z):
    return polygamma(0, z)

def trigamma(z):
    return polygamma(1, z)

def tetragamma(z):
    return polygamma(2, z)


Basic.singleton['zeta'] = Zeta
Basic.singleton['dirichlet_eta'] = DirichletEta
Basic.singleton['polygamma'] = PolyGamma
