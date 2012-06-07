from sympy import C, pi, I
from sympy.core import Dummy, sympify
from sympy.functions import legendre, assoc_legendre
from sympy.functions.elementary.miscellaneous import sqrt

Pl = legendre
Plm= assoc_legendre

_x = Dummy("x")

def Plmcos(l, m, th):
    """
    Plm(cos(th)).
    """
    l = sympify(l)
    m = sympify(m)
    sin = C.sin
    cos = C.cos
    P = Plm(l, m, _x).subs(_x, cos(th))
    # assume th in (0,pi) => sin(th) is nonnegative
    _sinth = Dummy("_sinth", nonnegative=True)
    P = P.subs(1-cos(th)**2, _sinth**2).subs(_sinth, sin(th))
    return P

def Ylm(l, m, theta, phi):
    """
    Spherical harmonics Ylm.

    Examples
    ========

    >>> from sympy import symbols, Ylm
    >>> theta, phi = symbols("theta phi")
    >>> Ylm(0, 0, theta, phi)
    1/(2*sqrt(pi))
    >>> Ylm(1, -1, theta, phi)
    sqrt(6)*exp(-I*phi)*sin(theta)/(4*sqrt(pi))
    >>> Ylm(1, 0, theta, phi)
    sqrt(3)*cos(theta)/(2*sqrt(pi))

    """
    l, m, theta, phi = [sympify(x) for x in (l, m, theta, phi)]
    factorial = C.factorial
    return sqrt((2*l+1)/(4*pi) * factorial(l-m)/factorial(l+m)) * \
            Plmcos(l, m, theta) * C.exp(I*m*phi)

def Ylm_c(l, m, theta, phi):
    """Conjugate spherical harmonics."""
    return (-1)**m * Ylm(l, -m, theta, phi)

def Zlm(l, m, th, ph):
    """
    Real spherical harmonics.
    """
    from sympy import simplify
    if m > 0:
        zz = C.NegativeOne()**m*(Ylm(l, m, th, ph) + Ylm_c(l, m, th, ph))/sqrt(2)
    elif m == 0:
        return Ylm(l, m, th, ph)
    else:
        zz = C.NegativeOne()**m*(Ylm(l, -m, th, ph) - Ylm_c(l, -m, th, ph))/(I*sqrt(2))

    zz = zz.expand(complex=True)
    zz = simplify(zz)
    return zz
