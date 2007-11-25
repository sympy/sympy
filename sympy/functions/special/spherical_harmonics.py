from sympy import Function, Basic, Rational, pi, I
from sympy.functions import legendre, assoc_legendre

Pl = legendre
Plm= assoc_legendre


def Plmcos(l, m, th):
    l = Basic.sympify(l)
    m = Basic.sympify(m)
    x = Basic.Symbol("x", dummy = True)
    sin = Basic.sin
    cos = Basic.cos
    P = Plm(l, m, x).subs(x, cos(th))
    P = P.subs(1-cos(th)**2, sin(th)**2)
    return P

def Ylm(l, m, theta, phi):
    l, m, theta, phi = [Basic.sympify(x) for x in (l, m, theta, phi)]
    factorial = Basic.Factorial
    return Basic.sqrt((2*l+1)/(4*pi) * factorial(l-m)/factorial(l+m)) * \
            Plmcos(l, m, theta) * Basic.exp(I*m*phi)

def Ylm_c(l, m, theta, phi):
    """Conjugate spherical harmonics."""
    return (-1)**m * Ylm(l, -m, theta, phi)

def Zlm(l, m, th, ph):
    from sympy import simplify
    sqrt = Basic.sqrt
    if m > 0:
        zz = (-1)**m*(Ylm(l, m, th, ph) + Ylm_c(l, m, th, ph))/sqrt(2)
    elif  m == 0:
        return  Ylm(l, m, th, ph)
    else:
        zz = (-1)**m*(Ylm(l, -m, th, ph) - Ylm_c(l, -m, th, ph))/(I*sqrt(2))

    zz = zz.expand(complex=True)
    zz = simplify(zz)
    return zz
