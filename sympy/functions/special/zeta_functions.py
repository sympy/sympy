
from sympy.core import Function, S, Basic, pi

###############################################################################
###################### HURWITZ GENERALIZED ZETA FUNCTION ######################
###############################################################################

class zeta(Function):

    nargs = (1, 2)

    @classmethod
    def canonize(cls, z, a=S.One):
        z, a = map(Basic.sympify, (z, a))

        if isinstance(a, Basic.Number):
            if a is S.NaN:
                return S.NaN
            elif a is S.Zero:
                return cls(z)

        if isinstance(z, Basic.Number):
            if z is S.NaN:
                return S.NaN
            elif z is S.Infinity:
                return S.One
            elif z is S.Zero:
                if a.is_negative:
                    return S.Half - a - 1
                else:
                    return S.Half - a
            elif z is S.One:
                return S.ComplexInfinity
            elif isinstance(z, Basic.Integer):
                if isinstance(a, Basic.Integer):
                    if z.is_negative:
                        zeta = (-1)**z * Basic.bernoulli(-z+1)/(-z+1)
                    elif z.is_even:
                        B, F = Basic.bernoulli(z), Basic.Factorial(z)
                        zeta = 2**(z-1) * abs(B) * pi**z / F

                    if a.is_negative:
                        return zeta + Basic.harmonic(abs(a), z)
                    else:
                        return zeta - Basic.harmonic(a-1, z)


class dirichlet_eta(Function):
    """
    Dirichlet eta function
    """
    nargs = 1

    @classmethod
    def canonize(cls, s):
        if s == 1:
            return Basic.log(2)
        else:
            return (1-2**(1-s)) * zeta(s)
