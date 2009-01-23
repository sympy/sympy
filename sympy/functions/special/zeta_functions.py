
from sympy.core import Function, S, C, Basic, sympify, pi
from sympy.utilities.decorator import deprecated

###############################################################################
###################### HURWITZ GENERALIZED ZETA FUNCTION ######################
###############################################################################

class zeta(Function):

    nargs = (1, 2)

    @classmethod
    @deprecated
    def canonize(cls, z, a=S.One):
        return cls.eval(z, a)

    @classmethod
    def eval(cls, z, a=S.One):
        z, a = map(sympify, (z, a))

        if a.is_Number:
            if a is S.NaN:
                return S.NaN
            elif a is S.Zero:
                return cls(z)

        if z.is_Number:
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
            elif z.is_Integer:
                if a.is_Integer:
                    if z.is_negative:
                        zeta = (-1)**z * C.bernoulli(-z+1)/(-z+1)
                    elif z.is_even:
                        B, F = C.bernoulli(z), C.Factorial(z)
                        zeta = 2**(z-1) * abs(B) * pi**z / F
                    else:
                        return

                    if a.is_negative:
                        return zeta + C.harmonic(abs(a), z)
                    else:
                        return zeta - C.harmonic(a-1, z)


class dirichlet_eta(Function):
    """
    Dirichlet eta function
    """
    nargs = 1

    @classmethod
    @deprecated
    def canonize(cls, s):
        return cls.eval(s)

    @classmethod
    def eval(cls, s):
        if s == 1:
            return C.log(2)
        else:
            return (1-2**(1-s)) * zeta(s)
