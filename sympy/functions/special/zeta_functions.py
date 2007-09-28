
from sympy.core import SingleValuedFunction, S, Basic, pi

###############################################################################
###################### HURWITZ GENERALIZED ZETA FUNCTION ######################
###############################################################################

class zeta(SingleValuedFunction):

    nofargs = (1, 2)

    @classmethod
    def _eval_apply(self, z, a=S.One):
        z, a = map(Basic.sympify, (z, a))

        if isinstance(a, Basic.Number):
            if isinstance(a, Basic.NaN):
                return S.NaN
            elif isinstance(a, Basic.Zero):
                return self(z)

        if isinstance(z, Basic.Number):
            if isinstance(z, Basic.NaN):
                return S.NaN
            elif isinstance(z, Basic.Infinity):
                return S.One
            elif isinstance(z, Basic.Zero):
                if a.is_negative:
                    return S.Half - a - 1
                else:
                    return S.Half - a
            elif isinstance(z, Basic.One):
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


class dirichlet_eta(SingleValuedFunction):
    """
    Dirichlet eta function
    """
    nofargs = 1

    @classmethod
    def _eval_apply(cls, s):
        if s == 1:
            return Basic.log(2)
        else:
            return (1-2**(1-s)) * zeta(s)
