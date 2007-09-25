
from sympy.core import *

###############################################################################
###################### HURWITZ GENERALIZED ZETA FUNCTION ######################
###############################################################################

class Zeta(DefinedFunction):

    nofargs = (1, 2)

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
                        zeta = (-1)**z * S.Bernoulli(-z+1)/(-z+1)
                    elif z.is_even:
                        B, F = S.Bernoulli(z), S.Factorial(z)
                        zeta = 2**(z-1) * abs(B) * pi**z / F

                    if a.is_negative:
                        return zeta + S.Harmonic(abs(a), z)
                    else:
                        return zeta - S.Harmonic(a-1, z)

class ApplyZeta(Apply):
    pass

Basic.singleton['zeta'] = Zeta

#####

class DirichletEta(DefinedFunction): # TBD
    """
    Dirichlet eta function
    """
    nofargs = 1

    def _eval_apply(self, s):
        if s == 1:
            return Basic.log(2)
        else:
            return (1-2**(1-s)) * S.Zeta(s)


Basic.singleton['dirichlet_eta'] = DirichletEta

