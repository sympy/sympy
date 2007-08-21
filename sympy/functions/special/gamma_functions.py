
from sympy.core.basic import Basic, S, cache_it, cache_it_immutable
from sympy.core.function import DefinedFunction, Apply, Lambda

###############################################################################
############################ COMPLETE GAMMA FUNCTION ##########################
###############################################################################

class Gamma(DefinedFunction):

    nofargs = 1

    def fdiff(self, argindex=1):
        if argindex == 1:
            x = Basic.Symbol('x', dummy=True)
            return Lambda(self(x)*S.PolyGamma(0, x), x)
        else:
            raise ArgumentIndexError(self, argindex)

    def _eval_apply(self, arg):
        arg = Basic.sympify(arg)

        if arg.is_integer:
            if arg.is_positive:
                return S.Factorial(arg-1)
            else:
                return #S.ComplexInfinity
        elif isinstance(arg, Basic.Number):
            if isinstance(arg, Basic.NaN):
                return S.NaN
            elif isinstance(arg, Basic.Infinity):
                return S.Infinity
            elif isinstance(arg, Basic.Rational):
                if arg.q == 2:
                    n = abs(arg.p) // arg.q

                    if arg.is_positive:
                        coeff = S.One
                    else:
                        n = n - 1

                        if n & 1 == 0:
                            coeff = S.One
                        else:
                            coeff = S.NegativeOne

                    for i in range(3, 2*n, 2):
                        coeff *= i

                    if arg.is_positive:
                        return coeff*S.Sqrt(S.Pi)/2**n
                    else:
                        return 2**n*S.Sqrt(S.Pi)/coeff
            elif isinstance(arg, Basic.Real):
                return

class ApplyGamma(Apply):

    def _eval_rewrite_as_rf(self, arg):
        arg = arg.expand(basic=True)

        if isinstance(arg, Basic.Add):
            coeff = terms = S.Zero

            for term in arg:
                if isinstance(term, Basic.Number):
                    coeff += term
                else:
                    terms += term

            if isinstance(coeff, Basic.Integer):
                expo = coeff
            elif isinstance(coeff, Basic.Rational):
                expo = coeff.p / coeff.q
                terms += Rational(coeff.q)
            else:
                return

            return S.Gamma(terms)*S.RisingFactorial(terms, expo)

    def _eval_is_real(self):
        return self.args[0].is_real

Basic.singleton['gamma'] = Gamma

###############################################################################
################## LOWER and UPPER INCOMPLETE GAMMA FUNCTIONS #################
###############################################################################

class LowerGamma(DefinedFunction):
    """Lower incomplete gamma function"""

    nofargs = 2

    def _eval_apply(self, a, x):
        if isinstance(a, Basic.Number):
            if isinstance(a, Basic.One):
                return S.One - S.Exp(-x)
            elif isinstance(a, Basic.Integer):
                b = a - 1

                if b.is_positive:
                    return b*self(b, x) - x**b * S.Exp(-x)

class ApplyLowerGamma(Apply):
    pass

    #def _eval_expand_func(self):
    #    b = self.args[0] - 1
    #
    #    if isinstance(a, Basic.Integer) and b.is_positive:
    #        return (b*LowerGamma(b, x) - x**b * exp(-x)).expand(func=True)

class UpperGamma(DefinedFunction):
    """Upper incomplete gamma function"""

    nofargs = 2

    def _eval_apply(self, a, x):
        if isinstance(x, Basic.Number):
            if isinstance(x, Basic.NaN):
                return S.NaN
            elif isinstance(x, Basic.Infinity):
                return S.Zero
            elif isinstance(x, Basic.Zero):
                return S.Gamma(a)

        if isinstance(a, Basic.Number):
            if isinstance(a, Basic.One):
                return S.Exp(-x)
            elif isinstance(a, Basic.Integer):
                b = a - 1

                if b.is_positive:
                    return b*self(b, x) + x**b * S.Exp(-x)

class ApplyUpperGamma(Apply):
    pass

    #def _eval_expand_func(self):
    #    b = self.args[0] - 1
    #
    #    if isinstance(a, Basic.Integer) and b.is_positive:
    #        return (b*UpperGamma(b, x) + x**b * exp(-x)).expand(func=True)

Basic.singleton['lowergamma'] = LowerGamma
Basic.singleton['uppergamma'] = UpperGamma
