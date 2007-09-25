
from sympy.core import Basic, S, DefinedFunction, Apply, Lambda, symbols

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

        if isinstance(arg, Basic.Number):
            if isinstance(arg, Basic.NaN):
                return S.NaN
            elif isinstance(arg, Basic.Infinity):
                return S.Infinity
            elif isinstance(arg, Basic.Integer):
                if arg.is_positive:
                    return S.Factorial(arg-1)
                else:
                    return S.ComplexInfinity
            elif isinstance(arg, Basic.Rational):
                if arg.q == 2:
                    n = abs(arg.p) / arg.q

                    if arg.is_positive:
                        k, coeff = n, S.One
                    else:
                        n = k = n + 1

                        if n & 1 == 0:
                            coeff = S.One
                        else:
                            coeff = S.NegativeOne

                    for i in range(3, 2*k, 2):
                        coeff *= i

                    if arg.is_positive:
                        return coeff*Basic.sqrt(S.Pi) / 2**n
                    else:
                        return 2**n*Basic.sqrt(S.Pi) / coeff

class ApplyGamma(Apply):

    def _eval_expand_func(self, *args):
        arg = self[0]._eval_expand_basic()

        if isinstance(arg, Basic.Add):
            for i, coeff in enumerate(arg[:]):
                if isinstance(arg[i], Basic.Number):
                    terms = Basic.Add(*(arg[:i] + arg[i+1:]))

                    if isinstance(coeff, Basic.Rational):
                        if coeff.q != 1:
                            terms += Basic.Rational(1, coeff.q)
                            coeff = Basic.Integer(int(coeff))
                    else:
                        continue

                    return S.Gamma(terms)*S.RisingFactorial(terms, coeff)

        return self.func(*self[:])

    def _eval_is_real(self):
        return self[0].is_real

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

class UpperGamma(DefinedFunction):
    """Upper incomplete gamma function"""

    nofargs = 2

    def fdiff(self, argindex=2):
        if argindex == 2:
            a, z = symbols('az', dummy=True)
            return Lambda(-S.Exp(-z)*z**(a-1), a, z)
        else:
            raise ArgumentIndexError(self, argindex)

    def _eval_apply(self, a, z):
        if isinstance(z, Basic.Number):
            if isinstance(z, Basic.NaN):
                return S.NaN
            elif isinstance(z, Basic.Infinity):
                return S.Zero
            elif isinstance(z, Basic.Zero):
                return S.Gamma(a)

        if isinstance(a, Basic.Number):
            if isinstance(a, Basic.One):
                return S.Exp(-z)
            elif isinstance(a, Basic.Integer):
                b = a - 1

                if b.is_positive:
                    return b*self(b, z) + z**b * S.Exp(-z)

class ApplyUpperGamma(Apply):
    pass

Basic.singleton['lowergamma'] = LowerGamma
Basic.singleton['uppergamma'] = UpperGamma

###############################################################################
########################### GAMMA RELATED FUNCTIONS ###########################
###############################################################################

class PolyGamma(DefinedFunction):

    nofargs = 2

    def fdiff(self, argindex=2):
        if argindex == 2:
            n, z = symbols('nz', dummy=True)
            return Lambda(S.PolyGamma(n+1, z), n, z)
        else:
            raise ArgumentIndexError(self, argindex)

    def _eval_apply(self, n, z):
        n, z = map(Basic.sympify, (n, z))

        if n.is_integer:
            if n.is_negative:
                return S.LogGamma(z)
            else:
                if isinstance(z, Basic.Number):
                    if isinstance(z, Basic.NaN):
                        return S.NaN
                    elif isinstance(z, Basic.Infinity):
                        if isinstance(n, Basic.Number):
                            if isinstance(n, Basic.Zero):
                                return S.Infinity
                            else:
                                return S.Zero
                    elif isinstance(z, Basic.Integer):
                        if z.is_nonpositive:
                            return S.ComplexInfinity
                        else:
                            if isinstance(n, Basic.Zero):
                                return -S.EulerGamma + S.Harmonic(z-1, 1)
                            elif n.is_odd:
                                return (-1)**(n+1)*S.Factorial(n)*S.Zeta(n+1, z)

class ApplyPolyGamma(Apply):

    def _eval_expand_func(self, *args):
        n, z = self[0], self[1].expand(func=True)

        if isinstance(n, Basic.Integer) and n.is_nonnegative:
            if isinstance(z, Basic.Add):
                coeff, factors = z.as_coeff_factors()

                if isinstance(coeff, Basic.Integer):
                    tail = Add(*[ z + i for i in xrange(0, int(coeff)) ])
                    return self(n, z-coeff) + (-1)**n*S.Factorial(n)*tail
            elif isinstance(z, Basic.Mul):
                coeff, terms = z.as_coeff_terms()

                if isinstance(coeff, Basic.Integer) and coeff.is_positive:
                    tail = [ self(n, z + i/coeff) for i in xrange(0, int(coeff)) ]

                    if isinstance(n, Basic.Zero):
                        return S.Log(coeff) + Add(*tail)/coeff**(n+1)
                    else:
                        return Add(*tail)/coeff**(n+1)

        return self(n, z)

    def _eval_rewrite_as_zeta(self, n, z):
        return (-1)**(n+1)*S.Factorial(n)*S.Zeta(n+1, z-1)

class LogGamma(DefinedFunction):

    nofargs = 1

    def _eval_apply(self, z):
        return

class ApplyLogGamma(Apply):
    pass

Basic.singleton['polygamma'] = PolyGamma
Basic.singleton['loggamma'] = LogGamma
