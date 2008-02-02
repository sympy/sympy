
from sympy.core import Basic, Add, S, C, sympify
from sympy.core.function import Function
from zeta_functions import zeta
from sympy.functions.elementary.exponential import log
from sympy.functions.elementary.miscellaneous import sqrt

###############################################################################
############################ COMPLETE GAMMA FUNCTION ##########################
###############################################################################

class gamma(Function):

    nargs = 1

    def fdiff(self, argindex=1):
        if argindex == 1:
            return gamma(self.args[0])*polygamma(0, self.args[0])
        else:
            raise ArgumentIndexError(self, argindex)

    @classmethod
    def canonize(cls, arg):
        arg = sympify(arg)

        if isinstance(arg, C.Number):
            if arg is S.NaN:
                return S.NaN
            elif arg is S.Infinity:
                return S.Infinity
            elif isinstance(arg, C.Integer):
                if arg.is_positive:
                    return C.Factorial(arg-1)
                else:
                    return S.ComplexInfinity
            elif isinstance(arg, C.Rational):
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
                        return coeff*sqrt(S.Pi) / 2**n
                    else:
                        return 2**n*sqrt(S.Pi) / coeff


    def _eval_expand_func(self, *args):
        arg = self.args[0]._eval_expand_basic()

        if isinstance(arg, C.Add):
            for i, coeff in enumerate(arg.args[:]):
                if isinstance(arg.args[i], C.Number):
                    terms = C.Add(*(arg.args[:i] + arg.args[i+1:]))

                    if isinstance(coeff, C.Rational):
                        if coeff.q != 1:
                            terms += C.Rational(1, coeff.q)
                            coeff = C.Integer(int(coeff))
                    else:
                        continue

                    return gamma(terms)*C.RisingFactorial(terms, coeff)

        return self.func(*self.args)

    def _eval_is_real(self):
        return self.args[0].is_real


###############################################################################
################## LOWER and UPPER INCOMPLETE GAMMA FUNCTIONS #################
###############################################################################

class lowergamma(Function):
    """Lower incomplete gamma function"""

    nargs = 2

    @classmethod
    def canonize(cls, a, x):
        if isinstance(a, C.Number):
            if a is S.One:
                return S.One - C.exp(-x)
            elif isinstance(a, C.Integer):
                b = a - 1

                if b.is_positive:
                    return b*cls(b, x) - x**b * C.exp(-x)


class uppergamma(Function):
    """Upper incomplete gamma function"""

    nargs = 2

    def fdiff(self, argindex=2):
        if argindex == 2:
            a, z = self[0:2]
            return -C.exp(-z)*z**(a-1)
        else:
            raise ArgumentIndexError(self, argindex)

    @classmethod
    def canonize(cls, a, z):
        if isinstance(z, C.Number):
            if z is S.NaN:
                return S.NaN
            elif z is S.Infinity:
                return S.Zero
            elif z is S.Zero:
                return gamma(a)

        if isinstance(a, C.Number):
            if a is S.One:
                return C.exp(-z)
            elif isinstance(a, C.Integer):
                b = a - 1

                if b.is_positive:
                    return b*cls(b, z) + z**b * C.exp(-z)



###############################################################################
########################### GAMMA RELATED FUNCTIONS ###########################
###############################################################################

class polygamma(Function):

    nargs = 2

    def fdiff(self, argindex=2):
        if argindex == 2:
            n, z = self.args[0:2]
            return polygamma(n+1, z)
        else:
            raise ArgumentIndexError(self, argindex)

    @classmethod
    def canonize(cls, n, z):
        n, z = map(sympify, (n, z))

        if n.is_integer:
            if n.is_negative:
                return loggamma(z)
            else:
                if isinstance(z, C.Number):
                    if z is S.NaN:
                        return S.NaN
                    elif z is S.Infinity:
                        if isinstance(n, C.Number):
                            if n is S.Zero:
                                return S.Infinity
                            else:
                                return S.Zero
                    elif isinstance(z, C.Integer):
                        if z.is_nonpositive:
                            return S.ComplexInfinity
                        else:
                            if n is S.Zero:
                                return -S.EulerGamma + C.harmonic(z-1, 1)
                            elif n.is_odd:
                                return (-1)**(n+1)*C.Factorial(n)*zeta(n+1, z)


    def _eval_expand_func(self, *args):
        n, z = self.args[0], self.args[1].expand(func=True)

        if isinstance(n, C.Integer) and n.is_nonnegative:
            if isinstance(z, C.Add):
                coeff, factors = z.as_coeff_factors()

                if isinstance(coeff, C.Integer):
                    tail = Add(*[ z + i for i in xrange(0, int(coeff)) ])
                    return polygamma(n, z-coeff) + (-1)**n*C.Factorial(n)*tail
            elif isinstance(z, C.Mul):
                coeff, terms = z.as_coeff_terms()

                if isinstance(coeff, C.Integer) and coeff.is_positive:
                    tail = [ polygamma(n, z + i/coeff) for i in xrange(0, int(coeff)) ]

                    if n is S.Zero:
                        return log(coeff) + Add(*tail)/coeff**(n+1)
                    else:
                        return Add(*tail)/coeff**(n+1)

        return polygamma(n, z)

    def _eval_rewrite_as_zeta(self, n, z):
        return (-1)**(n+1)*C.Factorial(n)*zeta(n+1, z-1)

class loggamma(Function):

    nargs = 1
