from sympy.core import Add, S, C, sympify, oo, pi
from sympy.core.function import Function, ArgumentIndexError
from zeta_functions import zeta
from sympy.functions.elementary.exponential import log
from sympy.functions.elementary.integers import floor
from sympy.functions.elementary.miscellaneous import sqrt
from sympy.functions.combinatorial.numbers import bernoulli

###############################################################################
############################ COMPLETE GAMMA FUNCTION ##########################
###############################################################################

class gamma(Function):
    """The gamma function returns a function which passes through the integral
    values of the factorial function, i.e. though defined in the complex plane,
    when n is an integer, gamma(n) = (n - 1)!

    Reference:
        http://en.wikipedia.org/wiki/Gamma_function
    """

    nargs = 1

    def fdiff(self, argindex=1):
        if argindex == 1:
            return gamma(self.args[0])*polygamma(0, self.args[0])
        else:
            raise ArgumentIndexError(self, argindex)

    @classmethod
    def eval(cls, arg):
        if arg.is_Number:
            if arg is S.NaN:
                return S.NaN
            elif arg is S.Infinity:
                return S.Infinity
            elif arg.is_Integer:
                if arg.is_positive:
                    return C.factorial(arg-1)
                else:
                    return S.ComplexInfinity
            elif arg.is_Rational:
                if arg.q == 2:
                    n = abs(arg.p) // arg.q

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


    def _eval_expand_func(self, deep=True, **hints):
        if deep:
            arg = self.args[0].expand(deep, **hints)
        else:
            arg = self.args[0]

        if arg.is_Add:
            coeff, tail = arg.as_coeff_add()
            if coeff and coeff.q != 1:
                tail = (C.Rational(1, coeff.q),) + tail
                coeff = floor(coeff)
            tail = arg._new_rawargs(*tail, **dict(reeval=False))
            return gamma(tail)*C.RisingFactorial(tail, coeff)

        return self.func(*self.args)

    def _eval_is_real(self):
        return self.args[0].is_real

    def _eval_rewrite_as_tractable(self, z):
        return C.exp(loggamma(z))


###############################################################################
################## LOWER and UPPER INCOMPLETE GAMMA FUNCTIONS #################
###############################################################################

class lowergamma(Function):
    """Lower incomplete gamma function"""

    nargs = 2

    @classmethod
    def eval(cls, a, x):
        if a.is_Number:
            if a is S.One:
                return S.One - C.exp(-x)
            elif a.is_Integer:
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
    def eval(cls, a, z):
        if z.is_Number:
            if z is S.NaN:
                return S.NaN
            elif z is S.Infinity:
                return S.Zero
            elif z is S.Zero:
                return gamma(a)

        if a.is_Number:
            if a is S.One:
                return C.exp(-z)
            elif a.is_Integer:
                b = a - 1

                if b.is_positive:
                    return b*cls(b, z) + z**b * C.exp(-z)



###############################################################################
########################### GAMMA RELATED FUNCTIONS ###########################
###############################################################################

class polygamma(Function):
    """The function polygamma(n, z) returns log(gamma(z)).diff(n + 1)

    Reference:
        http://en.wikipedia.org/wiki/Polygamma_function
    """

    nargs = 2

    def fdiff(self, argindex=2):
        if argindex == 2:
            n, z = self.args[0:2]
            return polygamma(n+1, z)
        else:
            raise ArgumentIndexError(self, argindex)

    def _eval_is_positive(self):
        if self.args[1].is_positive and self.args[0] > 0:
            return self.args[0].is_odd

    def _eval_is_negative(self):
        if self.args[1].is_positive and self.args[0] > 0:
            return self.args[0].is_even

    def _eval_is_real(self):
        return self.args[0].is_real

    def _eval_aseries(self, n, args0, x, logx):
        if args0[1] != oo or not \
           (self.args[0].is_Integer and self.args[0].is_nonnegative):
            return super(polygamma, self)._eval_aseries(n, args0, x, logx)
        z = self.args[1]
        N = self.args[0]

        if N == 0:
            # digamma function series
            # Abramowitz & Stegun, p. 259, 6.3.18
            r = log(z) - 1/(2*z)
            o = None
            if n < 2:
                o = C.Order(1/z, x)
            else:
                m = C.ceiling((n+1)/2)
                l = [bernoulli(2*k) / (2*k*z**(2*k)) for k in range(1, m)]
                r -= Add(*l)
                o = C.Order(1/z**(2*m), x)
            return r._eval_nseries(x, n, logx) + o
        else:
            # proper polygamma function
            # Abramowitz & Stegun, p. 260, 6.4.10
            # We return terms to order higher than O(x**n) on purpose
            # -- otherwise we would not be able to return any terms for
            #    quite a long time!
            fac = gamma(N)
            e0 = fac + N*fac/(2*z)
            m = C.ceiling((n+1)/2)
            for k in range(1, m):
                fac = fac*(2*k+N-1)*(2*k+N-2) / ((2*k)*(2*k-1))
                e0 += bernoulli(2*k)*fac/z**(2*k)
            o = C.Order(1/z**(2*m), x)
            if n == 0:
                o = C.Order(1/z, x)
            elif n == 1:
                o = C.Order(1/z**2, x)
            r = e0._eval_nseries(z, n, logx) + o
            return -1 * (-1/z)**N * r

    @classmethod
    def eval(cls, n, z):
        n, z = map(sympify, (n, z))

        if n.is_integer:
            if n.is_negative:
                return loggamma(z)
            else:
                if z.is_Number:
                    if z is S.NaN:
                        return S.NaN
                    elif z is S.Infinity:
                        if n.is_Number:
                            if n is S.Zero:
                                return S.Infinity
                            else:
                                return S.Zero
                    elif z.is_Integer:
                        if z.is_nonpositive:
                            return S.ComplexInfinity
                        else:
                            if n is S.Zero:
                                return -S.EulerGamma + C.harmonic(z-1, 1)
                            elif n.is_odd:
                                return (-1)**(n+1)*C.factorial(n)*zeta(n+1, z)


    def _eval_expand_func(self, deep=True, **hints):
        if deep:
            hints['func'] = False
            n = self.args[0].expand(deep, **hints)
            z = self.args[1].expand(deep, **hints)
        else:
            n, z = self.args[0], self.args[1].expand(deep, func=True)

        if n.is_Integer and n.is_nonnegative:
            if z.is_Add:
                coeff = z.args[0]
                if coeff.is_Integer:
                    e = -(n + 1)
                    if coeff > 0:
                        tail = Add(*[C.Pow(z - i, e)  for i in xrange(1, int(coeff) + 1)])
                    else:
                        tail = -Add(*[C.Pow(z + i, e)  for i in xrange(0, int(-coeff))])
                    return polygamma(n, z - coeff) + (-1)**n*C.factorial(n)*tail

            elif z.is_Mul:
                coeff, z = z.as_two_terms()
                if coeff.is_Integer and coeff.is_positive:
                    tail = [ polygamma(n, z + C.Rational(i, coeff)) for i in xrange(0, int(coeff)) ]
                    if n == 0:
                        return Add(*tail)/coeff + log(coeff)
                    else:
                        return Add(*tail)/coeff**(n+1)

        return polygamma(n, z)

    def _eval_rewrite_as_zeta(self, n, z):
        return (-1)**(n+1)*C.factorial(n)*zeta(n+1, z-1)

class loggamma(Function):

    nargs = 1

    def _eval_aseries(self, n, args0, x, logx):
        if args0[0] != oo:
            return super(loggamma, self)._eval_aseries(n, args0, x, logx)
        z = self.args[0]
        m = min(n, C.ceiling((n+S(1))/2))
        r = log(z)*(z-S(1)/2) - z + log(2*pi)/2
        l = [bernoulli(2*k) / (2*k*(2*k-1)*z**(2*k-1)) for k in range(1, m)]
        o = None
        if m == 0:
            o = C.Order(1, x)
        else:
            o = C.Order(1/z**(2*m-1), x)
        # It is very inefficient to first add the order and then do the nseries
        return (r + Add(*l))._eval_nseries(x, n, logx) + o

    def _eval_rewrite_as_intractable(self, z):
        return log(gamma(z))

    def _eval_is_real(self):
        return self.args[0].is_real

    def fdiff(self, argindex=1):
        if argindex == 1:
            return polygamma(0, self.args[0])
        else:
            raise ArgumentIndexError(self, argindex)

def digamma(x):
    return polygamma(0, x)

def trigamma(x):
    return polygamma(1, x)
