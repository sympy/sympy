
from basic import Basic, S, Memoizer
from numbers import Integer, Rational
from function import DefinedFunction, Lambda, Function, Apply

class IntegerSequence(DefinedFunction):
    """
    http://en.wikipedia.org/wiki/Category:Integer_sequences
    """

    pass

###############################################################################
########################### FACTORIALS and BINOMIAL ###########################
###############################################################################

"""
class Factorial(IntegerSequence):

    nofargs = 1

    def _eval_apply(self, arg):
        arg = Basic.sympify(arg)

        if isinstance(arg, Basic.Number):
            if isinstance(arg, Basic.Zero):
                return S.One
            elif isinstance(arg, Basic.Integer):
                if arg.is_negative:
                    return S.Zero
                else:
                    k, n = arg.p, 1

                    while k:
                        n *= k
                        k -= 1

                    return Integer(n)
"""

class Binomial(IntegerSequence):

    nofargs = 2

    def _eval_apply(self, r, k):
        r, k = map(Basic.sympify, (r, k))

        if isinstance(k, Basic.Integer):
            from sympy.specfun.factorials import factorial
            if k.is_negative:
                return S.Zero
            else:
                rk, result = r - k, []

                for i in xrange(1, k.p+1):
                    result.append(rk+i)

                numer = Basic.Mul(*result)
                denom = factorial(k)

                return (numer/denom).expand()

#Basic.singleton['factorial'] = Factorial
Basic.singleton['binomial'] = Binomial
