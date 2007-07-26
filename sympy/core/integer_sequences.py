
from basic import Basic, Atom
from function import DefinedFunction, Lambda, Function, Apply

class IntegerSequence(DefinedFunction):
    """
    http://en.wikipedia.org/wiki/Category:Integer_sequences
    """

    pass

#class Factorial(IntegerSequence):

#    nofargs = 1

#    def _eval_apply(self, arg):
#        if isinstance(arg, Basic.Zero): return Basic.One()
#        if isinstance(arg, Basic.Integer) and arg.is_positive:
#            r = arg.p
#            m = 1
#            while r:
#                m *= r
#                r -= 1
#            return Basic.Integer(m)
        # return Gamma(arg)

class Binomial(IntegerSequence):
    """
    http://en.wikipedia.org/wiki/Binomial_coefficient
    """

    nofargs = 2

    def _eval_apply(self, z, k):
        from sympy.modules.specfun.factorials import Factorial
        if isinstance(k, Basic.Integer):
            assert not k.is_negative,`k`
            l = []
            zk = z - k
            for n in xrange(1, k.p+1):
                l.append(zk+n)
            return (Basic.Mul(*l)/Factorial()(k)).expand()

#Basic.singleton['factorial'] = Factorial
Basic.singleton['binomial'] = Binomial

