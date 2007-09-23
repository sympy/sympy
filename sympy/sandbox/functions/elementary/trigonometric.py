from sympy.sandbox.core import Basic, Mul, Integer, Rational, sqrt, pi, I
from sympy.sandbox.core.function import SingleValuedFunction, FunctionSignature

def without(L, x):
    L = L[:]
    L.remove(x)
    return L


class sin(SingleValuedFunction):

    signature = FunctionSignature((Basic,), (Basic,))

    @classmethod
    def canonize(cls, (arg,), **options):
        #if arg == nan:
        #    return nan
        if arg == 0:
            return Integer(0)
        factors = arg.split('*')
        if I in factors:
            # Simplify sin(I*x)
            return I * Basic.sinh(Mul(*without(factors, I)))
        if pi in factors:
            # Simplify sin((p/q)*pi)
            c = Mul(*without(factors, pi))
            if c.is_Rational:
                cases = {1:Integer(0), 2:Integer(1), 3:sqrt(3)/2,
                    4:sqrt(2)/2, 6:Rational(1,2)}
                if c.q in cases:
                    return (-1)**((c.p//c.q)%2) * cases[c.q]
        if any(x.is_Rational and x.p < 0 for x in factors):
            return -sin(-arg)

