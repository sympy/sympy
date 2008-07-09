
from sympy.core.basic import Basic, S, C, sympify
from sympy.core.function import Lambda, Function

###############################################################################
######################### FLOOR and CEILING FUNCTIONS #########################
###############################################################################

class floor(Function):
    """Floor is a univariate function which returns the largest integer
       value not greater than its argument. However this implementaion
       generalizes floor to complex numbers.

       More information can be found in "Concrete mathematics" by Graham,
       pp. 87 or visit http://mathworld.wolfram.com/FloorFunction.html.

       >>> from sympy import *

       >>> floor(17)
       17

       >>> floor(Rational(23, 10))
       2

       >>> floor(2*E)
       5

       >>> floor(-Real(0.567))
       -1

       >>> floor(-I/2)
       -I

    """

    nargs = 1

    @classmethod
    def canonize(cls, arg):
        if arg.is_integer:
            return arg
        elif arg.is_Number:
            if arg is S.Infinity:
                return S.Infinity
            elif arg is S.NegativeInfinity:
                return S.NegativeInfinity
            elif arg is S.NaN:
                return S.NaN
            elif arg.is_Integer:
                return arg
            elif arg.is_Rational:
                return C.Integer(arg.p // arg.q)
            elif arg.is_Real:
                return C.Integer(int(arg.floor()))
        elif arg.is_NumberSymbol:
            return arg.approximation_interval(C.Integer)[0]
        elif arg is S.ImaginaryUnit:
            return S.ImaginaryUnit
        elif arg.is_Add:
            included, excluded = [], []

            for term in arg.args:
                coeff = term.as_coefficient(S.ImaginaryUnit)

                if coeff is not None and coeff.is_real:
                    excluded.append(cls(coeff)*S.ImaginaryUnit)
                elif term.is_real:
                    if term.is_integer:
                        excluded.append(term)
                    else:
                        included.append(term)
                else:
                    return

            if excluded:
                return cls(C.Add(*included)) + C.Add(*excluded)
        else:
            coeff, terms = arg.as_coeff_terms(S.ImaginaryUnit)

            if not terms and not arg.atoms(C.Symbol):
                if arg.is_negative:
                    return -ceiling(-arg)
                else:
                    return cls(arg.evalf())
            elif terms == ( S.ImaginaryUnit, ) and coeff.is_real:
                return cls(coeff)*S.ImaginaryUnit

    def _eval_is_bounded(self):
        return self.args[0].is_bounded

    def _eval_is_real(self):
        return self.args[0].is_real

    def _eval_is_integer(self):
        return self.args[0].is_real

class ceiling(Function):
    """Ceiling is a univariate function which returns the smallest integer
       value not less than its argument. Ceiling function is generalized
       in this implementation to complex numbers.

       More information can be found in "Concrete mathematics" by Graham,
       pp. 87 or visit http://mathworld.wolfram.com/CeilingFunction.html.

       >>> from sympy import *

       >>> ceiling(17)
       17

       >>> ceiling(Rational(23, 10))
       3

       >>> ceiling(2*E)
       6

       >>> ceiling(-Real(0.567))
       0

       >>> ceiling(I/2)
       I

    """

    nargs = 1

    @classmethod
    def canonize(cls, arg):
        if arg.is_integer:
            return arg
        elif arg.is_Number:
            if arg is S.Infinity:
                return S.Infinity
            elif arg is S.NegativeInfinity:
                return S.NegativeInfinity
            elif arg is S.NaN:
                return S.NaN
            elif arg.is_Integer:
                return arg
            elif arg.is_Rational:
                return C.Integer(arg.p // arg.q + 1)
            elif arg.is_Real:
                return C.Integer(int(arg.ceiling()))
        elif arg.is_NumberSymbol:
            return arg.approximation_interval(C.Integer)[1]
        elif arg is S.ImaginaryUnit:
            return S.ImaginaryUnit
        elif arg.is_Add:
            included, excluded = [], []

            for term in arg.args:
                coeff = term.as_coefficient(S.ImaginaryUnit)

                if coeff is not None and coeff.is_real:
                    excluded.append(cls(coeff)*S.ImaginaryUnit)
                elif term.is_real:
                    if term.is_integer:
                        excluded.append(term)
                    else:
                        included.append(term)
                else:
                    return

            if excluded:
                return cls(C.Add(*included)) + C.Add(*excluded)
        else:
            coeff, terms = arg.as_coeff_terms(S.ImaginaryUnit)

            if not terms and not arg.atoms(C.Symbol):
                if arg.is_negative:
                    return -floor(-arg)
                else:
                    return cls(arg.evalf())
            elif terms == ( S.ImaginaryUnit, ) and coeff.is_real:
                return cls(coeff)*S.ImaginaryUnit

    def _eval_is_bounded(self):
        return self.args[0].is_bounded

    def _eval_is_real(self):
        return self.args[0].is_real

    def _eval_is_integer(self):
        return self.args[0].is_real
