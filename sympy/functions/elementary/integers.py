
from sympy.core.basic import Basic, S, cache_it, cache_it_immutable
from sympy.core.function import DefinedFunction, Apply, Lambda

###############################################################################
######################### FLOOR and CEILING FUNCTIONS #########################
###############################################################################

class Floor(DefinedFunction):
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

    nofargs = 1

    def _eval_apply(self, arg):
        arg = Basic.sympify(arg)

        if arg.is_integer:
            return arg
        elif isinstance(arg, Basic.Number):
            if isinstance(arg, Basic.Infinity):
                return S.Infinity
            elif isinstance(arg, Basic.NegativeInfinity):
                return S.NegativeInfinity
            elif isinstance(arg, Basic.NaN):
                return S.NaN
            elif isinstance(arg, Basic.Integer):
                return arg
            elif isinstance(arg, Basic.Rational):
                return Basic.Integer(arg.p // arg.q)
            elif isinstance(arg, Basic.Real):
                return Basic.Integer(int(arg.floor()))
        elif isinstance(arg, Basic.NumberSymbol):
            return arg.approximation_interval(Basic.Integer)[0]
        elif isinstance(arg, Basic.ImaginaryUnit):
            return S.ImaginaryUnit
        elif isinstance(arg, Basic.Add):
            included, excluded = [], []

            for term in arg:
                coeff = term.as_coefficient(S.ImaginaryUnit)

                if coeff is not None and coeff.is_real:
                    excluded.append(self(coeff)*S.ImaginaryUnit)
                elif term.is_real:
                    if term.is_integer:
                        excluded.append(term)
                    else:
                        included.append(term)
                else:
                    return

            if excluded:
                return self(Basic.Add(*included)) + Basic.Add(*excluded)
        else:
            coeff, terms = arg.as_coeff_terms(S.ImaginaryUnit)

            if not terms and not arg.atoms(type=Basic.Symbol):
                if arg.is_negative:
                    return -S.Ceiling(-arg)
                else:
                    return self(arg.evalf())
            elif terms == [ S.ImaginaryUnit ] and coeff.is_real:
                return self(coeff)*S.ImaginaryUnit

class ApplyFloor(Apply):

    def _eval_is_bounded(self):
        return self.args[0].is_bounded

    def _eval_is_real(self):
        return self.args[0].is_real

    def _eval_is_integer(self):
        return self.args[0].is_real

class Ceiling(DefinedFunction):
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

    nofargs = 1

    def _eval_apply(self, arg):
        arg = Basic.sympify(arg)

        if arg.is_integer:
            return arg
        elif isinstance(arg, Basic.Number):
            if isinstance(arg, Basic.Infinity):
                return S.Infinity
            elif isinstance(arg, Basic.NegativeInfinity):
                return S.NegativeInfinity
            elif isinstance(arg, Basic.NaN):
                return S.NaN
            elif isinstance(arg, Basic.Integer):
                return arg
            elif isinstance(arg, Basic.Rational):
                return Basic.Integer(arg.p // arg.q + 1)
            elif isinstance(arg, Basic.Real):
                return Basic.Integer(int(arg.ceiling()))
        elif isinstance(arg, Basic.NumberSymbol):
            return arg.approximation_interval(Basic.Integer)[1]
        elif isinstance(arg, Basic.ImaginaryUnit):
            return S.ImaginaryUnit
        elif isinstance(arg, Basic.Add):
            included, excluded = [], []

            for term in arg:
                coeff = term.as_coefficient(S.ImaginaryUnit)

                if coeff is not None and coeff.is_real:
                    excluded.append(self(coeff)*S.ImaginaryUnit)
                elif term.is_real:
                    if term.is_integer:
                        excluded.append(term)
                    else:
                        included.append(term)
                else:
                    return

            if excluded:
                return self(Basic.Add(*included)) + Basic.Add(*excluded)
        else:
            coeff, terms = arg.as_coeff_terms(S.ImaginaryUnit)

            if not terms and not arg.atoms(type=Basic.Symbol):
                if arg.is_negative:
                    return -S.Floor(-arg)
                else:
                    return self(arg.evalf())
            elif terms == [ S.ImaginaryUnit ] and coeff.is_real:
                return self(coeff)*S.ImaginaryUnit

class ApplyCeiling(Apply):

    def _eval_is_bounded(self):
        return self.args[0].is_bounded

    def _eval_is_real(self):
        return self.args[0].is_real

    def _eval_is_integer(self):
        return self.args[0].is_real

Basic.singleton['floor'] = Floor
Basic.singleton['ceiling'] = Ceiling
