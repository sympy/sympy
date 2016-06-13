from __future__ import print_function, division

from sympy.core.function import Function, ArgumentIndexError
from sympy.core import S, sympify, oo, diff
from sympy.core.logic import fuzzy_not
from sympy.core.relational import Eq
from sympy.functions.elementary.complexes import im
from sympy.functions.elementary.piecewise import Piecewise
from sympy.functions.special.delta_functions import DiracDelta, Heaviside

###############################################################################
############################# SINGULARITY FUNCTION ############################
###############################################################################


class SingularityFunction(Function):
    """


    """

    is_real = True

    def fdiff(self, argindex=1):
        if argindex == 1:
            x = sympify(self.args[0])
            a = sympify(self.args[1])
            n = sympify(self.args[2])
            if n == 0 or n == -1:
                return self.func(x, a, n-1)
            elif n > 0:
                return n*self.func(x, a, n-1)
        else:
            raise ArgumentIndexError(self, argindex)

    @classmethod
    def eval(cls, variable, offset, exponent):
        """

        Examples
        ========

        """

        x = sympify(variable)
        a = sympify(offset)
        n = sympify(exponent)
        shift = (x - a)

        if fuzzy_not(im(shift).is_zero):
            raise ValueError("Singularity Functions are defined only for Real Numbers.")
        if fuzzy_not(im(n).is_zero):
            raise ValueError("Singularity Functions are not defined for imaginary exponents.")
        if shift is S.NaN or n is S.NaN:
            return S.NaN
        if (n + 2).is_negative:
            raise ValueError("Singularity Functions are not defined for exponents less than -2.")
        if shift.is_negative:
            return S.Zero
        if not n.is_negative and (shift.is_positive or shift.is_zero):
            return (x - a)**n
        if n == -1 or n == -2:
            if shift.is_negative or shift.is_positive:
                return S.Zero
            if shift.is_zero:
                return S.Infinity

    def _eval_rewrite_as_Piecewise(self, *args):
        """


        """
        x = self.args[0]
        a = self.args[1]
        n = self.args[2]

        if n == -1 or n == -2:
            return Piecewise((oo, Eq((x - a), 0)), (0, True))
        elif n == 0:
            return Piecewise((S(1), (x - a) > 0), (0, True))
        elif n > 0:
            return Piecewise(((x - a)**n, (x - a) > 0), (0, True))

    def _eval_rewrite_as_Heaviside(self, *args):
        """


        """
        x = self.args[0]
        a = self.args[1]
        n = sympify(self.args[2])

        if n == -2:
            return diff(Heaviside(x - a), x, 2)
        if n == -1:
            return diff(Heaviside(x - a), x, 1)
        if n.is_nonnegative:
            return (x - a)**n*Heaviside(x - a)

    _eval_rewrite_as_DiracDelta = _eval_rewrite_as_Heaviside
    _eval_rewrite_as_HeavisideDiracDelta = _eval_rewrite_as_Heaviside
