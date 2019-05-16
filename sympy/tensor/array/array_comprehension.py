from __future__ import print_function, division
import functools
from sympy.core.sympify import sympify
from sympy.core.expr import Expr
from sympy.core import Basic

class ArrayComprehension(Basic):
    """
    Generate a list comprehension
    """
    def __new__(cls, expr, *values, **assumptions):
        if any(len(l) != 3 or None for l in values):
            raise ValueError('ArrayComprehension requires values lower and upper bound'
                              'for the expression')
        cls.default_assumptions = assumptions
        obj = Expr.__new__(cls, **assumptions)
        obj.expr = expr
        obj.values = values
        arglist = [expr]
        arglist.extend(values)
        obj._args = tuple(arglist)
        return obj

    def doit(self):
        for index, value in enumerate(self.values):
            var, inf, sup = map(sympify, value)
            # Array will not ne expanded if there is a symbolic dimension
            if any(i.is_symbol for i in [inf, sup]):
                return self
            else:
                self.expr = self._expand_array([inf, sup])


    def _expand_array(self, limits):
        print(limits)
        return 0