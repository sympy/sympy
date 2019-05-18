from __future__ import print_function, division
import functools
from sympy.core.sympify import sympify
from sympy.core.expr import Expr
from sympy.core import Basic
from sympy.core.compatibility import Iterable
from sympy.tensor.array import MutableDenseNDimArray
from sympy import  Symbol

class ArrayComprehension(Basic):
    """
    Generate a list comprehension
    """
    def __new__(cls, expr, *bounds, **assumptions):
        if any(len(l) != 3 or None for l in bounds):
            raise ValueError('ArrayComprehension requires values lower and upper bound'
                              ' for the expression')
        cls.default_assumptions = assumptions
        obj = Expr.__new__(cls, **assumptions)
        obj.expr = expr
        if cls._check_bounds_validity(bounds):
            obj.bounds = bounds
        arglist = [expr]
        arglist.extend(obj.bounds)
        obj._args = tuple(arglist)
        return obj

    @classmethod
    def _check_bounds_validity(cls, bounds):
        return True

    def doit(self):
        expr = self.expr
        for index, value in enumerate(self.bounds):
            var, inf, sup = map(sympify, value)
            # Array will not ne expanded if there is a symbolic dimension
            if Basic(inf, sup).atoms(Symbol):
                return self
        arr = self._expand_array()
        return arr

    # Add/replace a boudary of variable
    def add_bound(self, bound):
        return 0

    def _expand_array(self):
        # To perform a subs at every element of the array.
        def _array_subs(arr, var, val):
            arr = MutableDenseNDimArray(arr)
            for i in range(len(arr)):
                index = arr._get_tuple_index(i)
                arr[index] = arr[index].subs(var, val)
            return arr.tolist()

        # Recursive function to perform subs at every variable according to its boundary
        def f(expr, bounds):
            if len(bounds)== 1:
                var, inf, sup = bounds[0]
                list_gen = [expr.subs(var, val) for val in range(inf, sup+1)]
                return list_gen
            else:
                var, inf, sup = bounds[0]
                list_gen = []
                res = f(expr, bounds[1:])
                for val in range(inf, sup+1):
                    list_gen.append(_array_subs(res, var, val))
                return list_gen

        return f(self.expr, self.bounds)
