from __future__ import print_function, division
import functools
from sympy.core.sympify import sympify
from sympy.core.expr import Expr
from sympy.core import Basic
from sympy.core.compatibility import Iterable
from sympy.tensor.array import MutableDenseNDimArray, ImmutableDenseNDimArray
from sympy import  Symbol
from sympy.core.sympify import sympify
from sympy.core.compatibility import Iterable
from sympy.core.numbers import Integer

class ArrayComprehension(Basic):
    """
    Generate a list comprehension
    If there is a symbolic dimension, for example, say [i for i in range(1, N)] where
    N is a Symbol, then the expression will not be expanded to an array. Otherwise,
    calling the doit() function will launch the expansion.

    Examples
    ========

    >>> from sympy.tensor.array import ArrayComprehension
    >>> from sympy.abc import i, j, k
    >>> a = ArrayComprehension(10*i+j, (i, 1, 4), (j, 1, 3))
    >>> a.doit()
    [[11, 12, 13], [21, 22, 23], [31, 32, 33], [41, 42, 43]]
    >>> b = ArrayComprehension(10*i+j, (i, 1, 4), (j, 1, k))
    >>> b.doit()
    ArrayComprehension(10*i + j, (i, 1, 4), (j, 1, k))

    """
    def __new__(cls, expr, *bounds, **assumptions):
        if any(len(l) != 3 or None for l in bounds):
            raise ValueError('ArrayComprehension requires values lower and upper bound'
                              ' for the expression')
        cls.default_assumptions = assumptions
        obj = Expr.__new__(cls, **assumptions)
        obj.expr = expr
        obj.bounds = cls._check_bounds_validity(expr, bounds)
        arglist = [expr]
        arglist.extend(obj.bounds)
        obj._args = tuple(arglist)
        return obj

    @classmethod
    def _check_bounds_validity(cls, expr, bounds):
        bounds = sympify(bounds)
        for var, inf, sup in bounds:
            if var not in expr.free_symbols:
                raise ValueError('Varialbe {} does not exist in expression'.format(var))
            if any(not isinstance(i, (Integer, Symbol)) for i in [inf, sup]):
                raise TypeError('Bounds should be an Integer or a Symbol')
            if isinstance(inf, Integer) and isinstance(sup, Integer) and inf > sup:
                raise ValueError('Lower bound should be inferior to upper bound')
        return bounds

    def doit(self):
        expr = self.expr
        for index, value in enumerate(self.bounds):
            var, inf, sup = map(sympify, value)
            # Array will not ne expanded if there is a symbolic dimension
            if Basic(inf, sup).atoms(Symbol):
                return self
        arr = self._expand_array()
        return arr

    # Substitute the variable with a value, so that the symbolic dimension can be expanded as well
    def subs(self, var, val):
        return 0

    def _expand_array(self):
        # To perform a subs at every element of the array.
        def _array_subs(arr, var, val):
            arr = MutableDenseNDimArray(arr)
            for i in range(len(arr)):
                index = arr._get_tuple_index(i)
                arr[index] = arr[index].subs(var, val)
            return arr.tolist()

        list_gen = self.expr
        for var, inf, sup in reversed(self.bounds):
            list_expr = list_gen
            list_gen = []
            for val in range(inf, sup+1):
                if not isinstance(list_expr, Iterable):
                    list_gen.append(list_expr.subs(var, val))
                else:
                    list_gen.append(_array_subs(list_expr, var, val))
        return ImmutableDenseNDimArray(list_gen)
