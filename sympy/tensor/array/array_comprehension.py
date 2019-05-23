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
    >>> a = ArrayComprehension(10*i + j, (i, 1, 4), (j, 1, 3))
    >>> a
    ArrayComprehension(10*i + j, (i, 1, 4), (j, 1, 3))
    >>> a.doit()
    [[11, 12, 13], [21, 22, 23], [31, 32, 33], [41, 42, 43]]
    >>> b = ArrayComprehension(10*i + j, (i, 1, 4), (j, 1, k))
    >>> b.doit()
    ArrayComprehension(10*i + j, (i, 1, 4), (j, 1, k))

    """
    def __new__(cls, expr, *bounds, **assumptions):
        if any(len(l) != 3 or None for l in bounds):
            raise ValueError('ArrayComprehension requires values lower and upper bound'
                              ' for the expression')
        cls.default_assumptions = assumptions
        arglist = [sympify(expr)]
        arglist.extend(cls._check_bounds_validity(expr, bounds))
        obj = Basic.__new__(cls, *arglist, **assumptions)
        obj._expr = arglist[0]
        obj._bounds = arglist[1:]
        return obj

    @property
    def expr(self):
        """
        Return the expression that will be expanded to an array

        Examples
        ========

        >>> from sympy.tensor.array import ArrayComprehension
        >>> from sympy.abc import i, j, k
        >>> a = ArrayComprehension(10*i + j, (i, 1, 4), (j, 1, 3))
        >>> a.expr
        10*i + j
        """
        return self._args[0]

    @property
    def bounds(self):
        """
        Return a list of the bounds that will be applied while expanding the array. Each
        bound contrains firstly the an variable (not necessarily to be component of the
        expression, e.g. an array of constant). Then the lower bound and the upper bound
        define the length of this expansion.

        Examples
        ========

        >>> from sympy.tensor.array import ArrayComprehension
        >>> from sympy.abc import i, j, k
        >>> a = ArrayComprehension(10*i + j, (i, 1, 4), (j, 1, 3))
        >>> a.bounds
        ((i, 1, 4), (j, 1, 3))
        """
        return self._args[1:]

    @classmethod
    def _check_bounds_validity(cls, expr, bounds):
        bounds = sympify(bounds)
        for var, inf, sup in bounds:
            if any(not isinstance(i, Expr) for i in [inf, sup]):
                raise TypeError('Bounds should be an Expression(combination of Integer and Symbol)')
            if isinstance(inf, Integer) and isinstance(sup, Integer) and (inf > sup) == True:
                raise ValueError('Lower bound should be inferior to upper bound')
        return bounds

    def doit(self):
        expr = self._expr
        for index, value in enumerate(self.bounds):
            var, inf, sup = map(sympify, value)
            # Array will not ne expanded if there is a symbolic dimension
            if Basic(inf, sup).atoms(Symbol):
                return self

        arr = self._expand_array()
        return arr

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
