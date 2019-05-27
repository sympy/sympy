from __future__ import print_function, division
import functools
from sympy.core.sympify import sympify
from sympy.core.expr import Expr
from sympy.core import Basic
from sympy.core.compatibility import Iterable
from sympy.tensor.array import MutableDenseNDimArray, ImmutableDenseNDimArray
from sympy import Symbol
from sympy.core.sympify import sympify
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
        arglist = [sympify(expr)]
        arglist.extend(cls._check_bounds_validity(expr, bounds))
        obj = Basic.__new__(cls, *arglist, **assumptions)
        obj._expr = obj._args[0]
        obj._bounds = obj._args[1:]
        obj._shape = cls._calculate_shape_from_bounds(obj._bounds)
        obj._rank = len(obj._shape)
        obj._loop_size = cls._calculate_loop_size(obj._shape)
        return obj

    @property
    def expr(self):
        """
        Returns the expression that will be expanded to an array

        Examples
        ========

        >>> from sympy.tensor.array import ArrayComprehension
        >>> from sympy.abc import i, j, k
        >>> a = ArrayComprehension(10*i + j, (i, 1, 4), (j, 1, 3))
        >>> a.expr
        10*i + j
        """
        return self._expr

    @property
    def bounds(self):
        """
        Returns a list of the bounds that will be applied while expanding the array. Each
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

    @property
    def free_symbols(self):
        """
        Returns a set of the free_symbols in the array. Variables appeared in the bounds
        are supposed to be excluded from the free symbol set.

        Examples
        ========

        >>> from sympy.tensor.array import ArrayComprehension
        >>> from sympy.abc import i, j, k
        >>> a = ArrayComprehension(10*i + j, (i, 1, 4), (j, 1, 3))
        >>> a.free_symbols
        set()
        >>> b = ArrayComprehension(10*i + j, (i, 1, 4), (j, 1, k+3))
        >>> b.free_symbols
        {k}
        """
        expr_free_sym = self.expr.free_symbols
        for var, inf, sup in self.bounds:
            expr_free_sym.discard(var)
            if len(inf.free_symbols) > 0:
                expr_free_sym = expr_free_sym.union(inf.free_symbols)
            if len(sup.free_symbols) > 0:
                expr_free_sym = expr_free_sym.union(sup.free_symbols)
        return expr_free_sym

    @property
    def shape(self):
        """
        Returns the shape of the expanded array, which can have symbols. Note that both
        the lower and the upper bounds are included while calculating the shape.

        Examples
        ========

        >>> from sympy.tensor.array import ArrayComprehension
        >>> from sympy.abc import i, j, k
        >>> a = ArrayComprehension(10*i + j, (i, 1, 4), (j, 1, 3))
        >>> a.shape
        (4, 3)
        >>> b = ArrayComprehension(10*i + j, (i, 1, 4), (j, 1, k+3))
        >>> b.shape
        (4, k + 3)
        """
        return self._shape

    def rank(self):
        """
        Returns the rank of the expanded array.

        Examples
        ========

        >>> from sympy.tensor.array import ArrayComprehension
        >>> from sympy.abc import i, j, k
        >>> a = ArrayComprehension(10*i + j, (i, 1, 4), (j, 1, 3))
        >>> a.rank()
        2
        """
        return self._rank

    def __len__(self):
        """
        Overload common function len().Returns the number of element in the expanded
        array. Note that symbolic length is not supported and will raise an error.

        Examples
        ========

        >>> from sympy.tensor.array import ArrayComprehension
        >>> from sympy.abc import i, j, k
        >>> a = ArrayComprehension(10*i + j, (i, 1, 4), (j, 1, 3))
        >>> len(a)
        12
        """
        if len(self._loop_size.free_symbols) != 0:
            raise ValueError('Symbolic length is not supported')
        return self._loop_size

    @classmethod
    def _check_bounds_validity(cls, expr, bounds):
        bounds = sympify(bounds)
        for var, inf, sup in bounds:
            if any(not isinstance(i, Expr) for i in [inf, sup]):
                raise TypeError('Bounds should be an Expression(combination of Integer and Symbol)')
            if (inf > sup) == True:
                raise ValueError('Lower bound should be inferior to upper bound')
            if var in inf.free_symbols or var in sup.free_symbols:
                raise ValueError('Variable should not be part of its bounds')
        return bounds

    @classmethod
    def _calculate_shape_from_bounds(cls, bounds):
        shape = []
        for var, inf, sup in bounds:
            shape.append(sup - inf + 1)
        return tuple(shape)

    @classmethod
    def _calculate_loop_size(cls, shape):
        if len(shape) == 0:
            return 0
        loop_size = 1
        for l in shape:
            loop_size = loop_size * l

        return loop_size

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
