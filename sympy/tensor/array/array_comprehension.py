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
    def __new__(cls, function, *symbols, **assumptions):
        if any(len(l) != 3 or None for l in symbols):
            raise ValueError('ArrayComprehension requires values lower and upper bound'
                              ' for the expression')
        arglist = [sympify(function)]
        arglist.extend(cls._check_limits_validity(function, symbols))
        obj = Basic.__new__(cls, *arglist, **assumptions)
        obj._function = obj._args[0]
        obj._limits = obj._args[1:]
        obj._shape = cls._calculate_shape_from_limits(obj._limits)
        obj._rank = len(obj._shape)
        obj._loop_size = cls._calculate_loop_size(obj._shape)
        return obj

    @property
    def function(self):
        """
        Return the function applied across limits

        Examples
        ========

        >>> from sympy.tensor.array import ArrayComprehension
        >>> from sympy.abc import i, j, k
        >>> a = ArrayComprehension(10*i + j, (i, 1, 4), (j, 1, 3))
        >>> a.function
        10*i + j
        """
        return self._function

    @property
    def limits(self):
        """
        Return a list of the limits that will be applied while expanding the array. Each
        bound contrains firstly the an variable (not necessarily to be component of the
        expression, e.g. an array of constant). Then the lower bound and the upper bound
        define the length of this expansion.

        Examples
        ========

        >>> from sympy.tensor.array import ArrayComprehension
        >>> from sympy.abc import i, j, k
        >>> a = ArrayComprehension(10*i + j, (i, 1, 4), (j, 1, 3))
        >>> a.limits
        ((i, 1, 4), (j, 1, 3))
        """
        return self._limits

    @property
    def free_symbols(self):
        """
        Return a set of the free_symbols in the array. Variables appeared in the bounds
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
        expr_free_sym = self._function.free_symbols
        for var, inf, sup in self._limits:
            expr_free_sym.discard(var)
            if len(inf.free_symbols) > 0:
                expr_free_sym = expr_free_sym.union(inf.free_symbols)
            if len(sup.free_symbols) > 0:
                expr_free_sym = expr_free_sym.union(sup.free_symbols)
        return expr_free_sym

    @property
    def variables(self):
        """
        Return a list of the variables in the limits

        Examples
        ========

        >>> from sympy.tensor.array import ArrayComprehension
        >>> from sympy.abc import i, j, k
        >>> a = ArrayComprehension(10*i + j, (i, 1, 4), (j, 1, 3))
        >>> a.variables
        [i, j]
        """
        return [l[0] for l in self._limits]

    @property
    def bound_symbols(self):
        """
        Return only variables that are dummy variables. Note that all variables are
        dummy variables since a limit without lower bound or upper bound is not accpted.
        """
        return [l[0] for l in self._limits if len(l) != 1]

    @property
    def shape(self):
        """
        Return the shape of the expanded array, which can have symbols. Note that both
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

    @property
    def is_numeric(self):
        """
        Return True if the expanded array is numeric, which means that there is not
        symbolic dimension.

        Examples
        ========

        >>> from sympy.tensor.array import ArrayComprehension
        >>> from sympy.abc import i, j, k
        >>> a = ArrayComprehension(10*i + j, (i, 1, 4), (j, 1, 3))
        >>> a.is_numeric
        True
        >>> b = ArrayComprehension(10*i + j, (i, 1, 4), (j, 1, k+3))
        >>> b.is_numeric
        False

        """
        for var, inf, sup in self._limits:
            if Basic(inf, sup).atoms(Symbol):
                return False
        return True

    def rank(self):
        """
        Return the rank of the expanded array.

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
    def _check_limits_validity(cls, function, limits):
        limits = sympify(limits)
        for var, inf, sup in limits:
            if any(not isinstance(i, Expr) for i in [inf, sup]):
                raise TypeError('Bounds should be an Expression(combination of Integer and Symbol)')
            if (inf > sup) == True:
                raise ValueError('Lower bound should be inferior to upper bound')
            if var in inf.free_symbols or var in sup.free_symbols:
                raise ValueError('Variable should not be part of its bounds')
        return limits

    @classmethod
    def _calculate_shape_from_limits(cls, limits):
        shape = []
        for var, inf, sup in limits:
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
        if not self.is_numeric:
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

        list_gen = self._function
        for var, inf, sup in reversed(self._limits):
            list_expr = list_gen
            list_gen = []
            for val in range(inf, sup+1):
                if not isinstance(list_expr, Iterable):
                    list_gen.append(list_expr.subs(var, val))
                else:
                    list_gen.append(_array_subs(list_expr, var, val))
        return ImmutableDenseNDimArray(list_gen)
