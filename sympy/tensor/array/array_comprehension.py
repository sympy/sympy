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
    Generates a list comprehension
    If there is a symbolic dimension, for example, say [i for i in range(1, N)] where
    N is a Symbol, then the expression will not be expanded to an array. Otherwise,
    calling the doit() function will launch the expansion.

    Examples
    ========

    >>> from sympy.tensor.array import ArrayComprehension
    >>> from sympy import symbols
    >>> i, j, k = symbols('i j k')
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
        return Basic.__new__(cls, *arglist, **assumptions)

    @property
    def function(self):
        """Use for obtaining the function applied across limits

        Returns
        =======

        Expr
            The applied function

        Examples
        ========

        >>> from sympy.tensor.array import ArrayComprehension
        >>> from sympy import symbols
        >>> i, j = symbols('i j')
        >>> a = ArrayComprehension(10*i + j, (i, 1, 4), (j, 1, 3))
        >>> a.function
        10*i + j
        """
        return self.args[0]

    @property
    def limits(self):
        """
        Used for obtaining the list of limits that will be applied
        while expanding the array

        Returns
        =======

        list
            The required list of the limits

        Examples
        ========

        >>> from sympy.tensor.array import ArrayComprehension
        >>> from sympy import symbols
        >>> i, j = symbols('i j')
        >>> a = ArrayComprehension(10*i + j, (i, 1, 4), (j, 1, 3))
        >>> a.limits
        ((i, 1, 4), (j, 1, 3))
        """
        return self.args[1:]

    @property
    def free_symbols(self):
        """
        Used for obtaining the set of the free_symbols in the array
        Variables appeared in the bounds are supposed to be excluded
        from the free symbol set.

        Returns
        =======

        set
            The set of the free_symbols

        Examples
        ========

        >>> from sympy.tensor.array import ArrayComprehension
        >>> from sympy import symbols
        >>> i, j, k = symbols('i j k')
        >>> a = ArrayComprehension(10*i + j, (i, 1, 4), (j, 1, 3))
        >>> a.free_symbols
        set()
        >>> b = ArrayComprehension(10*i + j, (i, 1, 4), (j, 1, k+3))
        >>> b.free_symbols
        {k}
        """
        expr_free_sym = self.function.free_symbols
        for var, inf, sup in self.limits:
            expr_free_sym.discard(var)
            if inf.free_symbols:
                expr_free_sym = expr_free_sym.union(inf.free_symbols)
            if sup.free_symbols:
                expr_free_sym = expr_free_sym.union(sup.free_symbols)
        return expr_free_sym

    @property
    def variables(self):
        """Used for obtaining the tuples of the variables in the limits

        Returns
        =======

        tuple
            The variables in the limits

        Examples
        ========

        >>> from sympy.tensor.array import ArrayComprehension
        >>> from sympy import symbols
        >>> i, j, k = symbols('i j k')
        >>> a = ArrayComprehension(10*i + j, (i, 1, 4), (j, 1, 3))
        >>> a.variables
        [i, j]
        """
        return [l[0] for l in self.limits]

    @property
    def bound_symbols(self):
        """Used to obtain the dummay variables

        Returns
        =======

        list
            The list of bound symbols

        Note
        ====

        Note that all variables are dummy variables since a limit without
        lower bound or upper bound is not accepted.
        """
        return [l[0] for l in self.limits if len(l) != 1]

    @property
    def shape(self):
        """

        Used for obtaining the shape of the expanded array, which may have symbols

        Returns
        =======

        tuple
            The shape of the expanded array.

        Note
        ====

        Both the lower and the upper bounds are included while
        calculating the shape.

        Examples
        ========

        >>> from sympy.tensor.array import ArrayComprehension
        >>> from sympy import symbols
        >>> i, j, k = symbols('i j k')
        >>> a = ArrayComprehension(10*i + j, (i, 1, 4), (j, 1, 3))
        >>> a.shape
        (4, 3)
        >>> b = ArrayComprehension(10*i + j, (i, 1, 4), (j, 1, k+3))
        >>> b.shape
        (4, k + 3)
        """
        return self._calculate_shape_from_limits(self.limits)

    @property
    def is_numeric(self):
        """
        Used to test if the array is numeric which means there is no symbolic
        dimension

        Returns
        =======

        Boolean
            True if the expanded array is numeric, False otherwise

        Examples
        ========

        >>> from sympy.tensor.array import ArrayComprehension
        >>> from sympy import symbols
        >>> i, j, k = symbols('i j k')
        >>> a = ArrayComprehension(10*i + j, (i, 1, 4), (j, 1, 3))
        >>> a.is_numeric
        True
        >>> b = ArrayComprehension(10*i + j, (i, 1, 4), (j, 1, k+3))
        >>> b.is_numeric
        False
        """
        for _, inf, sup in self.limits:
            if Basic(inf, sup).atoms(Symbol):
                return False
        return True

    def rank(self):
        """Used to obtain the rank of the expanded array

        Returns
        =======

        int
            The rank

        Examples
        ========

        >>> from sympy.tensor.array import ArrayComprehension
        >>> from sympy import symbols
        >>> i, j, k = symbols('i j k')
        >>> a = ArrayComprehension(10*i + j, (i, 1, 4), (j, 1, 3))
        >>> a.rank()
        2
        """
        return len(self.shape)

    def __len__(self):
        """
        Used to obtain the length of the expanded array which means the number
        of elements in the array.

        Returns
        =======

        int
            The required length

        Raises
        ======

        ValueError : When the length of the array is symbolic

        Examples
        ========

        >>> from sympy.tensor.array import ArrayComprehension
        >>> from sympy import symbols
        >>> i, j = symbols('i j')
        >>> a = ArrayComprehension(10*i + j, (i, 1, 4), (j, 1, 3))
        >>> len(a)
        12
        """
        loop_size = self._calculate_loop_size(self.shape)
        if loop_size.free_symbols:
            raise ValueError('Symbolic length is not supported')
        return loop_size

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
        for _, inf, sup in limits:
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

        list_gen = self.function
        for var, inf, sup in reversed(self.limits):
            list_expr = list_gen
            list_gen = []
            for val in range(inf, sup+1):
                if not isinstance(list_expr, Iterable):
                    list_gen.append(list_expr.subs(var, val))
                else:
                    list_gen.append(_array_subs(list_expr, var, val))
        return ImmutableDenseNDimArray(list_gen)
