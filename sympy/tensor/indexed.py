"""Module that defines indexed objects

    The classes IndexedBase, Indexed and Idx would represent a matrix element
    M[i, j] as in the following graph::

       1) The Indexed class represents the entire indexed object.
                  |
               ___|___
              '       '
               M[i, j]
              /   \__\______
              |             |
              |             |
              |     2) The Idx class represent indices and each Idx can
              |        optionally contain information about its range.
              |
        3) IndexedBase represents the `stem' of an indexed object, here `M'.
           The stem used by itself is usually taken to represent the entire
           array.

    There can be any number of indices on an Indexed object.  No transformation
    properties are implemented in these Base objects, but implicit contraction
    of repeated indices is supported.

    Note that the support for complicated (i.e. non-atomic) integer expressions
    as indices is limited.  (This should be improved in future releases.)

    Examples
    ========

    To express the above matrix element example you would write:

    >>> from sympy.tensor import IndexedBase, Idx
    >>> from sympy import symbols
    >>> M = IndexedBase('M')
    >>> i, j = map(Idx, ['i', 'j'])
    >>> M[i, j]
    M[i, j]

    Repreated indices in a product implies a summation, so to express a
    matrix-vector product in terms of Indexed objects:

    >>> x = IndexedBase('x')
    >>> M[i, j]*x[j]
    M[i, j]*x[j]

    If the indexed objects will be converted to component based arrays, e.g.
    with the code printers or the autowrap framework, you also need to provide
    (symbolic or numerical) dimensions.  This can be done by passing an
    optional shape parameter to IndexedBase upon construction:

    >>> dim1, dim2 = symbols('dim1 dim2', integer=True)
    >>> A = IndexedBase('A', shape=(dim1, 2*dim1, dim2))
    >>> A.shape
    (dim1, 2*dim1, dim2)
    >>> A[i, j, 3].shape
    (dim1, 2*dim1, dim2)

    If an IndexedBase object has no shape information, it is assumed that the
    array is as large as the ranges of it's indices:

    >>> n, m = symbols('n m', integer=True)
    >>> i = Idx('i', m)
    >>> j = Idx('j', n)
    >>> M[i, j].shape
    (m, n)
    >>> M[i, j].ranges
    [(0, m - 1), (0, n - 1)]

    The above can be compared with the following:

    >>> A[i, 2, j].shape
    (dim1, 2*dim1, dim2)
    >>> A[i, 2, j].ranges
    [(0, m - 1), None, (0, n - 1)]

    To analyze the structure of indexed expressions, you can use the methods
    get_indices() and get_contraction_structure():

    >>> from sympy.tensor import get_indices, get_contraction_structure
    >>> get_indices(A[i, j, j])
    (set([i]), {})
    >>> get_contraction_structure(A[i, j, j])
    {(j,): set([A[i, j, j]])}

    See the appropriate docstrings for a detailed explanation of the output.

"""

#   TODO:  (some ideas for improvement)
#
#   o test and guarantee numpy compatibility
#      - implement full support for broadcasting
#      - strided arrays
#
#   o more functions to analyze indexed expressions
#      - identify standard constructs, e.g matrix-vector product in a subexpression
#
#   o functions to generate component based arrays (numpy and sympy.Matrix)
#      - generate a single array directly from Indexed
#      - convert simple sub-expressions
#
#   o sophisticated indexing (possibly in subclasses to preserve simplicity)
#      - Idx with range smaller than dimension of Indexed
#      - Idx with stepsize != 1
#      - Idx with step determined by function call

from sympy.core import Expr, Basic, Tuple, Symbol, Integer, sympify, S
from sympy.core.compatibility import is_sequence

class IndexException(Exception):
    pass

class IndexedBase(Expr):
    """Represent the base or stem of an indexed object

    The IndexedBase class represent an array that contains elements. The main purpose
    of this class is to allow the convenient creation of objects of the Indexed
    class.  The __getitem__ method of IndexedBase returns an instance of
    Indexed.  Alone, without indices, the IndexedBase class can be used as a
    notation for e.g. matrix equations, resembling what you could do with the
    Symbol class.  But, the IndexedBase class adds functionality that is not
    available for Symbol instances:

      -  An IndexedBase object can optionally store shape information.  This can
         be used in to check array conformance and conditions for numpy
         broadcasting.  (TODO)
      -  An IndexedBase object implements syntactic sugar that allows easy symbolic
         representation of array operations, using implicit summation of
         repeated indices.
      -  The IndexedBase object symbolizes a mathematical structure equivalent
         to arrays, and is recognized as such for code generation and automatic
         compilation and wrapping.

    >>> from sympy.tensor import IndexedBase, Idx
    >>> from sympy import symbols
    >>> A = IndexedBase('A'); A
    A
    >>> type(A)
    <class 'sympy.tensor.indexed.IndexedBase'>

    When an IndexedBase object recieves indices, it returns an array with named
    axes, represented by an Indexed object:

    >>> i, j = symbols('i j', integer=True)
    >>> A[i, j, 2]
    A[i, j, 2]
    >>> type(A[i, j, 2])
    <class 'sympy.tensor.indexed.Indexed'>

    The IndexedBase constructor takes an optional shape argument.  If given,
    it overrides any shape information in the indices. (But not the index
    ranges!)

    >>> m, n, o, p = symbols('m n o p', integer=True)
    >>> i = Idx('i', m)
    >>> j = Idx('j', n)
    >>> A[i, j].shape
    (m, n)
    >>> B = IndexedBase('B', shape=(o, p))
    >>> B[i, j].shape
    (o, p)

    """
    is_commutative = False

    def __new__(cls, label, shape=None, **kw_args):
        if isinstance(label, basestring):
            label = Symbol(label)

        obj = Expr.__new__(cls, label, **kw_args)
        if is_sequence(shape):
            obj._shape = Tuple(*shape)
        else:
            obj._shape = shape
        return obj

    @property
    def args(self):
        if self._shape:
            return self._args + (self._shape,)
        else:
            return self._args

    def _hashable_content(self):
        return Expr._hashable_content(self) + (self._shape,)

    def __getitem__(self, indices, **kw_args):
        if is_sequence(indices):
            # Special case needed because M[*my_tuple] is a syntax error.
            if self.shape and len(self.shape) != len(indices):
                raise IndexException("Rank mismatch")
            return Indexed(self, *indices, **kw_args)
        else:
            if self.shape and len(self.shape) != 1:
                raise IndexException("Rank mismatch")
            return Indexed(self, indices, **kw_args)

    @property
    def shape(self):
        return self._shape

    @property
    def label(self):
        return self.args[0]

    def _sympystr(self, p):
        return p.doprint(self.label)


class Indexed(Expr):
    """Represents a mathematical object with indices.

    >>> from sympy.tensor import Indexed, IndexedBase, Idx
    >>> from sympy import symbols
    >>> i, j = map(Idx, ['i', 'j'])
    >>> Indexed('A', i, j)
    A[i, j]

    It is recommended that Indexed objects are created via IndexedBase:

    >>> A = IndexedBase('A')
    >>> Indexed('A', i, j) == A[i, j]
    True

    """
    is_commutative = False

    def __new__(cls, base, *args, **kw_args):
        if not args: raise IndexException("Indexed needs at least one index")
        if isinstance(base, (basestring, Symbol)):
            base = IndexedBase(base)
        elif not isinstance(base, IndexedBase):
            raise TypeError("Indexed expects string, Symbol or IndexedBase as base")
        return Expr.__new__(cls, base, *args, **kw_args)

    @property
    def base(self):
        return self.args[0]

    @property
    def indices(self):
        return self.args[1:]

    @property
    def rank(self):
        """returns the number of indices"""
        return len(self.args)-1

    @property
    def shape(self):
        """returns a list with dimensions of each index.

        Dimensions is a property of the array, not of the indices.  Still, if
        the IndexedBase does not define a shape attribute, it is assumed that
        the ranges of the indices correspond to the shape of the array.

        >>> from sympy.tensor.indexed import IndexedBase, Idx
        >>> from sympy import symbols
        >>> n, m = symbols('n m', integer=True)
        >>> i = Idx('i', m)
        >>> j = Idx('j', m)
        >>> A = IndexedBase('A', shape=(n, n))
        >>> B = IndexedBase('B')
        >>> A[i, j].shape
        (n, n)
        >>> B[i, j].shape
        (m, m)
        """
        if self.base.shape:
            return self.base.shape
        try:
            return Tuple(*[i.upper - i.lower + 1 for i in self.indices])
        except AttributeError:
            raise IndexException("Range is not defined for all indices in: %s" % self)
        except TypeError:
            raise IndexException("Shape cannot be inferred from Idx with undefined range: %s"%self)

    @property
    def ranges(self):
        """returns a list of tuples with lower and upper range of each index

        If an index does not define the data members upper and lower, the
        corresponding slot in the list contains ``None'' instead of a tuple.
        """
        ranges = []
        for i in self.indices:
            try:
                ranges.append(Tuple(i.lower, i.upper))
            except AttributeError:
                ranges.append(None)
        return ranges

    def _sympystr(self, p):
        indices = map(p.doprint, self.indices)
        return "%s[%s]" % (p.doprint(self.base), ", ".join(indices))


class Idx(Expr):
    """Represents an index, either symbolic or integer.

    There are a number of ways to create an Idx object.  The constructor
    takes two arguments:

    ``label``
        An integer or a symbol that labels the index.
    ``range``
        Optionally you can specify a range as either

            - Symbol or integer: This is interpreted as dimension. lower and
              upper ranges are set to 0 and range-1
            - tuple: This is interpreted as the lower and upper bounds in the
              range.

    Note that the Idx constructor is rather pedantic, and will not accept
    non-integer symbols.  The only exception is that you can use oo and -oo to
    specify an unbounded range.  For all other cases, both label and bounds
    must be declared as integers, in the sense that for a index label n,
    n.is_integer must return True.

    For convenience, if the label is given as a string, it is automatically
    converted to an integer symbol.  (Note that this conversion is not done for
    range or dimension arguments.)

    :Examples:

    >>> from sympy.tensor import IndexedBase, Idx
    >>> from sympy import symbols, oo
    >>> n, i, L, U = symbols('n i L U', integer=True)

    0) Construction from a string.  An integer symbol is created from the
    string.

    >>> Idx('qwerty')
    qwerty

    1) Both upper and lower bound specified

    >>> idx = Idx(i, (L, U)); idx
    i
    >>> idx.lower, idx.upper
    (L, U)

    2) Only dimension specified, lower bound defaults to 0

    >>> idx = Idx(i, n); idx.lower, idx.upper
    (0, n - 1)
    >>> idx = Idx(i, 4); idx.lower, idx.upper
    (0, 3)
    >>> idx = Idx(i, oo); idx.lower, idx.upper
    (0, oo)

    3) No bounds given, interpretation of this depends on context.

    >>> idx = Idx(i); idx.lower, idx.upper
    (None, None)

    4) for a literal integer instead of a symbolic label the bounds are still
    there:

    >>> idx = Idx(2, n); idx.lower, idx.upper
    (0, n - 1)

    """

    is_integer = True

    def __new__(cls, label, range=None, **kw_args):

        if isinstance(label, basestring):
            label = Symbol(label, integer=True)
        label, range = map(sympify, (label, range))

        if not label.is_integer:
            raise TypeError("Idx object requires an integer label")

        elif is_sequence(range):
            assert len(range) == 2, "Idx got range tuple with wrong length"
            for bound in range:
                if not (bound.is_integer or abs(bound) is S.Infinity):
                    raise TypeError("Idx object requires integer bounds")
            args = label, Tuple(*range)
        elif isinstance(range, Expr):
            if not (range.is_integer or range is S.Infinity):
                raise TypeError("Idx object requires an integer dimension")
            args = label, Tuple(S.Zero, range-S.One)
        elif range:
            raise TypeError("range must be ordered iterable or integer sympy expression")
        else:
            args = label,

        obj = Expr.__new__(cls, *args, **kw_args)
        return obj

    @property
    def label(self):
        """Returns the name/label of the index, or it's integer value"""
        return self.args[0]

    @property
    def lower(self):
        """Returns the lower bound of the index"""
        try:
            return self.args[1][0]
        except IndexError:
            return

    @property
    def upper(self):
        """Returns the upper bound of the index"""
        try:
            return self.args[1][1]
        except IndexError:
            return

    def _sympystr(self, p):
        return p.doprint(self.label)
