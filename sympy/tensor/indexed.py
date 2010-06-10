"""Module that defines indexed objects with arbitrary transformation properties.

    The Idx class represent indices and each Idx can optionally contain
    information about its range.  Transformation properties can be defined in
    subclasses.

    Examples
    ========

    To express a matrix-vector product in terms of Indexed objects:


    >>> from sympy.tensor import Indexed, Idx
    >>> from sympy import symbols, Eq
    >>> i,j,n,m,M,x,y = symbols('i j n m M x y')
    >>> M = Indexed(M, Idx(i, m), Idx(j, n))
    >>> x = Indexed(x, Idx(j, n))
    >>> y = Indexed(y, Idx(i, m))
    >>> Eq(y, M*x)
    y(i) == M(i, j)*x(j)

"""

from sympy import Expr, Basic, SymTuple, Symbol, Integer, sympify, S


class Indexed(Expr):
    """Represent an arbitrary object with indices, e.g. a symbolic array element.

    >>> from sympy.tensor import Indexed
    >>> from sympy.abc import a, x, y
    >>> Indexed(a, x, y)
    a(x, y)

    """

    @property
    def label(self):
        return self.args[0]

    @property
    def indices(self):
        return self.args[1:]

    @property
    def rank(self):
        """returns the number of indices"""
        return len(args)-1

    @property
    def dimensions(self):
        return [(i.lower, i.upper) for i in self.indices]

    def _sympystr(self, p):
        indices = map(p.doprint, self.indices)
        return "%s(%s)" % (p.doprint(self.label), ", ".join(indices))


class Idx(Basic):
    """Represents an index, either symbolic or integer.

    Optionally you can specify a range [default=0]
    .. Symbol, integer  --  interpreted as dimension, lower and upper ranges are
                            set to 0 and range-1
    .. tuple  --  interpreted as lower, upper elements in range.

    You can use oo or -oo to specify an unbounded range, and the default range
    (0,0) means that the index is fixed.

    >>> from sympy.tensor import Indexed, Idx
    >>> from sympy import symbols, oo
    >>> n, x, L, U = symbols('n x L U')

    1) Both upper and lower bound specified

    >>> idx = Idx(x, (L, U)); idx
    x
    >>> idx.lower, idx.upper
    (L, U)

    2) Only dimension specified, lower bound defaults to 0

    >>> idx = Idx(x, n); idx.lower, idx.upper
    (0, -1 + n)
    >>> idx = Idx(x, 4); idx.lower, idx.upper
    (0, 3)
    >>> idx = Idx(x, oo); idx.lower, idx.upper
    (0, oo)

    3) No bounds given, interpretation of this depends on context.

    >>> idx = Idx(x); idx.lower, idx.upper
    (None, None)


    """

    def __new__(cls, label, range=None, **kw_args):

        label, range = map(sympify, (label, range))

        if label.is_Number:
            obj = Basic.__new__(cls, label, SymTuple(S.Zero, S.Zero), **kw_args)
            return obj

        if isinstance(range, (tuple, list, SymTuple)):
            assert len(range) == 2, "Idx got range tuple with wrong length"
            args = label, SymTuple(*range)
        elif isinstance(range, (Symbol, Integer)) or range is S.Infinity:
            args = label, SymTuple(S.Zero, range-S.One)
        elif range:
            raise TypeError("range must be tuple, symbol or integer")
        else:
            args = label, SymTuple(None, None)

        obj = Basic.__new__(cls, *args, **kw_args)
        return obj

    @property
    def label(self):
        return self.args[0]

    @property
    def lower(self):
        return self.args[1][0]

    @property
    def upper(self):
        return self.args[1][1]

    def _sympystr(self, p):
        return p.doprint(self.label)

