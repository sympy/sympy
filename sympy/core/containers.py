"""Module for SymPy containers

    (SymPy objects that store other SymPy objects)

    The containers implemented in this module are subclassed to Basic.
    They are supposed to work seamlessly within the SymPy framework.
"""

from sympy.core.basic import Basic
from sympy.core.sympify import sympify, converter
from sympy.utilities.iterables import iterable

class Tuple(Basic):
    """
    Wrapper around the builtin tuple object

    The Tuple is a subclass of Basic, so that it works well in the
    SymPy framework.  The wrapped tuple is available as self.args, but
    you can also access elements or slices with [:] syntax.

    >>> from sympy import symbols
    >>> from sympy.core.containers import Tuple
    >>> a, b, c, d = symbols('a b c d')
    >>> Tuple(a, b, c)[1:]
    (b, c)
    >>> Tuple(a, b, c).subs(a, d)
    (d, b, c)

    """

    def __new__(cls, *args, **assumptions):
        args = [ sympify(arg) for arg in args ]
        obj = Basic.__new__(cls, *args, **assumptions)
        return obj

    def __getitem__(self,i):
        if isinstance(i,slice):
            indices = i.indices(len(self))
            return Tuple(*[self.args[i] for i in range(*indices)])
        return self.args[i]

    def __len__(self):
        return len(self.args)

    def __contains__(self, item):
        return item in self.args

    def __iter__(self):
        return iter(self.args)

    def __add__(self, other):
        if isinstance(other, Tuple):
            return Tuple(*(self.args + other.args))
        elif isinstance(other, tuple):
            return Tuple(*(self.args + other))
        else:
            return NotImplemented

    def __radd__(self, other):
        if isinstance(other, Tuple):
            return Tuple(*(other.args + self.args))
        elif isinstance(other, tuple):
            return Tuple(*(other + self.args))
        else:
            return NotImplemented

    def __eq__(self, other):
        if isinstance(other, Basic):
            return super(Tuple, self).__eq__(other)
        return self.args == other

    def __ne__(self, other):
        if isinstance(other, Basic):
            return super(Tuple, self).__ne__(other)
        return self.args != other

    def __hash__(self):
        return hash(self.args)

    def _to_mpmath(self, prec):
        return tuple([a._to_mpmath(prec) for a in self.args])

    def __lt__(self, other):
        return self.args < other.args

converter[tuple] = lambda tup: Tuple(*tup)

def tuple_wrapper(method):
    """
    Decorator that converts any tuple in the function arguments into a Tuple.

    The motivation for this is to provide simple user interfaces.  The user can
    call a function with regular tuples in the argument, and the wrapper will
    convert them to Tuples before handing them to the function.

    >>> from sympy.core.containers import tuple_wrapper, Tuple
    >>> def f(*args):
    ...    return args
    >>> g = tuple_wrapper(f)

    The decorated function g sees only the Tuple argument:

    >>> g(0, (1, 2), 3)
    (0, (1, 2), 3)

    """
    def wrap_tuples(*args, **kw_args):
        newargs=[]
        for arg in args:
            if type(arg) is tuple:
                newargs.append(Tuple(*arg))
            else:
                newargs.append(arg)
        return method(*newargs, **kw_args)
    return wrap_tuples

class Dict(Basic):
    """
    Wrapper around the builtin dict object

    The Dict is a subclass of Basic, so that it works well in the
    SymPy framework.  Because it is immutable, it may be included
    in sets, but its values must all be given at instantiation and
    cannot be changed afterwards.  Otherwise it behaves identically
    to the Python dict.

    >>> from sympy import S
    >>> from sympy.core.containers import Dict

    >>> D = Dict({1:'one', 2:'two'})
    >>> for key in D: print key, D[key]
    1 one
    2 two

    The args are sympified so the 1 and 2 are Integers and the values
    are Symbols. Queries automatically sympify args so the following work:

    >>> 1 in D
    True
    >>> D.has('one') # searches keys and values
    True
    >>> 'one' in D # not in the keys
    False
    >>> D[1]
    one

    """

    def __new__(cls, *args):
        if len(args)==1 and args[0].__class__ is dict:
            items = [Tuple(k, v) for k, v in args[0].items()]
        elif iterable(args) and all(len(arg) == 2 for arg in args):
            items = [Tuple(k, v) for k, v in args]
        else:
            raise TypeError('Pass Dict args as Dict((k1, v1), ...) or Dict({k1: v1, ...})')
        elements = frozenset(items)
        obj = Basic.__new__(cls, elements)
        obj.elements = elements
        obj._dict = dict(items) # In case Tuple decides it wants to sympify
        return obj

    def __getitem__(self, key):
        """x.__getitem__(y) <==> x[y]"""
        return self._dict[sympify(key)]

    def __setitem__(self, key, value):
        raise NotImplementedError("SymPy Dicts are Immutable")

    @property
    def args(self):
        return tuple(self.elements)

    def items(self):
        '''D.items() -> list of D's (key, value) pairs, as 2-tuples'''
        return self._dict.items()

    def keys(self):
        '''D.keys() -> list of D's keys'''
        return self._dict.keys()

    def values(self):
        '''D.values() -> list of D's values'''
        return self._dict.values()

    def __iter__(self):
        '''x.__iter__() <==> iter(x)'''
        return iter(self._dict)

    def __len__(self):
        '''x.__len__() <==> len(x)'''
        return self._dict.__len__()

    def get(self, key, default=None):
        '''D.get(k[,d]) -> D[k] if k in D, else d.  d defaults to None.'''
        return self._dict.get(sympify(key), default)

    def __contains__(self, key):
        '''D.__contains__(k) -> True if D has a key k, else False'''
        return sympify(key) in self._dict

class TableForm(Basic):
    """
    Allows to create nice table representation of data.

    Example::

        >>> from sympy import TableForm
        >>> t = TableForm([[5, 7], [4, 2], [10, 3]])

    Then use "print t" to print the table. You can use the SymPy's printing
    system to produce tables in any format (ascii, latex, html, ...).

    """

    def __new__(cls, data, headings=None, alignment="left"):
        """
        Creates a TableForm.

        Parameters:

            data ... 2D data to be put into the table
            headings ... gives the labels for entries in each dimension:
                         None ... no labels in any dimension
                         "automatic" ... gives successive integer labels
                         [[l1, l2, ...], ...] gives labels for each entry in
                             each dimension (can be None for some dimension)
            alignment ... "left", "center", "right" (alignment of the columns)

        Example:

        >>> from sympy import TableForm
        >>> t = TableForm([[5, 7], [4, 2], [10, 3]])

        """
        # We only support 2D data. Check the consistency:
        from sympy.core.numbers import Integer
        _w = Integer(len(data[0]))
        _h = Integer(len(data))
        for line in data:
            assert len(line) == _w
        _lines = Tuple(*data)

        if headings is None:
            _headings = [None, None]
        elif headings == "automatic":
            _headings = [range(1, _h + 1), range(1, _w + 1)]
        else:
            h1, h2 = headings
            if h1 == "automatic":
                h1 = range(1, _h + 1)
            if h2 == "automatic":
                h2 = range(1, _w + 1)
            _headings = [h1, h2]

        _alignment = alignment

        obj = Basic.__new__(cls, _w, _h, _lines)
        obj._w = _w
        obj._h = _h
        obj._lines = _lines
        obj._headings = _headings
        obj._alignment = _alignment
        return obj

    def as_str(self):
        """
        Returns the string representation of 'self'.

        Example:

        >>> from sympy import TableForm
        >>> t = TableForm([[5, 7], [4, 2], [10, 3]])
        >>> s = t.as_str()

        """
        column_widths = [0] * self._w
        lines = []
        for line in self._lines:
            new_line = []
            for i in range(self._w):
                # Format the item somehow if needed:
                s = str(line[i])
                w = len(s)
                if w > column_widths[i]:
                    column_widths[i] = w
                new_line.append(s)
            lines.append(new_line)

        # Check heading:
        if self._headings[1]:
            new_line = []
            for i in range(self._w):
                # Format the item somehow if needed:
                s = str(self._headings[1][i])
                w = len(s)
                if w > column_widths[i]:
                    column_widths[i] = w
                new_line.append(s)
            self._headings[1] = new_line

        format_str = ""
        for w in column_widths:
            if self._alignment == "left":
                align = "-"
            elif self._alignment == "right":
                align = ""
            else:
                raise NotImplementedError()
            format_str += "%" + align + str(w) + "s "
        format_str += "\n"

        if self._headings[0]:
            self._headings[0] = [str(x) for x in self._headings[0]]
            heading_width = max([len(x) for x in self._headings[0]])
            format_str = "%" + str(heading_width) + "s | " + format_str

        s = ""
        if self._headings[1]:
            d = self._headings[1]
            if self._headings[0]:
                d = [""] + d
            first_line = format_str % tuple(d)
            s += first_line
            s += "-" * (len(first_line) - 2) + "\n"
        for i, line in enumerate(lines):
            d = line
            if self._headings[0]:
                d = [self._headings[0][i]] + d
            s += format_str % tuple(d)
        return s
