from basic import Basic, Atom
from core import C
from sympify import sympify, _sympify, SympifyError
from singleton import S, Singleton
from expr import Expr, AtomicExpr
from cache import cacheit
from compatibility import any, all
from function import FunctionClass
from sympy.logic.boolalg import Boolean

import re

class Symbol(AtomicExpr, Boolean):
    """
    Assumptions::
       commutative = True

    You can override the default assumptions in the constructor::
       >>> from sympy import symbols
       >>> A,B = symbols('A,B', commutative = False)
       >>> bool(A*B != B*A)
       True
       >>> bool(A*B*2 == 2*A*B) == True # multiplication by scalars is commutative
       True

    """

    is_comparable = False

    __slots__ = ['is_commutative', 'name']

    is_Symbol = True

    def __new__(cls, name, commutative=True, **assumptions):
        """Symbols are identified by name and commutativity::

        >>> from sympy import Symbol
        >>> bool(Symbol("x") == Symbol("x")) == True
        True
        >>> bool(Symbol("x", real=True) == Symbol("x", real=False)) == True
        True
        >>> bool(Symbol("x", commutative=True) ==
        ...      Symbol("x", commutative=False)) == True
        False

        """

        if 'dummy' in assumptions:
            import warnings
            warnings.warn(
                    "\nThe syntax Symbol('x', dummy=True) is deprecated and will"
                    "\nbe dropped in a future version of Sympy. Please use Dummy()"
                    "\nor symbols(..., cls=Dummy) to create dummy symbols.",
                    DeprecationWarning)
            if assumptions.pop('dummy'):
                return Dummy(name, commutative, **assumptions)
        return Symbol.__xnew_cached_(cls, name, commutative, **assumptions)

    def __new_stage2__(cls, name, commutative=True, **assumptions):
        assert isinstance(name, str),`type(name)`
        obj = Expr.__new__(cls, **assumptions)
        obj.is_commutative = commutative
        obj.name = name
        return obj

    __xnew__       = staticmethod(__new_stage2__)            # never cached (e.g. dummy)
    __xnew_cached_ = staticmethod(cacheit(__new_stage2__))   # symbols are always cached

    def __getnewargs__(self):
        return (self.name, self.is_commutative)

    def _hashable_content(self):
        return (self.is_commutative, self.name)

    def as_dummy(self):
        return Dummy(self.name, self.is_commutative, **self.assumptions0)

    def __call__(self, *args):
        from function import Function
        return Function(self.name, nargs=len(args))(*args, **self.assumptions0)

    def as_real_imag(self, deep=True):
        return (C.re(self), C.im(self))

    def _eval_expand_complex(self, deep=True, **hints):
        re, im = self.as_real_imag()
        return re + im*S.ImaginaryUnit

    def _sage_(self):
        import sage.all as sage
        return sage.var(self.name)

    @property
    def is_number(self):
        return False

    @property
    def free_symbols(self):
        return set([self])

class Dummy(Symbol):
    """Dummy symbols are each unique, identified by an internal count index ::

    >>> from sympy import Dummy
    >>> bool(Dummy("x") == Dummy("x")) == True
    False

    If a name is not supplied then a string value of the count index will be
    used. This is useful when a temporary variable is needed and the name
    of the variable used in the expression is not important. ::
    >>> Dummy._count = 0 # /!\ this should generally not be changed; it is being
    >>> Dummy()          # used here to make sure that the doctest passes.
    _0

    """

    _count = 0

    __slots__ = ['dummy_index']

    is_Dummy = True

    def __new__(cls, name=None, commutative=True, **assumptions):
        if name is None:
            name = str(Dummy._count)

        obj = Symbol.__xnew__(cls, name, commutative=commutative, **assumptions)

        Dummy._count += 1
        obj.dummy_index = Dummy._count
        return obj

    def _hashable_content(self):
        return Symbol._hashable_content(self) + (self.dummy_index,)

class Wild(Symbol):
    """
    Wild() matches any expression but another Wild().
    """

    __slots__ = ['exclude', 'properties']

    is_Wild = True

    def __new__(cls, name, exclude=None, properties=None, **assumptions):
        if type(exclude) is list:
            exclude = tuple(exclude)
        if type(properties) is list:
            properties = tuple(properties)

        return Wild.__xnew__(cls, name, exclude, properties, **assumptions)

    def __getnewargs__(self):
        return (self.name, self.exclude, self.properties)

    @staticmethod
    @cacheit
    def __xnew__(cls, name, exclude, properties, **assumptions):
        obj = Symbol.__xnew__(cls, name, **assumptions)

        if exclude is None:
            obj.exclude = None
        else:
            obj.exclude = tuple([sympify(x) for x in exclude])
        if properties is None:
            obj.properties = None
        else:
            obj.properties = tuple(properties)
        return obj

    def _hashable_content(self):
        return (self.name, self.exclude, self.properties )

    # TODO add check against another Wild
    def matches(self, expr, repl_dict={}, evaluate=False):
        if self in repl_dict:
            if repl_dict[self] == expr:
                return repl_dict
            else:
                return None
        if self.exclude:
            for x in self.exclude:
                if x in expr:
                    return None
        if self.properties:
            for f in self.properties:
                if not f(expr):
                    return None
        repl_dict = repl_dict.copy()
        repl_dict[self] = expr
        return repl_dict

    def __call__(self, *args, **assumptions):
        from sympy.core.function import WildFunction
        return WildFunction(self.name, nargs=len(args))(*args, **assumptions)

class Pure(Expr):
    """A commutative singleton symbol different from all other symbols. """
    __metaclass__ = Singleton

    __slots__ = ['is_commutative', 'name']
    is_commutative, name = True, 'pure'

    is_Pure   = True

_re_var_range = re.compile(r"^(.*?)(\d*):(\d+)$")
_re_var_split = re.compile(r"\s|,")

def symbols(names, **args):
    """
    Transform strings into instances of :class:`Symbol` class.

    :func:`symbols` function returns a sequence of symbols with names taken
    from ``names`` argument, which can be a comma or whitespace delimited
    string, or a sequence of strings::

        >>> from sympy import symbols, Function

        >>> x, y, z = symbols('x,y,z')
        >>> a, b, c = symbols('a b c')

    The type of output is dependent on the properties of input arguments::

        >>> x = symbols('x')
        >>> (x,) = symbols('x,')

        >>> symbols(('a', 'b', 'c'))
        (a, b, c)
        >>> symbols(['a', 'b', 'c'])
        [a, b, c]
        >>> symbols(set(['a', 'b', 'c']))
        set([a, b, c])

    If an iterable container is needed set ``seq`` argument to ``True``::

        >>> symbols('x', seq=True)
        (x,)

    To cut on typing, range syntax is supported co create indexed symbols::

        >>> symbols('x:10')
        (x0, x1, x2, x3, x4, x5, x6, x7, x8, x9)

        >>> symbols('x5:10')
        (x5, x6, x7, x8, x9)

        >>> symbols('x5:10,y:5')
        (x5, x6, x7, x8, x9, y0, y1, y2, y3, y4)

        >>> symbols(('x5:10', 'y:5'))
        ((x5, x6, x7, x8, x9), (y0, y1, y2, y3, y4))

    All newly created symbols have assumptions set accordingly to ``args``::

        >>> a = symbols('a', integer=True)
        >>> a.is_integer
        True

        >>> x, y, z = symbols('x,y,z', real=True)
        >>> x.is_real and y.is_real and z.is_real
        True

    Despite its name, :func:`symbols` can create symbol--like objects of
    other type, for example instances of Function or Wild classes. To
    achieve this, set ``cls`` keyword argument to the desired type::

        >>> symbols('f,g,h', cls=Function)
        (f, g, h)

        >>> type(_[0])
        <class 'sympy.core.function.UndefinedFunction'>

    """
    result = []

    if isinstance(names, basestring):
        names = _re_var_split.split(names)

        cls = args.pop('cls', Symbol)
        seq = args.pop('seq', False)

        for name in names:
            if not name:
                continue

            match = _re_var_range.match(name)

            if match is not None:
                name, start, end = match.groups()

                if not start:
                    start = 0
                else:
                    start = int(start)

                for i in xrange(start, int(end)):
                    symbol = cls("%s%i" % (name, i), **args)
                    result.append(symbol)

                seq = True
            else:
                symbol = cls(name, **args)
                result.append(symbol)

        if not seq and len(result) <= 1:
            if not result:
                return None
            elif names[-1]:
                return result[0]

        return tuple(result)
    else:
        for name in names:
            syms = symbols(name, **args)

            if syms is not None:
                result.append(syms)

        return type(names)(result)

def var(names, **args):
    """
    Create symbols and inject them into the global namespace.

    This calls :func:`symbols` with the same arguments and puts the results
    into the *global* namespace. It's recommended not to use :func:`var` in
    library code, where :func:`symbols` has to be used::

        >>> from sympy import var

        >>> var('x')
        x
        >>> x
        x

        >>> var('a,ab,abc')
        (a, ab, abc)
        >>> abc
        abc

        >>> var('x,y', real=True)
        (x, y)
        >>> x.is_real and y.is_real
        True

    See :func:`symbol` documentation for more details on what kinds of
    arguments can be passed to :func:`var`.

    """
    def traverse(symbols, frame):
        """Recursively inject symbols to the global namespace. """
        for symbol in symbols:
            if isinstance(symbol, Basic):
                frame.f_globals[symbol.name] = symbol
            elif isinstance(symbol, FunctionClass):
                frame.f_globals[symbol.__name__] = symbol
            else:
                traverse(symbol, frame)

    from inspect import currentframe
    frame = currentframe().f_back

    try:
        syms = symbols(names, **args)

        if syms is not None:
            if isinstance(syms, Basic):
                frame.f_globals[syms.name] = syms
            elif isinstance(syms, FunctionClass):
                frame.f_globals[syms.__name__] = syms
            else:
                traverse(syms, frame)
    finally:
        del frame # break cyclic dependencies as stated in inspect docs

    return syms
