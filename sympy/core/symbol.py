from basic import Basic, Atom
from core import C
from sympify import sympify, _sympify, SympifyError
from singleton import S, Singleton
from expr import Expr, AtomicExpr
from cache import cacheit
from compatibility import any, all
from sympy.logic.boolalg import Boolean

import re

class Symbol(AtomicExpr, Boolean):
    """
    Assumptions::
       commutative = True

    You can override the default assumptions in the constructor::
       >>> from sympy import symbols
       >>> A,B = symbols('AB', commutative = False)
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

def symbols(*names, **kwargs):
    """
    Return a list of symbols with names taken from 'names'
    argument, which can be a string, then each character
    forms a separate symbol, or a sequence of strings.

    >>> from sympy import symbols
    >>> x, y, z = symbols('xyz')

    Please note that this syntax is deprecated and will be dropped in a
    future version of sympy. Use comma or whitespace separated characters
    instead. Currently the old behavior is standard, this can be changed
    using the 'each_char' keyword:

    >>> symbols('xyz', each_char=False)
    xyz

    All newly created symbols have assumptions set accordingly
    to 'kwargs'. Main intention behind this function is to
    simplify and shorten examples code in doc-strings.

    >>> a = symbols('a', integer=True)
    >>> a.is_integer
    True
    >>> xx, yy, zz = symbols('xx', 'yy', 'zz', real=True)
    >>> xx.is_real and yy.is_real and zz.is_real
    True

    """
    func = Symbol
    # when polys12 is in place remove from here...
    if 'cls' in kwargs and kwargs.pop('cls') is Dummy:
        func = Dummy
    # ... to here when polys12 is in place
    # use new behavior if space or comma in string
    if not 'each_char' in kwargs and len(names) == 1 and \
    isinstance(names[0], str) and (' ' in names[0] or ',' in names[0]):
        kwargs['each_char'] = False
    if not kwargs.pop("each_char", True):
        # the new way:
        s = names[0]
        if not isinstance(s, list):
            s = re.split('\s|,', s)
        res = []
        for t in s:
            # skip empty strings
            if not t:
                continue
            sym = func(t, **kwargs)
            res.append(sym)
        res = tuple(res)
        if len(res) == 0:   # var('')
            res = None
        elif len(res) == 1: # var('x')
            res = res[0]
                            # otherwise var('a b ...')
        return res
    else:
        # this is the old, deprecated behavior:
        if len(names) == 1:
            result = [ func(name, **kwargs) for name in names[0] ]
        else:
            result = [ func(name, **kwargs) for name in names ]
        if len(result) == 1:
            return result[0]
        else:
            return result

def var(*names, **kwargs):
    """
    Create symbols and inject them into global namespace.

    This calls symbols() with the same arguments and puts the results into
    global namespace. Unlike symbols(), it uses each_char=False by default
    for compatibility reasons.

    NOTE: The new variable is both returned and automatically injected into
    the parent's *global* namespace.  It's recommended not to use "var" in
    library code, it is better to use symbols() instead.

    >>> from sympy import var
    >>> var('m')
    m
    >>> var('n xx yy zz')
    (n, xx, yy, zz)
    >>> n
    n
    >>> var('x y', real=True)
    (x, y)
    >>> x.is_real and y.is_real
    True

    """
    import inspect
    frame = inspect.currentframe().f_back
    try:
        kwargs['each_char'] = False
        s = symbols(*names, **kwargs)
        if s is None:
            return s
        if isinstance(s, Symbol):
            s_list = [s]
        else:
            s_list = s
        for t in s_list:
            frame.f_globals[t.name] = t
        return s
    finally:
        # we should explicitly break cyclic dependencies as stated in inspect
        # doc
        del frame
