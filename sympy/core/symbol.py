
from basic import Basic, Atom, S, C, sympify
from methods import RelMeths, ArithMeths
from cache import cacheit

# from function import Function, WildFunction   /cyclic/

class Symbol(Atom, RelMeths, ArithMeths):
    """
    Assumptions::
       real = True
       commutative = True

    You can override the default assumptions in the constructor::
       >>> A,B = symbols('AB', commutative = False)
       >>> bool(A*B != B*A)
       True
       >>> bool(A*B*2 == 2*A*B) == True # multiplication by scalars is commutative
       True
    """

    is_comparable = False

    __slots__ = ['name']

    is_Symbol = True

    def __new__(cls, name, commutative=True, dummy=False,
                **assumptions):
        """if dummy == True, then this Symbol is totally unique, i.e.::

        >>> bool(Symbol("x") == Symbol("x")) == True
        True

        but with the dummy variable ::

        >>> bool(Symbol("x", dummy = True) == Symbol("x", dummy = True)) == True
        False

        """

        # XXX compatibility stuff
        if dummy==True:
            return Dummy(name, commutative=commutative, **assumptions)
        else:
            return Symbol.__xnew_cached_(cls, name, commutative, **assumptions)

    def __new_stage2__(cls, name, commutative=True, **assumptions):
        obj = Basic.__new__(cls,
                            commutative=commutative,
                            **assumptions)
        assert isinstance(name, str),`type(name)`
        obj.name = name
        return obj

    __xnew__       = staticmethod(__new_stage2__)            # never cached (e.g. dummy)
    __xnew_cached_ = staticmethod(cacheit(__new_stage2__))   # symbols are always cached

    def _hashable_content(self):
        return (self.name,)

    def tostr(self, level=0):
        return self.name

    def torepr(self):
        return '%s(%r)' % (self.__class__.__name__, self.name)

    def as_dummy(self):
        assumptions = self._assumptions.copy()
        return Dummy(self.name, **assumptions)

    def __call__(self, *args):
        assumptions = self._assumptions
        return Function(self.name, nargs=len(args))(*args, **assumptions)

    def _eval_integral(self, s):
        if self==s:
            return self**2/2
        return self*s

    def _eval_defined_integral(self, s, a, b):
        if self==s:
            return (b**2-a**2)/2
        return self*(b-a)

    def _eval_expand_complex(self, *args):
        return C.re(self) + C.im(self)*S.ImaginaryUnit

    def _sage_(self):
        import sage.all as sage
        return sage.var(self.name)

    @property
    def is_number(self):
        return False

class Dummy(Symbol):
    """Dummy Symbol

       use this through Symbol:

       >>> x1 = Symbol('x', dummy=True)
       >>> x2 = Symbol('x', dummy=True)
       >>> bool(x1 == x2)
       False
    """

    dummycount = 0

    __slots__ = ['dummy_index']

    def __new__(cls, name, commutative=True, **assumptions):
        obj = Symbol.__xnew__(cls, name, commutative=commutative, **assumptions)

        Dummy.dummycount += 1
        obj.dummy_index = Dummy.dummycount
        return obj

    def _hashable_content(self):
        return (self.name, self.dummy_index)

    def tostr(self, level=0):
        return '_' + self.name


class Temporary(Dummy):
    """
    Indexed dummy symbol.
    """

    __slots__ = []

    def __new__(cls, **assumptions):
        obj = Dummy.__new__(cls, 'T%i' % Dummy.dummycount, **assumptions)
        return obj


class Wild(Symbol):
    """
    Wild() matches any expression but another Wild().
    """

    __slots__ = ['exclude']

    def __new__(cls, name, exclude=None, **assumptions):
        if type(exclude) is list:
            exclude = tuple(exclude)

        return Wild.__xnew__(cls, name, exclude, **assumptions)

    @staticmethod
    @cacheit
    def __xnew__(cls, name, exclude, **assumptions):
        obj = Symbol.__xnew__(cls, name, **assumptions)

        if exclude is None:
            obj.exclude = None
        else:
            obj.exclude = tuple([sympify(x) for x in exclude])
        return obj

    def _hashable_content(self):
        return (self.name, self.exclude)

    # TODO add check against another Wild
    def matches(pattern, expr, repl_dict={}, evaluate=False):
        for p,v in repl_dict.items():
            if p==pattern:
                if v==expr: return repl_dict
                return None
        if pattern.exclude:
            for x in pattern.exclude:
                if x in expr:
                    return None
                #else:
                #    print x, expr, pattern, expr, pattern.exclude
        repl_dict = repl_dict.copy()
        repl_dict[pattern] = expr
        return repl_dict

    def __call__(self, *args, **assumptions):
        return WildFunction(self.name, nargs=len(args))(*args, **assumptions)

    def tostr(self, level=0):
        return self.name + '_'


def symbols(*names, **kwargs):
    """Returns a list of symbols with names taken from 'names'
       argument, which can be a string, then each character
       forms a separate symbol, or a sequence of strings.

       All newly created symbols have assumptions set accordingly
       to 'kwargs'. Main intention behind this function is to
       simplify and shorten examples code in doc-strings.

       >>> x, y, z = symbols('xyz')

       >>> a = symbols('a', integer=True)

       >>> a.is_integer
       True

       >>> xx, yy, zz = symbols('xx', 'yy', 'zz', real=True)

       >>> xx.is_real and yy.is_real and zz.is_real
       True

    """
    if len(names) == 1:
        result = [ Symbol(name, **kwargs) for name in names[0] ]
    else:
        result = [ Symbol(name, **kwargs) for name in names ]

    if len(result) == 1:
        return result[0]
    else:
        return result

def var(s):
    """
        Create a symbolic variable with the name *s*.

        INPUT:
            s -- a string, either a single variable name, or
                 a space separated list of variable names, or
                 a list of variable names.

        NOTE: The new variable is both returned and automatically injected into
        the parent's *global* namespace.  It's recommended not to use "var" in
        library code, it is better to use symbols() instead.

        EXAMPLES:
        We define some symbolic variables:
            >>> var('m')
            m
            >>> var('n xx yy zz')
            (n, xx, yy, zz)
            >>> n
            n

    """
    import inspect
    frame = inspect.currentframe().f_back

    try:
        if not isinstance(s, list):
            s = s.split(" ")

        res = []

        for t in s:
            # skip empty strings
            if not t:
                continue
            sym = Symbol(t)
            frame.f_globals[t] = sym
            res.append(sym)

        res = tuple(res)
        if len(res) == 0:   # var('')
            res = None
        elif len(res) == 1: # var('x')
            res = res[0]
                            # otherwise var('a b ...')
        return res

    finally:
        # we should explicitly break cyclic dependencies as stated in inspect
        # doc
        del frame

# /cyclic/
import basic as _
_.Symbol    = Symbol
_.Wild      = Wild
_.Temporary = Temporary
del _
