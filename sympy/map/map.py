from sympy.assumptions import ask, Q
from sympy.core import Expr, S
from sympy.core.decorators import (
    call_highest_priority, sympify_method_args, sympify_return
)
from sympy.core.sympify import _sympify
from sympy.core.symbol import Str
from sympy.sets import ProductSet

__all__ = [
    'Map', 'UndefinedMap', 'InverseMap', 'IdentityMap', 'AppliedMap'
]

@sympify_method_args
class Map(Expr):
    """
    Abstract base class for general mathematical mappings, which also serves
    as a constructor for undefined map classes.

    Explanation
    ===========

    Map is a binary relation over two sets that associates to every
    element of the first set exactly one element of the second set [1].
    The first set is called domain, and the second set is called
    codomain of the map.

    .. note::
       User must be aware that argument of `f(x)` can be interpreted to
       both `x` or `(x,)`. The former agrees with mathematical intuition,
       while the latter is mathematically rigorous (just as the argument
       of `f(x,y)` is `(x, y)`).

    Examples
    ========

    >>> from sympy import symbols, S
    >>> from sympy.map import Map
    >>> x, y = symbols('x y')

    Constructing undefined map:

    >>> f = Map('f', S.Reals, S.Reals**2)

    >>> f
    f
    >>> f.domain
    Reals
    >>> f.codomain
    ProductSet(Reals, Reals)
    >>> f.nargs
    1

    >>> f(x)
    f(x)
    >>> f(x).doit()
    f(x)

    Subclassing for defined map:

    >>> class G(Map):
    ...     name = 'g'
    ...     domain = S.Reals**2
    ...     codomain = S.Reals
    ...     def eval(self, tup):
    ...         return sum(tup)
    >>> g = G()

    >>> g.nargs
    2
    >>> g((x,y)).doit()
    x + y

    References
    ==========

    .. [1] https://en.wikipedia.org/wiki/Function_(mathematics)

    """

    # these attributes are designed to be overridden if needed
    domain = codomain = S.UniversalSet

    def __new__(cls, *args, **kwargs):
        if cls is Map:
            return UndefinedMap(*args, **kwargs)
        return super().__new__(cls, *args)

    @property
    def nargs(self):
        if isinstance(self.domain, ProductSet):
            return len(self.domain.args)
        else:
            return 1

    def eval(self, *args):
        return

    def __call__(self, *args, evaluate=False, **kwargs):
        return AppliedMap(self, *args, evaluate=evaluate)

    @sympify_return([('other', Expr)], NotImplemented)
    @call_highest_priority('__rmatmul__')
    def __matmul__(self, other):
        return CompositeMap(self, other, evaluate=True)

    @sympify_return([('other', Expr)], NotImplemented)
    @call_highest_priority('__matmul__')
    def __rmatmul__(self, other):
        return CompositeMap(other, self, evaluate=True)

    def inverse(self, evaluate=False):
        return InverseMap(self, evaluate=evaluate)
    inv = inverse

    def _eval_inverse(self):
        return

    def composite(self, other, evaluate=False):
        return CompositeMap(self, other, evaluate=evaluate)

    def _eval_composite(self, other):
        # define composition with other 'special' map here.
        # composition with inverse map or identity map need not
        # be defined here since CompositeMap deals with it.
        return

    def _eval_iteration(self, n):
        # define special result for positive n-th iteration
        # of self here.
        return

    def doit(self, **hints):
        deep = hints.get('deep', True)
        if deep:
            args = (a.doit(**hints) for a in self.args)
            return self.func(*args, evaluate=True)
        else:
            return self.func(*self.args, evaluate=True)

    def as_base_iternum(self):
        """
        Interprete *self* as IteratedMap, returning the
        base and number of iteration.

        """
        return self, S.One

class UndefinedMap(Map):
    """
    A class for undefined mappings.

    Parameters
    ==========

    name : str
        Name of the map

    domain, codomain : Set

    """

    def __new__(
        cls, name, domain=S.UniversalSet, codomain=S.UniversalSet,
        **kwargs
    ):
        if not isinstance(name, Str):
            name = Str(name)

        domain, codomain = _sympify(domain), _sympify(codomain)

        obj = super().__new__(cls, name, domain, codomain)
        return obj

    @property
    def name(self):
        return self.args[0]

    @property
    def domain(self):
        return self.args[1]

    @property
    def codomain(self):
        return self.args[2]

class InverseMap(Map):
    """
    A class for unevaluated inverse mappings.

    Parameters
    ==========

    mapping : Map

    evaluate : bool, optional
        If True, return the evaluated inverse function of *mapping*

    Example
    =======

    >>> from sympy.map import Map, InverseMap

    >>> class Exp(Map):
    ...     def _eval_inverse(self):
    ...         return Log()
    >>> exp = Exp()

    >>> class Log(Map):
    ...     def _eval_inverse(self):
    ...         return Exp()
    >>> log = Log()

    >>> exp.inv() == InverseMap(exp)
    True
    >>> exp.inv(evaluate=True) == InverseMap(exp, evaluate=True) == log
    True

    """
    def __new__(cls, mapping, evaluate=False, **kwargs):

        if ask(Q.invertible(mapping)) is False:
            raise TypeError("%s is not invertible." % mapping)

        if evaluate:
            obj = mapping._eval_inverse()
            if obj is not None:
                return obj
        return super().__new__(cls, mapping)

    @property
    def base(self):
        return self.args[0]

    @property
    def domain(self):
        return self.base.codomain

    @property
    def codomain(self):
        return self.base.domain

    def _eval_inverse(self):
        return self.base

class IdentityMap(Map):
    """
    General identity mapping.

    .. note::
       IdentityMap is unary. Pass tuple for n-dimensional argument.

    Explanation
    ===========

    Identity map is a map that always return its argument.

    Parameters
    ==========

    domain : Set

    Examples
    ========

    >>> from sympy.abc import x, y
    >>> from sympy.map import IdentityMap

    >>> I = IdentityMap()
    >>> I(x, evaluate=True)
    x
    >>> I((x, y), evaluate=True)
    (x, y)

    """
    def __new__(cls, domain=S.UniversalSet, **kwargs):
        return super().__new__(cls, domain)

    @property
    def domain(self):
        return self.args[0]

    @property
    def codomain(self):
        return self.domain

    def __call__(self, arg, **kwargs):
        return super().__call__(arg, **kwargs)

    def eval(self, x):
        return x

    def _eval_inverse(self):
        return self

class AppliedMap(Expr):
    """
    Unevaluated result of Map applied to arguments.

    .. note::
       If user include parameter in mapping definition (e.g. $a$ in $f(x) = a*x)$,
       `free_symbol` of its `AppliedMap` will return incoincident value with its
       evaluated result. Instead, design the mapping to be multivariate function
       (e.g. $f(x, a) = a*x$).
       User should not call this class directly. Instead, use `__call__` method
       of the map instance.

    Parameters
    ==========

    map : Map

    args
        arguments applied to *map*

    evaluate : bool, optional
        If True, returns evaluated application of *map* to *args*

    Examples
    ========

    >>> from sympy.abc import x
    >>> from sympy.map import Map, AppliedMap

    >>> class F(Map):
    ...     def eval(self, x):
    ...         return x + 2
    >>> f = F()

    >>> isinstance(f(x), AppliedMap)
    True
    >>> f(x, evaluate=True)
    x + 2

    """
    def __new__(cls, mapping, *args, evaluate=False, **kwargs):
        args = [_sympify(a) for a in args]

        if evaluate:

            # convert f(f(x)) to CompositeMap(f,f)(x)
            if len(args) == 1 and isinstance(args[0], cls):
                mapping = CompositeMap(mapping, args[0].map, evaluate=True)
                args = args[0].arguments

            result = mapping.eval(*args)
            if result is not None:
                return result

        return super().__new__(cls, mapping, *args)

    @property
    def map(self):
        return self.args[0]

    @property
    def arguments(self):
        return self.args[1:]

    def doit(self, **hints):
        deep = hints.get('deep', True)
        if deep:
            args = (a.doit(**hints) for a in self.args)
            return self.func(*args, evaluate=True)
        else:
            return self.func(*self.args, evaluate=True)

from .composite import CompositeMap, IteratedMap
