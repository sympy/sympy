from sympy.assumptions import ask, Q
from sympy.core import Expr, S, Tuple
from sympy.core.compatibility import iterable
from sympy.core.decorators import (
    call_highest_priority, sympify_method_args, sympify_return
)
from sympy.core.sympify import _sympify
from sympy.core.symbol import Str
from sympy.sets import ProductSet, FiniteSet

__all__ = [
    'Map', 'UndefinedMap', 'RestrictedMap', 'InverseMap', 'IdentityMap',
    'AppliedMap',
    'isappliedmap'
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
       User must be aware that argument of ``f(x)`` can be interpreted to
       both ``x`` or ``(x,)``. The former agrees with mathematical intuition,
       while the latter is mathematically rigorous (just as the argument
       of ``f(x,y)`` is ``(x, y)``).
       Due to the catch mentioned above, domain of the argument is not
       checked by default when map is applied to them.

    Examples
    ========

    >>> from sympy import symbols, S
    >>> from sympy.map import Map
    >>> x, y = symbols('x y', real=True)

    Constructing undefined map:

    >>> f = Map('f', S.Reals, S.Reals**2)

    >>> f
    f : Reals -> ProductSet(Reals, Reals)
    >>> f.domain
    Reals
    >>> f.codomain
    ProductSet(Reals, Reals)
    >>> f.arity
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
    ...     def eval(self, *args):
    ...         return sum(args)
    >>> g = G()

    >>> g.arity
    2
    >>> g(x,y).doit()
    x + y

    References
    ==========

    .. [1] https://en.wikipedia.org/wiki/Function_(mathematics)

    """

    # these attributes are designed to be overridden if needed
    # they can be overridden by both class attribute or instance attribute
    domain = codomain = S.UniversalSet

    def __new__(cls, *args, **kwargs):
        if cls is Map:
            return UndefinedMap(*args, **kwargs)
        return super().__new__(cls, *args)

    @property
    def arity(self):
        if self.domain in (FiniteSet(()), FiniteSet(S.EmptySet)):
            return 0
        if isinstance(self.domain, ProductSet):
            return len(self.domain.args)
        return 1
    nargs = arity

    def apply(self, *args, **kwargs):
        """
        Action of the operator on arguments. This method is always
        called when arguments are applied.

        Examples
        ========

        >>> from sympy import Map
        >>> class F(Map):
        ...     name = 'f'
        ...     def apply(self, x, **kwargs):
        ...         return x + 1
        >>> f = F()
    
        >>> f(1)
        2

        """
        evaluate = kwargs.get('evaluate', False)
        if evaluate:
            return self.eval(*args)        

    def eval(self, *args):
        """
        Evaluation of the operator on its arguments. This method
        is called when ``evaluate=True`` is passed.

        Examples
        ========

        >>> from sympy import Map
        >>> class F(Map):
        ...     name = 'f'
        ...     def eval(self, x):
        ...         return x + 1
        >>> f = F()

        >>> f(1)
        f(1)
        >>> f(1, evaluate=True)
        2

        """
        return

    def __call__(self, *args, evaluate=False, **kwargs):
        return AppliedMap(self, args, evaluate=evaluate, **kwargs)

    @sympify_return([('other', Expr)], NotImplemented)
    @call_highest_priority('__rmatmul__')
    def __matmul__(self, other):
        return composite_op(self, other, evaluate=True)

    @sympify_return([('other', Expr)], NotImplemented)
    @call_highest_priority('__matmul__')
    def __rmatmul__(self, other):
        return composite_op(other, self, evaluate=True)

    def inverse(self, evaluate=False):
        """
        Returns the inverse function of *self*.

        """
        return InverseMap(self, evaluate=evaluate)
    inv = inverse

    def iterate(self, n, evaluate=False):
        iter_op = composite_op.exponent_operator()
        return iter_op(self, n, evaluate=evaluate)

    def _eval_inverse(self):
        return

    def composite(self, other, evaluate=False):
        return composite_op(self, other, evaluate=evaluate)

    def _eval_composite(self, other):
        # define composition with other 'special' map here.
        # composition with inverse map or identity map need not
        # be defined here since CompositeMap deals with it.
        return

    def _eval_iterate(self, n):
        # define special result for positive n-th iteration
        # of self here.
        return

    def restrict(self, new_domain, evaluate=False):
        """
        Restrict the map to smaller domain.

        """
        return RestrictedMap(self, new_domain, evaluate=evaluate)

    def _eval_restrict(self, new_domain):
        return

    def is_restriction(self, other):
        """
        Returns True if *self* is restricted function of *other*.

        Explanation
        ===========

        Restricted function has smaller domain, but pertains the relation between
        domain and codomain. For the identification, ``_map_content`` method is
        referred to.

        Examples
        ========

        >>> from sympy import Map, S
        >>> class F(Map):
        ...     @property
        ...     def domain(self):
        ...         return self.args[0]
        ...     @property
        ...     def codomain(self):
        ...         return self.args[1]
        ...     def eval(self, x):
        ...         return x + 1
        >>> f1, f2 = F(S.Reals, S.Reals), F(S.Integers, S.Reals)

        >>> f1 == f2
        False
        >>> f2.is_restriction(f1)
        True
        >>> f1(1) == f2(1)
        False
        >>> f1(1, evaluate=True) == f2(1, evaluate=True)
        True

        """
        return self._eval_is_restriction(other)

    def _eval_is_restriction(self, other):
        return (
            self._map_content() == other._map_content() and
            self.arity == other.arity and
            self.codomain.is_subset(other.codomain) and
            self.domain.is_subset(other.domain)
        )

    def _map_content(self):
        # Unique fingerprint of map, independent of domain and codomain.
        # Used for restricted function identification
        return self.func

    def doit(self, **hints):
        deep = hints.get('deep', True)
        if deep:
            args = (a.doit(**hints) for a in self.args)
            return self.func(*args, evaluate=True)
        else:
            return self.func(*self.args, evaluate=True)

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

        obj = super().__new__(cls, name, domain, codomain, **kwargs)
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

    def _map_content(self):
        return self.func, self.name

class RestrictedMap(Map):
    r"""
    A class for general restricted function, whose domain is restricted to a subset of its
    original domain.

    Parameters
    ==========

    mapping : Map
        Original map before restriction.

    new_domain : Set
        Restricted domain.

    evaluate : bool, optional
        If True, return the evaluated restriction of *mapping*.

    Examples
    ========

    >>> from sympy import Map, S
    >>> f = Map('f', S.Reals, S.Integers)
    >>> f
    f : Reals -> Integers
    >>> f.restrict(S.Integers)
    RestrictedMap(f, Integers) : Integers -> Integers

    """
    def __new__(cls, mapping, new_domain, evaluate=False, **kwargs):

        new_domain = _sympify(new_domain)
        if not new_domain.is_subset(mapping.domain):
            raise TypeError("%s is not subset of %s's domain %s." % (new_domain, mapping, mapping.domain))

        if evaluate:

            if isinstance(mapping, RestrictedMap):
                return cls(mapping.base, new_domain, evaluate=True)

            result = mapping._eval_restrict(new_domain)
            if result is not None:
                return result
        return super().__new__(cls, mapping, new_domain, **kwargs)

    @property
    def base(self):
        return self.args[0]

    @property
    def domain(self):
        return self.args[1]

    @property
    def codomain(self):
        return self.base.codomain

    def eval(self, *args):
        return self.base(*args, evaluate=True)

    def _map_content(self):
        return self.base._map_content()

class InverseMap(Map):
    """
    A class for unevaluated inverse mappings.

    Parameters
    ==========

    mapping : Map
        Original map which will be inversed.

    evaluate : bool, optional
        If True, return the evaluated inverse function of *mapping*.

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
        return super().__new__(cls, mapping, **kwargs)

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

    def as_base_iternum(self):
        return self.base, S.NegativeOne

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
        return super().__new__(cls, domain, **kwargs)

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

    def _map_content(self):
        return self.func

class AppliedMap(Expr):
    """
    Unevaluated result of Map applied to arguments.

    .. note::
       If user include parameter in mapping definition (e.g. $a$ in $f(x) = a*x)$,
       ``free_symbol`` of its ``AppliedMap`` will return incoincident value with its
       evaluated result. Instead, design the mapping to be multivariate function
       (e.g. $f(x, a) = a*x$).
       This class is not designed to be constructed directly. Instead, use ``__call__`` method
       of the map instance.

    Parameters
    ==========

    map : Map

    args : tuple of arguments
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
    def __new__(cls, mapping, args, evaluate=False, **kwargs):
        kwargs.update(evaluate=evaluate)

        # consult mapping.apply
        result = mapping.apply(*args, **kwargs)

        if result is None:
            # generate AppliedMap class
            args = Tuple(*[_sympify(a) for a in args])
            result = super().__new__(cls, mapping, args)

        # check codomain
        if mapping.codomain.contains(result) == False:
                    raise TypeError(
                "%s is not in %s's codomain %s." % (result, mapping, mapping.codomain)
                )

        return result

    @property
    def map(self):
        return self.args[0]

    @property
    def arguments(self):
        return self.args[1]

    def _new_rawargs(self, *args, **kwargs):
        kwargs["evaluate"] = False
        return self.func(self.map, args, **kwargs)

    def _contained(self, other):
        # Let `f(x) in f.codomain` return True
        # see Set.contains and Set.__contains__ methods.
        if other.is_superset(self.map.codomain):
            return True
        # Do not return False; allow other to check as well.

    def doit(self, **hints):
        deep = hints.get('deep', True)
        if deep:
            args = (a.doit(**hints) for a in self.args)
            return self.func(*args, evaluate=True)
        else:
            return self.func(*self.args, evaluate=True)

    def as_base_exp(self, operator):
        if self.map.is_restriction(operator):
            return self.map._eval_as_base_exp(*self.arguments)
        return self, S.One

def isappliedmap(arg, maps):
    """
    Return ``True`` if *arg* is unevaluated applied result of
    *maps* (if *maps* is ``Map``), or one of *maps* (if *maps* is iterable)

    Parameters
    ==========

    arg : Expr

    maps : Map, or iterable of Map

    Examples
    ========

    >>> from sympy import S, Map, isappliedmap
    >>> f = Map('f')
    >>> g = Map('g')

    >>> isappliedmap(f(1), f)
    True
    >>> isappliedmap(f(1), g)
    False
    >>> isappliedmap(f(1), (f, g))
    True

    """
    if isinstance(arg, AppliedMap):
        arg_map = arg.map
        if not iterable(maps):
            # maps is map
            return arg.map.is_restriction(maps)
        else:
            for m in maps:
                if arg.map.is_restriction(m):
                    return True
    return False

from .composite import function_set, composite_op
