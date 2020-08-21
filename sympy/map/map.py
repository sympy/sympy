from functools import cmp_to_key
from sympy.assumptions import ask, Q
from sympy.core import Add, S, Basic, Expr, Tuple
from sympy.core.decorators import call_highest_priority, sympify_method_args
from sympy.core.symbol import Str, Dummy
from sympy.core.sympify import _sympify
from sympy.core.compatibility import iterable
from sympy.sets import FiniteSet, ProductSet, Interval

__all__ = [
    'Map', 'UndefinedMap',
    'RestrictedMap', 'InverseMap', 'IdentityMap', 'ConstantMap',
    'AppliedMap', 'isappliedmap',
]

@sympify_method_args
class Map(Basic):
    """
    Abstract base class for general mathematical maps, which also serves
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

    >>> from sympy import symbols, S, Add
    >>> from sympy.map import Map
    >>> x, y = symbols('x y', real=True)

    Constructing undefined map:

    >>> f = Map('f', S.Reals**2, S.Reals)

    >>> f
    f : ProductSet(Reals, Reals) -> Reals
    >>> f(x)
    f(x)
    >>> f.domain
    ProductSet(Reals, Reals)
    >>> f.codomain
    Reals

    Since domain of ``f`` is binary cartesian product, its arity is 2.

    >>> f.arity
    2

    Subclassing for defined map:

    >>> class G(Map):
    ...     name = 'g'
    ...     domain = S.Reals**2
    ...     codomain = S.Reals
    ...     def eval(self, *args):
    ...         return Add(*args, evaluate=True)
    >>> g = G()

    >>> g.arity
    2
    >>> g(x,y).doit()
    x + y

    References
    ==========

    .. [1] https://en.wikipedia.org/wiki/Function_(mathematics)

    """

    _op_priority = 11

    # since sympy's is_commutative stands for the element being commutative to multiplication,
    # we use these attributes for the operator being commutative and associative.
    is_commutative_map = None
    is_associative = None
    is_invertible = None

    # these attributes are designed to be overridden if needed
    # they can be overridden by both class attribute or instance attribute
    domain = codomain = S.UniversalSet

    def _corresponding_oldfunc(self):
        """
        Introduced for compatibility with Function

        See Also
        ========

        AppliedMap._allowed_superclasshook
        """
        return [Function]

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

    @property
    def range(self):
        result = self._eval_range()
        if result is not None:
            return result
        return self.codomain

    def _eval_range(self):
        return

    def doit(self, **hints):
        deep = hints.get('deep', True)
        if deep:
            args = (a.doit(**hints) for a in self.args)
            return self.func(*args, evaluate=True)
        else:
            return self.func(*self.args, evaluate=True)

    ### AppliedMap related

    def __call__(self, *args, evaluate=False, **kwargs):
        return AppliedMap(self, args, evaluate=evaluate, **kwargs)

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
        evaluate = kwargs.pop('evaluate', False)
        if evaluate:
            return self.eval(*args, **kwargs)

    def eval(self, *args, **kwargs):
        """
        Evaluation of the operator on its arguments. This method
        is called by ``apply`` when ``evaluate=True`` is passed.

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

    def process_args(self, args, **kwargs):
        """
        Process the arguments, i.e. flattening or commutating. This
        method is run when the function is not evaluated by ``apply``
        and ``eval``.

        """
        if ask(Q.associative(self)):
            args = self.flatten(args)
        if ask(Q.commutative_map(self)):
            args.sort(key=cmp_to_key(Basic.compare))
        return args

    def flatten(self, seq):
        new_seq = []
        for o in seq:
            if isappliedmap(o, self):
                new_seq.extend(o.arguments)
            else:
                new_seq.append(o)
        return new_seq

    ### InverseMap related

    def inverse(self, evaluate=False):
        """
        Returns the inverse function of *self*.

        """
        return InverseMap(self, evaluate=evaluate)
    inv = inverse

    def _eval_inverse(self):
        return

    ### RestrictedMap related

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

    ### mapop related

    def composite(self, other, **kwargs):
        return CompositeMap(self, other, **kwargs)

    def _eval_composite(self, other):
        return

    def iterate(self, n, **kwargs):
        return IteratedMap(self, n, **kwargs)

    def _eval_iterate(self, n):
        return

    def as_base_iternum(self):
        return self, S.One

    def as_coeff_Mul(self, *args, **kwargs):
        return S.One, self

    def as_base_exp(self):
        return self, S.One

    def _eval_power(self, other):
        return

    ### DerivativeFunction related

    def diff(self, *indices, evaluate=True):
        diffop = DiffOp(*indices)
        return diffop(self, evaluate=evaluate)

    def _eval_derivative_n_times(self, index, count):
        from sympy import Integer
        if isinstance(count, (int, Integer)):
            obj = self
            for i in range(count):
                obj2 = obj.fdiff(index)
                if obj == obj2 or obj2 is None:
                    break
                obj = obj2
            return obj2
        else:
            return None

    def fdiff(self, i=1):
        return None

    ### assumptions of AppliedMap of self

    def _applied_is_commutative(self, expr):
        return

    def _applied_is_rational(self, expr):
        return self.codomain.is_subset(S.Rationals)

    def _applied_is_algebraic(self, expr):
        return

    ### operator overriding

    def __pos__(self):
        return self

    def __neg__(self):
        return MapMul(S.NegativeOne, self, evaluate=True)

    @call_highest_priority('__radd__')
    def __add__(self, other):
        return MapAdd(self, other, evaluate=True)

    @call_highest_priority('__add__')
    def __radd__(self, other):
        return MapAdd(other, self, evaluate=True)

    @call_highest_priority('__rsub__')
    def __sub__(self, other):
        return MapAdd(self, -other, evaluate=True)

    @call_highest_priority('__sub__')
    def __rsub__(self, other):
        return MapAdd(other, -self, evaluate=True)

    @call_highest_priority('__rmul__')
    def __mul__(self, other):
        return MapMul(self, other, evaluate=True)

    @call_highest_priority('__mul__')
    def __rmul__(self, other):
        return MapMul(other, self, evaluate=True)

    @call_highest_priority('__rpow__')
    def __pow__(self, other):
        return MapPow(self, other, evaluate=True)

    @call_highest_priority('__pow__')
    def __rpow__(self, other):
        return MapPow(other, self, evaluate=True)

    @call_highest_priority('__rdiv__')
    def __div__(self, other):
        denom = MapPow(other, S.NegativeOne, evaluate=True)
        if isinstance(self, ConstantMap) and self.output is S.One:
            return denom
        else:
            return MapMul(self, denom, evaluate=True)
    __truediv__ = __div__

    @call_highest_priority('__div__')
    def __rdiv__(self, other):
        denom = MapPow(self, S.NegativeOne, evaluate=True)
        if isinstance(other, ConstantMap) and other.output is S.One:
            return denom
        else:
            return MapMul(denom, other, evaluate=True)
    __rtruediv__ = __rdiv__

    @call_highest_priority('__rmatmul__')
    def __matmul__(self, other):
        return CompositeMap(self, other, evaluate=True)

    @call_highest_priority('__matmul__')
    def __rmatmul__(self, other):
        return CompositeMap(other, self, evaluate=True)

class UndefinedMap(Map):
    """
    A class for undefined maps.

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

    base : Map
        Original map before restriction.

    new_domain : Set
        Restricted domain.

    evaluate : bool, optional
        If True, return the evaluated restriction of *base*.

    Examples
    ========

    >>> from sympy import Map, S
    >>> f = Map('f', S.Reals, S.Integers)
    >>> f
    f : Reals -> Integers
    >>> f.restrict(S.Integers)
    RestrictedMap(f, Integers) : Integers -> Integers

    """

    @property
    def is_commutative_map(self):
        return ask(Q.commutative_map(self.base))

    @property
    def is_associative(self):
        return ask(Q.associative(self.base))

    @property
    def invertible(self):
        return ask(Q.invertible(self.base))

    def __new__(cls, base, new_domain, evaluate=False, **kwargs):

        new_domain = _sympify(new_domain)
        if not new_domain.is_subset(base.domain):
            raise TypeError("%s is not subset of %s's domain %s." % (new_domain, base, base.domain))

        if evaluate:

            if isinstance(base, RestrictedMap):
                return cls(base.base, new_domain, evaluate=True)

            result = base._eval_restrict(new_domain)
            if result is not None:
                return result
        return super().__new__(cls, base, new_domain, **kwargs)

    @property
    def base(self):
        return self.args[0]

    @property
    def domain(self):
        return self.args[1]

    @property
    def codomain(self):
        return self.base.codomain

    def eval(self, *args, **kwargs):
        kwargs.update(evaluate=True)
        return self.base(*args, **kwargs)

    def _map_content(self):
        return self.base._map_content()

class InverseMap(Map):
    """
    A class for unevaluated inverse maps.

    Parameters
    ==========

    base : Map
        Original map which will be inversed.

    evaluate : bool, optional
        If True, return the evaluated inverse function of *base*.

    Example
    =======

    >>> from sympy.map import Map, InverseMap

    >>> class F(Map):
    ...     def _eval_inverse(self):
    ...         return G()
    >>> f = F()

    >>> class G(Map):
    ...     def _eval_inverse(self):
    ...         return F()
    >>> g = G()

    >>> f.inv() == InverseMap(f)
    True
    >>> f.inv(evaluate=True) == InverseMap(f, evaluate=True) == g
    True

    """
    # invertible function can be always inverted back
    is_invertible = True

    def __new__(cls, base, evaluate=False, **kwargs):

        if ask(Q.invertible(base)) is False:
            raise TypeError("%s is not invertible." % base)

        if evaluate:
            obj = base._eval_inverse()
            if obj is not None:
                return obj
        return super().__new__(cls, base, **kwargs)

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
    General identity map.

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

    is_commutative_map = False
    is_associative = False
    is_invertible = True

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

    def eval(self, x, **kwargs):
        return x

    def _eval_inverse(self):
        return self

    def fdiff(self, i=1):
        return ConstantMap(1, self.domain)

class ConstantMap(Map):
    """
    Map which always returns same output.

    Parameters
    ==========

    output : Output of the constant map

    domain : Set, optional
        Domain of the map. Default is ``S.UniversalSet``.

    Examples
    ========

    >>> from sympy import S, ConstantMap, Symbol
    >>> x = Symbol('x', real=True)

    >>> f = ConstantMap(1, S.Reals)
    >>> f
    1 : Reals -> FiniteSet(1)

    >>> f(x)
    1(x)
    >>> f(x, evaluate=True)
    1

    """

    is_commutative_map = True
    is_associative = True
    is_invertible = False

    def __new__(cls, output, domain=S.UniversalSet, **kwargs):
        output = _sympify(output)
        domain = _sympify(domain)
        return super().__new__(cls, output, domain)

    @property
    def output(self):
        return self.args[0]

    @property
    def domain(self):
        return self.args[1]

    @property
    def codomain(self):
        return FiniteSet(self.output)

    def eval(self, *args, **kwargs):
        return self.output

    def _map_content(self):
        return self.func, self.output

    def fdiff(self, i=1):
        return self.func(0, self.domain)

class AppliedMap(Expr):
    """
    Unevaluated result of Map applied to arguments.

    .. note::
       If user include parameter in map definition (e.g. $a$ in $f(x) = a*x)$,
       ``free_symbol`` of its ``AppliedMap`` will return incoincident value with its
       evaluated result. Instead, design the map to be multivariate function
       (e.g. $f(x, a) = a*x$).
       This class is not designed to be constructed directly. Instead, use ``__call__`` method
       of the map instance.

    Parameters
    ==========

    map : Map

    args : tuple of arguments
        Arguments applied to *map*

    evaluate : bool, optional
        If True, returns evaluated application of *map* to *args*.
        Default is False.

    Examples
    ========

    >>> from sympy.abc import x
    >>> from sympy.map import Map, AppliedMap

    >>> class F(Map):
    ...     def eval(self, *args):
    ...         return x + 2
    >>> f = F()

    >>> isinstance(f(x), AppliedMap)
    True
    >>> f(x, evaluate=True)
    x + 2

    """

    def _allowed_superclasshook(self):
        """
        Introduced for compatibility with Function

        Examples
        ========

        >>> from sympy import S, Sin, sin
        >>> from sympy.testing.pytest import warns_deprecated_sympy
        >>> from sympy.abc import x

        >>> isinstance(Sin(S.Reals)(x), sin)
        True

        """
        return self.map._corresponding_oldfunc()

    def __new__(cls, map, args, evaluate=False, **kwargs):
        kwargs.update(evaluate=evaluate)
        args = [_sympify(a) for a in args]

        # consult map.apply
        result = map.apply(*args, **kwargs)

        if result is None:

            # Even function cannot be evaluated, arguments still
            # might can be processed. (i.e. flattening, commutating..)
            if evaluate:
                args = map.process_args(args, **kwargs)

            # generate AppliedMap class
            args = Tuple(*args)
            result = super().__new__(cls, map, args)

            # Compatibility with old function
            result._args = _DeprecatedArgs(map, args)

        # check codomain
        if map.codomain.contains(result) == False:
                    raise TypeError(
                "%s is not in %s's codomain %s." % (result, map, map.codomain)
                )

        return result

    @property
    def map(self):
        # When _DeprecatedArgs is removed, change this to self.args[0]
        return self.args.args[0]

    @property
    def arguments(self):
        return self.args.args[1]

    @property
    def free_symbols(self):
        return self.arguments.free_symbols

    def denest(self, deep=True, **kwargs):
        """
        Flatten the nested structure, but do not evaluate it

        Examples
        ========

        >>> from sympy import Map, Add
        >>> class F(Map):
        ...     name = 'f'
        ...     is_associative = True
        ...     def eval(self, *args):
        ...         return Add(*args, evaluate=True)
        >>> f = F()
        >>> expr = f(1, f(2, f(3, 4)))

        >>> expr
        f(1, f(2, f(3, 4)))
        >>> expr.denest(deep=False)
        f(1, 2, f(3, 4))
        >>> expr.denest()
        f(1, 2, 3, 4)

        """
        if deep:
            arguments = []
            for a in self.arguments:
                if hasattr(a, 'denest'):
                    arguments.append(a.denest(deep=True, **kwargs))
                else:
                    arguments.append(a)
        else:
            arguments = self.arguments

        if ask(Q.associative(self.map)):
            arguments = self.map.flatten(arguments)

        return self._new_rawargs(*arguments)

    def _new_rawargs(self, *args, **kwargs):
        return self.func(self.map, args, **kwargs)

    def doit(self, **hints):
        """
        Evaluate *self*.

        Examples
        ========

        >>> from sympy import Map, Add
        >>> class F(Map):
        ...     name = 'f'
        ...     is_associative = True
        ...     def eval(self, *args):
        ...         return Add(*args, evaluate=True)
        >>> f = F()
        >>> expr = f(1, 2, f(3, f(4, 5)))

        >>> expr
        f(1, 2, f(3, f(4, 5)))
        >>> expr.doit(deep=False)
        3 + f(3, f(4, 5))
        >>> expr.doit()
        15

        """
        deep = hints.get('deep', True)
        if deep:
            self = self.denest(**hints)
            args = [a.doit(**hints) for a in self.args]
            return self.func(*args, evaluate=True)
        else:
            return self.func(*self.args, evaluate=True)

    def _eval_derivative(self, s):
        # f(x).diff(s) -> x.diff(s) * f.diff(1)(s)
        i = 0
        l = []
        for a in self.arguments:
            i += 1
            da = a.diff(s)
            if da.is_zero:
                continue
            deriv_func = self.map.diff(i)
            df = deriv_func(*self.arguments, evaluate=True)
            l.append(df * da)
        return Add(*l)

    ### assumptions

    def _eval_is_commutative(self):
        return self.map._applied_is_commutative(self)

    def _eval_is_rational(self):
        return self.map._applied_is_rational(self)

    def _eval_is_algebraic(self):
        return self.map._applied_is_algebraic(self)

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

    >>> from sympy import Map, isappliedmap
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
        if not iterable(maps):
            # maps is map
            return arg.map.is_restriction(maps)
        else:
            for m in maps:
                if arg.map.is_restriction(m):
                    return True
    return False

from .composite import CompositeMap, IteratedMap
from .mapop import MapAdd, MapMul, MapPow
from .derivative import DiffOp, DerivativeFunction

### Special container for compatibility

from sympy.utilities.exceptions import SymPyDeprecationWarning

class _DeprecatedArgs(Tuple):
    """
    In old function, arguments are stored in *args* and in many parts of sympy
    *args* used to access them. However, in ``AppliedMap``, *args* contains maps
    and arguments. It is *arguments* which really stores the arguments that are
    applied to the map.

    For backward compatibility, this class is introduced to match the behavior
    of two objects.

    Examples
    ========

    >>> from sympy.testing.pytest import warns_deprecated_sympy
    >>> from sympy import S, Sin, sin as old_sin
    >>> from sympy.abc import x
    >>> new_sin = Sin(S.Reals)

    >>> old_sin(x).args
    (x,)

    >>> new_sin(x).args
    (sin : Reals -> Interval(-1, 1), (x,))
    >>> new_sin(x).arguments
    (x,)

    >>> with warns_deprecated_sympy():
    ...     print(old_sin(x).args[0] == new_sin(x).args[0] == x)
    True

    """
    def __getitem__(self, index):
        SymPyDeprecationWarning(
            feature="Accessing arguments of AppliedMap with `args`",
            useinstead="arguments",
            issue=19953,
            deprecated_since_version="1.7"
        ).warn()
        return self.args[1][index]

### imported for Function.__instancecheck__

from sympy.core.function import ArgumentIndexError, Function, Derivative, Subs
