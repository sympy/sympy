from collections import defaultdict
from collections.abc import MutableMapping

from sympy.assumptions.ask import Q
from sympy.core import (Add, Mul, Pow, Integer, Number, NumberSymbol, Symbol)
from sympy.core.compatibility import iterable
from sympy.core.numbers import ImaginaryUnit
from sympy.functions.elementary.complexes import Abs
from sympy.logic.boolalg import (Boolean, Equivalent, And, Or, Implies)
from sympy.matrices.expressions import MatMul

# APIs here may be subject to change


### Helper classes ###

class ArgFactHandler:
    """
    Class to help apply boolean function to the arguments of expression.

    Parameters
    ==========

    args : tuple of Symbols.
        Placeholder symbols representing the arguments of *f*.

    f : Any Boolean expression.

    """
    def __init__(self, args, f):
        self.args = args
        self.f = f

    def __repr__(self):
        return "%s(%s, %s)" % (type(self).__name__, str(self.args), str(self.f))

    def apply(self, args):
        return self.f.subs(zip(self.args, args))


class AllArgs(ArgFactHandler):
    """
    Apply the function to all arguments of the expression.

    Examples
    ========

    >>> from sympy import Q
    >>> from sympy.assumptions.sathandlers import AllArgs
    >>> from sympy.abc import x, y
    >>> allargs = AllArgs(x, Q.negative(x) | Q.positive(x))
    >>> allargs(x*y)
    (Q.negative(x) | Q.positive(x)) & (Q.negative(y) | Q.positive(y))

    """
    def __init__(self, x, f):
        args = (x,)
        super().__init__(args, f)

    def __call__(self, expr):
        return And(*[self.apply((arg,)) for arg in expr.args])


class AnyArgs(ArgFactHandler):
    """
    Apply the function to any argument of the expression.

    Examples
    ========

    >>> from sympy import Q
    >>> from sympy.assumptions.sathandlers import AnyArgs
    >>> from sympy.abc import x, y
    >>> anyargs = AnyArgs(x, Q.negative(x) & Q.positive(x))
    >>> anyargs(x*y)
    (Q.negative(x) & Q.positive(x)) | (Q.negative(y) & Q.positive(y))

    """
    def __init__(self, x, f):
        args = (x,)
        super().__init__(args, f)

    def __call__(self, expr):
        return Or(*[self.apply((arg,)) for arg in expr.args])


class ExactlyOneArg(ArgFactHandler):
    """
    Apply the function to exactly one argument of the expression.

    Examples
    ========

    >>> from sympy import Q
    >>> from sympy.assumptions.sathandlers import ExactlyOneArg
    >>> from sympy.abc import x, y
    >>> onearg = ExactlyOneArg(x, Q.positive(x))
    >>> onearg(x*y)
    (Q.positive(x) & ~Q.positive(y)) | (Q.positive(y) & ~Q.positive(x))

    """
    def __init__(self, x, f):
        args = (x,)
        super().__init__(args, f)

    def __call__(self, expr):
        pred_args = [self.apply((arg,)) for arg in expr.args]
        res = Or(*[And(pred_args[i], *[~lit for lit in pred_args[:i] +
            pred_args[i+1:]]) for i in range(len(pred_args))])
        return res


### Fact registry ###

class FactRegistry(MutableMapping):
    """
    Base class for fact registries

    """
    def __init__(self, d=None):
        d = d or {}
        self.d = defaultdict(frozenset, d)
        super().__init__()

    def __setitem__(self, key, item):
        self.d[key] = frozenset(item)

    def __delitem__(self, key):
        del self.d[key]

    def __iter__(self):
        return self.d.__iter__()

    def __len__(self):
        return len(self.d)

    def __repr__(self):
        return repr(self.d)

    def __call__(self, expr):
        handlers = self[expr.func]
        return {h(expr) for h in handlers}

    def register(self, item):
        def _(func):
            self[item] |= {func}
            return func
        return _

    def register_many(self, *classes):
        def _(func):
            for cls in classes:
                self[cls] |= {func}
            return func
        return _


class ClassFactRegistry(FactRegistry):
    """
    Register handlers against classes.

    Explanation
    ===========

    ``register`` method registers the handler function for the class.

    ``registry[C]`` returns a set of handlers for class ``C``, or any
    of its superclasses.

    ``registry(expr)`` returns a set of facts for *expr*.

    Examples
    ========

    Here, we register the facts for ``Abs``.

    >>> from sympy import Abs, Q
    >>> from sympy.logic.boolalg import Equivalent
    >>> from sympy.assumptions.sathandlers import ClassFactRegistry
    >>> reg = ClassFactRegistry()
    >>> @reg.register(Abs)
    ... def f1(expr):
    ...     return Q.nonnegative(expr)
    >>> @reg.register(Abs)
    ... def f2(expr):
    ...     arg = expr.args[0]
    ...     return Equivalent(~Q.zero(arg), ~Q.zero(expr))

    Getting the item of the registry with class returns the registered
    handlers for the class.

    >>> reg[Abs] == frozenset({f1, f2})
    True

    Calling the registry with expression returns the defined facts for the
    expression.

    >>> from sympy.abc import x
    >>> reg(Abs(x))
    {Q.nonnegative(Abs(x)), Equivalent(~Q.zero(x), ~Q.zero(Abs(x)))}

    Multiple facts can be registered at once by letting handler return an
    iterable.

    >>> reg2 = ClassFactRegistry()
    >>> @reg2.register(Abs)
    ... def _(expr):
    ...     arg = expr.args[0]
    ...     return [Q.even(arg) >> Q.even(expr), Q.odd(arg) >> Q.odd(expr)]
    >>> reg2(Abs(x))
    {Implies(Q.even(x), Q.even(Abs(x))), Implies(Q.odd(x), Q.odd(Abs(x)))}

    """
    def __getitem__(self, key):
        ret = self.d[key]
        for k in self.d:
            if issubclass(key, k):
                ret |= self.d[k]
        return ret

    def __call__(self, expr):
        handlers = self[expr.func]
        ret = set()
        for h in handlers:
            val = h(expr)
            if isinstance(val, (bool, Boolean)):
                ret.add(val)
            elif iterable(val):
                ret.update(val)
            else:
                raise TypeError("Unexpected result from %s(%s) : %s" % (self, expr, val))
        return ret

class_fact_registry = ClassFactRegistry()



### Class fact registration ###

x = Symbol('x')

## Abs ##

@class_fact_registry.register(Abs)
def _(expr):
    arg = expr.args[0]
    return [Q.nonnegative(expr), 
            Equivalent(~Q.zero(arg), ~Q.zero(expr)),
            Q.even(arg) >> Q.even(expr),
            Q.odd(arg) >> Q.odd(expr),
            Q.integer(arg) >> Q.integer(expr),
            ]


### Add ##

@class_fact_registry.register(Add)
def _(expr):
    return [AllArgs(x, Q.positive(x))(expr) >> Q.positive(expr),
            AllArgs(x, Q.negative(x))(expr) >> Q.negative(expr),
            AllArgs(x, Q.real(x))(expr) >> Q.real(expr),
            AllArgs(x, Q.rational(x))(expr) >> Q.rational(expr),
            AllArgs(x, Q.integer(x))(expr) >> Q.integer(expr),
            ExactlyOneArg(x, ~Q.integer(x))(expr) >> ~Q.integer(expr),
            ]

@class_fact_registry.register(Add)
def _(expr):
    allargs_real = AllArgs(x, Q.real(x))(expr)
    onearg_irrational = ExactlyOneArg(x, Q.irrational(x))(expr)
    return Implies(allargs_real, Implies(onearg_irrational, Q.irrational(expr)))


### Mul ###

@class_fact_registry.register(Mul)
def _(expr):
    return [Equivalent(Q.zero(expr), AnyArgs(x, Q.zero(x))(expr)),
            AllArgs(x, Q.positive(x))(expr) >> Q.positive(expr),
            AllArgs(x, Q.real(x))(expr) >> Q.real(expr),
            AllArgs(x, Q.rational(x))(expr) >> Q.rational(expr),
            AllArgs(x, Q.integer(x))(expr) >> Q.integer(expr),
            ExactlyOneArg(x, ~Q.rational(x))(expr) >> ~Q.integer(expr),
            AllArgs(x, Q.commutative(x))(expr) >> Q.commutative(expr),
            ]

@class_fact_registry.register(Mul)
def _(expr):
    # Implicitly assumes Mul has more than one arg
    # Would be AllArgs(x, Q.prime(x) | Q.composite(x)) except 1 is composite
    # More advanced prime assumptions will require inequalities, as 1 provides
    # a corner case.
    allargs_prime = AllArgs(x, Q.prime(x))(expr)
    return Implies(allargs_prime, ~Q.prime(expr))

@class_fact_registry.register(Mul)
def _(expr):
    # General Case: Odd number of imaginary args implies mul is imaginary(To be implemented)
    allargs_imag_or_real = AllArgs(x, Q.imaginary(x) | Q.real(x))(expr)
    onearg_imaginary = ExactlyOneArg(x, Q.imaginary(x))(expr)
    return Implies(allargs_imag_or_real, Implies(onearg_imaginary, Q.imaginary(expr)))

@class_fact_registry.register(Mul)
def _(expr):
    allargs_real = AllArgs(x, Q.real(x))(expr)
    onearg_irrational = ExactlyOneArg(x, Q.irrational(x))(expr)
    return Implies(allargs_real, Implies(onearg_irrational, Q.irrational(expr)))

@class_fact_registry.register(Mul)
def _(expr):
    # Including the integer qualification means we don't need to add any facts
    # for odd, since the assumptions already know that every integer is
    # exactly one of even or odd.
    allargs_integer = AllArgs(x, Q.integer(x))(expr)
    anyargs_even = AnyArgs(x, Q.even(x))(expr)
    return Implies(allargs_integer, Equivalent(anyargs_even, Q.even(expr)))


### MatMul ###

@class_fact_registry.register(MatMul)
def _(expr):
    allargs_square = AllArgs(x, Q.square(x))(expr)
    allargs_invertible = AllArgs(x, Q.invertible(x))(expr)
    return Implies(allargs_square, Equivalent(Q.invertible(expr), allargs_invertible))


### Pow ###

@class_fact_registry.register(Pow)
def _(expr):
    base, exp = expr.base, expr.exp
    return [
        (Q.real(base) & Q.even(exp) & Q.nonnegative(exp)) >> Q.nonnegative(expr),
        (Q.nonnegative(base) & Q.odd(exp) & Q.nonnegative(exp)) >> Q.nonnegative(expr),
        (Q.nonpositive(base) & Q.odd(exp) & Q.nonnegative(exp)) >> Q.nonpositive(expr),
        Equivalent(Q.zero(expr), Q.zero(base) & Q.positive(exp))
    ]


### Integer ###

@class_fact_registry.register(Integer)
def _(expr):
    from sympy import isprime
    return [Equivalent(Q.prime(expr), isprime(expr)),
            Equivalent(Q.composite(expr), expr.is_composite),
    ]


### Numbers ###

_old_assump_getters = {
    Q.positive: lambda o: o.is_positive,
    Q.zero: lambda o: o.is_zero,
    Q.negative: lambda o: o.is_negative,
    Q.rational: lambda o: o.is_rational,
    Q.irrational: lambda o: o.is_irrational,
    Q.even: lambda o: o.is_even,
    Q.odd: lambda o: o.is_odd,
    Q.imaginary: lambda o: o.is_imaginary,
}

@class_fact_registry.register_many(Number, NumberSymbol, ImaginaryUnit)
def _(expr):
    ret = []
    for p, getter in _old_assump_getters.items():
        pred = p(expr)
        prop = getter(expr)
        if prop is not None:
            ret.append(Equivalent(pred, prop))
    return ret
