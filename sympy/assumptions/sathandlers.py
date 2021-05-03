from collections import defaultdict
from collections.abc import MutableMapping

from sympy.assumptions.ask import Q
from sympy.core import (Add, Mul, Pow, Number, NumberSymbol, Symbol)
from sympy.core.numbers import ImaginaryUnit
from sympy.functions.elementary.complexes import Abs
from sympy.logic.boolalg import (Equivalent, And, Or, Implies)
from sympy.matrices.expressions import MatMul

# APIs here may be subject to change


### Helper functions ###

def allargs(symbol, fact, expr):
    """
    Apply all arguments of the expression to the fact structure.

    Parameters
    ==========

    symbol : Symbol
        A placeholder symbol.

    fact : Boolean
        Resulting ``Boolean`` expression.

    expr : Expr

    Examples
    ========

    >>> from sympy import Q
    >>> from sympy.assumptions.sathandlers import allargs
    >>> from sympy.abc import x, y
    >>> allargs(x, Q.negative(x) | Q.positive(x), x*y)
    (Q.negative(x) | Q.positive(x)) & (Q.negative(y) | Q.positive(y))

    """
    return And(*[fact.subs(symbol, arg) for arg in expr.args])


def anyarg(symbol, fact, expr):
    """
    Apply any argument of the expression to the fact structure.

    Parameters
    ==========

    symbol : Symbol
        A placeholder symbol.

    fact : Boolean
        Resulting ``Boolean`` expression.

    expr : Expr

    Examples
    ========

    >>> from sympy import Q
    >>> from sympy.assumptions.sathandlers import anyarg
    >>> from sympy.abc import x, y
    >>> anyarg(x, Q.negative(x) & Q.positive(x), x*y)
    (Q.negative(x) & Q.positive(x)) | (Q.negative(y) & Q.positive(y))

    """
    return Or(*[fact.subs(symbol, arg) for arg in expr.args])


def exactlyonearg(symbol, fact, expr):
    """
    Apply exactly one argument of the expression to the fact structure.

    Parameters
    ==========

    symbol : Symbol
        A placeholder symbol.

    fact : Boolean
        Resulting ``Boolean`` expression.

    expr : Expr

    Examples
    ========

    >>> from sympy import Q
    >>> from sympy.assumptions.sathandlers import exactlyonearg
    >>> from sympy.abc import x, y
    >>> exactlyonearg(x, Q.positive(x), x*y)
    (Q.positive(x) & ~Q.positive(y)) | (Q.positive(y) & ~Q.positive(x))

    """
    pred_args = [fact.subs(symbol, arg) for arg in expr.args]
    res = Or(*[And(pred_args[i], *[~lit for lit in pred_args[:i] +
        pred_args[i+1:]]) for i in range(len(pred_args))])
    return res


### Fact registry ###

class ClassFactRegistry(MutableMapping):
    """
    Register handlers against classes.

    Explanation
    ===========

    ``register`` method registers the handler function for a class. Here,
    handler function should return a single fact.

    ``registry(expr)`` returns a dictionary of computed facts and their handler
    functions for *expr*.

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

    Calling the registry with expression returns a dictionary of computed
    facts and their handler functions for the expression.

    >>> from sympy.abc import x
    >>> reg(Abs(x)) # doctest: +SKIP
    {Q.nonnegative(Abs(x)): <function __main__.f1(expr)>,
     Equivalent(~Q.zero(x), ~Q.zero(Abs(x))): <function __main__.f2(expr)>}

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

    def register(self, cls):
        def _(func):
            self.d[cls] |= {func}
            return func
        return _

    def __getitem__(self, key):
        ret = self.d[key]
        for k in self.d:
            if issubclass(key, k):
                ret |= self.d[k]
        return ret

    def __call__(self, expr):
        handlers = self[expr.func]
        return {h(expr): h for h in handlers}

class_fact_registry = ClassFactRegistry()



### Class fact registration ###

x = Symbol('x')

## Abs ##

for fact in [
    lambda expr: Q.nonnegative(expr),
    lambda expr: Equivalent(~Q.zero(expr.args[0]), ~Q.zero(expr)),
    lambda expr: Q.even(expr.args[0]) >> Q.even(expr),
    lambda expr: Q.odd(expr.args[0]) >> Q.odd(expr),
    lambda expr: Q.integer(expr.args[0]) >> Q.integer(expr),
]:
    class_fact_registry.register(Abs)(fact)


### Add ##

for fact in [
    lambda expr: allargs(x, Q.positive(x), expr) >> Q.positive(expr),
    lambda expr: allargs(x, Q.negative(x), expr) >> Q.negative(expr),
    lambda expr: allargs(x, Q.real(x), expr) >> Q.real(expr),
    lambda expr: allargs(x, Q.rational(x), expr) >> Q.rational(expr),
    lambda expr: allargs(x, Q.integer(x), expr) >> Q.integer(expr),
    lambda expr: exactlyonearg(x, ~Q.integer(x), expr) >> ~Q.integer(expr)
]:
    class_fact_registry.register(Add)(fact)


@class_fact_registry.register(Add)
def _(expr):
    allargs_real = allargs(x, Q.real(x), expr)
    onearg_irrational = exactlyonearg(x, Q.irrational(x), expr)
    return Implies(allargs_real, Implies(onearg_irrational, Q.irrational(expr)))


### Mul ###

for fact in [
    lambda expr: Equivalent(Q.zero(expr), anyarg(x, Q.zero(x), expr)),
    lambda expr: allargs(x, Q.positive(x), expr) >> Q.positive(expr),
    lambda expr: allargs(x, Q.real(x), expr) >> Q.real(expr),
    lambda expr: allargs(x, Q.rational(x), expr) >> Q.rational(expr),
    lambda expr: allargs(x, Q.integer(x), expr) >> Q.integer(expr),
    lambda expr: exactlyonearg(x, ~Q.rational(x), expr) >> ~Q.integer(expr),
    lambda expr: allargs(x, Q.commutative(x), expr) >> Q.commutative(expr),
]:
    class_fact_registry.register(Mul)(fact)

@class_fact_registry.register(Mul)
def _(expr):
    # Implicitly assumes Mul has more than one arg
    # Would be allargs(x, Q.prime(x) | Q.composite(x)) except 1 is composite
    # More advanced prime assumptions will require inequalities, as 1 provides
    # a corner case.
    allargs_prime = allargs(x, Q.prime(x), expr)
    return Implies(allargs_prime, ~Q.prime(expr))

@class_fact_registry.register(Mul)
def _(expr):
    # General Case: Odd number of imaginary args implies mul is imaginary(To be implemented)
    allargs_imag_or_real = allargs(x, Q.imaginary(x) | Q.real(x), expr)
    onearg_imaginary = exactlyonearg(x, Q.imaginary(x), expr)
    return Implies(allargs_imag_or_real, Implies(onearg_imaginary, Q.imaginary(expr)))

@class_fact_registry.register(Mul)
def _(expr):
    allargs_real = allargs(x, Q.real(x), expr)
    onearg_irrational = exactlyonearg(x, Q.irrational(x), expr)
    return Implies(allargs_real, Implies(onearg_irrational, Q.irrational(expr)))

@class_fact_registry.register(Mul)
def _(expr):
    # Including the integer qualification means we don't need to add any facts
    # for odd, since the assumptions already know that every integer is
    # exactly one of even or odd.
    allargs_integer = allargs(x, Q.integer(x), expr)
    anyarg_even = anyarg(x, Q.even(x), expr)
    return Implies(allargs_integer, Equivalent(anyarg_even, Q.even(expr)))


### MatMul ###

@class_fact_registry.register(MatMul)
def _(expr):
    allargs_square = allargs(x, Q.square(x), expr)
    allargs_invertible = allargs(x, Q.invertible(x), expr)
    return Implies(allargs_square, Equivalent(Q.invertible(expr), allargs_invertible))


### Pow ###

for fact in [
    lambda expr: (Q.real(expr.base) & Q.even(expr.exp) & Q.nonnegative(expr.exp)) >> Q.nonnegative(expr),
    lambda expr: (Q.nonnegative(expr.base) & Q.odd(expr.exp) & Q.nonnegative(expr.exp)) >> Q.nonnegative(expr),
    lambda expr: (Q.nonpositive(expr.base) & Q.odd(expr.exp) & Q.nonnegative(expr.exp)) >> Q.nonpositive(expr),
    lambda expr: Equivalent(Q.zero(expr), Q.zero(expr.base) & Q.positive(expr.exp)),
]:
    class_fact_registry.register(Pow)(fact)


### Numbers ###

for cls in (Number, NumberSymbol, ImaginaryUnit):
    for fact in [
        lambda o: Equivalent(Q.positive(o), o.is_positive) if o.is_positive is not None else True,
        lambda o: Equivalent(Q.zero(o), o.is_zero) if o.is_zero is not None else True,
        lambda o: Equivalent(Q.negative(o), o.is_negative) if o.is_negative is not None else True,
        lambda o: Equivalent(Q.rational(o), o.is_rational) if o.is_rational is not None else True,
        lambda o: Equivalent(Q.irrational(o), o.is_irrational) if o.is_irrational is not None else True,
        lambda o: Equivalent(Q.even(o), o.is_even) if o.is_even is not None else True,
        lambda o: Equivalent(Q.odd(o), o.is_odd) if o.is_odd is not None else True,
        lambda o: Equivalent(Q.imaginary(o), o.is_imaginary) if o.is_imaginary is not None else True,
        lambda o: Equivalent(Q.prime(o), o.is_prime) if o.is_prime is not None else True,
        lambda o: Equivalent(Q.composite(o), o.is_composite) if o.is_composite is not None else True,
    ]:
        class_fact_registry.register(cls)(fact)
