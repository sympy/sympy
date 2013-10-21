from __future__ import print_function, division

from collections import MutableMapping, defaultdict

from sympy.core import (Add, Mul, Pow, Integer, Lambda, Dummy, Number,
    NumberSymbol,)
from sympy.core.numbers import ImaginaryUnit
from sympy.core.sympify import _sympify
from sympy.core.rules import Transform
from sympy.core.logic import fuzzy_or, fuzzy_and
from sympy.matrices.expressions import MatMul

from sympy.assumptions.ask import Q
from sympy.assumptions.assume import Predicate, AppliedPredicate
from sympy.logic.boolalg import (Equivalent, Implies, And, Or, BooleanFunction,
    _find_predicates)

# APIs here may be subject to change

# XXX: Better name?
class UnevaluatedOnFree(BooleanFunction):
    """
    Represents a Boolean function that remains unevaluated on free predicates

    This is intended to be a superclass of other classes, which define the
    behavior on singly applied predicates.

    A free predicate is a predicate that is not applied, or a combination
    thereof. For example, Q.zero or Or(Q.positive, Q.negative).

    A singly applied predicate is a free predicated applied everywhere to a
    single expression. For instance, Q.zero(x) and Or(Q.positive(x*y),
    Q.negative(x*y)) are singly applied, but Or(Q.positive(x), Q.negative(y))
    and Or(Q.positive, Q.negative(y)) are not.

    The boolean literals True and False are considered to be both free and
    singly applied.

    This class raises ValueError unless the input is a free predicate or a
    singly applied predicate.

    On a free predicate, this class remains unevaluated. On a singly applied
    predicate, the method apply() is called and returned, or the original
    expression returned if apply() returns None. When apply() is called,
    self.expr is set to the unique expression that the predicates are applied
    at.

    The typical usage is to create this class with free predicates and
    evaluate it using .rcall().

    """
    def __new__(cls, arg):
        # Mostly type checking here
        arg = _sympify(arg)
        predicates = arg.atoms(Predicate)
        applied_predicates = arg.atoms(AppliedPredicate)
        if predicates and applied_predicates:
            raise ValueError("arg must be either completely free or singly applied")
        if not applied_predicates:
            obj = BooleanFunction.__new__(cls, arg)
            obj.expr = None
            return obj
        predicate_args = set([pred.args[0] for pred in applied_predicates])
        if len(predicate_args) > 1:
            raise ValueError("The AppliedPredicates in arg must be applied to a single expression.")
        obj = BooleanFunction.__new__(cls, arg)
        obj.expr = predicate_args.pop()
        applied = obj.apply()
        if applied is None:
            return obj
        return applied

    def apply(self):
        return

class AllArgs(UnevaluatedOnFree):
    """
    Class representing vectorizing a predicate over all the .args of an
    expression

    See the docstring of UnevaluatedOnFree for more information on this
    class.

    The typical usage is to evaluate predicates with expressions using .rcall().

    Example
    =======

    >>> from sympy.assumptions.newhandlers import AllArgs
    >>> from sympy import symbols, Q
    >>> x, y = symbols('x y')
    >>> a = AllArgs(Q.positive | Q.negative)
    >>> a
    AllArgs(Or(Q.negative, Q.positive))
    >>> a.rcall(x*y)
    And(Or(Q.negative(x), Q.positive(x)), Or(Q.negative(y), Q.positive(y)))
    """

    def apply(self):
        return And(*[self.args[0].xreplace({self.expr: arg}) for arg in
            self.expr.args])


class AnyArgs(UnevaluatedOnFree):
    """
    Class representing vectorizing a predicate over any of the .args of an
    expression

    See the docstring of UnevaluatedOnFree for more information on this
    class.

    The typical usage is to evaluate predicates with expressions using .rcall().

    Example
    =======

    >>> from sympy.assumptions.newhandlers import AnyArgs
    >>> from sympy import symbols, Q
    >>> x, y = symbols('x y')
    >>> a = AnyArgs(Q.positive & Q.negative)
    >>> a
    AnyArgs(And(Q.negative, Q.positive))
    >>> a.rcall(x*y)
    Or(And(Q.negative(x), Q.positive(x)), And(Q.negative(y), Q.positive(y)))
    """

    def apply(self):
        return Or(*[self.args[0].xreplace({self.expr: arg}) for arg in
            self.expr.args])

# TODO: How can we avoid calling evalf on the same expression multiple times
# (once for each evalf fact)?

def _old_assump_replacer(obj):
    # Things to be careful of:
    # - real means real or infinite in the old assumptions.
    # - nonzero does not imply real in the old assumptions.
    # - finite means finite and not zero in the old assumptions.
    if not isinstance(obj, AppliedPredicate):
        return obj

    e = obj.args[0]
    ret = None

    if obj.func == Q.positive:
        ret = fuzzy_and([e.is_finite, e.is_positive])
    if obj.func == Q.zero:
        ret = e.is_zero
    if obj.func == Q.negative:
        ret = fuzzy_and([e.is_finite, e.is_negative])
    if obj.func == Q.nonpositive:
        ret = fuzzy_and([fuzzy_or([e.is_zero, e.is_finite]), e.is_nonpositive])
    if obj.func == Q.nonzero:
        ret = fuzzy_and([e.is_real, e.is_finite, e.is_nonzero])
    if obj.func == Q.nonnegative:
        ret = fuzzy_and([fuzzy_or([e.is_zero, e.is_finite]), e.is_nonnegative])
    if ret is None:
        return obj
    return ret

def evaluate_old_assump(pred):
    """
    Replace assumptions of expressions replaced with their values in the old
    assumptions (like Q.negative(-1) => True). Useful because some direct
    computations for numeric objects is defined most conveniently in the old
    assumptions.

    """
    return pred.xreplace(Transform(_old_assump_replacer))

class CheckOldAssump(UnevaluatedOnFree):
    def apply(self):
        return Equivalent(self.args[0], evaluate_old_assump(self.args[0]))

class CheckIsPrime(UnevaluatedOnFree):
    def apply(self):
        from sympy import isprime
        return Equivalent(self.args[0], isprime(self.expr))

class CustomLambda(object):
    """
    Interface to lambda with rcall

    Workaround until we get a better way to represent certain facts.
    """
    def __init__(self, lamda):
        self.lamda = lamda
    def rcall(self, *args):
        return self.lamda(*args)

class ClassFactRegistry(MutableMapping):
    """
    Register handlers against classes

    ``registry[C] = handler`` registers ``handler`` for class
    ``C``. ``registry[C]`` returns a set of handlers for class ``C``, or any
    of its superclasses.
    """
    def __init__(self, d=None):
        d = d or {}
        self.d = defaultdict(frozenset, d)
        super(ClassFactRegistry, self).__init__()

    def __setitem__(self, key, item):
        self.d[key] = frozenset(item)

    def __getitem__(self, key):
        ret = self.d[key]
        for k in self.d:
            if issubclass(key, k):
                ret |= self.d[k]
        return ret

    def __delitem__(self, key):
        del self.d[key]

    def __iter__(self):
        return self.d.__iter__()

    def __len__(self):
        return len(self.d)

    def __repr__(self):
        return repr(self.d)

fact_registry = ClassFactRegistry()

def register_fact(klass, fact, registry=fact_registry):
    registry[klass] |= set([fact])

for klass, fact in [
    (Mul, Equivalent(Q.zero, AnyArgs(Q.zero))),
    (MatMul, Implies(AllArgs(Q.square), Equivalent(Q.invertible, AllArgs(Q.invertible)))),
    (Add, Implies(AllArgs(Q.positive), Q.positive)),
    (Mul, Implies(AllArgs(Q.positive), Q.positive)),
    # This one can still be made easier to read. I think we need basic pattern
    # matching, so that we can just write Equivalent(Q.zero(x**y), Q.zero(x) & Q.positive(y))
    (Pow, CustomLambda(lambda power: Equivalent(Q.zero(power), Q.zero(power.base) & Q.positive(power.exp)))),
    (Integer, CheckIsPrime(Q.prime)),
    (Number, CheckOldAssump(Q.negative)),
    (Number, CheckOldAssump(Q.zero)),
    (Number, CheckOldAssump(Q.positive)),
    (Number, CheckOldAssump(Q.nonnegative)),
    (Number, CheckOldAssump(Q.nonzero)),
    (Number, CheckOldAssump(Q.nonpositive)),
    # For some reason NumberSymbol does not subclass Number
    (NumberSymbol, CheckOldAssump(Q.negative)),
    (NumberSymbol, CheckOldAssump(Q.zero)),
    (NumberSymbol, CheckOldAssump(Q.positive)),
    (NumberSymbol, CheckOldAssump(Q.nonnegative)),
    (NumberSymbol, CheckOldAssump(Q.nonzero)),
    (NumberSymbol, CheckOldAssump(Q.nonpositive)),
    (ImaginaryUnit, CheckOldAssump(Q.negative)),
    (ImaginaryUnit, CheckOldAssump(Q.zero)),
    (ImaginaryUnit, CheckOldAssump(Q.positive)),
    (ImaginaryUnit, CheckOldAssump(Q.nonnegative)),
    (ImaginaryUnit, CheckOldAssump(Q.nonzero)),
    (ImaginaryUnit, CheckOldAssump(Q.nonpositive)),
    ]:

    register_fact(klass, fact)
