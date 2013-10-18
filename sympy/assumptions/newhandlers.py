from __future__ import print_function, division

from collections import MutableMapping, defaultdict

from sympy.core import Add, Mul, Pow, Integer, Lambda, Dummy
from sympy.core.sympify import _sympify

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
    predicate, the method apply is called and returned. In that case,
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

class CheckIsPrime(UnevaluatedOnFree):
    def apply(self):
        from sympy import isprime
        return Equivalent(Q.prime(self.expr), isprime(self.expr))

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
    def __init__(self):
        self.d = defaultdict(frozenset)
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

def register_handler(klass, handler, registry=fact_registry):
    registry[klass] |= set([handler])

for klass, fact in [
    (Mul, Equivalent(Q.zero, AnyArgs(Q.zero))),
    (MatMul, Implies(AllArgs(Q.square), Equivalent(Q.invertible, AllArgs(Q.invertible)))),
    (Add, Implies(AllArgs(Q.positive), Q.positive)),
    # This one can still be made easier to read. I think we need basic pattern
    # matching, so that we can just write Equivalent(Q.zero(x**y), Q.zero(x) & Q.positive(y))
    (Pow, CustomLambda(lambda power: Equivalent(Q.zero(power), Q.zero(power.base) & Q.positive(power.exp)))),
    (Integer, CheckIsPrime(Q.prime)),
    ]:

    register_handler(klass, fact)
