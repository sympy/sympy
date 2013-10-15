from __future__ import print_function, division

from collections import MutableMapping, defaultdict

from sympy.core import Add, Mul, Pow
from sympy.core.sympify import _sympify

from sympy.matrices.expressions import MatMul

from sympy.assumptions.ask import Q
from sympy.assumptions.assume import Predicate, AppliedPredicate
from sympy.logic.boolalg import (Equivalent, Implies, And, Or, BooleanFunction,
    _find_predicates)

# APIs here may be subject to change

class Handler(object):
    def get_relevant_fact(self, key):
        return True

class ArgHandler(Handler):
    """
    A handler relationships among args

    This handles cases when there is some relationship between
    Q.assumption(expr) and Q.assumption(arg) for arg in expr.args.

    Subclasses should override the get_relationship method.

    Alternately, you can instantiate this class with the get_relationship
    function as the second argument.  That is, ``handler = ArgHandler(predicate, lambda
    key, keyed_args: relationship)`` is the same as::

        class MyHandler(ArgHandler):
            def get_relationship(self, key, keyed_args):
                return relationship

        handler = MyHandler(predicate)

    """

    def __init__(self, predicate, get_relationship=None):
        self.predicate = predicate
        if get_relationship:
            self.get_relationship = get_relationship

    def get_relevant_fact(self, key):
        expr = key.args[0]
        if key.func == self.predicate: # TODO: isinstance doesn't work here
            return self.get_relationship(key, list(map(self.predicate, expr.args)))
        return True

class EquivalentAnyArgs(ArgHandler):
    """
    Q.assumption(expr) iff any(Q.assumption(arg) for arg in expr.args)
    """
    def get_relationship(self, key, keyed_args):
        return Equivalent(key, Or(*keyed_args))

class EquivalentAllArgs(ArgHandler):
    """
    Q.assumption(expr) iff all(Q.assumption(arg) for arg in expr.args)
    """
    def get_relationship(self, key, keyed_args):
        return Equivalent(key, And(*keyed_args))

class AllArgsImplies(ArgHandler):
    """
    all(Q.assumption(arg) for arg in expr.args) implies Q.assumption(expr)
    (but the reverse implication does not hold)
    """
    def get_relationship(self, key, keyed_args):
        return Implies(And(*keyed_args), key)

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
        return obj.apply()

    def apply(self):
        return self

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

class ClassHandlerRegistry(MutableMapping):
    """
    Register handlers against classes

    ``registry[C] = handler`` registers ``handler`` for class
    ``C``. ``registry[C]`` returns a set of handlers for class ``C``, or any
    of its superclasses.
    """
    def __init__(self):
        self.d = defaultdict(frozenset)
        super(ClassHandlerRegistry, self).__init__()

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

handler_registry = ClassHandlerRegistry()

def register_handler(klass, handler, registry=handler_registry):
    registry[klass] |= set([handler])

for handler, klass in [
    (EquivalentAnyArgs(Q.zero), Mul),
    (EquivalentAllArgs(Q.invertible), MatMul),
    (AllArgsImplies(Q.positive), Add),
    (ArgHandler(Q.zero, lambda key, keyed_args: Implies(key, keyed_args[0])),
    Pow),
    (ArgHandler(Q.zero, lambda key, keyed_args:
        Implies(And(Q.zero(keyed_args[0].args[0]), Q.positive(keyed_args[1].args[0])), key)), Pow),
    ]:

    register_handler(klass, handler)
