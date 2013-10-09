from __future__ import print_function, division

from collections import MutableMapping, defaultdict

from sympy.core import Add, Mul, Pow
from sympy.matrices.expressions import MatMul

from sympy.assumptions.ask import Q
from sympy.logic.boolalg import Equivalent, Implies, And, Or

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

# TODO: Create a handler registry system

class class_handler_registry(MutableMapping):
    """
    Register handlers against classes

    ``registry[C] = handler`` registers ``handler`` for class
    ``C``. ``registry[C]`` returns a set of handlers for class ``C``, or any
    of its superclasses.
    """
    def __init__(self):
        self.d = defaultdict(frozenset)
        super(class_handler_registry, self).__init__()

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

handler_registry = class_handler_registry()

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
