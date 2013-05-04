#######################################
# Teach LogPy how to manipulate SymPy #
#######################################
from logpy.unify import reify, reify_isinstance_list, seq_registry
from sympy import Basic, Symbol, Number, Expr

from sympy.assumptions import AppliedPredicate, Predicate

Basic._as_tuple = lambda self: (self.func, ) + tuple(self.args)

Predicate._as_tuple = lambda self: (type(self), self.name, self.handlers)

AppliedPredicate._as_tuple = lambda self: (type(self), self.func, self.arg)

slot_classes = Symbol, Number
for slot in slot_classes:
    slot._as_tuple = lambda self: (
            type(self),) + tuple(getattr(self, a) for a in self.__slots__)

def from_tuple(tup):
    try:
        return tup[0](*tup[1:], evaluate=False)
    except TypeError:
        return tup[0](*tup[1:])

def from_tuple_simple(tup):
    return tup[0](*tup[1:])

Basic._from_tuple = staticmethod(from_tuple)

Predicate._from_tuple = staticmethod(from_tuple_simple)

AppliedPredicate._from_tuple = staticmethod(lambda (t, pred, arg): pred(arg))

for slot in slot_classes:
    slot._from_tuple = staticmethod(from_tuple_simple)

#######################
# Simplification code #
#######################

import logpy
from logpy.variables import variables
from logpy import var
from logpy.assoccomm import eq_assoccomm as eqac
from logpy.goals import goalify
from sympy import assuming, ask

# TODO - make sure this doesn't happen prematurely
asko = goalify(ask)

def refine_one(expr, *assumptions, **kwargs):
    reduces = kwargs['reduces']
    vars = kwargs['vars']
    with assuming(*assumptions):
        with variables(*vars):
            source, target, condition = var(), var(), var()
            result = logpy.run(1, target, (reduces, source, target, condition),
                                          (eqac, source, expr),
                                          (asko, condition, True))
    return result[0] if result else expr

###############
# Commutivity #
###############

from logpy.assoccomm import commutative, associative
from logpy import facts
from sympy import Add, Mul, MatAdd, MatMul

facts(commutative, [Add], [Mul], [MatAdd])
facts(associative, [Add], [Mul], [MatAdd], [MatMul])
