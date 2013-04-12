#######################################
# Teach LogPy how to manipulate SymPy #
#######################################
from logpy.unify import (reify, unify_seq, unify_isinstance_list,
                         reify_isinstance_list, reify_dispatch, seq_registry)
from sympy import Basic, Symbol, Integer, Rational, Dummy

from sympy.assumptions import AppliedPredicate
slot_classes = Symbol, Integer, Rational
def seq_Basic(x):
    return (x.__class__,) + x.args

def seq_Predicate(x):
    return (x.__class__, x.func, x.arg)

def seq_slot(x):
    return (type(x),) + tuple(getattr(x, a) for a in x.__slots__)

def reify_Basic(u, s):
    return u.func(*[reify(arg, s) for arg in u.args])

def reify_slot(u, s):
    return u.func(*[reify(getattr(u, a), s) for a in u.__slots__])

def build(tup):
    if isinstance(tup, tuple) and issubclass(tup[0], Basic):
        return tup[0](*map(build, tup[1:]))
    else:
        return tup

seq_registry.extend([(slot_classes, seq_slot),
                     (AppliedPredicate, seq_Predicate),
                     (Basic, seq_Basic)])

reify_dispatch[Dummy] = lambda u, s: u  # Dangerous
reify_isinstance_list.append((slot_classes, reify_slot))
reify_isinstance_list.append((Basic, reify_Basic))

#######################
# Simplification code #
#######################

import logpy
from logpy.variables import variables
from logpy import var, eq
from logpy.assoccomm import eq_assoccomm as eqac
from logpy.goals import goalify
from sympy import assuming, ask

# TODO - make sure this doesn't happen prematurely
asko = goalify(ask)

def refine_one(expr, *assumptions, **kwargs):
    reduces = kwargs.get('reduces')
    vars = kwargs.get('vars')
    with assuming(*assumptions):
        with variables(*vars):
            source, target, condition = var(), var(), var()
            result = logpy.run(1, target, (reduces, source, target, condition),
                                          (eqac, source, expr),
                                          (asko, condition, True))
    return build(result[0]) if result else expr

###############
# Commutivity #
###############

from logpy.assoccomm import commutative, associative
from logpy import facts
from sympy import Add, Mul, MatAdd, MatMul
facts(commutative, [Add], [Mul], [MatAdd])
facts(associative, [Add], [Mul], [MatAdd], [MatMul])
