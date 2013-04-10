#######################################
# Teach LogPy how to manipulate SymPy #
#######################################
from logpy.unify import (reify, unify_seq, unify_isinstance_list,
                         reify_isinstance_list, reify_dispatch)
from sympy import Basic, Symbol, Integer, Rational, Dummy

def unify_Basic(u, v, s):
    return unify_seq((type(u),) + u.args, (type(v),) + v.args, s)
def reify_Basic(u, s):
    return u.func(*[reify(arg, s) for arg in u.args])

from sympy.assumptions import AppliedPredicate
def unify_Predicate(u, v, s):
    return unify_seq((type(u), u.func, u.args[0]), (type(v), v.func, v.args[0]), s)

slot_classes = Symbol, Integer, Rational
def unify_slot(u, v, s):
    return unify_seq((type(u),) + tuple(getattr(u, a) for a in u.__slots__),
                     (type(v),) + tuple(getattr(v, a) for a in v.__slots__),
                     s)
def reify_slot(u, s):
    return u.func(*[reify(getattr(u, a), s) for a in u.__slots__])


reify_dispatch[Dummy] = lambda u, s: u  # Dangerous
unify_isinstance_list.append(((AppliedPredicate, AppliedPredicate), unify_Predicate))
unify_isinstance_list.append(((slot_classes, slot_classes), unify_slot))
unify_isinstance_list.append(((Basic, Basic), unify_Basic))
reify_isinstance_list.append((slot_classes, reify_slot))
reify_isinstance_list.append((Basic, reify_Basic))

#######################
# Simplification code #
#######################

import logpy
from logpy.variables import variables
from logpy import var, eq
from logpy.goals import goalify
from sympy import assuming, ask

asko = goalify(ask)

def refine_one(expr, *assumptions, **kwargs):
    reduces = kwargs.get('reduces')
    vars = kwargs.get('vars')
    with assuming(*assumptions):
        with variables(*vars):
            source, target, condition = var(), var(), var()
            result = logpy.run(1, target, (reduces, source, target, condition),
                                          (eq, source, expr),
                                          (asko, condition, True))
    return result[0] if result else expr
