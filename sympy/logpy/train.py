#######################################
# Teach LogPy how to manipulate SymPy #
#######################################

from sympy import Basic, Symbol, Number, Expr, Dummy
from sympy.assumptions import AppliedPredicate, Predicate

Basic._as_logpy = lambda self: (self.func, ) + tuple(self.args)

Predicate._as_logpy = lambda self: (type(self), self.name, self.handlers)

AppliedPredicate._as_logpy = lambda self: (type(self), self.func, self.arg)

Dummy._as_logpy = lambda self: (type(self), self.name, self.dummy_index)

slot_classes = Symbol, Number
for slot in slot_classes:
    slot._as_logpy = lambda self: (
            type(self),) + tuple(getattr(self, a) for a in self.__slots__)

def _from_logpy(tup):
    try:
        return tup[0](*tup[1:], evaluate=False)
    except TypeError:
        return tup[0](*tup[1:])

def _from_logpy_simple(tup):
    return tup[0](*tup[1:])

Basic._from_logpy = staticmethod(_from_logpy)

Predicate._from_logpy = staticmethod(_from_logpy_simple)

AppliedPredicate._from_logpy = staticmethod(lambda (t, pred, arg): pred(arg))

def dummy_from_logpy((t, name, idx)):
    obj = Dummy()
    obj.name = name
    obj.dummy_index = idx
    return obj

Dummy._from_logpy = staticmethod(dummy_from_logpy)

for slot in slot_classes:
    slot._from_logpy = staticmethod(_from_logpy_simple)

###############
# Commutivity #
###############

try:
    from logpy.assoccomm import commutative, associative
    from logpy import facts
    from sympy import Add, Mul, MatAdd, MatMul

    facts(commutative, [Add], [Mul], [MatAdd])
    facts(associative, [Add], [Mul], [MatAdd], [MatMul])
except ImportError:
    pass

