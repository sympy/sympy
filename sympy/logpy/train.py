#######################################
# Teach LogPy how to manipulate SymPy #
#######################################

from sympy import Basic, Symbol, Number, Expr, Dummy, Q
from sympy.assumptions import AppliedPredicate, Predicate
from term import termify

def new_eval_false(func, args):
    try:
        return func(*args, evaluate=False)
    except TypeError:
        return func(*args)

def new_simple(func, args):
    return func(*args)

Basic._term_op      = lambda self: self.func
Basic._term_args    = lambda self: self.args
Basic._term_new     = classmethod(new_eval_false)
Basic._term_isleaf  = lambda self: not self.args

Predicate._term_isleaf = lambda self: True
Predicate._term_new = new_simple

def dummy_new(cls, (name, idx)):
    obj = cls()
    obj.name = name
    obj.dummy_index = idx
    return obj

Dummy._term_args = lambda self: (self.name, self.dummy_index)
Dummy._term_new = classmethod(dummy_new)

slot_classes = Symbol, Number
for slot in slot_classes:
    termify(slot)
    slot._term_new = classmethod(lambda cls, args: cls(*args))

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

