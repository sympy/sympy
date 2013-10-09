from sympy import Mul, Basic, Q, Expr, And, symbols, Equivalent, Implies, Or

from sympy.assumptions.newhandlers import (handler_registry, AllArgsImplies,
    EquivalentAnyArgs, EquivalentAllArgs, class_handler_registry, ArgHandler)

x, y = symbols('x y')

def test_class_handler_registry():
    my_handler_registry = class_handler_registry()

    # The predicate doesn't matter here, so just use is_true
    all_args_implies_is_true = AllArgsImplies(Q.is_true)
    equiv_all_args_is_true = EquivalentAllArgs(Q.is_true)

    my_handler_registry[Mul] = set([all_args_implies_is_true])
    my_handler_registry[Expr] = set([equiv_all_args_is_true])

    assert my_handler_registry[Basic] == set()
    assert my_handler_registry[Expr] == set([equiv_all_args_is_true])
    assert my_handler_registry[Mul] == set([equiv_all_args_is_true, all_args_implies_is_true])

def test_ArgHandler():
    class AndArgHandler(ArgHandler):
        def get_relationship(self, key, mapped_args):
            return And(key, *mapped_args)

    handler = AndArgHandler(Q.zero)
    assert handler.get_relevant_fact(Q.zero(x*y)) == And(Q.zero(x*y),
    Q.zero(x), Q.zero(y))
    assert handler.get_relevant_fact(Q.positive(x*y)) == True

def test_EquivalentAllArgs():
    handler = EquivalentAllArgs(Q.invertible)
    assert handler.get_relevant_fact(Q.invertible(x*y)) == Equivalent(Q.invertible(x*y),
    And(Q.invertible(x), Q.invertible(y)))
    assert handler.get_relevant_fact(Q.zero(x*y)) == True

def test_EquivalentAnyArgs():
    handler = EquivalentAnyArgs(Q.zero)
    assert handler.get_relevant_fact(Q.zero(x*y)) == Equivalent(Q.zero(x*y),
    Or(Q.zero(x), Q.zero(y)))
    assert handler.get_relevant_fact(Q.integer(x*y)) == True

def test_AllArgsImplies():
    handler = AllArgsImplies(Q.positive)
    assert handler.get_relevant_fact(Q.positive(x + y)) == \
    Implies(And(Q.positive(x), Q.positive(y)), Q.positive(x + y))
    assert handler.get_relevant_fact(Q.zero(x + y)) == True
