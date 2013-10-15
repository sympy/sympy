from sympy import Mul, Basic, Q, Expr, And, symbols, Equivalent, Implies, Or

from sympy.assumptions.newhandlers import (fact_registry, AllArgsImplies,
    EquivalentAnyArgs, EquivalentAllArgs, ClassFactRegistry, ArgHandler,
    AllArgs, UnevaluatedOnFree, AnyArgs)

from sympy.utilities.pytest import raises

x, y = symbols('x y')

def test_class_handler_registry():
    my_handler_registry = ClassFactRegistry()

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
        def get_relationship(self, key, keyed_args):
            return And(key, *keyed_args)

    handler = AndArgHandler(Q.zero)
    assert handler.get_relevant_fact(Q.zero(x*y)) == And(Q.zero(x*y),
        Q.zero(x), Q.zero(y))
    assert handler.get_relevant_fact(Q.positive(x*y)) == True

    handler = ArgHandler(Q.positive, lambda key, keyed_args: And(key,
        *keyed_args))
    assert handler.get_relevant_fact(Q.positive(x*y)) == And(Q.positive(x*y),
        Q.positive(x), Q.positive(y))
    assert handler.get_relevant_fact(Q.zero(x*y)) == True

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

def test_UnevaluatedOnFree():
    a = UnevaluatedOnFree(Q.positive)
    b = UnevaluatedOnFree(Q.positive | Q.negative)
    c = UnevaluatedOnFree(Q.positive & ~Q.positive) # It shouldn't do any deduction
    assert a.rcall(x) == UnevaluatedOnFree(Q.positive(x))
    assert b.rcall(x) == UnevaluatedOnFree(Q.positive(x) | Q.negative(x))
    assert c.rcall(x) == UnevaluatedOnFree(Q.positive(x) & ~Q.positive(x))
    assert a.rcall(x).expr == x
    raises(ValueError, lambda: UnevaluatedOnFree(Q.positive(x) | Q.negative))
    raises(ValueError, lambda: UnevaluatedOnFree(Q.positive(x) |
        Q.negative(y)))

    class MyUnevaluatedOnFree(UnevaluatedOnFree):
        def apply(self):
            return self.args[0]

    a = MyUnevaluatedOnFree(Q.positive)
    b = MyUnevaluatedOnFree(Q.positive | Q.negative)
    c = MyUnevaluatedOnFree(Q.positive(x))
    d = MyUnevaluatedOnFree(Q.positive(x) | Q.negative(x))

    assert a.rcall(x) == c == Q.positive(x)
    assert b.rcall(x) == d == Q.positive(x) | Q.negative(x)

    raises(ValueError, lambda: MyUnevaluatedOnFree(Q.positive(x) | Q.negative(y)))

def test_AllArgs():
    a = AllArgs(Q.zero)
    b = AllArgs(Q.positive | Q.negative)
    assert a.rcall(x*y) == And(Q.zero(x), Q.zero(y))
    assert b.rcall(x*y) == And(Q.positive(x) | Q.negative(x), Q.positive(y) | Q.negative(y))

def test_AnyArgs():
    a = AnyArgs(Q.zero)
    b = AnyArgs(Q.positive & Q.negative)
    assert a.rcall(x*y) == Or(Q.zero(x), Q.zero(y))
    assert b.rcall(x*y) == Or(Q.positive(x) & Q.negative(x), Q.positive(y) & Q.negative(y))
