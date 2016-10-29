from __future__ import print_function, division

from sympy.core import oo, Tuple

from sympy.assumptions.assume import global_assumptions, AppliedPredicate
from sympy.logic.algorithms.dpll2 import dpll_satisfiable as satisfiable
from sympy.logic.boolalg import And
from sympy.assumptions.ask_generated import get_known_facts_cnf
from sympy.assumptions.sathandlers import fact_registry


def satask(proposition, assumptions=True, context=global_assumptions,
        use_known_facts=True, iterations=oo):
    relevant_facts = get_all_relevant_facts(proposition, assumptions, context,
        use_known_facts=use_known_facts, iterations=iterations)

    can_be_true = satisfiable(And(proposition, assumptions,
        relevant_facts, *context))
    can_be_false = satisfiable(And(~proposition, assumptions,
        relevant_facts, *context))

    if can_be_true and can_be_false:
        return None

    if can_be_true and not can_be_false:
        return True

    if not can_be_true and can_be_false:
        return False

    if not can_be_true and not can_be_false:
        # TODO: Run additional checks to see which combination of the
        # assumptions, global_assumptions, and relevant_facts are
        # inconsistent.
        raise ValueError("Inconsistent assumptions")


def _extract_exprs(proposition, assumptions, context):
    keys = proposition.atoms(AppliedPredicate)
    # XXX: We need this since True/False are not Basic
    keys |= Tuple(*assumptions).atoms(AppliedPredicate)
    if context:
        keys |= And(*context).atoms(AppliedPredicate)

    return {key.args[0] for key in keys}


def get_relevant_facts(relevant_facts, exprs, use_known_facts=True):
    newexprs = set()
    for expr in exprs:
        for fact in fact_registry[expr.func]:
            newfact = fact.rcall(expr)
            relevant_facts.add(newfact)
            newexprs |= set([key.args[0] for key in
                newfact.atoms(AppliedPredicate)])

    return relevant_facts, newexprs - exprs


def get_all_relevant_facts(proposition, assumptions=True,
        context=global_assumptions, use_known_facts=True, iterations=oo):
    # The relevant facts might introduce new keys, e.g., Q.zero(x*y) will
    # introduce the keys Q.zero(x) and Q.zero(y), so we need to run it until
    # we stop getting new things. Hopefully this strategy won't lead to an
    # infinite loop in the future.
    i = 0
    relevant_facts = set()
    exprs = _extract_exprs(proposition, And.make_args(assumptions), context)
    all_exprs = set()
    while exprs and i < iterations:
        all_exprs |= exprs
        (relevant_facts, exprs) = get_relevant_facts(
            relevant_facts, exprs, use_known_facts=use_known_facts)
        i += 1

    if use_known_facts:
        for expr in all_exprs:
            relevant_facts.add(get_known_facts_cnf().rcall(expr))

    return And(*relevant_facts)
