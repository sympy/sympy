from __future__ import print_function, division

from sympy.core import Basic, Mul, Add, Pow
from sympy.matrices.expressions import MatMul

from sympy.assumptions.assume import global_assumptions, AppliedPredicate
from sympy.logic.inference import satisfiable
from sympy.logic.boolalg import And, Implies, Equivalent, Or
from sympy.assumptions.ask import Q
from sympy.utilities.iterables import sift
from sympy.assumptions.ask_generated import known_facts_cnf
from sympy.assumptions.newhandlers import handler_registry

def newask(proposition, assumptions=True, context=global_assumptions, use_known_facts=True):
    relevant_facts = get_all_relevant_facts(proposition, assumptions, context,
        use_known_facts=use_known_facts)

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

def get_relevant_facts(proposition, assumptions=True,
    context=global_assumptions, use_known_facts=True):
    keys = proposition.atoms(AppliedPredicate)
    if isinstance(assumptions, Basic):
        # XXX: We need this since True/False are not Basic
        keys |= assumptions.atoms(AppliedPredicate)
    if context:
        keys |= And(*context).atoms(AppliedPredicate)

    predicates = set([i.args[0] for i in keys])

    relevant_facts = True

    if use_known_facts:
        for predicate in predicates:
            relevant_facts &= known_facts_cnf.rcall(predicate)

    for key in keys:
        expr = key.args[0]
        for handler in handler_registry[expr.func]:
            relevant_facts &= handler.get_relevant_fact(key)

    return relevant_facts

def get_all_relevant_facts(proposition, assumptions=True,
    context=global_assumptions, use_known_facts=True):
    # The relevant facts might introduce new keys, e.g., Q.zero(x*y) will
    # introduce the keys Q.zero(x) and Q.zero(y), so we need to run it until
    # we stop getting new things.  Hopefully this strategy won't lead to an
    # infinite loop in the future.
    relevant_facts = True
    old_relevant_facts = False
    while relevant_facts != old_relevant_facts:
        old_relevant_facts, relevant_facts = (relevant_facts,
            get_relevant_facts(proposition, assumptions & relevant_facts,
                context, use_known_facts=use_known_facts))

    return relevant_facts
