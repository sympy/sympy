from __future__ import print_function, division

from sympy.core import Basic

from sympy.assumptions.assume import global_assumptions, AppliedPredicate
from sympy.logic.inference import satisfiable
from sympy.logic.boolalg import And, Implies
from sympy.assumptions.ask import Q

def newask(proposition, assumptions=True, context=global_assumptions):
    # The relevant facts might introduce new keys, e.g., Q.zero(x*y) will
    # introduce the keys Q.zero(x) and Q.zero(y), so we need to run it until
    # we stop getting new things.  Hopefully this strategy won't lead to an
    # infinite loop in the future.
    relevant_facts = False
    old_relevant_facts = True
    while relevant_facts != old_relevant_facts:
        old_relevant_facts, relevant_facts = (relevant_facts,
            get_relevant_facts(proposition, assumptions & old_relevant_facts,
                context))

    # TODO: Can this be faster to do it in one pass using xor?
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

def get_relevant_facts(proposition, assumptions=True, context=global_assumptions):
    keys = proposition.atoms(AppliedPredicate)
    if isinstance(assumptions, Basic):
        # XXX: We need this since True/False are not Basic
        keys |= assumptions.atoms(AppliedPredicate)
    if context:
        keys |= And(*context).atoms(AppliedPredicate)

    relevant_facts = True

    # TODO: Write this in a more scalable and extendable way

    real_keys = [key for key in keys if key.func == Q.real]
    positive_keys = [key for key in keys if key.func == Q.positive]

    for key in real_keys:
        if Q.positive(key.args[0]) in positive_keys:
            relevant_facts &= Implies(Q.positive(key.args[0]), key)

    for key in positive_keys:
        if Q.real(key.args[0]) in real_keys:
            relevant_facts &= Implies(key, Q.real(key.args[0]))

    return relevant_facts
