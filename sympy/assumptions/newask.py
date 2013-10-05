from __future__ import print_function, division

from sympy.core import Basic, Mul

from sympy.assumptions.assume import global_assumptions, AppliedPredicate
from sympy.logic.inference import satisfiable
from sympy.logic.boolalg import And, Implies, Equivalent, Or
from sympy.assumptions.ask import Q

def newask(proposition, assumptions=True, context=global_assumptions):
    relevant_facts = get_all_relevant_facts(proposition, assumptions, context)

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
    zero_keys = [key for key in keys if key.func == Q.zero]
    nonzero_keys = [key for key in keys if key.func == Q.nonzero]

    # To keep things straight, for implications, only worry about the
    # Implies(key, Q.something(key.args[0])) fact.

    for key in positive_keys:
        relevant_facts &= Implies(key, Q.real(key.args[0]))

    for key in zero_keys:
        relevant_facts &= Equivalent(key, ~Q.nonzero(key.args[0]))
        relevant_facts &= Implies(key, ~Q.positive(key.args[0]))
        relevant_facts &= Implies(key, Q.real(key.args[0]))

        # Now for something interesting...
        if isinstance(key.args[0], Mul):
            relevant_facts &= Equivalent(key, Or(*[Q.zero(i) for i in
                key.args[0].args]))

    for key in nonzero_keys:
        relevant_facts &= Equivalent(key, ~Q.zero(key.args[0]))

    return relevant_facts

def get_all_relevant_facts(proposition, assumptions=True, context=global_assumptions):
    # The relevant facts might introduce new keys, e.g., Q.zero(x*y) will
    # introduce the keys Q.zero(x) and Q.zero(y), so we need to run it until
    # we stop getting new things.  Hopefully this strategy won't lead to an
    # infinite loop in the future.
    relevant_facts = True
    old_relevant_facts = False
    while relevant_facts != old_relevant_facts:
        old_relevant_facts, relevant_facts = (relevant_facts,
            get_relevant_facts(proposition, assumptions & relevant_facts,
                context))

    return relevant_facts
