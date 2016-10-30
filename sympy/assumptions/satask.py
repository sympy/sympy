from __future__ import print_function, division

from sympy.core import oo

from sympy.assumptions.assume import global_assumptions, AppliedPredicate
from sympy.logic.algorithms.dpll2 import KB, _satisfiable
from sympy.logic.boolalg import conjuncts, to_cnf
from sympy.assumptions.ask_generated import get_known_facts_cnf
from sympy.assumptions.sathandlers import fact_registry


class CNF(object):
    def __init__(self, clauses=None):
        if clauses is None:
            clauses = set()
        self.clauses = clauses

    def add(self, proposition):
        self.clauses |= conjuncts(to_cnf(proposition))

    def copy(self):
        return CNF(set(self.clauses))

    @classmethod
    def from_prop(cls, prop):
        res = CNF()
        res.add(prop)
        return res


def satask(proposition, assumptions=True, context=global_assumptions,
        use_known_facts=True, iterations=oo):
    ctx = CNF.from_prop(assumptions)
    for c in context:
        ctx.add(c)

    relevant_facts = get_all_relevant_facts(proposition, ctx,
        use_known_facts=use_known_facts, iterations=iterations)
    for c in relevant_facts:
        ctx.add(c)

    ctx2 = ctx.copy()
    ctx.add(proposition)
    can_be_true = _satisfiable(KB(ctx.clauses))
    ctx2.add(~proposition)
    can_be_false = _satisfiable(KB(ctx2.clauses))

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

def _get_exprs(prop):
    return {pred.args[0] for pred in prop.atoms(AppliedPredicate)}

def _extract_exprs(proposition, ctx):
    exprs = _get_exprs(proposition)
    for c in ctx.clauses:
        exprs |= _get_exprs(c)
    return exprs


def get_relevant_facts(exprs):
    newexprs = set()
    newfacts = set()
    for expr in exprs:
        for fact in fact_registry[expr.func]:
            newfact = fact.rcall(expr)
            newfacts.add(newfact)
            newexprs |= _get_exprs(newfact)

    return newfacts, newexprs


def get_all_relevant_facts(
        proposition, ctx, use_known_facts=True, iterations=oo):
    # The relevant facts might introduce new keys, e.g., Q.zero(x*y) will
    # introduce the keys Q.zero(x) and Q.zero(y), so we need to run it until
    # we stop getting new things. Hopefully this strategy won't lead to an
    # infinite loop in the future.
    i = 0
    relevant_facts = set()
    exprs = _extract_exprs(proposition, ctx)
    all_exprs = set()
    while exprs and i < iterations:
        all_exprs |= exprs
        (newfacts, newexprs) = get_relevant_facts(exprs)
        exprs = newexprs - all_exprs
        relevant_facts |= newfacts
        i += 1

    if use_known_facts:
        known_facts_CNF = CNF.from_prop(get_known_facts_cnf())
        for expr in all_exprs:
            for p in known_facts_CNF.clauses:
                relevant_facts.add(p.rcall(expr))

    return relevant_facts
