from __future__ import print_function, division

from sympy.core import oo

from sympy.assumptions.assume import global_assumptions, AppliedPredicate
from sympy.logic.algorithms.dpll2 import SATEncoding, EncodedCNF, _satisfiable
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

    def extend(self, props):
        for p in props:
            self.add(p)

    def copy(self):
        return CNF(set(self.clauses))

    @classmethod
    def from_prop(cls, prop):
        res = CNF()
        res.add(prop)
        return res

    def __iand__(self, other):
        self.clauses |= other.clauses
        return self

    def add_relevant_facts(self, proposition, iterations=oo):
        # The relevant facts might introduce new keys, e.g., Q.zero(x*y) will
        # introduce the keys Q.zero(x) and Q.zero(y), so we need to run it
        # until we stop getting new things. Hopefully this strategy won't lead
        # to an infinite loop in the future.
        i = 0
        exprs = _extract_exprs(proposition, self)
        all_exprs = set()
        while exprs and i < iterations:
            all_exprs |= exprs
            (newfacts, newexprs) = get_relevant_facts(exprs)
            exprs = newexprs - all_exprs
            self.extend(newfacts)
            i += 1
        return all_exprs

    def rcall(self, expr):
        clauses = set(p.rcall(expr) for p in self.clauses)
        return CNF(clauses)

    def satisfiable(self):
        return _satisfiable(EncodedCNF.from_cnf(self))


def satask(proposition, assumptions=True, context=global_assumptions,
        use_known_facts=True, iterations=oo):
    ctx = CNF.from_prop(assumptions)
    for c in context:
        ctx.add(c)

    exprs = ctx.add_relevant_facts(proposition, iterations=iterations)

    if use_known_facts:
        known_facts_CNF = CNF.from_prop(get_known_facts_cnf())
        kf_encoded = EncodedCNF.from_cnf(known_facts_CNF)

        def translate_literal(lit, delta):
            if lit > 0:
                return lit + delta
            else:
                return lit - delta

        def translate_data(data, delta):
            return [{translate_literal(i, delta) for i in clause} for clause in data]
        data = []
        symbols = []
        n_lit = len(kf_encoded.encoding.symbols)
        for i, expr in enumerate(exprs):
            symbols += [pred(expr) for pred in kf_encoded.encoding.symbols]
            data += translate_data(kf_encoded.data, i * n_lit)

        encoding = SATEncoding(symbols)
        sat_true = EncodedCNF(data, encoding)
    else:
        sat_true = EncodedCNF.from_cnf(ctx)

    for clause in ctx.clauses:
        sat_true.data.append(sat_true.encoding.encode(clause))

    sat_false = sat_true.copy()

    sat_true.add_prop(proposition)
    can_be_true = _satisfiable(sat_true)

    sat_false.add_prop(~proposition)
    can_be_false = _satisfiable(sat_false)

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
