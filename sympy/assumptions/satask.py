from sympy import Symbol, S
from sympy.assumptions.ask_generated import get_all_known_facts
from sympy.assumptions.assume import global_assumptions, AppliedPredicate
from sympy.assumptions.sathandlers import fact_registry
from sympy.core import oo
from sympy.logic.inference import satisfiable
from sympy.assumptions.cnf import CNF, EncodedCNF


def satask(proposition, assumptions=True, context=global_assumptions,
        use_known_facts=True, iterations=oo):
    props = CNF.from_prop(proposition)
    _props = CNF.from_prop(~proposition)
    if context:
        tmp = CNF()
        context = tmp.extend(context)
    assumptions = CNF.from_prop(assumptions)

    sat = get_all_relevant_facts(props, assumptions, context,
        use_known_facts=use_known_facts, iterations=iterations)
    if context:
        sat.add_from_cnf(context)
    sat.add_from_cnf(assumptions)

    return check_satisfiability(props, _props, sat)


def check_satisfiability(prop, _prop, factbase):
    sat_true = factbase.copy()
    sat_false = factbase.copy()
    sat_true.add_from_cnf(prop)
    sat_false.add_from_cnf(_prop)
    can_be_true = satisfiable(sat_true)
    can_be_false = satisfiable(sat_false)

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


def get_relevant_facts(proposition, assumptions=None,
    context=None, exprs=None,
    relevant_facts=None):

    newexprs = set()

    if not assumptions:
        assumptions = CNF({S.true})

    if not relevant_facts:
        relevant_facts = set()

    def find_symbols(pred):
        if isinstance(pred, CNF):
            symbols = set()
            for a in pred.all_predicates():
                symbols |= find_symbols(a)
            return symbols
        if isinstance(pred.args, AppliedPredicate):
            return {pred.args[0]}
        return pred.atoms(Symbol)

    if not exprs:
        req_keys = find_symbols(proposition)
        keys = proposition.all_predicates()
        # XXX: We need this since True/False are not Basic
        lkeys = set()
        lkeys |= assumptions.all_predicates()
        if context:
            lkeys |= context.all_predicates()

        lkeys = lkeys - {S.true, S.false}
        tmp_keys = None
        while tmp_keys != set():
            tmp = set()
            for l in lkeys:
                syms = find_symbols(l)
                if (syms & req_keys) != set():
                    tmp |= syms
            tmp_keys = tmp - req_keys
            req_keys |= tmp_keys
        keys |= {l for l in lkeys if find_symbols(l) & req_keys != set()}

        exprs = {key.args[0] if isinstance(key, AppliedPredicate) else key for key in keys}
        return exprs, relevant_facts

    for expr in exprs:
        for fact in fact_registry[expr.func]:
            cnf_fact = CNF.to_CNF(fact)
            newfact = cnf_fact.rcall(expr)
            relevant_facts = relevant_facts._and(newfact)
            newexprs |= {key.args[0] for key in newfact.all_predicates()
                             if isinstance(key, AppliedPredicate)}

    return newexprs - exprs, relevant_facts


def get_all_relevant_facts(proposition, assumptions=True,
    context=global_assumptions, use_known_facts=True, iterations=oo):
    # The relevant facts might introduce new keys, e.g., Q.zero(x*y) will
    # introduce the keys Q.zero(x) and Q.zero(y), so we need to run it until
    # we stop getting new things. Hopefully this strategy won't lead to an
    # infinite loop in the future.
    i = 0
    relevant_facts = CNF()
    exprs = None
    all_exprs = set()
    while exprs != set():
        exprs, relevant_facts = get_relevant_facts(proposition,
                assumptions, context, exprs=exprs,
                relevant_facts=relevant_facts)
        all_exprs |= exprs
        i += 1
        if i >= iterations:
            break

    if use_known_facts:
        known_facts_CNF = CNF()
        known_facts_CNF.add_clauses(get_all_known_facts())
        kf_encoded = EncodedCNF()
        kf_encoded.from_cnf(known_facts_CNF)

        def translate_literal(lit, delta):
            if lit > 0:
                return lit + delta
            else:
                return lit - delta

        def translate_data(data, delta):
            return [{translate_literal(i, delta) for i in clause} for clause in data]
        data = []
        symbols = []
        n_lit = len(kf_encoded.symbols)
        for i, expr in enumerate(all_exprs):
            symbols += [pred(expr) for pred in kf_encoded.symbols]
            data += translate_data(kf_encoded.data, i * n_lit)

        encoding = dict(list(zip(symbols, range(1, len(symbols)+1))))
        ctx = EncodedCNF(data, encoding)
    else:
        ctx = EncodedCNF()

    ctx.add_from_cnf(relevant_facts)

    return ctx
