from sympy.assumptions.assume import global_assumptions
from sympy.assumptions.cnf import CNF, EncodedCNF
from sympy.assumptions.ask import Q
from sympy.logic.inference import satisfiable
from sympy.logic.algorithms.lra_theory import UnhandledNumber, ALLOWED_PRED
from sympy.assumptions.assume import AppliedPredicate


def lra_satask(proposition, assumptions=True, context=global_assumptions):
    """
    Function to evaluate the proposition with assumptions using SAT algorithm
    in conjunction with an Linear Real Arithmetic theory solver.

    Used to handle inequalities. Should eventually be depreciated and combined
    into satask, but infinity handling and other things need to be implemented
    before that can happen.
    """
    props = CNF.from_prop(proposition)
    _props = CNF.from_prop(~proposition)

    cnf = CNF.from_prop(assumptions)
    assumptions = EncodedCNF()
    assumptions.from_cnf(cnf)

    context_cnf = CNF()
    if context:
        context_cnf = context_cnf.extend(context)

    assumptions.add_from_cnf(context_cnf)

    return check_satisfiability(props, _props, assumptions)


def check_satisfiability(prop, _prop, factbase):
    sat_true = factbase.copy()
    sat_false = factbase.copy()
    sat_true.add_from_cnf(prop)
    sat_false.add_from_cnf(_prop)

    sat_true = split_unequality(sat_true)
    sat_false = split_unequality(sat_false)

    for pred in sat_true.encoding.keys():
        if isinstance(pred, AppliedPredicate):
            if pred.function not in ALLOWED_PRED and pred.function:
                return None
            exprs = pred.arguments
            if any(expr.is_real is not True for expr in exprs):
                return None

    try:
        can_be_true = satisfiable(sat_true, use_lra_theory=True) is not False
        can_be_false = satisfiable(sat_false, use_lra_theory=True) is not False
    except UnhandledNumber:
        return None

    if can_be_true and can_be_false:
        return None

    if can_be_true and not can_be_false:
        return True

    if not can_be_true and can_be_false:
        return False



def split_unequality(enc_cnf):
    """
    Returns an encoded cnf without any Q.ne predicate.

    Converts every unequality into a disjunction of strict
    inequalities. For example, x != 3 would become
    x < 3 OR x > 3.

    Also converts all negated Q.ne predicates into
    equalities.
    """
    enc_cnf = enc_cnf.copy()
    cur_enc = len(enc_cnf.encoding.items()) + 1

    rev_encoding = {value: key for key, value in enc_cnf.encoding.items()}

    new_data = []
    for clause in enc_cnf.data:
        new_clause = []
        for lit in clause:
            prop = rev_encoding[abs(lit)]
            negated = lit < 0

            if not isinstance(prop, AppliedPredicate):
                new_clause.append(lit)
                continue

            if negated and prop.function == Q.eq:
                negated = False
                prop = Q.ne(*prop.arguments)

            if prop.function != Q.ne:
                new_clause.append(lit)
                continue

            if prop.function == Q.ne:
                arg1, arg2 = prop.arguments
                if negated:
                    new_prop = Q.eq(arg1, arg2)
                    if new_prop not in enc_cnf.encoding:
                        enc_cnf.encoding[new_prop] = cur_enc
                        cur_enc += 1

                    new_enc = enc_cnf.encoding[new_prop]
                    new_clause.append(new_enc)
                    continue
                else:
                    new_props = (Q.gt(arg1, arg2), Q.lt(arg1, arg2))
                    for new_prop in new_props:
                        if new_prop not in enc_cnf.encoding:
                            enc_cnf.encoding[new_prop] = cur_enc
                            cur_enc += 1

                        new_enc = enc_cnf.encoding[new_prop]
                        new_clause.append(new_enc)
                    continue

            if prop.function == Q.eq and negated:
                assert False



        new_data.append(new_clause)

    enc_cnf.data = new_data
    return reassign_encoding(enc_cnf)

def reassign_encoding(enc_cnf):
    enc_cnf = enc_cnf.copy()
    old_to_new = {}
    cur_enc = 1
    new_data = []
    for clause in enc_cnf.data:
        new_clause = []
        for lit in clause:
            if abs(lit) not in old_to_new:
                old_to_new[abs(lit)] = cur_enc
                cur_enc += 1
            sign = (lit > 0) - (lit < 0)
            assert sign != 0
            new_lit = old_to_new[abs(lit)]*sign
            new_clause.append(new_lit)
        new_data.append(new_clause)

    new_encoding = {}
    for prop, enc in enc_cnf.encoding.items():
        if enc in old_to_new:
            new_encoding[prop] = old_to_new[enc]

    enc_cnf.encoding = new_encoding
    enc_cnf.data = new_data
    return EncodedCNF(new_data, new_encoding)
