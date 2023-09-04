from sympy.assumptions.assume import global_assumptions
from sympy.assumptions.cnf import CNF, EncodedCNF
from sympy.assumptions.ask import Q
from sympy.logic.inference import satisfiable
from sympy.logic.algorithms.lra_theory import UnhandledInput, ALLOWED_PRED
from sympy.matrices.common import MatrixKind
from sympy.core.kind import NumberKind
from sympy.assumptions.assume import AppliedPredicate
from sympy.core.mul import Mul


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

    try:
        sat_true = _preprocess(sat_true)
        sat_false = _preprocess(sat_false)
    except UnhandledPred:
        return None

    for pred in sat_true.encoding.keys():
        if isinstance(pred, AppliedPredicate):
            if pred.function not in ALLOWED_PRED:
                return None
            exprs = pred.arguments
            for expr in exprs:
                if expr.kind == MatrixKind(NumberKind):
                    return None
                if expr.is_real is not True:
                    return None
                if isinstance(expr, Mul) and any(arg.is_real is not True for arg in expr.args):
                    return None

    try:
        can_be_true = satisfiable(sat_true, use_lra_theory=True) is not False
        can_be_false = satisfiable(sat_false, use_lra_theory=True) is not False
    except UnhandledInput:
        return None

    if can_be_true and can_be_false:
        return None

    if can_be_true and not can_be_false:
        return True

    if not can_be_true and can_be_false:
        return False


class UnhandledPred(Exception):
    pass


def _preprocess(enc_cnf):
    """
    Returns an encoded cnf with only Q.eq, Q.gt, Q.lt,
    Q.ge, and Q.le predicate.

    Converts every unequality into a disjunction of strict
    inequalities. For example, x != 3 would become
    x < 3 OR x > 3.

    Also converts all negated Q.ne predicates into
    equalities.
    """
    enc_cnf = enc_cnf.copy()
    cur_enc = 1
    rev_encoding = {value: key for key, value in enc_cnf.encoding.items()}

    new_encoding = {}
    new_data = []
    for clause in enc_cnf.data:
        new_clause = []
        for lit in clause:
            if lit == 0:
                new_clause.append(lit)
                new_encoding[lit] = False
                continue
            prop = rev_encoding[abs(lit)]
            negated = lit < 0
            sign = (lit > 0) - (lit < 0)

            if not isinstance(prop, AppliedPredicate):
                if prop not in new_encoding:
                    new_encoding[prop] = cur_enc
                    cur_enc += 1
                lit = new_encoding[prop]
                new_clause.append(sign*lit)
                continue

            prop = _pred_to_binrel(prop)

            if not (prop.function in ALLOWED_PRED or prop.function == Q.ne):
                raise UnhandledPred

            if negated and prop.function == Q.eq:
                negated = False
                prop = Q.ne(*prop.arguments)

            if prop.function == Q.ne:
                arg1, arg2 = prop.arguments
                if negated:
                    new_prop = Q.eq(arg1, arg2)
                    if new_prop not in new_encoding:
                        new_encoding[new_prop] = cur_enc
                        cur_enc += 1

                    new_enc = new_encoding[new_prop]
                    new_clause.append(new_enc)
                    continue
                else:
                    new_props = (Q.gt(arg1, arg2), Q.lt(arg1, arg2))
                    for new_prop in new_props:
                        if new_prop not in new_encoding:
                            new_encoding[new_prop] = cur_enc
                            cur_enc += 1

                        new_enc = new_encoding[new_prop]
                        new_clause.append(new_enc)
                    continue

            if prop.function == Q.eq and negated:
                assert False

            if prop not in new_encoding:
                new_encoding[prop] = cur_enc
                cur_enc += 1
            new_clause.append(new_encoding[prop]*sign)
        new_data.append(new_clause)

    assert len(new_encoding) >= cur_enc - 1

    enc_cnf = EncodedCNF(new_data, new_encoding)
    return enc_cnf


def _pred_to_binrel(pred):
    if pred.function in pred_to_pos_neg_zero:
        f = pred_to_pos_neg_zero[pred.function]
        pred = f(pred.arguments[0])

    if pred.function == Q.positive:
        pred = Q.gt(pred.arguments[0], 0)
    elif pred.function == Q.negative:
        pred = Q.lt(pred.arguments[0], 0)
    elif pred.function == Q.zero:
        pred = Q.eq(pred.arguments[0], 0)
    elif pred.function == Q.nonpositive:
        pred = Q.le(pred.arguments[0], 0)
    elif pred.function == Q.nonnegative:
        pred = Q.ge(pred.arguments[0], 0)
    elif pred.function == Q.nonzero:
        pred = Q.ne(pred.arguments[0], 0)

    return pred

pred_to_pos_neg_zero = {
    Q.prime : Q.positive,
    Q.composite : Q.positive,
    Q.extended_positive : Q.positive,
    Q.extended_negative : Q.negative,
    Q.extended_nonpositive : Q.nonpositive,
    Q.extended_negative : Q.negative,
    Q.extended_nonzero: Q.nonzero
}
