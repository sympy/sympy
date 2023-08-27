from sympy.assumptions.assume import global_assumptions
from sympy.assumptions.cnf import CNF, EncodedCNF
from sympy.assumptions.ask import _ask_single_fact, Q
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

    check_satisfiability(props, _props, assumptions)


def check_satisfiability(prop, _prop, factbase):
    sat_true = factbase.copy()
    sat_false = factbase.copy()
    sat_true.add_from_cnf(prop)
    sat_false.add_from_cnf(_prop)

    for pred in sat_true.encoding.keys():
        if isinstance(pred, AppliedPredicate):
            if pred.function not in ALLOWED_PRED:
                return None
            exprs = pred.arguments
            if any(expr.is_finite is not True for expr in exprs):
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
