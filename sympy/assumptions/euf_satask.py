from sympy.assumptions.assume import global_assumptions, AppliedPredicate
from sympy.assumptions.cnf import CNF, EncodedCNF
from sympy.assumptions.ask import Q
from sympy.logic.inference import satisfiable
from sympy.logic.algorithms.lra_theory import UnhandledInput, ALLOWED_PRED
from sympy.matrices.kind import MatrixKind
from sympy.core.kind import NumberKind
from sympy.core.symbol import Dummy
from sympy.utilities.iterables import numbered_symbols
from sympy import Lambda, Eq, Not
from sympy.core.mul import Mul
from sympy.core.singleton import S

# Allowed binary preds:
ALLOWED_BIN_PRED = {Q.eq, Q.ne}

def euf_satask(proposition, assumptions=True, context=global_assumptions):
    props = CNF.from_prop(proposition)
    _props = CNF.from_prop(~proposition)

    cnf = CNF.from_prop(assumptions)
    assumptions_enc = EncodedCNF()
    assumptions_enc.from_cnf(cnf)

    context_cnf = CNF()
    if context:
        context_cnf = context_cnf.extend(context)
    assumptions_enc.add_from_cnf(context_cnf)

    return check_satisfiability(props, _props, assumptions_enc)

def check_satisfiability(prop, _prop, factbase):
    sat_true = factbase.copy()
    sat_false = factbase.copy()
    sat_true.add_from_cnf(prop)
    sat_false.add_from_cnf(_prop)

    all_pred, all_exprs = get_all_pred_and_expr_from_enc_cnf(sat_true)

    for pred in all_pred:
        if isinstance(pred, AppliedPredicate):
            if len(pred.arguments) == 1:
                # unary predicate always allowed
                continue
            if pred.function not in ALLOWED_BIN_PRED:
                raise UnhandledInput(f"EUFSolver: {pred} not allowed binary predicate")
        else:
            raise UnhandledInput(f"EUFSolver: unsupported literal {pred}")

    for expr in all_exprs:
        if expr.kind == MatrixKind(NumberKind):
            raise UnhandledInput(f"EUFSolver: {expr} is of MatrixKind")
        if expr == S.NaN:
            raise UnhandledInput("EUFSolver: nan")

    # NOTE: For EUF, we don't do the numeric assumption expansion that LRA does.

    sat_true = _preprocess_euf(sat_true)
    sat_false = _preprocess_euf(sat_false)

    can_be_true = satisfiable(sat_true, use_euf_theory=True) is not False
    can_be_false = satisfiable(sat_false, use_euf_theory=True) is not False

    if can_be_true and can_be_false:
        return None
    if can_be_true and not can_be_false:
        return True
    if not can_be_true and can_be_false:
        return False
    raise ValueError("Inconsistent assumptions")

def _preprocess_euf(enc_cnf):
    """
    Rewrite EncodedCNF so only Q.eq/Q.ne at top level remain.
    Replace complex/unary predicates with dummy constants and equalities.
    """
    enc_cnf = enc_cnf.copy()
    rev_encoding = {v: k for k, v in enc_cnf.encoding.items()}

    new_enc = {}
    new_data = []
    dummies = numbered_symbols("c", Dummy)
    lambda_map = {}

    def lambda_for_term(term):
        lam = _make_curried_lambda(term)
        if lam not in lambda_map:
            dummy = next(dummies)
            lambda_map[lam] = dummy
            _ensure_in_encoding(Eq(lam, dummy))
        return lambda_map[lam]

    def _ensure_in_encoding(expr):
        if expr not in new_enc:
            new_enc[expr] = len(new_enc) + 1
        return new_enc[expr]

    def process_pred(pred):
        if isinstance(pred, AppliedPredicate):
            # Unary predicate
            if len(pred.arguments) == 1:
                arg = pred.arguments[0]
                arg_dummy = lambda_for_term(arg)
                lam_pred = _make_curried_lambda(arg_dummy, pred.function)
                if lam_pred not in lambda_map:
                    dummy = next(dummies)
                    lambda_map[lam_pred] = dummy
                    _ensure_in_encoding(Eq(lam_pred, dummy))
                return lambda_map[lam_pred]
            # Binary eq/ne
            elif pred.function in ALLOWED_BIN_PRED:
                left, right = pred.arguments
                left_dummy = lambda_for_term(left)
                right_dummy = lambda_for_term(right)
                if pred.function == Q.eq:
                    return Eq(left_dummy, right_dummy)
                else:
                    return Not(Eq(left_dummy, right_dummy))
            else:
                raise UnhandledInput(f"EUFSolver: not allowed {pred}")
        else:
            if hasattr(pred, "free_symbols") and pred.free_symbols:
                return lambda_for_term(pred)
            return pred

    for clause in enc_cnf.data:
        new_clause = []
        for lit in clause:
            if lit == 0:
                new_clause.append(0)
                continue
            pred = rev_encoding[abs(lit)]
            sign = 1 if lit > 0 else -1
            processed = process_pred(pred)
            lit_id = _ensure_in_encoding(processed)
            new_clause.append(sign * lit_id)
        new_data.append(new_clause)

    return EncodedCNF(new_data, new_enc)

def _make_curried_lambda(expr, func_or_pred=None):
    vars_sorted = sorted(expr.free_symbols, key=lambda s: s.name)
    lam_expr = func_or_pred(expr) if func_or_pred is not None else expr
    for v in reversed(vars_sorted):
        lam_expr = Lambda(v, lam_expr)
    return lam_expr

def get_all_pred_and_expr_from_enc_cnf(enc_cnf):
    all_exprs = set()
    all_pred = set()
    for pred in enc_cnf.encoding.keys():
        if isinstance(pred, AppliedPredicate):
            all_pred.add(pred)
            all_exprs.update(pred.arguments)
    return all_pred, all_exprs
