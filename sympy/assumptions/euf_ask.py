from sympy.assumptions.assume import global_assumptions, AppliedPredicate
from sympy.assumptions.cnf import CNF, EncodedCNF
from sympy.assumptions.ask import Q
from sympy.logic.inference import satisfiable
from sympy.logic.algorithms.euf_theory import EUFUnhandledInput
from sympy.matrices.kind import MatrixKind
from sympy.core.kind import NumberKind
from sympy.core.symbol import Dummy, Symbol
from sympy.utilities.iterables import numbered_symbols
from sympy import Lambda
from sympy.core.singleton import S

# Allowed binary preds
ALLOWED_BIN_PRED = {Q.eq, Q.ne}

def euf_ask(proposition, assumptions=True, context=global_assumptions):
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
                continue
            if pred.function not in ALLOWED_BIN_PRED:
                raise EUFUnhandledInput(f"EUFSolver: {pred} not allowed binary predicate")
        else:
            raise EUFUnhandledInput(f"EUFSolver: unsupported literal {pred}")

    for expr in all_exprs:
        if expr.kind == MatrixKind(NumberKind):
            raise EUFUnhandledInput(f"EUFSolver: {expr} is of MatrixKind")
        if expr == S.NaN:
            raise EUFUnhandledInput("EUFSolver: nan")

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
    Replace unary predicates with dummy constants and equalities.

    Example:
    Q.prime(x) becomes Q.eq(Lambda(x, Q.prime(x)), _c0)
    Q.eq(x, y) becomes Q.eq(_c1, _c2) where _c1 represents x and _c2 represents y
    """
    enc_cnf = enc_cnf.copy()
    rev_encoding = {v: k for k, v in enc_cnf.encoding.items()}
    new_enc, new_data = {}, []

    # Create numbered dummy constants for unary predicates
    dummies = numbered_symbols("_c", lambda name=None: Dummy(name, finite=True))
    lambda_map = {}

    def _ensure_in_encoding(expr):
        if expr not in new_enc:
            new_enc[expr] = len(new_enc) + 1
        return new_enc[expr]

    def lambda_for_term(term):
        """Create a Lambda object for a term and map it to a dummy constant"""
        lam = _make_curried_lambda(term)
        if lam not in lambda_map:
            dummy = next(dummies)
            lambda_map[lam] = dummy
            # Create equality: Lambda(term) = dummy
            _ensure_in_encoding(Q.eq(lam, dummy))
        return lam, lambda_map[lam]

    def process_pred(pred):
        if isinstance(pred, AppliedPredicate):
            # Unary predicate - convert to equality with dummy constant
            if len(pred.arguments) == 1:
                arg = pred.arguments[0]
                # Create Lambda(x, Q.prime(x)) and map it to a dummy
                lam_pred = _make_curried_lambda(arg, pred.function)
                if lam_pred not in lambda_map:
                    dummy_pred = next(dummies)
                    lambda_map[lam_pred] = dummy_pred
                    _ensure_in_encoding(Q.eq(lam_pred, dummy_pred))
                return Q.eq(lam_pred, lambda_map[lam_pred])

            # Binary eq/ne - keep as is but replace arguments with dummies
            elif pred.function in ALLOWED_BIN_PRED:
                if(isinstance(pred.arguments[0],Symbol)):
                    left_dummy = pred.arguments[0]
                else:
                    _, left_dummy = lambda_for_term(pred.arguments[0])
                if(isinstance(pred.arguments[1],Symbol)):
                    right_dummy = pred.arguments[1]
                else:
                    _, right_dummy = lambda_for_term(pred.arguments[1])
                if pred.function == Q.eq:
                    return Q.eq(left_dummy, right_dummy)
                else:  # Q.ne
                    return Q.ne(left_dummy, right_dummy)
            else:
                raise EUFUnhandledInput(f"EUFSolver: not allowed {pred}")
        else:
            # Handle non-predicate expressions
            if hasattr(pred, "free_symbols") and pred.free_symbols:
                _, dummy = lambda_for_term(pred)
                return dummy
            return pred

    # Process each clause
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
    """
    Create a curried Lambda function.

    If func_or_pred is provided, it's applied to expr first.
    Then variables are bound in reverse order.
    """
    if func_or_pred is not None:
        # Apply the predicate function to the expression
        lam_expr = func_or_pred(expr)
    else:
        if isinstance(expr, Symbol):
            # For a single symbol, create Lambda(x, x) directly
            return Lambda(expr, expr)
        else:
            lam_expr = Lambda(expr.args,expr)
            curr_lam_expr = lam_expr.curry()
            return curr_lam_expr

    # Get variables and sort them for consistent ordering
    vars_sorted = sorted(expr.free_symbols, key=lambda s: s.name)

    # Bind variables in reverse order (right-to-left currying)
    for v in reversed(vars_sorted):
        lam_expr = Lambda(v, lam_expr)

    return lam_expr

def get_all_pred_and_expr_from_enc_cnf(enc_cnf):
    all_exprs, all_pred = set(), set()
    for pred in enc_cnf.encoding.keys():
        if isinstance(pred, AppliedPredicate):
            all_pred.add(pred)
            all_exprs.update(pred.arguments)
    return all_pred, all_exprs

