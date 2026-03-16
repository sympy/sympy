from __future__ import annotations
from sympy.assumptions.assume import global_assumptions
from sympy.assumptions.cnf import CNF, EncodedCNF
from sympy.assumptions.ask import Q
from sympy.core.numbers import I
from sympy.logic.inference import satisfiable
from sympy.logic.algorithms.lra_theory import UnhandledInput, ALLOWED_PRED
from sympy.matrices.kind import MatrixKind
from sympy.core.kind import NumberKind
from sympy.assumptions.assume import AppliedPredicate
from sympy.core.mul import Mul
from sympy.core.singleton import S


# Some predicates such as Q.prime can't be handled by lra_satask. For example,
# (x > 0) & (x < 1) & Q.prime(x) is unsat but lra_satask would think it was sat.
# WHITE_LIST is a list of predicates that can always be handled.
WHITE_LIST = ALLOWED_PRED | {Q.positive, Q.negative, Q.zero, Q.nonzero, Q.nonpositive, Q.nonnegative,
    Q.negative_infinite, Q.positive_infinite}

REAL_IMPLYING_PREDICATES = WHITE_LIST | {Q.real}

EXTENED_REAL_IMPLYING_PREDICATES = REAL_IMPLYING_PREDICATES | {Q.extended_real, Q.extended_positive, Q.extended_negative, Q.extended_nonpositive,
    Q.extended_negative, Q.extended_nonzero}

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
    assumptions_encoded_cnf = EncodedCNF()
    assumptions_encoded_cnf.from_cnf(cnf)

    context_cnf = CNF()
    if context:
        context_cnf = context_cnf.extend(context)

    assumptions_encoded_cnf.add_from_cnf(context_cnf)


    # Use unit clauses from `assumptions_encoded_cnf` to
    # deduce some expressions are extened real.
    reverse_encoding = {value: key for key, value in assumptions_encoded_cnf.encoding.items()}
    real_exprs = set()
    for unit_clause in assumptions_encoded_cnf.data:
        # Check if `unit_clause` is a unit clause.
        if len(unit_clause) != 1:
            continue

        # Continue if negative because negative literals
        # cannot establish that an expression is entended real.
        (positive_literal,) = unit_clause
        if not (positive_literal > 0):
            continue

        applied_predicate = reverse_encoding[positive_literal]
        if not isinstance(applied_predicate, AppliedPredicate):
            continue

        if applied_predicate.function in REAL_IMPLYING_PREDICATES:
            real_exprs.update(applied_predicate.arguments)

    # SATIFSIABLE STARTED HERE

    factbase = assumptions_encoded_cnf
    sat_true = factbase.copy()
    sat_false = factbase.copy()
    sat_true.add_from_cnf(props)
    sat_false.add_from_cnf(_props)


    # Check old assumptions to:
    #   1. Deduce some expressions are extened real.
    #   2. Add relevant inequality assumptions as unit clauses
    #      to `assumptions_encoded_cnf`.
    all_pred, all_exprs = get_all_pred_and_expr_from_enc_cnf([sat_true])
    for expr in all_exprs:
        convertable_to_inequality_applied_predicate, is_real  = extract_assumption_from_old_assumption(expr)
        if is_real:
            real_exprs.add(expr)
        if convertable_to_inequality_applied_predicate is None:
            continue


        if convertable_to_inequality_applied_predicate not in assumptions_encoded_cnf.encoding:
            applied_predicate_encoding = len(assumptions_encoded_cnf.encoding) + 1
            assumptions_encoded_cnf.encoding[convertable_to_inequality_applied_predicate] = applied_predicate_encoding

        # Add relevant old assumptions to `assumptions`.
        applied_predicate_encoding = assumptions_encoded_cnf.encoding[convertable_to_inequality_applied_predicate]
        assumptions_encoded_cnf.data.append([applied_predicate_encoding])

    for expr in all_exprs:
        if expr not in real_exprs:
            raise UnhandledInput(f"LRASolver: {expr} must be real")


    for pred in all_pred:
        if pred.function not in EXTENED_REAL_IMPLYING_PREDICATES and pred.function != Q.ne:
            raise UnhandledInput(f"LRASolver: {pred} is an unhandled predicate")
    for expr in all_exprs:
        if expr.kind != NumberKind:
            raise UnhandledInput(f"LRASolver: Only scalar expresions are supported. {expr} must be of {NumberKind} but is of {expr.kind} instead.")
        if expr == S.NaN:
            raise UnhandledInput("LRASolver: nan")

    # check satisfiable

    # factbase = assumptions_encoded_cnf
    # sat_true = factbase.copy()
    # sat_false = factbase.copy()
    # sat_true.add_from_cnf(props)
    # sat_false.add_from_cnf(_props)


    # Preprocess sat_true and sat_false into encoded CNFs containing only
    # Q.eq, Q.gt, Q.lt, Q.ge, Q.le, Q.real, Q.extended_real predicates.
    # Converts every unequality into a disjunction of strict inequalities
    # (e.g. x != 3  =>  x < 3 OR x > 3) and converts negated Q.ne into
    # equalities.
    sat_true = _preprocess(sat_true)
    sat_false = _preprocess(sat_false)

    can_be_true = satisfiable(sat_true, use_lra_theory=True) is not False
    can_be_false = satisfiable(sat_false, use_lra_theory=True) is not False

    if can_be_true and can_be_false:
        return None

    if can_be_true and not can_be_false:
        return True

    if not can_be_true and can_be_false:
        return False

    if not can_be_true and not can_be_false:
        raise ValueError("Inconsistent assumptions")


def _preprocess(enc_cnf):
    rev_encoding = {value: key for key, value in enc_cnf.encoding.items()}
    new_encoding = {}
    new_data = []

    def get_enc(p):
        if p not in new_encoding:
            new_encoding[p] = len(new_encoding) + 1
        return new_encoding[p]

    for clause in enc_cnf.data:
        new_clause = []
        for lit in clause:
            if lit == 0:
                new_clause.append(lit)
                new_encoding[lit] = False # Preserve existing DIMACS zero-terminator behavior
                continue

            prop = rev_encoding[abs(lit)]
            is_negated = lit < 0
            sign = -1 if is_negated else 1

            prop = _pred_to_binrel(prop)

            if not isinstance(prop, AppliedPredicate):
                new_clause.append(sign * get_enc(prop))
                continue

            if (prop.function == Q.eq and is_negated) or (prop.function == Q.ne and not is_negated):
                # (x != y) -> (x > y) | (x < y)
                arg1, arg2 = prop.arguments
                new_clause.extend([get_enc(Q.gt(arg1, arg2)), get_enc(Q.lt(arg1, arg2))])
            elif prop.function == Q.ne and is_negated:
                # ~(x != y) -> x == y
                arg1, arg2 = prop.arguments
                new_clause.append(get_enc(Q.eq(arg1, arg2)))
            else:
                # Standard predicates
                assert prop.function in (Q.gt, Q.lt, Q.ge, Q.le, Q.eq)
                new_clause.append(sign * get_enc(prop))

        new_data.append(new_clause)

    return EncodedCNF(new_data, new_encoding)


def _pred_to_binrel(pred):
    if not isinstance(pred, AppliedPredicate):
        return pred

    if pred.function in pred_to_pos_neg_zero:
        f = pred_to_pos_neg_zero[pred.function]
        if f is False:
            return False
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
    Q.extended_positive: Q.positive,
    Q.extended_negative: Q.negative,
    Q.extended_nonpositive: Q.nonpositive,
    Q.extended_negative: Q.negative,
    Q.extended_nonzero: Q.nonzero,
    Q.negative_infinite: False,
    Q.positive_infinite: False
}


def get_all_pred_and_expr_from_enc_cnf(enc_cnfs):
    all_exprs = set()
    all_pred = set()
    for enc_cnf in enc_cnfs:
        for pred in enc_cnf.encoding.keys():
            if isinstance(pred, AppliedPredicate):
                all_pred.add(pred)
                all_exprs.update(pred.arguments)

    return all_pred, all_exprs


def extract_assumption_from_old_assumption(expr):
    """
    Extracts a single relevant assumption from the old assumption system
    if there are any. Otherwise, it either returns None or raises an exception
    if a dissallowed assumption is present.

    Gives a

    Example
    =======
    >>> from sympy.assumptions.lra_satask import extract_assumption_from_old_assumption
    >>> from sympy import symbols
    >>> x = symbols("x")
    >>> extract_assumption_from_old_assumption(x) is None
    True
    >>> y = symbols("y", positive=True)
    >>> extract_assumption_from_old_assumption(y)
    Q.positive(y)
    >>> extract_assumption_from_old_assumption(-y)
    Q.negative(-y)
    >>> y = symbols("y", integer=True)
    >>> extract_assumption_from_old_assumption(y) # raises exception
    """
    if not hasattr(expr, "free_symbols") or len(expr.free_symbols) == 0:
        if expr.is_real:
            return None, True
        return None, None

    # Test for I times imaginary variable. Such expressions are considered
    # real but aren't handled.
    if expr.has(I):
        raise UnhandledInput(f"LRASolver: {expr} must not contain I")

    if expr.is_integer == True and expr.is_zero != True:
        raise UnhandledInput(f"LRASolver: {expr} is an integer")
    if expr.is_integer == False:
        raise UnhandledInput(f"LRASolver: {expr} can't be an integer")
    if expr.is_rational == False:
        raise UnhandledInput(f"LRASolver: {expr} is irational")

    if expr.is_zero:
        return Q.zero(expr), True
    if expr.is_positive:
        return Q.positive(expr), True
    if expr.is_negative:
        return Q.negative(expr), True
    if expr.is_nonzero:
        return Q.nonzero(expr), True
    if expr.is_nonpositive:
        return Q.nonpositive(expr), True
    if expr.is_nonnegative:
        return Q.nonnegative(expr), True
    if expr.is_real:
        return None, True

    return None, None


def extract_pred_from_old_assum(all_exprs):
    """
    Returns a list of relevant new assumption predicate
    based on any old assumptions.

    Raises an UnhandledInput exception if any of the assumptions are
    unhandled.

    Ignored predicate:
    - commutative
    - complex
    - algebraic
    - transcendental
    - extended_real
    - real
    - all matrix predicate
    - rational
    - irrational

    Example
    =======
    >>> from sympy.assumptions.lra_satask import extract_pred_from_old_assum
    >>> from sympy import symbols
    >>> x, y = symbols("x y", positive=True)
    >>> extract_pred_from_old_assum([x, y, 2])
    [Q.positive(x), Q.positive(y)]
    """
    ret = []
    for expr in all_exprs:
        if not hasattr(expr, "free_symbols"):
            continue
        if len(expr.free_symbols) == 0:
            continue

        if expr.is_real is not True:
            raise UnhandledInput(f"LRASolver: {expr} must be real")
        # test for I times imaginary variable; such expressions are considered real
        if isinstance(expr, Mul) and any(arg.is_real is not True for arg in expr.args):
            raise UnhandledInput(f"LRASolver: {expr} must be real")

        if expr.is_integer == True and expr.is_zero != True:
            raise UnhandledInput(f"LRASolver: {expr} is an integer")
        if expr.is_integer == False:
            raise UnhandledInput(f"LRASolver: {expr} can't be an integer")
        if expr.is_rational == False:
            raise UnhandledInput(f"LRASolver: {expr} is irational")

        if expr.is_zero:
            ret.append(Q.zero(expr))
        elif expr.is_positive:
            ret.append(Q.positive(expr))
        elif expr.is_negative:
            ret.append(Q.negative(expr))
        elif expr.is_nonzero:
            ret.append(Q.nonzero(expr))
        elif expr.is_nonpositive:
            ret.append(Q.nonpositive(expr))
        elif expr.is_nonnegative:
            ret.append(Q.nonnegative(expr))

    return ret
