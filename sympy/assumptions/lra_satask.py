from __future__ import annotations
from sympy.assumptions.assume import global_assumptions
from sympy.assumptions.cnf import CNF, EncodedCNF
from sympy.assumptions.ask import Q
from sympy.core.add import Add
from sympy.core.numbers import I
from sympy.core.symbol import Symbol
from sympy.logic.inference import satisfiable
from sympy.logic.algorithms.lra_theory import UnhandledInput, ALLOWED_PRED
from sympy.matrices.kind import MatrixKind
from sympy.core.kind import NumberKind
from sympy.assumptions.assume import AppliedPredicate
from sympy.core.mul import Mul
from sympy.core.singleton import S


REAL_PREDICATES = {Q.real, Q.positive, Q.negative, Q.zero, Q.nonzero, Q.nonpositive, Q.nonnegative}
EXTENDED_REAL_PREDICATES = {Q.extended_real, Q.extended_positive, Q.extended_negative, Q.extended_nonpositive,
    Q.extended_nonnegative, Q.extended_nonzero, Q.positive_infinite, Q.negative_infinite, Q.gt, Q.lt, Q.ge, Q.le}


HANDLED_PREDICATES = EXTENED_REAL_IMPLYING_PREDICATES | {Q.ne}

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
    # try to deduce if expressions are real or extended real.
    reverse_encoding = {value: key for key, value in assumptions_encoded_cnf.encoding.items()}
    real_exprs = set()
    extended_real_exprs = set()
    for unit_clause in assumptions_encoded_cnf.data:
        # Check if `unit_clause` is a unit clause.
        if len(unit_clause) != 1:
            continue

        # Negated predicates cannot establish that
        # an expression is real (or extended real).
        (positive_literal,) = unit_clause
        if not (positive_literal > 0):
            continue

        applied_predicate = reverse_encoding[positive_literal]
        if not isinstance(applied_predicate, AppliedPredicate):
            continue

        if applied_predicate.function in REAL_PREDICATES:
            real_exprs.update(applied_predicate.arguments)
        if applied_predicate.function in EXTENDED_REAL_PREDICATES:
            extended_real_exprs.update(applied_predicate.arguments)


    # SATIFSIABLE STARTED HERE

    factbase = assumptions_encoded_cnf
    sat_true = factbase.copy()
    sat_false = factbase.copy()
    sat_true.add_from_cnf(props)
    sat_false.add_from_cnf(_props)


    # Check old assumptions to:
    #   1. Check if expressions are assumed to be real or extened real.
    #   2. Add relevant inequality assumptions to `assumptions_encoded_cnf`.
    all_pred, all_exprs = get_all_pred_and_expr_from_enc_cnf([sat_true])
    for expr in all_exprs:
        # If there are any predicates that can be represented as inequalities,
        # `relevant_applied_predicate` is set to the most specific of those predicates.
        relevant_applied_predicate, is_real, is_extended_real  = extract_assumption_from_old_assumption(expr)
        if is_real:
            real_exprs.add(expr)
        elif is_extended_real:
            extended_real_exprs.add(expr)

        if relevant_applied_predicate is None:
            continue


        if relevant_applied_predicate not in assumptions_encoded_cnf.encoding:
            applied_predicate_encoding = len(assumptions_encoded_cnf.encoding) + 1
            assumptions_encoded_cnf.encoding[relevant_applied_predicate] = applied_predicate_encoding

        # Add relevant old assumptions to `assumptions`.
        applied_predicate_encoding = assumptions_encoded_cnf.encoding[relevant_applied_predicate]
        assumptions_encoded_cnf.data.append([applied_predicate_encoding])


<<<<<<< HEAD
    all_symbols = set()
=======

    for pred in all_pred:
        if pred.function not in HANDLED_PREDICATES:
            raise UnhandledInput(f"LRASolver: {pred} is an unhandled predicate")
>>>>>>> 70850820c1 (Create HANDLED_PREDICATES for readablity)
    for expr in all_exprs:
        if expr.kind != NumberKind:
            raise UnhandledInput(f"LRASolver: Only scalar expresions are supported. {expr} must be of {NumberKind} but is of {expr.kind} instead.")
        if expr == S.NaN:
            raise UnhandledInput("LRASolver: nan")
        if expr not in real_exprs and expr not in extended_real_exprs:
            raise UnhandledInput(f"LRASolver: {expr} must be extended real.")

        # If there exist distinct extended real expresions (that may be infinite),
        # they must be variable disjoint from each other because extended real
        # arthmatic is not fully handled. For example, the query `ask(x > x + 1)`
        # is not handled because `x` and `x + 1` are not variable disjoint.
        # We have to be careful because `x > x + 1` may be False if x=oo.
        if expr in extended_real_exprs and expr not in real_exprs:
            if all_symbols.isdisjoint(expr.free_symbols):
                all_symbols.update(expr.free_symbols)
                continue

            raise UnhandledInput(f"LRASolver: Extended real expressions are not variable disjoint.")


    for pred in all_pred:
        if (pred.function not in REAL_PREDICATES and
            pred.function not in EXTENDED_REAL_PREDICATES and
            pred.function not in (Q.ne, Q.eq)):
                raise UnhandledInput(f"LRASolver: {pred} is an unhandled predicate")


    # check satisfiable

    # factbase = assumptions_encoded_cnf
    # sat_true = factbase.copy()
    # sat_false = factbase.copy()
    # sat_true.add_from_cnf(props)
    # sat_false.add_from_cnf(_props)


    # Preprocess sat_true and sat_false into encoded CNFs containing only
    # Q.eq, Q.gt, Q.lt, Q.ge, Q.le, Q.eq, and Q.real predicates. It does the following
    # transformations:
    # - Q.positive -> Q.gt, Q.negative -> Q.lt, Q.zero -> Q.eq, etc.
    # - Q.extended_positive -> Q.gt, Q.extended_negative -> Q.lt, etc.
    # - Q.ne, ~Q.eq, Q.nonzero -> disjunction of inequalities (e.g. x != 3  =>  x < 3 OR x > 3)
    # - Q.extended_real -> True
    # - Q.real -> True or Q.real depending on whether the expr is known to be real
    # and converts negated Q.ne into equalities.
    sat_true = _preprocess(sat_true, real_exprs)
    sat_false = _preprocess(sat_false, real_exprs)

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


def _preprocess(enc_cnf, real_exprs):
    """
    Examples
    ========
    >>> from sympy import symbols, Q
    >>> from sympy.assumptions.cnf import EncodedCNF
    >>> from sympy.assumptions.lra_satask import _preprocess
    >>> x, y = symbols('x y')

    >>> r = _preprocess(EncodedCNF([[1]], {Q.positive(x): 1}), {x}, set())
    >>> Q.gt(x, 0) in r.encoding and Q.positive(x) not in r.encoding
    True

    >>> r = _preprocess(EncodedCNF([[1]], {Q.nonzero(x): 1}), {x}, set())
    >>> Q.gt(x, 0) in r.encoding and Q.lt(x, 0) in r.encoding and len(r.data[0]) == 2
    True

    >>> r = _preprocess(EncodedCNF([[-1]], {Q.eq(x, y): 1}), {x, y}, set())
    >>> Q.gt(x, y) in r.encoding and Q.lt(x, y) in r.encoding and len(r.data[0]) == 2
    True

    >>> r = _preprocess(EncodedCNF([[-1]], {Q.ne(x, y): 1}), {x, y}, set())
    >>> Q.eq(x, y) in r.encoding and len(r.data[0]) == 1
    True

    >>> Q.gt(x, 0) in _preprocess(EncodedCNF([[1]], {Q.extended_positive(x): 1}), set(), {x}).encoding
    True

    >>> r = _preprocess(EncodedCNF([[1]], {Q.real(x): 1}), {x}, set())
    >>> True in r.encoding and Q.real(x) not in r.encoding
    True

    >>> Q.real(x) in _preprocess(EncodedCNF([[1]], {Q.real(x): 1}), set(), {x}).encoding
    True

    >>> r = _preprocess(EncodedCNF([[1], [2]], {Q.positive(x): 1, Q.negative(y): 2}), {x, y}, set())
    >>> Q.gt(x, 0) in r.encoding and Q.lt(y, 0) in r.encoding and len(r.data) == 2
    True
    """
    rev_encoding = {value: key for key, value in enc_cnf.encoding.items()}
    new_encoding = {}
    new_data = []

    def get_enc(p):
        if p not in new_encoding:
            new_encoding[p] = len(new_encoding) + 1
        return new_encoding[p]

    for clause in enc_cnf.data:
        new_clause = []
        skip_new_clause = False
        for lit in clause:
            if lit == 0:
                # TODO: Add tests for this and look into if it's actually needed.
                new_clause.append(lit)
                new_encoding[lit] = False # Preserve existing DIMACS zero-terminator behavior
                continue

            prop = rev_encoding[abs(lit)]
            is_negated = lit < 0
            sign = -1 if is_negated else 1

            if not isinstance(prop, AppliedPredicate):
                new_clause.append(sign * get_enc(prop))
                continue

            if prop.function == Q.real and prop.arguments[0] not in real_exprs:
                new_clause.append(sign * get_enc(prop))
                continue
            elif prop.function in (Q.real, Q.extended_real):
                prop = True

            if prop not in (True, False) and prop.function in (Q.positive_infinite, Q.negative_infinite):
                if prop.arguments[0] in real_exprs:
                    prop = False
                else:
                    raise UnhandledInput(f"LRASolver: Q.positive_infinite and Q.negative_infinite are not fully handled yet.")

            if prop in (True, False):
                new_clause.append(sign * get_enc(prop))
                continue

            # if prop is False:
            #     # This might result in new_clause to be empty
            #     # which is intended. In that case, the sat
            #     # solver will give unsat.
            #     continue
            # if prop is True:
            #     # This could result in new_data to be empty which
            #     # is intended. In that case, the sat solver will
            #     # give sat.
            #     skip_new_clause = True
            #     break

            prop = _pred_to_binrel(prop)

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
                if prop.function not in (Q.gt, Q.lt, Q.ge, Q.le, Q.eq):
                    pass
                assert prop.function in (Q.gt, Q.lt, Q.ge, Q.le, Q.eq)
                new_clause.append(sign * get_enc(prop))

        if skip_new_clause:
            continue
        new_data.append(new_clause)

    return EncodedCNF(new_data, new_encoding)


def _pred_to_binrel(pred):
    if not isinstance(pred, AppliedPredicate):
        return pred

    if pred.function in (Q.positive, Q.extended_positive):
        pred = Q.gt(pred.arguments[0], 0)
    elif pred.function in (Q.negative, Q.extended_negative):
        pred = Q.lt(pred.arguments[0], 0)
    elif pred.function == Q.zero:
        pred = Q.eq(pred.arguments[0], 0)
    elif pred.function in (Q.nonpositive, Q.extended_nonpositive):
        pred = Q.le(pred.arguments[0], 0)
    elif pred.function in (Q.nonnegative, Q.extended_nonnegative):
        pred = Q.ge(pred.arguments[0], 0)
    elif pred.function in (Q.nonzero, Q.extended_nonzero):
        pred = Q.ne(pred.arguments[0], 0)

    return pred


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
    # Test for I times imaginary variable. Such expressions are considered
    # real but aren't handled.
    if expr.has(I):
        raise UnhandledInput(f"LRASolver: {expr} must not contain I")

    if not hasattr(expr, "free_symbols") or len(expr.free_symbols) == 0:
        if expr.is_real:
            return None, True, True
        if expr.is_extended_real:
            return None, None, True
        return None, None, None

    if expr.is_integer == True and expr.is_zero != True:
        raise UnhandledInput(f"LRASolver: {expr} is an integer")
    if expr.is_integer == False:
        raise UnhandledInput(f"LRASolver: {expr} can't be an integer")
    if expr.is_rational == False:
        raise UnhandledInput(f"LRASolver: {expr} is irational")

    if expr.is_zero:
        return Q.zero(expr), True, True
    if expr.is_positive:
        return Q.positive(expr), True, True
    if expr.is_negative:
        return Q.negative(expr), True, True
    if expr.is_nonzero:
        return Q.nonzero(expr), True, True
    if expr.is_nonpositive:
        return Q.nonpositive(expr), True, True
    if expr.is_nonnegative:
        return Q.nonnegative(expr), True, True
    if expr.is_real:
        return None, True, True
    if expr.is_extended_real:
        return None, None, True

    return None, None, None
