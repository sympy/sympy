from __future__ import annotations
from sympy.assumptions.assume import global_assumptions
from sympy.assumptions.cnf import CNF, EncodedCNF
from sympy.assumptions.ask import Q
from sympy.logic.inference import satisfiable
from sympy.logic.algorithms.lra_theory import UnhandledInput, ALLOWED_PRED
from sympy.matrices.kind import MatrixKind
from sympy.core.kind import NumberKind
from sympy.assumptions.assume import AppliedPredicate
from sympy.core.mul import Mul
from sympy.core.singleton import S


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

# Some predicates such as Q.prime can't be handled by lra_satask.
# For example, (x > 0) & (x < 1) & Q.prime(x) is unsat but lra_satask would think it was sat.
# WHITE_LIST is a list of predicates that can always be handled.
WHITE_LIST = ALLOWED_PRED.keys() | {Q.positive, Q.negative, Q.zero, Q.nonzero, Q.nonpositive, Q.nonnegative,
                                    Q.extended_positive, Q.extended_negative, Q.extended_nonpositive,
                                    Q.extended_negative, Q.extended_nonzero, Q.negative_infinite,
                                    Q.positive_infinite}


def check_satisfiability(prop, _prop, factbase):
    sat_true = factbase.copy()
    sat_false = factbase.copy()
    sat_true.add_from_cnf(prop)
    sat_false.add_from_cnf(_prop)

    all_pred, all_exprs = get_all_pred_and_expr_from_enc_cnf(sat_true)

    # Predicates that are informational only and can be safely ignored by the
    # LRA solver. Q.real in particular was being raised as an unhandled input
    # even though it is documented as an "Ignored predicate" -- this caused
    # inconsistent behaviour for ask() depending on whether reality was
    # supplied as Q.real(x) or as symbols(..., real=True). See #27834.
    SILENT_IGNORE = {Q.real, Q.complex, Q.irrational, Q.rational, Q.algebraic,
                     Q.transcendental, Q.extended_real, Q.commutative}

    for pred in all_pred:
        if (pred.function not in WHITE_LIST and pred.function != Q.ne
                and pred.function not in SILENT_IGNORE):
            raise UnhandledInput(f"LRASolver: {pred} is an unhandled predicate")
    for expr in all_exprs:
        if expr.kind == MatrixKind(NumberKind):
            raise UnhandledInput(f"LRASolver: {expr} is of MatrixKind")
        if expr == S.NaN:
            raise UnhandledInput("LRASolver: nan")

    # Infer real variables from relational assumptions like Q.eq(x, 1) or
    # Q.lt(x, 0) -- see https://github.com/sympy/sympy/issues/27834.
    inferred_real = _infer_real_variables(factbase)
    # convert old assumptions into predicates and add them to sat_true and sat_false
    # also check for unhandled predicates
    for assm in extract_pred_from_old_assum(all_exprs, inferred_real=inferred_real):
        n = len(sat_true.encoding)
        if assm not in sat_true.encoding:
            sat_true.encoding[assm] = n+1
        sat_true.data.append([sat_true.encoding[assm]])

        n = len(sat_false.encoding)
        if assm not in sat_false.encoding:
            sat_false.encoding[assm] = n+1
        sat_false.data.append([sat_false.encoding[assm]])


    # Strip the silent-ignore predicates (Q.real etc.) from the encoding
    # because the downstream LRA solver asserts all predicates are
    # Q.eq/Q.lt/Q.gt/Q.le/Q.ge. We have already extracted the reality info
    # we need from them into ``inferred_real`` above.
    _strip_ignored_predicates(sat_true)
    _strip_ignored_predicates(sat_false)

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
    """
    Returns an encoded cnf with only Q.eq, Q.gt, Q.lt,
    Q.ge, and Q.le predicate.

    Converts every unequality into a disjunction of strict
    inequalities. For example, x != 3 would become
    x < 3 OR x > 3.

    Also converts all negated Q.ne predicates into
    equalities.
    """

    # loops through each literal in each clause
    # to construct a new, preprocessed encodedCNF

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

            prop = _pred_to_binrel(prop)

            if not isinstance(prop, AppliedPredicate):
                if prop not in new_encoding:
                    new_encoding[prop] = cur_enc
                    cur_enc += 1
                lit = new_encoding[prop]
                new_clause.append(sign*lit)
                continue


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

def get_all_pred_and_expr_from_enc_cnf(enc_cnf):
    all_exprs = set()
    all_pred = set()
    for pred in enc_cnf.encoding.keys():
        if isinstance(pred, AppliedPredicate):
            all_pred.add(pred)
            all_exprs.update(pred.arguments)

    return all_pred, all_exprs


def _infer_real_variables(enc_cnf):
    """Infer which variables must be real from relational assumptions.

    Walks the encoded factbase and returns the set of free Symbols that are
    forced to be real by the given assumptions. The rules are:

    - ``Q.eq(v, c)`` or ``Q.eq(c, v)`` where ``c`` is a real constant implies
      ``v`` is real (since ``c`` is real and equality is transitive with
      reals).
    - ``Q.lt`` / ``Q.gt`` / ``Q.le`` / ``Q.ge`` between a variable and a real
      constant also implies the variable is real, because strict/non-strict
      ordering is only defined on the reals.

    This lets ``ask`` work correctly when the user provides the realness of a
    variable indirectly through a relational, rather than via ``Q.real(v)``
    or ``symbols(..., real=True)`` -- see sympy/sympy#27834.

    Parameters
    ==========

    enc_cnf : EncodedCNF
        The encoded factbase of assumptions to scan.

    Returns
    =======

    set of Symbol
        The set of variables that are forced to be real by the assumptions.

    Examples
    ========

    >>> from sympy.assumptions.lra_satask import _infer_real_variables
    >>> from sympy.assumptions.cnf import CNF, EncodedCNF
    >>> from sympy import Q, symbols
    >>> y = symbols('y')
    >>> cnf = CNF.from_prop(Q.eq(y, 1))
    >>> enc = EncodedCNF()
    >>> enc.from_cnf(cnf)
    >>> _infer_real_variables(enc)
    {y}
    """
    inferred = set()
    for pred in enc_cnf.encoding.keys():
        if not isinstance(pred, AppliedPredicate):
            continue
        if pred.function not in (Q.eq, Q.lt, Q.gt, Q.le, Q.ge):
            continue
        args = pred.arguments
        if len(args) != 2:
            continue
        lhs, rhs = args
        # Equality / strict-or-non-strict ordering between a free symbol and
        # a real constant forces the symbol to be real.
        if (isinstance(lhs, AppliedPredicate)):
            continue
        if getattr(lhs, "is_real", None) is True and hasattr(rhs, "free_symbols"):
            inferred.update(rhs.free_symbols)
        if getattr(rhs, "is_real", None) is True and hasattr(lhs, "free_symbols"):
            inferred.update(lhs.free_symbols)
    return inferred


def _strip_ignored_predicates(enc_cnf):
    """Remove silently-ignored predicates (Q.real, Q.complex, ...) from the
    encoded CNF. The reality info from Q.real(v) should already have been
    extracted via :func:`_infer_real_variables` into the ``inferred_real`` set
    before this function is called.

    The LRA solver's preprocessing asserts that every literal refers to one
    of Q.eq/Q.lt/Q.gt/Q.le/Q.ge, so any other predicate has to be removed
    before the CNF is handed to it.
    """
    silent = {Q.real, Q.complex, Q.irrational, Q.rational, Q.algebraic,
              Q.transcendental, Q.extended_real, Q.commutative}
    ignored_codes = set()
    for pred, code in list(enc_cnf.encoding.items()):
        if isinstance(pred, AppliedPredicate) and pred.function in silent:
            ignored_codes.add(abs(code))
    if not ignored_codes:
        return
    # Remove from encoding
    for pred, code in list(enc_cnf.encoding.items()):
        if abs(code) in ignored_codes:
            del enc_cnf.encoding[pred]
    # Remove from data (clauses). If a clause becomes empty after removal, the
    # factbase is unsatisfiable -- but this should not happen for well-formed
    # inputs because the silent-ignore predicates are always positive and
    # we only strip clauses consisting of a single such literal.
    new_data = []
    for clause in enc_cnf.data:
        new_clause = [lit for lit in clause if abs(lit) not in ignored_codes]
        if new_clause:
            new_data.append(new_clause)
    enc_cnf.data = new_data


def extract_pred_from_old_assum(all_exprs, inferred_real=()):
    """
    Returns a list of relevant new assumption predicate
    based on any old assumptions.

    *inferred_real* is an optional iterable of free Symbols that have been
    forced to be real by the relational assumptions (see :func:`_infer_real_variables`).
    These symbols are accepted as real even if their ``is_real`` flag is not
    explicitly set -- this is what allows ``ask(Q.zero(x - 1), Q.eq(x, 1))``
    to return ``True`` for an unconstrained symbol ``x``.

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
    inferred_real_set = set(inferred_real)
    for expr in all_exprs:
        if not hasattr(expr, "free_symbols"):
            continue
        if len(expr.free_symbols) == 0:
            continue

        # A Symbol is treated as real if it is real by flag, OR if it was
        # forced to be real by a relational assumption (e.g. Q.eq(x, 1)).
        is_expr_real = expr.is_real is True or (
            expr.free_symbols and expr.free_symbols.issubset(inferred_real_set)
        )
        if not is_expr_real:
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
