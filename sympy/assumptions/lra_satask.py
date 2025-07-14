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
WHITE_LIST = ALLOWED_PRED | {Q.positive, Q.negative, Q.zero, Q.nonzero, Q.nonpositive, Q.nonnegative,
                                            Q.extended_positive, Q.extended_negative, Q.extended_nonpositive,
                                            Q.extended_negative, Q.extended_nonzero, Q.negative_infinite,
                                            Q.positive_infinite}


def check_satisfiability(prop, _prop, factbase):
    """
    Determines the satisfiability of a logical proposition within the linear real arithmetic (LRA) domain,
    given a base of known assumptions.

    This function evaluates both a logical proposition (`prop`) and its negation (`_prop`) to determine
    whether the proposition is always true, always false, or indeterminate under a given set of
    facts (`factbase`). It is an extension of the baseline `check_satisfiability()` in
    `sympy.assumptions.satask`, but adds domain-specific logic for linear real arithmetic (LRA),
    with stricter type enforcement and domain reasoning.

    Compared to `sympy.assumptions.satask.check_satisfiability`, this version introduces:

    1. Domain Restriction to LRA (Linear Real Arithmetic):
        - Only predicates and expressions involving real or rational scalar values are supported.
        - Types such as matrices, complex values, `NaN`, or symbolic constructs not known to be real
          are explicitly disallowed and will raise errors.
        - This ensures that satisfiability checking is compatible with solvers or theories designed
          for linear real arithmetic.

    2. Realness Enforcement:
        - The function ensures that all symbols involved in arithmetic comparisons are explicitly
          assumed to be real where required.

    3. Preprocessing and Implicit Assumption Injection:
        - Realness assumptions are introduced only under the following conditions:
            - Inequalities (Q.gt, Q.lt): If both sides are single symbolic variables and lack explicit domains, Q.real(x) and Q.real(y) are added.
            - Equalities (Q.eq): If one side is a single symbol known to be real and the other is a single symbol with unknown domain, realness is inferred. If one is real and the other explicitly not real, a ValueError is raised.
            - Explicit declarations: Symbols explicitly declared real (e.g., x.is_real == True) are recognized and used to avoid redundant inference.
        - This preprocessing helps the solver reason with more complete context, avoiding false
          negatives due to missing domain constraints.

    4. Early Rejection of Unsupported Inputs:
        - Before attempting satisfiability checks, the function scans the input for invalid constructs,
          such as:
            - `NaN` values
            - Matrix-like expressions
            - Unsupported predicate functions (anything not in the supported LRA predicate list)
        - If such constructs are found, the function raises `UnhandledInput` instead of continuing.
        - This guards against undefined or nonsensical logical evaluations.

    Parameters
    ----------
    prop : CNF
        The CNF form of the proposition to test for satisfiability.
    _prop : CNF
        The CNF form of the negation of the proposition.
    factbase : CNFEncoding
        A base of facts and assumptions to be included in the satisfiability checking.

    Returns
    -------
    True if the proposition can only be true.
    False if the proposition can only be false.
    None if the proposition can be both true and false under current assumptions.

    Raises
    ------
    ValueError
        If there are inconsistent assumptions (e.g., real vs. non-real variables).
    UnhandledInput
        If the expression contains unsupported constructs or is not entirely real/rational.
    """
    sat_true = factbase.copy()
    sat_false = factbase.copy()
    sat_true.add_from_cnf(prop)
    sat_false.add_from_cnf(_prop)

    all_exprs = set()
    all_pred = set()
    realness_preds = set()

    # Deduce realness with heuristics and collect all_exprs and all_pred
    for pred in sat_true.encoding.keys():
        if not isinstance(pred, AppliedPredicate):
            continue

        if pred.function is Q.eq:
            lhs, rhs = pred.arguments
            # heuristic: x = y & real(x) → add real(y)
            for a, b in [(lhs, rhs), (rhs, lhs)]:
                if a.is_real:
                    if b.is_real is None:
                        realness_preds.add(Q.real(b))
                    elif b.is_real is False:
                        raise ValueError("Inconsistent assumptions")
        if pred.function in (Q.gt, Q.lt):
            lhs, rhs = pred.arguments
            # heuristic: x > y → add real(x) and real(y), if they are symbols
            if lhs.is_Symbol:
                realness_preds.add(Q.real(lhs))
            if rhs.is_Symbol:
                realness_preds.add(Q.real(rhs))

        all_pred.add(pred)
        all_exprs.update(pred.arguments)

    # Ensure only supported predicates are used
    for pred in all_pred:
        if pred.function not in WHITE_LIST and pred.function != Q.ne:
            raise UnhandledInput(f"LRASolver: {pred} is an unhandled predicate")
    for expr in all_exprs:
        if expr.kind == MatrixKind(NumberKind):
            raise UnhandledInput(f"LRASolver: {expr} is of MatrixKind")
        if expr == S.NaN:
            raise UnhandledInput("LRASolver: nan")

    # convert old assumptions into predicates and add them to sat_true and sat_false
    # also check for unhandled predicates
    for assm in extract_pred_from_old_assum(all_exprs, realness_preds):
        n = len(sat_true.encoding)
        if assm not in sat_true.encoding:
            sat_true.encoding[assm] = n+1
        sat_true.data.append([sat_true.encoding[assm]])

        n = len(sat_false.encoding)
        if assm not in sat_false.encoding:
            sat_false.encoding[assm] = n+1
        sat_false.data.append([sat_false.encoding[assm]])

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

def extract_pred_from_old_assum(all_exprs, realness_preds=None):
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
            if expr.is_Symbol:
                if Q.real(expr) not in realness_preds:
                    raise UnhandledInput(f"LRASolver: {expr} must be real")
            else:
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
