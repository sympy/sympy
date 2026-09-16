from __future__ import annotations
from sympy.assumptions.cnf import CNF, EncodedCNF, Literal
from sympy.assumptions.ask import Q
from sympy.assumptions.reasoning_engine import ReasoningEngine
from sympy.logic.algorithms.lra_theory import UnhandledInput, ALLOWED_PRED
from sympy.matrices.kind import MatrixKind
from sympy.core.kind import NumberKind
from sympy.assumptions.assume import AppliedPredicate
from sympy.core.mul import Mul
from sympy.core.singleton import S


def lra_satask(proposition, assumptions=True):
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

    return check_satisfiability(props, _props, assumptions)

# Some predicates such as Q.prime can't be handled by lra_satask.
# For example, (x > 0) & (x < 1) & Q.prime(x) is unsat but lra_satask would think it was sat.
# WHITE_LIST is a list of predicates that can always be handled.
WHITE_LIST = ALLOWED_PRED.keys() | {Q.positive, Q.negative, Q.zero, Q.nonzero, Q.nonpositive, Q.nonnegative,
                                    Q.extended_positive, Q.extended_negative, Q.extended_nonpositive,
                                    Q.extended_negative, Q.extended_nonzero, Q.negative_infinite,
                                    Q.positive_infinite}

# Predicates that every real expression satisfies. Nothing reaches the theory
# solver unless all of its expressions are real, so these say nothing that is
# not already known and are handled by replacing them with True.
REAL_IMPLIED = {Q.real, Q.extended_real, Q.complex, Q.finite, Q.commutative,
                Q.hermitian}


def check_satisfiability(prop, _prop, factbase, known_real=frozenset()):
    """Answer *prop* with the LRA theory solver.

    *known_real* holds the expressions that are known to be real by something
    other than their old assumptions, such as the root level inference that
    ``satask`` does before handing over to this solver.
    """
    all_facts = factbase.copy()
    all_facts.add_from_cnf(prop)
    all_facts.add_from_cnf(_prop)
    all_pred, all_exprs = get_all_pred_and_expr_from_enc_cnf(all_facts)

    trivially_true = set()
    for pred in all_pred:
        if pred.function in REAL_IMPLIED:
            # What makes these true is that every argument is real, which is
            # either an old assumption or something *known_real* settled.
            if all(arg in known_real or getattr(arg, "is_real", None) is True
                   for arg in pred.arguments):
                trivially_true.add(pred)
                continue
            raise UnhandledInput(f"LRASolver: {pred} is an unhandled predicate")
        if pred.function not in WHITE_LIST and pred.function != Q.ne:
            raise UnhandledInput(f"LRASolver: {pred} is an unhandled predicate")
    for expr in all_exprs:
        if expr.kind == MatrixKind(NumberKind):
            raise UnhandledInput(f"LRASolver: {expr} is of MatrixKind")
        if expr == S.NaN:
            raise UnhandledInput("LRASolver: nan")

    factbase = factbase.copy()
    for assm in extract_pred_from_old_assum(all_exprs, known_real):
        factbase.add_prop(assm)

    preprocessed = EncodedCNF()
    preprocessed.from_cnf(_preprocess(factbase, trivially_true))
    engine = ReasoningEngine(preprocessed, use_lra_theory=True)
    sides = []
    for side in (prop, _prop):
        encoded = EncodedCNF()
        encoded.from_cnf(side)
        sides.append(_preprocess(encoded, trivially_true))
    query = engine.create_query(*sides)
    return engine.ask_query(query)


def _preprocess(enc_cnf, true_preds=frozenset()):
    """Return CNF containing only LRA predicates and Boolean constants.

    Replace *true_preds* with True, convert disequalities into disjunctions
    of strict inequalities, and convert negated disequalities into equalities.
    """
    rev_encoding = {value: key for key, value in enc_cnf.encoding.items()}
    clauses = set()
    for clause in enc_cnf.data:
        new_clause = set()
        for lit in clause:
            if lit == 0:
                new_clause.add(Literal(S.false))
                continue
            pred = rev_encoding[abs(lit)]
            pred = S.true if pred in true_preds else _pred_to_binrel(pred)
            negated = lit < 0
            if pred in (True, False):
                new_clause.add(Literal(S.true if bool(pred) != negated else S.false))
                continue
            if isinstance(pred, AppliedPredicate):
                if negated and pred.function == Q.eq:
                    pred = Q.ne(*pred.arguments)
                    negated = False
                if pred.function == Q.ne:
                    if negated:
                        new_clause.add(Literal(Q.eq(*pred.arguments)))
                    else:
                        new_clause.update((Literal(Q.gt(*pred.arguments)),
                                           Literal(Q.lt(*pred.arguments))))
                    continue
            new_clause.add(Literal(pred, negated))
        clauses.add(frozenset(new_clause))
    return CNF(clauses)


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

def extract_pred_from_old_assum(all_exprs, known_real=frozenset()):
    """
    Returns a list of relevant new assumption predicate
    based on any old assumptions.

    Raises an UnhandledInput exception if any of the assumptions are
    unhandled. An expression listed in *known_real* is taken to be real
    even if its old assumptions do not say so.

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

        if expr.is_real is not True and expr not in known_real:
            raise UnhandledInput(f"LRASolver: {expr} must be real")
        # test for I times imaginary variable; such expressions are considered real
        if isinstance(expr, Mul) and any(arg.is_real is not True
                                         and arg not in known_real
                                         for arg in expr.args):
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
