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


_SIGN_TO_BINREL = {
    Q.positive: Q.gt,
    Q.negative: Q.lt,
    Q.zero: Q.eq,
    Q.nonzero: Q.ne,
    Q.nonpositive: Q.le,
    Q.nonnegative: Q.ge,
    Q.extended_positive: Q.gt,
    Q.extended_negative: Q.lt,
    Q.extended_nonpositive: Q.le,
    Q.extended_nonzero: Q.ne,
}

REAL_IMPLIED = {Q.real, Q.extended_real, Q.complex, Q.finite, Q.commutative,
                Q.hermitian}


def lra_satask(proposition, assumptions=True):
    """Evaluate a proposition using SAT and linear real arithmetic."""
    return check_satisfiability(CNF.from_prop(proposition),
                               CNF.from_prop(~proposition),
                               CNF.from_prop(assumptions))


def check_satisfiability(prop, _prop, factbase, known_real=frozenset()):
    """Answer *prop* with the LRA theory solver, taking CNF inputs.

    *known_real* holds expressions established to be real independently of
    their old assumptions, such as by root level inference in ``satask``.
    """
    predicates = prop.all_predicates() | _prop.all_predicates() | factbase.all_predicates()
    replacements = {pred: _pred_to_binrel(pred, known_real) for pred in predicates}
    expressions = {arg for pred in predicates if isinstance(pred, AppliedPredicate)
                   for arg in pred.arguments}
    for expr in expressions:
        _validate_expression(expr, known_real)
    factbase = factbase.copy()
    for pred in extract_pred_from_old_assum(expressions):
        factbase.add(pred)
        replacements[pred] = _pred_to_binrel(pred, known_real)

    encoded = EncodedCNF()
    encoded.from_cnf(_preprocess(factbase, replacements))
    engine = ReasoningEngine(encoded, use_lra_theory=True)
    query = engine.create_query(_preprocess(prop, replacements),
                                _preprocess(_prop, replacements))
    return engine.ask_query(query)


def _preprocess(cnf, replacements):
    """Rewrite CNF literals without changing the clause structure."""
    return CNF({frozenset(new_lit for lit in clause
                          for new_lit in _rewrite_literal(lit, replacements[lit.lit]))
                for clause in cnf.clauses})


def _rewrite_literal(literal, pred):
    """Return the disjunction representing a converted literal for LRA."""
    negated = literal.is_Not
    if pred in (True, False):
        return (Literal(S.true if bool(pred) != negated else S.false),)
    if isinstance(pred, AppliedPredicate) and pred.function in (Q.eq, Q.ne):
        if (pred.function == Q.ne) != negated:
            return (Literal(Q.gt(*pred.arguments)), Literal(Q.lt(*pred.arguments)))
        return (Literal(Q.eq(*pred.arguments)),)
    return (Literal(pred, negated),)


def _is_real(expr, known_real):
    return expr in known_real or getattr(expr, "is_real", None) is True


def _pred_to_binrel(pred, known_real):
    """Validate a predicate and convert it to a relation or Boolean constant."""
    if not isinstance(pred, AppliedPredicate):
        return pred
    function = pred.function
    if function in REAL_IMPLIED and all(_is_real(arg, known_real) for arg in pred.arguments):
        return S.true
    if function in _SIGN_TO_BINREL:
        return _SIGN_TO_BINREL[function](pred.arguments[0], 0)
    if function in (Q.negative_infinite, Q.positive_infinite):
        return S.false
    if function in ALLOWED_PRED or function == Q.ne:
        return pred
    raise UnhandledInput(f"LRASolver: {pred} is an unhandled predicate")


def _validate_expression(expr, known_real):
    """Reject domains and expressions that real arithmetic cannot represent."""
    if getattr(expr, "kind", None) == MatrixKind(NumberKind):
        raise UnhandledInput(f"LRASolver: {expr} is of MatrixKind")
    if expr == S.NaN:
        raise UnhandledInput("LRASolver: nan")
    if not getattr(expr, "free_symbols", None):
        return
    if not _is_real(expr, known_real):
        raise UnhandledInput(f"LRASolver: {expr} must be real")
    if isinstance(expr, Mul) and not all(_is_real(arg, known_real) for arg in expr.args):
        raise UnhandledInput(f"LRASolver: {expr} must be real")

    if expr.is_integer == True and expr.is_zero != True:
        raise UnhandledInput(f"LRASolver: {expr} is an integer")
    if expr.is_integer == False:
        raise UnhandledInput(f"LRASolver: {expr} can't be an integer")
    if expr.is_rational == False:
        raise UnhandledInput(f"LRASolver: {expr} is irational")


def extract_pred_from_old_assum(all_exprs):
    """Extract sign facts from expressions already validated for LRA.

    Examples
    ========
    >>> from sympy.assumptions.lra_satask import extract_pred_from_old_assum
    >>> from sympy import symbols
    >>> x, y = symbols("x y", positive=True)
    >>> extract_pred_from_old_assum([x, y, 2])
    [Q.positive(x), Q.positive(y)]
    """
    ret = []
    for expr in all_exprs:
        if not getattr(expr, "free_symbols", None):
            continue
        for name in ("zero", "positive", "negative", "nonzero", "nonpositive", "nonnegative"):
            if getattr(expr, "is_" + name):
                ret.append(getattr(Q, name)(expr))
                break

    return ret
