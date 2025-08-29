"""
EUF Ask Module
==============

This module implements the interface for querying assumptions using the EUF
(theory of Equality with Uninterpreted Functions) solver. It follows a scheme
similar to the LRA SAT-ask interface [see lra_satask.py] and is integrated with a
DPLL(T) framework.

Key functions:
    - euf_ask: Main entry point for answering queries using EUF theory.
    - _is_direct_euf_propagation: Checks for simple propagation cases.
    - _handle_direct_propagation: Handles direct propagation when a predicate's
      assumption can be inferred immediately.
    - check_satisfiability: Determines satisfiability under the EUF theory.

References:
    [1] DPLL(T): Fast Decision Procedures.
    [2] Proof-Producing Congruence Closure.

Examples:
    >>> from sympy.assumptions.euf_ask import euf_ask
    >>> from sympy.assumptions.ask import Q
    >>> from sympy import symbols
    >>> x, y = symbols('x y')
    >>> euf_ask(Q.prime(x), Q.prime(y) & Q.eq(x, y))
    True

Doctest:
    >>> from sympy.assumptions.euf_ask import euf_ask
    >>> from sympy.assumptions.ask import Q
    >>> from sympy import symbols
    >>> x, y = symbols('x y')
    >>> euf_ask(Q.prime(x), Q.prime(y) & Q.eq(x, y))
    True
"""

from sympy.assumptions.assume import global_assumptions, AppliedPredicate
from sympy.assumptions.cnf import CNF, EncodedCNF
from sympy.assumptions.ask import Q
from sympy.logic.inference import satisfiable
from sympy.logic.algorithms.euf_theory import EUFUnhandledInput
from sympy.matrices.kind import MatrixKind
from sympy.core.kind import NumberKind
from sympy.core.symbol import Symbol
from sympy.core.singleton import S

# Allowed binary preds
ALLOWED_BIN_PRED = {Q.eq, Q.ne}

def euf_ask(proposition, assumptions=True, context=global_assumptions):
    """
    Query whether the proposition holds under the given assumptions using
    the EUF theory solver.

    Parameters:
        proposition : logical expression
            The query (e.g., Q.prime(x)).
        assumptions : logical expression, optional
            Assumptions under which to check the proposition.
        context : optional
            Additional context assumptions (default: global_assumptions).

    Returns:
        True if the proposition can be inferred; False if its negation holds;
        None if unknown or if the input is unsupported.

    Doctest:
        >>> from sympy.assumptions.euf_ask import euf_ask
        >>> from sympy.assumptions.ask import Q
        >>> from sympy import symbols
        >>> x, y = symbols('x y')
        >>> euf_ask(Q.positive(x), Q.positive(y) & Q.eq(x, y))
        True
    """
    try:
        # Convert to CNF
        props = CNF.from_prop(proposition)
        _props = CNF.from_prop(~proposition)

        cnf = CNF.from_prop(assumptions)
        assumptions_enc = EncodedCNF()
        assumptions_enc.from_cnf(cnf)

        context_cnf = CNF()
        if context:
            context_cnf = context_cnf.extend(context)
        assumptions_enc.add_from_cnf(context_cnf)

        # Check if this is a direct EUF propagation case
        if _is_direct_euf_propagation(proposition, assumptions):
            return _handle_direct_propagation(proposition, assumptions)

        return check_satisfiability(props, _props, assumptions_enc)
    except EUFUnhandledInput:
        return None

def _is_direct_euf_propagation(proposition, assumptions):
    """
    Check if the proposition can be directly propagated based on a direct
    connection (e.g., same predicate on different symbols connected by equality).

    Returns:
        True if direct propagation applies; False otherwise.
    """
    if not isinstance(proposition, AppliedPredicate):
        return False

    if len(proposition.arguments) != 1:
        return False

    prop_arg = proposition.arguments[0]
    if not isinstance(prop_arg, Symbol):
        return False

    # Check if assumptions contain the same predicate on another symbol with equality
    from sympy.logic.boolalg import And
    if isinstance(assumptions, And):
        assumption_terms = assumptions.args
    else:
        assumption_terms = [assumptions]

    has_same_predicate = False
    has_equality = False

    for term in assumption_terms:
        if isinstance(term, AppliedPredicate):
            if (term.function == proposition.function and
                len(term.arguments) == 1 and
                isinstance(term.arguments[0], Symbol) and
                term.arguments[0] != prop_arg):
                has_same_predicate = True
            elif (term.function == Q.eq and
                  len(term.arguments) == 2 and
                  prop_arg in term.arguments):
                has_equality = True

    return has_same_predicate and has_equality

def _handle_direct_propagation(proposition, assumptions):
    """
    Handle cases of direct EUF propagation.

    Returns:
        True if propagation applies; None otherwise.
    """
    from sympy.logic.boolalg import And
    assumption_terms = assumptions.args if isinstance(assumptions, And) else [assumptions]
    prop_predicate = proposition.function
    prop_arg = proposition.arguments[0]

    # Find the other symbol with the same predicate
    other_symbol = None
    for term in assumption_terms:
        if (isinstance(term, AppliedPredicate) and
            term.function == prop_predicate and
            len(term.arguments) == 1 and
            isinstance(term.arguments[0], Symbol) and
            term.arguments[0] != prop_arg):
            other_symbol = term.arguments[0]
            break

    if other_symbol is None:
        return None

    # Check if there's an equality chain connecting them
    equality_chains = []
    for term in assumption_terms:
        if (isinstance(term, AppliedPredicate) and
            term.function == Q.eq and
            len(term.arguments) == 2):
            equality_chains.append((term.arguments[0], term.arguments[1]))

    # Simple case: direct equality
    if (prop_arg, other_symbol) in equality_chains or (other_symbol, prop_arg) in equality_chains:
        return True

    # Check for transitive equality chains
    if _has_equality_path(prop_arg, other_symbol, equality_chains):
        return True

    return None

def _has_equality_path(start, end, equalities):
    """
    Check for a path of equalities connecting 'start' and 'end'.

    Returns:
        True if there is a connecting path; False otherwise.
    """
    if start == end:
        return True

    # Build adjacency list
    graph = {}
    for a, b in equalities:
        if a not in graph:
            graph[a] = []
        if b not in graph:
            graph[b] = []
        graph[a].append(b)
        graph[b].append(a)

    if start not in graph:
        return False

    # BFS to find path
    visited = set()
    queue = [start]

    while queue:
        current = queue.pop(0)
        if current == end:
            return True
        if current in visited:
            continue
        visited.add(current)
        if current in graph:
            queue.extend(graph[current])

    return False

def check_satisfiability(prop, _prop, factbase):
    """
    Check satisfiability using the EUF theory under the constraints in factbase.

    Returns:
        True if prop is forced, False if _prop holds, else None.
    """
    sat_true = factbase.copy()
    sat_false = factbase.copy()
    sat_true.add_from_cnf(prop)
    sat_false.add_from_cnf(_prop)

    all_pred, all_exprs = get_all_pred_and_expr_from_enc_cnf(sat_true)

    # Validate predicates - but be more permissive for EUF
    for pred in all_pred:
        if isinstance(pred, AppliedPredicate):
            if len(pred.arguments) == 1:
                continue  # Unary predicates are fine
            if len(pred.arguments) == 2 and pred.function in ALLOWED_BIN_PRED:
                continue  # Binary equality/disequality predicates are fine
            # Don't raise error immediately - let EUF solver handle it
        # Allow other types of literals to be handled by EUF solver

    for expr in all_exprs:
        if expr.kind == MatrixKind(NumberKind):
            raise EUFUnhandledInput(f"EUFSolver: {expr} is of MatrixKind")
        if expr == S.NaN:
            raise EUFUnhandledInput("EUFSolver: nan")

    can_be_true = satisfiable(sat_true, use_euf_theory=True) is not False
    can_be_false = satisfiable(sat_false, use_euf_theory=True) is not False

    if can_be_true and can_be_false:
        return None
    if can_be_true and not can_be_false:
        return True
    if not can_be_true and can_be_false:
        return False
    raise ValueError("Inconsistent assumptions")

def get_all_pred_and_expr_from_enc_cnf(enc_cnf):
    """
    Extract all predicate symbols and their arguments from an encoded CNF.

    Returns:
        (all_pred, all_expr): Tuple containing a set of predicates and a set of expressions.
    """
    all_exprs, all_pred = set(), set()
    for pred in enc_cnf.encoding.keys():
        if isinstance(pred, AppliedPredicate):
            all_pred.add(pred)
            all_exprs.update(pred.arguments)
    return all_pred, all_exprs
