"""
sympy.assumptions.euf
=====================

Interface for the Equality with Uninterpreted Functions (EUF) solver for SymPy's
assumptions system. This module provides the `euf_ask` function, which integrates
the EUF solver with the assumptions query framework, handling only unary predicates
and equality/disequality.

Classes and Functions
---------------------
- EUFUnhandledInput: Exception for unsupported or contradictory input.
- is_unary_predicate: Utility to check if a predicate is unary.
- euf_ask: Main interface for querying the EUF solver with assumptions.

Doctests
--------
>>> from sympy.assumptions.ask import Q
>>> from sympy.assumptions.euf import euf_ask
>>> from sympy import symbols
>>> x, y = symbols('x y')
>>> euf_ask(Q.prime(x), Q.prime(y) & Q.eq(x, y))
True
>>> euf_ask(Q.prime(x), Q.composite(y) & Q.eq(x, y))
False
>>> euf_ask(Q.integer(x), Q.even(y) & Q.eq(x, y))
True
>>> euf_ask(Q.even(x), Q.prime(y) & Q.eq(x, y))
>>> euf_ask(Q.eq(x, y), Q.eq(x, y))
True
>>> euf_ask(Q.ne(x, y), Q.eq(x, y))
False
>>> euf_ask(Q.eq(x, y), Q.ne(x, y))
False
"""

from sympy.assumptions.ask import Q
from sympy.logic.algorithms.euf_theory import EUFUnhandledInput, EUFSolver
from sympy.assumptions.cnf import CNF, EncodedCNF
from sympy.assumptions.assume import AppliedPredicate

WHITE_LIST = {Q.eq, Q.ne}

def is_unary_predicate(pred):
    """
    Determine if a predicate is unary.

    Parameters
    ----------
    pred : AppliedPredicate

    Returns
    -------
    bool
        True if the predicate takes one argument, False otherwise.

    Examples
    --------
    >>> from sympy.assumptions.ask import Q
    >>> from sympy.abc import x, y
    >>> from sympy.assumptions.euf import is_unary_predicate
    >>> is_unary_predicate(Q.prime(x))
    True
    >>> is_unary_predicate(Q.eq(x, y))
    False
    """
    return len(pred.arguments) == 1

def euf_ask(proposition, assumptions=True, context=None):
    """
    Evaluate a proposition with assumptions using the EUF solver.

    Only unary predicates and equality/disequality are supported.
    Returns True/False if determined, or None if undetermined.

    Parameters
    ----------
    proposition : AppliedPredicate
        The query to evaluate.
    assumptions : Boolean, optional
        The assumptions under which to evaluate the proposition.
    context : CNF, optional
        Additional context assumptions.

    Returns
    -------
    bool or None

    Raises
    ------
    EUFUnhandledInput
        If unsupported predicates are encountered.

    Examples
    --------
    >>> from sympy.assumptions.ask import Q
    >>> from sympy.assumptions.euf import euf_ask
    >>> from sympy import symbols
    >>> x, y = symbols('x y')
    >>> euf_ask(Q.prime(x), Q.prime(y) & Q.eq(x, y))
    True
    >>> euf_ask(Q.prime(x), Q.composite(y) & Q.eq(x, y))
    False
    >>> euf_ask(Q.integer(x), Q.even(y) & Q.eq(x, y))
    True
    """
    cnf = CNF.from_prop(assumptions)
    if context:
        cnf.extend(context)
    enc_cnf = EncodedCNF()
    enc_cnf.from_cnf(cnf)

    # Check for unsupported predicates
    for pred in enc_cnf.encoding:
        if isinstance(pred, AppliedPredicate):
            if is_unary_predicate(pred):
                pass
            elif pred.function in WHITE_LIST:
                pass
            else:
                raise EUFUnhandledInput(f"EUFSolver: {pred.function} is not supported")

    solver = EUFSolver.from_encoded_cnf(enc_cnf)

    if not solver.check_consistency():
        return False

    # Handle query
    if isinstance(proposition, AppliedPredicate):
        if proposition.function in [Q.eq, Q.ne] and len(proposition.arguments) == 2:
            a, b = proposition.arguments
            root_a = solver.find(solver.get_var(a))
            root_b = solver.find(solver.get_var(b))
            if proposition.function == Q.eq:
                if root_a == root_b:
                    return True
                return False
            else:
                if root_a == root_b:
                    return False
                return True
        elif len(proposition.arguments) == 1:
            pred, arg = proposition.function, proposition.arguments[0]
            res = solver.query(pred, arg)
            if res is not None:
                return res
    return None
