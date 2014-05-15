"""Inference in propositional logic"""
from __future__ import print_function, division

from sympy.logic.boolalg import And, Or, Not, Implies, Equivalent, \
    conjuncts, to_cnf
from sympy.core.basic import C
from sympy.core.sympify import sympify


def is_literal(expr):
    """
    Returns True if expr is a literal, else False.

    Examples
    ========

    >>> from sympy import Symbol, Or
    >>> from sympy.abc import A, B
    >>> from sympy.logic.inference import is_literal
    >>> is_literal(A)
    True
    >>> is_literal(~A)
    True
    >>> is_literal(Or(A, B))
    False

    """

    try:
        literal_symbol(expr)
        return True
    except (ValueError):
        return False


def literal_symbol(literal):
    """
    The symbol in this literal (without the negation).

    Examples
    ========

    >>> from sympy import Symbol
    >>> from sympy.abc import A
    >>> from sympy.logic.inference import literal_symbol
    >>> literal_symbol(A)
    A
    >>> literal_symbol(~A)
    A

    """

    if literal is True or literal is False:
        return literal
    try:
        if literal.is_Symbol:
            return literal
        if literal.is_Not:
            return literal_symbol(literal.args[0])
        else:
            raise ValueError
    except (AttributeError, ValueError):
        raise ValueError("Argument must be a boolean literal.")


def satisfiable(expr, algorithm="dpll2"):
    """
    Check satisfiability of a propositional sentence.
    Returns a model when it succeeds

    Examples:

    >>> from sympy.abc import A, B
    >>> from sympy.logic.inference import satisfiable
    >>> satisfiable(A & ~B)
    {A: True, B: False}
    >>> satisfiable(A & ~A)
    False

    """
    if expr is True:
        return {}
    if expr is False:
        return False
    expr = to_cnf(expr)
    if algorithm == "dpll":
        from sympy.logic.algorithms.dpll import dpll_satisfiable
        return dpll_satisfiable(expr)
    elif algorithm == "dpll2":
        from sympy.logic.algorithms.dpll2 import dpll_satisfiable
        return dpll_satisfiable(expr)
    raise NotImplementedError


def pl_true(expr, model={}):
    """
    Return True if the propositional logic expression is true in the model,
    and False if it is false. If the model does not specify the value for
    every proposition, this may return None to indicate 'not obvious';
    this may happen even when the expression is tautological.

    The model is implemented as a dict containing the pair symbol, boolean value.

    Examples
    ========

    >>> from sympy.abc import A, B
    >>> from sympy.logic.inference import pl_true
    >>> pl_true( A & B, {A: True, B : True})
    True

    """

    if isinstance(expr, bool):
        return expr

    expr = sympify(expr)

    if expr.is_Symbol:
        return model.get(expr)

    args = expr.args
    func = expr.func

    if func is Not:
        p = pl_true(args[0], model)
        if p is None:
            return None
        else:
            return not p
    elif func is Or:
        result = False
        for arg in args:
            p = pl_true(arg, model)
            if p is True:
                return True
            if p is None:
                result = None
        return result
    elif func is And:
        result = True
        for arg in args:
            p = pl_true(arg, model)
            if p is False:
                return False
            if p is None:
                result = None
        return result

    elif func is Implies:
        p, q = args
        return pl_true(Or(Not(p), q), model)

    elif func is Equivalent:
        p, q = args
        pt = pl_true(p, model)
        if pt is None:
            return None
        qt = pl_true(q, model)
        if qt is None:
            return None
        return pt == qt
    else:
        raise ValueError("Illegal operator in logic expression" + str(expr))


class KB(object):
    """Base class for all knowledge bases"""
    def __init__(self, sentence=None):
        self.clauses = []
        if sentence:
            self.tell(sentence)

    def tell(self, sentence):
        raise NotImplementedError

    def ask(self, query):
        raise NotImplementedError

    def retract(self, sentence):
        raise NotImplementedError


class PropKB(KB):
    """A KB for Propositional Logic.  Inefficient, with no indexing."""

    def tell(self, sentence):
        """Add the sentence's clauses to the KB

        Examples
        ========

        >>> from sympy.logic.inference import PropKB
        >>> from sympy.abc import x, y
        >>> l = PropKB()
        >>> l.clauses
        []

        >>> l.tell(x | y)
        >>> l.clauses
        [Or(x, y)]

        >>> l.tell(y)
        >>> l.clauses
        [Or(x, y), y]
        """
        for c in conjuncts(to_cnf(sentence)):
            if not c in self.clauses:
                self.clauses.append(c)

    def ask(self, query):
        """Checks if the query is true given the set of clauses.

        Examples
        ========

        >>> from sympy.logic.inference import PropKB
        >>> from sympy.abc import x, y
        >>> l = PropKB()
        >>> l.tell(x & ~y)
        >>> l.ask(x)
        True
        >>> l.ask(y)
        False
        """
        if len(self.clauses) == 0:
            return False
        from sympy.logic.algorithms.dpll import dpll
        query_conjuncts = self.clauses[:]
        query_conjuncts.extend(conjuncts(to_cnf(query)))
        s = set()
        for q in query_conjuncts:
            s = s.union(q.atoms(C.Symbol))
        return bool(dpll(query_conjuncts, list(s), {}))

    def retract(self, sentence):
        """Remove the sentence's clauses from the KB

        Examples
        ========

        >>> from sympy.logic.inference import PropKB
        >>> from sympy.abc import x, y
        >>> l = PropKB()
        >>> l.clauses
        []

        >>> l.tell(x | y)
        >>> l.clauses
        [Or(x, y)]

        >>> l.retract(x | y)
        >>> l.clauses
        []
        """
        for c in conjuncts(to_cnf(sentence)):
            if c in self.clauses:
                self.clauses.remove(c)
