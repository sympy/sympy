"""Inference in propositional logic"""
from __future__ import print_function, division

from sympy.logic.boolalg import And, Or, Not, Implies, Equivalent, \
    conjuncts, to_cnf, true
from sympy.core.basic import C
from sympy.core.compatibility import ordered
from sympy.core.sympify import sympify


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


def satisfiable(expr, algorithm="dpll2", all_models=False):
    """
    Check satisfiability of a propositional sentence.
    Returns a model when it succeeds.
    Returns {true: true} for trivially true expressions.

    If all_models is True then returns a Model object.
    Calling this object returns the next model or False if no more
    models are available. The object also supports iteration.


    Examples
    ========

    >>> from sympy.abc import A, B
    >>> from sympy.logic.inference import satisfiable
    >>> satisfiable(A & ~B)
    {A: True, B: False}
    >>> satisfiable(A & ~A)
    False
    >>> satisfiable(True)
    {True: True}

    """
    expr = to_cnf(expr)
    if algorithm == "dpll":
        from sympy.logic.algorithms.dpll import dpll_satisfiable
        return dpll_satisfiable(expr)
    elif algorithm == "dpll2":
        from sympy.logic.algorithms.dpll2 import dpll_satisfiable
        return dpll_satisfiable(expr, all_models)
    raise NotImplementedError


def valid(expr):
    """
    Check validity of a propositional sentence.
    A valid propositional sentence is True under every assignment.

    Examples
    ========

    >>> from sympy.abc import A, B
    >>> from sympy.logic.inference import valid
    >>> valid(A | ~A)
    True
    >>> valid(A | B)
    False

    References
    ==========

    .. [1] http://en.wikipedia.org/wiki/Validity

    """
    return not satisfiable(Not(expr))


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


def entails(expr, formula_set={}):
    """
    Check whether the given expr_set entail an expr.
    If formula_set is empty then it returns the validity of expr.

    Examples
    ========

    >>> from sympy.abc import A, B, C
    >>> from sympy.logic.inference import entails
    >>> entails(A, [A >> B, B >> C])
    False
    >>> entails(C, [A >> B, B >> C, A])
    True
    >>> entails(A >> B)
    False
    >>> entails(A >> (B >> A))
    True

    References
    ==========

    .. [1] http://en.wikipedia.org/wiki/Logical_consequence

    """
    formula_set = list(formula_set)
    formula_set.append(Not(expr))
    return not satisfiable(And(*formula_set))


class KB(object):
    """Base class for all knowledge bases"""
    def __init__(self, sentence=None):
        self.clauses_ = set()
        if sentence:
            self.tell(sentence)

    def tell(self, sentence):
        raise NotImplementedError

    def ask(self, query):
        raise NotImplementedError

    def retract(self, sentence):
        raise NotImplementedError

    @property
    def clauses(self):
        return list(ordered(self.clauses_))


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
        [y, Or(x, y)]
        """
        for c in conjuncts(to_cnf(sentence)):
            self.clauses_.add(c)

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
        return entails(query, self.clauses_)

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
            self.clauses_.discard(c)
