"""Inference in propositional logic"""
from __future__ import print_function, division

from sympy.logic.boolalg import And, Or, Not, Implies, Equivalent, \
    BooleanAtom, BooleanFunction, conjuncts, to_cnf, eliminate_implications
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


def satisfiable(expr, return_model=True, algorithm="dpll2"):
    """
    Check satisfiability of a propositional sentence.

    return_model:   Return a model when if expression is satisfiable.
                    Faster when set to False.
    algorithm:      Select algorithm to use for SAT solving - dpll/dpll2


    Examples
    ========

    >>> from sympy.abc import A, B
    >>> from sympy.logic.inference import satisfiable
    >>> satisfiable(A & ~B)
    {A: True, B: False}
    >>> satisfiable(A & ~A)
    False
    >>> satisfiable(A & ~B, return_model=False)
    True
    >>> satisfiable(A & ~A, return_model=False)
    False
    >>> 
    """

    if expr is True:
        return {}
    if expr is False:
        return False

    if return_model is True:
        expr = to_cnf(expr)
        if algorithm == "dpll":
            from sympy.logic.algorithms.dpll import dpll_satisfiable
            return dpll_satisfiable(expr)
        elif algorithm == "dpll2":
            from sympy.logic.algorithms.dpll2 import dpll_satisfiable
            return dpll_satisfiable(expr)
        raise ValueError("'algorithm' must be one of 'dpll', 'dpll2'")

    elif return_model is False:
        return semantic_tableaux(expr)

    else:
        raise ValueError("'return_model' must contain a boolean value")


def pl_true(expr, model={}, deep=False):
    """
    Returns the value of the expression under the given model.

    If the model does not specify the value for every proposition,
    this may return None to indicate 'not obvious'.

    model:  dict containing the symbol, boolean value pair.
    deep:   gives the value of the expression under partial
            interpretations correctly. May still return None.


    Examples
    ========

    >>> from sympy.abc import A, B, C
    >>> from sympy.logic.inference import pl_true
    >>> pl_true(A & B, {A: True, B: True})
    True
    >>> pl_true(A & B, {A: False})
    False
    >>> pl_true(A & B, {A: True})
    >>> pl_true(A >> (B >> A))
    True
    >>> pl_true(A & ~A)
    False
    >>> pl_true(A & B)
    >>> pl_true((C >> A) >> (B >> A), {C: True})
    >>> pl_true((C >> A) >> (B >> A), {C: True}, deep=True)
    True
    >>> pl_true(A & B & (~A | ~B), {A: True})
    >>> pl_true(A & B & (~A | ~B), {A: True}, deep=True)
    False
    >>> pl_true(A | B, {A: False}, deep=True)
    """

    expr = sympify(expr)

    if isinstance(expr, BooleanAtom):
        return bool(expr)

    if expr.is_Symbol:
        return model.get(expr)

    expr = eliminate_implications(expr)

    if deep:
        atoms = set()
        expr = _pl_interpretation(expr, atoms, model)
        if isinstance(expr, BooleanAtom):
            return bool(expr)
        for atom in atoms:
            model[atom] = True

        if pl_true(expr, model, deep=False):
            if satisfiable(Not(expr), return_model=False):
                return None
            else:
                return True
        else:
            if satisfiable(expr, return_model=False):
                return None
            else:
                return False

    if isinstance(expr, Not):
        p = pl_true(expr.args[0], model)
        if p is None:
            return None
        else:
            return not p

    elif isinstance(expr, Or):
        args = set(expr.args)
        gen = (arg for arg in args if is_literal(arg))
        for arg in gen:
            if Not(arg) in args:
                return True
        result = False
        for arg in args:
            p = pl_true(arg, model)
            if p is True:
                return True
            if p is None:
                result = None
        return result

    elif isinstance(expr, And):
        args = set(expr.args)
        gen = (arg for arg in args if is_literal(arg))
        for arg in gen:
            if Not(arg) in args:
                return False
        result = True
        for arg in args:
            p = pl_true(arg, model)
            if p is False:
                return False
            if p is None:
                result = None
        return result

    else:
        raise TypeError("Illegal operator %s" % expr.func)


def _pl_interpretation(expr, atoms, i={}):
    if expr.is_Atom:
        if isinstance(i.get(expr), bool):
            return i[expr]
        else:
            atoms.add(expr)
            return expr
    else:
        args = [_pl_interpretation(arg, atoms, i) for arg in expr.args]
        return expr.func(*args)


def semantic_tableaux(expr):
    """
    Checks satisfiability of a formula using a Semantic Tableaux.

    This method is much faster than a traditional SAT solver, however,
    it returns only True or False. To obtain a model use 'satisfiable'.

    Examples
    ========

    >>> from sympy.abc import A, B
    >>> from sympy.logic.inference import semantic_tableaux
    >>> semantic_tableaux(A | B)
    True
    >>> semantic_tableaux(A & ~A)
    False

    References
    ==========

    .. [1] http://en.wikipedia.org/wiki/Method_of_analytic_tableaux
    """

    expr = sympify(expr)
    if isinstance(expr, BooleanAtom):
        return bool(expr)
    if not isinstance(expr, BooleanFunction):
        raise TypeError("Illegal propositional expression '%s'" % expr)
    expr = eliminate_implications(expr)

    s = [set([expr])]

    while s:
        e = s.pop()

        clause = None
        flag = False

        for c in e:
            if c.is_Atom or c.is_Not :
                if Not(c) in e:
                    flag = True
            else:
                clause = c

        if flag:
            continue
        elif clause is None:
            return True
        else:
            e.remove(clause)

        if isinstance(clause, Or):
            for arg in clause.args[1:]:
                temp = set(e)
                temp.add(arg)
                s.append(temp)
            e.add(clause.args[0])

        elif isinstance(clause, And):
            e.update(clause.args)

        else:
            raise TypeError("Illegal operator %s" % expr.func)

        s.append(e)

    return False


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
