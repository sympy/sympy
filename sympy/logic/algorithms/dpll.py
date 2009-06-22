from sympy.core import Symbol
from sympy.logic.boolalg import conjuncts, to_cnf
from sympy.logic.inference import find_pure_symbol, find_unit_clause, pl_true

def dpll_satisfiable(expr):
    """Check satisfiability of a propositional sentence.
    It returns a model rather than True when it succeeds
    >>> from sympy import symbols
    >>> A, B = symbols('AB')
    >>> dpll_satisfiable(A & ~B)
    {A: True, B: False}
    >>> dpll_satisfiable(A & ~A)
    False

    References: Implemented as described in http://aima.cs.berkeley.edu/
    """
    clauses = conjuncts(to_cnf(expr))
    symbols = list(expr.atoms(Symbol))
    return dpll(clauses, symbols, {})

def dpll(clauses, symbols, model):
    """See if the clauses are true in a partial model.
    http://en.wikipedia.org/wiki/DPLL_algorithm
    """
    unknown_clauses = [] ## clauses with an unknown truth value
    for c in clauses:
        val =  pl_true(c, model)
        if val == False:
            return False
        if val != True:
            unknown_clauses.append(c)
    if not unknown_clauses:
        return model
    P, value = find_pure_symbol(symbols, unknown_clauses)
    if P:
        model_1 = model.copy()
        model_1.update({P: value})
        syms = [x for x in symbols if x != P]
        return dpll(clauses, syms, model_1)
    P, value = find_unit_clause(unknown_clauses, model)
    if P:
        model_1 = model.copy()
        model_1.update({P: value})
        syms = [x for x in symbols if x != P]
        return dpll(clauses, syms, model_1)
    P = symbols.pop()
    model_1, model_2 = model.copy(), model.copy()
    model_1.update({P: True})
    model_2.update({P: False})
    return (dpll(clauses, symbols, model_1) or
            dpll(clauses, symbols, model_2))
