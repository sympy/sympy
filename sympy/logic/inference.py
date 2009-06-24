"""Inference in propositional logic"""
from sympy.logic.boolalg import And, Or, Not, Implies, Equivalent, disjuncts, \
    to_cnf
from sympy.core import Symbol, sympify

def find_pure_symbol(symbols, unknown_clauses):
    """Find a symbol and its value if it appears only as a positive literal
    (or only as a negative) in clauses.
    >>> from sympy import symbols
    >>> A, B, C = symbols('ABC')
    >>> find_pure_symbol([A, B, C], [A|~B,~B|~C,C|A])
    (A, True)
    """
    for sym in symbols:
        found_pos, found_neg = False, False
        for c in unknown_clauses:
            if not found_pos and sym in disjuncts(c): found_pos = True
            if not found_neg and Not(sym) in disjuncts(c): found_neg = True
        if found_pos != found_neg: return sym, found_pos
    return None, None

def find_unit_clause(clauses, model):
    """A unit clause has only 1 variable that is not bound in the model.
    >>> from sympy import symbols
    >>> A, B, C = symbols('ABC')
    >>> find_unit_clause([A | B | C, B | ~C, A | ~B], {A:True})
    (B, False)
    """
    for clause in clauses:
        num_not_in_model = 0
        for literal in disjuncts(clause):
            sym = literal_symbol(literal)
            if sym not in model:
                num_not_in_model += 1
                P, value = sym, not (isinstance(literal, Not))
        if num_not_in_model == 1:
            return P, value
    return None, None

def literal_symbol(literal):
    """The symbol in this literal (without the negation).
    >>> from sympy import Symbol
    >>> A = Symbol('A')
    >>> literal_symbol(A)
    A
    >>> literal_symbol(~A)
    A
    """
    if isinstance(literal, Not):
        return literal.args[0]
    else:
        return literal

def satisfiable(expr, algorithm="dpll"):
    """Check satisfiability of a propositional sentence.
    Returns a model when it succeeds

    Examples
    >>> from sympy import symbols
    >>> A, B = symbols('AB')
    >>> satisfiable(A & ~B)
    {A: True, B: False}
    >>> satisfiable(A & ~A)
    False
    """
    expr = to_cnf(expr)
    if algorithm == "dpll":
        from sympy.logic.algorithms.dpll import dpll_satisfiable
        return dpll_satisfiable(expr)
    raise NotImplementedError

def pl_true(expr, model={}):
    """Return True if the propositional logic expression is true in the model,
    and False if it is false. If the model does not specify the value for
    every proposition, this may return None to indicate 'not obvious';
    this may happen even when the expression is tautological.

    The model is implemented as a dict containing the pair symbol, boolean value.

    Examples:
    >>> from sympy import symbols
    >>> A, B = symbols('AB')
    >>> pl_true( A & B, {A: True, B : True})
    True
    """

    if isinstance(expr, bool):
        return expr

    expr = sympify(expr)
    if expr.is_Atom:
        return model.get(expr)

    args = expr.args
    if isinstance(expr, Not):
        p = pl_true(args[0], model)
        if p is None: return None
        else: return not p
    elif isinstance(expr, Or):
        result = False
        for arg in args:
            p = pl_true(arg, model)
            if p == True: return True
            if p == None: result = None
        return result
    elif isinstance(expr, And):
        result = True
        for arg in args:
            p = pl_true(arg, model)
            if p == False: return False
            if p == None: result = None
        return result

    elif isinstance(expr, Implies):
        p, q = args
        return pl_true(Or(Not(p), q), model)

    elif isinstance(expr, Equivalent):
        p, q = args
        pt = pl_true(p, model)
        if pt == None:
            return None
        qt = pl_true(q, model)
        if qt == None:
            return None
        return pt == qt
    else:
        raise ValueError, "Illegal operator in logic expression" + str(expr)
