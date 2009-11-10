"""Implementation of DPLL algorithm

Further improvements: eliminate calls to pl_true, implement branching rules,
efficient unit propagation.

References:
  - http://en.wikipedia.org/wiki/DPLL_algorithm
  - http://bioinformatics.louisville.edu/ouyang/MingOuyangThesis.pdf
"""
from sympy.core import Symbol
from sympy.logic.boolalg import Or, Not, conjuncts, disjuncts, to_cnf, \
    to_int_repr
from sympy.logic.inference import pl_true, literal_symbol

def dpll_satisfiable(expr):
    """
    Check satisfiability of a propositional sentence.
    It returns a model rather than True when it succeeds
    >>> from sympy import symbols
    >>> from sympy.abc import A, B
    >>> from sympy.logic.algorithms.dpll import dpll_satisfiable
    >>> dpll_satisfiable(A & ~B)
    {A: True, B: False}
    >>> dpll_satisfiable(A & ~A)
    False

    """
    symbols = list(expr.atoms(Symbol))
    symbols_int_repr = range(1, len(symbols) + 1)
    clauses = conjuncts(to_cnf(expr))
    clauses_int_repr = to_int_repr(clauses, symbols)
    result = dpll_int_repr(clauses_int_repr, symbols_int_repr, {})
    if not result: return result
    output = {}
    for key in result:
        output.update({symbols[key-1]: result[key]})
    return output

def dpll(clauses, symbols, model):
    """
    Compute satisfiability in a partial model.
    Clauses is an array of conjuncts.

    >>> from sympy.abc import A, B, C
    >>> from sympy.logic.algorithms.dpll import dpll
    >>> dpll([A, B, C], [A, B], {C: False})
    False

    """
    # compute DP kernel
    P, value = find_unit_clause(clauses, model)
    while P:
        model.update({P: value})
        symbols.remove(P)
        if not value: P = ~P
        clauses = unit_propagate(clauses, P)
        P, value = find_unit_clause(clauses, model)
    P, value = find_pure_symbol(symbols, clauses)
    while P:
        model.update({P: value})
        symbols.remove(P)
        if not value: P = ~P
        clauses = unit_propagate(clauses, P)
        P, value = find_pure_symbol(symbols, clauses)
    # end DP kernel
    unknown_clauses = []
    for c in clauses:
        val =  pl_true(c, model)
        if val == False:
            return False
        if val != True:
            unknown_clauses.append(c)
    if not unknown_clauses:
        return model
    if not clauses: return model
    P = symbols.pop()
    model_copy = model.copy()
    model.update({P: True})
    model_copy.update({P: False})
    symbols_copy = symbols[:]
    return (dpll(clauses, symbols, model) or
            dpll(clauses, symbols_copy, model_copy))

def dpll_int_repr(clauses, symbols, model):
    """
    Compute satisfiability in a partial model.
    Arguments are expected to be in integer representation

    >>> from sympy.logic.algorithms.dpll import dpll_int_repr
    >>> dpll_int_repr([[1], [2], [3]], [1, 2], {3: False})
    False

    """
    # compute DP kernel
    P, value = find_unit_clause_int_repr(clauses, model)
    while P:
        model.update({P: value})
        symbols.remove(P)
        if not value: P = -P
        clauses = unit_propagate_int_repr(clauses, P)
        P, value = find_unit_clause_int_repr(clauses, model)
    P, value = find_pure_symbol_int_repr(symbols, clauses)
    while P:
        model.update({P: value})
        symbols.remove(P)
        if not value: P = -P
        clauses = unit_propagate_int_repr(clauses, P)
        P, value = find_pure_symbol_int_repr(symbols, clauses)
    # end DP kernel
    unknown_clauses = []
    for c in clauses:
        val =  pl_true_int_repr(c, model)
        if val == False:
            return False
        if val != True:
            unknown_clauses.append(c)
    if not unknown_clauses:
        return model
    if not clauses: return model
    P = symbols.pop()
    model_copy = model.copy()
    model.update({P: True})
    model_copy.update({P: False})
    symbols_copy = symbols[:]
    return (dpll_int_repr(clauses, symbols, model) or
            dpll_int_repr(clauses, symbols_copy, model_copy))

### helper methods for DPLL

def pl_true_int_repr(clause, model={}):
    """
    Lightweight version of pl_true.
    Argument clause represents the args of an Or clause. This is used
    inside dpll_int_repr, it is not meant to be used directly.

    >>> from sympy.logic.algorithms.dpll import pl_true_int_repr
    >>> pl_true_int_repr([1, 2], {1: False})
    >>> pl_true_int_repr([1, 2], {1: False, 2: False})
    False

    """
    if len(clause) == 1:
        c = clause[0]
        if c in model: return model.get(c)
        if -c in model: return not model.get(-c)
        return None
    result = False
    for l in clause:
        p = pl_true_int_repr([l], model)
        if p: return True
        if p is None: result = None
    return result

def unit_propagate(clauses, symbol):
    """
    Returns an equivalent set of clauses
    If a set of clauses contains the unit clause l, the other clauses are
    simplified by the application of the two following rules:

      1. every clause containing l is removed
      2. in every clause that contains ~l this literal is deleted

    Arguments are expected to be in CNF.

    >>> from sympy import symbols
    >>> from sympy.abc import A, B, C
    >>> from sympy.logic.algorithms.dpll import unit_propagate
    >>> unit_propagate([A | B, C | ~B, B], B)
    [C, B]

    """
    output = []
    for c in clauses:
        if not isinstance(c, Or):
            output.append(c)
            continue
        for arg in c.args:
            if arg == ~symbol:
                output.append(Or(*[x for x in c.args if x != ~symbol]))
                break
            if arg == symbol:
                break
        else:
            output.append(c)
    return output

def unit_propagate_int_repr(clauses, symbol):
    """
    Same as above, but arguments are expected to be in integer
    representation

    >>> from sympy.logic.algorithms.dpll import unit_propagate_int_repr
    >>> unit_propagate_int_repr([[1, 2], [3, -2], [2]], 2)
    [[3], [2]]

    """
    output = []
    for c in clauses:
        if len(c) == 1:
            output.append(c)
            continue
        for l in c:
            if l == symbol: break
            elif l == -symbol:
                output.append([x for x in c if x != l])
                break
        else: output.append(c)
    return output

def find_pure_symbol(symbols, unknown_clauses):
    """
    Find a symbol and its value if it appears only as a positive literal
    (or only as a negative) in clauses.

    >>> from sympy import symbols
    >>> from sympy.abc import A, B, C
    >>> from sympy.logic.algorithms.dpll import find_pure_symbol
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

def find_pure_symbol_int_repr(symbols, unknown_clauses):
    """
    Same as find_pure_symbol, but arguments are expected
    to be in integer representation

    >>> from sympy.logic.algorithms.dpll import find_pure_symbol_int_repr
    >>> find_pure_symbol_int_repr([1,2,3], [[1, -2], [-2, -3], [3, 1]])
    (1, True)

    """
    for sym in symbols:
        found_pos, found_neg = False, False
        for c in unknown_clauses:
            if not found_pos and sym in c: found_pos = True
            if not found_neg and -sym in c: found_neg = True
        if found_pos ^ found_neg: return sym, found_pos
    return None, None

def find_unit_clause(clauses, model):
    """
    A unit clause has only 1 variable that is not bound in the model.

    >>> from sympy import symbols
    >>> from sympy.abc import A, B, C
    >>> from sympy.logic.algorithms.dpll import find_unit_clause
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

def find_unit_clause_int_repr(clauses, model):
    """
    Same as find_unit_clause, but arguments are expected to be in
    integer representation.

    >>> from sympy.logic.algorithms.dpll import find_unit_clause_int_repr
    >>> find_unit_clause_int_repr([[1, 2, 3], [2, -3], [1, -2]], {1: True})
    (2, False)

    """
    for c in clauses:
        num_not_in_model = False
        for literal in c:
            sym = abs(literal)
            if sym not in model:
                if num_not_in_model:
                    num_not_in_model = False
                    break
                num_not_in_model = True
                P, value = sym, literal > 0
        if num_not_in_model:
            return P, value
    return None, None
