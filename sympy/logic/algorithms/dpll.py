"""Implementation of DPLL algorithm

Further improvements: eliminate calls to pl_true, implement branching rules,
efficient unit propagation.

References:
  - http://en.wikipedia.org/wiki/DPLL_algorithm
  - http://bioinformatics.louisville.edu/ouyang/MingOuyangThesis.pdf
"""
from sympy.core import Symbol
from sympy.logic.boolalg import Or, Not, conjuncts, disjuncts, to_cnf
from sympy.logic.inference import pl_true, literal_symbol

def dpll_satisfiable(expr):
    """Check satisfiability of a propositional sentence.
    It returns a model rather than True when it succeeds
    >>> from sympy import symbols
    >>> A, B = symbols('AB')
    >>> dpll_satisfiable(A & ~B)
    {A: True, B: False}
    >>> dpll_satisfiable(A & ~A)
    False
    """
    clauses = conjuncts(to_cnf(expr))
    symbols = list(expr.atoms(Symbol))
    return dpll(clauses, symbols, {})

def dpll(clauses, symbols, model):
    """Compute satisfiability in a partial model"""
    # compute DP kernel
    P, value = find_unit_clause(clauses, model)
    while P:
        model.update({P: value})
        assert P in model
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

### helper methods for DPLL

def unit_propagate(clauses, symbol):
    """Returns an equivalent set of clauses
    http://en.wikipedia.org/wiki/Unit_propagation
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
