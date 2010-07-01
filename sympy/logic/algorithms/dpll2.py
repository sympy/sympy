"""Implementation of DPLL algorithm

Features:
  - Clause learning
  - Watch literal scheme
  - VSIDS heuristic

References:
  - http://en.wikipedia.org/wiki/DPLL_algorithm
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
    symbols_int_repr = set(range(1, len(symbols) + 1))
    clauses = conjuncts(to_cnf(expr))
    clauses_int_repr = to_int_repr(clauses, symbols)

    solver = SATSolver(clauses_int_repr, symbols_int_repr, {})
    result = solver.find_model()

    if not result:
        return result
    output = {}
    for key in result:
        output.update({symbols[key-1]: result[key]})
    return output

class SATSolver(object):
    """
    Class for representing a SAT solver capable of
     finding a model to a boolean theory in conjunctive
     normal form.
    """

    def __init__(self, clauses, variables, var_settings, heuristic = 'vsids'):
        self.var_settings = var_settings
        self.heuristic = heuristic

        self.initialize_variables(variables)
        self.initialize_clauses(clauses)

        if 'vsids' == heuristic:
            self.calculate_heuristic = self.calculate_vsids
        else:
            raise NotImplementedError

    def initialize_variables(self, variables):
        """Set up the variable data structures needed."""
        pass

    def initialize_clauses(self, clauses):
        """Set up the clause data structures needed."""
        # For now just using sets of literals.
        self.clauses = clauses

    def find_model(self):
        """Main DPLL loop."""
        pass

    #########################
    #      Propagation      #
    #########################
    """
    Propagation methods should attempt to soundly simplify the boolean
      theory, and return True if any simplification occurred and False
      otherwise.
    """
    def simplify(self):
        """Iterate over the various forms of propagation to simplify the theory."""
        changed = True
        while changed:
            changed = False
            changed |= self.unit_prop()
            changed |= self.pure_literal()

    def unit_prop(self):
        """Perform unit propagation on the current theory."""
        return False

    def pure_literal(self):
        """Look for pure literals and assign them when found."""
        return False

    #########################
    #      Heuristics       #
    #########################
    def calculate_vsids(self):
        """
            VSIDS Heuristic
        """
        pass
