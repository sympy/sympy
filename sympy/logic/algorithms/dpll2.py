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
from heapq import heappush, heappop

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

    solver = SATSolver(clauses_int_repr, symbols_int_repr, set([]))
    result = solver.find_model()

    if not result:
        return result

    return dict((symbols[abs(lit) - 1], lit > 0) for lit in solver.var_settings)


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
            self.vsids_init()
            self.heur_calculate = self.vsids_calculate
            self.heur_lit_assigned = self.vsids_lit_assigned
            self.heur_lit_removed = self.vsids_lit_removed
        else:
            raise NotImplementedError

        # Create the base level
        self.levels = [Level(0)]
        self.current_level.varsettings = var_settings

    def initialize_variables(self, variables):
        """Set up the variable data structures needed."""
        self.appears_pos = [[] for i in range(len(variables) + 1)]
        self.appears_neg = [[] for i in range(len(variables) + 1)]
        self.variable_set = [False] * (len(variables) + 1)

    def initialize_clauses(self, clauses):
        """Set up the clause data structures needed."""
        self.clauses = []
        self.clause_sat = []
        for cls in clauses:
            self.clauses.append(list(cls))
            self.clause_sat.append(False)

        self.remaining_clauses = set(range(len(self.clauses)))

        for i in range(len(self.clauses)):
            for lit in self.clauses[i]:
                if lit < 0:
                    self.appears_neg[-lit].append(i)
                else:
                    self.appears_pos[lit].append(i)

    def find_model(self):
        """Main DPLL loop."""

        # We use this variable to keep track of if we should flip a
        #  variable setting in successive rounds
        flip_var = False

        # While the theory still has clauses remaining
        while not self.is_satisfied:
            if flip_var:
                # We have just backtracked and we are trying to opposite literal
                flip_var = False
                lit = self.current_level.decision

            else:
                # Pick a literal to set
                lit = self.heur_calculate()

                # Start the new decision level
                self.levels.append(Level(lit))

            # Assign the literal, updating the clauses it satisfies
            self.assign_literal(lit)

            # Simplify the theory
            self.simplify()
            # Check if we've made the theory unsat
            if self.is_unsatisfied:
                # TODO: Detect and add a learned clause here

                # We unroll all of the decisions until we can flip a literal
                while self.current_level.flipped:
                    self.undo()

                    # If we've unrolled all the way, the theory is unsat
                    if 1 == len(self.levels):
                        return False

                flip_lit = -self.current_level.decision
                self.undo()
                self.levels.append(Level(flip_lit, flipped = True))
                flip_var = True

        return True


    ########################
    #    Helper Methods    #
    ########################
    @property
    def is_satisfied(self):
        return 0 == len(self.remaining_clauses)

    @property
    def is_unsatisfied(self):
        for cls in self.remaining_clauses:
            if self.clause_unsatisfied(cls):
                return True
        return False

    @property
    def current_level(self):
        return self.levels[-1]

    def clause_satisfied(self, cls):
        return self.clause_sat[cls]

    def clause_unsatisfied(self, cls):
        for lit in self.clauses[cls]:
            if -lit not in self.var_settings:
                return False
        return True

    def fetch_clause(self, cls):
        return filter(lambda x: not self.variable_set[abs(x)], self.clauses[cls])

    def assign_literal(self, lit):
        self.var_settings.add(lit)
        self.current_level.var_settings.append(lit)
        self.heur_lit_assigned(lit)
        self.variable_set[abs(lit)] = True

        if lit < 0:
            for cls in self.appears_neg[-lit]:
                if not self.clause_sat[cls]:
                    self.clause_sat[cls] = True
                    self.remaining_clauses.remove(cls)
                    self.current_level.clauses_removed.append(cls)
        else:
            for cls in self.appears_pos[lit]:
                if not self.clause_sat[cls]:
                    self.clause_sat[cls] = True
                    self.remaining_clauses.remove(cls)
                    self.current_level.clauses_removed.append(cls)

    def undo(self):
        """
        Undo the changes of the most recent decision level.
        """
        # Undo the variable settings
        for lit in self.current_level.var_settings:
            self.var_settings.remove(lit)
            self.heur_lit_removed(lit)
            self.variable_set[abs(lit)] = False

        # Undo the clause satisfactions
        for cls in self.current_level.clauses_removed:
            self.clause_sat[cls] = False
            self.remaining_clauses.add(cls)

        # Pop the level off the stack
        self.levels.pop()

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
        for cls in self.remaining_clauses:
            full = self.fetch_clause(cls)
            if 1 == len(full):
                self.assign_literal(full[0])
                return True
        return False

    def pure_literal(self):
        """Look for pure literals and assign them when found."""
        return False

    #########################
    #      Heuristics       #
    #########################
    def vsids_init(self):
        self.lit_heap = []
        for var in range(1, len(self.appears_neg)):
            heappush(self.lit_heap, (-len(self.appears_neg[var]), -var))
            heappush(self.lit_heap, (-len(self.appears_pos[var]), var))

    def vsids_decay(self):
        # We divide everything in the heap by 2 for a decay factor
        #  Note: This doesn't change the heap property
        self.lit_heap = map(lambda x: x / 2.0, self.lit_heap)

    def vsids_calculate(self):
        """
            VSIDS Heuristic
        """
        while self.variable_set[abs(self.lit_heap[0][1])]:
            heappop(self.lit_heap)

        return heappop(self.lit_heap)[1]

    def vsids_lit_assigned(self, lit):
        pass

    def vsids_lit_removed(self, lit):
        var = abs(lit)
        heappush(self.lit_heap, (-len(self.appears_neg), -var))
        heappush(self.lit_heap, (-len(self.appears_pos), var))

class Level(object):
    """
    Represents a single level in the DPLL algorithm, and contains
    enough information for a sound backtracking procedure.
    """

    def __init__(self, decision, flipped = False):
        self.decision = decision
        self.clauses_removed = []
        self.var_settings = []
        self.flipped = flipped
