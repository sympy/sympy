"""Implementation of DPLL algorithm

Features:
  - Clause learning
  - Watch literal scheme
  - VSIDS heuristic

References:
  - http://en.wikipedia.org/wiki/DPLL_algorithm
"""
from sympy.core import Symbol
from sympy import Predicate
from sympy.logic.boolalg import conjuncts, to_cnf, to_int_repr
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
    symbols = list(expr.atoms(Symbol, Predicate))
    symbols_int_repr = set(range(1, len(symbols) + 1))
    clauses = conjuncts(to_cnf(expr))
    clauses_int_repr = to_int_repr(clauses, symbols)

    solver = SATSolver(clauses_int_repr, symbols_int_repr, set())
    result = solver.find_model()

    if not result:
        return result
    # Uncomment to confirm the solution is valid (hitting set for the clauses)
    #else:
        #for cls in clauses_int_repr:
            #assert solver.var_settings.intersection(cls)

    return dict((symbols[abs(lit) - 1], lit > 0) for lit in solver.var_settings)


class SATSolver(object):
    """
    Class for representing a SAT solver capable of
     finding a model to a boolean theory in conjunctive
     normal form.
    """

    def __init__(self, clauses, variables, var_settings, heuristic = 'vsids', \
                 clause_learning = 'none', INTERVAL = 500):
        self.var_settings = var_settings
        self.heuristic = heuristic
        self.is_unsatisfied = False
        self.unit_prop_queue = []
        self.update_functions = []
        self.INTERVAL = INTERVAL

        self.initialize_variables(variables)
        self.initialize_clauses(clauses)

        if 'vsids' == heuristic:
            self.vsids_init()
            self.heur_calculate = self.vsids_calculate
            self.heur_lit_assigned = self.vsids_lit_assigned
            self.heur_lit_unset = self.vsids_lit_unset
            self.heur_clause_added = self.vsids_clause_added

            # Note: Uncomment this if/when clause learning is enabled
            #self.update_functions.append(self.vsids_decay)

        else:
            raise NotImplementedError

        if 'simple' == clause_learning:
            self.add_learned_clause = self.simple_add_learned_clause
            self.compute_conflict = self.simple_compute_conflict
            self.update_functions.append(self.simple_clean_clauses)
        elif 'none' == clause_learning:
            self.add_learned_clause = lambda x: None
            self.compute_conflict = lambda: None
        else:
            raise NotImplementedError

        # Create the base level
        self.levels = [Level(0)]
        self.current_level.varsettings = var_settings

        # Keep stats
        self.num_decisions = 0
        self.num_learned_clauses = 0
        self.original_num_clauses = len(self.clauses)

    def initialize_variables(self, variables):
        """Set up the variable data structures needed."""
        self.sentinels = {}
        self.occurrence_count = {}
        for i in xrange(1, len(variables)+1):
            self.sentinels[i] = set()
            self.sentinels[-i] = set()
            self.occurrence_count[i] = 0
            self.occurrence_count[-i] = 0

        self.variable_set = [False] * (len(variables) + 1)

    def initialize_clauses(self, clauses):
        """Set up the clause data structures needed.

        For each clause, the following changes are made:
        - Unit clauses are queued for propagation right away.
        - Non-unit clauses have their first and last literals set as sentinels.
        - The number of clauses a literal appears in is computed.
        """
        self.clauses = []
        for cls in clauses:
            self.clauses.append(list(cls))

        for i in range(len(self.clauses)):

            # Handle the unit clauses
            if 1 == len(self.clauses[i]):
                self.unit_prop_queue.append(self.clauses[i][0])
                continue

            self.sentinels[self.clauses[i][0]].add(i)
            self.sentinels[self.clauses[i][-1]].add(i)

            for lit in self.clauses[i]:
                self.occurrence_count[lit] += 1


    def find_model(self):
        """Main DPLL loop.

        Variables are chosen successively, and assigned to be either
        True or False. If a solution is not found with this setting,
        the opposite is chosen and the search continues. The solver
        halts when every variable has a setting.
        """

        # We use this variable to keep track of if we should flip a
        #  variable setting in successive rounds
        flip_var = False

        # Check if unit prop says the theory is unsat right off the bat
        self.simplify()
        if self.is_unsatisfied:
            return False

        # While the theory still has clauses remaining
        while True:
            # Perform cleanup / fixup at regular intervals
            if self.num_decisions % self.INTERVAL == 0:
                for func in self.update_functions:
                    func()

            if flip_var:
                # We have just backtracked and we are trying to opposite literal
                flip_var = False
                lit = self.current_level.decision

            else:
                # Pick a literal to set
                lit = self.heur_calculate()
                self.num_decisions += 1

                # Stopping condition for a satisfying theory
                if 0 == lit:
                    return True

                # Start the new decision level
                self.levels.append(Level(lit))

            # Assign the literal, updating the clauses it satisfies
            self.assign_literal(lit)

            # Simplify the theory
            self.simplify()

            # Check if we've made the theory unsat
            if self.is_unsatisfied:

                self.is_unsatisfied = False

                # We unroll all of the decisions until we can flip a literal
                while self.current_level.flipped:
                    self.undo()

                    # If we've unrolled all the way, the theory is unsat
                    if 1 == len(self.levels):
                        return False

                # Detect and add a learned clause
                self.add_learned_clause(self.compute_conflict())

                # Try the opposite setting of the most recent decision
                flip_lit = -self.current_level.decision
                self.undo()
                self.levels.append(Level(flip_lit, flipped = True))
                flip_var = True


    ########################
    #    Helper Methods    #
    ########################
    @property
    def current_level(self):
        """The current decision level data structure"""
        return self.levels[-1]

    def clause_sat(self, cls):
        """Check if a clause is satisfied by the current variable setting."""
        for lit in self.clauses[cls]:
            if lit in self.var_settings:
                return True
        return False

    def is_sentinel(self, lit, cls):
        """Check if a literal is a sentinel of a given clause."""
        return cls in self.sentinels[lit]

    def assign_literal(self, lit):
        """Make a literal assignment.

        The literal assignment must be recorded as part of the current
        decision level. Additionally, if the literal is marked as a
        sentinel of any clause, then a new sentinel must be chosen. If
        this is not possible, then unit propagation is triggered and
        another literal is added to the queue to be set in the future.
        """
        self.var_settings.add(lit)
        self.current_level.var_settings.add(lit)
        self.variable_set[abs(lit)] = True
        self.heur_lit_assigned(lit)

        sentinel_list = list(self.sentinels[-lit])

        for cls in sentinel_list:
            if not self.clause_sat(cls):
                other_sentinel = None
                for newlit in self.clauses[cls]:
                    if newlit != -lit:
                        if self.is_sentinel(newlit, cls):
                            other_sentinel = newlit
                        elif not self.variable_set[abs(newlit)]:
                            self.sentinels[-lit].remove(cls)
                            self.sentinels[newlit].add(cls)
                            other_sentinel = None
                            break

                # Check if no sentinel update exists
                if other_sentinel:
                    self.unit_prop_queue.append(other_sentinel)

    def undo(self):
        """
        Undo the changes of the most recent decision level.
        """
        # Undo the variable settings
        for lit in self.current_level.var_settings:
            self.var_settings.remove(lit)
            self.heur_lit_unset(lit)
            self.variable_set[abs(lit)] = False

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
        result = len(self.unit_prop_queue) > 0

        while self.unit_prop_queue:
            next_lit = self.unit_prop_queue.pop()
            if -next_lit in self.var_settings:
                self.is_unsatisfied = True
                self.unit_prop_queue = []
                return False
            else:
                self.assign_literal(next_lit)

        return result

    def pure_literal(self):
        """Look for pure literals and assign them when found."""
        return False

    #########################
    #      Heuristics       #
    #########################
    def vsids_init(self):
        """Initialize the data structures needed for the VSIDS heuristic."""
        self.lit_heap = []
        self.lit_scores = {}
        for var in range(1, len(self.variable_set)):
            self.lit_scores[var] = -float(self.occurrence_count[var])
            self.lit_scores[-var] = -float(self.occurrence_count[-var])
            heappush(self.lit_heap, (self.lit_scores[var], var))
            heappush(self.lit_heap, (self.lit_scores[-var], -var))

    def vsids_decay(self):
        """Decay the VSIDS scores for every literal."""
        # We divide every literal score by 2 for a decay factor
        #  Note: This doesn't change the heap property
        for lit in self.lit_scores.keys():
            self.lit_scores[lit] /= 2.0

    def vsids_calculate(self):
        """
            VSIDS Heuristic Calculation
        """
        if len(self.lit_heap) == 0:
            return 0

        # Clean out the front of the heap as long the variables are set
        while self.variable_set[abs(self.lit_heap[0][1])]:
            heappop(self.lit_heap)
            if len(self.lit_heap) == 0:
                return 0

        return heappop(self.lit_heap)[1]

    def vsids_lit_assigned(self, lit):
        """Handle the assignment of a literal for the VSIDS heuristic."""
        pass

    def vsids_lit_unset(self, lit):
        """Handle the unsetting of a literal for the VSIDS heuristic."""
        var = abs(lit)
        heappush(self.lit_heap, (self.lit_scores[var], var))
        heappush(self.lit_heap, (self.lit_scores[-var], -var))

    def vsids_clause_added(self, cls):
        """Handle the addition of a new clause for the VSIDS heuristic."""
        self.num_learned_clauses += 1
        for lit in cls:
            self.lit_scores[lit] += 1


    ########################
    #   Clause Learning    #
    ########################

    def simple_add_learned_clause(self, cls):
        """Add a new clause to the theory."""

        cls_num = len(self.clauses)
        self.clauses.append(cls)

        for lit in cls:
            self.occurrence_count[lit] += 1

        self.sentinels[cls[0]].add(cls_num)
        self.sentinels[cls[-1]].add(cls_num)

        self.heur_clause_added(cls)

    def simple_compute_conflict(self):
        """ Build a clause representing the fact that at least one decision made
        so far is wrong.
        """
        return [-(level.decision) for level in self.levels[1:]]

    def simple_clean_clauses(self):
        """Clean up learned clauses."""
        pass

class Level(object):
    """
    Represents a single level in the DPLL algorithm, and contains
    enough information for a sound backtracking procedure.
    """

    def __init__(self, decision, flipped = False):
        self.decision = decision
        self.var_settings = set()
        self.flipped = flipped
