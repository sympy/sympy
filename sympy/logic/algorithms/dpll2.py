"""Implementation of DPLL algorithm

Features:
  - Clause learning
  - Watch literal scheme
  - VSIDS heuristic
  - IPASIR style interface for incremental solving

References:
  - https://en.wikipedia.org/wiki/DPLL_algorithm
  - https://satcompetition.github.io/2020/track_incremental.html
"""
from __future__ import annotations

from collections import defaultdict
from copy import deepcopy
from enum import Enum
from heapq import heappush, heappop

from sympy.core.sorting import ordered
from sympy.assumptions.cnf import EncodedCNF

from sympy.logic.algorithms.lra_theory import LRASolver


class IpasirStatus(Enum):
    # Status codes used by the IPASIR style interface of SATSolver. The values are
    # the ones mandated by the IPASIR standard for ``ipasir_solve``.
    UNKNOWN = 0
    SATISFIABLE = 10
    UNSATISFIABLE = 20


def dpll_satisfiable(expr, all_models=False, use_lra_theory=False):
    """
    Check satisfiability of a propositional sentence.
    It returns a model rather than True when it succeeds.
    Returns a generator of all models if all_models is True.

    Examples
    ========

    >>> from sympy.abc import A, B
    >>> from sympy.logic.algorithms.dpll2 import dpll_satisfiable
    >>> dpll_satisfiable(A & ~B)
    {A: True, B: False}
    >>> dpll_satisfiable(A & ~A)
    False

    """
    if not isinstance(expr, EncodedCNF):
        exprs = EncodedCNF()
        exprs.add_prop(expr)
        expr = exprs

    # Return UNSAT when False (encoded as 0) is present in the CNF
    if {0} in expr.data:
        if all_models:
            return (f for f in [False])
        return False

    if use_lra_theory:
        lra, immediate_conflicts = LRASolver.from_encoded_cnf(expr)
    else:
        lra = None
        immediate_conflicts = []
    solver = SATSolver(expr.data + immediate_conflicts, expr.variables, set(), expr.symbols, lra_theory=lra)
    models = solver._find_model()

    if all_models:
        return _all_models(models)

    try:
        return next(models)
    except StopIteration:
        return False

    # Uncomment to confirm the solution is valid (hitting set for the clauses)
    #else:
        #for cls in clauses_int_repr:
            #assert solver.var_settings.intersection(cls)


def _all_models(models):
    satisfiable = False
    try:
        while True:
            yield next(models)
            satisfiable = True
    except StopIteration:
        if not satisfiable:
            yield False


class SATSolver:
    """
    Class for representing a SAT solver capable of
     finding a model to a boolean theory in conjunctive
     normal form.
    """

    def __init__(self, clauses, variables, var_settings, symbols=None,
                heuristic='vsids', clause_learning='none', INTERVAL=500,
                 lra_theory = None):

        self.var_settings = var_settings
        self.heuristic = heuristic
        self.is_unsatisfied = False
        self._unit_prop_queue = []
        self.update_functions = []
        self.INTERVAL = INTERVAL

        if symbols is None:
            self.symbols = list(ordered(variables))
        else:
            self.symbols = symbols

        self._initialize_variables(variables)
        self._initialize_clauses(clauses)

        if 'vsids' == heuristic:
            self._vsids_init()
            self.heur_calculate = self._vsids_calculate
            self.heur_lit_assigned = self._vsids_lit_assigned
            self.heur_lit_unset = self._vsids_lit_unset
            self.heur_clause_added = self._vsids_clause_added

            # Note: Uncomment this if/when clause learning is enabled
            #self.update_functions.append(self._vsids_decay)

        else:
            raise NotImplementedError

        if 'simple' == clause_learning:
            self.add_learned_clause = self._simple_add_learned_clause
            self.compute_conflict = self._simple_compute_conflict
            self.update_functions.append(self._simple_clean_clauses)
        elif 'none' == clause_learning:
            self.add_learned_clause = lambda x: None
            self.compute_conflict = lambda: None
        else:
            raise NotImplementedError

        self.lra = lra_theory

        # Create the base level
        self.levels = []
        self._create_level(0)
        self._current_level.var_settings = set(var_settings)
        if self.lra and self._current_level.var_settings:
            raise NotImplementedError("A non-empty var_settings is not "
                                      "supported when using the LRA theory.")

        # Keep stats
        self.num_decisions = 0
        self.num_learned_clauses = 0
        self.original_num_clauses = len(self.clauses)

        self.lra = lra_theory

        # State of the IPASIR style interface
        self._status = IpasirStatus.UNKNOWN
        self._models = None
        self._clause_buffer = []
        self._assumptions = []

    def _initialize_variables(self, variables):
        """Set up the variable data structures needed."""
        self.sentinels = defaultdict(set)
        self.occurrence_count = defaultdict(int)
        self.variable_set = [False] * (len(variables) + 1)

    def _initialize_clauses(self, clauses):
        """Set up the clause data structures needed.

        For each clause, the following changes are made:
        - Unit clauses are queued for propagation right away.
        - Non-unit clauses have their first and last literals set as sentinels.
        - The number of clauses a literal appears in is computed.
        """
        self.clauses = [list(clause) for clause in clauses]

        for i, clause in enumerate(self.clauses):

            # Handle the unit clauses
            if 1 == len(clause):
                self._unit_prop_queue.append(clause[0])
                continue

            self.sentinels[clause[0]].add(i)
            self.sentinels[clause[-1]].add(i)

            for lit in clause:
                self.occurrence_count[lit] += 1

    def _find_model(self):
        """
        Main DPLL loop. Returns a generator of models.

        Variables are chosen successively, and assigned to be either
        True or False. If a solution is not found with this setting,
        the opposite is chosen and the search continues. The solver
        halts when every variable has a setting.

        Examples
        ========

        >>> from sympy.logic.algorithms.dpll2 import SATSolver
        >>> l = SATSolver([{2, -3}, {1}, {3, -3}, {2, -2},
        ... {3, -2}], {1, 2, 3}, set())
        >>> list(l._find_model())
        [{1: True, 2: False, 3: False}, {1: True, 2: True, 3: True}]

        >>> from sympy.abc import A, B, C
        >>> l = SATSolver([{2, -3}, {1}, {3, -3}, {2, -2},
        ... {3, -2}], {1, 2, 3}, set(), [A, B, C])
        >>> list(l._find_model())
        [{A: True, B: False, C: False}, {A: True, B: True, C: True}]

        """

        # We use this variable to keep track of if we should flip a
        #  variable setting in successive rounds
        flip_var = False

        # Check if unit prop says the theory is unsat right off the bat
        self._simplify()
        if self.is_unsatisfied:
            return

        # The assumptions are decisions the search cannot take back. Failing
        # one rules out models, not the clauses themselves.
        for assumed_lit in self._assumptions:
            if assumed_lit in self.var_settings:
                continue

            if -assumed_lit not in self.var_settings:
                self.levels.append(Level(assumed_lit))
                self._assign_literal(assumed_lit)
                self._simplify()
                if not self.is_unsatisfied:
                    continue
                self.is_unsatisfied = False

            while len(self.levels) > 1:
                self._undo()
            return

        # Undoing an assumption would answer a different question.
        assumption_level = len(self.levels)

        # While the theory still has clauses remaining
        while True:
            # Perform cleanup / fixup at regular intervals
            if self.num_decisions % self.INTERVAL == 0:
                for func in self.update_functions:
                    func()

            if flip_var:
                # We have just backtracked and we are trying to opposite literal
                flip_var = False
                lit = self._current_level.decision

            else:
                # Pick a literal to set
                lit = self.heur_calculate()
                self.num_decisions += 1

                # Stopping condition for a satisfying theory
                if 0 == lit:
                    res = None
                    if self.lra:
                        res = self.lra.check()
                    if res is None or res[0]:
                        yield {self.symbols[abs(lit) - 1]:
                                    lit > 0 for lit in self.var_settings}
                    else:
                        self._simple_add_learned_clause(res[1])

                        # Backtrack until reaching a level with one of the conflict causing literals.
                        inconsistent_literals = [-lit for lit in res[1]]
                        while True:
                            if len(self.levels) == assumption_level:
                                # If theory-inconsistent literals were set right off the bat
                                # at level 0, the formula is unsat.
                                return

                            if any(inconsistent_lit in self._current_level.var_settings for inconsistent_lit in inconsistent_literals):
                                break
                            self._undo()

                    # To find the next model after yield, or after adding a conflict clause,
                    # simulate a conflict and backtrack to the most recent unflipped decision.
                    while self._current_level.flipped:
                        self._undo()
                    if len(self.levels) == assumption_level:
                        return
                    flip_lit = -self._current_level.decision
                    self._undo()
                    self._create_level(flip_lit, flipped=True)
                    flip_var = True
                    continue

                # Start the new decision level
                self._create_level(lit)

            # Assign the literal, updating the clauses it satisfies
            conflict = self._assign_literal(lit)
            if conflict is not None:
                self.is_unsatisfied = True
                self._simple_add_learned_clause(conflict)
                self._unit_prop_queue = []
            else:
                # simplify the theory
                self._simplify()

            # Check if we've made the theory unsat
            if self.is_unsatisfied:

                self.is_unsatisfied = False

                # We unroll all of the decisions until we can flip a literal
                while self._current_level.flipped:
                    self._undo()

                    # If we've unrolled all the way, the theory is unsat
                    if assumption_level == len(self.levels):
                        return

                # The literal to flip would be an assumption, not ours to flip.
                if assumption_level == len(self.levels):
                    return

                # Detect and add a learned clause
                self.add_learned_clause(self.compute_conflict())

                # Try the opposite setting of the most recent decision
                flip_lit = -self._current_level.decision
                self._undo()
                self._create_level(flip_lit, flipped=True)
                flip_var = True

    ###############################
    #    IPASIR Style Interface   #
    ###############################

    # The following code is not implemented fully, there are
    # a lot of things that can be added to the Interface. This is
    # the most minimalistic version of it. Add new required features as TODO.

    """
    A subset of the IPASIR standard for incremental SAT solving, using the
    names that CaDiCaL gives them in its C++ API. Only the parts needed to
    inspect the root level before searching, to keep solving after new
    clauses have been added and to solve under assumptions are implemented so
    far; ``failed`` is not supported yet.

    # https://github.com/arminbiere/cadical/blob/master/src/cadical.hpp
    """
    def propagate(self):
        """Propagate the unit clauses at the root level, deciding nothing.

        Returns ``UNSATISFIABLE`` on a conflict, ``SATISFIABLE`` if it leaves
        no variable unassigned, and ``UNKNOWN`` otherwise.

        A conflict the LRA theory finds is reported too, but never a model.

        TODO: IPASIR propagates at any decision level, while this is limited
        to the root.

        Examples
        ========

        >>> from sympy.logic.algorithms.dpll2 import SATSolver, IpasirStatus
        >>> l = SATSolver([{1}, {-1, 2}], {1, 2}, set())
        >>> l.propagate() == IpasirStatus.SATISFIABLE
        True
        >>> l.fixed(2)
        1

        """
        if len(self.levels) > 1:
            raise ValueError("propagate() can only be used at the root level.")

        self._simplify()
        if self.is_unsatisfied:
            self._status = IpasirStatus.UNSATISFIABLE
        elif self.lra is None and all(self.variable_set[1:]):
            # Nothing is left to decide on, so the assignments are a model.
            self._status = IpasirStatus.SATISFIABLE

        return self._status

    def fixed(self, lit):
        """Return 1 if *lit* is implied by the clauses, -1 if ``-lit`` is
        implied, and 0 if neither is known yet.

        TODO: IPASIR answers this at any decision level, which needs the root
        level assignments to be kept apart from those a decision implies.

        """
        if len(self.levels) > 1:
            raise ValueError("fixed() is only defined at the root level.")

        if lit in self.var_settings:
            return 1
        if -lit in self.var_settings:
            return -1
        return 0

    def solve(self):
        """Search for a model, reusing the work ``propagate()`` already did,
        and return ``SATISFIABLE`` or ``UNSATISFIABLE``.

        """
        if self._models is not None:
            raise ValueError("solve() can only be called again once new "
                "clauses have been added with add() or clause(), as "
                "restarting the same search is not implemented yet.")

        while len(self.levels) > 1:
            self._undo()

        assumed = bool(self._assumptions)

        self._models = self._find_model()
        if next(self._models, None) is None:
            self._status = IpasirStatus.UNSATISFIABLE
        else:
            self._status = IpasirStatus.SATISFIABLE

        # Asking again without them is a different question, not a restart.
        if assumed:
            self._assumptions = []
            self._models = None

        return self._status

    def val(self, lit):
        """Return *lit* if it is true in the model found by ``solve()``,
        ``-lit`` if it is false there, and 0 if the model does not assign it.

        """
        if self._status != IpasirStatus.SATISFIABLE:
            raise ValueError("val() is only defined once solve() has returned "
                "SATISFIABLE.")

        if lit in self.var_settings:
            return lit
        if -lit in self.var_settings:
            return -lit
        return 0

    def assume(self, lit):
        """Constrain the next call to ``solve()`` with *lit*, which it drops
        again afterwards, letting one solver answer several questions.

        TODO: ``failed()``, which reports the assumptions a search could not
        satisfy, is not implemented yet.

        Examples
        ========

        >>> from sympy.logic.algorithms.dpll2 import SATSolver, IpasirStatus
        >>> l = SATSolver([{1, 2}, {-1, -2}], {1, 2}, set())
        >>> l.assume(1)
        >>> l.solve() == IpasirStatus.SATISFIABLE
        True
        >>> l.val(2)
        -2

        """
        if lit == 0 or abs(lit) >= len(self.variable_set):
            raise ValueError(f"{lit} is not a literal of one of the variables "
                "the solver was created with.")

        while len(self.levels) > 1:
            self._undo()

        self._models = None
        if self._status == IpasirStatus.SATISFIABLE:
            self._status = IpasirStatus.UNKNOWN

        self._assumptions.append(lit)

    def add(self, lit):
        """Add *lit* to the clause being built, or add that clause to the
        solver when *lit* is 0.

        The search restarts from the root level, so ``solve()`` may be called
        again, and an unsatisfiable solver stays unsatisfiable.

        TODO: only literals of the variables the solver was created with can
        be added, as there is no way to introduce a new variable yet.

        Examples
        ========

        >>> from sympy.logic.algorithms.dpll2 import SATSolver, IpasirStatus
        >>> l = SATSolver([{1, 2}], {1, 2}, set())
        >>> l.solve() == IpasirStatus.SATISFIABLE
        True
        >>> l.add(-1)
        >>> l.add(0)
        >>> l.solve() == IpasirStatus.SATISFIABLE
        True
        >>> l.val(1)
        -1

        """
        if lit != 0:
            if abs(lit) >= len(self.variable_set):
                raise ValueError("%s is not a literal of one of the variables "
                    "the solver was created with." % lit)

            self._clause_buffer.append(lit)
            return

        clause_to_add = self._clause_buffer
        self._clause_buffer = []

        # The decisions were made without this clause, so the search restarts.
        while len(self.levels) > 1:
            self._undo()

        self._models = None
        if self._status == IpasirStatus.SATISFIABLE:
            self._status = IpasirStatus.UNKNOWN

        clause_num = len(self.clauses)
        self.clauses.append(clause_to_add)

        for clause_lit in clause_to_add:
            self.occurrence_count[clause_lit] += 1

        # Only a literal that is not already false can be a sentinel. With
        # fewer than two of those the clause is satisfied, unit, or false, and
        # none of those needs to be watched.
        unassigned = [clause_lit for clause_lit in clause_to_add
                      if not self.variable_set[abs(clause_lit)]]

        if len(unassigned) > 1:
            self.sentinels[unassigned[0]].add(clause_num)
            self.sentinels[unassigned[-1]].add(clause_num)
        elif not any(clause_lit in self.var_settings
                     for clause_lit in clause_to_add):
            if unassigned:
                self._unit_prop_queue.append(unassigned[0])
            else:
                self.is_unsatisfied = True
                self._status = IpasirStatus.UNSATISFIABLE

    def clause(self, *lits):
        """Add the clause made up of *lits*, given one by one or as a single
        iterable, which covers the ``clause`` overloads of CaDiCaL.

        Without any literal it adds the empty clause, which is false.

        """
        if len(lits) == 1 and not isinstance(lits[0], int):
            lits = lits[0]

        for lit in lits:
            self.add(lit)

        # Zero is never a literal, which is why IPASIR uses it to mark the
        # end of a clause rather than passing a length around.
        self.add(0)

    def copy(self):
        """Return an independent solver with the same clauses and state, so
        that adding clauses to it or searching with it changes nothing here.

        Unlike ``copy`` in CaDiCaL this returns a new solver and copies the
        state of the search along with the formula.

        Examples
        ========

        >>> from sympy.logic.algorithms.dpll2 import SATSolver, IpasirStatus
        >>> l = SATSolver([{1}, {-1, 2}], {1, 2}, set())
        >>> temporary = l.copy()
        >>> temporary.clause(-2)
        >>> temporary.solve() == IpasirStatus.UNSATISFIABLE
        True
        >>> l.solve() == IpasirStatus.SATISFIABLE
        True

        """
        # A generator cannot be copied, and the symbols are only ever read.
        models, symbols = self._models, self.symbols
        self._models, self.symbols = None, None

        try:
            other = deepcopy(self)
        finally:
            self._models, self.symbols = models, symbols

        other.symbols = symbols

        return other

    ########################
    #    Helper Methods    #
    ########################
    @property
    def _current_level(self):
        """The current decision level data structure

        Examples
        ========

        >>> from sympy.logic.algorithms.dpll2 import SATSolver
        >>> l = SATSolver([{1}, {2}], {1, 2}, set())
        >>> next(l._find_model())
        {1: True, 2: True}
        >>> l._current_level.decision
        0
        >>> l._current_level.flipped
        False
        >>> l._current_level.var_settings
        {1, 2}

        """
        return self.levels[-1]

    def _clause_sat(self, cls):
        """Check if a clause is satisfied by the current variable setting.

        Examples
        ========

        >>> from sympy.logic.algorithms.dpll2 import SATSolver
        >>> l = SATSolver([{1}, {-1}], {1}, set())
        >>> try:
        ...     next(l._find_model())
        ... except StopIteration:
        ...     pass
        >>> l._clause_sat(0)
        False
        >>> l._clause_sat(1)
        True

        """
        for lit in self.clauses[cls]:
            if lit in self.var_settings:
                return True
        return False

    def _is_sentinel(self, lit, cls):
        """Check if a literal is a sentinel of a given clause.

        Examples
        ========

        >>> from sympy.logic.algorithms.dpll2 import SATSolver
        >>> l = SATSolver([{2, -3}, {1}, {3, -3}, {2, -2},
        ... {3, -2}], {1, 2, 3}, set())
        >>> next(l._find_model())
        {1: True, 2: False, 3: False}
        >>> l._is_sentinel(2, 3)
        True
        >>> l._is_sentinel(-3, 1)
        False

        """
        return cls in self.sentinels[lit]

    def _assign_literal(self, lit):
        """Make a literal assignment.

        The literal assignment must be recorded as part of the current
        decision level. Additionally, if the literal is marked as a
        sentinel of any clause, then a new sentinel must be chosen. If
        this is not possible, then unit propagation is triggered and
        another literal is added to the queue to be set in the future.

        Examples
        ========

        >>> from sympy.logic.algorithms.dpll2 import SATSolver
        >>> l = SATSolver([{2, -3}, {1}, {3, -3}, {2, -2},
        ... {3, -2}], {1, 2, 3}, set())
        >>> next(l._find_model())
        {1: True, 2: False, 3: False}
        >>> l.var_settings
        {-3, -2, 1}

        >>> l = SATSolver([{2, -3}, {1}, {3, -3}, {2, -2},
        ... {3, -2}], {1, 2, 3}, set())
        >>> l._assign_literal(-1)
        >>> try:
        ...     next(l._find_model())
        ... except StopIteration:
        ...     pass
        >>> l.var_settings
        {-1}

        """
        self.var_settings.add(lit)
        self._current_level.var_settings.add(lit)
        self.variable_set[abs(lit)] = True
        self.heur_lit_assigned(lit)

        conflict = None
        if self.lra:
            res = self.lra.assert_lit(lit)
            if res and res[0] is False:
                conflict = res[1]

        sentinel_list = list(self.sentinels[-lit])

        for cls in sentinel_list:
            if not self._clause_sat(cls):
                other_sentinel = None
                for newlit in self.clauses[cls]:
                    if newlit != -lit:
                        if self._is_sentinel(newlit, cls):
                            other_sentinel = newlit
                        elif not self.variable_set[abs(newlit)]:
                            self.sentinels[-lit].remove(cls)
                            self.sentinels[newlit].add(cls)
                            other_sentinel = None
                            break

                # Check if no sentinel update exists
                if other_sentinel:
                    self._unit_prop_queue.append(other_sentinel)

        return conflict

    def _create_level(self, lit, flipped=False):
        """
        Start a new decision level for ``lit``.

        If a theory solver is being used it is told to start a new level too,
        so that the bounds asserted while this level is current can all be
        undone together when `_undo` pops the level.
        """
        if self.lra:
            self.lra.push_level()
        self.levels.append(Level(lit, flipped=flipped))

    def _undo(self):
        """
        _undo the changes of the most recent decision level.

        Examples
        ========

        >>> from sympy.logic.algorithms.dpll2 import SATSolver
        >>> l = SATSolver([{2, -3}, {1}, {3, -3}, {2, -2},
        ... {3, -2}], {1, 2, 3}, set())
        >>> next(l._find_model())
        {1: True, 2: False, 3: False}
        >>> level = l._current_level
        >>> level.decision, level.var_settings, level.flipped
        (-3, {-3, -2}, False)
        >>> l._undo()
        >>> level = l._current_level
        >>> level.decision, level.var_settings, level.flipped
        (0, {1}, False)

        """
        if len(self.levels) == 1:
            raise IndexError("Cannot pop base decision level 0.")

        # Undo the variable settings
        for lit in self._current_level.var_settings:
            self.var_settings.remove(lit)
            self.heur_lit_unset(lit)
            self.variable_set[abs(lit)] = False

        if self.lra:
            self.lra.pop_level()

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
    def _simplify(self):
        """Iterate over the various forms of propagation to simplify the theory.

        Examples
        ========

        >>> from sympy.logic.algorithms.dpll2 import SATSolver
        >>> l = SATSolver([{2, -3}, {1}, {3, -3}, {2, -2},
        ... {3, -2}], {1, 2, 3}, set())
        >>> l.variable_set
        [False, False, False, False]
        >>> l.sentinels
        {-3: {0, 2}, -2: {3, 4}, 2: {0, 3}, 3: {2, 4}}

        >>> l._simplify()

        >>> l.variable_set
        [False, True, False, False]
        >>> l.sentinels
        {-3: {0, 2}, -2: {3, 4}, -1: set(), 2: {0, 3},
        ...3: {2, 4}}

        """
        changed = True
        while changed:
            changed = False
            changed |= self._unit_prop()
            changed |= self._pure_literal()

    def _unit_prop(self):
        """Perform unit propagation on the current theory."""
        result = len(self._unit_prop_queue) > 0
        while self._unit_prop_queue:
            next_lit = self._unit_prop_queue.pop()
            if -next_lit in self.var_settings:
                self.is_unsatisfied = True
            else:
                conflict = self._assign_literal(next_lit)
                if conflict is not None:
                    self.is_unsatisfied = True
                    self._simple_add_learned_clause(conflict)

            if self.is_unsatisfied:
                self._unit_prop_queue = []
                return False

        return result

    def _pure_literal(self):
        """Look for pure literals and assign them when found."""
        return False

    #########################
    #      Heuristics       #
    #########################
    def _vsids_init(self):
        """Initialize the data structures needed for the VSIDS heuristic."""
        self.lit_heap = []
        self.lit_scores = {}

        for var in range(1, len(self.variable_set)):
            self.lit_scores[var] = float(-self.occurrence_count[var])
            self.lit_scores[-var] = float(-self.occurrence_count[-var])
            heappush(self.lit_heap, (self.lit_scores[var], var))
            heappush(self.lit_heap, (self.lit_scores[-var], -var))

    def _vsids_decay(self):
        """Decay the VSIDS scores for every literal.

        Examples
        ========

        >>> from sympy.logic.algorithms.dpll2 import SATSolver
        >>> l = SATSolver([{2, -3}, {1}, {3, -3}, {2, -2},
        ... {3, -2}], {1, 2, 3}, set())

        >>> l.lit_scores
        {-3: -2.0, -2: -2.0, -1: 0.0, 1: 0.0, 2: -2.0, 3: -2.0}

        >>> l._vsids_decay()

        >>> l.lit_scores
        {-3: -1.0, -2: -1.0, -1: 0.0, 1: 0.0, 2: -1.0, 3: -1.0}

        """
        # We divide every literal score by 2 for a decay factor
        #  Note: This doesn't change the heap property
        for lit in self.lit_scores.keys():
            self.lit_scores[lit] /= 2.0

    def _vsids_calculate(self):
        """
            VSIDS Heuristic Calculation

        Examples
        ========

        >>> from sympy.logic.algorithms.dpll2 import SATSolver
        >>> l = SATSolver([{2, -3}, {1}, {3, -3}, {2, -2},
        ... {3, -2}], {1, 2, 3}, set())

        >>> l.lit_heap
        [(-2.0, -3), (-2.0, 2), (-2.0, -2), (0.0, 1), (-2.0, 3), (0.0, -1)]

        >>> l._vsids_calculate()
        -3

        >>> l.lit_heap
        [(-2.0, -2), (-2.0, 2), (0.0, -1), (0.0, 1), (-2.0, 3)]

        """
        if len(self.lit_heap) == 0:
            return 0

        # Clean out the front of the heap as long the variables are set
        while self.variable_set[abs(self.lit_heap[0][1])]:
            heappop(self.lit_heap)
            if len(self.lit_heap) == 0:
                return 0

        return heappop(self.lit_heap)[1]

    def _vsids_lit_assigned(self, lit):
        """Handle the assignment of a literal for the VSIDS heuristic."""
        pass

    def _vsids_lit_unset(self, lit):
        """Handle the unsetting of a literal for the VSIDS heuristic.

        Examples
        ========

        >>> from sympy.logic.algorithms.dpll2 import SATSolver
        >>> l = SATSolver([{2, -3}, {1}, {3, -3}, {2, -2},
        ... {3, -2}], {1, 2, 3}, set())
        >>> l.lit_heap
        [(-2.0, -3), (-2.0, 2), (-2.0, -2), (0.0, 1), (-2.0, 3), (0.0, -1)]

        >>> l._vsids_lit_unset(2)

        >>> l.lit_heap
        [(-2.0, -3), (-2.0, -2), (-2.0, -2), (-2.0, 2), (-2.0, 3), (0.0, -1),
        ...(-2.0, 2), (0.0, 1)]

        """
        var = abs(lit)
        heappush(self.lit_heap, (self.lit_scores[var], var))
        heappush(self.lit_heap, (self.lit_scores[-var], -var))

    def _vsids_clause_added(self, cls):
        """Handle the addition of a new clause for the VSIDS heuristic.

        Examples
        ========

        >>> from sympy.logic.algorithms.dpll2 import SATSolver
        >>> l = SATSolver([{2, -3}, {1}, {3, -3}, {2, -2},
        ... {3, -2}], {1, 2, 3}, set())

        >>> l.num_learned_clauses
        0
        >>> l.lit_scores
        {-3: -2.0, -2: -2.0, -1: 0.0, 1: 0.0, 2: -2.0, 3: -2.0}

        >>> l._vsids_clause_added({2, -3})

        >>> l.num_learned_clauses
        1
        >>> l.lit_scores
        {-3: -1.0, -2: -2.0, -1: 0.0, 1: 0.0, 2: -1.0, 3: -2.0}

        """
        self.num_learned_clauses += 1
        for lit in cls:
            self.lit_scores[lit] += 1

    ########################
    #   Clause Learning    #
    ########################
    def _simple_add_learned_clause(self, cls):
        """Add a new clause to the theory.

        Examples
        ========

        >>> from sympy.logic.algorithms.dpll2 import SATSolver
        >>> l = SATSolver([{2, -3}, {1}, {3, -3}, {2, -2},
        ... {3, -2}], {1, 2, 3}, set())

        >>> l.num_learned_clauses
        0
        >>> l.clauses
        [[2, -3], [1], [3, -3], [2, -2], [3, -2]]
        >>> l.sentinels
        {-3: {0, 2}, -2: {3, 4}, 2: {0, 3}, 3: {2, 4}}

        >>> l._simple_add_learned_clause([3])

        >>> l.clauses
        [[2, -3], [1], [3, -3], [2, -2], [3, -2], [3]]
        >>> l.sentinels
        {-3: {0, 2}, -2: {3, 4}, 2: {0, 3}, 3: {2, 4, 5}}

        """
        cls_num = len(self.clauses)
        self.clauses.append(cls)

        for lit in cls:
            self.occurrence_count[lit] += 1

        self.sentinels[cls[0]].add(cls_num)
        self.sentinels[cls[-1]].add(cls_num)

        self.heur_clause_added(cls)

    def _simple_compute_conflict(self):
        """ Build a clause representing the fact that at least one decision made
        so far is wrong.

        Examples
        ========

        >>> from sympy.logic.algorithms.dpll2 import SATSolver
        >>> l = SATSolver([{2, -3}, {1}, {3, -3}, {2, -2},
        ... {3, -2}], {1, 2, 3}, set())
        >>> next(l._find_model())
        {1: True, 2: False, 3: False}
        >>> l._simple_compute_conflict()
        [3]

        """
        return [-(level.decision) for level in self.levels[1:]]

    def _simple_clean_clauses(self):
        """Clean up learned clauses."""
        pass


class Level:
    """
    Represents a single level in the DPLL algorithm, and contains
    enough information for a sound backtracking procedure.
    """

    def __init__(self, decision, flipped=False):
        self.decision = decision
        self.var_settings = set()
        self.flipped = flipped

    def __repr__(self):
        return "<Level decision=%s, flipped=%s, var_settings=%s>" % (
            self.decision, self.flipped, self.var_settings)
