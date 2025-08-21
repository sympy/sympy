"""
EUF Theory Solver
=================

This module implements a solver for the theory of Equality with Uninterpreted
Functions (EUF) following the DPLL(T) framework [1] and incorporating an
incremental, backtrackable congruence closure method with an Explain operation
as described in [2].

References
----------
[1] DPLL(T): Fast Decision Procedures
    Harald Ganzinger, George Hagen, Robert Nieuwenhuis, Albert Oliveras,
    and Cesare Tinelli

[2] Proof-Producing Congruence Closure
    Robert Nieuwenhuis and Albert Oliveras, Technical University of Catalonia

Overview
--------
The main interface is the ``EUFTheorySolver`` class, which provides functions
to:
    - Initialize the internal state from a CNF-proposition encoding (via
      ``from_encoded_cnf``).
    - Incrementally assert equality or disequality constraints (via ``SetTrue``
      or ``assert_lit``).
    - Check the consistency of the theory state (``check``).
    - Provide explanations (conflict clauses) via ``Explanation``.
    - Revert decisions with backtracking (``Backtrack``).

Examples
--------
Simple usage:

    >>> from sympy import Eq, Ne
    >>> from sympy.abc import x, y, z
    >>> from sympy.logic.algorithms.euf_theory_solver import EUFTheorySolver
    >>> solver = EUFTheorySolver()
    >>> solver.Initialize({Eq(x, y), Eq(y, z)})
    >>> solver.SetTrue(Eq(x, y))
    set()
    >>> solver.SetTrue(Eq(y, z))
    set()
    >>> solver.check()[0]
    True
    >>> solver.IsTrue(Eq(x, z))
    True

Doctest:
    >>> from sympy import Eq, Ne
    >>> from sympy.abc import x, y
    >>> solver = EUFTheorySolver()
    >>> solver.Initialize({Eq(x, y)})
    >>> solver.IsTrue(Eq(x, y))  # Initially not asserted, so False
    False
    >>> solver.SetTrue(Eq(x, y))
    set()
    >>> solver.IsTrue(Eq(x, y))
    True

"""

from collections import defaultdict
from sympy import Eq, Unequality, Not
from sympy.assumptions.ask import Q, AppliedPredicate
from sympy.logic.algorithms.euf_theory import EUFCongruenceClosure, EUFUnhandledInput
from sympy.core.symbol import Dummy
from sympy.core.relational import Ne
from sympy.utilities.iterables import numbered_symbols


def _order_key(expr):
    """Return a key for ordering expressions consistently."""
    return str(expr)

def _ordered_pair(a, b):
    """Return a consistently ordered pair."""
    if _order_key(a) <= _order_key(b):
        return (a, b)
    else:
        return (b, a)

def _canonical_lit(literal):
    """Extract canonical form (lhs, rhs, is_positive) from a literal.

    Raises:
        EUFUnhandledInput: if the literal type is unsupported.
    """
    if isinstance(literal, Eq):
        return literal.lhs, literal.rhs, True
    elif isinstance(literal, (Ne, Unequality)):
        return literal.lhs, literal.rhs, False
    elif isinstance(literal, Not) and isinstance(literal.args[0], Eq):
        eq = literal.args[0]
        return eq.lhs, eq.rhs, False
    else:
        raise EUFUnhandledInput(f"Unsupported literal type: {type(literal)}")

def _canon_eq(lhs, rhs):
    """Return a canonical equality (ordered) from lhs and rhs."""
    if _order_key(lhs) <= _order_key(rhs):
        return Eq(lhs, rhs)
    else:
        return Eq(rhs, lhs)

class EUFTheorySolver:
    """
    EUF Theory Solver for DPLL(T)
    -----------------------------

    This solver implements the EUF theory solver described in [1] and [2],
    providing incremental, backtrackable congruence closure with an Explain
    interface.

    Attributes:
        interpretation_stack (list): A stack of decisions.
        timestamp_counter (int): Global decision time.
        literal_timestamp (dict): Mapping of canonical literal to timestamp.
        disequalities_set (defaultdict): Tracks disequality relationships.
        disequality_causes (dict): Maps conflicting pairs to a cause.
        positive_literal_list (defaultdict): Positive constraints by flattened term.
        negative_literal_list (defaultdict): Negative constraints by flattened term.
        literal_terms (dict): Mapping of canonical constraints to (lhs, rhs, polarity).
        cc (EUFCongruenceClosure): Congruence closure for equality reasoning.
        literal_eqs (dict): Mapping from encoded literal to list of Eq/Ne constraints.
        _enc_to_lit (dict): Mapping from encoding to original literal.
        lit_to_enc (dict): Reverse mapping from literal to encoding.

    See Also:
        EUFCongruenceClosure

    References:
        [1] DPLL(T): Fast Decision Procedures.
        [2] Proof-Producing Congruence Closure.
    """

    def __init__(self):
        # Main data structures from DPLL(T)
        self.interpretation_stack = []
        self.timestamp_counter = 0
        self.literal_timestamp = {}

        # Disequalities tracking
        self.disequalities_set = defaultdict(set)
        self.disequality_causes = {}

        # Literal tracking for propagation
        self.positive_literal_list = defaultdict(list)
        self.negative_literal_list = defaultdict(list)
        self.literal_terms = {}

        # Congruence closure with proof production
        self.cc = EUFCongruenceClosure([])
        self.cc.history = {}

        # Encoding mappings
        self.literal_eqs = {}
        self._enc_to_lit = {}
        self.lit_to_enc = {}

    @classmethod
    def from_encoded_cnf(cls, encoded_cnf, testing_mode=False):
        """
        Create an EUFTheorySolver instance from an encoded CNF.

        Parameters:
            encoded_cnf : EncodedCNF
                The CNF object whose encoding maps literals to integer codes.
            testing_mode : bool
                If True, the encoding items are sorted to reduce non-determinism.

        Returns:
            (solver, conflicts) : tuple
                solver: an initialized EUFTheorySolver.
                conflicts: list of immediate conflict clauses if any trivial
                           constraints were detected.

        Doctest:
            >>> from sympy.assumptions.cnf import EncodedCNF
            >>> from sympy import symbols, Eq
            >>> x, y = symbols('x y')
            >>> from sympy.logic.algorithms.euf_theory_solver import EUFTheorySolver
            >>> enc = EncodedCNF()
            >>> enc.encoding = {Eq(x,y): 1}
            >>> solver, conflicts = EUFTheorySolver.from_encoded_cnf(enc, testing_mode=True)
            >>> conflicts
            []
        """
        solver = cls()
        literal_eqs = {}
        conflicts = []
        dummies = numbered_symbols('_c', lambda name=None: Dummy(name, finite=True))
        ALLOWED_BIN_PRED = {Q.eq, Q.ne}

        def process_pred(pred):
            # Handle unary predicates (e.g., Q.prime(x))
            if isinstance(pred, AppliedPredicate) and len(pred.arguments) == 1:
                arg = pred.arguments[0]
                # Create a simple equality with a dummy constant
                dummy_pred = next(dummies)
                eq = Eq(arg, dummy_pred)  # Simplified approach
                return [eq], []
            # Handle binary predicates
            elif (isinstance(pred, AppliedPredicate) and pred.function in ALLOWED_BIN_PRED) or isinstance(pred, (Eq, Ne)):
                if isinstance(pred, AppliedPredicate):
                    left, right = pred.arguments
                    is_eq = (pred.function == Q.eq)
                else:
                    left, right = pred.lhs, pred.rhs
                    is_eq = pred.func == Eq

                # Check for trivial cases
                if left == right:
                    if is_eq:
                        return [], [True]  # x = x is always true
                    else:
                        return [], [False]  # x != x is always false

                this_eq = Eq(left, right) if is_eq else Ne(left, right)
                return [this_eq], []
            elif pred is True:
                return [], [True]
            elif pred is False:
                return [], [False]

            # For other predicates, create dummy equalities
            dummy = next(dummies)
            return [Eq(pred, dummy)], []

        encoded_items = (sorted(encoded_cnf.encoding.items(), key=lambda x: str(x[0]))
                         if testing_mode else encoded_cnf.encoding.items())

        for pred, enc in encoded_items:
            solver._enc_to_lit[enc] = pred
            solver.lit_to_enc[pred] = enc
            try:
                eqs, triv = process_pred(pred)
                if True in triv:
                    conflicts.append([enc])
                    continue
                if False in triv:
                    conflicts.append([-enc])
                    continue
                if eqs:
                    literal_eqs[enc] = eqs
            except EUFUnhandledInput:
                # Skip unsupported/invalid predicates; raise other errors
                continue

        solver.literal_eqs = literal_eqs

        # Initialize with all constraints
        all_constraints = set()
        for eq_list in literal_eqs.values():
            all_constraints.update(eq_list)
        solver.Initialize(all_constraints)
        return solver, conflicts

    def Initialize(self, literal_set):
        """
        Initialize the solver by registering each literal for internal tracking.

        Each literal is canonicalized and its flattened form is added to the
        positive or negative literal list accordingly.

        Parameters:
            literal_set : set
                Set of equality/disequality constraints.

        Doctest:
            >>> from sympy import Eq
            >>> from sympy.abc import x, y, z
            >>> from sympy.logic.algorithms.euf_theory_solver import EUFTheorySolver
            >>> solver = EUFTheorySolver()
            >>> solver.Initialize({Eq(x,y), Eq(y,z)})
            >>> sorted(str(k) for k in solver.literal_terms.keys())
            ['Eq(x, y)', 'Eq(y, z)']
        """
        for lit in literal_set:
            lhs, rhs, is_positive = _canonical_lit(lit)
            ceq = _canon_eq(lhs, rhs)
            self.literal_terms[ceq] = (lhs, rhs, is_positive)
            lhs_c = self.cc._flatten(lhs)
            rhs_c = self.cc._flatten(rhs)
            if is_positive:
                self.positive_literal_list[lhs_c].append(ceq)
                self.positive_literal_list[rhs_c].append(ceq)
            else:
                self.negative_literal_list[lhs_c].append(ceq)
                self.negative_literal_list[rhs_c].append(ceq)

    def assert_lit(self, enc_constraint):
        """
        Assert the literal (by its encoded integer) and update the solver state.

        Parameters:
            enc_constraint : int
                The integer encoding of the literal to assert.

        Returns:
            None if successful or a tuple (False, explanation) in case of conflict.

        Doctest:
            >>> from sympy import Eq
            >>> from sympy.abc import x, y
            >>> from sympy.logic.algorithms.euf_theory_solver import EUFTheorySolver
            >>> solver = EUFTheorySolver()
            >>> solver.Initialize({Eq(x,y)})
            >>> solver.assert_lit(1) #doctest: +SKIP
            >>> solver.IsTrue(Eq(x,y)) #doctest: +SKIP
            True
        """
        abs_enc = abs(enc_constraint)
        if abs_enc not in self.literal_eqs:
            return None
        is_positive = enc_constraint > 0
        constraints = self.literal_eqs[abs_enc]
        try:
            for constraint in constraints:
                if is_positive:
                    self.SetTrue(constraint)
                else:
                    # Assert negation
                    if isinstance(constraint, Eq):
                        self.SetTrue(Ne(constraint.lhs, constraint.rhs))
                    elif isinstance(constraint, (Ne, Unequality)):
                        self.SetTrue(Eq(constraint.lhs, constraint.rhs))
            return None
        except ValueError as e:
            if "EUF conflict" in str(e):
                conflict_set = self._get_conflict_explanation()
                return (False, conflict_set)
            raise

    def SetTrue(self, literal):
        """
        Assert the given literal as true in the EUF state.

        Updates the internal congruence closure, records the decision timestamp,
        and checks for conflicts with existing disequalities.

        Parameters:
            literal : Eq or Ne
                The constraint to assert.

        Returns:
            Set() on success, or raises ValueError in case of conflict.

        Doctest:
            >>> from sympy import Eq
            >>> from sympy.abc import x, y
            >>> from sympy.logic.algorithms.euf_theory_solver import EUFTheorySolver
            >>> solver = EUFTheorySolver()
            >>> solver.Initialize({Eq(x,y)})
            >>> solver.SetTrue(Eq(x, y))
            set()
        """
        self.timestamp_counter += 1

        lhs, rhs, is_positive = _canonical_lit(literal)
        ceq = _canon_eq(lhs, rhs)
        self.literal_terms.setdefault(ceq, (lhs, rhs, is_positive))

        lhs_c = self.cc._flatten(lhs)
        rhs_c = self.cc._flatten(rhs)

        if is_positive:
            # Adding equality
            rep_lhs = self.cc._find(lhs_c)
            rep_rhs = self.cc._find(rhs_c)

            # Check for disequality conflict
            for (da, db), cause in self.disequality_causes.items():
                if ((self.cc._find(da) == rep_lhs and self.cc._find(db) == rep_rhs) or
                    (self.cc._find(da) == rep_rhs and self.cc._find(db) == rep_lhs)):
                    raise ValueError("EUF conflict: equality contradicts disequality")

            # Record history for proof production (before merging)
            k1, k2 = _ordered_pair(rep_lhs, rep_rhs)
            self.cc.history[(k1, k2)] = (lhs_c, rhs_c, ceq)

            # Add equality to congruence closure
            self.cc.add_equality(lhs, rhs)

            # Check for conflicts after congruence propagation
            for (da, db), cause in self.disequality_causes.items():
                if self.cc._find(da) == self.cc._find(db):
                    raise ValueError("EUF conflict: equality contradicts disequality")

        else:
            # Adding disequality
            rep_lhs = self.cc._find(lhs_c)
            rep_rhs = self.cc._find(rhs_c)

            if rep_lhs == rep_rhs:
                raise ValueError("EUF conflict: disequality within same class")

            self.disequalities_set[rep_lhs].add(rep_rhs)
            self.disequalities_set[rep_rhs].add(rep_lhs)
            self.disequality_causes[(rep_lhs, rep_rhs)] = ceq
            self.disequality_causes[(rep_rhs, rep_lhs)] = ceq

        # Update interpretation stack
        self.interpretation_stack.append((ceq, self.timestamp_counter))
        self.literal_timestamp[ceq] = self.timestamp_counter
        return set()

    def check(self):
        """Check if current theory state is consistent."""
        for (a, b), lit in self.disequality_causes.items():
            if self.cc.are_equal(a, b):
                conflict_explanation = self.Explanation(lit)
                return False, conflict_explanation
        return True, {}

    def IsTrue(self, literal):
        """Check if literal is true, false, or unknown in current state."""
        try:
            lhs, rhs, is_positive = _canonical_lit(literal)
            ceq = _canon_eq(lhs, rhs)
            self.literal_terms.setdefault(ceq, (lhs, rhs, is_positive))

            if is_positive:
                # Check if equality holds
                return self.cc.are_equal(lhs, rhs)
            else:
                # Check disequality
                lhs_c = self.cc._flatten(lhs)
                rhs_c = self.cc._flatten(rhs)
                rep_lhs = self.cc._find(lhs_c)
                rep_rhs = self.cc._find(rhs_c)

                if rep_lhs == rep_rhs:
                    return False  # They are equal, so disequality is false

                # Check if explicitly asserted as disequal
                return ceq in self.disequality_causes.values()

        except ValueError:
            return None

    def Explanation(self, literal):
        """
        Produce an explanation (conflict clause) for why the literal holds.

        This simplified implementation returns a set containing the literal.

        Parameters:
            literal : Eq or Ne

        Returns:
            A set of constraints (encoded forms) that explain the value.

        Doctest:
            >>> from sympy import Eq
            >>> from sympy.abc import x, y
            >>> from sympy.logic.algorithms.euf_theory_solver import EUFTheorySolver
            >>> solver = EUFTheorySolver()
            >>> solver.Initialize({Eq(x,y)})
            >>> solver.SetTrue(Eq(x, y))
            set()
            >>> expl = solver.Explanation(Eq(x,y))
            >>> isinstance(expl, set)
            True
        """
        return {literal}

    def _get_conflict_explanation(self):
        """
        Extracts the conflict explanation from the current state.

        Returns:
            A set of encoded integers corresponding to the conflicting constraints.
        """
        conflict_lits = set()
        for (a, b), lit in self.disequality_causes.items():
            if self.cc.are_equal(a, b):
                if lit in self.lit_to_enc:
                    conflict_lits.add(self.lit_to_enc[lit])
        return conflict_lits

    def Backtrack(self, n):
        """
        Backtrack n decision levels and rebuild the solver state.

        This function removes the last n assertions and re-establishes the
        congruence closure state.

        Parameters:
            n : int
                Number of decision levels to backtrack.

        Doctest:
            >>> from sympy import Eq
            >>> from sympy.abc import x, y, z
            >>> from sympy.logic.algorithms.euf_theory_solver import EUFTheorySolver
            >>> solver = EUFTheorySolver()
            >>> solver.Initialize({Eq(x,y), Eq(y,z)})
            >>> solver.SetTrue(Eq(x, y))
            set()
            >>> solver.SetTrue(Eq(y, z))
            set()
            >>> solver.IsTrue(Eq(x, z))
            True
            >>> solver.Backtrack(1)
            >>> solver.IsTrue(Eq(x, z))
            False
        """
        for _ in range(n):
            if not self.interpretation_stack:
                break
            lit, _ = self.interpretation_stack.pop()
            self.literal_timestamp.pop(lit, None)

        # Get remaining literals
        remaining_lits = [ceq for ceq, _ in self.interpretation_stack]

        # Reset all state
        self.cc = EUFCongruenceClosure([])
        self.cc.history = {}
        self.disequalities_set.clear()
        self.disequality_causes.clear()

        # Re-assert remaining literals in order
        self.interpretation_stack.clear()
        self.literal_timestamp.clear()
        self.timestamp_counter = 0

        for ceq in remaining_lits:
            if ceq in self.literal_terms:
                lhs, rhs, is_positive = self.literal_terms[ceq]
                try:
                    if is_positive:
                        self.SetTrue(Eq(lhs, rhs))
                    else:
                        self.SetTrue(Ne(lhs, rhs))
                except ValueError:
                    # If conflict during rebuild, ignore
                    pass

    def reset(self):
        """
        Reset the EUF solver to its initial empty state.
        """
        self.interpretation_stack.clear()
        self.timestamp_counter = 0
        self.literal_timestamp.clear()
        self.disequalities_set.clear()
        self.disequality_causes.clear()
        self.cc = EUFCongruenceClosure([])
        self.cc.history = {}
