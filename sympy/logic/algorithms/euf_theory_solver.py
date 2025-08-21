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
    """Extract canonical form: (lhs, rhs, is_positive)."""
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
    """Return canonical equality (ordered)."""
    if _order_key(lhs) <= _order_key(rhs):
        return Eq(lhs, rhs)
    else:
        return Eq(rhs, lhs)

class EUFTheorySolver:
    """
    EUF Theory Solver implementing DPLL(T) as described in:
    "DPLL(T): Fast Decision Procedures" by Ganzinger et al.
    """

    def __init__(self):
        # Main data structures from DPLL(T) paper
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

        # Congruence closure with proof support
        self.cc = EUFCongruenceClosure([])
        self.cc.history = {}

        # Encoding mappings
        self.literal_eqs = {}
        self._enc_to_lit = {}
        self.lit_to_enc = {}

    @classmethod
    def from_encoded_cnf(cls, encoded_cnf, testing_mode=False):
        """
        Create EUF solver from encoded CNF.
        Returns (solver, conflicts) where conflicts are immediate contradictions.
        """
        solver = cls()
        literal_eqs = {}
        conflicts = []
        dummies = numbered_symbols('_c', lambda name=None: Dummy(name, finite=True))
        ALLOWED_BIN_PRED = {Q.eq, Q.ne}

        def process_pred(pred):
            # Handle unary predicates: Q.prime(x), etc.
            if isinstance(pred, AppliedPredicate) and len(pred.arguments) == 1:
                arg = pred.arguments[0]
                # Create a simple equality with a dummy constant
                dummy_pred = next(dummies)
                eq = Eq(arg, dummy_pred)  # Simplified approach
                return [eq], []

            # Binary predicates Q.eq/Q.ne and Eq/Ne
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

        encoded_items = (
            sorted(encoded_cnf.encoding.items(), key=lambda x: str(x[0]))
            if testing_mode else encoded_cnf.encoding.items()
        )

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
        """Initialize solver with a set of literals for proper tracking."""
        for lit in literal_set:
                lhs, rhs, is_positive = _canonical_lit(lit)
                ceq = _canon_eq(lhs, rhs)
                self.literal_terms[ceq] = (lhs, rhs, is_positive)

                # Flatten and add to appropriate literal lists
                lhs_c = self.cc._flatten(lhs)
                rhs_c = self.cc._flatten(rhs)

                if is_positive:
                    self.positive_literal_list[lhs_c].append(ceq)
                    self.positive_literal_list[rhs_c].append(ceq)
                else:
                    self.negative_literal_list[lhs_c].append(ceq)
                    self.negative_literal_list[rhs_c].append(ceq)

    def assert_lit(self, enc_constraint):
        """Assert a literal representing a constraint and update internal state."""
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
        """Assert literal to be true in the theory."""
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
        """Produce explanation for why literal is true/false."""
        return {literal}

    def _get_conflict_explanation(self):
        """Extract conflict explanation from current state."""
        conflict_lits = set()
        for (a, b), lit in self.disequality_causes.items():
            if self.cc.are_equal(a, b):
                if lit in self.lit_to_enc:
                    conflict_lits.add(self.lit_to_enc[lit])
        return conflict_lits

    def Backtrack(self, n):
        """Backtrack n decision levels and rebuild state."""
        # Remove last n entries from interpretation stack
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
        """Reset solver to initial state."""
        self.interpretation_stack.clear()
        self.timestamp_counter = 0
        self.literal_timestamp.clear()
        self.disequalities_set.clear()
        self.disequality_causes.clear()
        self.cc = EUFCongruenceClosure([])
        self.cc.history = {}
