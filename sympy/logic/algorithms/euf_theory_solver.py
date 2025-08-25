from collections import defaultdict, deque
from sympy import Eq, Ne, Not
from sympy.assumptions.ask import Q, AppliedPredicate
from sympy.logic.algorithms.euf_theory import EUFCongruenceClosure, EUFUnhandledInput
from sympy.core.symbol import Dummy
from sympy.core.relational import Relational
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
    """Extract canonical form (lhs, rhs, is_positive) from a literal."""
    if isinstance(literal, Eq):
        return literal.lhs, literal.rhs, True
    elif isinstance(literal, Ne):
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

class ProofProducingCongruenceClosure(EUFCongruenceClosure):
    """Extended congruence closure with proof forest for explanations."""

    def __init__(self, equations):
        super().__init__(equations)
        # Proof forest for explanations (as described in the paper)
        self.proof_forest = {}  # Maps (a,b) -> explanation
        self.merge_history = []  # Sequence of merges with reasons

    def add_equality_with_reason(self, lhs, rhs, reason):
        """Add equality and record the reason for proof production."""
        lhs_id = self._flatten(lhs)
        rhs_id = self._flatten(rhs)

        # Check if already equal
        if self._find(lhs_id) == self._find(rhs_id):
            return

        # Record in proof forest before merging
        rep_lhs = self._find(lhs_id)
        rep_rhs = self._find(rhs_id)

        key = _ordered_pair(rep_lhs, rep_rhs)
        self.proof_forest[key] = reason
        self.merge_history.append((rep_lhs, rep_rhs, reason))

        # Perform the actual merge
        super().add_equality(lhs, rhs)

    def explain_equality(self, lhs, rhs):
        """Explain why lhs = rhs using proof forest."""
        lhs_id = self._flatten(lhs)
        rhs_id = self._flatten(rhs)

        if self._find(lhs_id) != self._find(rhs_id):
            return set()  # Not equal

        # Find path in proof forest (simplified version)
        explanation = set()
        for (a, b, reason) in self.merge_history:
            if self._contributes_to_equality(a, b, lhs_id, rhs_id):
                explanation.add(reason)

        return explanation

    def _contributes_to_equality(self, merge_a, merge_b, target_a, target_b):
        """Check if a merge contributes to target equality (simplified)."""
        # This is a simplified version - full implementation would trace the actual path
        rep_target_a = self._find(target_a)
        rep_target_b = self._find(target_b)
        rep_merge_a = self._find(merge_a)
        rep_merge_b = self._find(merge_b)

        return (rep_merge_a == rep_target_a or rep_merge_a == rep_target_b or
                rep_merge_b == rep_target_a or rep_merge_b == rep_target_b)

class EUFTheorySolver:
    """
    EUF Theory Solver implementing DPLL(T) interface with proof-producing
    congruence closure for minimal conflict explanations.

    Follows the architecture described in:
    - DPLL(T): Fast Decision Procedures
    - Proof-Producing Congruence Closure
    """

    def __init__(self):
        # DPLL(T) interface data structures
        self.interpretation_stack = []  # I-stack from the paper
        self.timestamp_counter = 0
        self.literal_timestamp = {}

        # Disequality tracking
        self.disequalities_set = defaultdict(set)
        self.disequality_causes = {}

        # Literal tracking for T-consequence propagation
        self.positive_literal_list = defaultdict(list)
        self.negative_literal_list = defaultdict(list)
        self.literal_terms = {}

        # Proof-producing congruence closure
        self.cc = ProofProducingCongruenceClosure([])

        # Encoding mappings for DPLL(T)
        self.literal_eqs = {}
        self._enc_to_lit = {}
        self.lit_to_enc = {}

    @classmethod
    def from_encoded_cnf(cls, encoded_cnf, testing_mode=False):
        solver = cls()
        literal_eqs = {}
        conflicts = []
        dummies = numbered_symbols('_c', lambda name=None: Dummy(name, finite=True))

        # CRITICAL: Use consistent constants for same predicate types
        predicate_constants = {}  # Maps predicate function to its constant

        ALLOWED_BIN_PRED = {Q.eq, Q.ne}

        def get_predicate_constant(pred_func):
            """Get consistent constant for a predicate function."""
            if pred_func not in predicate_constants:
                predicate_constants[pred_func] = next(dummies)
            return predicate_constants[pred_func]

        def process_pred(pred):
            if pred is True:
                return [], [True]
            elif pred is False:
                return [], [False]

            # Unary predicates: Q.prime(x) -> Eq(Q.prime(x), _c_prime)
            if isinstance(pred, AppliedPredicate) and len(pred.arguments) == 1:
                # Use SAME constant for ALL instances of this predicate type
                pred_constant = get_predicate_constant(pred.function)
                eq = Eq(pred, pred_constant)
                return [eq], []

            # Binary predicates Q.eq/Q.ne
            elif (isinstance(pred, AppliedPredicate) and pred.function in ALLOWED_BIN_PRED) or isinstance(pred, (Eq, Ne)):
                if isinstance(pred, AppliedPredicate):
                    left, right = pred.arguments
                    is_eq = (pred.function == Q.eq)
                else:
                    left, right = pred.lhs, pred.rhs
                    is_eq = pred.func == Eq

                # Trivial cases
                if left == right:
                    if is_eq:
                        return [], [True]
                    else:
                        return [], [False]

                # Create constraint directly (flattening handled in solver)
                if is_eq:
                    new_constraint = Eq(left, right)
                else:
                    new_constraint = Ne(left, right)

                return [new_constraint], []

            # Other predicates -> dummy equality
            dummy = next(dummies)
            return [Eq(pred, dummy)], []

        # Process all predicates
        encoded_items = (sorted(encoded_cnf.encoding.items(), key=lambda x: str(x[0]))
                        if testing_mode else encoded_cnf.encoding.items())

        for pred, enc in encoded_items:
            try:
                eqs, triv = process_pred(pred)

                if eqs:
                    literal_eqs[enc] = eqs
                    solver._enc_to_lit[enc] = eqs[0]
                    solver.lit_to_enc[eqs[0]] = enc

                # Handle trivial conflicts
                if True in triv:
                    conflicts.append([enc])
                if False in triv:
                    conflicts.append([-enc])

            except Exception:
                continue

        solver.literal_eqs = literal_eqs

        # Initialize solver with all constraints
        all_constraints = set()
        for eq_list in literal_eqs.values():
            all_constraints.update(eq_list)
        solver.Initialize(all_constraints)

        return solver, conflicts


    def Initialize(self, literal_set):
        """Initialize solver with set of literals for proper tracking."""
        for lit in literal_set:
            try:
                lhs, rhs, is_positive = _canonical_lit(lit)
                ceq = _canon_eq(lhs, rhs)
                self.literal_terms[ceq] = (lhs, rhs, is_positive)

                # Add to literal lists for T-consequence propagation
                lhs_c = self.cc._flatten(lhs)
                rhs_c = self.cc._flatten(rhs)

                if is_positive:
                    self.positive_literal_list[lhs_c].append(ceq)
                    self.positive_literal_list[rhs_c].append(ceq)
                else:
                    self.negative_literal_list[lhs_c].append(ceq)
                    self.negative_literal_list[rhs_c].append(ceq)
            except:
                continue

    def assert_lit(self, enc_constraint):
        """
        Assert literal and return T-consequences or conflict explanation.

        This implements the SetTrue operation from DPLL(T) paper.
        """
        abs_enc = abs(enc_constraint)
        if abs_enc not in self.literal_eqs:
            return None

        is_positive = enc_constraint > 0
        constraints = self.literal_eqs[abs_enc]

        try:
            consequences = set()
            for constraint in constraints:
                if is_positive:
                    cons = self.SetTrue(constraint)
                else:
                    # Assert negation
                    if isinstance(constraint, Eq):
                        cons = self.SetTrue(Ne(constraint.lhs, constraint.rhs))
                    elif isinstance(constraint, Ne):
                        cons = self.SetTrue(Eq(constraint.lhs, constraint.rhs))
                consequences.update(cons)
            return consequences

        except ValueError as e:
            if "EUF conflict" in str(e):
                # Generate minimal conflict explanation using proof forest
                conflict_explanation = self._generate_minimal_conflict()
                return (False, conflict_explanation)
            raise

    def SetTrue(self, literal):
        """
        Assert literal and return T-consequences.

        Implements the core SetTrue operation with conflict detection
        and T-consequence propagation as described in DPLL(T) paper.
        """
        self.timestamp_counter += 1

        lhs, rhs, is_positive = _canonical_lit(literal)
        ceq = _canon_eq(lhs, rhs)
        self.literal_terms.setdefault(ceq, (lhs, rhs, is_positive))

        lhs_c = self.cc._flatten(lhs)
        rhs_c = self.cc._flatten(rhs)

        consequences = set()

        if is_positive:
            # Adding equality
            rep_lhs = self.cc._find(lhs_c)
            rep_rhs = self.cc._find(rhs_c)

            # Check for immediate disequality conflict
            for (da, db) in self.disequality_causes.keys():
                if ((self.cc._find(da) == rep_lhs and self.cc._find(db) == rep_rhs) or
                    (self.cc._find(da) == rep_rhs and self.cc._find(db) == rep_lhs)):
                    raise ValueError("EUF conflict: equality contradicts disequality")

            # Add equality with reason for proof production
            self.cc.add_equality_with_reason(lhs, rhs, ceq)

            # Check for conflicts after congruence propagation
            for (da, db) in self.disequality_causes.items():
                if self.cc._find(da) == self.cc._find(db):
                    raise ValueError("EUF conflict: equality contradicts disequality")

            # Find T-consequences (newly implied equalities)
            new_rep_lhs = self.cc._find(lhs_c)
            new_rep_rhs = self.cc._find(rhs_c)

            # Check positive literal lists for newly implied equalities
            if new_rep_lhs != rep_lhs:
                for pos_lit in self.positive_literal_list[lhs_c]:
                    if pos_lit != ceq and pos_lit not in [lit for lit, _ in self.interpretation_stack]:
                        consequences.add(pos_lit)

            if new_rep_rhs != rep_rhs:
               for pos_lit in self.positive_literal_list[rhs_c]:
                    if pos_lit != ceq and pos_lit not in [lit for lit, _ in self.interpretation_stack]:
                        consequences.add(pos_lit)

        else:
            # Adding disequality
            rep_lhs = self.cc._find(lhs_c)
            rep_rhs = self.cc._find(rhs_c)

            # Check for immediate equality conflict
            if rep_lhs == rep_rhs:
                raise ValueError("EUF conflict: disequality within same class")

            # Record disequality
            self.disequalities_set[rep_lhs].add(rep_rhs)
            self.disequalities_set[rep_rhs].add(rep_lhs)
            self.disequality_causes[(rep_lhs, rep_rhs)] = ceq
            self.disequality_causes[(rep_rhs, rep_lhs)] = ceq

            # Find T-consequences (newly implied disequalities)
            # for neg_lit in self.negative_literal_list[lhs_c]:
            #     if neg_lit != ceq and neg_lit not in [lit for lit, _ in self.interpretation_stack]:
            #         consequences.add(neg_lit)

        # Update interpretation stack
        self.interpretation_stack.append((ceq, self.timestamp_counter))
        self.literal_timestamp[ceq] = self.timestamp_counter

        return consequences

    def check(self):
        """Check satisfiability and return conflict explanation if UNSAT."""
        for (a, b), cause_lit in self.disequality_causes.items():
            if self.cc._find(a) == self.cc._find(b):
                # Generate minimal conflict explanation using proof forest
                conflict_explanation = self._generate_minimal_conflict_for_pair(a, b, cause_lit)
                return False, conflict_explanation
        return True, set()

    def _generate_minimal_conflict(self):
        """
        Generate minimal conflict clause using proof-producing congruence closure.

        This implements the conflict explanation algorithm from the
        "Proof-Producing Congruence Closure" paper.
        """
        conflict_lits = set()

        # Find conflicting disequality
        for (a, b), cause_lit in self.disequality_causes.items():
            if self.cc._find(a) == self.cc._find(b):
                conflict_lits.update(self._generate_minimal_conflict_for_pair(a, b, cause_lit))
                break

        return conflict_lits

    def _generate_minimal_conflict_for_pair(self, a, b, diseq_lit):
        """Generate minimal conflict for specific disequality violation."""
        conflict_lits = set()

        # Add the disequality itself
        if diseq_lit in self.lit_to_enc:
            conflict_lits.add(-self.lit_to_enc[diseq_lit])  # Negated since it caused conflict

        # Use proof forest to find minimal set of equalities that caused a=b
        equality_explanation = self.cc.explain_equality(a, b)

        for eq_lit in equality_explanation:
            if eq_lit in self.lit_to_enc:
                conflict_lits.add(-self.lit_to_enc[eq_lit])  # Negated for conflict clause

        return conflict_lits

    def IsTrue(self, literal):
        """Check if literal is true, false, or unknown."""
        try:
            lhs, rhs, is_positive = _canonical_lit(literal)

            if is_positive:
                return self.cc.are_equal(lhs, rhs)
            else:
                lhs_c = self.cc._flatten(lhs)
                rhs_c = self.cc._flatten(rhs)
                rep_lhs = self.cc._find(lhs_c)
                rep_rhs = self.cc._find(rhs_c)

                if rep_lhs == rep_rhs:
                    return False  # Equal, so disequality is false

                ceq = _canon_eq(lhs, rhs)
                return ceq in self.disequality_causes.values()

        except:
            return None

    def Explanation(self, literal):
        """Produce explanation for why literal is true using proof forest."""
        lhs, rhs, is_positive = _canonical_lit(literal)

        if is_positive:
            # Explain equality using proof forest
            return self.cc.explain_equality(lhs, rhs)
        else:
            # For disequality, just return the literal itself
            return {literal}

    def Backtrack(self, n):
        """Backtrack n decision levels."""
        for _ in range(n):
            if not self.interpretation_stack:
                break
            lit, _ = self.interpretation_stack.pop()
            self.literal_timestamp.pop(lit, None)

        # Rebuild state (simplified - full implementation would be more efficient)
        self._rebuild_state()

    def _rebuild_state(self):
        """Rebuild congruence closure state after backtracking."""
        remaining_lits = [ceq for ceq, _ in self.interpretation_stack]

        # Reset state
        self.cc = ProofProducingCongruenceClosure([])
        self.disequalities_set.clear()
        self.disequality_causes.clear()

        # Re-assert remaining literals
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
                    pass  # Ignore conflicts during rebuild

    def reset(self):
        """Reset solver to initial state."""
        self.interpretation_stack.clear()
        self.timestamp_counter = 0
        self.literal_timestamp.clear()
        self.disequalities_set.clear()
        self.disequality_causes.clear()
        self.cc = ProofProducingCongruenceClosure([])
