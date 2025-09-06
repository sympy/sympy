from collections import defaultdict, deque
from sympy import Eq, Ne, Not
from sympy.assumptions.ask import Q
from sympy.assumptions.assume import AppliedPredicate
from sympy.logic.algorithms.euf_theory import EUFCongruenceClosure, EUFUnhandledInput
from sympy.core.symbol import Dummy
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

class EUFDisequalityContradictionException(Exception):
    """
    Raised when a disequality which is already asserted is being contradicted.
    """

class EUFEqualityContradictionException(Exception):
    """
    Raised when an equality which is already asserted is being contradicted.
    """

class ProofProducingCongruenceClosure(EUFCongruenceClosure):
    """Extended congruence closure with proof forest for explanations."""

    def __init__(self, equations):
        super().__init__(equations)
        # Proof forest for explanations (as described in the paper)
        self.proof_forest = {}  # Maps (rep_a, rep_b) -> reason (literal that caused merge)
        self.direct_explanations = {}  # Maps (a,b) -> literal that directly asserted a=b
        self.merge_sequence = []  # Track order of merges for explanation

        # New: Map flattened term to equality literals that involve it
        self.equalities_for_term = defaultdict(set)

    def add_equality_with_reason(self, lhs, rhs, reason):
        """Add equality and record reason for proof production."""
        lhs_id = self._flatten(lhs)
        rhs_id = self._flatten(rhs)

        rep_lhs = self._find(lhs_id)
        rep_rhs = self._find(rhs_id)

        # Skip if already equal
        if rep_lhs == rep_rhs:
            return

        # Record direct assertion for flattened terms
        key = _ordered_pair(lhs_id, rhs_id)
        self.direct_explanations[key] = reason

        # Record in proof forest and merge sequence
        proof_key = _ordered_pair(rep_lhs, rep_rhs)
        self.proof_forest[proof_key] = reason
        self.merge_sequence.append((rep_lhs, rep_rhs, reason, lhs_id, rhs_id))

        # Record reason for each flattened term involved
        self.equalities_for_term[lhs_id].add(reason)
        self.equalities_for_term[rhs_id].add(reason)

        # Perform the actual merge using base class method
        super().add_equality(lhs, rhs)

    def explain_equality(self, a, b):
        """Explain why a = b using DFS over equality literals connecting flattened terms."""
        a_id = self._flatten(a)
        b_id = self._flatten(b)

        if a_id == b_id:
            return set()

        visited = set()
        path = []
        explanation = set()

        def dfs(current):
            if current == b_id:
                return True
            visited.add(current)
            for eq in self.equalities_for_term[current]:
                other = eq.rhs if eq.lhs == current else eq.lhs
                other_id = self._flatten(other)
                if other_id not in visited:
                    path.append(eq)
                    if dfs(other_id):
                        return True
                    path.pop()
            return False

        found = dfs(a_id)
        if found:
            explanation = set(path)
        else:
            explanation = set()
        return explanation

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

        # Conflict tracking for explanation generation
        self.current_conflict_diseq = None
        self.current_conflict_eq = None
        self.current_conflict_eq_explanation = None

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

                # Create constraint directly
                if is_eq:
                    new_constraint = Eq(left, right)
                else:
                    new_constraint = Ne(left, right)

                return [new_constraint], []

            else:
                raise EUFUnhandledInput(f"Unsupported predicate: {pred}")

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
                    if isinstance(eqs[0], Eq):
                        opp_equation = Ne(eqs[0].lhs, eqs[0].rhs)
                        solver._enc_to_lit[-enc] = opp_equation
                        solver.lit_to_enc[opp_equation] = -enc
                    if isinstance(eqs[0], Ne):
                        opp_equation = Eq(eqs[0].lhs, eqs[0].rhs)
                        solver._enc_to_lit[-enc] = opp_equation
                        solver.lit_to_enc[opp_equation] = -enc

                # Handle trivial conflicts
                if True in triv:
                    conflicts.append([enc])
                if False in triv:
                    conflicts.append([-enc])

            except EUFUnhandledInput:
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


    def explain_disequality(self, a, b):
        """
        Explain why a != b by finding the source disequality and equality explanations.
        Returns set of literals that explain the disequality.
        """
        explanation = set()

        a_id = self.cc._flatten(a)
        b_id = self.cc._flatten(b)
        rep_a = self.cc._find(a_id)
        rep_b = self.cc._find(b_id)

        # Find the source disequality
        for diseq_pair, diseq_cause in self.disequality_causes.items():
            da, db = diseq_pair
            rep_da = self.cc._find(da)
            rep_db = self.cc._find(db)

            if ((rep_da == rep_a and rep_db == rep_b) or
                (rep_da == rep_b and rep_db == rep_a)):
                # Add the source disequality to explanation
                explanation.add(diseq_cause)

                # Now explain why da = a and db = b (or da = b and db = a)
                if rep_da == rep_a and rep_db == rep_b:
                    # Explain da = a and db = b
                    if da != a_id:
                        explanation.update(self.cc.explain_equality(da, a))
                    if db != b_id:
                        explanation.update(self.cc.explain_equality(db, b))
                else:  # rep_da == rep_b and rep_db == rep_a
                    # Explain da = b and db = a
                    if da != b_id:
                        explanation.update(self.cc.explain_equality(da, b))
                    if db != a_id:
                        explanation.update(self.cc.explain_equality(db, a))
                break

        return explanation

    def _generate_minimal_conflict(self):
        """
        Generate minimal conflict clause based on the type of contradiction.
        """
        conflict_clause = set()

        if hasattr(self, '_current_exception_type'):
            if self._current_exception_type == 'disequality_contradiction':
                # EUFDisequalityContradictionException case
                # We have: source disequality + conflicting equality

                # Add the conflicting equality (negated)
                if self.current_conflict_eq is not None:
                    eq_encoding = self.lit_to_enc.get(self.current_conflict_eq)
                    if eq_encoding is not None:
                        conflict_clause.add(-eq_encoding)

                # Get explanation for the disequality that's being violated
                if self.current_conflict_eq is not None:
                    lhs = self.current_conflict_eq.lhs
                    rhs = self.current_conflict_eq.rhs
                    diseq_explanation = self.explain_disequality(lhs, rhs)

                    # Add all literals from disequality explanation (negated)
                    for exp_lit in diseq_explanation:
                        exp_encoding = self.lit_to_enc.get(exp_lit)
                        if exp_encoding is not None:
                            conflict_clause.add(-exp_encoding)

            elif self._current_exception_type == 'equality_contradiction':
                # EUFEqualityContradictionException case
                # We're trying to assert disequality but terms are already equal

                # Add the conflicting disequality (negated)
                if self.current_conflict_diseq is not None:
                    diseq_encoding = self.lit_to_enc.get(self.current_conflict_diseq)
                    if diseq_encoding is not None:
                        conflict_clause.add(-diseq_encoding)
                    else :
                        new_diseq = Ne(self.current_conflict_diseq.rhs, self.current_conflict_diseq.lhs)
                        diseq_encoding = self.lit_to_enc.get(new_diseq)
                        if diseq_encoding is not None:
                            conflict_clause.add(-diseq_encoding)

                # Get explanation for why the terms are already equal
                if self.current_conflict_diseq is not None:
                    lhs = self.current_conflict_diseq.lhs
                    rhs = self.current_conflict_diseq.rhs
                    eq_explanation = self.cc.explain_equality(lhs, rhs)

                    # Add all literals from equality explanation (negated)
                    for exp_lit in eq_explanation:
                        exp_encoding = self.lit_to_enc.get(exp_lit)
                        if exp_encoding is not None:
                            conflict_clause.add(-exp_encoding)
                        else:
                            new_eq = Eq(exp_lit.rhs, exp_lit.lhs)
                            exp_encoding = self.lit_to_enc.get(new_eq)
                            if exp_encoding is not None:
                                conflict_clause.add(-exp_encoding)

        return conflict_clause

    def SetTrue(self, literal):
        """
        Enhanced SetTrue with proper exception type tracking.
        """
        self.timestamp_counter += 1
        lhs, rhs, is_positive = _canonical_lit(literal)
        ceq = _canon_eq(lhs, rhs)
        self.literal_terms.setdefault(ceq, (lhs, rhs, is_positive))

        lhs_c = self.cc._flatten(lhs)
        rhs_c = self.cc._flatten(rhs)
        consequences = set()

        if is_positive:
            # Adding equality a = b
            rep_lhs = self.cc._find(lhs_c)
            rep_rhs = self.cc._find(rhs_c)

            # Check for immediate disequality conflict BEFORE adding equality
            for diseq_pair, diseq_cause in self.disequality_causes.items():
                da, db = diseq_pair
                if ((self.cc._find(da) == rep_lhs and self.cc._find(db) == rep_rhs) or
                    (self.cc._find(da) == rep_rhs and self.cc._find(db) == rep_lhs)):
                    # Found conflict - equality contradicts disequality
                    self.current_conflict_diseq = diseq_cause
                    self.current_conflict_eq = ceq
                    self._current_exception_type = 'disequality_contradiction'
                    raise EUFDisequalityContradictionException("EUF conflict: Equality contradicts Disequality")

            # Add equality with reason for proof production
            if rep_lhs != rep_rhs:
                self.cc.add_equality_with_reason(lhs, rhs, ceq)

            # Check for conflicts after congruence propagation
            for diseq_pair, diseq_cause in list(self.disequality_causes.items()):
                da, db = diseq_pair
                if self.cc._find(da) == self.cc._find(db):
                    # Congruence closure made diseq terms equal
                    self.current_conflict_diseq = diseq_cause
                    self.current_conflict_eq = ceq
                    self._current_exception_type = 'disequality_contradiction'
                    raise EUFDisequalityContradictionException("EUF conflict: Congruence Closure violated Disequality")

        else:
            ceq = Ne(ceq.lhs, ceq.rhs)

            # Adding disequality a != b
            rep_lhs = self.cc._find(lhs_c)
            rep_rhs = self.cc._find(rhs_c)

            # Check for immediate equality conflict
            if rep_lhs == rep_rhs:
                # Terms are already equal - disequality contradicts equality
                self.current_conflict_diseq = ceq  # The disequality we're trying to assert
                self.current_conflict_eq = None
                self._current_exception_type = 'equality_contradiction'
                raise EUFEqualityContradictionException("EUF conflict: disequality within same equivalence class")

            lhs, rhs, is_positive = _canonical_lit(literal)
            if self.cc.are_equal(lhs,rhs):
                # Terms are already equal - disequality contradicts equality
                self.current_conflict_diseq = ceq  # The disequality we're trying to assert
                self.current_conflict_eq = None
                self._current_exception_type = 'equality_contradiction'
                raise EUFEqualityContradictionException("EUF conflict: disequality within same equivalence class")

            # Record disequality
            diseq_key = _ordered_pair(rep_lhs, rep_rhs)
            self.disequalities_set[rep_lhs].add(rep_rhs)
            self.disequalities_set[rep_rhs].add(rep_lhs)
            self.disequality_causes[diseq_key] = ceq

        # Update interpretation stack
        self.interpretation_stack.append((literal, self.timestamp_counter))
        self.literal_timestamp[literal] = self.timestamp_counter

        return consequences

    def assert_lit(self, enc_constraint):
        """
        Enhanced assert_lit with proper conflict handling.
        """
        abs_enc = abs(enc_constraint)

        if abs_enc not in self.literal_eqs:
            return True, set()

        is_positive = enc_constraint > 0
        constraint = self.literal_eqs[abs_enc][0]

        # Clear previous conflict state
        self.current_conflict_diseq = None
        self.current_conflict_eq = None
        self._current_exception_type = None

        try:
            if is_positive:
                consequences = self.SetTrue(constraint)
            else:
                if isinstance(constraint, Eq):
                    new_eq = Ne(constraint.lhs, constraint.rhs)
                    consequences = self.SetTrue(new_eq)
                else:  # Ne
                    new_eq = Eq(constraint.lhs, constraint.rhs)
                    consequences = self.SetTrue(new_eq)

            if consequences is None:
                consequences = set()
            return True, consequences

        except (EUFDisequalityContradictionException, EUFEqualityContradictionException):
            # Generate minimal conflict clause based on exception type
            conflict_explanation = self._generate_minimal_conflict()
            return False, conflict_explanation


    def IsTrue(self, literal):
        """Check if literal is true, false, or unknown."""
        lhs, rhs, is_positive = _canonical_lit(literal)
        lhs_c = self.cc._flatten(lhs)
        rhs_c = self.cc._flatten(rhs)
        rep_lhs = self.cc._find(lhs_c)
        rep_rhs = self.cc._find(rhs_c)
        if is_positive:
            if self.cc.are_equal(lhs, rhs) is True:
                return True
            else:
                for diseq_pair in self.disequality_causes.keys():
                    da, db = diseq_pair
                    poss1 = self.cc._find(da) == rep_lhs and self.cc._find(db) == rep_rhs
                    poss2 = self.cc._find(da) == rep_rhs and self.cc._find(db) == rep_lhs
                    if poss1 or poss2:
                        return False

        else:
            if self.cc.are_equal(lhs,rhs):
                return False # Equal, so disequality is false

            for diseq_pair in self.disequality_causes.keys():
                da, db = diseq_pair
                poss1 = self.cc._find(da) == rep_lhs and self.cc._find(db) == rep_rhs
                poss2 = self.cc._find(da) == rep_rhs and self.cc._find(db) == rep_lhs
                if poss1 or poss2:
                    return True


    def Explanation(self, literal):
        """Produce explanation for why literal is true using proof forest."""
        lhs, rhs, is_positive = _canonical_lit(literal)
        rep_lhs = self.cc._find(lhs)
        rep_rhs = self.cc._find(rhs)
        if is_positive:
            # Explain equality using proof forest
            return self.cc.explain_equality(lhs, rhs)
        else:
            explanation = set()
            eq1_explanation = set()
            eq2_explanation = set()
            for diseq_pair, diseq_cause in self.disequality_causes.items():
                da, db = diseq_pair
                poss1 = self.cc._find(da) == rep_lhs and self.cc._find(db) == rep_rhs
                poss2 = self.cc._find(da) == rep_rhs and self.cc._find(db) == rep_lhs
                if poss1:
                    eq1_explanation = self.cc.explain_equality(da,rep_lhs)
                    eq2_explanation = self.cc.explain_equality(db,rep_rhs)
                if poss2:
                    eq1_explanation = self.cc.explain_equality(da,rep_rhs)
                    eq2_explanation = self.cc.explain_equality(db,rep_lhs)
                if poss1 or poss2:
                    neq = Ne(da,db)
                    explanation.add(self.lit_to_enc[neq])
                    for eq in eq1_explanation:
                        explanation.add(eq)
                    for eq in eq2_explanation:
                        explanation.add(eq)
        return explanation

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

        # Reset conflict tracking
        self.current_conflict_diseq = None
        self.current_conflict_eq = None
        self.current_conflict_eq_explanation = None
        self._current_exception_type = None
