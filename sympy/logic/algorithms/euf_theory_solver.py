from collections import defaultdict
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
    """Extended congruence closure with proper function application and negation tracking."""

    def __init__(self, equations):
        """Initialize proof-producing congruence closure with enhanced function tracking."""
        super().__init__(equations)

        # Existing proof structures
        self.proof_forest = {}
        self.direct_explanations = {}
        self.merge_sequence = []
        self.equalities_for_term = defaultdict(set)

        # Enhanced function tracking
        self.function_applications = defaultdict(list)  # Maps function_name -> [(args, func_app, reason, is_positive), ...]
        self.constant_equalities = {}  # Maps func_app -> constant
        self.function_negations = defaultdict(list)  # Maps function_name -> [(args, func_app, reason), ...] for negated apps

    def add_equality_with_reason(self, lhs, rhs, reason):
        """Add equality constraint and record the reason for proof production."""
        lhs_id = self._flatten(lhs)
        rhs_id = self._flatten(rhs)

        rep_lhs = self._find(lhs_id)
        rep_rhs = self._find(rhs_id)

        if rep_lhs == rep_rhs:
            return

        # Track function applications (both positive and negative)
        self._track_function_applications(lhs, rhs, reason)

        # Existing logic
        key = _ordered_pair(lhs_id, rhs_id)
        self.direct_explanations[key] = reason

        proof_key = _ordered_pair(rep_lhs, rep_rhs)
        self.proof_forest[proof_key] = reason

        self.merge_sequence.append((rep_lhs, rep_rhs, reason, lhs_id, rhs_id))

        self.equalities_for_term[lhs_id].add(reason)
        self.equalities_for_term[rhs_id].add(reason)

        super().add_equality(lhs, rhs)

    def add_disequality_with_reason(self, lhs, rhs, reason):
        """Add disequality and track negated function applications."""
        # Check if this is a negated function application
        if self._is_function_application(lhs) and self._is_constant(rhs):
            func_name = lhs.func.__name__
            args = lhs.args
            self.function_negations[func_name].append((args, lhs, reason))
        elif self._is_function_application(rhs) and self._is_constant(lhs):
            func_name = rhs.func.__name__
            args = rhs.args
            self.function_negations[func_name].append((args, rhs, reason))

    def _track_function_applications(self, lhs, rhs, reason):
        """Track function applications for congruence explanations."""
        func_app, constant = None, None

        if self._is_function_application(lhs) and not self._is_function_application(rhs):
            func_app, constant = lhs, rhs
        elif self._is_function_application(rhs) and not self._is_function_application(lhs):
            func_app, constant = rhs, lhs

        if func_app and constant:
            if hasattr(func_app, 'func') and hasattr(func_app.func, '__name__'):
                func_name = func_app.func.__name__
                args = func_app.args

                # Store positive function application
                self.function_applications[func_name].append((args, func_app, reason, True))
                self.constant_equalities[func_app] = constant

    def _is_function_application(self, expr):
        """Check if expression is a function application."""
        return (hasattr(expr, 'func') and
                hasattr(expr.func, '__name__') and
                hasattr(expr, 'args') and
                len(expr.args) > 0)

    def _is_constant(self, expr):
        """Check if expression is a constant (like _c0, _c1, etc.)."""
        return (hasattr(expr, 'name') and
                isinstance(expr.name, str) and
                (expr.name.startswith('_c') | expr.name.startswith('c')))

    def explain_equality(self, a, b):
        """Generate explanation for why terms a and b are equal with proper function congruence."""
        a_id = self._flatten(a)
        b_id = self._flatten(b)

        if a_id == b_id:
            return set()

        # Check for function congruence explanation first
        func_explanation = self._explain_function_congruence(a, b)
        if func_explanation:
            return func_explanation

        # Fall back to DFS
        return self._explain_via_dfs(a, b)

    def _explain_function_congruence(self, a, b):
        """Explain equality via function congruence with proper argument tracking."""
        # Case: function application equals constant
        if self._is_function_application(a) and self._is_constant(b):
            return self._explain_func_app_equals_constant(a, b)
        elif self._is_function_application(b) and self._is_constant(a):
            return self._explain_func_app_equals_constant(b, a)

        # Case: both are function applications of same function
        if (self._is_function_application(a) and self._is_function_application(b) and
            hasattr(a, 'func') and hasattr(b, 'func') and a.func == b.func):
            return self._explain_same_function_applications(a, b)

        return None

    def _explain_func_app_equals_constant(self, func_app, constant):
        """Explain why func_app = constant using congruence with another function application."""
        if not hasattr(func_app, 'func') or not hasattr(func_app.func, '__name__'):
            return None

        func_name = func_app.func.__name__
        target_args = func_app.args

        # Look for another function application with same function name and same constant
        for stored_args, stored_func_app, stored_reason, is_positive in self.function_applications[func_name]:
            if stored_func_app == func_app:
                continue  # Skip self

            # Check if it equals the same constant
            stored_constant = self.constant_equalities.get(stored_func_app)
            if stored_constant != constant:
                continue

            # Check arguments have same arity
            if len(stored_args) != len(target_args):
                continue

            # Collect argument equality explanations
            args_explanation = set()
            all_args_equal = True

            for stored_arg, target_arg in zip(stored_args, target_args):
                if self._find(self._flatten(stored_arg)) == self._find(self._flatten(target_arg)):
                    # Get explanation for why these arguments are equal
                    if stored_arg != target_arg:  # Not trivially equal
                        arg_explanation = self._explain_via_dfs(stored_arg, target_arg)
                        args_explanation.update(arg_explanation)
                else:
                    all_args_equal = False
                    break

            if all_args_equal:
                # Found valid congruence explanation
                explanation = {stored_reason}  # Original function assertion
                explanation.update(args_explanation)  # All argument equalities
                return explanation

        # Also check against negated function applications (for conflicts)
        for neg_args, neg_func_app, neg_reason in self.function_negations[func_name]:
            if len(neg_args) != len(target_args):
                continue

            # Check if arguments are equal
            args_explanation = set()
            all_args_equal = True

            for neg_arg, target_arg in zip(neg_args, target_args):
                if self._find(self._flatten(neg_arg)) == self._find(self._flatten(target_arg)):
                    if neg_arg != target_arg:
                        arg_explanation = self._explain_via_dfs(neg_arg, target_arg)
                        args_explanation.update(arg_explanation)
                else:
                    all_args_equal = False
                    break

            if all_args_equal:
                # This would be a conflict case, but for explanation purposes
                # we return the negation reason and argument equalities
                explanation = {neg_reason}  # Negated function assertion
                explanation.update(args_explanation)  # All argument equalities
                return explanation

        return None

    def _explain_same_function_applications(self, func_app1, func_app2):
        """Explain why f(a1, ..., an) = f(b1, ..., bn) using argument equalities."""
        args1 = func_app1.args
        args2 = func_app2.args

        if len(args1) != len(args2):
            return None

        explanation = set()
        for arg1, arg2 in zip(args1, args2):
            if self._find(self._flatten(arg1)) == self._find(self._flatten(arg2)):
                if arg1 != arg2:  # Not trivially equal
                    arg_explanation = self._explain_via_dfs(arg1, arg2)
                    explanation.update(arg_explanation)
            else:
                # Arguments not equal, cannot explain via congruence
                return None

        return explanation

    def _explain_via_dfs(self, a, b):
        """DFS-based explanation for equality."""
        a_id = self._flatten(a)
        b_id = self._flatten(b)

        if a_id == b_id:
            return set()

        visited = set()
        path = []

        def dfs(current):
            if current == b_id:
                return True

            visited.add(current)

            for eq in self.equalities_for_term[current]:
                if hasattr(eq, 'lhs') and hasattr(eq, 'rhs'):
                    if self._flatten(eq.lhs) == current:
                        other = eq.rhs
                    elif self._flatten(eq.rhs) == current:
                        other = eq.lhs
                    else:
                        continue

                    other_id = self._flatten(other)

                    if other_id not in visited:
                        path.append(eq)
                        if dfs(other_id):
                            return True
                        path.pop()

            return False

        found = dfs(a_id)
        return set(path) if found else set()



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

                lhs_diseq = diseq_cause.lhs
                rhs_diseq = diseq_cause.rhs
                if self.cc._find(lhs_diseq) == rep_a and self.cc._find(rhs_diseq) == rep_b:
                    explanation.update(self.cc.explain_equality(lhs_diseq, a))
                    explanation.update(self.cc.explain_equality(rhs_diseq, b))
                else:
                    explanation.update(self.cc.explain_equality(lhs_diseq, b))
                    explanation.update(self.cc.explain_equality(rhs_diseq, a))
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

            # Record disequality
            diseq_key = _ordered_pair(rep_lhs, rep_rhs)
            self.disequalities_set[rep_lhs].add(rep_rhs)
            self.disequalities_set[rep_rhs].add(rep_lhs)
            self.disequality_causes[diseq_key] = ceq

            # Track negated function applications
            self.cc.add_disequality_with_reason(lhs, rhs, ceq)

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
