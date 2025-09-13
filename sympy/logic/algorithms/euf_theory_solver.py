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
    """
    Extended congruence closure with proof forest for explanations and function congruence support.

    This class implements a proof-producing variant of congruence closure for Equality with
    Uninterpreted Functions (EUF) theory. It extends the basic congruence closure algorithm
    to track reasons (proof certificates) for why terms become equal, enabling the generation
    of minimal explanations for equalities and conflicts in DPLL(T) solvers.

    Key Features:
    - Proof forest: Maps representative pairs to the equality literal that caused their merge
    - Direct explanations: Maps term pairs to literals that directly asserted their equality
    - Function congruence: Tracks function applications and their reasoning
    - Explanation generation: Uses DFS to find minimal sets of equalities explaining why a=b

    Attributes:
        proof_forest (dict): Maps (rep_a, rep_b) -> equality_literal that caused merge
        direct_explanations (dict): Maps (term_a, term_b) -> equality_literal for direct assertions
        merge_sequence (list): Chronological log of all merges with reasons
        equalities_for_term (defaultdict): Maps flattened_term -> {equality_literals}
        function_applications (defaultdict): Maps function_name -> [(args, func_app, reason, is_positive)]
        constant_equalities (dict): Maps func_app -> constant it equals
        function_negations (defaultdict): Maps function_name -> [(args, func_app, reason)] for negated apps
    """

    def __init__(self, equations):
        """
        Initialize proof-producing congruence closure with enhanced function tracking.

        Parameters:
            equations (list): Initial set of equality constraints to process

        Note:
            This extends the base EUFCongruenceClosure with proof tracking capabilities.
            All data structures are initialized to empty states.
        """
        super().__init__(equations)

        # Proof tracking structures
        self.proof_forest = {}  # (rep_a, rep_b) -> equality_literal
        self.direct_explanations = {}  # (term_a, term_b) -> equality_literal
        self.merge_sequence = []  # [(rep_a, rep_b, reason, orig_a, orig_b), ...]
        self.equalities_for_term = defaultdict(set)  # term -> {equality_literals}

        # Function application tracking
        self.function_applications = defaultdict(list)  # func_name -> [(args, func_app, reason, is_positive)]
        self.constant_equalities = {}  # func_app -> constant
        self.function_negations = defaultdict(list)  # func_name -> [(args, func_app, reason)]

    def add_equality_with_reason(self, lhs, rhs, reason):
        """
        Add equality constraint and record the reason for proof production.

        This method extends the base congruence closure to track WHY the equality
        holds by maintaining proof certificates for later explanation generation.

        Parameters:
            lhs (sympy.Expr): Left-hand side of the equality
            rhs (sympy.Expr): Right-hand side of the equality
            reason (sympy.Eq): The equality literal that justifies this constraint

        Returns:
            None

        Side Effects:
            - Updates proof_forest with merge reasoning
            - Updates direct_explanations for term pairs
            - Tracks function applications if present
            - Calls parent class add_equality for actual union-find merge

        Complexity:
            O(alpha(n)) where alpha is inverse Ackermann function
        """
        lhs_id = self._flatten(lhs)
        rhs_id = self._flatten(rhs)

        rep_lhs = self._find(lhs_id)
        rep_rhs = self._find(rhs_id)

        # Early termination if already equal
        if rep_lhs == rep_rhs:
            return

        # Track function applications for congruence reasoning
        self._track_function_applications(lhs, rhs, reason)

        # Record direct explanation between original flattened terms
        key = _ordered_pair(lhs_id, rhs_id)
        self.direct_explanations[key] = reason

        # Record in proof forest between current representatives
        proof_key = _ordered_pair(rep_lhs, rep_rhs)
        self.proof_forest[proof_key] = reason

        # Log chronological merge sequence
        self.merge_sequence.append((rep_lhs, rep_rhs, reason, lhs_id, rhs_id))

        # Update explanation graph edges
        self.equalities_for_term[lhs_id].add(reason)
        self.equalities_for_term[rhs_id].add(reason)

        # Perform actual union-find merge
        super().add_equality(lhs, rhs)

    def add_disequality_with_reason(self, lhs, rhs, reason):
        """
        Add disequality constraint and track negated function applications.

        This method tracks disequalities that involve function applications,
        which is crucial for generating proper conflict explanations when
        function congruence is violated.

        Parameters:
            lhs (sympy.Expr): Left-hand side of the disequality
            rhs (sympy.Expr): Right-hand side of the disequality
            reason (sympy.Ne): The disequality literal that justifies this constraint

        Returns:
            None

        Side Effects:
            - Updates function_negations if a function application is involved
        """
        # Check if this involves a negated function application
        if self._is_function_application(lhs) and self._is_constant(rhs):
            func_name = lhs.func.__name__
            args = lhs.args
            self.function_negations[func_name].append((args, lhs, reason))
        elif self._is_function_application(rhs) and self._is_constant(lhs):
            func_name = rhs.func.__name__
            args = rhs.args
            self.function_negations[func_name].append((args, rhs, reason))

    def _track_function_applications(self, lhs, rhs, reason):
        """
        Track function applications for congruence explanations.

        Identifies and stores function applications that equal constants,
        which enables function congruence reasoning in explanations.

        Parameters:
            lhs (sympy.Expr): Left-hand side of equality
            rhs (sympy.Expr): Right-hand side of equality
            reason (sympy.Eq): Reason for this equality

        Returns:
            None

        Side Effects:
            - Updates function_applications if f(args) = constant pattern found
            - Updates constant_equalities mapping
        """
        func_app, constant = None, None

        # Identify function application = constant pattern
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
        """
        Check if expression is a function application.

        Parameters:
            expr (sympy.Expr): Expression to check

        Returns:
            bool: True if expr is a function application (e.g., f(x), g(a,b)), False otherwise.
        """
        return (hasattr(expr, 'func') and
                hasattr(expr.func, '__name__') and
                hasattr(expr, 'args') and
                len(expr.args) > 0)

    def _is_constant(self, expr):
        """
        Check if expression is a constant symbol (like _c0, _c1, c1, etc.).

        Parameters:
            expr (sympy.Expr): Expression to check

        Returns:
            bool: True if expr appears to be a constant symbol, False otherwise

        Note:
            This heuristically identifies constants by checking for names starting
            with '_c' or 'c', which is the convention used in EUF solvers.
        """
        return (hasattr(expr, 'name') and
                isinstance(expr.name, str) and
                (expr.name.startswith('_c') or expr.name.startswith('c')))

    def explain_equality(self, a, b):
        """
        Generate minimal explanation for why terms a and b are equal.

        This method produces a set of equality literals that together constitute
        a proof/explanation for why a = b holds in the current congruence closure.
        It first tries function congruence reasoning, then falls back to DFS.

        Parameters:
            a (sympy.Expr): First term to explain equality for
            b (sympy.Expr): Second term to explain equality for

        Returns:
            set: Set of sympy.Eq literals that explain why a = b.
                 Empty set if a and b are identical or not equal.

        Algorithm:
            1. Check if terms are identical (return empty set)
            2. Try function congruence explanation first
            3. Fall back to DFS-based explanation through equality graph

        Complexity:
            O(V + E) where V = number of terms, E = number of equalities
        """
        a_id = self._flatten(a)
        b_id = self._flatten(b)

        # Base case: terms are syntactically identical
        if a_id == b_id:
            return set()

        # Try function congruence explanation first (more specific)
        func_explanation = self._explain_function_congruence(a, b)
        if func_explanation:
            return func_explanation

        # Fall back to general DFS explanation
        return self._explain_via_dfs(a, b)

    def _explain_function_congruence(self, a, b):
        """
        Explain equality via function congruence with proper argument tracking.

        Handles three cases:
        1. Function application equals constant: f(x) = c
        2. Constant equals function application: c = f(x)
        3. Two function applications of same function: f(x) = f(y)

        Parameters:
            a (sympy.Expr): First term
            b (sympy.Expr): Second term

        Returns:
            set or None: Set of equality literals explaining the congruence,
                        or None if no function congruence explanation applies
        """
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
        """
        Explain why func_app = constant using congruence with another function application.

        Searches for another application of the same function that equals the same constant,
        then explains via argument equalities and the original function assertion.

        Parameters:
            func_app (sympy.Expr): Function application (e.g., f(x, y))
            constant (sympy.Expr): Constant term (e.g., _c0)

        Returns:
            set or None: Set of equality literals explaining the congruence:
                        {original_func_assertion} U {argument_equalities}
                        Returns None if no valid congruence explanation found

        Algorithm:
            1. Find other applications of same function with same constant value
            2. Check if all corresponding arguments are equal
            3. Collect explanations for argument equalities
            4. Return original assertion + argument equality explanations

        Example:
            Given: f(a) = c, a = b
            Query: Why f(b) = c?
            Returns: {Eq(f(a), c), Eq(a, b)}
        """
        if not hasattr(func_app, 'func') or not hasattr(func_app.func, '__name__'):
            return None

        func_name = func_app.func.__name__
        target_args = func_app.args

        # Search positive function applications
        for stored_args, stored_func_app, stored_reason, is_positive in self.function_applications[func_name]:
            if stored_func_app == func_app:
                continue  # Skip self-reference

            # Check if stored application equals same constant
            stored_constant = self.constant_equalities.get(stored_func_app)
            if stored_constant != constant:
                continue

            # Check argument arity compatibility
            if len(stored_args) != len(target_args):
                continue

            # Collect argument equality explanations
            args_explanation = set()
            all_args_equal = True

            for stored_arg, target_arg in zip(stored_args, target_args):
                if self._find(self._flatten(stored_arg)) == self._find(self._flatten(target_arg)):
                    if stored_arg != target_arg:  # Non-trivial equality
                        arg_explanation = self._explain_via_dfs(stored_arg, target_arg)
                        args_explanation.update(arg_explanation)
                else:
                    all_args_equal = False
                    break

            if all_args_equal:
                # Found valid congruence explanation
                explanation = {stored_reason}  # Original f(args) = constant
                explanation.update(args_explanation)  # All arg1 = arg2 explanations
                return explanation

        # Also check negated function applications (for conflict scenarios)
        for neg_args, neg_func_app, neg_reason in self.function_negations[func_name]:
            if len(neg_args) != len(target_args):
                continue

            # Check argument equality for negated case
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
                # Explanation involves negated assertion (conflict case)
                explanation = {neg_reason}  # Negated function assertion
                explanation.update(args_explanation)  # Argument equalities
                return explanation

        return None

    def _explain_same_function_applications(self, func_app1, func_app2):
        """
        Explain why f(a1, ..., an) = f(b1, ..., bn) using argument equalities.

        For two applications of the same function to be equal, all corresponding
        arguments must be equal (function congruence property).

        Parameters:
            func_app1 (sympy.Expr): First function application
            func_app2 (sympy.Expr): Second function application

        Returns:
            set or None: Set of equality literals explaining argument equalities,
                        or None if arguments are not all equal

        Precondition:
            func_app1.func == func_app2.func (same function)

        Example:
            Given: a = b, c = d
            Query: Why f(a, c) = f(b, d)?
            Returns: {Eq(a, b), Eq(c, d)}
        """
        args1 = func_app1.args
        args2 = func_app2.args

        if len(args1) != len(args2):
            return None

        explanation = set()
        for arg1, arg2 in zip(args1, args2):
            if self._find(self._flatten(arg1)) == self._find(self._flatten(arg2)):
                if arg1 != arg2:  # Non-trivial argument equality
                    arg_explanation = self._explain_via_dfs(arg1, arg2)
                    explanation.update(arg_explanation)
            else:
                # Arguments not equal - cannot explain via congruence
                return None

        return explanation

    def _explain_via_dfs(self, a, b):
        """
        Generate equality explanation using depth-first search on equality graph.

        This is the fallback explanation method that traverses the equality graph
        (where nodes are terms and edges are equality literals) to find a path
        from term a to term b.

        Parameters:
            a (sympy.Expr): Source term
            b (sympy.Expr): Target term

        Returns:
            set: Set of equality literals forming a path from a to b.
                 Empty set if no path exists or terms are identical.

        Algorithm:
            Uses DFS with backtracking to find path in equality graph.
            Graph construction: equalities_for_term[term] -> {equality_literals}

        Complexity:
            O(V + E) where V = number of terms, E = number of equality literals

        Note:
            This method provides a general explanation mechanism when function
            congruence doesn't apply. It may not produce the most minimal explanation.
        """
        a_id = self._flatten(a)
        b_id = self._flatten(b)

        if a_id == b_id:
            return set()

        visited = set()
        path = []  # Current path of equality literals

        def dfs(current):
            """
            Internal DFS helper function.

            Parameters:
                current: Current term ID in the search

            Returns:
                bool: True if path to b_id found, False otherwise

            Side Effects:
                Modifies 'path' list with current search path
            """
            if current == b_id:
                return True

            visited.add(current)

            # Explore all equality literals involving current term
            for eq in self.equalities_for_term[current]:
                if hasattr(eq, 'lhs') and hasattr(eq, 'rhs'):
                    # Determine other term in this equality
                    if self._flatten(eq.lhs) == current:
                        other = eq.rhs
                    elif self._flatten(eq.rhs) == current:
                        other = eq.lhs
                    else:
                        continue  # This equality doesn't involve current term

                    other_id = self._flatten(other)

                    if other_id not in visited:
                        path.append(eq)
                        if dfs(other_id):
                            return True
                        path.pop()  # Backtrack

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
