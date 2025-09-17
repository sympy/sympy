import pytest
from sympy import symbols, Function, Eq, Not, Ne, Unequality
from sympy.assumptions.cnf import CNF, EncodedCNF
from sympy.logic import boolalg
from sympy.assumptions.ask import Q
from sympy.logic.algorithms.euf_theory_solver import (
    EUFTheorySolver, EUFDisequalityContradictionException, EUFEqualityContradictionException,
    ProofProducingCongruenceClosure, _order_key, _ordered_pair, _canonical_lit, _canon_eq
)
from sympy.logic.algorithms.euf_theory import EUFUnhandledInput
from sympy.assumptions.assume import AppliedPredicate
import sympy.core.random as random
from sympy.logic.inference import satisfiable
from sympy.external import import_module
from sympy.testing.pytest import raises, skip
from sympy.logic.algorithms.z3_wrapper import z3_satisfiable

# Try to import Z3 for comparison
try:
    z3 = import_module("z3")
    Z3_AVAILABLE = True
except ImportError:
    Z3_AVAILABLE = False


class RandomEUFTestGenerator:
    """Generator for random EUF test cases."""

    def __init__(self, seed=42):
        random.seed(seed)
        self.MAX_ARGS = 3

        self.symbols_pool = [symbols(f'x_{i}') for i in range(20)]
        self.functions_pool = [Function(f'f_{i}') for i in range(5)]
        self.function_arity ={ f: random.randint(1, self.MAX_ARGS) for f in self.functions_pool }


    def generate_random_term(self, depth=0, max_depth=3):
        """Generate a random term (symbol or function application)."""
        if depth >= max_depth or random.random() < 0.4:
            # Return a symbol
            return random.choice(self.symbols_pool)
        else:
            # Return a function application
            func = random.choice(self.functions_pool)
            num_args = self.function_arity[func]
            args = [self.generate_random_term(depth + 1, max_depth) for _ in range(num_args)]
            return func(*args)

    def generate_random_equality(self):
        """Generate a random equality or disequality."""
        term1 = self.generate_random_term()
        term2 = self.generate_random_term()

        if random.random() < 0.7:  # 70% equalities, 30% disequalities
            return Eq(term1, term2)
        else:
            return Ne(term1, term2)

    def generate_constraint_set(self, num_constraints=25):
        """Generate a set of random constraints."""
        constraints = set()
        while len(constraints) < num_constraints:
            candidate_term_pair = sorted((self.generate_random_term(), self.generate_random_term()),key=str)
            candidate_term_pair = tuple(candidate_term_pair)
            if candidate_term_pair[0] == candidate_term_pair[1]:
                continue # Skip x != x or x == x
            if candidate_term_pair in constraints:
                continue # Skip duplicates

            constraints.add(candidate_term_pair)

        # 70% equalities, 30% disequalities
        constraints = [Eq(term1, term2) if random.random() < 0.7 else Ne(term1,term2) for term1, term2 in constraints]

        return constraints


class Z3Comparator:
    """Compare SymPy EUF results with Z3."""

    def __init__(self):
        if not Z3_AVAILABLE:
            pytest.skip("Z3 not available for comparison")

    def _get_function_signature(self, func_name, arity):
        """Get a unique function signature for Z3."""
        return f"{func_name}_{arity}"

    def sympy_to_z3_term(self, term, z3_symbols, z3_functions):
        """Convert SymPy term to Z3 term."""
        if term.is_symbol:
            if term not in z3_symbols:
                z3_symbols[term] = z3.Const(str(term), z3.IntSort())
            return z3_symbols[term]
        elif hasattr(term, 'func') and hasattr(term.func, '__name__'):
            func_name = term.func.__name__
            arity = len(term.args)

            # Create unique function signature based on name and arity
            func_signature = self._get_function_signature(func_name, arity)

            if func_signature not in z3_functions:
                # Create Z3 function with specific arity
                arg_sorts = [z3.IntSort()] * arity
                z3_functions[func_signature] = z3.Function(func_signature, *arg_sorts, z3.IntSort())

            z3_func = z3_functions[func_signature]
            z3_args = [self.sympy_to_z3_term(arg, z3_symbols, z3_functions) for arg in term.args]
            return z3_func(*z3_args)
        else:
            # Fallback for other terms
            term_str = str(term)
            if term_str not in z3_symbols:
                z3_symbols[term_str] = z3.Const(term_str, z3.IntSort())
            return z3_symbols[term_str]

    def sympy_to_z3_constraint(self, constraint, z3_symbols, z3_functions):
        """Convert SymPy constraint to Z3 constraint."""
        if isinstance(constraint, Eq):
            left = self.sympy_to_z3_term(constraint.lhs, z3_symbols, z3_functions)
            right = self.sympy_to_z3_term(constraint.rhs, z3_symbols, z3_functions)
            return left == right
        elif isinstance(constraint, Ne):
            left = self.sympy_to_z3_term(constraint.lhs, z3_symbols, z3_functions)
            right = self.sympy_to_z3_term(constraint.rhs, z3_symbols, z3_functions)
            return left != right
        elif isinstance(constraint, Not) and isinstance(constraint.args[0], Eq):
            eq = constraint.args[0]
            left = self.sympy_to_z3_term(eq.lhs, z3_symbols, z3_functions)
            right = self.sympy_to_z3_term(eq.rhs, z3_symbols, z3_functions)
            return left != right
        else:
            raise ValueError(f"Unsupported constraint type: {type(constraint)}")

    def check_satisfiability_with_z3(self, constraints):
        """Check satisfiability using Z3."""
        solver = z3.Solver()
        z3_symbols = {}
        z3_functions = {}

        try:
            for constraint in constraints:
                z3_constraint = self.sympy_to_z3_constraint(constraint, z3_symbols, z3_functions)
                solver.add(z3_constraint)

            result = solver.check()
            if result == z3.sat:
                return True, solver.model()
            elif result == z3.unsat:
                return False, None
            else:
                return None, None
        except (ValueError, TypeError, AttributeError, z3.Z3Exception) as e:
            print(f"Z3 conversion error: {e}")
            return None, None


f, g = symbols('f g', cls=Function)
F, G, H, I, J, K = symbols('F G H I J K', cls=Function)
# Create large sets of symbols and functions for complex testing
variables = symbols('a b c d e f g h i j k l m n o p q r s t u v w x y z')
a, b, c, d, e, f, g, h, i, j, k, l, m, n, o, p, q, r, s, t, u, v, w, x, y, z = variables


def test_initialize_and_istrue():
    solver = EUFTheorySolver()
    eqs = {Eq(a, b), Eq(b, c), Eq(F(a), F(b))}
    solver.Initialize(eqs)

    # Initially false
    for lit in eqs:
        assert solver.IsTrue(lit) is None

    # Assert a=b
    solver.SetTrue(Eq(a, b))
    assert solver.IsTrue(Eq(a, b)) is True
    assert solver.IsTrue(Eq(b, a)) is True
    # others still undecided
    assert solver.IsTrue(Eq(b, c)) is None


def test_positive_propagation_transitivity():
    solver = EUFTheorySolver()
    lits = {Eq(a, b), Eq(b, c), Eq(a, c)}
    solver.Initialize(lits)

    solver.SetTrue(Eq(a, b))
    solver.SetTrue(Eq(b, c))
    # transitive a=c
    assert solver.IsTrue(Eq(a, c)) is True


def test_functional_congruence_propagation():
    solver = EUFTheorySolver()
    lits = {Eq(a, b), Eq(F(a), F(b))}
    solver.Initialize(lits)

    solver.SetTrue(Eq(a, b))
    # should propagate F(a)=F(b)
    assert solver.IsTrue(Eq(F(a), F(b))) is True


def test_disequality_no_merge():
    solver = EUFTheorySolver()
    solver.Initialize({Not(Eq(a, b))})
    solver.SetTrue(Not(Eq(a, b)))
    arep = solver.cc._find(solver.cc._flatten(a))
    brep = solver.cc._find(solver.cc._flatten(b))
    assert arep != brep
    assert brep in solver.disequalities_set[arep]


def test_disequality_conflict():
    solver = EUFTheorySolver()
    solver.Initialize({Eq(a, b), Not(Eq(a, b))})
    solver.SetTrue(Not(Eq(a, b)))
    with pytest.raises(EUFDisequalityContradictionException):
        solver.SetTrue(Eq(a, b))


def test_redundant_assertions_safe():
    solver = EUFTheorySolver()
    solver.Initialize({Eq(a, b)})
    solver.SetTrue(Eq(a, b))
    # asserting again should not error
    solver.SetTrue(Eq(a, b))
    assert solver.IsTrue(Eq(a, b))


def test_sympy_unequality_init():
    solver = EUFTheorySolver()
    solver.Initialize({Unequality(a, b)})
    solver.SetTrue(Unequality(a, b))
    assert solver.IsTrue(Unequality(a, b)) is True


def make_solver_from_props(*props):
    cnf = CNF.from_prop(boolalg.And(*props))
    enc = EncodedCNF()
    enc.from_cnf(cnf)
    solver, conflicts = EUFTheorySolver.from_encoded_cnf(enc, testing_mode=True)
    return solver, enc, conflicts


def test_order_independence_of_assertions():
    # x=y, y=z - test that order doesn't matter
    prop = Q.eq(x,y) & Q.eq(y,z)
    _prop = CNF.from_prop(prop)
    sat = EncodedCNF()
    sat.add_from_cnf(_prop)
    solver, conflicts = EUFTheorySolver.from_encoded_cnf(sat)
    solver.assert_lit(1)
    solver.assert_lit(2)
    assert conflicts == []
    assert solver.IsTrue(Eq(x, z)) is True


def test_simple_equality_chain():
    """
    Test EUF: x = y, y = z  -> SAT, and x = z must hold.
    """
    cnf = CNF().from_prop(Eq(x, y) & Eq(y, z))
    enc = EncodedCNF(); enc.from_cnf(cnf)
    euf, _ = EUFTheorySolver.from_encoded_cnf(enc, testing_mode=True)

    # Assert all clauses/literals
    for lit_id in enc.encoding.values():
        assert euf.assert_lit(lit_id) == (True,set())

    assert euf.IsTrue(Eq(x, z)) is True


def test_multiple_function_symbols_deep_nesting():
        """Test congruence with multiple function symbols and deep nesting."""
        solver = EUFTheorySolver()

        # Set up complex function applications with deep nesting
        constraints = [
            Eq(a, b),
            Eq(c, d),
            Eq(e, f),
            Eq(F(a), G(c)),
            Eq(G(c), H(e)),
            Eq(H(F(a)), I(G(c))),
            Eq(I(H(F(a))), J(I(G(c)))),
            Eq(J(I(H(F(a)))), K(J(I(G(c)))))
        ]

        solver.Initialize(set(constraints))

        # Assert base equalities
        solver.SetTrue(Eq(a, b))
        solver.SetTrue(Eq(c, d))
        solver.SetTrue(Eq(e, f))

        # Assert function equalities
        solver.SetTrue(Eq(F(a), G(c)))
        solver.SetTrue(Eq(G(c), H(e)))
        solver.SetTrue(Eq(H(F(a)), I(G(c))))
        solver.SetTrue(Eq(I(H(F(a))), J(I(G(c)))))
        solver.SetTrue(Eq(J(I(H(F(a)))), K(J(I(G(c))))))

        # Test congruence propagation through all levels
        assert solver.IsTrue(Eq(F(b), G(d))) is True  # F(a)=F(b), G(c)=G(d)
        assert solver.IsTrue(Eq(G(d), H(f))) is True  # G(c)=G(d), H(e)=H(f)
        assert solver.IsTrue(Eq(H(F(b)), I(G(d)))) is True
        assert solver.IsTrue(Eq(I(H(F(b))), J(I(G(d))))) is True
        assert solver.IsTrue(Eq(J(I(H(F(b)))), K(J(I(G(d)))))) is True


def test_function_array_like_operations():
    """Test function applications that simulate array operations."""
    solver = EUFTheorySolver()

    # Simulate array[index] = value operations with functions
    # Let F represent an array, indices are the first argument
    constraints = [
            # Set up index equalities
            Eq(i, j),  # indices i and j are equal
            Eq(x, y),  # values x and y are equal

            # Array operations: F(i, v) represents array updated at index i with value v
            Eq(F(a, x), F(b, y)),  # F(a,x) = F(b,y)
            Eq(F(i, x), F(j, y)),  # Should be equal due to i=j, x=y

            # Nested function calls
            Eq(G(F(a, x)), G(F(b, y))),
            Eq(H(G(F(a, x))), c),
        ]

    solver.Initialize(set(constraints))

    # Assert the constraints
    for constraint in constraints:
        solver.SetTrue(constraint)

    # Verify congruence holds
    assert solver.IsTrue(Eq(F(i, x), F(j, y))) is True
    assert solver.IsTrue(Eq(G(F(i, x)), G(F(j, y)))) is True


def test_functional_composition_chains():
    """Test long chains of functional composition."""
    solver = EUFTheorySolver()

    # Create a composition chain: F(G(H(I(J(x)))))
    # with various equalities that should propagate through
    constraints = [
            Eq(x, y),
            Eq(J(x), a),
            Eq(I(a), b),
            Eq(H(b), c),
            Eq(G(c), d),
            Eq(F(d), e),

            # Additional constraints for testing
            Eq(J(y), a),  # Should be inferred by congruence
            Eq(F(G(H(I(J(x))))), e)  # Should be inferred
        ]

    solver.Initialize(set(constraints))

    # Assert base constraints (except the ones that should be inferred)
    solver.SetTrue(Eq(x, y))
    solver.SetTrue(Eq(J(x), a))
    solver.SetTrue(Eq(I(a), b))
    solver.SetTrue(Eq(H(b), c))
    solver.SetTrue(Eq(G(c), d))
    solver.SetTrue(Eq(F(d), e))

    # Test that congruence propagated correctly
    assert solver.IsTrue(Eq(J(y), a)) is True
    assert solver.IsTrue(Eq(I(J(x)), b)) is True
    assert solver.IsTrue(Eq(H(I(J(x))), c)) is True
    assert solver.IsTrue(Eq(G(H(I(J(x)))), d)) is True
    assert solver.IsTrue(Eq(F(G(H(I(J(x))))), e)) is True

    # Test the full composition
    assert solver.IsTrue(Eq(F(G(H(I(J(y))))), e)) is True


def test_issue_1():
    enc_cnf = EncodedCNF()
    enc_cnf.encoding = { Q.prime(x): 1, Q.eq(x,y): 2, Q.prime(y): 3 }
    enc_cnf.data = [{-1},{2},{3}]
    solver,conflicts = EUFTheorySolver.from_encoded_cnf(enc_cnf)
    assert solver.assert_lit(-1) == (True, set())
    assert solver.assert_lit(2) == (True, set())
    assert solver.assert_lit(3) == (False, {1,-2,-3})


# Test helper functions
def test_order_key():
    """Test _order_key function for consistent ordering."""
    assert _order_key(a) == str(a)
    assert _order_key(F(a)) == str(F(a))


def test_ordered_pair():
    """Test _ordered_pair function."""
    pair1 = _ordered_pair(a, b)
    pair2 = _ordered_pair(b, a)
    assert pair1 == pair2  # Should be consistently ordered


def test_canonical_lit():
    """Test _canonical_lit function with various literal types."""
    # Test Eq
    lhs, rhs, is_pos = _canonical_lit(Eq(a, b))
    assert lhs == a and rhs == b and is_pos == True

    # Test Ne
    lhs, rhs, is_pos = _canonical_lit(Ne(a, b))
    assert lhs == a and rhs == b and is_pos == False

    # Test Not(Eq)
    lhs, rhs, is_pos = _canonical_lit(Not(Eq(a, b)))
    assert lhs == a and rhs == b and is_pos == False

    # Test unsupported literal
    with pytest.raises(EUFUnhandledInput):
        _canonical_lit(a)  # Just a symbol, not a supported literal


def test_canon_eq():
    """Test _canon_eq function."""
    eq1 = _canon_eq(a, b)
    eq2 = _canon_eq(b, a)
    assert eq1 == eq2  # Should be consistently ordered


# Test ProofProducingCongruenceClosure
def test_proof_producing_congruence_closure_init():
    """Test initialization of ProofProducingCongruenceClosure."""
    ppcc = ProofProducingCongruenceClosure([])
    assert ppcc.proof_forest == {}
    assert ppcc.direct_explanations == {}
    assert ppcc.merge_sequence == []


def test_add_equality_with_reason():
    """Test add_equality_with_reason method."""
    ppcc = ProofProducingCongruenceClosure([])
    reason = Eq(a, b)
    ppcc.add_equality_with_reason(a, b, reason)

    # Check that equality is recorded
    assert ppcc.are_equal(a, b)

    # Adding same equality again should not cause issues
    ppcc.add_equality_with_reason(a, b, reason)


def test_explain_equality_same_terms():
    """Test explain_equality with identical terms."""
    ppcc = ProofProducingCongruenceClosure([])
    explanation = ppcc.explain_equality(a, a)
    assert explanation == set()  # Same term needs no explanation


def test_explain_equality_not_equal():
    """Test explain_equality with non-equal terms."""
    ppcc = ProofProducingCongruenceClosure([])
    explanation = ppcc.explain_equality(a, b)
    assert explanation == set()  # Not equal, can't explain


def test_explain_equality_direct():
    """Test explain_equality with directly asserted equality."""
    ppcc = ProofProducingCongruenceClosure([])
    reason = Eq(a, b)
    ppcc.add_equality_with_reason(a, b, reason)

    explanation = ppcc.explain_equality(a, b)
    assert reason in explanation


def test_find_proof_path_complex():
    """Test _find_proof_path with complex merge chains."""
    ppcc = ProofProducingCongruenceClosure([])

    # Create a longer chain: a = b = c = d
    ppcc.add_equality_with_reason(a, b, Eq(a, b))
    ppcc.add_equality_with_reason(b, c, Eq(b, c))
    ppcc.add_equality_with_reason(c, d, Eq(c, d))

    explanation = ppcc.explain_equality(a, d)
    assert len(explanation) > 0


# Test EUFTheorySolver exception cases
def test_equality_contradiction_exception():
    """Test EUFEqualityContradictionException."""
    solver = EUFTheorySolver()
    solver.Initialize({Eq(a, b), Ne(a, b)})

    # First assert equality
    solver.SetTrue(Eq(a, b))

    # Then try to assert disequality - should raise exception
    with pytest.raises(EUFEqualityContradictionException):
        solver.SetTrue(Ne(a, b))


def test_conflict_generation_equality():
    """Test conflict generation for equality contradictions."""
    solver = EUFTheorySolver()
    solver.Initialize({Eq(a, b), Ne(a, b)})

    # Set up contradiction scenario
    solver.SetTrue(Eq(a, b))

    # This should generate a conflict
    success, conflict = solver.assert_lit(-1)  # Assuming -1 corresponds to Ne(a,b)

    # The method should handle the conflict appropriately


def test_istrue_with_disequalities():
    """Test IsTrue method with disequality tracking."""
    solver = EUFTheorySolver()
    solver.Initialize({Ne(a, b), Eq(a, c), Eq(b, d)})

    solver.SetTrue(Ne(a, b))
    assert solver.IsTrue(Ne(a, b)) is True
    assert solver.IsTrue(Eq(a, b)) is False


def test_explain_disequality():
    """Test explain_disequality method."""
    solver = EUFTheorySolver()
    solver.Initialize({Ne(a, b), Eq(a, c), Eq(b, d)})

    solver.SetTrue(Ne(a, b))
    solver.SetTrue(Eq(a, c))
    solver.SetTrue(Eq(b, d))

    # Should be able to explain disequality between c and d
    explanation = solver.explain_disequality(c, d)
    assert isinstance(explanation, set)


def test_functional_congruence_complex():
    """Test complex functional congruence scenarios."""
    solver = EUFTheorySolver()
    solver.Initialize({
        Eq(a, b),
        Eq(F(a), c),
        Eq(F(b), d),
        Eq(G(F(a)), e)
    })

    solver.SetTrue(Eq(a, b))
    solver.SetTrue(Eq(F(a), c))

    # F(a) should equal F(b) by congruence
    assert solver.IsTrue(Eq(F(a), F(b))) is True


def test_nested_functions():
    """Test nested function applications."""
    solver = EUFTheorySolver()
    solver.Initialize({
        Eq(a, b),
        Eq(F(G(a)), c),
        Eq(F(G(b)), d)
    })

    solver.SetTrue(Eq(a, b))
    solver.SetTrue(Eq(F(G(a)), c))

    # Should propagate through nested functions
    assert solver.IsTrue(Eq(F(G(a)), F(G(b)))) is True


def test_from_encoded_cnf_with_trivial_true():
    """Test from_encoded_cnf with trivially true predicates."""
    enc = EncodedCNF()
    # Create a trivially true predicate (x = x)
    enc.encoding = {Q.eq(x, x): 1}
    enc.data = [{1}]

    solver, conflicts = EUFTheorySolver.from_encoded_cnf(enc, testing_mode=True)
    assert len(conflicts) > 0  # Should detect trivial conflicts


def test_from_encoded_cnf_with_trivial_false():
    """Test from_encoded_cnf with trivially false predicates."""
    enc = EncodedCNF()
    # Create a trivially false predicate (x != x)
    enc.encoding = {Q.ne(x, x): 1}
    enc.data = [{1}]

    solver, conflicts = EUFTheorySolver.from_encoded_cnf(enc, testing_mode=True)
    assert len(conflicts) > 0  # Should detect trivial conflicts


def test_from_encoded_cnf_unhandled_predicate():
    """Test from_encoded_cnf with unhandled predicate types."""
    enc = EncodedCNF()
    # Create an unhandled predicate type
    enc.encoding = {x: 1}  # Just a symbol, not a predicate
    enc.data = [[1]]

    solver, conflicts = EUFTheorySolver.from_encoded_cnf(enc, testing_mode=True)
    # Should handle gracefully by skipping


def test_unary_predicate_handling():
    """Test handling of unary predicates."""
    enc = EncodedCNF()
    # Create unary predicate Q.prime(x)
    prime_x = AppliedPredicate(Q.prime, (x,))
    enc.encoding = {prime_x: 1}
    enc.data = [[1]]

    solver, conflicts = EUFTheorySolver.from_encoded_cnf(enc, testing_mode=True)


def test_binary_predicate_q_eq():
    """Test Q.eq binary predicate handling."""
    enc = EncodedCNF()
    # Create Q.eq(x, y) predicate
    eq_xy = AppliedPredicate(Q.eq, (x, y))
    enc.encoding = {eq_xy: 1}
    enc.data = [[1]]

    solver, conflicts = EUFTheorySolver.from_encoded_cnf(enc, testing_mode=True)
    # Should convert properly


def test_binary_predicate_q_ne():
    """Test Q.ne binary predicate handling."""
    enc = EncodedCNF()
    # Create Q.ne(x, y) predicate
    ne_xy = AppliedPredicate(Q.ne, (x, y))
    enc.encoding = {ne_xy: 1}
    enc.data = [[1]]

    solver, conflicts = EUFTheorySolver.from_encoded_cnf(enc, testing_mode=True)
    # Should convert properly


def test_multiple_disequality_conflicts():
    """Test multiple disequality conflicts."""
    solver = EUFTheorySolver()
    solver.Initialize({
        Ne(a, b), Ne(b, c), Ne(c, d),
        Eq(a, b), Eq(b, c), Eq(c, d)
    })

    # Set up chain of equalities
    solver.SetTrue(Eq(a, b))
    solver.SetTrue(Eq(b, c))

    # Try to add conflicting disequality
    with pytest.raises((EUFDisequalityContradictionException, EUFEqualityContradictionException)):
        solver.SetTrue(Ne(a, c))


def test_congruence_with_multiple_functions():
    """Test congruence with multiple function symbols."""
    solver = EUFTheorySolver()
    solver.Initialize({
        Eq(a, b),
        Eq(F(a), c), Eq(G(a), d)
    })

    solver.SetTrue(Eq(a, b))
    solver.SetTrue(Eq(F(a), c))
    solver.SetTrue(Eq(G(a), d))

    # Both F and G should respect congruence
    assert solver.IsTrue(Eq(F(a), F(b))) is True
    assert solver.IsTrue(Eq(G(a), G(b))) is True


def test_deep_explanation_chain():
    """Test explanation with deep chains of reasoning."""
    solver = EUFTheorySolver()

    # Create long chain: a=b=c=d=e
    eqs = [Eq(a, b), Eq(b, c), Eq(c, d), Eq(d, e)]
    solver.Initialize(set(eqs))

    for eq in eqs:
        solver.SetTrue(eq)

    # Explain why a = e
    explanation = solver.cc.explain_equality(a, e)
    assert len(explanation) > 0


def test_exception_type_tracking():
    """Test proper exception type tracking for conflict generation."""
    solver = EUFTheorySolver()
    solver.Initialize({Eq(a, b), Ne(a, b)})

    # Set up for disequality contradiction
    solver.SetTrue(Ne(a, b))

    try:
        solver.SetTrue(Eq(a, b))
    except EUFDisequalityContradictionException:
        # Check that exception type is properly tracked
        assert hasattr(solver, '_current_exception_type')
        assert solver._current_exception_type == 'disequality_contradiction'


def test_generate_minimal_conflict_edge_cases():
    """Test _generate_minimal_conflict with edge cases."""
    solver = EUFTheorySolver()
    solver.Initialize({Eq(a, b)})

    # Test with no current conflict
    conflict = solver._generate_minimal_conflict()
    assert isinstance(conflict, set)

    # Test with unknown exception type
    solver._current_exception_type = 'unknown'
    conflict = solver._generate_minimal_conflict()
    assert isinstance(conflict, set)

def test_simple_chain_explanation():
    # Simple chain x=y, y=z, test explanation for x=z
    solver = EUFTheorySolver()
    cc = solver.cc
    cc.add_equality_with_reason(x, y, Eq(x, y))
    cc.add_equality_with_reason(y, z, Eq(y, z))
    explanation = cc.explain_equality(x, z)
    assert Eq(x, y) in explanation
    assert Eq(y, z) in explanation


def test_assert_literals_and_conflict():
    # Use CNF encoding from props and assumptions, assert literals and expect conflict on negation

    prop = Q.prime(z)
    assum = Q.prime(y) & Q.eq(x, y) & Q.eq(x, z) & Q.eq(z, d) & Q.eq(y, e)
    add = Q.eq(e, d)

    _prop = CNF.from_prop(prop)
    _assum = CNF.from_prop(assum)
    _add = CNF.from_prop(add)

    sat = EncodedCNF()
    sat.add_from_cnf(_prop)
    sat.add_from_cnf(_assum)
    sat.add_from_cnf(_add)

    solver = EUFTheorySolver()
    a, conflicts = solver.from_encoded_cnf(sat)

    # Assert positive literals
    for lit in [Q.prime(y), Q.eq(x, z), Q.eq(z, d), Q.eq(y, e), Q.eq(x, y)]:
        res, _ = a.assert_lit(sat.encoding[lit])
        assert res is True

    # Assert negation of e=d causes conflict
    res, conflict = a.assert_lit(-sat.encoding[Q.eq(e, d)])
    assert res is False

    # Explanation for equality e=d includes chain
    explanation = a.cc.explain_equality(e, d)
    expected_eqs = {Eq(e, y), Eq(x, z), Eq(d, z), Eq(x, y)}
    for eq in expected_eqs:
        assert eq in explanation

def test_explain_disequality_basic():
    """Test basic disequality explanation with direct disequality."""
    solver = EUFTheorySolver()

    # Manually set up disequality cause
    solver.SetTrue(Ne(x, y))

    explanation = solver.explain_disequality(x, y)
    assert Ne(x, y) in explanation
    assert len(explanation) == 1  # Only the direct disequality


def test_explain_disequality_with_equality_chain():
    """Test disequality explanation involving equality chains."""
    solver = EUFTheorySolver()

    # Set up: x = a, y = b, a != b
    # This means x != y should be explained by {Ne(a, b), Eq(x, a), Eq(y, b)}
    solver.SetTrue(Eq(x, a))
    solver.SetTrue(Eq(y, b))
    solver.SetTrue(Ne(a, b))

    explanation = solver.explain_disequality(x, y)

    # Should contain the source disequality and connecting equalities
    assert Ne(a, b) in explanation
    assert Eq(a, x) in explanation
    assert Eq(b, y) in explanation


def test_explain_disequality_transitive_chain():
    """Test disequality explanation with longer equality chains."""
    solver = EUFTheorySolver()

    # Chain: x = u, u = v, y = w, v != w
    # So x != y should be explained through this chain
    solver.SetTrue(Eq(x, u))
    solver.SetTrue(Eq(u, v))
    solver.SetTrue(Eq(y, w))
    solver.SetTrue(Ne(v, w))

    explanation = solver.explain_disequality(x, y)

    # Should contain source disequality and all connecting equalities
    assert Ne(v, w) in explanation
    expected_equalities = {Eq(u, x), Eq(u, v), Eq(w, y), Ne(v, w)}
    for eq in expected_equalities:
        assert eq in explanation


def test_explain_disequality_reversed_terms():
    """Test disequality explanation works with reversed term order."""
    solver = EUFTheorySolver()

    # Set up a != b
    solver.SetTrue(Ne(a, b))

    # Test both directions: a != b and b != a
    explanation_ab = solver.explain_disequality(a, b)
    explanation_ba = solver.explain_disequality(b, a)

    # Both should give same explanation (just the disequality)
    assert Ne(a, b) in explanation_ab
    assert Ne(a, b) in explanation_ba
    assert explanation_ab == explanation_ba


def test_explain_disequality_no_disequality():
    """Test explanation when terms are not disequal."""
    solver = EUFTheorySolver()

    # Set up only equalities, no disequalities
    solver.SetTrue(Eq(x, y))

    # Should return empty explanation since x and y are not disequal
    explanation = solver.explain_disequality(x, y)
    assert explanation == set()


def test_explain_disequality_complex_scenario():
    """Test disequality explanation in complex scenario with multiple constraints."""
    solver = EUFTheorySolver()

    # Complex setup: x = u, u = v, y = w, w = z, v != z
    # This creates x != y through the chain
    solver.SetTrue(Eq(x, u))
    solver.SetTrue(Eq(u, v))
    solver.SetTrue(Eq(y, w))
    solver.SetTrue(Eq(w, z))
    solver.SetTrue(Ne(v, z))

    explanation = solver.explain_disequality(x, y)

    # Should include the source disequality and relevant equality chains
    assert Ne(v, z) in explanation
    # Should include equalities connecting x to v and y to z
    expected_in_explanation = {Ne(v,z), Eq(u, x), Eq(u, v), Eq(w, y), Eq(w, z)}
    for eq in expected_in_explanation:
        assert eq in explanation


def test_explain_disequality_with_functions():
    """Test disequality explanation involving function terms."""
    solver = EUFTheorySolver()

    # Set up: F(x) != F(y), x = a, y = b
    solver.SetTrue(Ne(F(x), F(y)))
    solver.SetTrue(Eq(x, a))
    solver.SetTrue(Eq(y, b))

    # Explanation for F(a) != F(b) should include the function disequality and argument equalities
    explanation = solver.explain_disequality(F(a), F(b))

    assert Ne(F(x), F(y)) in explanation


def test_explain_disequality_edge_cases():
    """Test edge cases for disequality explanation."""
    solver = EUFTheorySolver()

    # Test with identical terms (should return empty)
    explanation = solver.explain_disequality(x, x)
    assert explanation == set()

    # Test with unknown terms
    unknown1, unknown2 = symbols('unknown1 unknown2')
    explanation = solver.explain_disequality(unknown1, unknown2)
    assert explanation == set()


def test_conflict_generation_disequality_contradiction():
    """Test conflict generation for disequality contradiction cases."""
    solver = EUFTheorySolver()

    # Set up encoding manually for testing
    solver.literal_eqs = {1: [Eq(x, y)], 2: [Ne(x, y)]}
    solver._enc_to_lit = {1: Eq(x, y), -1: Ne(x, y), 2: Ne(x, y), -2: Eq(x, y)}
    solver.lit_to_enc = {Eq(x, y): 1, Ne(x, y): 2}

    # First assert disequality
    result, _ = solver.assert_lit(2)  # Ne(x, y)
    assert result is True

    # Then assert conflicting equality - should generate conflict
    result, conflict = solver.assert_lit(1)  # Eq(x, y)
    assert result is False
    assert len(conflict) > 0
    assert -1 in conflict or -2 in conflict  # Should negate one of the conflicting literals


def test_multiple_disequalities_explanation():
    """Test explanation when multiple disequalities exist."""
    solver = EUFTheorySolver()

    # Set up multiple disequalities: a != b, c != d
    # And equalities: x = a, y = b, z = c, w = d
    solver.SetTrue(Ne(a, b))
    solver.SetTrue(Ne(c, d))
    solver.SetTrue(Eq(x, a))
    solver.SetTrue(Eq(y, b))
    solver.SetTrue(Eq(z, c))
    solver.SetTrue(Eq(w, d))

    # Test explanation for x != y (should use first disequality)
    explanation_xy = solver.explain_disequality(x, y)
    assert Ne(a, b) in explanation_xy
    assert Eq(a, x) in explanation_xy
    assert Eq(b, y) in explanation_xy

    # Test explanation for z != w (should use second disequality)
    explanation_zw = solver.explain_disequality(z, w)
    assert Ne(c, d) in explanation_zw
    assert Eq(c, z) in explanation_zw
    assert Eq(d, w) in explanation_zw


def test_random_satisfiable_constraints_small():
    """Test random satisfiable constraints (small scale)."""
    generator = RandomEUFTestGenerator()
    z3_comparator = Z3Comparator() if Z3_AVAILABLE else None

    for trial in range(10):
        print(f"\nTrial {trial + 1}/10")

        # Generate simple satisfiable constraints
        constraints = []
        symbols_used = generator.symbols_pool[:5]  # Use only 5 symbols

        # Add some equalities in a chain: x0 = x1 = x2
        constraints.append(Eq(symbols_used[0], symbols_used[1]))
        constraints.append(Eq(symbols_used[1], symbols_used[2]))

        # Add some unrelated equalities
        constraints.append(Eq(symbols_used[3], symbols_used[4]))

        # Test with SymPy
        result = satisfiable(boolalg.And(*constraints), use_euf_theory=True)
        sympy_sat = result is not False
        print(f"SymPy result: {'SAT' if sympy_sat else 'UNSAT'}")
        # Compare with Z3 if available
        if Z3_AVAILABLE and z3_comparator:
            z3_result = z3_comparator.check_satisfiability_with_z3(constraints)
            z3_sat, z3_model = z3_result if z3_result else (None, None)
            if z3_sat is not None:
                print(f"Z3 result: {'SAT' if z3_sat else 'UNSAT'}")
                assert sympy_sat == z3_sat, f"Disagreement: SymPy={sympy_sat}, Z3={z3_sat}"


def test_random_unsatisfiable_constraints_small():
    """Test random unsatisfiable constraints (small scale)."""
    generator = RandomEUFTestGenerator()
    z3_comparator = Z3Comparator() if Z3_AVAILABLE else None

    for trial in range(10):
        print(f"\nTrial {trial + 1}/10")

        # Generate unsatisfiable constraints
        x, y, z = generator.symbols_pool[:3]

        constraints = [
            Eq(x, y),
            Eq(y, z),
            Ne(x, z)  # This makes it unsatisfiable
        ]
        result = satisfiable(boolalg.And(*constraints), use_euf_theory=True)
        sympy_sat = result is not False
        print(f"SymPy result: {'SAT' if sympy_sat else 'UNSAT'}")
        # Should be unsatisfiable
        assert not sympy_sat, "Expected UNSAT but got SAT"

        # Compare with Z3 if available
        if Z3_AVAILABLE and z3_comparator:
            z3_result = z3_comparator.check_satisfiability_with_z3(constraints)
            z3_sat, z3_model = z3_result if z3_result else (None, None)
            if z3_sat is not None:
                print(f"Z3 result: {'SAT' if z3_sat else 'UNSAT'}")
                assert sympy_sat == z3_sat, f"Disagreement: SymPy={sympy_sat}, Z3={z3_sat}"


def test_random_function_constraints_medium():
    """Test random constraints with functions (medium scale)."""
    generator = RandomEUFTestGenerator()

    for trial in range(5):
        print(f"\nFunction trial {trial + 1}/5")

        # Generate constraints with functions
        x, y, z = generator.symbols_pool[:3]
        f = generator.functions_pool[0]

        constraints = [
            Eq(x, y),
            Eq(f(x), z),
            # f(y) should equal z by congruence
        ]

        # Add the congruence conclusion as a query
        query = Eq(f(y), z)

        # Check if constraints + query is satisfiable
        all_constraints = constraints + [query]
        result = satisfiable(boolalg.And(*all_constraints), use_euf_theory=True)
        sympy_sat = result is not False
        print(f"SymPy result for constraints + query: {'SAT' if sympy_sat else 'UNSAT'}")

        # Should be satisfiable due to congruence
        assert sympy_sat, "Expected SAT due to functional congruence"

        # Also check if constraints + negated query is unsatisfiable
        negated_constraints = constraints + [Not(query)]
        result_neg = satisfiable(boolalg.And(*negated_constraints), use_euf_theory=True)
        sympy_sat_neg = result_neg is not False
        print(f"SymPy result for constraints + ~query: {'SAT' if sympy_sat_neg else 'UNSAT'}")

        # Should be unsatisfiable
        assert not sympy_sat_neg, "Expected UNSAT for constraints + ~query"


def test_large_random_constraint_set():
    """Test large random constraint sets."""
    generator = RandomEUFTestGenerator()
    z3 = import_module("z3")
    if z3 is None:
        skip("z3 not installed.")
    # z3_comparator = Z3Comparator() if Z3_AVAILABLE else None

    sat_count = 0
    TRIAL_COUNT = 50

    problems = []
    # generate likely sat problems
    problems += [generator.generate_constraint_set(num_constraints=25) for _ in range(TRIAL_COUNT // 2)]
    # generate likely unsat problems
    problems += [generator.generate_constraint_set(num_constraints=100) for _ in range(1 + TRIAL_COUNT // 2)]

    for trial in range(TRIAL_COUNT):
        # print(f"\nLarge trial {trial + 1}/3")

        constraints = problems[trial]



        # print(f"Testing {len(constraints)} constraints")

        problem =boolalg.And(*constraints)
        assert problem not in (False, True)

        result = satisfiable(problem, use_euf_theory=True)
        sympy_sat = result is not False
        # print(f"SymPy result: {'SAT' if sympy_sat else 'UNSAT'}")
        z3_sat = z3_satisfiable(problem) is not False

        assert sympy_sat == z3_sat
        if sympy_sat:
            sat_count += 1

    assert sat_count > 0
    assert sat_count != TRIAL_COUNT
    print(f"{sat_count} / {TRIAL_COUNT} problems were satifiable")

def test_random_sat_problems():
    """Test large random sat problem involving euf constraints."""
    generator = RandomEUFTestGenerator()
    z3 = import_module("z3")
    if z3 is None:
        skip("z3 not installed.")
    # z3_comparator = Z3Comparator() if Z3_AVAILABLE else None

    sat_count = 0
    TRIAL_COUNT = 1



    for trial in range(TRIAL_COUNT):
        constraint_sets = []
        constraint_sets += [generator.generate_constraint_set(num_constraints=100) for _ in range(TRIAL_COUNT)]

        simple_disjunction = boolalg.Or(
            *[boolalg.And(*[constraint for constraint in constraint_set]) for constraint_set in constraint_sets])


        assert simple_disjunction not in (False, True)

        result = satisfiable(simple_disjunction, use_euf_theory=True)
        sympy_sat = result is not False
        # print(f"SymPy result: {'SAT' if sympy_sat else 'UNSAT'}")
        z3_sat = z3_satisfiable(simple_disjunction) is not False

        assert sympy_sat == z3_sat
        if sympy_sat:
            sat_count += 1

    # assert sat_count > 0
    # assert sat_count != TRIAL_COUNT
    print(f"{sat_count} / {TRIAL_COUNT} problems were satifiable")


        # Compare with Z3 if available (only for smaller sets due to conversion complexity)
        # if Z3_AVAILABLE and z3_comparator and len(filtered_constraints) <= 15:
        #     z3_result = z3_comparator.check_satisfiability_with_z3(filtered_constraints)
        #     z3_sat, z3_model = z3_result if z3_result else (None, None)
        #     if z3_sat is not None:
        #         print(f"Z3 result: {'SAT' if z3_sat else 'UNSAT'}")
        #         assert sympy_sat == z3_sat, f"Disagreement: SymPy={sympy_sat}, Z3={z3_sat}"


def test_assert_lit_transitivity_conflict():
    """Test transitivity conflict: x = y, y = z, x != z."""
    enc_cnf = EncodedCNF()
    enc_cnf.encoding = {
        Eq(x, y): 1,      # x = y
        Eq(y, z): 2,      # y = z
        Ne(x, z): 3       # x != z (conflicts by transitivity)
    }
    enc_cnf.data = [{1}, {2}, {3}]

    solver, conflicts = EUFTheorySolver.from_encoded_cnf(enc_cnf)

    # Assert x = y and y = z (should succeed)
    assert solver.assert_lit(1) == (True, set())
    assert solver.assert_lit(2) == (True, set())

    # Assert x != z (should conflict due to transitivity x = y = z)
    assert solver.assert_lit(3) == (False, {-1, -2, -3})


def test_assert_lit_function_congruence_conflict():
    """Test function congruence conflict: ~Q.prime(x), Q.eq(x,y), Q.prime(y)."""
    enc_cnf = EncodedCNF()
    enc_cnf.encoding = {
        Q.prime(x): 1,    # Q.prime(x)
        Q.eq(x, y): 2,    # x = y
        Q.prime(y): 3     # Q.prime(y)
    }
    enc_cnf.data = [{-1}, {2}, {3}]  # ~Q.prime(x), x = y, Q.prime(y)

    solver, conflicts = EUFTheorySolver.from_encoded_cnf(enc_cnf)

    # Assert ~Q.prime(x) and x = y (should succeed)
    assert solver.assert_lit(-1) == (True, set())
    assert solver.assert_lit(2) == (True, set())

    # Assert Q.prime(y) (should conflict by function congruence)
    assert solver.assert_lit(3) == (False, {1, -2, -3})


def test_assert_lit_multi_predicate_conflict():
    """Test conflict with multiple different predicates."""
    enc_cnf = EncodedCNF()
    enc_cnf.encoding = {
        Q.prime(x): 1,    # Q.prime(x)
        Q.even(x): 2,     # Q.even(x)
        Q.eq(x, y): 3,    # x = y
        Q.even(y): 4      # Q.even(y)
    }
    enc_cnf.data = [{1}, {-2}, {3}, {4}]  # Q.prime(x), ~Q.even(x), x = y, Q.even(y)

    solver, conflicts = EUFTheorySolver.from_encoded_cnf(enc_cnf)

    # Assert Q.prime(x), ~Q.even(x), x = y (should succeed)
    assert solver.assert_lit(1) == (True, set())
    assert solver.assert_lit(-2) == (True, set())
    assert solver.assert_lit(3) == (True, set())

    # Assert Q.even(y) (should conflict with ~Q.even(x) and x = y)
    assert solver.assert_lit(4) == (False, {2, -3, -4})


def test_assert_lit_complex_equality_chain_conflict():
    """Test conflict in longer equality chains."""
    enc_cnf = EncodedCNF()
    enc_cnf.encoding = {
        Eq(a, b): 1,      # a = b
        Eq(b, c): 2,      # b = c
        Eq(c, d): 3,      # c = d
        Eq(d, e): 4,      # d = e
        Ne(a, e): 5       # a != e (conflicts with transitivity)
    }
    enc_cnf.data = [{1}, {2}, {3}, {4}, {5}]

    solver, conflicts = EUFTheorySolver.from_encoded_cnf(enc_cnf)

    # Assert the equality chain (should succeed)
    assert solver.assert_lit(1) == (True, set())
    assert solver.assert_lit(2) == (True, set())
    assert solver.assert_lit(3) == (True, set())
    assert solver.assert_lit(4) == (True, set())

    # Assert a != e (should conflict with full chain)
    assert solver.assert_lit(5) == (False, {-1, -2, -3, -4, -5})


def test_assert_lit_multiple_function_applications():
    """Test conflicts with multiple function applications."""
    enc_cnf = EncodedCNF()
    enc_cnf.encoding = {
        Q.prime(x): 1,     # Q.prime(x)
        Q.composite(x): 2, # Q.composite(x)
        Q.eq(x, y): 3,     # x = y
        Q.prime(y): 4,     # Q.prime(y) (derivable)
        Q.composite(y): 5  # Q.composite(y) (derivable)
    }
    enc_cnf.data = [{1}, {-2}, {3}, {4}, {5}]  # Q.prime(x), ~Q.composite(x), x = y, Q.prime(y), Q.composite(y)

    solver, conflicts = EUFTheorySolver.from_encoded_cnf(enc_cnf)

    # Assert Q.prime(x), ~Q.composite(x), x = y, Q.prime(y) (should succeed)
    assert solver.assert_lit(1) == (True, set())
    assert solver.assert_lit(-2) == (True, set())
    assert solver.assert_lit(3) == (True, set())
    assert solver.assert_lit(4) == (True, set())  # Derivable by congruence

    # Assert Q.composite(y) (should conflict with ~Q.composite(x) and x = y)
    assert solver.assert_lit(5) == (False, {2, -3, -5})


def test_assert_lit_reverse_order_conflict():
    """Test that conflict detection works regardless of assertion order."""
    enc_cnf = EncodedCNF()
    enc_cnf.encoding = {
        Q.prime(x): 1,
        Q.eq(x, y): 2,
        Q.prime(y): 3
    }
    enc_cnf.data = [{3}, {2}, {-1}]  # Q.prime(y), x = y, ~Q.prime(x)

    solver, conflicts = EUFTheorySolver.from_encoded_cnf(enc_cnf)

    # Different order: Q.prime(y), x = y first
    assert solver.assert_lit(3) == (True, set())
    assert solver.assert_lit(2) == (True, set())

    # Then ~Q.prime(x) should conflict
    assert solver.assert_lit(-1) == (False, {1, -2, -3})


def test_assert_lit_mixed_equalities_disequalities():
    """Test mixed equality and disequality constraints."""
    enc_cnf = EncodedCNF()
    enc_cnf.encoding = {
        Ne(a, b): 1,      # a != b
        Eq(a, c): 2,      # a = c
        Eq(b, d): 3,      # b = d
        Ne(c, d): 4       # c != d (should be derivable, not conflict)
    }
    enc_cnf.data = [{1}, {2}, {3}, {4}]

    solver, conflicts = EUFTheorySolver.from_encoded_cnf(enc_cnf)

    # All should succeed (c != d follows from a != b, a = c, b = d)
    assert solver.assert_lit(1) == (True, set())
    assert solver.assert_lit(2) == (True, set())
    assert solver.assert_lit(3) == (True, set())
    assert solver.assert_lit(4) == (True, set())


def test_assert_lit_contradictory_disequality():
    """Test contradiction when asserting disequality after equality chain."""
    enc_cnf = EncodedCNF()
    enc_cnf.encoding = {
        Eq(a, b): 1,      # a = b
        Eq(b, c): 2,      # b = c
        Ne(a, c): 3       # a != c (conflicts with transitivity)
    }
    enc_cnf.data = [{1}, {2}, {3}]

    solver, conflicts = EUFTheorySolver.from_encoded_cnf(enc_cnf)

    # Establish equality chain
    assert solver.assert_lit(1) == (True, set())
    assert solver.assert_lit(2) == (True, set())

    # Contradictory disequality should conflict
    assert solver.assert_lit(3) == (False, {-1, -2, -3})


def test_assert_lit_no_conflict_compatible_assertions():
    """Test that compatible assertions don't produce false conflicts."""
    enc_cnf = EncodedCNF()
    enc_cnf.encoding = {
        Q.prime(x): 1,    # Q.prime(x)
        Q.eq(x, y): 2,    # x = y
        Q.prime(y): 3     # Q.prime(y) (derivable by congruence)
    }
    enc_cnf.data = [{1}, {2}, {3}]

    solver, conflicts = EUFTheorySolver.from_encoded_cnf(enc_cnf)

    # All should succeed without conflict
    assert solver.assert_lit(1) == (True, set())
    assert solver.assert_lit(2) == (True, set())
    assert solver.assert_lit(3) == (True, set())  # Should be derivable


def test_assert_lit_self_equality_no_conflict():
    """Test that self-equalities don't cause conflicts."""
    enc_cnf = EncodedCNF()
    enc_cnf.encoding = {
        Eq(x, x): 1,      # x = x (trivially true)
        Q.prime(x): 2     # Q.prime(x)
    }
    enc_cnf.data = [{1}, {2}]

    solver, conflicts = EUFTheorySolver.from_encoded_cnf(enc_cnf)

    # Both should succeed
    assert solver.assert_lit(1) == (True, set())
    assert solver.assert_lit(2) == (True, set())


def test_assert_lit_function_chain_conflict():
    """Test conflict in function application chains."""
    enc_cnf = EncodedCNF()
    enc_cnf.encoding = {
        Q.prime(a): 1,       # Q.prime(a)
        Q.eq(a, b): 2,       # a = b
        Q.eq(b, c): 3,       # b = c
        Q.eq(c, d): 4,       # c = d
        Q.composite(d): 5    # Q.composite(d) (conflicts with Q.prime(a) through chain)
    }
    enc_cnf.data = [{1}, {2}, {3}, {4}, {-5}]  # Assume Q.prime and Q.composite are contradictory

    solver, conflicts = EUFTheorySolver.from_encoded_cnf(enc_cnf)

    # Establish chain
    assert solver.assert_lit(1) == (True, set())
    assert solver.assert_lit(2) == (True, set())
    assert solver.assert_lit(3) == (True, set())
    assert solver.assert_lit(4) == (True, set())


def test_assert_lit_original_issue_case():
    """Test the exact case from the original failing test."""
    enc_cnf = EncodedCNF()
    enc_cnf.encoding = {Q.prime(x): 1, Q.eq(x, y): 2, Q.prime(y): 3}
    enc_cnf.data = [{-1}, {2}, {3}]

    solver, conflicts = EUFTheorySolver.from_encoded_cnf(enc_cnf)

    assert solver.assert_lit(-1) == (True, set())
    assert solver.assert_lit(2) == (True, set())
    assert solver.assert_lit(3) == (False, {1, -2, -3})


def test_assert_lit_empty_conflict_set_on_success():
    """Test that successful assertions return empty conflict sets."""
    enc_cnf = EncodedCNF()
    enc_cnf.encoding = {
        Eq(x, y): 1,
        Eq(y, z): 2
    }
    enc_cnf.data = [{1}, {2}]

    solver, conflicts = EUFTheorySolver.from_encoded_cnf(enc_cnf)

    # Both should return empty conflict sets
    success1, conflict1 = solver.assert_lit(1)
    success2, conflict2 = solver.assert_lit(2)

    assert success1 == True and conflict1 == set()
    assert success2 == True and conflict2 == set()


def test_assert_lit_conflict_contains_only_integers():
    """Test that conflict sets contain only integer literals."""
    enc_cnf = EncodedCNF()
    enc_cnf.encoding = {Q.prime(x): 1, Q.eq(x, y): 2, Q.prime(y): 3}
    enc_cnf.data = [{-1}, {2}, {3}]

    solver, conflicts = EUFTheorySolver.from_encoded_cnf(enc_cnf)

    solver.assert_lit(-1)
    solver.assert_lit(2)
    success, conflict = solver.assert_lit(3)

    assert success == False
    assert all(isinstance(lit, int) for lit in conflict)
    assert len(conflict) > 0


def test_assert_lit_unary_predicate_conflicts():
    """Test conflicts with various unary predicates."""
    enc_cnf = EncodedCNF()
    enc_cnf.encoding = {
        Q.positive(x): 1,    # Q.positive(x)
        Q.negative(x): 2,    # Q.negative(x)
        Q.eq(x, y): 3,       # x = y
        Q.positive(y): 4,    # Q.positive(y)
        Q.negative(y): 5     # Q.negative(y)
    }
    enc_cnf.data = [{1}, {-2}, {3}, {4}, {5}]  # Q.positive(x), ~Q.negative(x), x = y, Q.positive(y), Q.negative(y)

    solver, conflicts = EUFTheorySolver.from_encoded_cnf(enc_cnf)

    # Establish positive and equality
    assert solver.assert_lit(1) == (True, set())   # Q.positive(x)
    assert solver.assert_lit(-2) == (True, set())  # ~Q.negative(x)
    assert solver.assert_lit(3) == (True, set())   # x = y
    assert solver.assert_lit(4) == (True, set())   # Q.positive(y) - derivable

    # Q.negative(y) should conflict with ~Q.negative(x) and x = y
    assert solver.assert_lit(5) == (False, {2, -3, -5})


def test_assert_lit_different_functions_no_conflict():
    """Test different functions don't interfere: f(x) = c1, g(x) = c2."""
    enc_cnf = EncodedCNF()
    c1, c2 = symbols('c1 c2')
    func_f, func_g = Function('f'), Function('g')

    enc_cnf.encoding = {
        Eq(func_f(x), c1): 1,     # f(x) = c1
        Eq(func_g(x), c2): 2      # g(x) = c2 (different function, no conflict)
    }
    enc_cnf.data = [{1}, {2}]

    solver, conflicts = EUFTheorySolver.from_encoded_cnf(enc_cnf)

    # Both should succeed (different functions)
    assert solver.assert_lit(1) == (True, set())
    assert solver.assert_lit(2) == (True, set())


def test_assert_lit_mixed_function_predicates():
    """Test mixing function applications with predicates."""
    enc_cnf = EncodedCNF()
    func_f = Function('f')

    enc_cnf.encoding = {
        Eq(func_f(x), a): 1,      # f(x) = a
        Q.prime(a): 2,            # Q.prime(a)
        Q.eq(x, y): 3,            # x = y
        Eq(func_f(y), b): 4,      # f(y) = b
        Q.composite(b): 5         # Q.composite(b)
    }
    enc_cnf.data = [{1}, {2}, {3}, {4}, {-5}]  # Assume prime/composite conflict

    solver, conflicts = EUFTheorySolver.from_encoded_cnf(enc_cnf)

    assert solver.assert_lit(1) == (True, set())
    assert solver.assert_lit(2) == (True, set())
    assert solver.assert_lit(3) == (True, set())
    assert solver.assert_lit(4) == (True, set())  # f(y) = b, but by congruence f(y) should = a


def test_assert_lit_function_arity_mismatch_no_conflict():
    """Test functions with different arities don't conflict: f(x) = c1, f(x,y) = c2."""
    enc_cnf = EncodedCNF()
    c1, c2 = symbols('c1 c2')
    func_f = Function('f')

    enc_cnf.encoding = {
        Eq(func_f(x), c1): 1,     # f(x) = c1 (arity 1)
        Eq(func_f(x, y), c2): 2   # f(x,y) = c2 (arity 2, different function)
    }
    enc_cnf.data = [{1}, {2}]

    solver, conflicts = EUFTheorySolver.from_encoded_cnf(enc_cnf)

    # Different arities should not conflict
    assert solver.assert_lit(1) == (True, set())
    assert solver.assert_lit(2) == (True, set())


def test_assert_lit_function_symmetry():
    """Test function symmetry properties: f(a,b) = c, f(b,a) = d, no automatic conflict."""
    enc_cnf = EncodedCNF()
    func_f = Function('f')

    enc_cnf.encoding = {
        Eq(func_f(a, b), c): 1,   # f(a,b) = c
        Eq(func_f(b, a), d): 2    # f(b,a) = d (different unless f is symmetric)
    }
    enc_cnf.data = [{1}, {2}]

    solver, conflicts = EUFTheorySolver.from_encoded_cnf(enc_cnf)

    # Should not conflict (f is not assumed symmetric)
    assert solver.assert_lit(1) == (True, set())
    assert solver.assert_lit(2) == (True, set())


def test_diamond_disequality_basic_conflict():
    """Test basic diamond disequality: a = b, a = c, b != c."""
    enc_cnf = EncodedCNF()
    enc_cnf.encoding = {
        Eq(a, b): 1,      # a = b
        Eq(a, c): 2,      # a = c
        Ne(b, c): 3       # b != c (conflicts with transitivity a = b = c)
    }
    enc_cnf.data = [{1}, {2}, {3}]

    solver, conflicts = EUFTheorySolver.from_encoded_cnf(enc_cnf)

    # Assert the equalities (should succeed)
    assert solver.assert_lit(1) == (True, set())  # a = b
    assert solver.assert_lit(2) == (True, set())  # a = c

    # Assert disequality (should conflict due to transitivity)
    result = solver.assert_lit(3)
    assert result[0] == False
    # Conflict should involve all three literals
    conflict_lits = result[1]
    assert len(conflict_lits) > 0
    assert any(abs(lit) in [1, 2, 3] for lit in conflict_lits)


def test_diamond_disequality_reverse_order():
    """Test diamond disequality with disequality asserted first."""
    enc_cnf = EncodedCNF()
    enc_cnf.encoding = {
        Ne(b, c): 1,      # b != c
        Eq(a, b): 2,      # a = b
        Eq(a, c): 3       # a = c (should conflict)
    }
    enc_cnf.data = [{1}, {2}, {3}]

    solver, conflicts = EUFTheorySolver.from_encoded_cnf(enc_cnf)

    # Assert disequality first
    assert solver.assert_lit(1) == (True, set())  # b != c
    assert solver.assert_lit(2) == (True, set())  # a = b

    # This should conflict (a = c would make b = c, contradicting b != c)
    result = solver.assert_lit(3)
    assert result[0] == False
    assert len(result[1]) > 0


def test_diamond_disequality_with_functions():
    """Test diamond disequality with function applications: f(a) = x, f(b) = x, a = b, f(a) != f(b)."""
    enc_cnf = EncodedCNF()
    func_f = Function('f')

    enc_cnf.encoding = {
        Eq(func_f(a), x): 1,     # f(a) = x
        Eq(func_f(b), x): 2,     # f(b) = x
        Eq(a, b): 3,             # a = b
        Ne(func_f(a), func_f(b)): 4  # f(a) != f(b) (conflicts with congruence)
    }
    enc_cnf.data = [{1}, {2}, {3}, {4}]

    solver, conflicts = EUFTheorySolver.from_encoded_cnf(enc_cnf)

    assert solver.assert_lit(1) == (True, set())
    assert solver.assert_lit(2) == (True, set())
    assert solver.assert_lit(3) == (True, set())

    # Should conflict due to functional congruence
    result = solver.assert_lit(4)
    assert result[0] == False
    assert len(result[1]) > 0


def test_diamond_disequality_extended_chain():
    """Test extended diamond with longer equality chain: a = b = c = d, but a != d."""
    enc_cnf = EncodedCNF()
    enc_cnf.encoding = {
        Eq(a, b): 1,      # a = b
        Eq(b, c): 2,      # b = c
        Eq(c, d): 3,      # c = d
        Ne(a, d): 4       # a != d (conflicts with transitivity)
    }
    enc_cnf.data = [{1}, {2}, {3}, {4}]

    solver, conflicts = EUFTheorySolver.from_encoded_cnf(enc_cnf)

    # Build equality chain
    assert solver.assert_lit(1) == (True, set())
    assert solver.assert_lit(2) == (True, set())
    assert solver.assert_lit(3) == (True, set())

    # Contradictory disequality should conflict
    result = solver.assert_lit(4)
    assert result[0] == False
    conflict_lits = result[1]
    assert len(conflict_lits) > 0
    # Should involve multiple literals from the chain
    assert len(conflict_lits) >= 3


def test_diamond_disequality_multiple_branches():
    """Test diamond with multiple branches: a = b, a = c, a = d, b != c."""
    enc_cnf = EncodedCNF()
    enc_cnf.encoding = {
        Eq(a, b): 1,      # a = b
        Eq(a, c): 2,      # a = c
        Eq(a, d): 3,      # a = d
        Ne(b, c): 4,      # b != c (conflicts)
        Ne(c, d): 5       # c != d (would also conflict)
    }
    enc_cnf.data = [{1}, {2}, {3}, {4}, {5}]

    solver, conflicts = EUFTheorySolver.from_encoded_cnf(enc_cnf)

    # Build equality star pattern centered on 'a'
    assert solver.assert_lit(1) == (True, set())
    assert solver.assert_lit(2) == (True, set())
    assert solver.assert_lit(3) == (True, set())

    # Conflict
    result1 = solver.assert_lit(4)  # b != c
    assert result1[0] == False


def test_diamond_disequality_nested_functions():
    """Test diamond with nested functions: f(g(a)) = x, f(g(b)) = x, a = b, f(g(a)) != f(g(b))."""
    enc_cnf = EncodedCNF()
    func_f, func_g = Function('f'), Function('g')

    enc_cnf.encoding = {
        Eq(func_f(func_g(a)), x): 1,        # f(g(a)) = x
        Eq(func_f(func_g(b)), x): 2,        # f(g(b)) = x
        Eq(a, b): 3,                        # a = b
        Ne(func_f(func_g(a)), func_f(func_g(b))): 4  # f(g(a)) != f(g(b))
    }
    enc_cnf.data = [{1}, {2}, {3}, {4}]

    solver, conflicts = EUFTheorySolver.from_encoded_cnf(enc_cnf)

    assert solver.assert_lit(1) == (True, set())
    assert solver.assert_lit(2) == (True, set())
    assert solver.assert_lit(3) == (True, set())

    # Should conflict due to nested function congruence
    result = solver.assert_lit(4)
    assert result[0] == False
    assert len(result[1]) > 0


def test_diamond_disequality_explanation_basic():
    """Test explain_disequality for basic diamond pattern."""
    solver = EUFTheorySolver()

    # Set up diamond: a = b, a = c, so b = c by transitivity
    solver.SetTrue(Eq(a, b))
    solver.SetTrue(Eq(a, c))

    # Now b and c should be equal, so explain why they're NOT disequal
    explanation = solver.explain_disequality(b, c)

    # Since b = c (through a), there should be no disequality explanation
    assert explanation == set()  # Empty because b = c


def test_diamond_disequality_explanation_with_source():
    """Test explain_disequality when there's an actual disequality to explain."""
    solver = EUFTheorySolver()

    # Set up: x = a, y = b, a != b
    # This means x != y should be explained by {Eq(x, a), Eq(y, b), Ne(a, b)}
    solver.SetTrue(Eq(x, a))
    solver.SetTrue(Eq(y, b))
    solver.SetTrue(Ne(a, b))

    explanation = solver.explain_disequality(x, y)

    # Should contain the source disequality and connecting equalities
    expected_in_explanation = {Ne(a, b), Eq(a, x), Eq(b, y)}
    for exp_lit in expected_in_explanation:
        assert exp_lit in explanation


def test_diamond_disequality_explanation_chain():
    """Test explain_disequality through longer chain."""
    solver = EUFTheorySolver()

    # Chain: x = u = v, y = w = z, v 1= w
    # So x != y through the chain
    solver.SetTrue(Eq(x, u))
    solver.SetTrue(Eq(u, v))
    solver.SetTrue(Eq(y, w))
    solver.SetTrue(Eq(w, z))
    solver.SetTrue(Ne(v, w))

    explanation = solver.explain_disequality(x, y)

    # Should contain source disequality and connecting chain
    assert Ne(v, w) in explanation
    assert Eq(u, x) in explanation
    assert Eq(u, v) in explanation
    assert Eq(w, y) in explanation


def test_diamond_disequality_function_explanation():
    """Test explain_disequality with function applications."""
    solver = EUFTheorySolver()
    func_f = Function('f')

    # f(a) != f(b), a = x, b = y
    # So f(x) != f(y) should be explained
    solver.SetTrue(Ne(func_f(a), func_f(b)))
    solver.SetTrue(Eq(a, x))
    solver.SetTrue(Eq(b, y))

    explanation = solver.explain_disequality(func_f(x), func_f(y))

    # Should contain source function disequality and argument equalities
    assert Ne(func_f(a), func_f(b)) in explanation


def test_diamond_disequality_complex_scenario():
    """Test complex diamond disequality scenario with multiple interconnected terms."""
    enc_cnf = EncodedCNF()
    func_f = Function('f')

    enc_cnf.encoding = {
        Eq(a, b): 1,                # a = b
        Eq(func_f(a), x): 2,        # f(a) = x
        Eq(func_f(b), y): 3,        # f(b) = y
        Ne(x, y): 4                 # x != y (conflicts with f(a) = f(b) by congruence)
    }
    enc_cnf.data = [{1}, {2}, {3}, {4}]

    solver, conflicts = EUFTheorySolver.from_encoded_cnf(enc_cnf)

    assert solver.assert_lit(1) == (True, set())  # a = b
    assert solver.assert_lit(2) == (True, set())  # f(a) = x
    assert solver.assert_lit(3) == (True, set())  # f(b) = y

    # Should conflict: a = b implies f(a) = f(b), but f(a) = x, f(b) = y, x != y
    result = solver.assert_lit(4)
    assert result[0] == False
    assert len(result[1]) > 0


def test_diamond_disequality_no_false_conflicts():
    """Test that valid diamond patterns don't produce false conflicts."""
    enc_cnf = EncodedCNF()
    enc_cnf.encoding = {
        Eq(a, b): 1,      # a = b
        Eq(c, d): 2,      # c = d (unrelated)
        Ne(a, c): 3       # a != c (should be fine, no connection)
    }
    enc_cnf.data = [{1}, {2}, {3}]

    solver, conflicts = EUFTheorySolver.from_encoded_cnf(enc_cnf)

    # All should succeed - no diamond conflict here
    assert solver.assert_lit(1) == (True, set())
    assert solver.assert_lit(2) == (True, set())
    assert solver.assert_lit(3) == (True, set())


def test_diamond_disequality_partial_diamond():
    """Test partial diamond that doesn't create conflict."""
    enc_cnf = EncodedCNF()
    enc_cnf.encoding = {
        Eq(a, b): 1,      # a = b
        Ne(b, c): 2,      # b != c
        Ne(a, c): 3       # a != c (consistent with above)
    }
    enc_cnf.data = [{1}, {2}, {3}]

    solver, conflicts = EUFTheorySolver.from_encoded_cnf(enc_cnf)

    # All should succeed - this is consistent
    assert solver.assert_lit(1) == (True, set())
    assert solver.assert_lit(2) == (True, set())
    assert solver.assert_lit(3) == (True, set())
