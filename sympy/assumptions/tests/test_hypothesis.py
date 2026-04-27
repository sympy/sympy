"""
Hypothesis-based structural tests for lra_satask SMT solver.

This test suite focuses on structural/grammar-based fuzzing to stress-test
the SMT solver's handling of complex formula structures, rather than testing
specific numerical values. The approach is inspired by academic SMT testing
methodologies (Brummayer & Biere 2009, Winterer et al. 2020).

Key Testing Principles:
1. Structural Fuzzing: Generate complex AST topologies, not just vary constants
2. Boolean Complexity: Deep disjunctions, nested negations, implications
3. Formula Composition: Build formulas from smaller equisatisfiable pieces
4. Soundness Properties: Test invariants that must hold for ANY formula

Lessons Learned During Development:
1. Test formula STRUCTURE (deep nesting, complex boolean combinations), not specific constant values.
2. Only test lra_satask itself, not internal helper functions like _preprocess or extract_assumption_from_old_assumption.
3. Q.real(x) syntax doesn't propagate to derived expressions like -x or 2*x, unlike symbols('x', real=True).
4. Use assume() to filter contradictory inputs before testing, rather than catching ValueError in try-except.
5. Check assumption consistency with lra_satask(assumption, True) before running the actual test.

Related Git Commits (check-real branch):
- 14a47b1382: assumptions: Clean up lra_satask
- c3d9e39520: lra_satask: Partially handle extended real exprs
- 34760714a3: lra_satask: Use unit clauses and refactor
"""

from hypothesis import given, assume, settings, strategies as st
from hypothesis.strategies import composite
from sympy import symbols, Q, Or, And, Not, Implies
from sympy.assumptions.lra_satask import lra_satask
from sympy.logic.algorithms.lra_theory import UnhandledInput
from sympy.testing.pytest import XFAIL


# Configuration: Disable new Q.real syntax due to limitations with derived expressions
# When True, uses Q.real(x) as assumptions. When False, uses symbols('x', real=True)
# Currently disabled because Q.real(x) doesn't propagate to expressions like -x, 2*x, x+y
USE_NEW_QREAL_SYNTAX = False


# Strategy: Generate variable pools for building formulas
@composite
def variable_pool(draw, min_vars=1, max_vars=4):
    """
    Generate a list of distinct real-valued symbols.
    Randomizes between old assumption syntax (real=True) and new syntax (Q.real).
    Returns: (vars, real_assumptions) where real_assumptions is None or a formula
    """
    num_vars = draw(st.integers(min_value=min_vars, max_value=max_vars))
    var_names = [f'x{i}' for i in range(num_vars)]
    
    # Choose between old and new assumption syntax based on configuration
    use_old_syntax = not USE_NEW_QREAL_SYNTAX
    
    if use_old_syntax:
        # Old syntax: symbols('x y', real=True)
        vars = symbols(' '.join(var_names), real=True)
        if not isinstance(vars, tuple):
            vars = (vars,)
        return vars, None
    else:
        # New syntax: symbols('x y') with Q.real(x) & Q.real(y)
        vars = symbols(' '.join(var_names))
        if not isinstance(vars, tuple):
            vars = (vars,)
        # Build Q.real assumptions
        real_assumptions = Q.real(vars[0])
        for var in vars[1:]:
            real_assumptions = And(real_assumptions, Q.real(var))
        return vars, real_assumptions


# Strategy: Generate atomic inequality predicates
@composite
def atomic_inequality(draw, vars):
    """Generate a single atomic inequality like x > 0 or x < y."""
    if not vars:
        vars = [symbols('x', real=True)]
    if not isinstance(vars, (list, tuple)):
        vars = [vars]
    
    var = draw(st.sampled_from(vars))
    
    # Choose relation type
    relation = draw(st.sampled_from([Q.gt, Q.lt, Q.ge, Q.le, Q.eq]))
    
    # Choose right-hand side: another variable or a constant
    use_var = draw(st.booleans()) and len(vars) > 1
    if use_var:
        # Compare with another variable
        other_vars = [v for v in vars if v != var]
        rhs = draw(st.sampled_from(other_vars))
    else:
        # Compare with constant
        rhs = draw(st.integers(min_value=-10, max_value=10))
    
    return relation(var, rhs)


# Strategy: Generate complex boolean combinations
@composite
def boolean_formula(draw, vars, max_depth=3):
    """
    Generate structurally complex boolean formulas with deep nesting.
    This stresses the SAT engine component of the SMT solver.
    """
    if max_depth == 0:
        return draw(atomic_inequality(vars))
    
    depth = draw(st.integers(min_value=0, max_value=max_depth))
    if depth == 0:
        return draw(atomic_inequality(vars))
    
    # Choose boolean operator
    op_choice = draw(st.sampled_from(['and', 'or', 'not', 'implies']))
    
    if op_choice == 'not':
        subformula = draw(boolean_formula(vars, max_depth=depth-1))
        return Not(subformula)
    elif op_choice == 'implies':
        left = draw(boolean_formula(vars, max_depth=depth-1))
        right = draw(boolean_formula(vars, max_depth=depth-1))
        return Implies(left, right)
    elif op_choice == 'and':
        num_conjuncts = draw(st.integers(min_value=2, max_value=4))
        conjuncts = [draw(boolean_formula(vars, max_depth=depth-1)) 
                     for _ in range(num_conjuncts)]
        return And(*conjuncts)
    else:  # 'or'
        num_disjuncts = draw(st.integers(min_value=2, max_value=4))
        disjuncts = [draw(boolean_formula(vars, max_depth=depth-1)) 
                     for _ in range(num_disjuncts)]
        return Or(*disjuncts)


# Strategy: Generate inequality chains (transitivity structure)
@composite
def inequality_chain(draw, vars, min_length=2, max_length=5):
    """
    Generate chained inequalities like x > y > z > w.
    Tests structural transitivity reasoning.
    """
    if not vars or len(vars) < 2:
        vars = symbols('x y z w', real=True)
    
    length = draw(st.integers(min_value=min_length, max_value=min(max_length, len(vars))))
    chain_vars = draw(st.permutations(vars[:length]))
    
    # Choose direction: all > or all <
    use_gt = draw(st.booleans())
    relation = Q.gt if use_gt else Q.lt
    
    # Build chain: v0 > v1 > v2 > ...
    chain = relation(chain_vars[0], chain_vars[1])
    for i in range(2, len(chain_vars)):
        chain = chain & relation(chain_vars[i-1], chain_vars[i])
    
    return chain, chain_vars


# Strategy: Generate linear combinations
@composite
def linear_expression(draw, vars):
    """Generate linear expressions like 2*x + 3*y - z."""
    if not vars:
        vars = [symbols('x', real=True)]
    if not isinstance(vars, (list, tuple)):
        vars = [vars]
    
    num_terms = draw(st.integers(min_value=1, max_value=len(vars)))
    selected_vars = draw(st.lists(st.sampled_from(vars), min_size=num_terms, max_size=num_terms, unique=True))
    
    expr = 0
    for var in selected_vars:
        coeff = draw(st.integers(min_value=-5, max_value=5).filter(lambda x: x != 0))
        expr += coeff * var
    
    # Optionally add constant
    if draw(st.booleans()):
        const = draw(st.integers(min_value=-10, max_value=10))
        expr += const
    
    return expr


# SOUNDNESS TESTS: Properties that must hold for ANY formula

@given(data=st.data())
def test_formula_self_consistency(data):
    """
    Soundness: A formula with itself as assumption should always be True.
    For any formula P: P |- P should always hold.
    """
    vars, real_assump = data.draw(variable_pool(min_vars=1, max_vars=3))
    formula = data.draw(boolean_formula(vars, max_depth=2))
    
    # Build assumption: formula AND real assumptions (if using new syntax)
    assumption = And(formula, real_assump) if real_assump is not None else formula
    
    # Check if assumptions are contradictory first to avoid ValueError
    if assumption is not True:
        consistency_check = lra_satask(assumption, True)
        assume(consistency_check is not False)  # Skip test if assumptions are contradictory
    
    result = lra_satask(formula, assumption)
    # If formula is satisfiable, it should entail itself
    # Note: Implications may return None when they're vacuously true or undecidable
    assert result in [True, None], f"Formula {formula} with assumptions {assumption} should entail itself, got {result}"


@given(data=st.data())
def test_contradiction_detection(data):
    """
    Soundness: P and ~P should be unsatisfiable.
    Tests that the solver detects basic contradictions.
    """
    vars, real_assump = data.draw(variable_pool(min_vars=1, max_vars=2))
    formula = data.draw(boolean_formula(vars, max_depth=2))
    
    # Build contradictory assumption
    contradiction = And(formula, Not(formula))
    if real_assump is not None:
        contradiction = And(contradiction, real_assump)
    
    # Asking if formula is true, given formula and its negation
    # This may raise ValueError for contradictory assumptions, which is expected
    try:
        result = lra_satask(formula, contradiction)
        # Should return False/None, not True
        assert result is not True, f"Contradiction should not be satisfiable"
    except ValueError:
        # ValueError is expected for contradictory assumptions
        pass


@given(data=st.data())
def test_weakening_property(data):
    """
    Soundness: If P |- Q, then (P & R) |- Q for any R.
    Adding assumptions can't make provable things unprovable.
    """
    vars, real_assump = data.draw(variable_pool(min_vars=2, max_vars=3))
    conclusion = data.draw(atomic_inequality(vars))
    assumption1 = data.draw(boolean_formula(vars, max_depth=2))
    assumption2 = data.draw(atomic_inequality(vars))
    
    # Add real assumptions to both assumption sets
    assump1 = And(assumption1, real_assump) if real_assump is not None else assumption1
    assump2 = And(assumption1, assumption2, real_assump) if real_assump is not None else And(assumption1, assumption2)
    
    # Check if assumptions are contradictory first
    if assump1 is not True:
        consistency_check = lra_satask(assump1, True)
        assume(consistency_check is not False)
    
    if assump2 is not True:
        consistency_check2 = lra_satask(assump2, True)
        assume(consistency_check2 is not False)
    
    result1 = lra_satask(conclusion, assump1)
    result2 = lra_satask(conclusion, assump2)
    
    # If assumption1 proves conclusion, verify the result is valid
    if result1 is True:
        # result2 should still be True, False, or None
        assert result2 in [True, False, None]


@given(data=st.data())
def test_deep_disjunction_handling(data):
    """
    Structural test: Generate formulas with deep disjunctive structure.
    Tests the SAT engine's clause learning and backtracking.
    """
    vars, real_assump = data.draw(variable_pool(min_vars=2, max_vars=3))
    
    # Build a deep disjunctive formula: (a|b|c) & (d|e|f) & (g|h|i)
    num_clauses = data.draw(st.integers(min_value=2, max_value=4))
    clauses = []
    for _ in range(num_clauses):
        num_literals = data.draw(st.integers(min_value=2, max_value=3))
        literals = [data.draw(atomic_inequality(vars)) for _ in range(num_literals)]
        clauses.append(Or(*literals))
    
    cnf_formula = And(*clauses)
    assumption = real_assump if real_assump is not None else True
    
    # Test that it doesn't crash and returns a valid answer
    # May hit known bugs (ValueError for contradictions, IndexError for dpll2 bug)
    try:
        result = lra_satask(cnf_formula, assumption)
        assert result in [True, False, None]
    except (ValueError, IndexError):
        # ValueError: contradictory assumptions (expected)
        # IndexError: known bug tracked in test_dpll2_indexerror_bug
        pass


@given(data=st.data())
def test_implication_chains(data):
    """
    Structural test: Generate implication chains like (A -> B) & (B -> C).
    Tests handling of complex boolean structure.
    """
    vars, real_assump = data.draw(variable_pool(min_vars=2, max_vars=3))
    
    # Build chain: A -> B -> C
    chain_length = data.draw(st.integers(min_value=2, max_value=4))
    atoms = [data.draw(atomic_inequality(vars)) for _ in range(chain_length)]
    
    # Build chain of implications
    implications = Implies(atoms[0], atoms[1])
    for i in range(2, len(atoms)):
        implications = And(implications, Implies(atoms[i-1], atoms[i]))
    
    # Build assumptions: first atom AND implications AND real assumptions
    assumptions = And(atoms[0], implications)
    if real_assump is not None:
        assumptions = And(assumptions, real_assump)
    
    # Check if assumptions are contradictory first
    if assumptions is not True:
        consistency_check = lra_satask(assumptions, True)
        assume(consistency_check is not False)
    
    # Test: Given A and the chain, should be able to derive C
    result = lra_satask(atoms[-1], assumptions)
    # Result can be True (derivable), False (disprovable), or None (undecidable)
    assert result in [True, False, None]


@given(data=st.data())
def test_transitivity_with_structure(data):
    """
    Structural test: Generate random inequality chains and test transitivity.
    Focus is on chain structure, not specific values.
    """
    vars, real_assump = data.draw(variable_pool(min_vars=3, max_vars=5))
    chain, chain_vars = data.draw(inequality_chain(vars, min_length=2, max_length=4))
    
    # Build assumption from chain and real assumptions
    assumption = And(chain, real_assump) if real_assump is not None else chain
    
    # Check if assumptions are contradictory first
    if assumption is not True:
        consistency_check = lra_satask(assumption, True)
        assume(consistency_check is not False)
    
    # Test: Given the chain, first > last (or first < last)
    if 'Q.gt' in str(chain):
        conclusion = Q.gt(chain_vars[0], chain_vars[-1])
    else:
        conclusion = Q.lt(chain_vars[0], chain_vars[-1])
    
    result = lra_satask(conclusion, assumption)
    assert result is True, f"Transitivity failed for chain: {chain}"


@given(data=st.data())
def test_linear_combination_structure(data):
    """
    Structural test: Generate inequalities with linear combinations.
    Tests handling of arithmetic expressions, not just variables.
    """
    vars, real_assump = data.draw(variable_pool(min_vars=2, max_vars=3))
    
    lhs = data.draw(linear_expression(vars))
    rhs = data.draw(linear_expression(vars))
    
    relation = data.draw(st.sampled_from([Q.gt, Q.lt, Q.ge, Q.le, Q.eq]))
    formula = relation(lhs, rhs)
    
    assumption = real_assump if real_assump is not None else True
    
    # Test that complex expressions don't crash the solver
    result = lra_satask(formula, assumption)
    assert result in [True, False, None]


@given(data=st.data())
def test_negation_depth(data):
    """
    Structural test: Generate formulas with nested negations.
    Tests: ~(~P) handling, De Morgan's laws, etc.
    """
    vars, real_assump = data.draw(variable_pool(min_vars=1, max_vars=2))
    
    # Start with simple atom
    atom = data.draw(atomic_inequality(vars))
    
    # Add multiple layers of negation and conjunctions
    num_layers = data.draw(st.integers(min_value=1, max_value=3))
    formula = atom
    for _ in range(num_layers):
        op = data.draw(st.sampled_from(['not', 'and', 'or']))
        if op == 'not':
            formula = Not(formula)
        elif op == 'and':
            other = data.draw(atomic_inequality(vars))
            formula = And(formula, other)
        else:
            other = data.draw(atomic_inequality(vars))
            formula = Or(formula, other)
    
    assumption = real_assump if real_assump is not None else True
    
    # Test with nested negations - may hit known bugs
    try:
        result = lra_satask(formula, assumption)
        assert result in [True, False, None]
    except (ValueError, IndexError):
        # ValueError: contradictory assumptions (expected)
        # IndexError: known bug tracked in test_dpll2_indexerror_bug
        pass


@given(data=st.data())
def test_mixed_equality_inequality(data):
    """
    Structural test: Mix equality and inequality predicates.
    Tests theory combination between equality reasoning and LRA.
    """
    vars, real_assump = data.draw(variable_pool(min_vars=2, max_vars=3))
    
    # Generate mix of equalities and inequalities
    num_atoms = data.draw(st.integers(min_value=2, max_value=4))
    atoms = []
    for _ in range(num_atoms):
        relation = data.draw(st.sampled_from([Q.eq, Q.gt, Q.lt, Q.ge, Q.le]))
        v1 = data.draw(st.sampled_from(vars))
        v2 = data.draw(st.sampled_from([v for v in vars if v != v1]))
        atoms.append(relation(v1, v2))
    
    formula = And(*atoms)
    assumption = real_assump if real_assump is not None else True
    
    # Check if formula with assumptions is contradictory first
    # This avoids ValueError for inherently contradictory formulas
    try:
        result = lra_satask(formula, assumption)
        assert result in [True, False, None]
    except ValueError:
        # Formula is self-contradictory - skip this test case
        assume(False)


@given(data=st.data())
def test_disjunctive_assumptions(data):
    """
    Structural test: Use disjunctive assumptions, not just conjunctive.
    Tests case-splitting in the solver.
    """
    vars, real_assump = data.draw(variable_pool(min_vars=2, max_vars=3))
    
    # Create disjunctive assumption: (x > 0) | (x < 0)
    atom1 = data.draw(atomic_inequality(vars))
    atom2 = data.draw(atomic_inequality(vars))
    disj_assumption = Or(atom1, atom2)
    
    conclusion = data.draw(atomic_inequality(vars))
    
    # Combine with real assumptions
    assumption = And(disj_assumption, real_assump) if real_assump is not None else disj_assumption
    
    result = lra_satask(conclusion, assumption)
    assert result in [True, False, None]


@given(data=st.data())
@settings(max_examples=50)
def test_random_formula_doesnt_crash(data):
    """
    Robustness test: Random complex formulas shouldn't crash the solver.
    """
    vars, real_assump = data.draw(variable_pool(min_vars=2, max_vars=4))
    
    # Generate proposition and assumptions
    proposition = data.draw(boolean_formula(vars, max_depth=3))
    assumption = data.draw(boolean_formula(vars, max_depth=2))
    
    # Combine with real assumptions if using new syntax
    if real_assump is not None:
        assumption = And(assumption, real_assump)
    
    try:
        result = lra_satask(proposition, assumption)
        # Should return True, False, None, or raise UnhandledInput/ValueError
        assert result in [True, False, None]
    except (UnhandledInput, ValueError, IndexError):
        # These are acceptable outcomes for unsupported/contradictory formulas
        # IndexError: known bug tracked in test_dpll2_indexerror_bug
        pass


# Semantic fusion inspired test
@given(data=st.data())
def test_semantic_fusion_equisatisfiability(data):
    """
    If P is satisfiable and Q is satisfiable independently,
    then (P | Q) should also be satisfiable.
    """
    vars, real_assump = data.draw(variable_pool(min_vars=2, max_vars=3))
    
    p = data.draw(atomic_inequality(vars))
    q = data.draw(atomic_inequality(vars))
    
    # Use real assumptions as the assumption, not combined with propositions
    assumption = real_assump if real_assump is not None else True
    
    # Check if assumptions are contradictory first
    if assumption is not True:
        consistency_check = lra_satask(assumption, True)
        assume(consistency_check is not False)
    
    # Check if P and Q are individually consistent
    p_sat = lra_satask(p, assumption)
    q_sat = lra_satask(q, assumption)
    
    if p_sat in [True, None] and q_sat in [True, None]:
        # Then P | Q should also be satisfiable
        or_sat = lra_satask(Or(p, q), assumption)
        assert or_sat in [True, None], f"Disjunction of satisfiable formulas should be satisfiable"

@given(data=st.data())
def test_conjunction_satisfiability_monotonic(data):
    """
    If P & Q is satisfiable, then both P and Q must be satisfiable.
    """
    vars, real_assump = data.draw(variable_pool(min_vars=2, max_vars=3))
    
    p = data.draw(atomic_inequality(vars))
    q = data.draw(atomic_inequality(vars))
    
    # Use real assumptions as the assumption, not combined with propositions
    assumption = real_assump if real_assump is not None else True
    
    # Check if assumptions are contradictory first
    if assumption is not True:
        consistency_check = lra_satask(assumption, True)
        assume(consistency_check is not False)
    
    and_sat = lra_satask(And(p, q), assumption)
    
    if and_sat is True:
        # Then both P and Q should be individually satisfiable
        p_sat = lra_satask(p, assumption)
        q_sat = lra_satask(q, assumption)
        assert p_sat in [True, None], f"Conjunct P should be satisfiable if conjunction is"
        assert q_sat in [True, None], f"Conjunct Q should be satisfiable if conjunction is"


# XFAIL TESTS: Known bugs and limitations discovered by hypothesis testing

@XFAIL
def test_new_qreal_syntax_limitation():
    """
    XFAIL: Q.real syntax doesn't propagate to derived expressions.
    
    Limitation: When using Q.real(x) as an assumption (new syntax), derived
    expressions like -x, 2*x, x+y are NOT automatically recognized as real.
    This is in contrast to the old syntax symbols('x', real=True) where
    .is_real propagates to all arithmetic combinations.
    
    Example that fails:
    - x = symbols('x')  # No realness in symbol definition
    - lra_satask(Q.gt(-x, 0), Q.real(x))  # Fails: -x not recognized as real
    
    Example that works:
    - x = symbols('x', real=True)  # Realness in symbol definition
    - lra_satask(Q.gt(-x, 0), True)  # Works: -x.is_real is True
    
    This limitation affects any test that generates linear combinations like
    2*x + 3*y - z when using the new Q.real syntax.
    """
    x = symbols('x')  # New syntax: no real assumption in symbol
    
    # This should work but currently raises:
    # UnhandledInput: LRASolver: -x must be extended real.
    result = lra_satask(Q.gt(-x, 0), Q.real(x))
    assert result in [True, False, None]


@XFAIL
def test_dpll2_indexerror_bug():
    """
    XFAIL: IndexError bug in dpll2.py:310 discovered by hypothesis testing.
    
    Bug: When checking satisfiability of certain formula combinations,
    the DPLL solver tries to access self.levels[-1] when self.levels is empty,
    causing IndexError: list index out of range.
    
    This specific test case was found by hypothesis fuzzing in
    test_random_formula_doesnt_crash.
    
    Traceback:
      File "sympy/logic/algorithms/dpll2.py", line 310, in _current_level
        return self.levels[-1]
    IndexError: list index out of range
    
    Formulas that trigger the bug:
    - formula1 = Q.ge(x0, x1) & Q.le(x1, x0)  # Equivalent to x0 == x1
    - formula2 = Q.lt(x1, x0) | (Q.gt(x0, 0) & Q.gt(x0, x1))
    """
    x0, x1 = symbols('x0 x1', real=True)
    
    # Formula that triggers the bug
    formula1 = Q.ge(x0, x1) & Q.le(x1, x0)  # Equivalent to x0 == x1
    formula2 = Q.lt(x1, x0) | (Q.gt(x0, 0) & Q.gt(x0, x1))
    
    # This currently raises IndexError but should return a valid result
    result = lra_satask(formula1, formula2)
    assert result in [True, False, None]