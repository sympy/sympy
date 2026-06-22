from __future__ import annotations
from sympy.core.numbers import Rational, I, oo
from sympy.core.relational import Eq
from sympy.core.symbol import symbols
from sympy.core.singleton import S
from sympy.matrices.dense import Matrix
from sympy.matrices.dense import randMatrix
from sympy.assumptions.ask import Q
from sympy.logic.boolalg import And
from sympy.abc import x, y, z
from sympy.assumptions.cnf import CNF, EncodedCNF
from sympy.functions.elementary.trigonometric import cos
from sympy.external import import_module

from sympy.logic.algorithms.lra_theory import LRASolver, UnhandledInput, LRARational, HANDLE_NEGATION
from sympy.core.random import random, choice, randint
from sympy.core.sympify import sympify
from sympy.ntheory.generate import randprime
from sympy.core.relational import StrictLessThan, StrictGreaterThan
import itertools

from sympy.testing.pytest import raises, XFAIL, skip

def make_random_problem(num_variables=2, num_constraints=2, sparsity=.1, rational=True,
                        disable_strict = False, disable_nonstrict=False, disable_equality=False):
    def rand(sparsity=sparsity):
        if random() < sparsity:
            return sympify(0)
        if rational:
            int1, int2 = [randprime(0, 50) for _ in range(2)]
            return Rational(int1, int2) * choice([-1, 1])
        else:
            return randint(1, 10) * choice([-1, 1])

    variables = symbols('x1:%s' % (num_variables + 1))
    constraints = []
    for _ in range(num_constraints):
        lhs, rhs = sum(rand() * x for x in variables), rand(sparsity=0) # sparsity=0  bc of bug with smtlib_code
        options = []
        if not disable_equality:
            options += [Eq(lhs, rhs)]
        if not disable_nonstrict:
            options += [lhs <= rhs, lhs >= rhs]
        if not disable_strict:
            options += [lhs < rhs, lhs > rhs]

        constraints.append(choice(options))

    return constraints

def check_if_satisfiable_with_z3(constraints):
    from sympy.external.importtools import import_module
    from sympy.printing.smtlib import smtlib_code
    from sympy.logic.boolalg import And
    boolean_formula = And(*constraints)
    z3 = import_module("z3")
    if z3:
        smtlib_string = smtlib_code(boolean_formula)
        s = z3.Solver()
        s.from_string(smtlib_string)
        res = str(s.check())
        if res == 'sat':
            return True
        elif res == 'unsat':
            return False
        else:
            raise ValueError(f"z3 was not able to check the satisfiability of {boolean_formula}")

def find_rational_assignment(constr, assignment, iter=20):
    eps = sympify(1)

    for _ in range(iter):
        assign = {key: val[0] + val[1]*eps for key, val in assignment.items()}
        try:
            for cons in constr:
                assert cons.subs(assign) == True
            return assign
        except AssertionError:
            eps = eps/2

    return None

def boolean_formula_to_encoded_cnf(bf):
    cnf = CNF.from_prop(bf)
    enc = EncodedCNF()
    enc.from_cnf(cnf)
    return enc


def test_from_encoded_cnf():
    s1, s2 = symbols("s1 s2")

    # Test preprocessing
    # Example is from section 3 of paper.
    phi = (x >= 0) & ((x + y <= 2) | (x + 2 * y - z >= 6)) & (Eq(x + y, 2) | (x + 2 * y - z > 4))
    enc = boolean_formula_to_encoded_cnf(phi)
    lra, _ = LRASolver.from_encoded_cnf(enc, testing_mode=True)
    assert lra.A.shape == (2, 5)
    assert str(lra.slack) == '[_s1, _s2]'
    assert str(lra.nonslack) == '[x, y, z]'
    assert lra.A == Matrix([[ 1,  1, 0, -1,  0],
                            [-1, -2, 1,  0, -1]])
    actual = {tuple(sorted((str(b.var), b.bound, b.upper, b.strict) for b in bs)) for bs in lra.atom_id_to_boundaries.values()}
    expected = {
        (('_s1', 2, False, False), ('_s1', 2, True, False)), # Eq(x + y, 2)
        (('_s1', 2, True, False),),                          # x + y <= 2
        (('_s2', -4, True, True),),                          # x + 2*y - z > 4 -> _s2 < -4
        (('_s2', -6, True, False),),                         # x + 2*y - z >= 6 -> _s2 <= -6
        (('x', 0, False, False),)                            # x >= 0
    }
    assert actual == expected


def test_problem():
    cons = [-2 * x - 2 * y >= 7, -9 * y >= 7, -6 * y >= 5]
    cnf = CNF().from_prop(And(*cons))
    enc = EncodedCNF()
    enc.from_cnf(cnf)
    lra, _ = LRASolver.from_encoded_cnf(enc)
    lra.assert_lit(1)
    lra.assert_lit(2)
    lra.assert_lit(3)
    is_sat, assignment = lra.check()
    assert is_sat is True


def test_random_problems():
    z3 = import_module("z3")
    if z3 is None:
        skip("z3 is not installed")

    special_cases = []; x1, x2, x3 = symbols("x1 x2 x3")
    special_cases.append([x1 - 3 * x2 <= -5, 6 * x1 + 4 * x2 <= 0, -7 * x1 + 3 * x2 <= 3])
    special_cases.append([-3 * x1 >= 3, Eq(4 * x1, -1)])
    special_cases.append([-4 * x1 < 4, 6 * x1 <= -6])
    special_cases.append([-3 * x2 >= 7, 6 * x1 <= -5, -3 * x2 <= -4])
    special_cases.append([x + y >= 2, x + y <= 1])
    special_cases.append([x >= 0, x + y <= 2, x + 2 * y - z >= 6])  # from paper example
    special_cases.append([-2 * x1 - 2 * x2 >= 7, -9 * x1 >= 7, -6 * x1 >= 5])
    special_cases.append([2 * x1 > -3, -9 * x1 < -6, 9 * x1 <= 6])
    special_cases.append([-2*x1 < -4, 9*x1 > -9])
    special_cases.append([-6*x1 >= -1, -8*x1 + x2 >= 5, -8*x1 + 7*x2 < 4, x1 > 7])
    special_cases.append([Eq(x1, 2), Eq(5*x1, -2), Eq(-7*x2, -6), Eq(9*x1 + 10*x2, 9)])
    special_cases.append([Eq(3*x1, 6), Eq(x1 - 8*x2, -9), Eq(-7*x1 + 5*x2, 3), Eq(3*x2, 7)])
    special_cases.append([-4*x1 < 4, 6*x1 <= -6])
    special_cases.append([-3*x1 + 8*x2 >= -8, -10*x2 > 9, 8*x1 - 4*x2 < 8, 10*x1 - 9*x2 >= -9])
    special_cases.append([x1 + 5*x2 >= -6, 9*x1 - 3*x2 >= -9, 6*x1 + 6*x2 < -10, -3*x1 + 3*x2 < -7])
    special_cases.append([-9*x1 < 7, -5*x1 - 7*x2 < -1, 3*x1 + 7*x2 > 1, -6*x1 - 6*x2 > 9])
    special_cases.append([9*x1 - 6*x2 >= -7, 9*x1 + 4*x2 < -8, -7*x2 <= 1, 10*x2 <= -7])

    feasible_count = 0
    for i in range(50):
        if i % 8 == 0:
            constraints = make_random_problem(num_variables=1, num_constraints=2, rational=False)
        elif i % 8 == 1:
            constraints = make_random_problem(num_variables=2, num_constraints=4, rational=False, disable_equality=True,
                                              disable_nonstrict=True)
        elif i % 8 == 2:
            constraints = make_random_problem(num_variables=2, num_constraints=4, rational=False, disable_strict=True)
        elif i % 8 == 3:
            constraints = make_random_problem(num_variables=3, num_constraints=12, rational=False)
        else:
            constraints = make_random_problem(num_variables=3, num_constraints=6, rational=False)

        if i < len(special_cases):
            constraints = special_cases[i]

        if False in constraints or True in constraints:
            continue

        phi = And(*constraints)
        if phi == False:
            continue
        cnf = CNF.from_prop(phi); enc = EncodedCNF()
        enc.from_cnf(cnf)
        assert all(0 not in clause for clause in enc.data)

        lra, _ = LRASolver.from_encoded_cnf(enc, testing_mode=True)
        s_subs = lra.s_subs

        lra.run_checks = True
        s_subs_rev = {value: key for key, value in s_subs.items()}
        lits = {lit for clause in enc.data for lit in clause}

        bounds = [(lra.atom_id_to_boundaries[l], l) for l in lits if l in lra.atom_id_to_boundaries]
        bounds = sorted(bounds, key=lambda x: (str(x[0][0].var), x[0][0].bound, str(x[0][0].upper))) # to remove nondeterminism

        for b, l in bounds:
            res = lra.assert_lit(l)
            if res and res[0] == False:
                feasible = res
                break
        else:
            feasible = lra.check()

        if feasible[0] == True:
            feasible_count += 1
            assert check_if_satisfiable_with_z3(constraints) is True
            cons_funcs = [cons.func for cons in constraints]
            assignment = feasible[1]
            assignment = {key.var : value for key, value in assignment.items()}
            if not (StrictLessThan in cons_funcs or StrictGreaterThan in cons_funcs):
                assignment = {key: value[0] for key, value in assignment.items()}
                for cons in constraints:
                    assert cons.subs(assignment) == True

            else:
                rat_assignment = find_rational_assignment(constraints, assignment)
                assert rat_assignment is not None
        else:
            assert check_if_satisfiable_with_z3(constraints) is False

            conflict = feasible[1]
            assert len(conflict) >= 2
            def get_expr(bs):
                if len(bs) == 2:
                    return Eq(bs[0].var.var, bs[0].bound)
                return bs[0].get_inequality()

            conflict = {get_expr(lra.atom_id_to_boundaries[abs(l)]) for l in conflict}
            conflict = {clause.subs(s_subs_rev) for clause in conflict}
            assert check_if_satisfiable_with_z3(conflict) is False

            # check that conflict clause is probably minimal
            for subset in itertools.combinations(conflict, len(conflict)-1):
                assert check_if_satisfiable_with_z3(subset) is True


@XFAIL
def test_pos_neg_zero():
    bf = Q.positive(x) & Q.negative(x) & Q.zero(y)
    enc = boolean_formula_to_encoded_cnf(bf)
    lra, _ = LRASolver.from_encoded_cnf(enc, testing_mode=True)
    for lit in enc.encoding.values():
        if lra.assert_lit(lit) is not None:
            break
    assert len(lra.atom_id_to_boundaries) == 3
    assert lra.check()[0] == False

    bf = Q.positive(x) & Q.lt(x, -1)
    enc = boolean_formula_to_encoded_cnf(bf)
    lra, _ = LRASolver.from_encoded_cnf(enc, testing_mode=True)
    for lit in enc.encoding.values():
        if lra.assert_lit(lit) is not None:
            break
    assert len(lra.atom_id_to_boundaries) == 2
    assert lra.check()[0] == False

    bf = Q.positive(x) & Q.zero(x)
    enc = boolean_formula_to_encoded_cnf(bf)
    lra, _ = LRASolver.from_encoded_cnf(enc, testing_mode=True)
    for lit in enc.encoding.values():
        if lra.assert_lit(lit) is not None:
            break
    assert len(lra.atom_id_to_boundaries) == 2
    assert lra.check()[0] == False

    bf = Q.positive(x) & Q.zero(y)
    enc = boolean_formula_to_encoded_cnf(bf)
    lra, _ = LRASolver.from_encoded_cnf(enc, testing_mode=True)
    for lit in enc.encoding.values():
        if lra.assert_lit(lit) is not None:
            break
    assert len(lra.atom_id_to_boundaries) == 2
    assert lra.check()[0] == True


@XFAIL
def test_pos_neg_infinite():
    bf = Q.positive_infinite(x) & Q.lt(x, 10000000) & Q.positive_infinite(y)
    enc = boolean_formula_to_encoded_cnf(bf)
    lra, _ = LRASolver.from_encoded_cnf(enc, testing_mode=True)
    for lit in enc.encoding.values():
        if lra.assert_lit(lit) is not None:
            break
    assert len(lra.atom_id_to_boundaries) == 3
    assert lra.check()[0] == False

    bf = Q.positive_infinite(x) & Q.gt(x, 10000000) & Q.positive_infinite(y)
    enc = boolean_formula_to_encoded_cnf(bf)
    lra, _ = LRASolver.from_encoded_cnf(enc, testing_mode=True)
    for lit in enc.encoding.values():
        if lra.assert_lit(lit) is not None:
            break
    assert len(lra.atom_id_to_boundaries) == 3
    assert lra.check()[0] == True

    bf = Q.positive_infinite(x) & Q.negative_infinite(x)
    enc = boolean_formula_to_encoded_cnf(bf)
    lra, _ = LRASolver.from_encoded_cnf(enc, testing_mode=True)
    for lit in enc.encoding.values():
        if lra.assert_lit(lit) is not None:
            break
    assert len(lra.atom_id_to_boundaries) == 2
    assert lra.check()[0] == False


def test_binrel_evaluation():
    bf = Q.gt(3, 2)
    enc = boolean_formula_to_encoded_cnf(bf)
    lra, conflicts = LRASolver.from_encoded_cnf(enc, testing_mode=True)
    assert len(lra.atom_id_to_boundaries) == 0
    assert conflicts == [[1]]

    bf = Q.lt(3, 2)
    enc = boolean_formula_to_encoded_cnf(bf)
    lra, conflicts = LRASolver.from_encoded_cnf(enc, testing_mode=True)
    assert len(lra.atom_id_to_boundaries) == 0
    assert conflicts == [[-1]]


def test_negation():
    assert HANDLE_NEGATION is True
    bf = Q.gt(x, 1) & ~Q.gt(x, 0)
    enc = boolean_formula_to_encoded_cnf(bf)
    lra, _ = LRASolver.from_encoded_cnf(enc, testing_mode=True)
    conflict = None
    for clause in enc.data:
        for lit in clause:
            res = lra.assert_lit(lit)
            if res is not None:
                conflict = res[1]
    assert len(lra.atom_id_to_boundaries) == 2
    assert conflict is not None
    assert sorted(conflict) in [[-1, 2], [-2, 1]]

    bf = ~Q.gt(x, 1) & ~Q.lt(x, 0)
    enc = boolean_formula_to_encoded_cnf(bf)
    lra, _ = LRASolver.from_encoded_cnf(enc, testing_mode=True)
    conflict_found = False
    for clause in enc.data:
        for lit in clause:
            if lra.assert_lit(lit) is not None:
                conflict_found = True
                break
    assert len(lra.atom_id_to_boundaries) == 2
    assert conflict_found is False
    assert lra.check()[0] is True

    bf = ~Q.gt(x, 0) & ~Q.lt(x, 1)
    enc = boolean_formula_to_encoded_cnf(bf)
    lra, _ = LRASolver.from_encoded_cnf(enc, testing_mode=True)
    conflict_found = False
    for clause in enc.data:
        for lit in clause:
            if lra.assert_lit(lit) is not None:
                conflict_found = True
                break
    assert len(lra.atom_id_to_boundaries) == 2
    assert conflict_found is True

    bf = ~Q.gt(x, 0) & ~Q.le(x, 0)
    enc = boolean_formula_to_encoded_cnf(bf)
    lra, _ = LRASolver.from_encoded_cnf(enc, testing_mode=True)
    conflict_found = False
    for clause in enc.data:
        for lit in clause:
            if lra.assert_lit(lit) is not None:
                conflict_found = True
                break

    assert len(lra.atom_id_to_boundaries) == 2
    assert conflict_found is True

    bf = ~Q.le(x+y, 2) & ~Q.ge(x-y, 2) & ~Q.ge(y, 0)
    enc = boolean_formula_to_encoded_cnf(bf)
    lra, _ = LRASolver.from_encoded_cnf(enc, testing_mode=True)
    conflict_found = False
    for clause in enc.data:
        for lit in clause:
            if lra.assert_lit(lit) is not None:
                conflict_found = True
                break
    assert conflict_found is False
    assert len(lra.atom_id_to_boundaries) == 3
    is_sat, conflict = lra.check()
    assert is_sat is False
    assert len(conflict) == 3
    assert all(i > 0 for i in conflict)


def test_unhandled_input():
    nan = S.NaN
    bf = Q.gt(3, nan) & Q.gt(x, nan)
    enc = boolean_formula_to_encoded_cnf(bf)
    raises(ValueError, lambda: LRASolver.from_encoded_cnf(enc, testing_mode=True))

    bf = Q.gt(3, I) & Q.gt(x, I)
    enc = boolean_formula_to_encoded_cnf(bf)
    raises(UnhandledInput, lambda: LRASolver.from_encoded_cnf(enc, testing_mode=True))

    bf = Q.gt(3, float("inf")) & Q.gt(x, float("inf"))
    enc = boolean_formula_to_encoded_cnf(bf)
    raises(UnhandledInput, lambda: LRASolver.from_encoded_cnf(enc, testing_mode=True))

    bf = Q.gt(3, oo) & Q.gt(x, oo)
    enc = boolean_formula_to_encoded_cnf(bf)
    raises(UnhandledInput, lambda: LRASolver.from_encoded_cnf(enc, testing_mode=True))

    # test non-linearity
    bf = Q.gt(x**2 + x, 2)
    enc = boolean_formula_to_encoded_cnf(bf)
    raises(UnhandledInput, lambda: LRASolver.from_encoded_cnf(enc, testing_mode=True))

    bf = Q.gt(cos(x) + x, 2)
    enc = boolean_formula_to_encoded_cnf(bf)
    raises(UnhandledInput, lambda: LRASolver.from_encoded_cnf(enc, testing_mode=True))


@XFAIL
def test_infinite_strict_inequalities():
    # Extensive testing of the interaction between strict inequalities
    # and constraints containing infinity is needed because
    # the paper's rule for strict inequalities don't work when
    # infinite numbers are allowed. Using the paper's rules you
    # can end up with situations where oo + delta > oo is considered
    # True when oo + delta should be equal to oo.
    # See https://math.stackexchange.com/questions/4757069/can-this-method-of-converting-strict-inequalities-to-equisatisfiable-nonstrict-i
    bf = (-x - y >= -float("inf")) & (x > 0) & (y >= float("inf"))
    enc = boolean_formula_to_encoded_cnf(bf)
    lra, _ = LRASolver.from_encoded_cnf(enc, testing_mode=True)
    for lit in sorted(enc.encoding.values()):
        if lra.assert_lit(lit) is not None:
            break
    assert len(lra.atom_id_to_boundaries) == 3
    assert lra.check()[0] == True


def test_pivot():
    for _ in range(10):
        m = randMatrix(5)
        rref = m.rref()
        for _ in range(5):
            i, j = randint(0, 4), randint(0, 4)
            if m[i, j] != 0:
                assert LRASolver._pivot(m, i, j).rref() == rref


def test_reset_bounds():
    """
    Tests that reset_bounds properly resets all state variables to their default values.
    """
    # Test solver behavior after reset
    bf = Q.ge(x, 1) & Q.lt(x, 1)
    enc = boolean_formula_to_encoded_cnf(bf)
    lra, _ = LRASolver.from_encoded_cnf(enc, testing_mode=True)

    conflict_found = False
    for clause in enc.data:
        for lit in clause:
            if lra.assert_lit(lit) is not None:
                conflict_found = True
                break

    assert conflict_found is True
    lra.reset_bounds()
    assert lra.check()[0] is True

    # Test individual state variable resets
    bf = Q.ge(x, 0) & Q.le(x, 1)
    enc = boolean_formula_to_encoded_cnf(bf)
    lra, _ = LRASolver.from_encoded_cnf(enc, testing_mode=True)

    state_variables = [
        ('lower', LRARational(10, 0), LRARational(-float("inf"), 0)),
        ('upper', LRARational(10, 0), LRARational(float("inf"), 0)),
        ('lower_literal', 5, None),
        ('upper_literal', -5, None),
        ('assign', LRARational(10, 0), LRARational(0, 0))
    ]

    for attr_name, test_value, expected_reset_value in state_variables:
        for var in lra.all_var:
            setattr(var, attr_name, test_value)

        lra.reset_bounds()

        for var in lra.all_var:
            actual_value = getattr(var, attr_name)
            assert actual_value == expected_reset_value, \
                f"Failed to reset {attr_name}: expected {expected_reset_value}, got {actual_value}"


def test_empty_cnf():
    cnf = CNF()
    enc = EncodedCNF()
    enc.from_cnf(cnf)
    lra, conflict = LRASolver.from_encoded_cnf(enc)
    assert len(conflict) == 0
    assert lra.check() == (True, {})


def test_example_from_paper():
    # Example from the section 4.6 of the paper.
    # https://link.springer.com/chapter/10.1007/11817963_11
    enc = EncodedCNF()
    cons = [
        x <= -4,
        x >= -8,
        -x + y <= 1,
        x + y >= -3
    ]
    for con in cons:
        enc.add_prop(con)

    lra, _ = LRASolver.from_encoded_cnf(enc, testing_mode=True)

    # Extracts the variables stored in the solver
    var_x = next(v for v in lra.all_var if str(v.var) == 'x')
    var_y = next(v for v in lra.all_var if str(v.var) == 'y')
    # var_s1 is a slack variable which corresponds for -x + y <= 1
    # var_s2 is a slack variable which corresponds for -x - y <= 3
    _s1 = lra.s_subs[-x + y]
    _s2 = lra.s_subs[-x - y]
    var_s1 = next(v for v in lra.all_var if v.var == _s1)
    var_s2 = next(v for v in lra.all_var if v.var == _s2)

    # State A_0
    assert var_x.assign == LRARational(0, 0)
    assert var_y.assign == LRARational(0, 0)
    assert var_s1.assign == LRARational(0, 0)
    assert var_s2.assign == LRARational(0, 0)

    # Assert x <= -4
    lra.assert_lit(1)
    is_sat, _ = lra.check()
    assert is_sat is True

    # State A_1
    assert var_x.upper == LRARational(-4, 0)
    assert var_x.assign == LRARational(-4, 0)
    assert var_y.assign == LRARational(0, 0)
    assert var_s1.assign == LRARational(4, 0)
    assert var_s2.assign == LRARational(4, 0)

    # Assert x >= -8
    lra.assert_lit(2)
    is_sat, _ = lra.check()
    assert is_sat is True

    # State A_2
    assert var_x.lower == LRARational(-8, 0)
    assert var_x.upper == LRARational(-4, 0)
    assert var_x.assign == LRARational(-4, 0)
    assert var_y.assign == LRARational(0, 0)
    assert var_s1.assign == LRARational(4, 0)
    assert var_s2.assign == LRARational(4, 0)

    # Asserts -x + y <= 1
    # Check is invoked to pivot s1 and y
    # y's range is (-inf, -3]
    lra.assert_lit(3)
    is_sat, _ = lra.check()
    assert is_sat is True

    # State A_3
    assert var_s1.upper == LRARational(1, 0)
    assert var_x.lower == LRARational(-8, 0)
    assert var_x.upper == LRARational(-4, 0)
    assert var_x.assign == LRARational(-4, 0)
    assert var_y.assign == LRARational(-3, 0)
    assert var_s1.assign == LRARational(1, 0)
    assert var_s2.assign == LRARational(7, 0)

    # Assert -x - y <= 3 (s2)
    # s2 and s1 are conflicting assertions as for both to be true
    # x >= -2 but x's range is [-8, -4]
    res = lra.assert_lit(4)
    assert res is None
    is_sat, _ = lra.check()
    assert is_sat is False

    # Backtrack to remove the conflicted assertion
    lra.backtrack()
    is_sat, _ = lra.check()
    assert is_sat is True

    # State A_3 after backtracking
    assert var_s1.upper == LRARational(1, 0)
    assert var_x.lower == LRARational(-8, 0)
    assert var_x.upper == LRARational(-4, 0)
    assert var_x.assign == LRARational(-4, 0)
    assert var_y.assign == LRARational(-3, 0)
    assert var_s1.assign == LRARational(1, 0)
    assert var_s2.assign == LRARational(7, 0)


def test_backtracking_single_variable():
    # This test is for checking the correctness over a single variable
    cons = [x >= -8, x <= -4, x >= -2]
    enc = EncodedCNF()
    for con in cons:
        enc.add_prop(con)
    lra, _ = LRASolver.from_encoded_cnf(enc, testing_mode=True)

    # Assert x in [-8, -4]
    lra.assert_lit(1)
    lra.assert_lit(2)
    is_sat, _ = lra.check()
    assert is_sat is True

    # Asserts x >= -2
    res = lra.assert_lit(3)
    # This directly contradicts x <= -4 which `assert_lit` catches instantly
    assert res is not None
    is_sat, _ = res
    assert is_sat is False

    lra.backtrack()
    is_sat, _ = lra.check()
    assert is_sat is True


def test_backtracking_multiple_variables():
    # This test is for checking the correctness over multiple variables
    enc = EncodedCNF()
    cons = [2*x + 3*y <= 12, x >= 3, y >= 3]
    for con in cons:
        enc.add_prop(con)

    lra, _ = LRASolver.from_encoded_cnf(enc, testing_mode=True)
    lra.assert_lit(1)
    lra.assert_lit(2)
    is_sat, _ = lra.check()
    assert is_sat is True

    # If 2x + 3y <= 12 and x >= 3, then y <= 2
    # We are asserting y >= 3 which is wrong
    res = lra.assert_lit(3)
    assert res is None
    is_sat, _ = lra.check()
    assert is_sat is False

    # backtracking to remove the faulty assert
    lra.backtrack()
    is_sat, _ = lra.check()
    assert is_sat is True


def test_backtracking_single_variable_multiple_backtracks():
    # This test is for checking correctness over multiple backtracking
    # Range of x should be [0, 2]
    enc = EncodedCNF()
    cons = [x <= 10, x >= 0, x >= 5, x <= 2]
    for con in cons:
        enc.add_prop(con)
    lra, _ = LRASolver.from_encoded_cnf(enc, testing_mode=True)

    # Setting 5 <= x <= 10
    lra.assert_lit(1)
    lra.assert_lit(2)
    lra.assert_lit(3)

    # x <= 2 is impossible while x >= 5 is present
    res = lra.assert_lit(4)
    assert res is not None
    is_sat, _ = res
    assert is_sat is False

    # First backtrack: Undo x <= 2 to resolve the conflict
    lra.backtrack()
    is_sat, _ = lra.check()
    assert is_sat is True

    # Second backtrack: Erase the x >= 5 constraint.
    # This widens the valid domain back to [0, 10],
    # allowing to accept the previously conflicting x <= 2 rule.
    lra.backtrack()

    # Setting 0 <= x <= 2
    lra.assert_lit(4)
    is_sat, _ = lra.check()
    assert is_sat is True


def test_backtracking_multiple_variables_multiple_backtracks():
    # This test is for checking correctness over multiple backtracking
    # for multiple variables, we need to constraint 6 to be True
    enc = EncodedCNF()
    cons = [
        x <= 10,
        x >= 0,
        y >= 0,
        x >= 5,
        y >= 5,
        x + y <= 4
    ]

    for con in cons:
        enc.add_prop(con)
    lra, _ = LRASolver.from_encoded_cnf(enc, testing_mode=True)

    # Establish the base valid state (x in [0, 10], y>=0)
    lra.assert_lit(1)
    lra.assert_lit(2)
    lra.assert_lit(3)
    is_sat, _ = lra.check()
    assert is_sat is True

    # Now x >= 5 and y >= 5, x + y >= 10
    lra.assert_lit(4)
    lra.assert_lit(5)
    is_sat, _ = lra.check()
    assert is_sat is True

    # x + y <= 4 is mathematically impossible when x>=5 and y>=5
    res = lra.assert_lit(6)
    assert res is None
    is_sat, _ = lra.check()
    assert is_sat is False

    # First backtrack: Pop the conflicting rule (Rule 6)
    lra.backtrack()
    is_sat, _ = lra.check()
    assert is_sat is True

    # Second and Third backtrack: pop the restrictive constraints
    # This restores the domain to x >= 0, y >= 0
    lra.backtrack()
    lra.backtrack()

    # x + y <= 4 is mathematically possible now
    lra.assert_lit(6)
    is_sat, _ = lra.check()
    assert is_sat is True


def test_backtracking_empty_history():
    enc = EncodedCNF()
    lra, _ = LRASolver.from_encoded_cnf(enc, testing_mode=True)

    raises(ValueError, lambda: lra.backtrack())
