from sympy.core.function import Function
from sympy.core.numbers import Rational, pi, I, oo
from sympy.core.relational import Eq
from sympy.core.symbol import symbols
from sympy.core.singleton import S
from sympy.matrices.dense import Matrix
from sympy.matrices.dense import randMatrix
from sympy.assumptions.ask import Q
from sympy.logic.boolalg import And
from sympy.abc import x, y, z, a
from sympy.assumptions.cnf import CNF, EncodedCNF

from sympy.logic.algorithms.lra_theory import LRASolver, UnhandledNumber

from sympy.core.random import random, choice, randint
from sympy.core.sympify import sympify
from sympy.ntheory.generate import randprime

from sympy.testing.pytest import raises, XFAIL
import time

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
    assert str(lra.nonslack) == '[_x1, _x2, _x3]'
    assert lra.A == Matrix([[ 1,  1, 0, -1,  0],
                            [-1, -2, 1,  0, -1]])
    assert {(str(b.var), b.bound, b.upper, b.equality, b.strict) for b in lra.enc_to_boundry.values()} == {('_s1', 2, None, True, False),
    ('_s1', 2, True, False, False),
    ('_s2', -4, True, False, True),
    ('_s2', -6, True, False, False),
    ('_x1', 0, False, False, False)}


def test_random_problems():
    from sympy.core.relational import StrictLessThan, StrictGreaterThan
    import itertools

    special_cases = []; x1, x2, x3 = symbols("x1 x2 x3")
    #special_cases.append([x1 - 3 * x2 <= -5, 6 * x1 + 4 * x2 <= 0, -7 * x1 + 3 * x2 <= 3]) bug with smtlib_code
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

    check_time = 0.0
    from_encoded_time = 0.0
    assert_time = 0.0

    feasible_count = 0
    for i in range(300):
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
            constraints = special_cases[i-1]

        # constraints = make_random_problem(num_variables=2, num_constraints=4, rational=False, disable_strict=False,
        #                                  disable_nonstrict=False, disable_equality=False)


        if False in constraints or True in constraints:
            continue

        phi = And(*constraints)
        if phi == False:
            continue
        cnf = CNF.from_prop(phi); enc = EncodedCNF()
        enc.from_cnf(cnf)
        assert all(0 not in clause for clause in enc.data)

        start = time.time()
        lra, _ = LRASolver.from_encoded_cnf(enc, testing_mode=True)
        x_subs = lra.x_subs
        s_subs = lra.s_subs
        end = time.time()
        from_encoded_time += end - start

        lra.run_checks = True
        s_subs_rev = {value: key for key, value in s_subs.items()}
        lits = {lit for clause in enc.data for lit in clause}

        bounds = [(lra.enc_to_boundry[l], l) for l in lits if l in lra.enc_to_boundry]
        bounds = sorted(bounds, key=lambda x: (str(x[0].var), x[0].bound, str(x[0].upper))) # to remove nondeterminism

        start = time.time()
        for b, l in bounds:
            if lra.result and lra.result[0] == False:
                break
            #print("var:", b.var, "bound:", b.bound, "upper:", b.upper, "strict:", b.strict)
            lra.assert_lit(l)
            #print(lra.assign, lra.lower, lra.upper)

        end = time.time()
        assert_time += end - start

        start = time.time()
        feasible = lra.check()
        end = time.time()
        check_time += end - start

        if feasible[0] == True:
            feasible_count += 1
            assert check_if_satisfiable_with_z3(constraints) is True
            cons_funcs = [cons.func for cons in constraints]
            assignment = feasible[1]
            assignment = {key.var : value for key, value in assignment.items()}
            constraints = [constr.subs(x_subs) for constr in constraints]
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
            conflict = {lra.enc_to_boundry[-l].get_inequality() for l in conflict}
            conflict = {clause.subs(s_subs_rev) for clause in conflict}
            assert check_if_satisfiable_with_z3(conflict) is False

            # check that conflict clause is probably minimal
            for subset in itertools.combinations(conflict, len(conflict)-1):
                assert check_if_satisfiable_with_z3(subset) is True
    print("check", check_time)
    print("from_encoded", from_encoded_time)
    print("assert", assert_time)
    print("total", check_time+from_encoded_time+assert_time)
    print("sat count", feasible_count)

@XFAIL
def test_pos_neg_zero():
    bf = Q.positive(x) & Q.negative(x) & Q.zero(y)
    enc = boolean_formula_to_encoded_cnf(bf)
    lra, _ = LRASolver.from_encoded_cnf(enc, testing_mode=True)
    for lit in enc.encoding.values():
        if lra.assert_lit(lit) is not None:
            break
    assert len(lra.enc_to_boundry) == 3
    assert lra.check()[0] == False

    bf = Q.positive(x) & Q.lt(x, -1)
    enc = boolean_formula_to_encoded_cnf(bf)
    lra, _ = LRASolver.from_encoded_cnf(enc, testing_mode=True)
    for lit in enc.encoding.values():
        if lra.assert_lit(lit) is not None:
            break
    assert len(lra.enc_to_boundry) == 2
    assert lra.check()[0] == False

    bf = Q.positive(x) & Q.zero(x)
    enc = boolean_formula_to_encoded_cnf(bf)
    lra, _ = LRASolver.from_encoded_cnf(enc, testing_mode=True)
    for lit in enc.encoding.values():
        if lra.assert_lit(lit) is not None:
            break
    assert len(lra.enc_to_boundry) == 2
    assert lra.check()[0] == False

    bf = Q.positive(x) & Q.zero(y)
    enc = boolean_formula_to_encoded_cnf(bf)
    lra, _ = LRASolver.from_encoded_cnf(enc, testing_mode=True)
    for lit in enc.encoding.values():
        if lra.assert_lit(lit) is not None:
            break
    assert len(lra.enc_to_boundry) == 2
    assert lra.check()[0] == True


@XFAIL
def test_pos_neg_infinite():
    bf = Q.positive_infinite(x) & Q.lt(x, 10000000) & Q.positive_infinite(y)
    enc = boolean_formula_to_encoded_cnf(bf)
    lra, _ = LRASolver.from_encoded_cnf(enc, testing_mode=True)
    for lit in enc.encoding.values():
        if lra.assert_lit(lit) is not None:
            break
    assert len(lra.enc_to_boundry) == 3
    assert lra.check()[0] == False

    bf = Q.positive_infinite(x) & Q.gt(x, 10000000) & Q.positive_infinite(y)
    enc = boolean_formula_to_encoded_cnf(bf)
    lra, _ = LRASolver.from_encoded_cnf(enc, testing_mode=True)
    for lit in enc.encoding.values():
        if lra.assert_lit(lit) is not None:
            break
    assert len(lra.enc_to_boundry) == 3
    assert lra.check()[0] == True

    bf = Q.positive_infinite(x) & Q.negative_infinite(x)
    enc = boolean_formula_to_encoded_cnf(bf)
    lra, _ = LRASolver.from_encoded_cnf(enc, testing_mode=True)
    for lit in enc.encoding.values():
        if lra.assert_lit(lit) is not None:
            break
    assert len(lra.enc_to_boundry) == 2
    assert lra.check()[0] == False


def test_binrel_evaluation():
    bf = Q.gt(3, 2)
    enc = boolean_formula_to_encoded_cnf(bf)
    lra, conflicts = LRASolver.from_encoded_cnf(enc, testing_mode=True)
    assert len(lra.enc_to_boundry) == 0
    assert conflicts == [[1]]

    bf = Q.lt(3, 2)
    enc = boolean_formula_to_encoded_cnf(bf)
    lra, conflicts = LRASolver.from_encoded_cnf(enc, testing_mode=True)
    assert len(lra.enc_to_boundry) == 0
    assert conflicts == [[-1]]


@XFAIL
def test_negation_handling():
    # note that SymPy assumes variables are complex by default
    # thus Not(x > 0) implies (x <= 0) OR (x has imaginary component)

    bf = ~Q.gt(x, 0) & Q.gt(x, 1)
    enc = boolean_formula_to_encoded_cnf(bf)
    lra, _ = LRASolver.from_encoded_cnf(enc)
    for lit in enc.encoding.values():
        if lra.assert_lit(lit) is not None:
            break
    assert lra.check()[0] is False

    bf = ~Q.gt(x, 0) & ~Q.le(x, 0)
    enc = boolean_formula_to_encoded_cnf(bf)
    lra, _ = LRASolver.from_encoded_cnf(enc)
    for lit in enc.encoding.values():
        if lra.assert_lit(lit) is not None:
            break
    assert lra.check()[0] is True

    bf = ~Q.gt(x, 0) & ~Q.le(x, 0) & Q.lt(x, 100)
    enc = boolean_formula_to_encoded_cnf(bf)
    lra, _ = LRASolver.from_encoded_cnf(enc)
    for lit in enc.encoding.values():
        if lra.assert_lit(lit) is not None:
            break
    res = lra.check()
    assert res[0] is False
    assert len(res[1]) == 3


def test_unhandled_input():
    nan = S.NaN
    bf = Q.gt(3, nan) & Q.gt(x, nan)
    enc = boolean_formula_to_encoded_cnf(bf)
    raises(ValueError, lambda: LRASolver.from_encoded_cnf(enc, testing_mode=True))

    bf = Q.gt(3, I) & Q.gt(x, I)
    enc = boolean_formula_to_encoded_cnf(bf)
    raises(UnhandledNumber, lambda: LRASolver.from_encoded_cnf(enc, testing_mode=True))

    bf = Q.gt(3, float("inf")) & Q.gt(x, float("inf"))
    enc = boolean_formula_to_encoded_cnf(bf)
    raises(UnhandledNumber, lambda: LRASolver.from_encoded_cnf(enc, testing_mode=True))

    bf = Q.gt(3, oo) & Q.gt(x, oo)
    enc = boolean_formula_to_encoded_cnf(bf)
    raises(UnhandledNumber, lambda: LRASolver.from_encoded_cnf(enc, testing_mode=True))

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
    assert len(lra.enc_to_boundry) == 3
    assert lra.check()[0] == True



def test_pivot():
    for _ in range(10):
        m = randMatrix(5)
        rref = m.rref()
        for _ in range(5):
            i, j = randint(0, 4), randint(0, 4)
            if m[i, j] != 0:
                assert LRASolver._pivot(m, i, j).rref() == rref