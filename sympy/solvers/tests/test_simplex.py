from sympy.core.numbers import Rational
from sympy.core.relational import Relational, Eq, Ne
from sympy.core.symbol import symbols
from sympy.core.sympify import sympify
from sympy.core.random import random, choice
from sympy.functions.elementary.miscellaneous import sqrt
from sympy.ntheory.generate import randprime
from sympy.assumptions.ask import Q
from sympy.matrices.dense import Matrix, eye
from sympy.solvers.solveset import linear_eq_to_matrix
from sympy.solvers.simplex import (lp,
    UnboundedLPError, InfeasibleLPError)

from sympy.external.importtools import import_module

from sympy.testing.pytest import raises

from sympy.abc import x, y, z


def test_lp():
    np = import_module("numpy")
    scipy = import_module("scipy")

    def get_results_with_scipy(objective, constraints, variables):
        if scipy is not None and np is not None:
            from sympy.solvers.inequalities import _np
            nonpos, rep, xx = _np(constraints, [])
            assert not rep  # only testing nonneg variables
            C, _D = linear_eq_to_matrix(objective, *variables)
            A, B = linear_eq_to_matrix(nonpos, *variables)
            assert _D[0] == 0  # scipy only deals with D = 0



            A_sci = Matrix([[A], [-eye(len(variables))]])
            B_sci = Matrix([[B], [Matrix([0] * len(variables))]])
            C_sci = C
            A_sci = np.array(A_sci.tolist())
            B_sci = np.array(B_sci.tolist())
            C_sci = np.array(C_sci.tolist())
            res = scipy.optimize.linprog(C_sci, A_ub=A_sci, b_ub=B_sci)
            return res

    r1 = y + 2*z <= 3
    r2 = -x - 3*z <= -2
    r3 = 2*x + y + 7*z <= 5
    constraints = [r1, r2, r3]
    objective = -x - y - 5 * z
    variables = [x, y, z]
    optimum, argmax = lp(max, objective, constraints)
    if scipy is not None and np is not None:
        scipy_res = get_results_with_scipy(objective, constraints, variables)
        assert optimum.evalf() == sympify(-scipy_res.fun)
    assert objective.subs(argmax) == optimum
    for constr in constraints:
        assert constr.subs(argmax) == True

    r1 = x - y + 2*z <= 3
    r2 = -x + 2*y - 3*z <= -2
    r3 = 2*x + y - 7*z <= -5
    constraints = [r1, r2, r3]
    objective = -x - y - 5*z
    variables = [x, y, z]
    optimum, argmax = lp(max, objective, constraints, [])
    if scipy is not None and np is not None:
        scipy_res = get_results_with_scipy(objective, constraints, variables)
        assert optimum.evalf() == sympify(-scipy_res.fun)
    assert objective.subs(argmax) == optimum
    for constr in constraints:
        assert constr.subs(argmax) == True

    r1 = x - y + 2*z <= -4
    r2 = -x + 2*y - 3*z <= 8
    r3 = 2*x + y - 7*z <= 10
    constraints = [r1, r2, r3]
    const = 2
    objective = -x-y-5*z+const # has constant term
    variables = [x, y, z]
    optimum, argmax = lp(max, objective, constraints, [])
    if scipy is not None and np is not None:
        scipy_res = get_results_with_scipy(objective-const, constraints, variables)
        assert optimum.evalf() == (sympify(-scipy_res.fun)+const)
    assert objective.subs(argmax) == optimum
    for constr in constraints:
        assert constr.subs(argmax) == True

    # Section 4 Problem 1 from
    # http://web.tecnico.ulisboa.pt/mcasquilho/acad/or/ftp/FergusonUCLA_LP.pdf
    # answer on page 55
    x1, x2, x3, x4 = symbols('x1 x2 x3 x4')
    r1 = x1 - x2 - 2*x3 - x4 <= 4
    r2 = 2*x1 + x3 -4*x4 <= 2
    r3 = -2*x1 + x2 + x4 <= 1
    optimum, argmax = lp(max,
        x1 - 2*x2 - 3*x3 - x4, [r1, r2, r3], [])
    assert optimum == 4
    assert list(argmax.values()) == [7, 0, 0, 3]

    # equality
    r1 = Eq(x, y)
    r2 = Eq(y, z)
    r3 = z <= 3
    constraints = [r1, r2, r3]
    objective = x + y + z  # if equalities are removed then objective could be x
    variables = [x, y, z]
    optimum, argmax = lp(max, objective, constraints, [])
    if scipy is not None and np is not None:
        scipy_res = get_results_with_scipy(objective, constraints, variables)
        assert optimum.evalf() == sympify(-scipy_res.fun)
    assert objective.subs(argmax) == optimum
    for constr in constraints:
        assert constr.subs(argmax) == True

    # Binary predicate
    r1 = Q.ge(x, y)
    r2 = Q.le(x, 4)
    constraints = [r1, r2]
    objective = x + y
    variables = [x, y, z]
    optimum, argmax = lp(max, objective, constraints, [])
    if scipy is not None and np is not None:
        scipy_res = get_results_with_scipy(objective, constraints, variables)
        assert optimum.evalf() == sympify(-scipy_res.fun)
    assert optimum == 8
    assert [x.subs(argmax), y.subs(argmax)] == [4, 4]

    # input contains Floats
    r1 = x - y + 2.0*z <= -4
    r2 = -x + 2*y - 3.0*z <= 8
    r3 = 2*x + y - 7*z <= 10
    constraints = [r1, r2, r3]
    objective = -x-y-5*z
    variables = [x, y, z]
    optimum, argmax = lp(max, objective, constraints, [])
    if scipy is not None and np is not None:
        scipy_res = get_results_with_scipy(objective, constraints, variables)
        assert optimum.evalf() == sympify(-scipy_res.fun)
    assert objective.subs(argmax) == optimum
    for constr in constraints:
        assert constr.subs(argmax) == True

    # input contains non-float or non-Rational
    r1 = x - y + sqrt(2) * z <= -4
    r2 = -x + 2*y - 3*z <= 8
    r3 = 2*x + y - 7*z <= 10
    raises(TypeError, lambda: lp(max, -x-y-5*z, [r1, r2, r3], []))

    r1 = x >= 0
    raises(UnboundedLPError, lambda: lp(max, x, [r1], []))
    r2 = x <= -1
    raises(InfeasibleLPError, lambda: lp(max, x, [r1, r2], []))

    # strict inequalities are not allowed
    r1 = x > 0
    raises(TypeError, lambda: lp(max, x, [r1], []))

    # not equals not allowed
    r1 = Ne(x, 0)
    raises(TypeError, lambda: lp(max, x, [r1], []))

    def make_random_problem(num_variables=2, num_constraints=2, sparsity=.1):
        def rand():
            if random() < sparsity:
                return sympify(0)
            int1, int2 = [randprime(0, 200) for _ in range(2)]
            return Rational(int1, int2)*choice([-1, 1])
        variables = symbols('x1:%s' % (num_variables + 1))
        constraints = [(sum(rand()*x for x in variables) <= rand())
                       for _ in range(num_constraints)]
        objective = sum(rand() * x for x in variables)
        return objective, constraints, variables

    # testing random problems
    if scipy is not None and np is not None:
        for _ in range(50):
            objective, constraints, variables = make_random_problem()
            constraints = [c for c in constraints if isinstance(c, Relational)] # in case c auto simplifies to True or False
            if len(constraints) == 0:
                continue

            # check lp maximization
            # scipy minimizes, so negative objective for it
            scipy_res = get_results_with_scipy(-objective, constraints, variables)
            if scipy_res.status == 0:
                optimum, argmax = lp(max, objective, constraints, [])
                scipy_op = -scipy_res.fun  # negated to give actual max
                assert abs(optimum.evalf() - scipy_op) < .1**10
                assert objective.subs(argmax) == optimum
                for constr in constraints:
                    assert constr.subs(argmax) == True
            elif scipy_res.status == 2:
                # scipy: problem is infeasible
                raises(InfeasibleLPError,
                       lambda: lp(max, objective, constraints, []))
            elif scipy_res.status == 3:
                # scipy: problem is unbounded
                raises(UnboundedLPError,
                       lambda: lp(max, objective, constraints, []))
            else:
                # scipy: either iteration limit reached or numerical difficulties
                pass

            # check lp minimization
            scipy_res = get_results_with_scipy(objective, constraints, variables)
            if scipy_res.status == 0:
                optimum, argmax = lp(min, objective, constraints, [])
                scipy_op = scipy_res.fun
                assert abs(optimum.evalf() - scipy_op) < .1**10
                assert objective.subs(argmax) == optimum
                for constr in constraints:
                    assert constr.subs(argmax) == True
            elif scipy_res.status == 2:
                # scipy: problem is infeasible
                raises(InfeasibleLPError,
                       lambda: lp(max, objective, constraints, []))
            elif scipy_res.status == 3:
                # scipy: problem is unbounded
                raises(UnboundedLPError,
                       lambda: lp(max, objective, constraints, []))
            else:
                # scipy: either iteration limit reached or numerical difficulties
                pass
