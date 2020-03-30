#
# The main tests for the code in single.py are currently located in
# sympy/solvers/tests/test_ode.py
#
r"""
This File contains test functions for the individual hints used for solving ODEs.

Examples of each solver will be returned by _get_examples_ode_sol_name_of_solver.

Examples should have a key 'XFAIL' which stores the list of hints if they are
expected to fail for that hint.

Functions that are for internal use:

1) _ode_solver_test(ode_examples) - It takes dictionary of examples returned by
   _get_examples method and tests them with their respective hints.

2) _test_particular_example(our_hint, example_name) - It tests the ODE example corresponding
   to the hint provided.

3) _test_all_hints(runxfail=False) - It is used to test all the examples with all the hints
  currently implemented. It calls _test_for_particular_hint() which outputs whether the given
  hint functions properly if it classifies the ODE example.
  If runxfail flag is set to True then it will only test the examples which are expected to fail.

  Everytime the ODE of partiular solver are added then _test_all_hints() is to execuetd to find
  the possible failures of different solver hints.

"""
from sympy import (cos, Derivative, diff,
    Eq, exp, log, pi, Rational, sin,
    sqrt, symbols, Ei, erfi)

from sympy.core import Function, Symbol
from sympy.functions import airyai, airybi, besselj, bessely
from sympy.integrals.risch import NonElementaryIntegral
from sympy.solvers.ode import classify_ode, dsolve
from sympy.solvers.ode.ode import allhints, _remove_redundant_solutions
from sympy.solvers.ode.single import (FirstLinear, ODEMatchError,
    SingleODEProblem, SingleODESolver)

from sympy.solvers.ode.subscheck import checkodesol

from sympy.testing.pytest import XFAIL, raises
import traceback


x = Symbol('x')
y = Symbol('y')
f = Function('f')
g = Function('g')
C1, C2, C3 = symbols('C1:4')


hint_message = """\
Hint did not match the example {example}.

The ODE is:
{eq}.

The expected hint was
{our_hint}\
"""

expected_sol_message = """\
Different solution found from dsolve for example {example}.

The ODE is:
{eq}

The expected solution was
{sol}

What dsolve returned is:
{dsolve_sol}\
"""

checkodesol_msg = """\
solution found is not correct for example {example}.

The ODE is:
{eq}\
"""

dsol_incorrect_msg = """\
solution returned by dsolve is incorrect when using {hint}.

The ODE is:
{eq}

The expected solution was
{sol}

what dsolve returned is:
{dsolve_sol}\
"""

exception_msg = """\
dsolve raised exception : {e}

when using {hint} for the example {example}

You can test this with:

from sympy.solvers.ode.tests.test_single import _test_particular_example

_test_particular_example('{hint}', '{example}')

The ODE is:
{eq}

\
"""

def _test_for_particular_hint(our_hint, ode_example, xfail):
    eq = ode_example['eq']
    expected_sol = ode_example['sol']
    example = ode_example['example_name']
    xpass = True
    if our_hint in classify_ode(eq):
        try:
            dsolve_sol = dsolve(eq, hint=our_hint)
            expected_checkodesol = [(True, 0) for i in range(len(expected_sol))]
            if len(expected_sol) == 1:
                expected_checkodesol = (True, 0)

            if checkodesol(eq, dsolve_sol) != expected_checkodesol:
                message = dsol_incorrect_msg.format(hint=our_hint, eq=eq, sol=expected_sol,dsolve_sol=dsolve_sol)
                raise AssertionError(message)
        except Exception as e:
            print(exception_msg.format(e=str(e), hint=our_hint, example=example, eq=eq))
            traceback.print_exc()
            xpass = False
    return xpass


def _ode_solver_test(ode_examples):
    our_hint = ode_examples['hint']
    for example in ode_examples['examples']:
        eq = ode_examples['examples'][example]['eq']
        sol = ode_examples['examples'][example]['sol']
        if our_hint not in classify_ode(eq):
            message = hint_message.format(example=example, eq=eq, our_hint=our_hint)
            raise AssertionError(message)

        dsolve_sol = dsolve(eq,hint=our_hint)
        if dsolve_sol not in sol:
            message = expected_sol_message.format(example=example, eq=eq, sol=sol, dsolve_sol=dsolve_sol)
            raise AssertionError(message)

        expected_checkodesol = [(True, 0) for i in range(len(sol))]
        if checkodesol(eq, sol) != expected_checkodesol:
            message = checkodesol.format(example=example, eq=eq)
            raise AssertionError(message)


def _test_all_hints(runxfail=False):
    all_hints = list(allhints)
    all_examples = _get_all_examples()
    for our_hint in all_hints:
        for ode_example in all_examples:
            xfail = our_hint in ode_example['XFAIL']
            if not runxfail:
                xpass = _test_for_particular_hint(our_hint, ode_example, xfail)
                if xpass and xfail:
                    print(ode_example,"is now working for hint",our_hint)
            else:
                if xfail:
                    xpass = _test_for_particular_hint(our_hint, ode_example, xfail)
                    if xpass:
                        print(ode_example,"is now working for hint",our_hint)


def _test_particular_example(our_hint, example_name):
    all_examples = _get_all_examples()
    for example in all_examples:
        if example['example_name'] == example_name:
            eq = example['eq']
            dsolve(eq, hint=our_hint)


def test_SingleODESolver():
    # Test that not implemented methods give NotImplementedError
    # Subclasses should override these methods.
    problem = SingleODEProblem(f(x).diff(x), f(x), x)
    solver = SingleODESolver(problem)
    raises(NotImplementedError, lambda: solver.matches())
    raises(NotImplementedError, lambda: solver.get_general_solution())
    raises(NotImplementedError, lambda: solver._matches())
    raises(NotImplementedError, lambda: solver._get_general_solution())

    # This ODE can not be solved by the FirstLinear solver. Here we test that
    # it does not match and the asking for a general solution gives
    # ODEMatchError

    problem = SingleODEProblem(f(x).diff(x) + f(x)*f(x), f(x), x)

    solver = FirstLinear(problem)
    raises(ODEMatchError, lambda: solver.get_general_solution())

    solver = FirstLinear(problem)
    assert solver.matches() is False

    #These are just test for order of ODE

    problem = SingleODEProblem(f(x).diff(x) + f(x), f(x), x)
    assert problem.order == 1

    problem = SingleODEProblem(f(x).diff(x,4) + f(x).diff(x,2) - f(x).diff(x,3), f(x), x)
    assert problem.order == 4


def test_nth_algebraic():
    eqn = Eq(Derivative(f(x), x), Derivative(g(x), x))
    sol = Eq(f(x), C1 + g(x))
    assert checkodesol(eqn, sol, order=1, solve_for_func=False)[0]
    assert sol == dsolve(eqn, f(x), hint='nth_algebraic'), dsolve(eqn, f(x), hint='nth_algebraic')
    assert sol == dsolve(eqn, f(x))

    eqn = (diff(f(x)) - x)*(diff(f(x)) + x)
    sol = [Eq(f(x), C1 - x**2/2), Eq(f(x), C1 + x**2/2)]
    assert checkodesol(eqn, sol, order=1, solve_for_func=False)[0]
    assert set(sol) == set(dsolve(eqn, f(x), hint='nth_algebraic'))
    assert set(sol) == set(dsolve(eqn, f(x)))

    eqn = (1 - sin(f(x))) * f(x).diff(x)
    sol = Eq(f(x), C1)
    assert checkodesol(eqn, sol, order=1, solve_for_func=False)[0]
    assert sol == dsolve(eqn, f(x), hint='nth_algebraic')
    assert sol == dsolve(eqn, f(x))

    M, m, r, t = symbols('M m r t')
    phi = Function('phi')
    eqn = Eq(-M * phi(t).diff(t),
             Rational(3, 2) * m * r**2 * phi(t).diff(t) * phi(t).diff(t,t))
    solns = [Eq(phi(t), C1), Eq(phi(t), C1 + C2*t - M*t**2/(3*m*r**2))]
    assert checkodesol(eqn, solns[0], order=2, solve_for_func=False)[0]
    assert checkodesol(eqn, solns[1], order=2, solve_for_func=False)[0]
    assert set(solns) == set(dsolve(eqn, phi(t), hint='nth_algebraic'))
    assert set(solns) == set(dsolve(eqn, phi(t)))

    eqn = f(x) * f(x).diff(x) * f(x).diff(x, x)
    sol = Eq(f(x), C1 + C2*x)
    assert checkodesol(eqn, sol, order=1, solve_for_func=False)[0]
    assert sol == dsolve(eqn, f(x), hint='nth_algebraic')
    assert sol == dsolve(eqn, f(x))

    eqn = f(x) * f(x).diff(x) * f(x).diff(x, x) * (f(x) - 1)
    sol = Eq(f(x), C1 + C2*x)
    assert checkodesol(eqn, sol, order=1, solve_for_func=False)[0]
    assert sol == dsolve(eqn, f(x), hint='nth_algebraic')
    assert sol == dsolve(eqn, f(x))

    eqn = f(x) * f(x).diff(x) * f(x).diff(x, x) * (f(x) - 1) * (f(x).diff(x) - x)
    solns = [Eq(f(x), C1 + x**2/2), Eq(f(x), C1 + C2*x)]
    assert checkodesol(eqn, solns[0], order=2, solve_for_func=False)[0]
    assert checkodesol(eqn, solns[1], order=2, solve_for_func=False)[0]
    assert set(solns) == set(dsolve(eqn, f(x), hint='nth_algebraic'))
    assert set(solns) == set(dsolve(eqn, f(x)))


def test_nth_algebraic_issue15999():
    eqn = f(x).diff(x) - C1
    sol = Eq(f(x), C1*x + C2) # Correct solution
    assert checkodesol(eqn, sol, order=1, solve_for_func=False) == (True, 0)
    assert dsolve(eqn, f(x), hint='nth_algebraic') == sol
    assert dsolve(eqn, f(x)) == sol


def test_nth_algebraic_redundant_solutions():
    # This one has a redundant solution that should be removed
    eqn = f(x)*f(x).diff(x)
    soln = Eq(f(x), C1)
    assert checkodesol(eqn, soln, order=1, solve_for_func=False)[0]
    assert soln == dsolve(eqn, f(x), hint='nth_algebraic')
    assert soln == dsolve(eqn, f(x))

    # This has two integral solutions and no algebraic solutions
    eqn = (diff(f(x)) - x)*(diff(f(x)) + x)
    sol = [Eq(f(x), C1 - x**2/2), Eq(f(x), C1 + x**2/2)]
    assert all(c[0] for c in checkodesol(eqn, sol, order=1, solve_for_func=False))
    assert set(sol) == set(dsolve(eqn, f(x), hint='nth_algebraic'))
    assert set(sol) == set(dsolve(eqn, f(x)))

    eqn = f(x) + f(x)*f(x).diff(x)
    solns = [Eq(f(x), 0),
             Eq(f(x), C1 - x)]
    assert all(c[0] for c in checkodesol(eqn, solns, order=1, solve_for_func=False))
    assert set(solns) == set(dsolve(eqn, f(x)))

    solns = [Eq(f(x), exp(x)),
             Eq(f(x), C1*exp(C2*x))]
    solns_final =  _remove_redundant_solutions(eqn, solns, 2, x)
    assert solns_final == [Eq(f(x), C1*exp(C2*x))]

    # This one needs a substitution f' = g.
    eqn = -exp(x) + (x*Derivative(f(x), (x, 2)) + Derivative(f(x), x))/x
    sol = Eq(f(x), C1 + C2*log(x) + exp(x) - Ei(x))
    assert checkodesol(eqn, sol, order=2, solve_for_func=False)[0]
    assert sol == dsolve(eqn, f(x))


#
# These tests can be combined with the above test if they get fixed
# so that dsolve actually works in all these cases.
#


# prep = True breaks this
def test_nth_algebraic_noprep1():
    eqn = Derivative(x*f(x), x, x, x)
    sol = Eq(f(x), (C1 + C2*x + C3*x**2) / x)
    assert checkodesol(eqn, sol, order=3, solve_for_func=False)[0]
    assert sol == dsolve(eqn, f(x), prep=False, hint='nth_algebraic')


@XFAIL
def test_nth_algebraic_prep1():
    eqn = Derivative(x*f(x), x, x, x)
    sol = Eq(f(x), (C1 + C2*x + C3*x**2) / x)
    assert checkodesol(eqn, sol, order=3, solve_for_func=False)[0]
    assert sol == dsolve(eqn, f(x), prep=True, hint='nth_algebraic')
    assert sol == dsolve(eqn, f(x))


# prep = True breaks this
def test_nth_algebraic_noprep2():
    eqn = Eq(Derivative(x*Derivative(f(x), x), x)/x, exp(x))
    sol = Eq(f(x), C1 + C2*log(x) + exp(x) - Ei(x))
    assert checkodesol(eqn, sol, order=2, solve_for_func=False)[0]
    assert sol == dsolve(eqn, f(x), prep=False, hint='nth_algebraic')


@XFAIL
def test_nth_algebraic_prep2():
    eqn = Eq(Derivative(x*Derivative(f(x), x), x)/x, exp(x))
    sol = Eq(f(x), C1 + C2*log(x) + exp(x) - Ei(x))
    assert checkodesol(eqn, sol, order=2, solve_for_func=False)[0]
    assert sol == dsolve(eqn, f(x), prep=True, hint='nth_algebraic')
    assert sol == dsolve(eqn, f(x))


def test_factorable():

    eq = f(x) + f(x)*f(x).diff(x)
    sols = [Eq(f(x), C1 - x), Eq(f(x), 0)]
    assert set(sols) == set(dsolve(eq, f(x), hint='factorable'))
    assert checkodesol(eq, sols) == 2*[(True, 0)]

    eq = f(x)*(f(x).diff(x)+f(x)*x+2)
    sols = [Eq(f(x), (C1 - sqrt(2)*sqrt(pi)*erfi(sqrt(2)*x/2))
            *exp(-x**2/2)), Eq(f(x), 0)]
    assert set(sols) == set(dsolve(eq, f(x), hint='factorable'))
    assert checkodesol(eq, sols) == 2*[(True, 0)]

    eq = (f(x).diff(x)+f(x)*x**2)*(f(x).diff(x, 2) + x*f(x))
    sols = [Eq(f(x), C1*airyai(-x) + C2*airybi(-x)),
            Eq(f(x), C1*exp(-x**3/3))]
    assert set(sols) == set(dsolve(eq, f(x), hint='factorable'))
    assert checkodesol(eq, sols[1]) == (True, 0)

    eq = (f(x).diff(x)+f(x)*x**2)*(f(x).diff(x, 2) + f(x))
    sols = [Eq(f(x), C1*exp(-x**3/3)), Eq(f(x), C1*sin(x) + C2*cos(x))]
    assert set(sols) == set(dsolve(eq, f(x), hint='factorable'))
    assert checkodesol(eq, sols) == 2*[(True, 0)]

    eq = (f(x).diff(x)**2-1)*(f(x).diff(x)**2-4)
    sols = [Eq(f(x), C1 - x), Eq(f(x), C1 + x), Eq(f(x), C1 + 2*x), Eq(f(x), C1 - 2*x)]
    assert set(sols) == set(dsolve(eq, f(x), hint='factorable'))
    assert checkodesol(eq, sols) == 4*[(True, 0)]

    eq = (f(x).diff(x, 2)-exp(f(x)))*f(x).diff(x)
    sol = Eq(f(x), C1)
    assert sol == dsolve(eq, f(x), hint='factorable')
    assert checkodesol(eq, sol) == (True, 0)

    eq = (f(x).diff(x)**2-1)*(f(x)*f(x).diff(x)-1)
    sol = [Eq(f(x), C1 - x), Eq(f(x), -sqrt(C1 + 2*x)),
           Eq(f(x), sqrt(C1 + 2*x)), Eq(f(x), C1 + x)]
    assert set(sol) == set(dsolve(eq, f(x), hint='factorable'))
    assert checkodesol(eq, sol) == 4*[(True, 0)]

    eq = Derivative(f(x), x)**4 - 2*Derivative(f(x), x)**2 + 1
    sol = [Eq(f(x), C1 - x), Eq(f(x), C1 + x)]
    assert set(sol) == set(dsolve(eq, f(x), hint='factorable'))
    assert checkodesol(eq, sol) == 2*[(True, 0)]

    eq = f(x)**2*Derivative(f(x), x)**6 - 2*f(x)**2*Derivative(f(x),
         x)**4 + f(x)**2*Derivative(f(x), x)**2 - 2*f(x)*Derivative(f(x),
         x)**5 + 4*f(x)*Derivative(f(x), x)**3 - 2*f(x)*Derivative(f(x),
         x) + Derivative(f(x), x)**4 - 2*Derivative(f(x), x)**2 + 1
    sol = [Eq(f(x), C1 - x), Eq(f(x), -sqrt(C1 + 2*x)),
           Eq(f(x), sqrt(C1 + 2*x)), Eq(f(x), C1 + x)]
    assert set(sol) == set(dsolve(eq, f(x), hint='factorable'))
    assert checkodesol(eq, sol) == 4*[(True, 0)]

    eq = (f(x).diff(x, 2)-exp(f(x)))*(f(x).diff(x, 2)+exp(f(x)))
    raises(NotImplementedError, lambda: dsolve(eq, hint = 'factorable'))

    eq = x**4*f(x)**2 + 2*x**4*f(x)*Derivative(f(x), (x, 2)) + x**4*Derivative(f(x),
         (x, 2))**2  + 2*x**3*f(x)*Derivative(f(x), x) + 2*x**3*Derivative(f(x),
         x)*Derivative(f(x), (x, 2)) - 7*x**2*f(x)**2 - 7*x**2*f(x)*Derivative(f(x),
         (x, 2)) + x**2*Derivative(f(x), x)**2 - 7*x*f(x)*Derivative(f(x), x) + 12*f(x)**2

    sol = [Eq(f(x), C1*besselj(2, x) + C2*bessely(2, x)), Eq(f(x), C1*besselj(sqrt(3),
           x) + C2*bessely(sqrt(3), x))]

    assert set(sol) == set(dsolve(eq, f(x), hint='factorable'))
    assert checkodesol(eq, sol) == 2*[(True, 0)]


def test_issue_15889():
    eq = exp(f(x).diff(x))-f(x)**2
    sol = Eq(NonElementaryIntegral(1/log(y**2), (y, f(x))), C1 + x)
    assert sol.dummy_eq(dsolve(eq))
    assert checkodesol(eq, sol) == (True, 0)

    eq = f(x).diff(x)**2 - f(x)**3
    sol = Eq(f(x), 4/(C1**2 - 2*C1*x + x**2))
    assert sol == dsolve(eq)
    assert checkodesol(eq, sol) == (True, 0)

    eq = f(x).diff(x)**2 - f(x)
    sol = Eq(f(x), C1**2/4 - C1*x/2 + x**2/4)
    assert sol == dsolve(eq)
    assert checkodesol(eq, sol) == (True, 0)

    eq = f(x).diff(x)**2 - f(x)**2
    sol = [Eq(f(x), C1*exp(x)), Eq(f(x), C1*exp(-x))]
    assert sol == dsolve(eq)
    assert checkodesol(eq, sol) == 2*[(True, 0)]

    eq = f(x).diff(x)**2 - f(x)**3
    sol = Eq(f(x), 4/(C1**2 - 2*C1*x + x**2))
    assert sol == dsolve(eq)
    assert checkodesol(eq, sol) == (True, 0)


def test_Riccati_special_minus2():
    # Type: Riccati special alpha = -2, a*dy/dx + b*y**2 + c*y/x +d/x**2
    eq = 2*f(x).diff(x) + f(x)**2 - f(x)/x + 3*x**(-2)
    sol = dsolve(eq, f(x), hint='Riccati_special_minus2')
    assert checkodesol(eq, sol, order=1, solve_for_func=False)[0]


def test_Bernoulli():
    # Type: Bernoulli, f'(x) + p(x)*f(x) == q(x)*f(x)**n
    eq = Eq(x*f(x).diff(x) + f(x) - f(x)**2, 0)
    sol = dsolve(eq, f(x), hint='Bernoulli')
    assert sol == Eq(f(x), 1/(C1*x + 1))
    assert checkodesol(eq, sol, order=1, solve_for_func=False)[0]

    eq = f(x).diff(x) - y*f(x)
    sol = dsolve(eq, hint='Bernoulli')
    assert sol == Eq(f(x), C1*exp(x*y))
    assert checkodesol(eq, sol)[0]

    eq = f(x)*f(x).diff(x) - 1
    sol = dsolve(eq,hint='Bernoulli')
    assert sol == [Eq(f(x), -sqrt(C1 + 2*x)), Eq(f(x), sqrt(C1 + 2*x))]
    assert checkodesol(eq, sol) == [(True, 0), (True, 0)]


def test_nth_order_linear_euler_eq_homogeneous():
    x, t, a, b, c = symbols('x t a b c')
    y = Function('y')
    our_hint = "nth_linear_euler_eq_homogeneous"

    eq = diff(f(t), t, 4)*t**4 - 13*diff(f(t), t, 2)*t**2 + 36*f(t)
    assert our_hint in classify_ode(eq)

    eq = a*y(t) + b*t*diff(y(t), t) + c*t**2*diff(y(t), t, 2)
    assert our_hint in classify_ode(eq)

    _ode_solver_test(_get_examples_ode_sol_euler_homogeneous())


def test_nth_order_linear_euler_eq_nonhomogeneous_undetermined_coefficients():
    x, t = symbols('x t')
    a, b, c, d = symbols('a b c d', integer=True)
    our_hint = "nth_linear_euler_eq_nonhomogeneous_undetermined_coefficients"

    eq = x**4*diff(f(x), x, 4) - 13*x**2*diff(f(x), x, 2) + 36*f(x) + x
    assert our_hint in classify_ode(eq, f(x))

    eq = a*x**2*diff(f(x), x, 2) + b*x*diff(f(x), x) + c*f(x) + d*log(x)
    assert our_hint in classify_ode(eq, f(x))

    _ode_solver_test(_get_examples_ode_sol_euler_undetermined_coeff())


def test_nth_order_linear_euler_eq_nonhomogeneous_variation_of_parameters():
    x, t = symbols('x, t')
    a, b, c, d = symbols('a, b, c, d', integer=True)
    our_hint = "nth_linear_euler_eq_nonhomogeneous_variation_of_parameters"

    eq = Eq(x**2*diff(f(x),x,2) - 8*x*diff(f(x),x) + 12*f(x), x**2)
    assert our_hint in classify_ode(eq, f(x))

    eq = Eq(a*x**3*diff(f(x),x,3) + b*x**2*diff(f(x),x,2) + c*x*diff(f(x),x) + d*f(x), x*log(x))
    assert our_hint in classify_ode(eq, f(x))

    _ode_solver_test(_get_examples_ode_sol_euler_var_para())


def _get_examples_ode_sol_euler_homogeneous():
    return {
            'hint': "nth_linear_euler_eq_homogeneous",
            'func': f(x),
            'examples':{
    'euler_hom_01': {
        'eq': Eq(-3*diff(f(x), x)*x + 2*x**2*diff(f(x), x, x), 0),
        'sol': [Eq(f(x), C1 + C2*x**Rational(5, 2))],
    },

    'euler_hom_02': {
        'eq': Eq(3*f(x) - 5*diff(f(x), x)*x + 2*x**2*diff(f(x), x, x), 0),
        'sol': [Eq(f(x), C1*sqrt(x) + C2*x**3)]
    },

    'euler_hom_03': {
        'eq': Eq(4*f(x) + 5*diff(f(x), x)*x + x**2*diff(f(x), x, x), 0),
        'sol': [Eq(f(x), (C1 + C2*log(x))/x**2)]
    },

    'euler_hom_04': {
        'eq': Eq(6*f(x) - 6*diff(f(x), x)*x + 1*x**2*diff(f(x), x, x) + x**3*diff(f(x), x, x, x), 0),
        'sol': [Eq(f(x), C1/x**2 + C2*x + C3*x**3)]
    },

    'euler_hom_05': {
        'eq': Eq(-125*f(x) + 61*diff(f(x), x)*x - 12*x**2*diff(f(x), x, x) + x**3*diff(f(x), x, x, x), 0),
        'sol': [Eq(f(x), x**5*(C1 + C2*log(x) + C3*log(x)**2))]
    },

    'euler_hom_06': {
        'eq': x**2*diff(f(x), x, 2) + x*diff(f(x), x) - 9*f(x),
        'sol': [Eq(f(x), C1*x**-3 + C2*x**3)]
    },

    'euler_hom_07': {
        'eq': sin(x)*x**2*f(x).diff(x, 2) + sin(x)*x*f(x).diff(x) + sin(x)*f(x),
        'sol': [Eq(f(x), C1*sin(log(x)) + C2*cos(log(x)))],
        'XFAIL': ['2nd_power_series_regular']
    },
    }
    }


def _get_examples_ode_sol_euler_undetermined_coeff():
    return {
            'hint': "nth_linear_euler_eq_nonhomogeneous_undetermined_coefficients",
            'func': f(x),
            'examples':{
    'euler_undet_01': {
        'eq': Eq(x**2*diff(f(x), x, x) + x*diff(f(x), x), 1),
        'sol': [Eq(f(x), C1 + C2*log(x) + log(x)**2/2)]
    },

    'euler_undet_02': {
        'eq': Eq(x**2*diff(f(x), x, x) - 2*x*diff(f(x), x) + 2*f(x), x**3),
        'sol': [Eq(f(x), x*(C1 + C2*x + Rational(1, 2)*x**2))]
    },

    'euler_undet_03': {
        'eq': Eq(x**2*diff(f(x), x, x) - x*diff(f(x), x) - 3*f(x), log(x)/x),
        'sol': [Eq(f(x), (C1 + C2*x**4 - log(x)**2/8 - log(x)/16)/x)]
    },

    'euler_undet_04': {
        'eq': Eq(x**2*diff(f(x), x, x) + 3*x*diff(f(x), x) - 8*f(x), log(x)**3 - log(x)),
        'sol': [Eq(f(x), C1/x**4 + C2*x**2 - Rational(1,8)*log(x)**3 - Rational(3,32)*log(x)**2 - Rational(1,64)*log(x) - Rational(7, 256))]
    },

    'euler_undet_05': {
        'eq': Eq(x**3*diff(f(x), x, x, x) - 3*x**2*diff(f(x), x, x) + 6*x*diff(f(x), x) - 6*f(x), log(x)),
        'sol': [Eq(f(x), C1*x + C2*x**2 + C3*x**3 - Rational(1, 6)*log(x) - Rational(11, 36))]
    },
    }
    }


def _get_examples_ode_sol_euler_var_para():
    return {
            'hint': "nth_linear_euler_eq_nonhomogeneous_variation_of_parameters",
            'func': f(x),
            'examples':{
    'euler_var_01': {
        'eq': Eq(x**2*Derivative(f(x), x, x) - 2*x*Derivative(f(x), x) + 2*f(x), x**4),
        'sol': [Eq(f(x), x*(C1 + C2*x + x**3/6))]
    },

    'euler_var_02': {
        'eq': Eq(3*x**2*diff(f(x), x, x) + 6*x*diff(f(x), x) - 6*f(x), x**3*exp(x)),
        'sol': [Eq(f(x), C1/x**2 + C2*x + x*exp(x)/3 - 4*exp(x)/3 + 8*exp(x)/(3*x) - 8*exp(x)/(3*x**2))]
    },

    'euler_var_03': {
        'eq': Eq(x**2*Derivative(f(x), x, x) - 2*x*Derivative(f(x), x) + 2*f(x), x**4*exp(x)),
        'sol':  [Eq(f(x), x*(C1 + C2*x + x*exp(x) - 2*exp(x)))]
    },

    'euler_var_04': {
        'eq': x**2*Derivative(f(x), x, x) - 2*x*Derivative(f(x), x) + 2*f(x) - log(x),
        'sol': [Eq(f(x), C1*x + C2*x**2 + log(x)/2 + Rational(3, 4))]
    },

    'euler_var_05': {
        'eq': -exp(x) + (x*Derivative(f(x), (x, 2)) + Derivative(f(x), x))/x,
        'sol': [Eq(f(x), C1 + C2*log(x) + exp(x) - Ei(x))]
    },
    }
    }


def _get_all_examples():
    all_solvers = [_get_examples_ode_sol_euler_homogeneous(), _get_examples_ode_sol_euler_undetermined_coeff(), _get_examples_ode_sol_euler_var_para()]
    all_examples = []
    for solver in all_solvers:
        for example in solver['examples']:
            temp = {
                'hint': solver['hint'],
                'func': solver['examples'][example].get('func',solver['func']),
                'eq': solver['examples'][example]['eq'],
                'sol': solver['examples'][example]['sol'],
                'XFAIL': solver['examples'][example].get('XFAIL',[]),
                'example_name': example,
            }
            all_examples.append(temp)
    return all_examples
