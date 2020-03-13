#
# The main tests for the code in single.py are currently located in
# sympy/solvers/tests/test_ode.py
#
from sympy import (cos, Derivative, diff,
    Eq, exp, log, pi, Rational, sin,
    sqrt, symbols, Ei, erfi)

from sympy.core import Function, Symbol
from sympy.functions import airyai, airybi, besselj, bessely
from sympy.integrals.risch import NonElementaryIntegral
from sympy.solvers.ode import dsolve
from sympy.solvers.ode.ode import _remove_redundant_solutions
from sympy.solvers.ode.single import (FirstLinear, ODEMatchError,
    SingleODEProblem, SingleODESolver)

from sympy.solvers.ode.subscheck import checkodesol

from sympy.testing.pytest import XFAIL, raises


x = Symbol('x')
y = Symbol('y')
f = Function('f')
g = Function('g')
C1, C2, C3 = symbols('C1:4')


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
    assert sol == Eq(f(x), 1/(x*(C1 + 1/x)))
    assert checkodesol(eq, sol, order=1, solve_for_func=False)[0]
