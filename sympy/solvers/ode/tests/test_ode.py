from sympy import (acos, acosh, asinh, atan, cos, Derivative, diff,
    Dummy, Eq, Ne, exp, Function, I, Integral, LambertW, log, O, pi,
    Rational, rootof, S, sin, sqrt, Subs, Symbol, tan, asin, sinh,
    Piecewise, symbols, Poly, sec, Ei, re, im, atan2, collect, hyper, simplify, integrate)
from sympy.solvers.ode import (classify_ode,
    homogeneous_order, infinitesimals, checkinfsol,
    dsolve)

from sympy.solvers.ode.subscheck import checkodesol, checksysodesol
from sympy.solvers.ode.ode import (_linear_coeff_match,
    _undetermined_coefficients_match, classify_sysode,
    constant_renumber, constantsimp, get_numbered_constants, solve_ics)

from sympy.functions import airyai, airybi, besselj, bessely

from sympy.solvers.deutils import ode_order
from sympy.testing.pytest import XFAIL, skip, raises, slow, ON_TRAVIS, SKIP
from sympy.utilities.misc import filldedent

C0, C1, C2, C3, C4, C5, C6, C7, C8, C9, C10 = symbols('C0:11')
u, x, y, z = symbols('u,x:z', real=True)
f = Function('f')
g = Function('g')
h = Function('h')

# Note: the tests below may fail (but still be correct) if ODE solver,
# the integral engine, solve(), or even simplify() changes. Also, in
# differently formatted solutions, the arbitrary constants might not be
# equal.  Using specific hints in tests can help to avoid this.

# Tests of order higher than 1 should run the solutions through
# constant_renumber because it will normalize it (constant_renumber causes
# dsolve() to return different results on different machines)


def test_get_numbered_constants():
    with raises(ValueError):
        get_numbered_constants(None)


def test_dsolve_system():

    eqs = [f(x).diff(x, 2), g(x).diff(x)]
    with raises(ValueError):
        dsolve(eqs) # NotImplementedError would be better

    eqs = [f(x).diff(x) - x, f(x).diff(x) + x]
    with raises(ValueError):
        # Could also be NotImplementedError. f(x)=0 is a solution...
        dsolve(eqs)

    eqs = [f(x, y).diff(x)]
    with raises(ValueError):
        dsolve(eqs)

    eqs = [f(x, y).diff(x)+g(x).diff(x), g(x).diff(x)]
    with raises(ValueError):
        dsolve(eqs)


def test_dsolve_all_hint():
    eq = f(x).diff(x)
    output = dsolve(eq, hint='all')

    # Match the Dummy variables:
    sol1 = output['separable_Integral']
    _y = sol1.lhs.args[1][0]
    sol1 = output['1st_homogeneous_coeff_subs_dep_div_indep_Integral']
    _u1 = sol1.rhs.args[1].args[1][0]

    expected = {'Bernoulli_Integral': Eq(f(x), C1 + Integral(0, x)),
        '1st_homogeneous_coeff_best': Eq(f(x), C1),
        'Bernoulli': Eq(f(x), C1),
        'nth_algebraic': Eq(f(x), C1),
        'nth_linear_euler_eq_homogeneous': Eq(f(x), C1),
        'nth_linear_constant_coeff_homogeneous': Eq(f(x), C1),
        'separable': Eq(f(x), C1),
        '1st_homogeneous_coeff_subs_indep_div_dep': Eq(f(x), C1),
        'nth_algebraic_Integral': Eq(f(x), C1),
        '1st_linear': Eq(f(x), C1),
        '1st_linear_Integral': Eq(f(x), C1 + Integral(0, x)),
        'lie_group': Eq(f(x), C1),
        '1st_homogeneous_coeff_subs_dep_div_indep': Eq(f(x), C1),
        '1st_homogeneous_coeff_subs_dep_div_indep_Integral': Eq(log(x), C1 + Integral(-1/_u1, (_u1, f(x)/x))),
        '1st_power_series': Eq(f(x), C1),
        'separable_Integral': Eq(Integral(1, (_y, f(x))), C1 + Integral(0, x)),
        '1st_homogeneous_coeff_subs_indep_div_dep_Integral': Eq(f(x), C1),
        'best': Eq(f(x), C1),
        'best_hint': 'nth_algebraic',
        'default': 'nth_algebraic',
        'order': 1}
    assert output == expected

    assert dsolve(eq, hint='best') == Eq(f(x), C1)


def test_dsolve_ics():
    # Maybe this should just use one of the solutions instead of raising...
    with raises(NotImplementedError):
        dsolve(f(x).diff(x) - sqrt(f(x)), ics={f(1):1})


@slow
@XFAIL
def test_nonlinear_3eq_order1_type1():
    if ON_TRAVIS:
        skip("Too slow for travis.")
    a, b, c = symbols('a b c')

    eqs = [
        a * f(x).diff(x) - (b - c) * g(x) * h(x),
        b * g(x).diff(x) - (c - a) * h(x) * f(x),
        c * h(x).diff(x) - (a - b) * f(x) * g(x),
    ]

    assert dsolve(eqs) # NotImplementedError


def test_dsolve_euler_rootof():
    eq = x**6 * f(x).diff(x, 6) - x*f(x).diff(x) + f(x)
    sol = Eq(f(x),
          C1*x
        + C2*x**rootof(x**5 - 14*x**4 + 71*x**3 - 154*x**2 + 120*x - 1, 0)
        + C3*x**rootof(x**5 - 14*x**4 + 71*x**3 - 154*x**2 + 120*x - 1, 1)
        + C4*x**rootof(x**5 - 14*x**4 + 71*x**3 - 154*x**2 + 120*x - 1, 2)
        + C5*x**rootof(x**5 - 14*x**4 + 71*x**3 - 154*x**2 + 120*x - 1, 3)
        + C6*x**rootof(x**5 - 14*x**4 + 71*x**3 - 154*x**2 + 120*x - 1, 4)
        )
    assert dsolve(eq) == sol


def test_linear_2eq_order1_type2_noninvertible():
    # a*d - b*c == 0
    eqs = [Eq(diff(f(x), x), f(x) + g(x) + 5),
           Eq(diff(g(x), x), f(x) + g(x) + 7)]
    sol = [Eq(f(x), C1*exp(2*x) + C2 - x - 3), Eq(g(x), C1*exp(2*x) - C2 + x - 3)]
    assert dsolve(eqs) == sol
    assert checksysodesol(eqs, sol) == (True, [0, 0])


def test_linear_2eq_order1_type2_fixme():
    # There is a FIXME comment about this in the code that handles this case.
    # The answer returned is currently incorrect as reported by checksysodesol
    # below in the XFAIL test below...

    # a*d - b*c == 0 and a + b*c/a = 0
    eqs = [Eq(diff(f(x), x),  f(x) + g(x) + 5),
           Eq(diff(g(x), x), -f(x) - g(x) + 7)]
    sol = [Eq(f(x), C1 + C2*(x + 1) + 12*x**2 + 5*x), Eq(g(x), -C1 - C2*x - 12*x**2 + 7*x)]
    assert dsolve(eqs) == sol
    # FIXME: checked in XFAIL test_linear_2eq_order1_type2_fixme_check below


@XFAIL
def test_linear_2eq_order1_type2_fixme_check():
    # See test_linear_2eq_order1_type2_fixme above
    eqs = [Eq(diff(f(x), x),  f(x) + g(x) + 5),
           Eq(diff(g(x), x), -f(x) - g(x) + 7)]
    sol = [Eq(f(x), C1 + C2*(x + 1) + 12*x**2 + 5*x), Eq(g(x), -C1 - C2*x - 12*x**2 + 7*x)]
    assert checksysodesol(eqs, sol) == (True, [0, 0])


def test_linear_2eq_order1_type4():
    eqs = [Eq(diff(f(x), x),  f(x) + x*g(x)),
           Eq(diff(g(x), x),-x*f(x) + g(x))]
    sol = [Eq(f(x), (C1*cos(x**2/2) + C2*sin(x**2/2))*exp(x)),
           Eq(g(x), (-C1*sin(x**2/2) + C2*cos(x**2/2))*exp(x))]
    dsolve_sol = dsolve(eqs)
    # FIXME: This should probably be fixed so that this happens in the solver:
    dsolve_sol = [s.doit() for s in sol]
    assert dsolve_sol == sol
    assert checksysodesol(eqs, sol) == (True, [0, 0])


def test_linear_2eq_order1_type4_broken():
    eqs = [Eq(f(x).diff(x), f(x) + x*g(x)),
           Eq(g(x).diff(x), x*f(x) - g(x))]
    # FIXME: This is not the correct solution:
    sol = [Eq(f(x), (C1*sin(x) + C2*cos(x))*exp(x**2/2)),
           Eq(g(x), (C1*cos(x) - C2*sin(x))*exp(x**2/2))]
    dsolve_sol = dsolve(eqs)
    # FIXME: This should probably be fixed so that this happens in the solver:
    dsolve_sol = [s.doit() for s in sol]
    assert dsolve_sol == sol
    # FIXME: Checked in XFAIL test_linear_2eq_order1_type4_broken_check below


@XFAIL
def test_linear_2eq_order1_type4_broken_check():
    # see test_linear_2eq_order1_type4_broken above
    eqs = [Eq(f(x).diff(x), f(x) + x*g(x)),
           Eq(g(x).diff(x), x*f(x) - g(x))]
    # FIXME: This is not the correct solution:
    sol = [Eq(f(x), (C1*sin(x) + C2*cos(x))*exp(x**2/2)),
           Eq(g(x), (C1*cos(x) - C2*sin(x))*exp(x**2/2))]
    assert checksysodesol(eqs, sol) == (True, [0, 0])


def test_linear_2eq_order1_type5():
    eqs = [Eq(diff(f(x), x),  x*f(x) + x**2*g(x)),
           Eq(diff(g(x), x),  2*x**2*f(x) + (x + 3*x**2)*g(x))]
    sol = [
        Eq(f(x), (C1*(-sqrt(17)/4 - Rational(3,4))*exp((Rational(3,2) - sqrt(17)/2)*Integral(x**2, x)) + C2*(Rational(-3,4) + sqrt(17)/4) *
                exp((Rational(3,2) + sqrt(17)/2)*Integral(x**2, x)))*exp(Integral(x, x))),
        Eq(g(x), (C1*exp((Rational(3,2) - sqrt(17)/2)*Integral(x**2, x)) + C2*exp((Rational(3,2) + sqrt(17)/2)*Integral(x**2, x))) *
                exp(Integral(x, x)))]

    assert dsolve(eqs) == sol
    assert checksysodesol(eqs, sol) == (True, [0, 0])


def test_linear_2eq_order1_type6_path1_broken():
    eqs = [Eq(diff(f(x), x), f(x) + x*g(x)),
           Eq(diff(g(x), x), 2*(1 + 2/x)*f(x) + 2*(x - 1/x) * g(x))]
    # FIXME: This is not the correct solution:
    sol = [
        Eq(f(x), (C1 + Integral(C2*x*exp(-2*Integral(1/x, x))*exp(Integral(-2*x - 1, x)), x))*exp(-Integral(-2*x - 1, x))),
        Eq(g(x), C1*exp(-2*Integral(1/x, x))
            + 2*(C1 + Integral(C2*x*exp(-2*Integral(1/x, x))*exp(Integral(-2*x - 1, x)), x))*exp(-Integral(-2*x - 1, x)))
        ]
    dsolve_sol = dsolve(eqs)
    # FIXME: Comparing solutions with == doesn't work in this case...
    assert [ds.lhs for ds in dsolve_sol] == [f(x), g(x)]
    assert [ds.rhs.equals(ss.rhs) for ds, ss in zip(dsolve_sol, sol)]
    # FIXME: checked in XFAIL test_linear_2eq_order1_type6_path1_broken_check below


@XFAIL
def test_linear_2eq_order1_type6_path1_broken_check():
    # See test_linear_2eq_order1_type6_path1_broken above
    eqs = [Eq(diff(f(x), x), f(x) + x*g(x)),
           Eq(diff(g(x), x), 2*(1 + 2/x)*f(x) + 2*(x - 1/x) * g(x))]
    # FIXME: This is not the correct solution:
    sol = [
        Eq(f(x), (C1 + Integral(C2*x*exp(-2*Integral(1/x, x))*exp(Integral(-2*x - 1, x)), x))*exp(-Integral(-2*x - 1, x))),
        Eq(g(x), C1*exp(-2*Integral(1/x, x))
            + 2*(C1 + Integral(C2*x*exp(-2*Integral(1/x, x))*exp(Integral(-2*x - 1, x)), x))*exp(-Integral(-2*x - 1, x)))
        ]
    assert checksysodesol(eqs, sol) == (True, [0, 0])  # XFAIL


def test_linear_2eq_order1_type6_path2_broken():
    # This is the reverse of the equations above and should also be handled by
    # type6.
    eqs = [Eq(diff(g(x), x), 2*(1 + 2/x)*g(x) + 2*(x - 1/x) * f(x)),
           Eq(diff(f(x), x), g(x) + x*f(x))]
    # FIXME: This is not the correct solution:
    sol = [
        Eq(g(x), C1*exp(-2*Integral(1/x, x))
             + 2*(C1 + Integral(-C2*exp(-2*Integral(1/x, x))*exp(Integral(-2*x - 1, x)), x))*exp(-Integral(-2*x - 1, x))),
        Eq(f(x), (C1 + Integral(-C2*exp(-2*Integral(1/x, x))*exp(Integral(-2*x - 1, x)), x))*exp(-Integral(-2*x - 1, x)))
    ]
    dsolve_sol = dsolve(eqs)
    # Comparing solutions with == doesn't work in this case...
    assert [ds.lhs for ds in dsolve_sol] == [g(x), f(x)]
    assert [ds.rhs.equals(ss.rhs) for ds, ss in zip(dsolve_sol, sol)]
    # FIXME: checked in XFAIL test_linear_2eq_order1_type6_path2_broken_check below


@XFAIL
def test_linear_2eq_order1_type6_path2_broken_check():
    # See test_linear_2eq_order1_type6_path2_broken above
    eqs = [Eq(diff(g(x), x), 2*(1 + 2/x)*g(x) + 2*(x - 1/x) * f(x)),
           Eq(diff(f(x), x), g(x) + x*f(x))]
    sol = [
        Eq(g(x), C1*exp(-2*Integral(1/x, x))
             + 2*(C1 + Integral(-C2*exp(-2*Integral(1/x, x))*exp(Integral(-2*x - 1, x)), x))*exp(-Integral(-2*x - 1, x))),
        Eq(f(x), (C1 + Integral(-C2*exp(-2*Integral(1/x, x))*exp(Integral(-2*x - 1, x)), x))*exp(-Integral(-2*x - 1, x)))
    ]
    assert checksysodesol(eqs, sol) == (True, [0, 0])  # XFAIL


def test_nth_euler_imroot():
    eq = x**2 * f(x).diff(x, 2) + x * f(x).diff(x) + 4 * f(x) - 1/x
    sol = Eq(f(x), C1*sin(2*log(x)) + C2*cos(2*log(x)) + 1/(5*x))
    dsolve_sol = dsolve(eq, hint='nth_linear_euler_eq_nonhomogeneous_variation_of_parameters')
    assert dsolve_sol == sol
    assert checkodesol(eq, sol, order=2, solve_for_func=False)[0]


def test_constant_coeff_circular_atan2():
    eq = f(x).diff(x, x) + y*f(x)
    sol = Eq(f(x), C1*exp(-x*sqrt(-y)) + C2*exp(x*sqrt(-y)))
    assert dsolve(eq) == sol
    assert checkodesol(eq, sol, order=2, solve_for_func=False)[0]


def test_linear_2eq_order2_type1_broken1():
    eqs = [Eq(f(x).diff(x, 2), 2*f(x) + g(x)),
           Eq(g(x).diff(x, 2), -f(x))]
    # FIXME: This is not the correct solution:
    sol = [
        Eq(f(x), 2*C1*(x + 2)*exp(x) + 2*C2*(x + 2)*exp(-x) + 2*C3*x*exp(x) + 2*C4*x*exp(-x)),
        Eq(g(x), -2*C1*x*exp(x) - 2*C2*x*exp(-x) + C3*(-2*x + 4)*exp(x) + C4*(-2*x - 4)*exp(-x))
    ]
    assert dsolve(eqs) == sol
    # FIXME: checked in XFAIL test_linear_2eq_order2_type1_broken1_check below


@XFAIL
def test_linear_2eq_order2_type1_broken1_check():
    # See test_linear_2eq_order2_type1_broken1 above
    eqs = [Eq(f(x).diff(x, 2), 2*f(x) + g(x)),
           Eq(g(x).diff(x, 2), -f(x))]
    # This is the returned solution but it isn't correct:
    sol = [
        Eq(f(x), 2*C1*(x + 2)*exp(x) + 2*C2*(x + 2)*exp(-x) + 2*C3*x*exp(x) + 2*C4*x*exp(-x)),
        Eq(g(x), -2*C1*x*exp(x) - 2*C2*x*exp(-x) + C3*(-2*x + 4)*exp(x) + C4*(-2*x - 4)*exp(-x))
    ]
    assert checksysodesol(eqs, sol) == (True, [0, 0])


@XFAIL
def test_linear_2eq_order2_type1_broken2():
    eqs = [Eq(f(x).diff(x, 2), 0),
            Eq(g(x).diff(x, 2), f(x))]
    sol = [
        Eq(f(x), C1 + C2*x),
        Eq(g(x), C4 + C3*x + C2*x**3/6 + C1*x**2/2)
    ]
    assert dsolve(eqs) == sol # UnboundLocalError


def test_linear_2eq_order2_type1_broken2_check():
    eqs = [Eq(f(x).diff(x, 2), 0),
            Eq(g(x).diff(x, 2), f(x))]
    sol = [
        Eq(f(x), C1 + C2*x),
        Eq(g(x), C4 + C3*x + C2*x**3/6 + C1*x**2/2)
    ]
    assert checksysodesol(eqs, sol) == (True, [0, 0])


def test_linear_2eq_order2_type1():
    eqs = [Eq(f(x).diff(x, 2), 2*f(x)),
           Eq(g(x).diff(x, 2), -f(x) + 2*g(x))]
    sol = [
        Eq(f(x), 2*sqrt(2)*C1*exp(sqrt(2)*x) + 2*sqrt(2)*C2*exp(-sqrt(2)*x)),
        Eq(g(x), -C1*x*exp(sqrt(2)*x) + C2*x*exp(-sqrt(2)*x) + C3*exp(sqrt(2)*x) + C4*exp(-sqrt(2)*x))
    ]
    assert dsolve(eqs) == sol
    assert checksysodesol(eqs, sol) == (True, [0, 0])

    eqs = [Eq(f(x).diff(x, 2), 2*f(x) + g(x)),
           Eq(g(x).diff(x, 2), + 2*g(x))]
    sol = [
        Eq(f(x), C1*x*exp(sqrt(2)*x) - C2*x*exp(-sqrt(2)*x) + C3*exp(sqrt(2)*x) + C4*exp(-sqrt(2)*x)),
        Eq(g(x), 2*sqrt(2)*C1*exp(sqrt(2)*x) + 2*sqrt(2)*C2*exp(-sqrt(2)*x))
    ]
    assert dsolve(eqs) == sol
    assert checksysodesol(eqs, sol) == (True, [0, 0])

    eqs = [Eq(f(x).diff(x, 2), f(x)),
           Eq(g(x).diff(x, 2), f(x))]
    sol = [Eq(f(x), C1*exp(x) + C2*exp(-x)),
           Eq(g(x), C1*exp(x) + C2*exp(-x) - C3*x - C4)]
    assert dsolve(eqs) == sol
    assert checksysodesol(eqs, sol) == (True, [0, 0])

    eqs = [Eq(f(x).diff(x, 2), f(x) + g(x)),
           Eq(g(x).diff(x, 2), -f(x) - g(x))]
    sol = [Eq(f(x), C1*x**3 + C2*x**2 + C3*x + C4),
           Eq(g(x), -C1*x**3 + 6*C1*x - C2*x**2 + 2*C2 - C3*x - C4)]
    assert dsolve(eqs) == sol
    assert checksysodesol(eqs, sol) == (True, [0, 0])


def test_linear_2eq_order2_type2():
    eqs = [Eq(f(x).diff(x, 2), f(x) + g(x) + 1),
           Eq(g(x).diff(x, 2), f(x) + g(x) + 1)]
    sol = [Eq(f(x), C1*exp(sqrt(2)*x) + C2*exp(-sqrt(2)*x) + C3*x + C4 - S.Half),
           Eq(g(x), C1*exp(sqrt(2)*x) + C2*exp(-sqrt(2)*x) - C3*x - C4 - S.Half)]
    assert dsolve(eqs) == sol
    assert checksysodesol(eqs, sol) == (True, [0, 0])

    eqs = [Eq(f(x).diff(x, 2), f(x) + g(x) + 1),
           Eq(g(x).diff(x, 2), -f(x) - g(x) + 1)]
    sol = [Eq(f(x), C1*x**3 + C2*x**2 + C3*x + C4 + x**4/12 + x**2/2),
           Eq(g(x), -C1*x**3 + 6*C1*x - C2*x**2 + 2*C2 - C3*x - C4 - x**4/12 + x**2/2)]
    assert dsolve(eqs) == sol
    assert checksysodesol(eqs, sol) == (True, [0, 0])


def test_linear_2eq_order2_type4_broken():
    Ca, Cb, Ra, Rb = symbols('Ca, Cb, Ra, Rb')
    eq = [f(x).diff(x, 2) + 2*f(x).diff(x) + f(x) + g(x) - 2*exp(I*x),
          g(x).diff(x, 2) + 2*g(x).diff(x) + f(x) + g(x) - 2*exp(I*x)]
    # FIXME: This is not the correct solution:
    # Solution returned with Ca, Ra etc symbols is clearly incorrect:
    sol = [
        Eq(f(x), C1 + C2*exp(2*x) + C3*exp(x*(1 + sqrt(3))) + C4*exp(x*(-sqrt(3) + 1)) + (I*Ca + Ra)*exp(I*x)),
        Eq(g(x), -C1 - 3*C2*exp(2*x) + C3*(-3*sqrt(3) - 4 + (1 + sqrt(3))**2)*exp(x*(1 + sqrt(3)))
            + C4*(-4 + (-sqrt(3) + 1)**2 + 3*sqrt(3))*exp(x*(-sqrt(3) + 1)) + (I*Cb + Rb)*exp(I*x))
    ]
    dsolve_sol = dsolve(eq)
    assert dsolve_sol == sol
    # FIXME: checked in XFAIL test_linear_2eq_order2_type4_broken_check below


@XFAIL
def test_linear_2eq_order2_type4_broken_check():
    # See test_linear_2eq_order2_type4_broken above
    Ca, Cb, Ra, Rb = symbols('Ca, Cb, Ra, Rb')
    eq = [f(x).diff(x, 2) + 2*f(x).diff(x) + f(x) + g(x) - 2*exp(I*x),
          g(x).diff(x, 2) + 2*g(x).diff(x) + f(x) + g(x) - 2*exp(I*x)]
    # Solution returned with Ca, Ra etc symbols is clearly incorrect:
    sol = [
        Eq(f(x), C1 + C2*exp(2*x) + C3*exp(x*(1 + sqrt(3))) + C4*exp(x*(-sqrt(3) + 1)) + (I*Ca + Ra)*exp(I*x)),
        Eq(g(x), -C1 - 3*C2*exp(2*x) + C3*(-3*sqrt(3) - 4 + (1 + sqrt(3))**2)*exp(x*(1 + sqrt(3)))
            + C4*(-4 + (-sqrt(3) + 1)**2 + 3*sqrt(3))*exp(x*(-sqrt(3) + 1)) + (I*Cb + Rb)*exp(I*x))
    ]
    assert checksysodesol(eq, sol) == (True, [0, 0]) # Fails here


def test_linear_2eq_order2_type5():
    eqs = [Eq(f(x).diff(x, 2), 2*(x*g(x).diff(x) - g(x))),
           Eq(g(x).diff(x, 2),-2*(x*f(x).diff(x) - f(x)))]
    sol = [Eq(f(x), C3*x + x*Integral((2*C1*cos(x**2) + 2*C2*sin(x**2))/x**2, x)),
           Eq(g(x), C4*x + x*Integral((-2*C1*sin(x**2) + 2*C2*cos(x**2))/x**2, x))]
    assert dsolve(eqs) == sol
    assert checksysodesol(eqs, sol) == (True, [0, 0])


def test_linear_2eq_order2_type8():
    eqs = [Eq(f(x).diff(x, 2), 2/x *(x*g(x).diff(x) - g(x))),
           Eq(g(x).diff(x, 2),-2/x *(x*f(x).diff(x) - f(x)))]
    # FIXME: This is what is returned but it does not seem correct:
    sol = [Eq(f(x), C3*x + x*Integral((-C1*cos(Integral(-2, x)) - C2*sin(Integral(-2, x)))/x**2, x)),
           Eq(g(x), C4*x + x*Integral((-C1*sin(Integral(-2, x)) + C2*cos(Integral(-2, x)))/x**2, x))]
    assert dsolve(eqs) == sol
    assert checksysodesol(eqs, sol) == (True, [0, 0]) # Fails here


@XFAIL
def test_nonlinear_3eq_order1_type4():
    eqs = [
        Eq(f(x).diff(x), (2*h(x)*g(x) - 3*g(x)*h(x))),
        Eq(g(x).diff(x), (4*f(x)*h(x) - 2*h(x)*f(x))),
        Eq(h(x).diff(x), (3*g(x)*f(x) - 4*f(x)*g(x))),
    ]
    dsolve(eqs)  # KeyError when matching
    # sol = ?
    # assert dsolve_sol == sol
    # assert checksysodesol(eqs, dsolve_sol) == (True, [0, 0, 0])


@slow
@XFAIL
def test_nonlinear_3eq_order1_type3():
    if ON_TRAVIS:
        skip("Too slow for travis.")
    eqs = [
        Eq(f(x).diff(x), (2*f(x)**2 - 3        )),
        Eq(g(x).diff(x), (4         - 2*h(x)   )),
        Eq(h(x).diff(x), (3*h(x)    - 4*f(x)**2)),
    ]
    dsolve(eqs)  # Not sure if this finishes...
    # sol = ?
    # assert dsolve_sol == sol
    # assert checksysodesol(eqs, dsolve_sol) == (True, [0, 0, 0])


@XFAIL
def test_nonlinear_3eq_order1_type5():
    eqs = [
        Eq(f(x).diff(x), f(x)*(2*f(x) - 3*g(x))),
        Eq(g(x).diff(x), g(x)*(4*g(x) - 2*h(x))),
        Eq(h(x).diff(x), h(x)*(3*h(x) - 4*f(x))),
    ]
    dsolve(eqs)  # KeyError
    # sol = ?
    # assert dsolve_sol == sol
    # assert checksysodesol(eqs, dsolve_sol) == (True, [0, 0, 0])


def test_linear_2eq_order1():
    x, y, z = symbols('x, y, z', cls=Function)
    k, l, m, n = symbols('k, l, m, n', Integer=True)
    t = Symbol('t')
    x0, y0 = symbols('x0, y0', cls=Function)

    eq1 = (Eq(diff(x(t),t), x(t) + y(t) + 9), Eq(diff(y(t),t), 2*x(t) + 5*y(t) + 23))
    sol1 = [Eq(x(t), C1*exp(t*(sqrt(6) + 3)) + C2*exp(t*(-sqrt(6) + 3)) - Rational(22, 3)), \
    Eq(y(t), C1*(2 + sqrt(6))*exp(t*(sqrt(6) + 3)) + C2*(-sqrt(6) + 2)*exp(t*(-sqrt(6) + 3)) - Rational(5, 3))]
    assert checksysodesol(eq1, sol1) == (True, [0, 0])

    eq2 = (Eq(diff(x(t),t), x(t) + y(t) + 81), Eq(diff(y(t),t), -2*x(t) + y(t) + 23))
    sol2 = [Eq(x(t), (C1*cos(sqrt(2)*t) + C2*sin(sqrt(2)*t))*exp(t) - Rational(58, 3)), \
    Eq(y(t), (-sqrt(2)*C1*sin(sqrt(2)*t) + sqrt(2)*C2*cos(sqrt(2)*t))*exp(t) - Rational(185, 3))]
    assert checksysodesol(eq2, sol2) == (True, [0, 0])

    eq3 = (Eq(diff(x(t),t), 5*t*x(t) + 2*y(t)), Eq(diff(y(t),t), 2*x(t) + 5*t*y(t)))
    sol3 = [Eq(x(t), (C1*exp(2*t) + C2*exp(-2*t))*exp(Rational(5, 2)*t**2)), \
    Eq(y(t), (C1*exp(2*t) - C2*exp(-2*t))*exp(Rational(5, 2)*t**2))]
    assert checksysodesol(eq3, sol3) == (True, [0, 0])

    eq4 = (Eq(diff(x(t),t), 5*t*x(t) + t**2*y(t)), Eq(diff(y(t),t), -t**2*x(t) + 5*t*y(t)))
    sol4 = [Eq(x(t), (C1*cos((t**3)/3) + C2*sin((t**3)/3))*exp(Rational(5, 2)*t**2)), \
    Eq(y(t), (-C1*sin((t**3)/3) + C2*cos((t**3)/3))*exp(Rational(5, 2)*t**2))]
    assert checksysodesol(eq4, sol4) == (True, [0, 0])

    eq5 = (Eq(diff(x(t),t), 5*t*x(t) + t**2*y(t)), Eq(diff(y(t),t), -t**2*x(t) + (5*t+9*t**2)*y(t)))
    sol5 = [Eq(x(t), (C1*exp((sqrt(77)/2 + Rational(9, 2))*(t**3)/3) + \
    C2*exp((-sqrt(77)/2 + Rational(9, 2))*(t**3)/3))*exp(Rational(5, 2)*t**2)), \
    Eq(y(t), (C1*(sqrt(77)/2 + Rational(9, 2))*exp((sqrt(77)/2 + Rational(9, 2))*(t**3)/3) + \
    C2*(-sqrt(77)/2 + Rational(9, 2))*exp((-sqrt(77)/2 + Rational(9, 2))*(t**3)/3))*exp(Rational(5, 2)*t**2))]
    assert checksysodesol(eq5, sol5) == (True, [0, 0])

    eq6 = (Eq(diff(x(t),t), 5*t*x(t) + t**2*y(t)), Eq(diff(y(t),t), (1-t**2)*x(t) + (5*t+9*t**2)*y(t)))
    sol6 = [Eq(x(t), C1*x0(t) + C2*x0(t)*Integral(t**2*exp(Integral(5*t, t))*exp(Integral(9*t**2 + 5*t, t))/x0(t)**2, t)), \
    Eq(y(t), C1*y0(t) + C2*(y0(t)*Integral(t**2*exp(Integral(5*t, t))*exp(Integral(9*t**2 + 5*t, t))/x0(t)**2, t) + \
    exp(Integral(5*t, t))*exp(Integral(9*t**2 + 5*t, t))/x0(t)))]
    s = dsolve(eq6)
    assert s == sol6   # too complicated to test with subs and simplify
    # assert checksysodesol(eq10, sol10) == (True, [0, 0])  # this one fails


def test_linear_2eq_order1_nonhomog_linear():
    e = [Eq(diff(f(x), x), f(x) + g(x) + 5*x),
         Eq(diff(g(x), x), f(x) - g(x))]
    raises(NotImplementedError, lambda: dsolve(e))


def test_linear_2eq_order1_nonhomog():
    # Note: once implemented, add some tests esp. with resonance
    e = [Eq(diff(f(x), x), f(x) + exp(x)),
         Eq(diff(g(x), x), f(x) + g(x) + x*exp(x))]
    raises(NotImplementedError, lambda: dsolve(e))


def test_linear_2eq_order1_type2_degen():
    e = [Eq(diff(f(x), x), f(x) + 5),
         Eq(diff(g(x), x), f(x) + 7)]
    s1 = [Eq(f(x), C1*exp(x) - 5), Eq(g(x), C1*exp(x) - C2 + 2*x - 5)]
    assert checksysodesol(e, s1) == (True, [0, 0])


def test_dsolve_linear_2eq_order1_diag_triangular():
    e = [Eq(diff(f(x), x), f(x)),
         Eq(diff(g(x), x), g(x))]
    s1 = [Eq(f(x), C1*exp(x)), Eq(g(x), C2*exp(x))]
    assert checksysodesol(e, s1) == (True, [0, 0])

    e = [Eq(diff(f(x), x), 2*f(x)),
         Eq(diff(g(x), x), 3*f(x) + 7*g(x))]
    s1 = [Eq(f(x), -5*C2*exp(2*x)),
          Eq(g(x), 5*C1*exp(7*x) + 3*C2*exp(2*x))]
    assert checksysodesol(e, s1) == (True, [0, 0])


def test_sysode_linear_2eq_order1_type1_D_lt_0():
    e = [Eq(diff(f(x), x), -9*I*f(x) - 4*g(x)),
         Eq(diff(g(x), x), -4*I*g(x))]
    s1 = [Eq(f(x), -4*C1*exp(-4*I*x) - 4*C2*exp(-9*I*x)), \
    Eq(g(x), 5*I*C1*exp(-4*I*x))]
    assert checksysodesol(e, s1) == (True, [0, 0])



def test_sysode_linear_2eq_order1_type1_D_lt_0_b_eq_0():
    e = [Eq(diff(f(x), x), -9*I*f(x)),
         Eq(diff(g(x), x), -4*I*g(x))]
    s1 = [Eq(f(x), -5*I*C2*exp(-9*I*x)), Eq(g(x), 5*I*C1*exp(-4*I*x))]
    assert checksysodesol(e, s1) == (True, [0, 0])



def test_sysode_linear_2eq_order1_many_zeros():
    t = Symbol('t')
    corner_cases = [(0, 0, 0, 0), (1, 0, 0, 0), (0, 1, 0, 0),
                    (0, 0, 1, 0), (0, 0, 0, 1), (1, 0, 0, I),
                    (I, 0, 0, -I), (0, I, 0, 0), (0, I, I, 0)]
    s1 = [[Eq(f(t), C1), Eq(g(t), C2)],
          [Eq(f(t), C1*exp(t)), Eq(g(t), -C2)],
          [Eq(f(t), C1 + C2*t), Eq(g(t), C2)],
          [Eq(f(t), C2), Eq(g(t), C1 + C2*t)],
          [Eq(f(t), -C2), Eq(g(t), C1*exp(t))],
          [Eq(f(t), C1*(1 - I)*exp(t)), Eq(g(t), C2*(-1 + I)*exp(I*t))],
          [Eq(f(t), 2*I*C1*exp(I*t)), Eq(g(t), -2*I*C2*exp(-I*t))],
          [Eq(f(t), I*C1 + I*C2*t), Eq(g(t), C2)],
          [Eq(f(t), I*C1*exp(I*t) + I*C2*exp(-I*t)), \
           Eq(g(t), I*C1*exp(I*t) - I*C2*exp(-I*t))]
         ]
    for r, sol in zip(corner_cases, s1):
        eq = [Eq(diff(f(t), t), r[0]*f(t) + r[1]*g(t)),
              Eq(diff(g(t), t), r[2]*f(t) + r[3]*g(t))]
        assert checksysodesol(eq, sol) == (True, [0, 0])


def test_dsolve_linsystem_symbol_piecewise():

    # example from https://groups.google.com/d/msg/sympy/xmzoqW6tWaE/sf0bgQrlCgAJ
    i, r1, c1, r2, c2, t = symbols('i, r1, c1, r2, c2, t')
    x1 = Function('x1')
    x2 = Function('x2')
    eq1 = r1*c1*Derivative(x1(t), t) + x1(t) - x2(t) - r1*i
    eq2 = r2*c1*Derivative(x1(t), t) + r2*c2*Derivative(x2(t), t) + x2(t) - r2*i
    sol = dsolve((eq1, eq2))
    # FIXME: assert checksysodesol(eq, sol) == (True, [0, 0])
    # Remove line below when checksysodesol works
    assert all(s.has(Piecewise) for s in sol)


@slow
def test_linear_2eq_order2():
    x, y, z = symbols('x, y, z', cls=Function)
    k, l, m, n = symbols('k, l, m, n', Integer=True)
    t, l = symbols('t, l')
    x0, y0 = symbols('x0, y0', cls=Function)

    eq1 = (Eq(diff(x(t),t,t), 5*x(t) + 43*y(t)), Eq(diff(y(t),t,t), x(t) + 9*y(t)))
    sol1 = [Eq(x(t), 43*C1*exp(t*rootof(l**4 - 14*l**2 + 2, 0)) + 43*C2*exp(t*rootof(l**4 - 14*l**2 + 2, 1)) + \
    43*C3*exp(t*rootof(l**4 - 14*l**2 + 2, 2)) + 43*C4*exp(t*rootof(l**4 - 14*l**2 + 2, 3))), \
    Eq(y(t), C1*(rootof(l**4 - 14*l**2 + 2, 0)**2 - 5)*exp(t*rootof(l**4 - 14*l**2 + 2, 0)) + \
    C2*(rootof(l**4 - 14*l**2 + 2, 1)**2 - 5)*exp(t*rootof(l**4 - 14*l**2 + 2, 1)) + \
    C3*(rootof(l**4 - 14*l**2 + 2, 2)**2 - 5)*exp(t*rootof(l**4 - 14*l**2 + 2, 2)) + \
    C4*(rootof(l**4 - 14*l**2 + 2, 3)**2 - 5)*exp(t*rootof(l**4 - 14*l**2 + 2, 3)))]
    assert dsolve(eq1) == sol1
    # FIXME: assert checksysodesol(eq1, sol1) == (True, [0, 0])  # this one fails

    eq2 = (Eq(diff(x(t),t,t), 8*x(t)+3*y(t)+31), Eq(diff(y(t),t,t), 9*x(t)+7*y(t)+12))
    sol2 = [Eq(x(t), 3*C1*exp(t*rootof(l**4 - 15*l**2 + 29, 0)) + 3*C2*exp(t*rootof(l**4 - 15*l**2 + 29, 1)) + \
    3*C3*exp(t*rootof(l**4 - 15*l**2 + 29, 2)) + 3*C4*exp(t*rootof(l**4 - 15*l**2 + 29, 3)) - Rational(181, 29)), \
    Eq(y(t), C1*(rootof(l**4 - 15*l**2 + 29, 0)**2 - 8)*exp(t*rootof(l**4 - 15*l**2 + 29, 0)) + \
    C2*(rootof(l**4 - 15*l**2 + 29, 1)**2 - 8)*exp(t*rootof(l**4 - 15*l**2 + 29, 1)) + \
    C3*(rootof(l**4 - 15*l**2 + 29, 2)**2 - 8)*exp(t*rootof(l**4 - 15*l**2 + 29, 2)) + \
    C4*(rootof(l**4 - 15*l**2 + 29, 3)**2 - 8)*exp(t*rootof(l**4 - 15*l**2 + 29, 3)) + Rational(183, 29))]
    assert dsolve(eq2) == sol2
    # FIXME: assert checksysodesol(eq2, sol2) == (True, [0, 0])  # this one fails

    eq3 = (Eq(diff(x(t),t,t) - 9*diff(y(t),t) + 7*x(t),0), Eq(diff(y(t),t,t) + 9*diff(x(t),t) + 7*y(t),0))
    sol3 = [Eq(x(t), C1*cos(t*(Rational(9, 2) + sqrt(109)/2)) + C2*sin(t*(Rational(9, 2) + sqrt(109)/2)) + C3*cos(t*(-sqrt(109)/2 + Rational(9, 2))) + \
    C4*sin(t*(-sqrt(109)/2 + Rational(9, 2)))), Eq(y(t), -C1*sin(t*(Rational(9, 2) + sqrt(109)/2)) + C2*cos(t*(Rational(9, 2) + sqrt(109)/2)) - \
    C3*sin(t*(-sqrt(109)/2 + Rational(9, 2))) + C4*cos(t*(-sqrt(109)/2 + Rational(9, 2))))]
    assert dsolve(eq3) == sol3
    assert checksysodesol(eq3, sol3) == (True, [0, 0])

    eq4 = (Eq(diff(x(t),t,t), 9*t*diff(y(t),t)-9*y(t)), Eq(diff(y(t),t,t),7*t*diff(x(t),t)-7*x(t)))
    sol4 = [Eq(x(t), C3*t + t*Integral((9*C1*exp(3*sqrt(7)*t**2/2) + 9*C2*exp(-3*sqrt(7)*t**2/2))/t**2, t)), \
    Eq(y(t), C4*t + t*Integral((3*sqrt(7)*C1*exp(3*sqrt(7)*t**2/2) - 3*sqrt(7)*C2*exp(-3*sqrt(7)*t**2/2))/t**2, t))]
    assert dsolve(eq4) == sol4
    assert checksysodesol(eq4, sol4) == (True, [0, 0])

    eq5 = (Eq(diff(x(t),t,t), (log(t)+t**2)*diff(x(t),t)+(log(t)+t**2)*3*diff(y(t),t)), Eq(diff(y(t),t,t), \
    (log(t)+t**2)*2*diff(x(t),t)+(log(t)+t**2)*9*diff(y(t),t)))
    sol5 = [Eq(x(t), -sqrt(22)*(C1*Integral(exp((-sqrt(22) + 5)*Integral(t**2 + log(t), t)), t) + C2 - \
    C3*Integral(exp((sqrt(22) + 5)*Integral(t**2 + log(t), t)), t) - C4 - \
    (sqrt(22) + 5)*(C1*Integral(exp((-sqrt(22) + 5)*Integral(t**2 + log(t), t)), t) + C2) + \
    (-sqrt(22) + 5)*(C3*Integral(exp((sqrt(22) + 5)*Integral(t**2 + log(t), t)), t) + C4))/88), \
    Eq(y(t), -sqrt(22)*(C1*Integral(exp((-sqrt(22) + 5)*Integral(t**2 + log(t), t)), t) + \
    C2 - C3*Integral(exp((sqrt(22) + 5)*Integral(t**2 + log(t), t)), t) - C4)/44)]
    assert dsolve(eq5) == sol5
    assert checksysodesol(eq5, sol5) == (True, [0, 0])

    eq6 = (Eq(diff(x(t),t,t), log(t)*t*diff(y(t),t) - log(t)*y(t)), Eq(diff(y(t),t,t), log(t)*t*diff(x(t),t) - log(t)*x(t)))
    sol6 = [Eq(x(t), C3*t + t*Integral((C1*exp(Integral(t*log(t), t)) + \
    C2*exp(-Integral(t*log(t), t)))/t**2, t)), Eq(y(t), C4*t + t*Integral((C1*exp(Integral(t*log(t), t)) - \
    C2*exp(-Integral(t*log(t), t)))/t**2, t))]
    assert dsolve(eq6) == sol6
    assert checksysodesol(eq6, sol6) == (True, [0, 0])

    eq7 = (Eq(diff(x(t),t,t), log(t)*(t*diff(x(t),t) - x(t)) + exp(t)*(t*diff(y(t),t) - y(t))), \
    Eq(diff(y(t),t,t), (t**2)*(t*diff(x(t),t) - x(t)) + (t)*(t*diff(y(t),t) - y(t))))
    sol7 = [Eq(x(t), C3*t + t*Integral((C1*x0(t) + C2*x0(t)*Integral(t*exp(t)*exp(Integral(t**2, t))*\
    exp(Integral(t*log(t), t))/x0(t)**2, t))/t**2, t)), Eq(y(t), C4*t + t*Integral((C1*y0(t) + \
    C2*(y0(t)*Integral(t*exp(t)*exp(Integral(t**2, t))*exp(Integral(t*log(t), t))/x0(t)**2, t) + \
    exp(Integral(t**2, t))*exp(Integral(t*log(t), t))/x0(t)))/t**2, t))]
    assert dsolve(eq7) == sol7
    # FIXME: assert checksysodesol(eq7, sol7) == (True, [0, 0])

    eq8 = (Eq(diff(x(t),t,t), t*(4*x(t) + 9*y(t))), Eq(diff(y(t),t,t), t*(12*x(t) - 6*y(t))))
    sol8 = [Eq(x(t), -sqrt(133)*(-4*C1*airyai(t*(-1 + sqrt(133))**(S(1)/3)) + 4*C1*airyai(-t*(1 + \
    sqrt(133))**(S(1)/3)) - 4*C2*airybi(t*(-1 + sqrt(133))**(S(1)/3)) + 4*C2*airybi(-t*(1 + sqrt(133))**(S(1)/3)) +\
    (-sqrt(133) - 1)*(C1*airyai(t*(-1 + sqrt(133))**(S(1)/3)) + C2*airybi(t*(-1 + sqrt(133))**(S(1)/3))) - (-1 +\
    sqrt(133))*(C1*airyai(-t*(1 + sqrt(133))**(S(1)/3)) + C2*airybi(-t*(1 + sqrt(133))**(S(1)/3))))/3192), \
    Eq(y(t), -sqrt(133)*(-C1*airyai(t*(-1 + sqrt(133))**(S(1)/3)) + C1*airyai(-t*(1 + sqrt(133))**(S(1)/3)) -\
    C2*airybi(t*(-1 + sqrt(133))**(S(1)/3)) + C2*airybi(-t*(1 + sqrt(133))**(S(1)/3)))/266)]
    assert dsolve(eq8) == sol8
    assert checksysodesol(eq8, sol8) == (True, [0, 0])
    assert filldedent(dsolve(eq8)) == filldedent('''
        [Eq(x(t), -sqrt(133)*(-4*C1*airyai(t*(-1 + sqrt(133))**(1/3)) +
        4*C1*airyai(-t*(1 + sqrt(133))**(1/3)) - 4*C2*airybi(t*(-1 +
        sqrt(133))**(1/3)) + 4*C2*airybi(-t*(1 + sqrt(133))**(1/3)) +
        (-sqrt(133) - 1)*(C1*airyai(t*(-1 + sqrt(133))**(1/3)) +
        C2*airybi(t*(-1 + sqrt(133))**(1/3))) - (-1 +
        sqrt(133))*(C1*airyai(-t*(1 + sqrt(133))**(1/3)) + C2*airybi(-t*(1 +
        sqrt(133))**(1/3))))/3192), Eq(y(t), -sqrt(133)*(-C1*airyai(t*(-1 +
        sqrt(133))**(1/3)) + C1*airyai(-t*(1 + sqrt(133))**(1/3)) -
        C2*airybi(t*(-1 + sqrt(133))**(1/3)) + C2*airybi(-t*(1 +
        sqrt(133))**(1/3)))/266)]''')
    assert checksysodesol(eq8, sol8) == (True, [0, 0])

    eq9 = (Eq(diff(x(t),t,t), t*(4*diff(x(t),t) + 9*diff(y(t),t))), Eq(diff(y(t),t,t), t*(12*diff(x(t),t) - 6*diff(y(t),t))))
    sol9 = [Eq(x(t), -sqrt(133)*(4*C1*Integral(exp((-sqrt(133) - 1)*Integral(t, t)), t) + 4*C2 - \
    4*C3*Integral(exp((-1 + sqrt(133))*Integral(t, t)), t) - 4*C4 - (-1 + sqrt(133))*(C1*Integral(exp((-sqrt(133) - \
    1)*Integral(t, t)), t) + C2) + (-sqrt(133) - 1)*(C3*Integral(exp((-1 + sqrt(133))*Integral(t, t)), t) + \
    C4))/3192), Eq(y(t), -sqrt(133)*(C1*Integral(exp((-sqrt(133) - 1)*Integral(t, t)), t) + C2 - \
    C3*Integral(exp((-1 + sqrt(133))*Integral(t, t)), t) - C4)/266)]
    assert dsolve(eq9) == sol9
    assert checksysodesol(eq9, sol9) == (True, [0, 0])

    eq10 = (t**2*diff(x(t),t,t) + 3*t*diff(x(t),t) + 4*t*diff(y(t),t) + 12*x(t) + 9*y(t), \
    t**2*diff(y(t),t,t) + 2*t*diff(x(t),t) - 5*t*diff(y(t),t) + 15*x(t) + 8*y(t))
    sol10 = [Eq(x(t), -C1*(-2*sqrt(-346/(3*(Rational(4333, 4) + 5*sqrt(70771857)/36)**Rational(1, 3)) + 4 + 2*(Rational(4333, 4) + \
    5*sqrt(70771857)/36)**Rational(1, 3)) + 13 + 2*sqrt(-284/sqrt(-346/(3*(Rational(4333, 4) + 5*sqrt(70771857)/36)**Rational(1, 3)) + \
    4 + 2*(Rational(4333, 4) + 5*sqrt(70771857)/36)**Rational(1, 3)) - 2*(Rational(4333, 4) + 5*sqrt(70771857)/36)**Rational(1, 3) + 8 + \
    346/(3*(Rational(4333, 4) + 5*sqrt(70771857)/36)**Rational(1, 3))))*exp((-sqrt(-346/(3*(Rational(4333, 4) + 5*sqrt(70771857)/36)**Rational(1, 3)) + \
    4 + 2*(Rational(4333, 4) + 5*sqrt(70771857)/36)**Rational(1, 3))/2 + 1 + sqrt(-284/sqrt(-346/(3*(Rational(4333, 4) + \
    5*sqrt(70771857)/36)**Rational(1, 3)) + 4 + 2*(Rational(4333, 4) + 5*sqrt(70771857)/36)**Rational(1, 3)) - 2*(Rational(4333, 4) + \
    5*sqrt(70771857)/36)**Rational(1, 3) + 8 + 346/(3*(Rational(4333, 4) + 5*sqrt(70771857)/36)**Rational(1, 3)))/2)*log(t)) - \
    C2*(-2*sqrt(-346/(3*(Rational(4333, 4) + 5*sqrt(70771857)/36)**Rational(1, 3)) + 4 + 2*(Rational(4333, 4) + 5*sqrt(70771857)/36)**Rational(1, 3)) + \
    13 - 2*sqrt(-284/sqrt(-346/(3*(Rational(4333, 4) + 5*sqrt(70771857)/36)**Rational(1, 3)) + 4 + 2*(Rational(4333, 4) + \
    5*sqrt(70771857)/36)**Rational(1, 3)) - 2*(Rational(4333, 4) + 5*sqrt(70771857)/36)**Rational(1, 3) + 8 + 346/(3*(Rational(4333, 4) + \
    5*sqrt(70771857)/36)**Rational(1, 3))))*exp((-sqrt(-346/(3*(Rational(4333, 4) + 5*sqrt(70771857)/36)**Rational(1, 3)) + 4 + \
    2*(Rational(4333, 4) + 5*sqrt(70771857)/36)**Rational(1, 3))/2 + 1 - sqrt(-284/sqrt(-346/(3*(Rational(4333, 4) + 5*sqrt(70771857)/36)**Rational(1, 3)) + \
    4 + 2*(Rational(4333, 4) + 5*sqrt(70771857)/36)**Rational(1, 3)) - 2*(Rational(4333, 4) + 5*sqrt(70771857)/36)**Rational(1, 3) + 8 + 346/(3*(Rational(4333, 4) + \
    5*sqrt(70771857)/36)**Rational(1, 3)))/2)*log(t)) - C3*t**(1 + sqrt(-346/(3*(Rational(4333, 4) + 5*sqrt(70771857)/36)**Rational(1, 3)) + 4 + \
    2*(Rational(4333, 4) + 5*sqrt(70771857)/36)**Rational(1, 3))/2 + sqrt(-2*(Rational(4333, 4) + 5*sqrt(70771857)/36)**Rational(1, 3) + 8 + 346/(3*(Rational(4333, 4) + \
    5*sqrt(70771857)/36)**Rational(1, 3)) + 284/sqrt(-346/(3*(Rational(4333, 4) + 5*sqrt(70771857)/36)**Rational(1, 3)) + 4 + 2*(Rational(4333, 4) + \
    5*sqrt(70771857)/36)**Rational(1, 3)))/2)*(2*sqrt(-346/(3*(Rational(4333, 4) + 5*sqrt(70771857)/36)**Rational(1, 3)) + 4 + 2*(Rational(4333, 4) + \
    5*sqrt(70771857)/36)**Rational(1, 3)) + 13 + 2*sqrt(-2*(Rational(4333, 4) + 5*sqrt(70771857)/36)**Rational(1, 3) + 8 + 346/(3*(Rational(4333, 4) + \
    5*sqrt(70771857)/36)**Rational(1, 3)) + 284/sqrt(-346/(3*(Rational(4333, 4) + 5*sqrt(70771857)/36)**Rational(1, 3)) + 4 + 2*(Rational(4333, 4) + \
    5*sqrt(70771857)/36)**Rational(1, 3)))) - C4*t**(-sqrt(-2*(Rational(4333, 4) + 5*sqrt(70771857)/36)**Rational(1, 3) + 8 + 346/(3*(Rational(4333, 4) + \
    5*sqrt(70771857)/36)**Rational(1, 3)) + 284/sqrt(-346/(3*(Rational(4333, 4) + 5*sqrt(70771857)/36)**Rational(1, 3)) + 4 + 2*(Rational(4333, 4) + \
    5*sqrt(70771857)/36)**Rational(1, 3)))/2 + 1 + sqrt(-346/(3*(Rational(4333, 4) + 5*sqrt(70771857)/36)**Rational(1, 3)) + 4 + 2*(Rational(4333, 4) + \
    5*sqrt(70771857)/36)**Rational(1, 3))/2)*(-2*sqrt(-2*(Rational(4333, 4) + 5*sqrt(70771857)/36)**Rational(1, 3) + 8 + 346/(3*(Rational(4333, 4) + \
    5*sqrt(70771857)/36)**Rational(1, 3)) + 284/sqrt(-346/(3*(Rational(4333, 4) + 5*sqrt(70771857)/36)**Rational(1, 3)) + 4 + 2*(Rational(4333, 4) + \
    5*sqrt(70771857)/36)**Rational(1, 3))) + 2*sqrt(-346/(3*(Rational(4333, 4) + 5*sqrt(70771857)/36)**Rational(1, 3)) + 4 + 2*(Rational(4333, 4) + \
    5*sqrt(70771857)/36)**Rational(1, 3)) + 13)), Eq(y(t), C1*(-sqrt(-346/(3*(Rational(4333, 4) + 5*sqrt(70771857)/36)**Rational(1, 3)) + 4 + \
    2*(Rational(4333, 4) + 5*sqrt(70771857)/36)**Rational(1, 3)) + 14 + (-sqrt(-346/(3*(Rational(4333, 4) + 5*sqrt(70771857)/36)**Rational(1, 3)) + 4 + \
    2*(Rational(4333, 4) + 5*sqrt(70771857)/36)**Rational(1, 3))/2 + 1 + sqrt(-284/sqrt(-346/(3*(Rational(4333, 4) + 5*sqrt(70771857)/36)**Rational(1, 3)) + \
    4 + 2*(Rational(4333, 4) + 5*sqrt(70771857)/36)**Rational(1, 3)) - 2*(Rational(4333, 4) + 5*sqrt(70771857)/36)**Rational(1, 3) + 8 + 346/(3*(Rational(4333, 4) + \
    5*sqrt(70771857)/36)**Rational(1, 3)))/2)**2 + sqrt(-284/sqrt(-346/(3*(Rational(4333, 4) + 5*sqrt(70771857)/36)**Rational(1, 3)) + 4 + \
    2*(Rational(4333, 4) + 5*sqrt(70771857)/36)**Rational(1, 3)) - 2*(Rational(4333, 4) + 5*sqrt(70771857)/36)**Rational(1, 3) + 8 + 346/(3*(Rational(4333, 4) + \
    5*sqrt(70771857)/36)**Rational(1, 3))))*exp((-sqrt(-346/(3*(Rational(4333, 4) + 5*sqrt(70771857)/36)**Rational(1, 3)) + 4 + 2*(Rational(4333, 4) + \
    5*sqrt(70771857)/36)**Rational(1, 3))/2 + 1 + sqrt(-284/sqrt(-346/(3*(Rational(4333, 4) + 5*sqrt(70771857)/36)**Rational(1, 3)) + 4 + \
    2*(Rational(4333, 4) + 5*sqrt(70771857)/36)**Rational(1, 3)) - 2*(Rational(4333, 4) + 5*sqrt(70771857)/36)**Rational(1, 3) + 8 + 346/(3*(Rational(4333, 4) + \
    5*sqrt(70771857)/36)**Rational(1, 3)))/2)*log(t)) + C2*(-sqrt(-346/(3*(Rational(4333, 4) + 5*sqrt(70771857)/36)**Rational(1, 3)) + 4 + \
    2*(Rational(4333, 4) + 5*sqrt(70771857)/36)**Rational(1, 3)) + 14 - sqrt(-284/sqrt(-346/(3*(Rational(4333, 4) + 5*sqrt(70771857)/36)**Rational(1, 3)) + \
    4 + 2*(Rational(4333, 4) + 5*sqrt(70771857)/36)**Rational(1, 3)) - 2*(Rational(4333, 4) + 5*sqrt(70771857)/36)**Rational(1, 3) + 8 + 346/(3*(Rational(4333, 4) + \
    5*sqrt(70771857)/36)**Rational(1, 3))) + (-sqrt(-346/(3*(Rational(4333, 4) + 5*sqrt(70771857)/36)**Rational(1, 3)) + 4 + 2*(Rational(4333, 4) + \
    5*sqrt(70771857)/36)**Rational(1, 3))/2 + 1 - sqrt(-284/sqrt(-346/(3*(Rational(4333, 4) + 5*sqrt(70771857)/36)**Rational(1, 3)) + 4 + \
    2*(Rational(4333, 4) + 5*sqrt(70771857)/36)**Rational(1, 3)) - 2*(Rational(4333, 4) + 5*sqrt(70771857)/36)**Rational(1, 3) + 8 + 346/(3*(Rational(4333, 4) + \
    5*sqrt(70771857)/36)**Rational(1, 3)))/2)**2)*exp((-sqrt(-346/(3*(Rational(4333, 4) + 5*sqrt(70771857)/36)**Rational(1, 3)) + 4 + 2*(Rational(4333, 4) + \
    5*sqrt(70771857)/36)**Rational(1, 3))/2 + 1 - sqrt(-284/sqrt(-346/(3*(Rational(4333, 4) + 5*sqrt(70771857)/36)**Rational(1, 3)) + 4 + \
    2*(Rational(4333, 4) + 5*sqrt(70771857)/36)**Rational(1, 3)) - 2*(Rational(4333, 4) + 5*sqrt(70771857)/36)**Rational(1, 3) + 8 + 346/(3*(Rational(4333, 4) + \
    5*sqrt(70771857)/36)**Rational(1, 3)))/2)*log(t)) + C3*t**(1 + sqrt(-346/(3*(Rational(4333, 4) + 5*sqrt(70771857)/36)**Rational(1, 3)) + 4 + \
    2*(Rational(4333, 4) + 5*sqrt(70771857)/36)**Rational(1, 3))/2 + sqrt(-2*(Rational(4333, 4) + 5*sqrt(70771857)/36)**Rational(1, 3) + 8 + 346/(3*(Rational(4333, 4) + \
    5*sqrt(70771857)/36)**Rational(1, 3)) + 284/sqrt(-346/(3*(Rational(4333, 4) + 5*sqrt(70771857)/36)**Rational(1, 3)) + 4 + 2*(Rational(4333, 4) + \
    5*sqrt(70771857)/36)**Rational(1, 3)))/2)*(sqrt(-346/(3*(Rational(4333, 4) + 5*sqrt(70771857)/36)**Rational(1, 3)) + 4 + 2*(Rational(4333, 4) + \
    5*sqrt(70771857)/36)**Rational(1, 3)) + sqrt(-2*(Rational(4333, 4) + 5*sqrt(70771857)/36)**Rational(1, 3) + 8 + 346/(3*(Rational(4333, 4) + \
    5*sqrt(70771857)/36)**Rational(1, 3)) + 284/sqrt(-346/(3*(Rational(4333, 4) + 5*sqrt(70771857)/36)**Rational(1, 3)) + 4 + 2*(Rational(4333, 4) + \
    5*sqrt(70771857)/36)**Rational(1, 3))) + 14 + (1 + sqrt(-346/(3*(Rational(4333, 4) + 5*sqrt(70771857)/36)**Rational(1, 3)) + 4 + 2*(Rational(4333, 4) + \
    5*sqrt(70771857)/36)**Rational(1, 3))/2 + sqrt(-2*(Rational(4333, 4) + 5*sqrt(70771857)/36)**Rational(1, 3) + 8 + 346/(3*(Rational(4333, 4) + \
    5*sqrt(70771857)/36)**Rational(1, 3)) + 284/sqrt(-346/(3*(Rational(4333, 4) + 5*sqrt(70771857)/36)**Rational(1, 3)) + 4 + 2*(Rational(4333, 4) + \
    5*sqrt(70771857)/36)**Rational(1, 3)))/2)**2) + C4*t**(-sqrt(-2*(Rational(4333, 4) + 5*sqrt(70771857)/36)**Rational(1, 3) + 8 + \
    346/(3*(Rational(4333, 4) + 5*sqrt(70771857)/36)**Rational(1, 3)) + 284/sqrt(-346/(3*(Rational(4333, 4) + 5*sqrt(70771857)/36)**Rational(1, 3)) + \
    4 + 2*(Rational(4333, 4) + 5*sqrt(70771857)/36)**Rational(1, 3)))/2 + 1 + sqrt(-346/(3*(Rational(4333, 4) + 5*sqrt(70771857)/36)**Rational(1, 3)) + \
    4 + 2*(Rational(4333, 4) + 5*sqrt(70771857)/36)**Rational(1, 3))/2)*(-sqrt(-2*(Rational(4333, 4) + 5*sqrt(70771857)/36)**Rational(1, 3) + \
    8 + 346/(3*(Rational(4333, 4) + 5*sqrt(70771857)/36)**Rational(1, 3)) + 284/sqrt(-346/(3*(Rational(4333, 4) + 5*sqrt(70771857)/36)**Rational(1, 3)) + \
    4 + 2*(Rational(4333, 4) + 5*sqrt(70771857)/36)**Rational(1, 3))) + (-sqrt(-2*(Rational(4333, 4) + 5*sqrt(70771857)/36)**Rational(1, 3) + 8 + \
    346/(3*(Rational(4333, 4) + 5*sqrt(70771857)/36)**Rational(1, 3)) + 284/sqrt(-346/(3*(Rational(4333, 4) + 5*sqrt(70771857)/36)**Rational(1, 3)) + \
    4 + 2*(Rational(4333, 4) + 5*sqrt(70771857)/36)**Rational(1, 3)))/2 + 1 + sqrt(-346/(3*(Rational(4333, 4) + 5*sqrt(70771857)/36)**Rational(1, 3)) + \
    4 + 2*(Rational(4333, 4) + 5*sqrt(70771857)/36)**Rational(1, 3))/2)**2 + sqrt(-346/(3*(Rational(4333, 4) + \
    5*sqrt(70771857)/36)**Rational(1, 3)) + 4 + 2*(Rational(4333, 4) + 5*sqrt(70771857)/36)**Rational(1, 3)) + 14))]
    assert dsolve(eq10) == sol10
    # FIXME: assert checksysodesol(eq10, sol10) == (True, [0, 0]) # this hangs or at least takes a while...


def test_linear_3eq_order1():
    x, y, z = symbols('x, y, z', cls=Function)
    t = Symbol('t')

    f = t**3 + log(t)
    g = t**2 + sin(t)
    eq1 = (Eq(diff(x(t),t),(4*f+g)*x(t)-f*y(t)-2*f*z(t)), Eq(diff(y(t),t),2*f*x(t)+(f+g)*y(t)-2*f*z(t)), Eq(diff(z(t),t),5*f*x(t)+f*y(t)+(-3*f+g)*z(t)))
    # sol4 = [Eq(x(t), (C1*exp(-2*Integral(t**3 + log(t), t)) + C2*(sqrt(3)*sin(sqrt(3)*Integral(t**3 + log(t), t))/6 \
    # + cos(sqrt(3)*Integral(t**3 + log(t), t))/2) + C3*(sin(sqrt(3)*Integral(t**3 + log(t), t))/2 - \
    # sqrt(3)*cos(sqrt(3)*Integral(t**3 + log(t), t))/6))*exp(Integral(-t**2 - sin(t), t))), Eq(y(t), \
    # (C2*(sqrt(3)*sin(sqrt(3)*Integral(t**3 + log(t), t))/6 + cos(sqrt(3)*Integral(t**3 + log(t), t))/2) + \
    # C3*(sin(sqrt(3)*Integral(t**3 + log(t), t))/2 - sqrt(3)*cos(sqrt(3)*Integral(t**3 + log(t), t))/6))*\
    # exp(Integral(-t**2 - sin(t), t))), Eq(z(t), (C1*exp(-2*Integral(t**3 + log(t), t)) + C2*cos(sqrt(3)*\
    # Integral(t**3 + log(t), t)) + C3*sin(sqrt(3)*Integral(t**3 + log(t), t)))*exp(Integral(-t**2 - sin(t), t)))]
    sol1 = [Eq(x(t), (C1*exp(-2*Integral(t**3 + log(t), t)) + C2*(sqrt(3)*sin(sqrt(3)*Integral(t**3 + log(t), t))/6 \
            + cos(sqrt(3)*Integral(t**3 + log(t), t))/2) + C3*(-sin(sqrt(3)*Integral(t**3 + log(t), t))/2 \
            + sqrt(3)*cos(sqrt(3)*Integral(t**3 + log(t), t))/6))*exp(Integral(-t**2 - sin(t), t))),
            Eq(y(t), (C2*(sqrt(3)*sin(sqrt(3)*Integral(t**3 + log(t), t))/6 + cos(sqrt(3)* \
            Integral(t**3 + log(t), t))/2) + C3*(-sin(sqrt(3)*Integral(t**3 + log(t), t))/2 \
            + sqrt(3)*cos(sqrt(3)*Integral(t**3 + log(t), t))/6))*exp(Integral(-t**2 - sin(t), t))),
            Eq(z(t), (C1*exp(-2*Integral(t**3 + log(t), t)) + C2*cos(sqrt(3)*Integral(t**3 + log(t), t)) - \
            C3*sin(sqrt(3)*Integral(t**3 + log(t), t)))*exp(Integral(-t**2 - sin(t), t)))]
    dsolve_sol = dsolve(eq1)
    # dsolve_sol = [eq.subs(C3, -C3) for eq in dsolve_sol]
    assert all(simplify(s1.rhs - ds1.rhs) == 0 for s1, ds1 in zip(sol1, dsolve_sol))


def test_linear_3eq_order1_nonhomog():
    e = [Eq(diff(f(x), x), -9*f(x) - 4*g(x)),
         Eq(diff(g(x), x), -4*g(x)),
         Eq(diff(h(x), x), h(x) + exp(x))]
    raises(NotImplementedError, lambda: dsolve(e))


def test_nonlinear_2eq_order1():
    x, y, z = symbols('x, y, z', cls=Function)
    t = Symbol('t')
    eq1 = (Eq(diff(x(t),t),x(t)*y(t)**3), Eq(diff(y(t),t),y(t)**5))
    sol1 = [
        Eq(x(t), C1*exp((-1/(4*C2 + 4*t))**(Rational(-1, 4)))),
        Eq(y(t), -(-1/(4*C2 + 4*t))**Rational(1, 4)),
        Eq(x(t), C1*exp(-1/(-1/(4*C2 + 4*t))**Rational(1, 4))),
        Eq(y(t), (-1/(4*C2 + 4*t))**Rational(1, 4)),
        Eq(x(t), C1*exp(-I/(-1/(4*C2 + 4*t))**Rational(1, 4))),
        Eq(y(t), -I*(-1/(4*C2 + 4*t))**Rational(1, 4)),
        Eq(x(t), C1*exp(I/(-1/(4*C2 + 4*t))**Rational(1, 4))),
        Eq(y(t), I*(-1/(4*C2 + 4*t))**Rational(1, 4))]
    assert dsolve(eq1) == sol1
    assert checksysodesol(eq1, sol1) == (True, [0, 0])

    eq2 = (Eq(diff(x(t),t), exp(3*x(t))*y(t)**3),Eq(diff(y(t),t), y(t)**5))
    sol2 = [
        Eq(x(t), -log(C1 - 3/(-1/(4*C2 + 4*t))**Rational(1, 4))/3),
        Eq(y(t), -(-1/(4*C2 + 4*t))**Rational(1, 4)),
        Eq(x(t), -log(C1 + 3/(-1/(4*C2 + 4*t))**Rational(1, 4))/3),
        Eq(y(t), (-1/(4*C2 + 4*t))**Rational(1, 4)),
        Eq(x(t), -log(C1 + 3*I/(-1/(4*C2 + 4*t))**Rational(1, 4))/3),
        Eq(y(t), -I*(-1/(4*C2 + 4*t))**Rational(1, 4)),
        Eq(x(t), -log(C1 - 3*I/(-1/(4*C2 + 4*t))**Rational(1, 4))/3),
        Eq(y(t), I*(-1/(4*C2 + 4*t))**Rational(1, 4))]
    assert dsolve(eq2) == sol2
    assert checksysodesol(eq2, sol2) == (True, [0, 0])

    eq3 = (Eq(diff(x(t),t), y(t)*x(t)), Eq(diff(y(t),t), x(t)**3))
    tt = Rational(2, 3)
    sol3 = [
        Eq(x(t), 6**tt/(6*(-sinh(sqrt(C1)*(C2 + t)/2)/sqrt(C1))**tt)),
        Eq(y(t), sqrt(C1 + C1/sinh(sqrt(C1)*(C2 + t)/2)**2)/3)]
    assert dsolve(eq3) == sol3
    # FIXME: assert checksysodesol(eq3, sol3) == (True, [0, 0])

    eq4 = (Eq(diff(x(t),t),x(t)*y(t)*sin(t)**2), Eq(diff(y(t),t),y(t)**2*sin(t)**2))
    sol4 = set([Eq(x(t), -2*exp(C1)/(C2*exp(C1) + t - sin(2*t)/2)), Eq(y(t), -2/(C1 + t - sin(2*t)/2))])
    assert dsolve(eq4) == sol4
    # FIXME: assert checksysodesol(eq4, sol4) == (True, [0, 0])

    eq5 = (Eq(x(t),t*diff(x(t),t)+diff(x(t),t)*diff(y(t),t)), Eq(y(t),t*diff(y(t),t)+diff(y(t),t)**2))
    sol5 = set([Eq(x(t), C1*C2 + C1*t), Eq(y(t), C2**2 + C2*t)])
    assert dsolve(eq5) == sol5
    assert checksysodesol(eq5, sol5) == (True, [0, 0])

    eq6 = (Eq(diff(x(t),t),x(t)**2*y(t)**3), Eq(diff(y(t),t),y(t)**5))
    sol6 = [
        Eq(x(t), 1/(C1 - 1/(-1/(4*C2 + 4*t))**Rational(1, 4))),
        Eq(y(t), -(-1/(4*C2 + 4*t))**Rational(1, 4)),
        Eq(x(t), 1/(C1 + (-1/(4*C2 + 4*t))**(Rational(-1, 4)))),
        Eq(y(t), (-1/(4*C2 + 4*t))**Rational(1, 4)),
        Eq(x(t), 1/(C1 + I/(-1/(4*C2 + 4*t))**Rational(1, 4))),
        Eq(y(t), -I*(-1/(4*C2 + 4*t))**Rational(1, 4)),
        Eq(x(t), 1/(C1 - I/(-1/(4*C2 + 4*t))**Rational(1, 4))),
        Eq(y(t), I*(-1/(4*C2 + 4*t))**Rational(1, 4))]
    assert dsolve(eq6) == sol6
    assert checksysodesol(eq6, sol6) == (True, [0, 0])


@slow
def test_nonlinear_3eq_order1():
    x, y, z = symbols('x, y, z', cls=Function)
    t, u = symbols('t u')
    eq1 = (4*diff(x(t),t) + 2*y(t)*z(t), 3*diff(y(t),t) - z(t)*x(t), 5*diff(z(t),t) - x(t)*y(t))
    sol1 = [Eq(4*Integral(1/(sqrt(-4*u**2 - 3*C1 + C2)*sqrt(-4*u**2 + 5*C1 - C2)), (u, x(t))),
        C3 - sqrt(15)*t/15), Eq(3*Integral(1/(sqrt(-6*u**2 - C1 + 5*C2)*sqrt(3*u**2 + C1 - 4*C2)),
        (u, y(t))), C3 + sqrt(5)*t/10), Eq(5*Integral(1/(sqrt(-10*u**2 - 3*C1 + C2)*
        sqrt(5*u**2 + 4*C1 - C2)), (u, z(t))), C3 + sqrt(3)*t/6)]
    assert [i.dummy_eq(j) for i, j in zip(dsolve(eq1), sol1)]
    # FIXME: assert checksysodesol(eq1, sol1) == (True, [0, 0, 0])

    eq2 = (4*diff(x(t),t) + 2*y(t)*z(t)*sin(t), 3*diff(y(t),t) - z(t)*x(t)*sin(t), 5*diff(z(t),t) - x(t)*y(t)*sin(t))
    sol2 = [Eq(3*Integral(1/(sqrt(-6*u**2 - C1 + 5*C2)*sqrt(3*u**2 + C1 - 4*C2)), (u, x(t))), C3 +
        sqrt(5)*cos(t)/10), Eq(4*Integral(1/(sqrt(-4*u**2 - 3*C1 + C2)*sqrt(-4*u**2 + 5*C1 - C2)),
        (u, y(t))), C3 - sqrt(15)*cos(t)/15), Eq(5*Integral(1/(sqrt(-10*u**2 - 3*C1 + C2)*
        sqrt(5*u**2 + 4*C1 - C2)), (u, z(t))), C3 + sqrt(3)*cos(t)/6)]
    assert [i.dummy_eq(j) for i, j in zip(dsolve(eq2), sol2)]
    # FIXME: assert checksysodesol(eq2, sol2) == (True, [0, 0, 0])


@slow
def test_dsolve_options():
    eq = x*f(x).diff(x) + f(x)
    a = dsolve(eq, hint='all')
    b = dsolve(eq, hint='all', simplify=False)
    c = dsolve(eq, hint='all_Integral')
    keys = ['1st_exact', '1st_exact_Integral', '1st_homogeneous_coeff_best',
        '1st_homogeneous_coeff_subs_dep_div_indep',
        '1st_homogeneous_coeff_subs_dep_div_indep_Integral',
        '1st_homogeneous_coeff_subs_indep_div_dep',
        '1st_homogeneous_coeff_subs_indep_div_dep_Integral', '1st_linear',
        '1st_linear_Integral', 'Bernoulli', 'Bernoulli_Integral',
        'almost_linear', 'almost_linear_Integral', 'best', 'best_hint',
        'default', 'lie_group',
        'nth_linear_euler_eq_homogeneous', 'order',
        'separable', 'separable_Integral']
    Integral_keys = ['1st_exact_Integral',
        '1st_homogeneous_coeff_subs_dep_div_indep_Integral',
        '1st_homogeneous_coeff_subs_indep_div_dep_Integral', '1st_linear_Integral',
        'Bernoulli_Integral', 'almost_linear_Integral', 'best', 'best_hint', 'default',
        'nth_linear_euler_eq_homogeneous',
        'order', 'separable_Integral']
    assert sorted(a.keys()) == keys
    assert a['order'] == ode_order(eq, f(x))
    assert a['best'] == Eq(f(x), C1/x)
    assert dsolve(eq, hint='best') == Eq(f(x), C1/x)
    assert a['default'] == 'separable'
    assert a['best_hint'] == 'separable'
    assert not a['1st_exact'].has(Integral)
    assert not a['separable'].has(Integral)
    assert not a['1st_homogeneous_coeff_best'].has(Integral)
    assert not a['1st_homogeneous_coeff_subs_dep_div_indep'].has(Integral)
    assert not a['1st_homogeneous_coeff_subs_indep_div_dep'].has(Integral)
    assert not a['1st_linear'].has(Integral)
    assert a['1st_linear_Integral'].has(Integral)
    assert a['1st_exact_Integral'].has(Integral)
    assert a['1st_homogeneous_coeff_subs_dep_div_indep_Integral'].has(Integral)
    assert a['1st_homogeneous_coeff_subs_indep_div_dep_Integral'].has(Integral)
    assert a['separable_Integral'].has(Integral)
    assert sorted(b.keys()) == keys
    assert b['order'] == ode_order(eq, f(x))
    assert b['best'] == Eq(f(x), C1/x)
    assert dsolve(eq, hint='best', simplify=False) == Eq(f(x), C1/x)
    assert b['default'] == 'separable'
    assert b['best_hint'] == '1st_linear'
    assert a['separable'] != b['separable']
    assert a['1st_homogeneous_coeff_subs_dep_div_indep'] != \
        b['1st_homogeneous_coeff_subs_dep_div_indep']
    assert a['1st_homogeneous_coeff_subs_indep_div_dep'] != \
        b['1st_homogeneous_coeff_subs_indep_div_dep']
    assert not b['1st_exact'].has(Integral)
    assert not b['separable'].has(Integral)
    assert not b['1st_homogeneous_coeff_best'].has(Integral)
    assert not b['1st_homogeneous_coeff_subs_dep_div_indep'].has(Integral)
    assert not b['1st_homogeneous_coeff_subs_indep_div_dep'].has(Integral)
    assert not b['1st_linear'].has(Integral)
    assert b['1st_linear_Integral'].has(Integral)
    assert b['1st_exact_Integral'].has(Integral)
    assert b['1st_homogeneous_coeff_subs_dep_div_indep_Integral'].has(Integral)
    assert b['1st_homogeneous_coeff_subs_indep_div_dep_Integral'].has(Integral)
    assert b['separable_Integral'].has(Integral)
    assert sorted(c.keys()) == Integral_keys
    raises(ValueError, lambda: dsolve(eq, hint='notarealhint'))
    raises(ValueError, lambda: dsolve(eq, hint='Liouville'))
    assert dsolve(f(x).diff(x) - 1/f(x)**2, hint='all')['best'] == \
        dsolve(f(x).diff(x) - 1/f(x)**2, hint='best')
    assert dsolve(f(x) + f(x).diff(x) + sin(x).diff(x) + 1, f(x),
                  hint="1st_linear_Integral") == \
        Eq(f(x), (C1 + Integral((-sin(x).diff(x) - 1)*
                exp(Integral(1, x)), x))*exp(-Integral(1, x)))


def test_classify_ode():
    assert classify_ode(f(x).diff(x, 2), f(x)) == \
        (
        'nth_algebraic',
        'nth_linear_constant_coeff_homogeneous',
        'nth_linear_euler_eq_homogeneous',
        'Liouville',
        '2nd_power_series_ordinary',
        'nth_algebraic_Integral',
        'Liouville_Integral',
        )
    assert classify_ode(f(x), f(x)) == ('nth_algebraic', 'nth_algebraic_Integral')
    assert classify_ode(Eq(f(x).diff(x), 0), f(x)) == (
        'nth_algebraic',
        'separable',
        '1st_linear',
        'Bernoulli',
        '1st_homogeneous_coeff_best',
        '1st_homogeneous_coeff_subs_indep_div_dep',
        '1st_homogeneous_coeff_subs_dep_div_indep',
        '1st_power_series', 'lie_group',
        'nth_linear_constant_coeff_homogeneous',
        'nth_linear_euler_eq_homogeneous',
        'nth_algebraic_Integral',
        'separable_Integral',
        '1st_linear_Integral',
        'Bernoulli_Integral',
        '1st_homogeneous_coeff_subs_indep_div_dep_Integral',
        '1st_homogeneous_coeff_subs_dep_div_indep_Integral')
    assert classify_ode(f(x).diff(x)**2, f(x)) == ('factorable',
         'nth_algebraic',
         'separable',
         '1st_linear',
         'Bernoulli',
         '1st_homogeneous_coeff_best',
         '1st_homogeneous_coeff_subs_indep_div_dep',
         '1st_homogeneous_coeff_subs_dep_div_indep',
         '1st_power_series',
         'lie_group',
         'nth_linear_constant_coeff_homogeneous',
         'nth_linear_euler_eq_homogeneous',
         'nth_algebraic_Integral',
         'separable_Integral',
         '1st_linear_Integral',
         'Bernoulli_Integral',
         '1st_homogeneous_coeff_subs_indep_div_dep_Integral',
         '1st_homogeneous_coeff_subs_dep_div_indep_Integral')
    # issue 4749: f(x) should be cleared from highest derivative before classifying
    a = classify_ode(Eq(f(x).diff(x) + f(x), x), f(x))
    b = classify_ode(f(x).diff(x)*f(x) + f(x)*f(x) - x*f(x), f(x))
    c = classify_ode(f(x).diff(x)/f(x) + f(x)/f(x) - x/f(x), f(x))
    assert a == ('1st_linear',
        'Bernoulli',
        'almost_linear',
        '1st_power_series', "lie_group",
        'nth_linear_constant_coeff_undetermined_coefficients',
        'nth_linear_constant_coeff_variation_of_parameters',
        '1st_linear_Integral',
        'Bernoulli_Integral',
        'almost_linear_Integral',
        'nth_linear_constant_coeff_variation_of_parameters_Integral')
    assert b == ('factorable',
         '1st_linear',
         'Bernoulli',
         '1st_power_series',
         'lie_group',
         'nth_linear_constant_coeff_undetermined_coefficients',
         'nth_linear_constant_coeff_variation_of_parameters',
         '1st_linear_Integral',
         'Bernoulli_Integral',
         'nth_linear_constant_coeff_variation_of_parameters_Integral')
    assert c == ('1st_linear',
         'Bernoulli',
         '1st_power_series',
         'lie_group',
         'nth_linear_constant_coeff_undetermined_coefficients',
         'nth_linear_constant_coeff_variation_of_parameters',
         '1st_linear_Integral',
         'Bernoulli_Integral',
         'nth_linear_constant_coeff_variation_of_parameters_Integral')

    assert classify_ode(
        2*x*f(x)*f(x).diff(x) + (1 + x)*f(x)**2 - exp(x), f(x)
    ) == ('Bernoulli', 'almost_linear', 'lie_group',
        'Bernoulli_Integral', 'almost_linear_Integral')
    assert 'Riccati_special_minus2' in \
        classify_ode(2*f(x).diff(x) + f(x)**2 - f(x)/x + 3*x**(-2), f(x))
    raises(ValueError, lambda: classify_ode(x + f(x, y).diff(x).diff(
        y), f(x, y)))
    # issue 5176
    k = Symbol('k')
    assert classify_ode(f(x).diff(x)/(k*f(x) + k*x*f(x)) + 2*f(x)/(k*f(x) +
        k*x*f(x)) + x*f(x).diff(x)/(k*f(x) + k*x*f(x)) + z, f(x)) == \
        ('separable', '1st_exact', '1st_linear', 'Bernoulli',
        '1st_power_series', 'lie_group', 'separable_Integral', '1st_exact_Integral',
        '1st_linear_Integral', 'Bernoulli_Integral')
    # preprocessing
    ans = ('nth_algebraic', 'separable', '1st_exact', '1st_linear', 'Bernoulli',
        '1st_homogeneous_coeff_best',
        '1st_homogeneous_coeff_subs_indep_div_dep',
        '1st_homogeneous_coeff_subs_dep_div_indep',
        '1st_power_series', 'lie_group',
        'nth_linear_constant_coeff_undetermined_coefficients',
        'nth_linear_euler_eq_nonhomogeneous_undetermined_coefficients',
        'nth_linear_constant_coeff_variation_of_parameters',
        'nth_linear_euler_eq_nonhomogeneous_variation_of_parameters',
        'nth_algebraic_Integral',
        'separable_Integral', '1st_exact_Integral',
        '1st_linear_Integral',
        'Bernoulli_Integral',
        '1st_homogeneous_coeff_subs_indep_div_dep_Integral',
        '1st_homogeneous_coeff_subs_dep_div_indep_Integral',
        'nth_linear_constant_coeff_variation_of_parameters_Integral',
        'nth_linear_euler_eq_nonhomogeneous_variation_of_parameters_Integral')
    #     w/o f(x) given
    assert classify_ode(diff(f(x) + x, x) + diff(f(x), x)) == ans
    #     w/ f(x) and prep=True
    assert classify_ode(diff(f(x) + x, x) + diff(f(x), x), f(x),
                        prep=True) == ans

    assert classify_ode(Eq(2*x**3*f(x).diff(x), 0), f(x)) == \
        ('factorable', 'nth_algebraic', 'separable', '1st_linear',
         'Bernoulli', '1st_power_series',
         'lie_group', 'nth_linear_euler_eq_homogeneous',
         'nth_algebraic_Integral', 'separable_Integral',
         '1st_linear_Integral', 'Bernoulli_Integral')


    assert classify_ode(Eq(2*f(x)**3*f(x).diff(x), 0), f(x)) == \
        ('factorable', 'nth_algebraic', 'separable', '1st_linear', 'Bernoulli',
         '1st_power_series', 'lie_group', 'nth_algebraic_Integral',
         'separable_Integral', '1st_linear_Integral', 'Bernoulli_Integral')
    # test issue 13864
    assert classify_ode(Eq(diff(f(x), x) - f(x)**x, 0), f(x)) == \
        ('1st_power_series', 'lie_group')
    assert isinstance(classify_ode(Eq(f(x), 5), f(x), dict=True), dict)


def test_classify_ode_ics():
    # Dummy
    eq = f(x).diff(x, x) - f(x)

    # Not f(0) or f'(0)
    ics = {x: 1}
    raises(ValueError, lambda: classify_ode(eq, f(x), ics=ics))


    ############################
    # f(0) type (AppliedUndef) #
    ############################


    # Wrong function
    ics = {g(0): 1}
    raises(ValueError, lambda: classify_ode(eq, f(x), ics=ics))

    # Contains x
    ics = {f(x): 1}
    raises(ValueError, lambda: classify_ode(eq, f(x), ics=ics))

    # Too many args
    ics = {f(0, 0): 1}
    raises(ValueError, lambda: classify_ode(eq, f(x), ics=ics))

    # point contains f
    # XXX: Should be NotImplementedError
    ics = {f(0): f(1)}
    raises(ValueError, lambda: classify_ode(eq, f(x), ics=ics))

    # Does not raise
    ics = {f(0): 1}
    classify_ode(eq, f(x), ics=ics)


    #####################
    # f'(0) type (Subs) #
    #####################

    # Wrong function
    ics = {g(x).diff(x).subs(x, 0): 1}
    raises(ValueError, lambda: classify_ode(eq, f(x), ics=ics))

    # Contains x
    ics = {f(y).diff(y).subs(y, x): 1}
    raises(ValueError, lambda: classify_ode(eq, f(x), ics=ics))

    # Wrong variable
    ics = {f(y).diff(y).subs(y, 0): 1}
    raises(ValueError, lambda: classify_ode(eq, f(x), ics=ics))

    # Too many args
    ics = {f(x, y).diff(x).subs(x, 0): 1}
    raises(ValueError, lambda: classify_ode(eq, f(x), ics=ics))

    # Derivative wrt wrong vars
    ics = {Derivative(f(x), x, y).subs(x, 0): 1}
    raises(ValueError, lambda: classify_ode(eq, f(x), ics=ics))

    # point contains f
    # XXX: Should be NotImplementedError
    ics = {f(x).diff(x).subs(x, 0): f(0)}
    raises(ValueError, lambda: classify_ode(eq, f(x), ics=ics))

    # Does not raise
    ics = {f(x).diff(x).subs(x, 0): 1}
    classify_ode(eq, f(x), ics=ics)

    ###########################
    # f'(y) type (Derivative) #
    ###########################

    # Wrong function
    ics = {g(x).diff(x).subs(x, y): 1}
    raises(ValueError, lambda: classify_ode(eq, f(x), ics=ics))

    # Contains x
    ics = {f(y).diff(y).subs(y, x): 1}
    raises(ValueError, lambda: classify_ode(eq, f(x), ics=ics))

    # Too many args
    ics = {f(x, y).diff(x).subs(x, y): 1}
    raises(ValueError, lambda: classify_ode(eq, f(x), ics=ics))

    # Derivative wrt wrong vars
    ics = {Derivative(f(x), x, z).subs(x, y): 1}
    raises(ValueError, lambda: classify_ode(eq, f(x), ics=ics))

    # point contains f
    # XXX: Should be NotImplementedError
    ics = {f(x).diff(x).subs(x, y): f(0)}
    raises(ValueError, lambda: classify_ode(eq, f(x), ics=ics))

    # Does not raise
    ics = {f(x).diff(x).subs(x, y): 1}
    classify_ode(eq, f(x), ics=ics)

def test_classify_sysode():
    from sympy.matrices.dense import Matrix
    # Here x is assumed to be x(t) and y as y(t) for simplicity.
    # Similarly diff(x,t) and diff(y,y) is assumed to be x1 and y1 respectively.
    k, l, m, n = symbols('k, l, m, n', Integer=True)
    k1, k2, k3, l1, l2, l3, m1, m2, m3 = symbols('k1, k2, k3, l1, l2, l3, m1, m2, m3', Integer=True)
    P, Q, R, p, q, r = symbols('P, Q, R, p, q, r', cls=Function)
    P1, P2, P3, Q1, Q2, R1, R2 = symbols('P1, P2, P3, Q1, Q2, R1, R2', cls=Function)
    x, y, z = symbols('x, y, z', cls=Function)
    t = symbols('t')
    x1 = diff(x(t),t) ; y1 = diff(y(t),t) ; z1 = diff(z(t),t)
    x2 = diff(x(t),t,t) ; y2 = diff(y(t),t,t)

    eq1 = (Eq(diff(x(t),t), 5*t*x(t) + 2*y(t)), Eq(diff(y(t),t), 2*x(t) + 5*t*y(t)))
    sol1 = {'no_of_equation': 2, 'func_coeff': {(0, x(t), 0): -5*t, (1, x(t), 1): 0, (0, x(t), 1): 1, \
    (1, y(t), 0): -5*t, (1, x(t), 0): -2, (0, y(t), 1): 0, (0, y(t), 0): -2, (1, y(t), 1): 1}, \
    'type_of_equation': 'type3', 'func': [x(t), y(t)], 'is_linear': True, 'eq': [-5*t*x(t) - 2*y(t) + \
    Derivative(x(t), t), -5*t*y(t) - 2*x(t) + Derivative(y(t), t)], 'order': {y(t): 1, x(t): 1}}
    assert classify_sysode(eq1) == sol1

    eq2 = (Eq(x2, k*x(t) - l*y1), Eq(y2, l*x1 + k*y(t)))
    sol2 = {'order': {y(t): 2, x(t): 2}, 'type_of_equation': 'type3', 'is_linear': True, 'eq': \
    [-k*x(t) + l*Derivative(y(t), t) + Derivative(x(t), t, t), -k*y(t) - l*Derivative(x(t), t) + \
    Derivative(y(t), t, t)], 'no_of_equation': 2, 'func_coeff': {(0, y(t), 0): 0, (0, x(t), 2): 1, \
    (1, y(t), 1): 0, (1, y(t), 2): 1, (1, x(t), 2): 0, (0, y(t), 2): 0, (0, x(t), 0): -k, (1, x(t), 1): \
    -l, (0, x(t), 1): 0, (0, y(t), 1): l, (1, x(t), 0): 0, (1, y(t), 0): -k}, 'func': [x(t), y(t)]}
    assert classify_sysode(eq2) == sol2

    eq3 = (Eq(x2+4*x1+3*y1+9*x(t)+7*y(t), 11*exp(I*t)), Eq(y2+5*x1+8*y1+3*x(t)+12*y(t), 2*exp(I*t)))
    sol3 = {'no_of_equation': 2, 'func_coeff': {(1, x(t), 2): 0, (0, y(t), 2): 0, (0, x(t), 0): 9, \
    (1, x(t), 1): 5, (0, x(t), 1): 4, (0, y(t), 1): 3, (1, x(t), 0): 3, (1, y(t), 0): 12, (0, y(t), 0): 7, \
    (0, x(t), 2): 1, (1, y(t), 2): 1, (1, y(t), 1): 8}, 'type_of_equation': 'type4', 'func': [x(t), y(t)], \
    'is_linear': True, 'eq': [9*x(t) + 7*y(t) - 11*exp(I*t) + 4*Derivative(x(t), t) + 3*Derivative(y(t), t) + \
    Derivative(x(t), t, t), 3*x(t) + 12*y(t) - 2*exp(I*t) + 5*Derivative(x(t), t) + 8*Derivative(y(t), t) + \
    Derivative(y(t), t, t)], 'order': {y(t): 2, x(t): 2}}
    assert classify_sysode(eq3) == sol3

    eq4 = (Eq((4*t**2 + 7*t + 1)**2*x2, 5*x(t) + 35*y(t)), Eq((4*t**2 + 7*t + 1)**2*y2, x(t) + 9*y(t)))
    sol4 = {'no_of_equation': 2, 'func_coeff': {(1, x(t), 2): 0, (0, y(t), 2): 0, (0, x(t), 0): -5, \
    (1, x(t), 1): 0, (0, x(t), 1): 0, (0, y(t), 1): 0, (1, x(t), 0): -1, (1, y(t), 0): -9, (0, y(t), 0): -35, \
    (0, x(t), 2): 16*t**4 + 56*t**3 + 57*t**2 + 14*t + 1, (1, y(t), 2): 16*t**4 + 56*t**3 + 57*t**2 + 14*t + 1, \
    (1, y(t), 1): 0}, 'type_of_equation': 'type10', 'func': [x(t), y(t)], 'is_linear': True, \
    'eq': [(4*t**2 + 7*t + 1)**2*Derivative(x(t), t, t) - 5*x(t) - 35*y(t), (4*t**2 + 7*t + 1)**2*Derivative(y(t), t, t)\
    - x(t) - 9*y(t)], 'order': {y(t): 2, x(t): 2}}
    assert classify_sysode(eq4) == sol4

    eq5 = (Eq(diff(x(t),t), x(t) + y(t) + 9), Eq(diff(y(t),t), 2*x(t) + 5*y(t) + 23))
    sol5 = {'no_of_equation': 2, 'func_coeff': {(0, x(t), 0): -1, (1, x(t), 1): 0, (0, x(t), 1): 1, (1, y(t), 0): -5, \
    (1, x(t), 0): -2, (0, y(t), 1): 0, (0, y(t), 0): -1, (1, y(t), 1): 1}, 'type_of_equation': 'type2', \
    'func': [x(t), y(t)], 'is_linear': True, 'eq': [-x(t) - y(t) + Derivative(x(t), t) - 9, -2*x(t) - 5*y(t) + \
    Derivative(y(t), t) - 23], 'order': {y(t): 1, x(t): 1}}
    assert classify_sysode(eq5) == sol5

    eq6 = (Eq(x1, exp(k*x(t))*P(x(t),y(t))), Eq(y1,r(y(t))*P(x(t),y(t))))
    sol6 = {'no_of_equation': 2, 'func_coeff': {(0, x(t), 0): 0, (1, x(t), 1): 0, (0, x(t), 1): 1, (1, y(t), 0): 0, \
    (1, x(t), 0): 0, (0, y(t), 1): 0, (0, y(t), 0): 0, (1, y(t), 1): 1}, 'type_of_equation': 'type2', 'func': \
    [x(t), y(t)], 'is_linear': False, 'eq': [-P(x(t), y(t))*exp(k*x(t)) + Derivative(x(t), t), -P(x(t), \
    y(t))*r(y(t)) + Derivative(y(t), t)], 'order': {y(t): 1, x(t): 1}}
    assert classify_sysode(eq6) == sol6

    eq7 = (Eq(x1, x(t)**2+y(t)/x(t)), Eq(y1, x(t)/y(t)))
    sol7 = {'no_of_equation': 2, 'func_coeff': {(0, x(t), 0): 0, (1, x(t), 1): 0, (0, x(t), 1): 1, (1, y(t), 0): 0, \
    (1, x(t), 0): -1/y(t), (0, y(t), 1): 0, (0, y(t), 0): -1/x(t), (1, y(t), 1): 1}, 'type_of_equation': 'type3', \
    'func': [x(t), y(t)], 'is_linear': False, 'eq': [-x(t)**2 + Derivative(x(t), t) - y(t)/x(t), -x(t)/y(t) + \
    Derivative(y(t), t)], 'order': {y(t): 1, x(t): 1}}
    assert classify_sysode(eq7) == sol7

    eq8 = (Eq(x1, P1(x(t))*Q1(y(t))*R(x(t),y(t),t)), Eq(y1, P1(x(t))*Q1(y(t))*R(x(t),y(t),t)))
    sol8 = {'func': [x(t), y(t)], 'is_linear': False, 'type_of_equation': 'type4', 'eq': \
    [-P1(x(t))*Q1(y(t))*R(x(t), y(t), t) + Derivative(x(t), t), -P1(x(t))*Q1(y(t))*R(x(t), y(t), t) + \
    Derivative(y(t), t)], 'func_coeff': {(0, y(t), 1): 0, (1, y(t), 1): 1, (1, x(t), 1): 0, (0, y(t), 0): 0, \
    (1, x(t), 0): 0, (0, x(t), 0): 0, (1, y(t), 0): 0, (0, x(t), 1): 1}, 'order': {y(t): 1, x(t): 1}, 'no_of_equation': 2}
    assert classify_sysode(eq8) == sol8

    eq9 = (Eq(x1,3*y(t)-11*z(t)),Eq(y1,7*z(t)-3*x(t)),Eq(z1,11*x(t)-7*y(t)))
    sol9 = {'no_of_equation': 3, 'eq': [-3*y(t) + 11*z(t) + Derivative(x(t), t), 3*x(t) - 7*z(t) + Derivative(y(t), t),
            -11*x(t) + 7*y(t) + Derivative(z(t), t)], 'func': [x(t), y(t), z(t)], 'order': {x(t): 1, y(t): 1, z(t): 1},
            'is_linear': True, 'is_constant': True, 'is_homogeneous': True, 'func_coeff': Matrix([
            [  0, -3, 11], [  3,  0, -7], [-11,  7,  0]]), 'type_of_equation': 'type1'}
    assert classify_sysode(eq9) == sol9

    eq10 = (x2 + log(t)*(t*x1 - x(t)) + exp(t)*(t*y1 - y(t)), y2 + (t**2)*(t*x1 - x(t)) + (t)*(t*y1 - y(t)))
    sol10 = {'no_of_equation': 2, 'func_coeff': {(1, x(t), 2): 0, (0, y(t), 2): 0, (0, x(t), 0): -log(t), \
    (1, x(t), 1): t**3, (0, x(t), 1): t*log(t), (0, y(t), 1): t*exp(t), (1, x(t), 0): -t**2, (1, y(t), 0): -t, \
    (0, y(t), 0): -exp(t), (0, x(t), 2): 1, (1, y(t), 2): 1, (1, y(t), 1): t**2}, 'type_of_equation': 'type11', \
    'func': [x(t), y(t)], 'is_linear': True, 'eq': [(t*Derivative(x(t), t) - x(t))*log(t) + (t*Derivative(y(t), t) - \
    y(t))*exp(t) + Derivative(x(t), t, t), t**2*(t*Derivative(x(t), t) - x(t)) + t*(t*Derivative(y(t), t) - y(t)) \
    + Derivative(y(t), t, t)], 'order': {y(t): 2, x(t): 2}}
    assert classify_sysode(eq10) == sol10

    eq11 = (Eq(x1,x(t)*y(t)**3), Eq(y1,y(t)**5))
    sol11 = {'no_of_equation': 2, 'func_coeff': {(0, x(t), 0): -y(t)**3, (1, x(t), 1): 0, (0, x(t), 1): 1, \
    (1, y(t), 0): 0, (1, x(t), 0): 0, (0, y(t), 1): 0, (0, y(t), 0): 0, (1, y(t), 1): 1}, 'type_of_equation': \
    'type1', 'func': [x(t), y(t)], 'is_linear': False, 'eq': [-x(t)*y(t)**3 + Derivative(x(t), t), \
    -y(t)**5 + Derivative(y(t), t)], 'order': {y(t): 1, x(t): 1}}
    assert classify_sysode(eq11) == sol11

    eq12 = (Eq(x1, y(t)), Eq(y1, x(t)))
    sol12 = {'no_of_equation': 2, 'eq': [-y(t) + Derivative(x(t), t), -x(t) + Derivative(y(t), t)],
            'func': [x(t), y(t)], 'order': {x(t): 1, y(t): 1}, 'is_linear': True, 'is_constant': True,
            'is_homogeneous': True, 'func_coeff': Matrix([[ 0, -1], [-1,  0]]), 'type_of_equation': 'type1'}
    assert classify_sysode(eq12) == sol12

    eq13 = (Eq(x1,x(t)*y(t)*sin(t)**2), Eq(y1,y(t)**2*sin(t)**2))
    sol13 = {'no_of_equation': 2, 'func_coeff': {(0, x(t), 0): -y(t)*sin(t)**2, (1, x(t), 1): 0, (0, x(t), 1): 1, \
    (1, y(t), 0): 0, (1, x(t), 0): 0, (0, y(t), 1): 0, (0, y(t), 0): -x(t)*sin(t)**2, (1, y(t), 1): 1}, \
    'type_of_equation': 'type4', 'func': [x(t), y(t)], 'is_linear': False, 'eq': [-x(t)*y(t)*sin(t)**2 + \
    Derivative(x(t), t), -y(t)**2*sin(t)**2 + Derivative(y(t), t)], 'order': {y(t): 1, x(t): 1}}
    assert classify_sysode(eq13) == sol13

    eq14 = (Eq(x1, 21*x(t)), Eq(y1, 17*x(t)+3*y(t)), Eq(z1, 5*x(t)+7*y(t)+9*z(t)))
    sol14 = {'no_of_equation': 3, 'eq': [-21*x(t) + Derivative(x(t), t), -17*x(t) - 3*y(t) + Derivative(y(t), t),
            -5*x(t) - 7*y(t) - 9*z(t) + Derivative(z(t), t)], 'func': [x(t), y(t), z(t)], 'order': {x(t): 1, y(t): 1, z(t): 1},
            'is_linear': True, 'is_constant': True, 'is_homogeneous': True, 'func_coeff': Matrix([
            [-21,  0,  0], [-17, -3,  0], [ -5, -7, -9]]), 'type_of_equation': 'type1'}
    assert classify_sysode(eq14) == sol14

    eq15 = (Eq(x1,4*x(t)+5*y(t)+2*z(t)),Eq(y1,x(t)+13*y(t)+9*z(t)),Eq(z1,32*x(t)+41*y(t)+11*z(t)))
    sol15 = {'no_of_equation': 3, 'eq': [-4*x(t) - 5*y(t) - 2*z(t) + Derivative(x(t), t),
            -x(t) - 13*y(t) - 9*z(t) + Derivative(y(t), t), -32*x(t) - 41*y(t) - 11*z(t) + Derivative(z(t), t)],
            'func': [x(t), y(t), z(t)], 'order': {x(t): 1, y(t): 1, z(t): 1}, 'is_linear': True,
            'is_constant': True, 'is_homogeneous': True, 'func_coeff': Matrix([[ -4,  -5,  -2],
            [ -1, -13,  -9], [-32, -41, -11]]), 'type_of_equation': 'type1'}
    assert classify_sysode(eq15) == sol15

    eq16 = (Eq(3*x1,4*5*(y(t)-z(t))),Eq(4*y1,3*5*(z(t)-x(t))),Eq(5*z1,3*4*(x(t)-y(t))))
    sol16 = {'no_of_equation': 3, 'eq': [-20*y(t) + 20*z(t) + 3*Derivative(x(t), t), 15*x(t) - 15*z(t) + 4*Derivative(y(t), t),
            -12*x(t) + 12*y(t) + 5*Derivative(z(t), t)], 'func': [x(t), y(t), z(t)], 'order': {x(t): 1, y(t): 1, z(t): 1},
            'is_linear': True, 'is_constant': True, 'is_homogeneous': True, 'func_coeff': Matrix([
            [    0, Rational(-20,3),  Rational(20,3)], [Rational(15,4),     0, Rational(-15,4)],
            [Rational(-12,5),  Rational(12,5),     0]]), 'type_of_equation': 'type1'}
    assert classify_sysode(eq16) == sol16

    # issue 8193: funcs parameter for classify_sysode has to actually work
    assert classify_sysode(eq1, funcs=[x(t), y(t)]) == sol1


def test_solve_ics():
    # Basic tests that things work from dsolve.
    assert dsolve(f(x).diff(x) - 1/f(x), f(x), ics={f(1): 2}) == \
        Eq(f(x), sqrt(2 * x + 2))
    assert dsolve(f(x).diff(x) - f(x), f(x), ics={f(0): 1}) == Eq(f(x), exp(x))
    assert dsolve(f(x).diff(x) - f(x), f(x), ics={f(x).diff(x).subs(x, 0): 1}) == Eq(f(x), exp(x))
    assert dsolve(f(x).diff(x, x) + f(x), f(x), ics={f(0): 1,
        f(x).diff(x).subs(x, 0): 1}) == Eq(f(x), sin(x) + cos(x))
    assert dsolve([f(x).diff(x) - f(x) + g(x), g(x).diff(x) - g(x) - f(x)],
        [f(x), g(x)], ics={f(0): 1, g(0): 0}) == [Eq(f(x), exp(x)*cos(x)),
            Eq(g(x), exp(x)*sin(x))]

    # Test cases where dsolve returns two solutions.
    eq = (x**2*f(x)**2 - x).diff(x)
    assert dsolve(eq, f(x), ics={f(1): 0}) == [Eq(f(x),
        -sqrt(x - 1)/x), Eq(f(x), sqrt(x - 1)/x)]
    assert dsolve(eq, f(x), ics={f(x).diff(x).subs(x, 1): 0}) == [Eq(f(x),
        -sqrt(x - S.Half)/x), Eq(f(x), sqrt(x - S.Half)/x)]

    eq = cos(f(x)) - (x*sin(f(x)) - f(x)**2)*f(x).diff(x)
    assert dsolve(eq, f(x),
        ics={f(0):1}, hint='1st_exact', simplify=False) == Eq(x*cos(f(x)) + f(x)**3/3, Rational(1, 3))
    assert dsolve(eq, f(x),
        ics={f(0):1}, hint='1st_exact', simplify=True) == Eq(x*cos(f(x)) + f(x)**3/3, Rational(1, 3))

    assert solve_ics([Eq(f(x), C1*exp(x))], [f(x)], [C1], {f(0): 1}) == {C1: 1}
    assert solve_ics([Eq(f(x), C1*sin(x) + C2*cos(x))], [f(x)], [C1, C2],
        {f(0): 1, f(pi/2): 1}) == {C1: 1, C2: 1}

    assert solve_ics([Eq(f(x), C1*sin(x) + C2*cos(x))], [f(x)], [C1, C2],
        {f(0): 1, f(x).diff(x).subs(x, 0): 1}) == {C1: 1, C2: 1}

    assert solve_ics([Eq(f(x), C1*sin(x) + C2*cos(x))], [f(x)], [C1, C2], {f(0): 1}) == \
        {C2: 1}

    # Some more complicated tests Refer to PR #16098

    assert set(dsolve(f(x).diff(x)*(f(x).diff(x, 2)-x), ics={f(0):0, f(x).diff(x).subs(x, 1):0})) == \
        {Eq(f(x), 0), Eq(f(x), x ** 3 / 6 - x / 2)}
    assert set(dsolve(f(x).diff(x)*(f(x).diff(x, 2)-x), ics={f(0):0})) == \
        {Eq(f(x), 0), Eq(f(x), C2*x + x**3/6)}

    K, r, f0 = symbols('K r f0')
    sol = Eq(f(x), K*f0*exp(r*x)/((-K + f0)*(f0*exp(r*x)/(-K + f0) - 1)))
    assert (dsolve(Eq(f(x).diff(x), r * f(x) * (1 - f(x) / K)), f(x), ics={f(0): f0})) == sol


    #Order dependent issues Refer to PR #16098
    assert set(dsolve(f(x).diff(x)*(f(x).diff(x, 2)-x), ics={f(x).diff(x).subs(x,0):0, f(0):0})) == \
        {Eq(f(x), 0), Eq(f(x), x ** 3 / 6)}
    assert set(dsolve(f(x).diff(x)*(f(x).diff(x, 2)-x), ics={f(0):0, f(x).diff(x).subs(x,0):0})) == \
        {Eq(f(x), 0), Eq(f(x), x ** 3 / 6)}

    # XXX: Ought to be ValueError
    raises(ValueError, lambda: solve_ics([Eq(f(x), C1*sin(x) + C2*cos(x))], [f(x)], [C1, C2], {f(0): 1, f(pi): 1}))

    # Degenerate case. f'(0) is identically 0.
    raises(ValueError, lambda: solve_ics([Eq(f(x), sqrt(C1 - x**2))], [f(x)], [C1], {f(x).diff(x).subs(x, 0): 0}))

    EI, q, L = symbols('EI q L')

    # eq = Eq(EI*diff(f(x), x, 4), q)
    sols = [Eq(f(x), C1 + C2*x + C3*x**2 + C4*x**3 + q*x**4/(24*EI))]
    funcs = [f(x)]
    constants = [C1, C2, C3, C4]
    # Test both cases, Derivative (the default from f(x).diff(x).subs(x, L)),
    # and Subs
    ics1 = {f(0): 0,
        f(x).diff(x).subs(x, 0): 0,
        f(L).diff(L, 2): 0,
        f(L).diff(L, 3): 0}
    ics2 = {f(0): 0,
        f(x).diff(x).subs(x, 0): 0,
        Subs(f(x).diff(x, 2), x, L): 0,
        Subs(f(x).diff(x, 3), x, L): 0}

    solved_constants1 = solve_ics(sols, funcs, constants, ics1)
    solved_constants2 = solve_ics(sols, funcs, constants, ics2)
    assert solved_constants1 == solved_constants2 == {
        C1: 0,
        C2: 0,
        C3: L**2*q/(4*EI),
        C4: -L*q/(6*EI)}


def test_ode_order():
    f = Function('f')
    g = Function('g')
    x = Symbol('x')
    assert ode_order(3*x*exp(f(x)), f(x)) == 0
    assert ode_order(x*diff(f(x), x) + 3*x*f(x) - sin(x)/x, f(x)) == 1
    assert ode_order(x**2*f(x).diff(x, x) + x*diff(f(x), x) - f(x), f(x)) == 2
    assert ode_order(diff(x*exp(f(x)), x, x), f(x)) == 2
    assert ode_order(diff(x*diff(x*exp(f(x)), x, x), x), f(x)) == 3
    assert ode_order(diff(f(x), x, x), g(x)) == 0
    assert ode_order(diff(f(x), x, x)*diff(g(x), x), f(x)) == 2
    assert ode_order(diff(f(x), x, x)*diff(g(x), x), g(x)) == 1
    assert ode_order(diff(x*diff(x*exp(f(x)), x, x), x), g(x)) == 0
    # issue 5835: ode_order has to also work for unevaluated derivatives
    # (ie, without using doit()).
    assert ode_order(Derivative(x*f(x), x), f(x)) == 1
    assert ode_order(x*sin(Derivative(x*f(x)**2, x, x)), f(x)) == 2
    assert ode_order(Derivative(x*Derivative(x*exp(f(x)), x, x), x), g(x)) == 0
    assert ode_order(Derivative(f(x), x, x), g(x)) == 0
    assert ode_order(Derivative(x*exp(f(x)), x, x), f(x)) == 2
    assert ode_order(Derivative(f(x), x, x)*Derivative(g(x), x), g(x)) == 1
    assert ode_order(Derivative(x*Derivative(f(x), x, x), x), f(x)) == 3
    assert ode_order(
        x*sin(Derivative(x*Derivative(f(x), x)**2, x, x)), f(x)) == 3


# In all tests below, checkodesol has the order option set to prevent
# superfluous calls to ode_order(), and the solve_for_func flag set to False
# because dsolve() already tries to solve for the function, unless the
# simplify=False option is set.
def test_old_ode_tests():
    # These are simple tests from the old ode module
    eq1 = Eq(f(x).diff(x), 0)
    eq2 = Eq(3*f(x).diff(x) - 5, 0)
    eq3 = Eq(3*f(x).diff(x), 5)
    eq4 = Eq(9*f(x).diff(x, x) + f(x), 0)
    eq5 = Eq(9*f(x).diff(x, x), f(x))
    # Type: a(x)f'(x)+b(x)*f(x)+c(x)=0
    eq6 = Eq(x**2*f(x).diff(x) + 3*x*f(x) - sin(x)/x, 0)
    eq7 = Eq(f(x).diff(x, x) - 3*diff(f(x), x) + 2*f(x), 0)
    # Type: 2nd order, constant coefficients (two real different roots)
    eq8 = Eq(f(x).diff(x, x) - 4*diff(f(x), x) + 4*f(x), 0)
    # Type: 2nd order, constant coefficients (two real equal roots)
    eq9 = Eq(f(x).diff(x, x) + 2*diff(f(x), x) + 3*f(x), 0)
    # Type: 2nd order, constant coefficients (two complex roots)
    eq10 = Eq(3*f(x).diff(x) - 1, 0)
    eq11 = Eq(x*f(x).diff(x) - 1, 0)
    sol1 = Eq(f(x), C1)
    sol2 = Eq(f(x), C1 + x*Rational(5, 3))
    sol3 = Eq(f(x), C1 + x*Rational(5, 3))
    sol4 = Eq(f(x), C1*sin(x/3) + C2*cos(x/3))
    sol5 = Eq(f(x), C1*exp(-x/3) + C2*exp(x/3))
    sol6 = Eq(f(x), (C1 - cos(x))/x**3)
    sol7 = Eq(f(x), (C1 + C2*exp(x))*exp(x))
    sol8 = Eq(f(x), (C1 + C2*x)*exp(2*x))
    sol9 = Eq(f(x), (C1*sin(x*sqrt(2)) + C2*cos(x*sqrt(2)))*exp(-x))
    sol10 = Eq(f(x), C1 + x/3)
    sol11 = Eq(f(x), C1 + log(x))
    assert dsolve(eq1) == sol1
    assert dsolve(eq1.lhs) == sol1
    assert dsolve(eq2) == sol2
    assert dsolve(eq3) == sol3
    assert dsolve(eq4) == sol4
    assert dsolve(eq5) == sol5
    assert dsolve(eq6) == sol6
    assert dsolve(eq7) == sol7
    assert dsolve(eq8) == sol8
    assert dsolve(eq9) == sol9
    assert dsolve(eq10) == sol10
    assert dsolve(eq11) == sol11
    assert checkodesol(eq1, sol1, order=1, solve_for_func=False)[0]
    assert checkodesol(eq2, sol2, order=1, solve_for_func=False)[0]
    assert checkodesol(eq3, sol3, order=1, solve_for_func=False)[0]
    assert checkodesol(eq4, sol4, order=2, solve_for_func=False)[0]
    assert checkodesol(eq5, sol5, order=2, solve_for_func=False)[0]
    assert checkodesol(eq6, sol6, order=1, solve_for_func=False)[0]
    assert checkodesol(eq7, sol7, order=2, solve_for_func=False)[0]
    assert checkodesol(eq8, sol8, order=2, solve_for_func=False)[0]
    assert checkodesol(eq9, sol9, order=2, solve_for_func=False)[0]
    assert checkodesol(eq10, sol10, order=1, solve_for_func=False)[0]
    assert checkodesol(eq11, sol11, order=1, solve_for_func=False)[0]


@slow
def test_1st_exact1():
    # Type: Exact differential equation, p(x,f) + q(x,f)*f' == 0,
    # where dp/df == dq/dx
    eq1 = sin(x)*cos(f(x)) + cos(x)*sin(f(x))*f(x).diff(x)
    eq2 = (2*x*f(x) + 1)/f(x) + (f(x) - x)/f(x)**2*f(x).diff(x)
    eq3 = 2*x + f(x)*cos(x) + (2*f(x) + sin(x) - sin(f(x)))*f(x).diff(x)
    eq4 = cos(f(x)) - (x*sin(f(x)) - f(x)**2)*f(x).diff(x)
    eq5 = 2*x*f(x) + (x**2 + f(x)**2)*f(x).diff(x)
    sol1 = [Eq(f(x), -acos(C1/cos(x)) + 2*pi), Eq(f(x), acos(C1/cos(x)))]
    sol2 = Eq(f(x), exp(C1 - x**2 + LambertW(-x*exp(-C1 + x**2))))
    sol2b = Eq(log(f(x)) + x/f(x) + x**2, C1)
    sol3 = Eq(f(x)*sin(x) + cos(f(x)) + x**2 + f(x)**2, C1)
    sol4 = Eq(x*cos(f(x)) + f(x)**3/3, C1)
    sol5 = Eq(x**2*f(x) + f(x)**3/3, C1)
    assert dsolve(eq1, f(x), hint='1st_exact') == sol1
    assert dsolve(eq2, f(x), hint='1st_exact') == sol2
    assert dsolve(eq3, f(x), hint='1st_exact') == sol3
    assert dsolve(eq4, hint='1st_exact') == sol4
    assert dsolve(eq5, hint='1st_exact', simplify=False) == sol5
    assert checkodesol(eq1, sol1, order=1, solve_for_func=False)[0]
    # issue 5080 blocks the testing of this solution
    # FIXME: assert checkodesol(eq2, sol2, order=1, solve_for_func=False)[0]
    assert checkodesol(eq2, sol2b, order=1, solve_for_func=False)[0]
    assert checkodesol(eq3, sol3, order=1, solve_for_func=False)[0]
    assert checkodesol(eq4, sol4, order=1, solve_for_func=False)[0]
    assert checkodesol(eq5, sol5, order=1, solve_for_func=False)[0]


@slow
@XFAIL
def test_1st_exact2_broken():
    """
    This is an exact equation that fails under the exact engine. It is caught
    by first order homogeneous albeit with a much contorted solution.  The
    exact engine fails because of a poorly simplified integral of q(0,y)dy,
    where q is the function multiplying f'.  The solutions should be
    Eq(sqrt(x**2+f(x)**2)**3+y**3, C1).  The equation below is
    equivalent, but it is so complex that checkodesol fails, and takes a long
    time to do so.
    """
    if ON_TRAVIS:
        skip("Too slow for travis.")
    eq = (x*sqrt(x**2 + f(x)**2) - (x**2*f(x)/(f(x) -
          sqrt(x**2 + f(x)**2)))*f(x).diff(x))
    sol = Eq(log(x),
        C1 - 9*sqrt(1 + f(x)**2/x**2)*asinh(f(x)/x)/(-27*f(x)/x +
        27*sqrt(1 + f(x)**2/x**2)) - 9*sqrt(1 + f(x)**2/x**2)*
        log(1 - sqrt(1 + f(x)**2/x**2)*f(x)/x + 2*f(x)**2/x**2)/
        (-27*f(x)/x + 27*sqrt(1 + f(x)**2/x**2)) +
        9*asinh(f(x)/x)*f(x)/(x*(-27*f(x)/x + 27*sqrt(1 + f(x)**2/x**2))) +
        9*f(x)*log(1 - sqrt(1 + f(x)**2/x**2)*f(x)/x + 2*f(x)**2/x**2)/
        (x*(-27*f(x)/x + 27*sqrt(1 + f(x)**2/x**2))))
    assert dsolve(eq) == sol # Slow
    # FIXME: Checked in test_1st_exact2_broken_check below


@slow
def test_1st_exact2_broken_check():
    # See test_1st_exact2_broken above
    eq = (x*sqrt(x**2 + f(x)**2) - (x**2*f(x)/(f(x) -
          sqrt(x**2 + f(x)**2)))*f(x).diff(x))
    sol = Eq(log(x),
        C1 - 9*sqrt(1 + f(x)**2/x**2)*asinh(f(x)/x)/(-27*f(x)/x +
        27*sqrt(1 + f(x)**2/x**2)) - 9*sqrt(1 + f(x)**2/x**2)*
        log(1 - sqrt(1 + f(x)**2/x**2)*f(x)/x + 2*f(x)**2/x**2)/
        (-27*f(x)/x + 27*sqrt(1 + f(x)**2/x**2)) +
        9*asinh(f(x)/x)*f(x)/(x*(-27*f(x)/x + 27*sqrt(1 + f(x)**2/x**2))) +
        9*f(x)*log(1 - sqrt(1 + f(x)**2/x**2)*f(x)/x + 2*f(x)**2/x**2)/
        (x*(-27*f(x)/x + 27*sqrt(1 + f(x)**2/x**2))))
    assert checkodesol(eq, sol, order=1, solve_for_func=False)[0]


def test_separable1():
    # test_separable1-5 are from Ordinary Differential Equations, Tenenbaum and
    # Pollard, pg. 55
    eq1 = f(x).diff(x) - f(x)
    eq2 = x*f(x).diff(x) - f(x)
    eq3 = f(x).diff(x) + sin(x)
    eq4 = f(x)**2 + 1 - (x**2 + 1)*f(x).diff(x)
    eq5 = f(x).diff(x)/tan(x) - f(x) - 2
    eq6 = f(x).diff(x) * (1 - sin(f(x))) - 1
    sol1 = Eq(f(x), C1*exp(x))
    sol2 = Eq(f(x), C1*x)
    sol3 = Eq(f(x), C1 + cos(x))
    sol4 = Eq(f(x), tan(C1 + atan(x)))
    sol5 = Eq(f(x), C1/cos(x) - 2)
    sol6 = Eq(-x + f(x) + cos(f(x)), C1)
    assert dsolve(eq1, hint='separable') == sol1
    assert dsolve(eq2, hint='separable') == sol2
    assert dsolve(eq3, hint='separable') == sol3
    assert dsolve(eq4, hint='separable') == sol4
    assert dsolve(eq5, hint='separable') == sol5
    assert dsolve(eq6, hint='separable') == sol6
    assert checkodesol(eq1, sol1, order=1, solve_for_func=False)[0]
    assert checkodesol(eq2, sol2, order=1, solve_for_func=False)[0]
    assert checkodesol(eq3, sol3, order=1, solve_for_func=False)[0]
    assert checkodesol(eq4, sol4, order=1, solve_for_func=False)[0]
    assert checkodesol(eq5, sol5, order=1, solve_for_func=False)[0]
    assert checkodesol(eq6, sol6, order=1, solve_for_func=False)[0]


@slow
def test_separable2():
    a = Symbol('a')
    eq6 = f(x)*x**2*f(x).diff(x) - f(x)**3 - 2*x**2*f(x).diff(x)
    eq7 = f(x)**2 - 1 - (2*f(x) + x*f(x))*f(x).diff(x)
    eq8 = x*log(x)*f(x).diff(x) + sqrt(1 + f(x)**2)
    eq9 = exp(x + 1)*tan(f(x)) + cos(f(x))*f(x).diff(x)
    eq10 = (x*cos(f(x)) + x**2*sin(f(x))*f(x).diff(x) -
            a**2*sin(f(x))*f(x).diff(x))
    sol6 = Eq(Integral((u - 2)/u**3, (u, f(x))),
        C1 + Integral(x**(-2), x))
    sol7 = Eq(-log(-1 + f(x)**2)/2, C1 - log(2 + x))
    sol8 = Eq(asinh(f(x)), C1 - log(log(x)))
    # integrate cannot handle the integral on the lhs (cos/tan)
    sol9 = Eq(Integral(cos(u)/tan(u), (u, f(x))),
        C1 + Integral(-exp(1)*exp(x), x))
    sol10 = Eq(-log(cos(f(x))), C1 - log(- a**2 + x**2)/2)
    assert dsolve(eq6, hint='separable_Integral').dummy_eq(sol6)
    assert dsolve(eq7, hint='separable', simplify=False) == sol7
    assert dsolve(eq8, hint='separable', simplify=False) == sol8
    assert dsolve(eq9, hint='separable_Integral').dummy_eq(sol9)
    assert dsolve(eq10, hint='separable', simplify=False) == sol10
    assert checkodesol(eq6, sol6, order=1, solve_for_func=False)[0]
    assert checkodesol(eq7, sol7, order=1, solve_for_func=False)[0]
    assert checkodesol(eq8, sol8, order=1, solve_for_func=False)[0]
    assert checkodesol(eq9, sol9, order=1, solve_for_func=False)[0]
    assert checkodesol(eq10, sol10, order=1, solve_for_func=False)[0]


def test_separable3():
    eq11 = f(x).diff(x) - f(x)*tan(x)
    eq12 = (x - 1)*cos(f(x))*f(x).diff(x) - 2*x*sin(f(x))
    eq13 = f(x).diff(x) - f(x)*log(f(x))/tan(x)
    sol11 = Eq(f(x), C1/cos(x))
    sol12 = Eq(log(sin(f(x))), C1 + 2*x + 2*log(x - 1))
    sol13 = Eq(log(log(f(x))), C1 + log(sin(x)))
    assert dsolve(eq11, hint='separable') == sol11
    assert dsolve(eq12, hint='separable', simplify=False) == sol12
    assert dsolve(eq13, hint='separable', simplify=False) == sol13
    assert checkodesol(eq11, sol11, order=1, solve_for_func=False)[0]
    assert checkodesol(eq12, sol12, order=1, solve_for_func=False)[0]
    assert checkodesol(eq13, sol13, order=1, solve_for_func=False)[0]


def test_separable4():
    # This has a slow integral (1/((1 + y**2)*atan(y))), so we isolate it.
    eq14 = x*f(x).diff(x) + (1 + f(x)**2)*atan(f(x))
    sol14 = Eq(log(atan(f(x))), C1 - log(x))
    assert dsolve(eq14, hint='separable', simplify=False) == sol14
    assert checkodesol(eq14, sol14, order=1, solve_for_func=False)[0]


def test_separable5():
    eq15 = f(x).diff(x) + x*(f(x) + 1)
    eq16 = exp(f(x)**2)*(x**2 + 2*x + 1) + (x*f(x) + f(x))*f(x).diff(x)
    eq17 = f(x).diff(x) + f(x)
    eq18 = sin(x)*cos(2*f(x)) + cos(x)*sin(2*f(x))*f(x).diff(x)
    eq19 = (1 - x)*f(x).diff(x) - x*(f(x) + 1)
    eq20 = f(x)*diff(f(x), x) + x - 3*x*f(x)**2
    eq21 = f(x).diff(x) - exp(x + f(x))
    sol15 = Eq(f(x), -1 + C1*exp(-x**2/2))
    sol16 = Eq(-exp(-f(x)**2)/2, C1 - x - x**2/2)
    sol17 = Eq(f(x), C1*exp(-x))
    sol18 = Eq(-log(cos(2*f(x)))/2, C1 + log(cos(x)))
    sol19 = Eq(f(x), (C1*exp(-x) - x + 1)/(x - 1))
    sol20 = Eq(log(-1 + 3*f(x)**2)/6, C1 + x**2/2)
    sol21 = Eq(-exp(-f(x)), C1 + exp(x))
    assert dsolve(eq15, hint='separable') == sol15
    assert dsolve(eq16, hint='separable', simplify=False) == sol16
    assert dsolve(eq17, hint='separable') == sol17
    assert dsolve(eq18, hint='separable', simplify=False) == sol18
    assert dsolve(eq19, hint='separable') == sol19
    assert dsolve(eq20, hint='separable', simplify=False) == sol20
    assert dsolve(eq21, hint='separable', simplify=False) == sol21
    assert checkodesol(eq15, sol15, order=1, solve_for_func=False)[0]
    assert checkodesol(eq16, sol16, order=1, solve_for_func=False)[0]
    assert checkodesol(eq17, sol17, order=1, solve_for_func=False)[0]
    assert checkodesol(eq18, sol18, order=1, solve_for_func=False)[0]
    assert checkodesol(eq19, sol19, order=1, solve_for_func=False)[0]
    assert checkodesol(eq20, sol20, order=1, solve_for_func=False)[0]
    assert checkodesol(eq21, sol21, order=1, solve_for_func=False)[0]


def test_separable_1_5_checkodesol():
    eq12 = (x - 1)*cos(f(x))*f(x).diff(x) - 2*x*sin(f(x))
    sol12 = Eq(-log(1 - cos(f(x))**2)/2, C1 - 2*x - 2*log(1 - x))
    assert checkodesol(eq12, sol12, order=1, solve_for_func=False)[0]


def test_homogeneous_order():
    assert homogeneous_order(exp(y/x) + tan(y/x), x, y) == 0
    assert homogeneous_order(x**2 + sin(x)*cos(y), x, y) is None
    assert homogeneous_order(x - y - x*sin(y/x), x, y) == 1
    assert homogeneous_order((x*y + sqrt(x**4 + y**4) + x**2*(log(x) - log(y)))/
        (pi*x**Rational(2, 3)*sqrt(y)**3), x, y) == Rational(-1, 6)
    assert homogeneous_order(y/x*cos(y/x) - x/y*sin(y/x) + cos(y/x), x, y) == 0
    assert homogeneous_order(f(x), x, f(x)) == 1
    assert homogeneous_order(f(x)**2, x, f(x)) == 2
    assert homogeneous_order(x*y*z, x, y) == 2
    assert homogeneous_order(x*y*z, x, y, z) == 3
    assert homogeneous_order(x**2*f(x)/sqrt(x**2 + f(x)**2), f(x)) is None
    assert homogeneous_order(f(x, y)**2, x, f(x, y), y) == 2
    assert homogeneous_order(f(x, y)**2, x, f(x), y) is None
    assert homogeneous_order(f(x, y)**2, x, f(x, y)) is None
    assert homogeneous_order(f(y, x)**2, x, y, f(x, y)) is None
    assert homogeneous_order(f(y), f(x), x) is None
    assert homogeneous_order(-f(x)/x + 1/sin(f(x)/ x), f(x), x) == 0
    assert homogeneous_order(log(1/y) + log(x**2), x, y) is None
    assert homogeneous_order(log(1/y) + log(x), x, y) == 0
    assert homogeneous_order(log(x/y), x, y) == 0
    assert homogeneous_order(2*log(1/y) + 2*log(x), x, y) == 0
    a = Symbol('a')
    assert homogeneous_order(a*log(1/y) + a*log(x), x, y) == 0
    assert homogeneous_order(f(x).diff(x), x, y) is None
    assert homogeneous_order(-f(x).diff(x) + x, x, y) is None
    assert homogeneous_order(O(x), x, y) is None
    assert homogeneous_order(x + O(x**2), x, y) is None
    assert homogeneous_order(x**pi, x) == pi
    assert homogeneous_order(x**x, x) is None
    raises(ValueError, lambda: homogeneous_order(x*y))


@slow
def test_1st_homogeneous_coeff_ode():
    # Type: First order homogeneous, y'=f(y/x)
    eq1 = f(x)/x*cos(f(x)/x) - (x/f(x)*sin(f(x)/x) + cos(f(x)/x))*f(x).diff(x)
    eq2 = x*f(x).diff(x) - f(x) - x*sin(f(x)/x)
    eq3 = f(x) + (x*log(f(x)/x) - 2*x)*diff(f(x), x)
    eq4 = 2*f(x)*exp(x/f(x)) + f(x)*f(x).diff(x) - 2*x*exp(x/f(x))*f(x).diff(x)
    eq5 = 2*x**2*f(x) + f(x)**3 + (x*f(x)**2 - 2*x**3)*f(x).diff(x)
    eq6 = x*exp(f(x)/x) - f(x)*sin(f(x)/x) + x*sin(f(x)/x)*f(x).diff(x)
    eq7 = (x + sqrt(f(x)**2 - x*f(x)))*f(x).diff(x) - f(x)
    eq8 = x + f(x) - (x - f(x))*f(x).diff(x)

    sol1 = Eq(log(x), C1 - log(f(x)*sin(f(x)/x)/x))
    sol2 = Eq(log(x), log(C1) + log(cos(f(x)/x) - 1)/2 - log(cos(f(x)/x) + 1)/2)
    sol3 = Eq(f(x), -exp(C1)*LambertW(-x*exp(-C1 + 1)))
    sol4 = Eq(log(f(x)), C1 - 2*exp(x/f(x)))
    sol5 = Eq(f(x), exp(2*C1 + LambertW(-2*x**4*exp(-4*C1))/2)/x)
    sol6 = Eq(log(x), C1 + exp(-f(x)/x)*sin(f(x)/x)/2 + exp(-f(x)/x)*cos(f(x)/x)/2)
    sol7 = Eq(log(f(x)), C1 - 2*sqrt(-x/f(x) + 1))
    sol8 = Eq(log(x), C1 - log(sqrt(1 + f(x)**2/x**2)) + atan(f(x)/x))

    # indep_div_dep actually has a simpler solution for eq2,
    # but it runs too slow
    assert dsolve(eq1, hint='1st_homogeneous_coeff_subs_dep_div_indep') == sol1
    assert dsolve(eq2, hint='1st_homogeneous_coeff_subs_dep_div_indep', simplify=False) == sol2
    assert dsolve(eq3, hint='1st_homogeneous_coeff_best') == sol3
    assert dsolve(eq4, hint='1st_homogeneous_coeff_best') == sol4
    assert dsolve(eq5, hint='1st_homogeneous_coeff_best') == sol5
    assert dsolve(eq6, hint='1st_homogeneous_coeff_subs_dep_div_indep') == sol6
    assert dsolve(eq7, hint='1st_homogeneous_coeff_best') == sol7
    assert dsolve(eq8, hint='1st_homogeneous_coeff_best') == sol8

    # FIXME: sol3 and sol5 don't work with checkodesol (because of LambertW?)
    # previous code was testing with these other solutions:
    sol3b = Eq(-f(x)/(1 + log(x/f(x))), C1)
    sol5b = Eq(log(C1*x*sqrt(1/x)*sqrt(f(x))) + x**2/(2*f(x)**2), 0)
    assert checkodesol(eq1, sol1, order=1, solve_for_func=False)[0]
    assert checkodesol(eq2, sol2, order=1, solve_for_func=False)[0]
    assert checkodesol(eq3, sol3b, order=1, solve_for_func=False)[0]
    assert checkodesol(eq4, sol4, order=1, solve_for_func=False)[0]
    assert checkodesol(eq5, sol5b, order=1, solve_for_func=False)[0]
    assert checkodesol(eq6, sol6, order=1, solve_for_func=False)[0]
    assert checkodesol(eq8, sol8, order=1, solve_for_func=False)[0]


def test_1st_homogeneous_coeff_ode_check2():
    eq2 = x*f(x).diff(x) - f(x) - x*sin(f(x)/x)
    sol2 = Eq(x/tan(f(x)/(2*x)), C1)
    assert checkodesol(eq2, sol2, order=1, solve_for_func=False)[0]


def test_1st_homogeneous_coeff_ode_check3():
    eq3 = f(x) + (x*log(f(x)/x) - 2*x)*diff(f(x), x)
    # This solution is correct:
    sol3 = Eq(f(x), -exp(C1)*LambertW(-x*exp(1 - C1)))
    assert dsolve(eq3) == sol3
    # FIXME: Checked in test_1st_homogeneous_coeff_ode_check3_check below

    # Alternate form:
    sol3a = Eq(f(x), x*exp(1 - LambertW(C1*x)))
    assert checkodesol(eq3, sol3a, solve_for_func=True)[0]


@XFAIL
def test_1st_homogeneous_coeff_ode_check3_check():
    # See test_1st_homogeneous_coeff_ode_check3 above
    eq3 = f(x) + (x*log(f(x)/x) - 2*x)*diff(f(x), x)
    sol3 = Eq(f(x), -exp(C1)*LambertW(-x*exp(1 - C1)))
    assert checkodesol(eq3, sol3) == (True, 0)  # XFAIL


def test_1st_homogeneous_coeff_ode_check7():
    eq7 = (x + sqrt(f(x)**2 - x*f(x)))*f(x).diff(x) - f(x)
    sol7 = Eq(log(f(x)), C1 - 2*sqrt(-x/f(x) + 1))
    assert dsolve(eq7) == sol7
    assert checkodesol(eq7, sol7, order=1, solve_for_func=False) == (True, 0)


def test_1st_homogeneous_coeff_ode2():
    eq1 = f(x).diff(x) - f(x)/x + 1/sin(f(x)/x)
    eq2 = x**2 + f(x)**2 - 2*x*f(x)*f(x).diff(x)
    eq3 = x*exp(f(x)/x) + f(x) - x*f(x).diff(x)
    sol1 = [Eq(f(x), x*(-acos(C1 + log(x)) + 2*pi)), Eq(f(x), x*acos(C1 + log(x)))]
    sol2 = Eq(log(f(x)), log(C1) + log(x/f(x)) - log(x**2/f(x)**2 - 1))
    sol3 = Eq(f(x), log((1/(C1 - log(x)))**x))
    # specific hints are applied for speed reasons
    assert dsolve(eq1, hint='1st_homogeneous_coeff_subs_dep_div_indep') == sol1
    assert dsolve(eq2, hint='1st_homogeneous_coeff_best', simplify=False) == sol2
    assert dsolve(eq3, hint='1st_homogeneous_coeff_subs_dep_div_indep') == sol3

    # FIXME: sol3 doesn't work with checkodesol (because of **x?)
    # previous code was testing with this other solution:
    sol3b = Eq(f(x), log(log(C1/x)**(-x)))
    assert checkodesol(eq1, sol1, order=1, solve_for_func=False)[0]
    assert checkodesol(eq2, sol2, order=1, solve_for_func=False)[0]
    assert checkodesol(eq3, sol3b, order=1, solve_for_func=False)[0]


def test_1st_homogeneous_coeff_ode_check9():
    _u2 = Dummy('u2')
    __a = Dummy('a')
    eq9 = f(x)**2 + (x*sqrt(f(x)**2 - x**2) - x*f(x))*f(x).diff(x)
    sol9 = Eq(-Integral(-1/(-(1 - sqrt(1 - _u2**2))*_u2 + _u2), (_u2, __a,
        x/f(x))) + log(C1*f(x)), 0)
    assert checkodesol(eq9, sol9, order=1, solve_for_func=False)[0]


def test_1st_homogeneous_coeff_ode3():
    # The standard integration engine cannot handle one of the integrals
    # involved (see issue 4551).  meijerg code comes up with an answer, but in
    # unconventional form.
    # checkodesol fails for this equation, so its test is in
    # test_1st_homogeneous_coeff_ode_check9 above. It has to compare string
    # expressions because u2 is a dummy variable.
    eq = f(x)**2 + (x*sqrt(f(x)**2 - x**2) - x*f(x))*f(x).diff(x)
    sol = Eq(log(f(x)), C1 + Piecewise(
            (acosh(f(x)/x), abs(f(x)**2)/x**2 > 1),
            (-I*asin(f(x)/x), True)))
    assert dsolve(eq, hint='1st_homogeneous_coeff_subs_indep_div_dep') == sol


def test_1st_homogeneous_coeff_corner_case():
    eq1 = f(x).diff(x) - f(x)/x
    c1 = classify_ode(eq1, f(x))
    eq2 = x*f(x).diff(x) - f(x)
    c2 = classify_ode(eq2, f(x))
    sdi = "1st_homogeneous_coeff_subs_dep_div_indep"
    sid = "1st_homogeneous_coeff_subs_indep_div_dep"
    assert sid not in c1 and sdi not in c1
    assert sid not in c2 and sdi not in c2


@slow
def test_nth_linear_constant_coeff_homogeneous():
    # From Exercise 20, in Ordinary Differential Equations,
    #                      Tenenbaum and Pollard, pg. 220
    a = Symbol('a', positive=True)
    k = Symbol('k', real=True)
    eq1 = f(x).diff(x, 2) + 2*f(x).diff(x)
    eq2 = f(x).diff(x, 2) - 3*f(x).diff(x) + 2*f(x)
    eq3 = f(x).diff(x, 2) - f(x)
    eq4 = f(x).diff(x, 3) + f(x).diff(x, 2) - 6*f(x).diff(x)
    eq5 = 6*f(x).diff(x, 2) - 11*f(x).diff(x) + 4*f(x)
    eq6 = Eq(f(x).diff(x, 2) + 2*f(x).diff(x) - f(x), 0)
    eq7 = diff(f(x), x, 3) + diff(f(x), x, 2) - 10*diff(f(x), x) - 6*f(x)
    eq8 = f(x).diff(x, 4) - f(x).diff(x, 3) - 4*f(x).diff(x, 2) + \
        4*f(x).diff(x)
    eq9 = f(x).diff(x, 4) + 4*f(x).diff(x, 3) + f(x).diff(x, 2) - \
        4*f(x).diff(x) - 2*f(x)
    eq10 = f(x).diff(x, 4) - a**2*f(x)
    eq11 = f(x).diff(x, 2) - 2*k*f(x).diff(x) - 2*f(x)
    eq12 = f(x).diff(x, 2) + 4*k*f(x).diff(x) - 12*k**2*f(x)
    eq13 = f(x).diff(x, 4)
    eq14 = f(x).diff(x, 2) + 4*f(x).diff(x) + 4*f(x)
    eq15 = 3*f(x).diff(x, 3) + 5*f(x).diff(x, 2) + f(x).diff(x) - f(x)
    eq16 = f(x).diff(x, 3) - 6*f(x).diff(x, 2) + 12*f(x).diff(x) - 8*f(x)
    eq17 = f(x).diff(x, 2) - 2*a*f(x).diff(x) + a**2*f(x)
    eq18 = f(x).diff(x, 4) + 3*f(x).diff(x, 3)
    eq19 = f(x).diff(x, 4) - 2*f(x).diff(x, 2)
    eq20 = f(x).diff(x, 4) + 2*f(x).diff(x, 3) - 11*f(x).diff(x, 2) - \
        12*f(x).diff(x) + 36*f(x)
    eq21 = 36*f(x).diff(x, 4) - 37*f(x).diff(x, 2) + 4*f(x).diff(x) + 5*f(x)
    eq22 = f(x).diff(x, 4) - 8*f(x).diff(x, 2) + 16*f(x)
    eq23 = f(x).diff(x, 2) - 2*f(x).diff(x) + 5*f(x)
    eq24 = f(x).diff(x, 2) - f(x).diff(x) + f(x)
    eq25 = f(x).diff(x, 4) + 5*f(x).diff(x, 2) + 6*f(x)
    eq26 = f(x).diff(x, 2) - 4*f(x).diff(x) + 20*f(x)
    eq27 = f(x).diff(x, 4) + 4*f(x).diff(x, 2) + 4*f(x)
    eq28 = f(x).diff(x, 3) + 8*f(x)
    eq29 = f(x).diff(x, 4) + 4*f(x).diff(x, 2)
    eq30 = f(x).diff(x, 5) + 2*f(x).diff(x, 3) + f(x).diff(x)
    eq31 = f(x).diff(x, 4) + f(x).diff(x, 2) + f(x)
    eq32 = f(x).diff(x, 4) + 4*f(x).diff(x, 2) + f(x)
    sol1 = Eq(f(x), C1 + C2*exp(-2*x))
    sol2 = Eq(f(x), (C1 + C2*exp(x))*exp(x))
    sol3 = Eq(f(x), C1*exp(x) + C2*exp(-x))
    sol4 = Eq(f(x), C1 + C2*exp(-3*x) + C3*exp(2*x))
    sol5 = Eq(f(x), C1*exp(x/2) + C2*exp(x*Rational(4, 3)))
    sol6 = Eq(f(x), C1*exp(x*(-1 + sqrt(2))) + C2*exp(x*(-sqrt(2) - 1)))
    sol7 = Eq(f(x), C3*exp(3*x) + (C1*exp(-sqrt(2)*x) + C2*exp(sqrt(2)*x))*exp(-2*x))
    sol8 = Eq(f(x), C1 + C2*exp(x) + C3*exp(-2*x) + C4*exp(2*x))
    sol9 = Eq(f(x), C3*exp(-x) + C4*exp(x) + (C1*exp(-sqrt(2)*x) + C2*exp(sqrt(2)*x))*exp(-2*x))
    sol10 = Eq(f(x),
        C1*sin(x*sqrt(a)) + C2*cos(x*sqrt(a)) + C3*exp(x*sqrt(a)) +
        C4*exp(-x*sqrt(a)))
    sol11 = Eq(f(x),
        C1*exp(x*(k - sqrt(k**2 + 2))) + C2*exp(x*(k + sqrt(k**2 + 2))))
    sol12 = Eq(f(x), C1*exp(-6*k*x) + C2*exp(2*k*x))
    sol13 = Eq(f(x), C1 + C2*x + C3*x**2 + C4*x**3)
    sol14 = Eq(f(x), (C1 + C2*x)*exp(-2*x))
    sol15 = Eq(f(x), (C1 + C2*x)*exp(-x) + C3*exp(x/3))
    sol16 = Eq(f(x), (C1 + x*(C2 + C3*x))*exp(2*x))
    sol17 = Eq(f(x), (C1 + C2*x)*exp(a*x))
    sol18 = Eq(f(x), C1 + C2*x + C3*x**2 + C4*exp(-3*x))
    sol19 = Eq(f(x), C1 + C2*x + C3*exp(x*sqrt(2)) + C4*exp(-x*sqrt(2)))
    sol20 = Eq(f(x), (C1 + C2*x)*exp(-3*x) + (C3 + C4*x)*exp(2*x))
    sol21 = Eq(f(x), C1*exp(x/2) + C2*exp(-x) + C3*exp(-x/3) + C4*exp(x*Rational(5, 6)))
    sol22 = Eq(f(x), (C1 + C2*x)*exp(-2*x) + (C3 + C4*x)*exp(2*x))
    sol23 = Eq(f(x), (C1*sin(2*x) + C2*cos(2*x))*exp(x))
    sol24 = Eq(f(x), (C1*sin(x*sqrt(3)/2) + C2*cos(x*sqrt(3)/2))*exp(x/2))
    sol25 = Eq(f(x),
        C1*cos(x*sqrt(3)) + C2*sin(x*sqrt(3)) + C3*sin(x*sqrt(2)) +
        C4*cos(x*sqrt(2)))
    sol26 = Eq(f(x), (C1*sin(4*x) + C2*cos(4*x))*exp(2*x))
    sol27 = Eq(f(x), (C1 + C2*x)*sin(x*sqrt(2)) + (C3 + C4*x)*cos(x*sqrt(2)))
    sol28 = Eq(f(x),
        (C1*sin(x*sqrt(3)) + C2*cos(x*sqrt(3)))*exp(x) + C3*exp(-2*x))
    sol29 = Eq(f(x), C1 + C2*sin(2*x) + C3*cos(2*x) + C4*x)
    sol30 = Eq(f(x), C1 + (C2 + C3*x)*sin(x) + (C4 + C5*x)*cos(x))
    sol31 = Eq(f(x), (C1*sin(sqrt(3)*x/2) + C2*cos(sqrt(3)*x/2))/sqrt(exp(x))
                   + (C3*sin(sqrt(3)*x/2) + C4*cos(sqrt(3)*x/2))*sqrt(exp(x)))
    sol32 = Eq(f(x), C1*sin(x*sqrt(-sqrt(3) + 2)) + C2*sin(x*sqrt(sqrt(3) + 2))
                   + C3*cos(x*sqrt(-sqrt(3) + 2)) + C4*cos(x*sqrt(sqrt(3) + 2)))
    sol1s = constant_renumber(sol1)
    sol2s = constant_renumber(sol2)
    sol3s = constant_renumber(sol3)
    sol4s = constant_renumber(sol4)
    sol5s = constant_renumber(sol5)
    sol6s = constant_renumber(sol6)
    sol7s = constant_renumber(sol7)
    sol8s = constant_renumber(sol8)
    sol9s = constant_renumber(sol9)
    sol10s = constant_renumber(sol10)
    sol11s = constant_renumber(sol11)
    sol12s = constant_renumber(sol12)
    sol13s = constant_renumber(sol13)
    sol14s = constant_renumber(sol14)
    sol15s = constant_renumber(sol15)
    sol16s = constant_renumber(sol16)
    sol17s = constant_renumber(sol17)
    sol18s = constant_renumber(sol18)
    sol19s = constant_renumber(sol19)
    sol20s = constant_renumber(sol20)
    sol21s = constant_renumber(sol21)
    sol22s = constant_renumber(sol22)
    sol23s = constant_renumber(sol23)
    sol24s = constant_renumber(sol24)
    sol25s = constant_renumber(sol25)
    sol26s = constant_renumber(sol26)
    sol27s = constant_renumber(sol27)
    sol28s = constant_renumber(sol28)
    sol29s = constant_renumber(sol29)
    sol30s = constant_renumber(sol30)
    assert dsolve(eq1) in (sol1, sol1s)
    assert dsolve(eq2) in (sol2, sol2s)
    assert dsolve(eq3) in (sol3, sol3s)
    assert dsolve(eq4) in (sol4, sol4s)
    assert dsolve(eq5) in (sol5, sol5s)
    assert dsolve(eq6) in (sol6, sol6s)
    got = dsolve(eq7)
    assert got in (sol7, sol7s), got
    assert dsolve(eq8) in (sol8, sol8s)
    got = dsolve(eq9)
    assert got in (sol9, sol9s), got
    assert dsolve(eq10) in (sol10, sol10s)
    assert dsolve(eq11) in (sol11, sol11s)
    assert dsolve(eq12) in (sol12, sol12s)
    assert dsolve(eq13) in (sol13, sol13s)
    assert dsolve(eq14) in (sol14, sol14s)
    assert dsolve(eq15) in (sol15, sol15s)
    got = dsolve(eq16)
    assert got in (sol16, sol16s), got
    assert dsolve(eq17) in (sol17, sol17s)
    assert dsolve(eq18) in (sol18, sol18s)
    assert dsolve(eq19) in (sol19, sol19s)
    assert dsolve(eq20) in (sol20, sol20s)
    assert dsolve(eq21) in (sol21, sol21s)
    assert dsolve(eq22) in (sol22, sol22s)
    assert dsolve(eq23) in (sol23, sol23s)
    assert dsolve(eq24) in (sol24, sol24s)
    assert dsolve(eq25) in (sol25, sol25s)
    assert dsolve(eq26) in (sol26, sol26s)
    assert dsolve(eq27) in (sol27, sol27s)
    assert dsolve(eq28) in (sol28, sol28s)
    assert dsolve(eq29) in (sol29, sol29s)
    assert dsolve(eq30) in (sol30, sol30s)
    assert dsolve(eq31) in (sol31,)
    assert dsolve(eq32) in (sol32,)
    assert checkodesol(eq1, sol1, order=2, solve_for_func=False)[0]
    assert checkodesol(eq2, sol2, order=2, solve_for_func=False)[0]
    assert checkodesol(eq3, sol3, order=2, solve_for_func=False)[0]
    assert checkodesol(eq4, sol4, order=3, solve_for_func=False)[0]
    assert checkodesol(eq5, sol5, order=2, solve_for_func=False)[0]
    assert checkodesol(eq6, sol6, order=2, solve_for_func=False)[0]
    assert checkodesol(eq7, sol7, order=3, solve_for_func=False)[0]
    assert checkodesol(eq8, sol8, order=4, solve_for_func=False)[0]
    assert checkodesol(eq9, sol9, order=4, solve_for_func=False)[0]
    assert checkodesol(eq10, sol10, order=4, solve_for_func=False)[0]
    assert checkodesol(eq11, sol11, order=2, solve_for_func=False)[0]
    assert checkodesol(eq12, sol12, order=2, solve_for_func=False)[0]
    assert checkodesol(eq13, sol13, order=4, solve_for_func=False)[0]
    assert checkodesol(eq14, sol14, order=2, solve_for_func=False)[0]
    assert checkodesol(eq15, sol15, order=3, solve_for_func=False)[0]
    assert checkodesol(eq16, sol16, order=3, solve_for_func=False)[0]
    assert checkodesol(eq17, sol17, order=2, solve_for_func=False)[0]
    assert checkodesol(eq18, sol18, order=4, solve_for_func=False)[0]
    assert checkodesol(eq19, sol19, order=4, solve_for_func=False)[0]
    assert checkodesol(eq20, sol20, order=4, solve_for_func=False)[0]
    assert checkodesol(eq21, sol21, order=4, solve_for_func=False)[0]
    assert checkodesol(eq22, sol22, order=4, solve_for_func=False)[0]
    assert checkodesol(eq23, sol23, order=2, solve_for_func=False)[0]
    assert checkodesol(eq24, sol24, order=2, solve_for_func=False)[0]
    assert checkodesol(eq25, sol25, order=4, solve_for_func=False)[0]
    assert checkodesol(eq26, sol26, order=2, solve_for_func=False)[0]
    assert checkodesol(eq27, sol27, order=4, solve_for_func=False)[0]
    assert checkodesol(eq28, sol28, order=3, solve_for_func=False)[0]
    assert checkodesol(eq29, sol29, order=4, solve_for_func=False)[0]
    assert checkodesol(eq30, sol30, order=5, solve_for_func=False)[0]
    assert checkodesol(eq31, sol31, order=4, solve_for_func=False)[0]
    assert checkodesol(eq32, sol32, order=4, solve_for_func=False)[0]

    # Issue #15237
    eqn = Derivative(x*f(x), x, x, x)
    hint = 'nth_linear_constant_coeff_homogeneous'
    raises(ValueError, lambda: dsolve(eqn, f(x), hint, prep=True))
    raises(ValueError, lambda: dsolve(eqn, f(x), hint, prep=False))


def test_nth_linear_constant_coeff_homogeneous_rootof():
    # One real root, two complex conjugate pairs
    eq = f(x).diff(x, 5) + 11*f(x).diff(x) - 2*f(x)
    r1, r2, r3, r4, r5 = [rootof(x**5 + 11*x - 2, n) for n in range(5)]
    sol = Eq(f(x),
            C5*exp(r1*x)
            + exp(re(r2)*x) * (C1*sin(im(r2)*x) + C2*cos(im(r2)*x))
            + exp(re(r4)*x) * (C3*sin(im(r4)*x) + C4*cos(im(r4)*x))
            )
    assert dsolve(eq) == sol
    # FIXME: assert checkodesol(eq, sol) == (True, [0])  # Hangs...

    # Three real roots, one complex conjugate pair
    eq = f(x).diff(x,5) - 3*f(x).diff(x) + f(x)
    r1, r2, r3, r4, r5 = [rootof(x**5 - 3*x + 1, n) for n in range(5)]
    sol = Eq(f(x),
            C3*exp(r1*x) + C4*exp(r2*x) + C5*exp(r3*x)
            + exp(re(r4)*x) * (C1*sin(im(r4)*x) + C2*cos(im(r4)*x))
            )
    assert dsolve(eq) == sol
    # FIXME: assert checkodesol(eq, sol) == (True, [0])  # Hangs...

    # Five distinct real roots
    eq = f(x).diff(x,5) - 100*f(x).diff(x,3) + 1000*f(x).diff(x) + f(x)
    r1, r2, r3, r4, r5 = [rootof(x**5 - 100*x**3 + 1000*x + 1, n) for n in range(5)]
    sol = Eq(f(x), C1*exp(r1*x) + C2*exp(r2*x) + C3*exp(r3*x) + C4*exp(r4*x) + C5*exp(r5*x))
    assert dsolve(eq) == sol
    # FIXME: assert checkodesol(eq, sol) == (True, [0])  # Hangs...

    # Rational root and unsolvable quintic
    eq = f(x).diff(x, 6) - 6*f(x).diff(x, 5) + 5*f(x).diff(x, 4) + 10*f(x).diff(x) - 50 * f(x)
    r2, r3, r4, r5, r6 = [rootof(x**5 - x**4 + 10, n) for n in range(5)]
    sol = Eq(f(x),
          C5*exp(5*x)
        + C6*exp(x*r2)
        + exp(re(r3)*x) * (C1*sin(im(r3)*x) + C2*cos(im(r3)*x))
        + exp(re(r5)*x) * (C3*sin(im(r5)*x) + C4*cos(im(r5)*x))
            )
    assert dsolve(eq) == sol
    # FIXME: assert checkodesol(eq, sol) == (True, [0])  # Hangs...

    # Five double roots (this is (x**5 - x + 1)**2)
    eq = f(x).diff(x, 10) - 2*f(x).diff(x, 6) + 2*f(x).diff(x, 5) + f(x).diff(x, 2) - 2*f(x).diff(x, 1) + f(x)
    r1, r2, r3, r4, r5 = [rootof(x**5 - x + 1, n) for n in range(5)]
    sol = Eq(f(x), (C1 + C2*x)*exp(x*r1) + (C10*sin(x*im(r4)) + C7*x*sin(x*im(r4)) + (
        C8 + C9*x)*cos(x*im(r4)))*exp(x*re(r4)) + (C3*x*sin(x*im(r2)) + C6*sin(x*im(r2)
        ) + (C4 + C5*x)*cos(x*im(r2)))*exp(x*re(r2)))
    got = dsolve(eq)
    assert sol == got, got
    # FIXME: assert checkodesol(eq, sol) == (True, [0])  # Hangs...


def test_nth_linear_constant_coeff_homogeneous_irrational():
    our_hint='nth_linear_constant_coeff_homogeneous'

    eq = Eq(sqrt(2) * f(x).diff(x,x,x) + f(x).diff(x), 0)
    sol = Eq(f(x), C1 + C2*sin(2**Rational(3, 4)*x/2) + C3*cos(2**Rational(3, 4)*x/2))
    assert our_hint in classify_ode(eq)
    assert dsolve(eq, f(x), hint=our_hint) == sol
    assert dsolve(eq, f(x)) == sol
    assert checkodesol(eq, sol, order=3, solve_for_func=False)[0]

    E = exp(1)
    eq = Eq(E * f(x).diff(x,x,x) + f(x).diff(x), 0)
    sol = Eq(f(x), C1 + C2*sin(x/sqrt(E)) + C3*cos(x/sqrt(E)))
    assert our_hint in classify_ode(eq)
    assert dsolve(eq, f(x), hint=our_hint) == sol
    assert dsolve(eq, f(x)) == sol
    assert checkodesol(eq, sol, order=3, solve_for_func=False)[0]

    eq = Eq(pi * f(x).diff(x,x,x) + f(x).diff(x), 0)
    sol = Eq(f(x), C1 + C2*sin(x/sqrt(pi)) + C3*cos(x/sqrt(pi)))
    assert our_hint in classify_ode(eq)
    assert dsolve(eq, f(x), hint=our_hint) == sol
    assert dsolve(eq, f(x)) == sol
    assert checkodesol(eq, sol, order=3, solve_for_func=False)[0]

    eq = Eq(I * f(x).diff(x,x,x) + f(x).diff(x), 0)
    sol = Eq(f(x), C1 + C2*exp(-sqrt(I)*x) + C3*exp(sqrt(I)*x))
    assert our_hint in classify_ode(eq)
    assert dsolve(eq, f(x), hint=our_hint) == sol
    assert dsolve(eq, f(x)) == sol
    assert checkodesol(eq, sol, order=3, solve_for_func=False)[0]


@XFAIL
@slow
def test_nth_linear_constant_coeff_homogeneous_rootof_sol():
    # See https://github.com/sympy/sympy/issues/15753
    if ON_TRAVIS:
        skip("Too slow for travis.")
    eq = f(x).diff(x, 5) + 11*f(x).diff(x) - 2*f(x)
    sol = Eq(f(x),
        C1*exp(x*rootof(x**5 + 11*x - 2, 0)) +
        C2*exp(x*rootof(x**5 + 11*x - 2, 1)) +
        C3*exp(x*rootof(x**5 + 11*x - 2, 2)) +
        C4*exp(x*rootof(x**5 + 11*x - 2, 3)) +
        C5*exp(x*rootof(x**5 + 11*x - 2, 4)))
    assert checkodesol(eq, sol, order=5, solve_for_func=False)[0]


@XFAIL
def test_noncircularized_real_imaginary_parts():
    # If this passes, lines numbered 3878-3882 (at the time of this commit)
    # of sympy/solvers/ode.py for nth_linear_constant_coeff_homogeneous
    # should be removed.
    y = sqrt(1+x)
    i, r = im(y), re(y)
    assert not (i.has(atan2) and r.has(atan2))


def test_collect_respecting_exponentials():
    # If this test passes, lines 1306-1311 (at the time of this commit)
    # of sympy/solvers/ode.py should be removed.
    sol = 1 + exp(x/2)
    assert sol == collect( sol, exp(x/3))


def test_undetermined_coefficients_match():
    assert _undetermined_coefficients_match(g(x), x) == {'test': False}
    assert _undetermined_coefficients_match(sin(2*x + sqrt(5)), x) == \
        {'test': True, 'trialset':
            set([cos(2*x + sqrt(5)), sin(2*x + sqrt(5))])}
    assert _undetermined_coefficients_match(sin(x)*cos(x), x) == \
        {'test': False}
    s = set([cos(x), x*cos(x), x**2*cos(x), x**2*sin(x), x*sin(x), sin(x)])
    assert _undetermined_coefficients_match(sin(x)*(x**2 + x + 1), x) == \
        {'test': True, 'trialset': s}
    assert _undetermined_coefficients_match(
        sin(x)*x**2 + sin(x)*x + sin(x), x) == {'test': True, 'trialset': s}
    assert _undetermined_coefficients_match(
        exp(2*x)*sin(x)*(x**2 + x + 1), x
    ) == {
        'test': True, 'trialset': set([exp(2*x)*sin(x), x**2*exp(2*x)*sin(x),
        cos(x)*exp(2*x), x**2*cos(x)*exp(2*x), x*cos(x)*exp(2*x),
        x*exp(2*x)*sin(x)])}
    assert _undetermined_coefficients_match(1/sin(x), x) == {'test': False}
    assert _undetermined_coefficients_match(log(x), x) == {'test': False}
    assert _undetermined_coefficients_match(2**(x)*(x**2 + x + 1), x) == \
        {'test': True, 'trialset': set([2**x, x*2**x, x**2*2**x])}
    assert _undetermined_coefficients_match(x**y, x) == {'test': False}
    assert _undetermined_coefficients_match(exp(x)*exp(2*x + 1), x) == \
        {'test': True, 'trialset': set([exp(1 + 3*x)])}
    assert _undetermined_coefficients_match(sin(x)*(x**2 + x + 1), x) == \
        {'test': True, 'trialset': set([x*cos(x), x*sin(x), x**2*cos(x),
        x**2*sin(x), cos(x), sin(x)])}
    assert _undetermined_coefficients_match(sin(x)*(x + sin(x)), x) == \
        {'test': False}
    assert _undetermined_coefficients_match(sin(x)*(x + sin(2*x)), x) == \
        {'test': False}
    assert _undetermined_coefficients_match(sin(x)*tan(x), x) == \
        {'test': False}
    assert _undetermined_coefficients_match(
        x**2*sin(x)*exp(x) + x*sin(x) + x, x
    ) == {
        'test': True, 'trialset': set([x**2*cos(x)*exp(x), x, cos(x), S.One,
        exp(x)*sin(x), sin(x), x*exp(x)*sin(x), x*cos(x), x*cos(x)*exp(x),
        x*sin(x), cos(x)*exp(x), x**2*exp(x)*sin(x)])}
    assert _undetermined_coefficients_match(4*x*sin(x - 2), x) == {
        'trialset': set([x*cos(x - 2), x*sin(x - 2), cos(x - 2), sin(x - 2)]),
        'test': True,
    }
    assert _undetermined_coefficients_match(2**x*x, x) == \
        {'test': True, 'trialset': set([2**x, x*2**x])}
    assert _undetermined_coefficients_match(2**x*exp(2*x), x) == \
        {'test': True, 'trialset': set([2**x*exp(2*x)])}
    assert _undetermined_coefficients_match(exp(-x)/x, x) == \
        {'test': False}
    # Below are from Ordinary Differential Equations,
    #                Tenenbaum and Pollard, pg. 231
    assert _undetermined_coefficients_match(S(4), x) == \
        {'test': True, 'trialset': set([S.One])}
    assert _undetermined_coefficients_match(12*exp(x), x) == \
        {'test': True, 'trialset': set([exp(x)])}
    assert _undetermined_coefficients_match(exp(I*x), x) == \
        {'test': True, 'trialset': set([exp(I*x)])}
    assert _undetermined_coefficients_match(sin(x), x) == \
        {'test': True, 'trialset': set([cos(x), sin(x)])}
    assert _undetermined_coefficients_match(cos(x), x) == \
        {'test': True, 'trialset': set([cos(x), sin(x)])}
    assert _undetermined_coefficients_match(8 + 6*exp(x) + 2*sin(x), x) == \
        {'test': True, 'trialset': set([S.One, cos(x), sin(x), exp(x)])}
    assert _undetermined_coefficients_match(x**2, x) == \
        {'test': True, 'trialset': set([S.One, x, x**2])}
    assert _undetermined_coefficients_match(9*x*exp(x) + exp(-x), x) == \
        {'test': True, 'trialset': set([x*exp(x), exp(x), exp(-x)])}
    assert _undetermined_coefficients_match(2*exp(2*x)*sin(x), x) == \
        {'test': True, 'trialset': set([exp(2*x)*sin(x), cos(x)*exp(2*x)])}
    assert _undetermined_coefficients_match(x - sin(x), x) == \
        {'test': True, 'trialset': set([S.One, x, cos(x), sin(x)])}
    assert _undetermined_coefficients_match(x**2 + 2*x, x) == \
        {'test': True, 'trialset': set([S.One, x, x**2])}
    assert _undetermined_coefficients_match(4*x*sin(x), x) == \
        {'test': True, 'trialset': set([x*cos(x), x*sin(x), cos(x), sin(x)])}
    assert _undetermined_coefficients_match(x*sin(2*x), x) == \
        {'test': True, 'trialset':
            set([x*cos(2*x), x*sin(2*x), cos(2*x), sin(2*x)])}
    assert _undetermined_coefficients_match(x**2*exp(-x), x) == \
        {'test': True, 'trialset': set([x*exp(-x), x**2*exp(-x), exp(-x)])}
    assert _undetermined_coefficients_match(2*exp(-x) - x**2*exp(-x), x) == \
        {'test': True, 'trialset': set([x*exp(-x), x**2*exp(-x), exp(-x)])}
    assert _undetermined_coefficients_match(exp(-2*x) + x**2, x) == \
        {'test': True, 'trialset': set([S.One, x, x**2, exp(-2*x)])}
    assert _undetermined_coefficients_match(x*exp(-x), x) == \
        {'test': True, 'trialset': set([x*exp(-x), exp(-x)])}
    assert _undetermined_coefficients_match(x + exp(2*x), x) == \
        {'test': True, 'trialset': set([S.One, x, exp(2*x)])}
    assert _undetermined_coefficients_match(sin(x) + exp(-x), x) == \
        {'test': True, 'trialset': set([cos(x), sin(x), exp(-x)])}
    assert _undetermined_coefficients_match(exp(x), x) == \
        {'test': True, 'trialset': set([exp(x)])}
    # converted from sin(x)**2
    assert _undetermined_coefficients_match(S.Half - cos(2*x)/2, x) == \
        {'test': True, 'trialset': set([S.One, cos(2*x), sin(2*x)])}
    # converted from exp(2*x)*sin(x)**2
    assert _undetermined_coefficients_match(
        exp(2*x)*(S.Half + cos(2*x)/2), x
    ) == {
        'test': True, 'trialset': set([exp(2*x)*sin(2*x), cos(2*x)*exp(2*x),
        exp(2*x)])}
    assert _undetermined_coefficients_match(2*x + sin(x) + cos(x), x) == \
        {'test': True, 'trialset': set([S.One, x, cos(x), sin(x)])}
    # converted from sin(2*x)*sin(x)
    assert _undetermined_coefficients_match(cos(x)/2 - cos(3*x)/2, x) == \
        {'test': True, 'trialset': set([cos(x), cos(3*x), sin(x), sin(3*x)])}
    assert _undetermined_coefficients_match(cos(x**2), x) == {'test': False}
    assert _undetermined_coefficients_match(2**(x**2), x) == {'test': False}


@slow
def test_nth_linear_constant_coeff_undetermined_coefficients():
    hint = 'nth_linear_constant_coeff_undetermined_coefficients'
    g = exp(-x)
    f2 = f(x).diff(x, 2)
    c = 3*f(x).diff(x, 3) + 5*f2 + f(x).diff(x) - f(x) - x
    eq1 = c - x*g
    eq2 = c - g
    # 3-27 below are from Ordinary Differential Equations,
    #                     Tenenbaum and Pollard, pg. 231
    eq3 = f2 + 3*f(x).diff(x) + 2*f(x) - 4
    eq4 = f2 + 3*f(x).diff(x) + 2*f(x) - 12*exp(x)
    eq5 = f2 + 3*f(x).diff(x) + 2*f(x) - exp(I*x)
    eq6 = f2 + 3*f(x).diff(x) + 2*f(x) - sin(x)
    eq7 = f2 + 3*f(x).diff(x) + 2*f(x) - cos(x)
    eq8 = f2 + 3*f(x).diff(x) + 2*f(x) - (8 + 6*exp(x) + 2*sin(x))
    eq9 = f2 + f(x).diff(x) + f(x) - x**2
    eq10 = f2 - 2*f(x).diff(x) - 8*f(x) - 9*x*exp(x) - 10*exp(-x)
    eq11 = f2 - 3*f(x).diff(x) - 2*exp(2*x)*sin(x)
    eq12 = f(x).diff(x, 4) - 2*f2 + f(x) - x + sin(x)
    eq13 = f2 + f(x).diff(x) - x**2 - 2*x
    eq14 = f2 + f(x).diff(x) - x - sin(2*x)
    eq15 = f2 + f(x) - 4*x*sin(x)
    eq16 = f2 + 4*f(x) - x*sin(2*x)
    eq17 = f2 + 2*f(x).diff(x) + f(x) - x**2*exp(-x)
    eq18 = f(x).diff(x, 3) + 3*f2 + 3*f(x).diff(x) + f(x) - 2*exp(-x) + \
        x**2*exp(-x)
    eq19 = f2 + 3*f(x).diff(x) + 2*f(x) - exp(-2*x) - x**2
    eq20 = f2 - 3*f(x).diff(x) + 2*f(x) - x*exp(-x)
    eq21 = f2 + f(x).diff(x) - 6*f(x) - x - exp(2*x)
    eq22 = f2 + f(x) - sin(x) - exp(-x)
    eq23 = f(x).diff(x, 3) - 3*f2 + 3*f(x).diff(x) - f(x) - exp(x)
    # sin(x)**2
    eq24 = f2 + f(x) - S.Half - cos(2*x)/2
    # exp(2*x)*sin(x)**2
    eq25 = f(x).diff(x, 3) - f(x).diff(x) - exp(2*x)*(S.Half - cos(2*x)/2)
    eq26 = (f(x).diff(x, 5) + 2*f(x).diff(x, 3) + f(x).diff(x) - 2*x -
        sin(x) - cos(x))
    # sin(2*x)*sin(x), skip 3127 for now, match bug
    eq27 = f2 + f(x) - cos(x)/2 + cos(3*x)/2
    eq28 = f(x).diff(x) - 1
    sol1 = Eq(f(x), C3*exp(x/3) - x + (C1 + x*(C2 - x**2/24 - 3*x/32))*exp(-x) - 1)
    sol2 = Eq(f(x), C3*exp(x/3) - x + (C1 + x*(C2 - x/8))*exp(-x) - 1)
    sol3 = Eq(f(x), 2 + C1*exp(-x) + C2*exp(-2*x))
    sol4 = Eq(f(x), 2*exp(x) + C1*exp(-x) + C2*exp(-2*x))
    sol5 = Eq(f(x), C1*exp(-2*x) + C2*exp(-x) + exp(I*x)/10 - 3*I*exp(I*x)/10)
    sol6 = Eq(f(x), -3*cos(x)/10 + sin(x)/10 + C1*exp(-x) + C2*exp(-2*x))
    sol7 = Eq(f(x), cos(x)/10 + 3*sin(x)/10 + C1*exp(-x) + C2*exp(-2*x))
    sol8 = Eq(f(x),
        4 - 3*cos(x)/5 + sin(x)/5 + exp(x) + C1*exp(-x) + C2*exp(-2*x))
    sol9 = Eq(f(x),
        -2*x + x**2 + (C1*sin(x*sqrt(3)/2) + C2*cos(x*sqrt(3)/2))*exp(-x/2))
    sol10 = Eq(f(x), -x*exp(x) - 2*exp(-x) + C1*exp(-2*x) + C2*exp(4*x))
    sol11 = Eq(f(x), C1 + C2*exp(3*x) + (-3*sin(x) - cos(x))*exp(2*x)/5)
    sol12 = Eq(f(x), x - sin(x)/4 + (C1 + C2*x)*exp(-x) + (C3 + C4*x)*exp(x))
    sol13 = Eq(f(x), C1 + x**3/3 + C2*exp(-x))
    sol14 = Eq(f(x), C1 - x - sin(2*x)/5 - cos(2*x)/10 + x**2/2 + C2*exp(-x))
    sol15 = Eq(f(x), (C1 + x)*sin(x) + (C2 - x**2)*cos(x))
    sol16 = Eq(f(x), (C1 + x/16)*sin(2*x) + (C2 - x**2/8)*cos(2*x))
    sol17 = Eq(f(x), (C1 + x*(C2 + x**3/12))*exp(-x))
    sol18 = Eq(f(x), (C1 + x*(C2 + x*(C3 - x**3/60 + x/3)))*exp(-x))
    sol19 = Eq(f(x), Rational(7, 4) - x*Rational(3, 2) + x**2/2 + C1*exp(-x) + (C2 - x)*exp(-2*x))
    sol20 = Eq(f(x), C1*exp(x) + C2*exp(2*x) + (6*x + 5)*exp(-x)/36)
    sol21 = Eq(f(x), Rational(-1, 36) - x/6 + C1*exp(-3*x) + (C2 + x/5)*exp(2*x))
    sol22 = Eq(f(x), C1*sin(x) + (C2 - x/2)*cos(x) + exp(-x)/2)
    sol23 = Eq(f(x), (C1 + x*(C2 + x*(C3 + x/6)))*exp(x))
    sol24 = Eq(f(x), S.Half - cos(2*x)/6 + C1*sin(x) + C2*cos(x))
    sol25 = Eq(f(x), C1 + C2*exp(-x) + C3*exp(x) +
               (-21*sin(2*x) + 27*cos(2*x) + 130)*exp(2*x)/1560)
    sol26 = Eq(f(x), C1 + x**2 + (C2 + x*(C3 - x/8))*sin(x) + (C4 + x*(C5 + x/8))*cos(x))
    sol27 = Eq(f(x), cos(3*x)/16 + C1*cos(x) + (C2 + x/4)*sin(x))
    sol28 = Eq(f(x), C1 + x)
    sol1s = constant_renumber(sol1)
    sol2s = constant_renumber(sol2)
    sol3s = constant_renumber(sol3)
    sol4s = constant_renumber(sol4)
    sol5s = constant_renumber(sol5)
    sol6s = constant_renumber(sol6)
    sol7s = constant_renumber(sol7)
    sol8s = constant_renumber(sol8)
    sol9s = constant_renumber(sol9)
    sol10s = constant_renumber(sol10)
    sol11s = constant_renumber(sol11)
    sol12s = constant_renumber(sol12)
    sol13s = constant_renumber(sol13)
    sol14s = constant_renumber(sol14)
    sol15s = constant_renumber(sol15)
    sol16s = constant_renumber(sol16)
    sol17s = constant_renumber(sol17)
    sol18s = constant_renumber(sol18)
    sol19s = constant_renumber(sol19)
    sol20s = constant_renumber(sol20)
    sol21s = constant_renumber(sol21)
    sol22s = constant_renumber(sol22)
    sol23s = constant_renumber(sol23)
    sol24s = constant_renumber(sol24)
    sol25s = constant_renumber(sol25)
    sol26s = constant_renumber(sol26)
    sol27s = constant_renumber(sol27)
    got = dsolve(eq1, hint=hint)
    assert got in (sol1, sol1s), got
    got = dsolve(eq2, hint=hint)
    assert got in (sol2, sol2s), got
    assert dsolve(eq3, hint=hint) in (sol3, sol3s)
    assert dsolve(eq4, hint=hint) in (sol4, sol4s)
    assert dsolve(eq5, hint=hint) in (sol5, sol5s)
    assert dsolve(eq6, hint=hint) in (sol6, sol6s)
    assert dsolve(eq7, hint=hint) in (sol7, sol7s)
    assert dsolve(eq8, hint=hint) in (sol8, sol8s)
    assert dsolve(eq9, hint=hint) in (sol9, sol9s)
    assert dsolve(eq10, hint=hint) in (sol10, sol10s)
    assert dsolve(eq11, hint=hint) in (sol11, sol11s)
    assert dsolve(eq12, hint=hint) in (sol12, sol12s)
    assert dsolve(eq13, hint=hint) in (sol13, sol13s)
    assert dsolve(eq14, hint=hint) in (sol14, sol14s)
    assert dsolve(eq15, hint=hint) in (sol15, sol15s)
    assert dsolve(eq16, hint=hint) in (sol16, sol16s)
    got = dsolve(eq17, hint=hint)
    assert got in (sol17, sol17s), got
    got = dsolve(eq18, hint=hint)
    assert got in (sol18, sol18s), got
    assert dsolve(eq19, hint=hint) in (sol19, sol19s)
    assert dsolve(eq20, hint=hint) in (sol20, sol20s)
    assert dsolve(eq21, hint=hint) in (sol21, sol21s)
    assert dsolve(eq22, hint=hint) in (sol22, sol22s)
    got = dsolve(eq23, hint=hint)
    assert got in (sol23, sol23s), got
    assert dsolve(eq24, hint=hint) in (sol24, sol24s)
    got = dsolve(eq25, hint=hint)
    assert got in (sol25, sol25s), got
    got = dsolve(eq26, hint=hint)
    assert got in (sol26, sol26s), got
    assert dsolve(eq27, hint=hint) in (sol27, sol27s)
    assert dsolve(eq28, hint=hint) == sol28
    assert checkodesol(eq1, sol1, order=3, solve_for_func=False)[0]
    assert checkodesol(eq2, sol2, order=3, solve_for_func=False)[0]
    assert checkodesol(eq3, sol3, order=2, solve_for_func=False)[0]
    assert checkodesol(eq4, sol4, order=2, solve_for_func=False)[0]
    assert checkodesol(eq5, sol5, order=2, solve_for_func=False)[0]
    assert checkodesol(eq6, sol6, order=2, solve_for_func=False)[0]
    assert checkodesol(eq7, sol7, order=2, solve_for_func=False)[0]
    assert checkodesol(eq8, sol8, order=2, solve_for_func=False)[0]
    assert checkodesol(eq9, sol9, order=2, solve_for_func=False)[0]
    assert checkodesol(eq10, sol10, order=2, solve_for_func=False)[0]
    assert checkodesol(eq11, sol11, order=2, solve_for_func=False)[0]
    assert checkodesol(eq12, sol12, order=4, solve_for_func=False)[0]
    assert checkodesol(eq13, sol13, order=2, solve_for_func=False)[0]
    assert checkodesol(eq14, sol14, order=2, solve_for_func=False)[0]
    assert checkodesol(eq15, sol15, order=2, solve_for_func=False)[0]
    assert checkodesol(eq16, sol16, order=2, solve_for_func=False)[0]
    assert checkodesol(eq17, sol17, order=2, solve_for_func=False)[0]
    assert checkodesol(eq18, sol18, order=3, solve_for_func=False)[0]
    assert checkodesol(eq19, sol19, order=2, solve_for_func=False)[0]
    assert checkodesol(eq20, sol20, order=2, solve_for_func=False)[0]
    assert checkodesol(eq21, sol21, order=2, solve_for_func=False)[0]
    assert checkodesol(eq22, sol22, order=2, solve_for_func=False)[0]
    assert checkodesol(eq23, sol23, order=3, solve_for_func=False)[0]
    assert checkodesol(eq24, sol24, order=2, solve_for_func=False)[0]
    assert checkodesol(eq25, sol25, order=3, solve_for_func=False)[0]
    assert checkodesol(eq26, sol26, order=5, solve_for_func=False)[0]
    assert checkodesol(eq27, sol27, order=2, solve_for_func=False)[0]
    assert checkodesol(eq28, sol28, order=1, solve_for_func=False)[0]


def test_issue_12623():
    t = symbols("t")
    u = symbols("u",cls=Function)
    R, L, C, E_0, alpha = symbols("R L C E_0 alpha",positive=True)
    omega = Symbol('omega')

    eqRLC_1 = Eq( u(t).diff(t,t) + R /L*u(t).diff(t) + 1/(L*C)*u(t), alpha)
    sol_1 = Eq(u(t), C*L*alpha + C1*exp(t*(-R - sqrt(C*R**2 - 4*L)/sqrt(C))/(2*L)) + C2*exp(t*(-R + sqrt(C*R**2 - 4*L)/sqrt(C))/(2*L)))
    assert dsolve(eqRLC_1) == sol_1
    assert checkodesol(eqRLC_1, sol_1) == (True, 0)

    eqRLC_2 = Eq( L*C*u(t).diff(t,t) + R*C*u(t).diff(t) + u(t), E_0*exp(I*omega*t) )
    sol_2 = Eq(u(t), C1*exp(t*(-R - sqrt(C*R**2 - 4*L)/sqrt(C))/(2*L)) + C2*exp(t*(-R + sqrt(C*R**2 - 4*L)/sqrt(C))/(2*L)) + E_0*exp(I*omega*t)/(-C*L*omega**2 + I*C*R*omega + 1))
    assert dsolve(eqRLC_2) == sol_2
    assert checkodesol(eqRLC_2, sol_2) == (True, 0)
    #issue-https://github.com/sympy/sympy/issues/12623


def test_issue_5787():
    # This test case is to show the classification of imaginary constants under
    # nth_linear_constant_coeff_undetermined_coefficients
    eq = Eq(diff(f(x), x), I*f(x) + S.Half - I)
    our_hint = 'nth_linear_constant_coeff_undetermined_coefficients'
    assert our_hint in classify_ode(eq)


def test_nth_linear_constant_coeff_undetermined_coefficients_imaginary_exp():
    # Equivalent to eq26 in
    # test_nth_linear_constant_coeff_undetermined_coefficients above. This
    # previously failed because the algorithm for undetermined coefficients
    # didn't know to multiply exp(I*x) by sufficient x because it is linearly
    # dependent on sin(x) and cos(x).
    hint = 'nth_linear_constant_coeff_undetermined_coefficients'
    eq26a = f(x).diff(x, 5) + 2*f(x).diff(x, 3) + f(x).diff(x) - 2*x - exp(I*x)
    sol26 = Eq(f(x), C1 + x**2*(I*exp(I*x)/8 + 1) + (C2 + C3*x)*sin(x) + (C4 + C5*x)*cos(x))
    assert dsolve(eq26a, hint=hint) == sol26
    assert checkodesol(eq26a, sol26) == (True, 0)


@slow
def test_nth_linear_constant_coeff_variation_of_parameters():
    hint = 'nth_linear_constant_coeff_variation_of_parameters'
    g = exp(-x)
    f2 = f(x).diff(x, 2)
    c = 3*f(x).diff(x, 3) + 5*f2 + f(x).diff(x) - f(x) - x
    eq1 = c - x*g
    eq2 = c - g
    eq3 = f(x).diff(x) - 1
    eq4 = f2 + 3*f(x).diff(x) + 2*f(x) - 4
    eq5 = f2 + 3*f(x).diff(x) + 2*f(x) - 12*exp(x)
    eq6 = f2 - 2*f(x).diff(x) - 8*f(x) - 9*x*exp(x) - 10*exp(-x)
    eq7 = f2 + 2*f(x).diff(x) + f(x) - x**2*exp(-x)
    eq8 = f2 - 3*f(x).diff(x) + 2*f(x) - x*exp(-x)
    eq9 = f(x).diff(x, 3) - 3*f2 + 3*f(x).diff(x) - f(x) - exp(x)
    eq10 = f2 + 2*f(x).diff(x) + f(x) - exp(-x)/x
    eq11 = f2 + f(x) - 1/sin(x)*1/cos(x)
    eq12 = f(x).diff(x, 4) - 1/x
    sol1 = Eq(f(x), C3*exp(x/3) - x + (C1 + x*(C2 - x**2/24 - 3*x/32))*exp(-x) - 1)
    sol2 = Eq(f(x), C3*exp(x/3) - x + (C1 + x*(C2 - x/8))*exp(-x) - 1)
    sol3 = Eq(f(x), C1 + x)
    sol4 = Eq(f(x), 2 + C1*exp(-x) + C2*exp(-2*x))
    sol5 = Eq(f(x), 2*exp(x) + C1*exp(-x) + C2*exp(-2*x))
    sol6 = Eq(f(x), -x*exp(x) - 2*exp(-x) + C1*exp(-2*x) + C2*exp(4*x))
    sol7 = Eq(f(x), (C1 + x*(C2 + x**3/12))*exp(-x))
    sol8 = Eq(f(x), C1*exp(x) + C2*exp(2*x) + (6*x + 5)*exp(-x)/36)
    sol9 = Eq(f(x), (C1 + x*(C2 + x*(C3 + x/6)))*exp(x))
    sol10 = Eq(f(x), (C1 + x*(C2 + log(x)))*exp(-x))
    sol11 = Eq(f(x), (C1 + log(sin(x) - 1)/2 - log(sin(x) + 1)/2
        )*cos(x) + (C2 + log(cos(x) - 1)/2 - log(cos(x) + 1)/2)*sin(x))
    sol12 = Eq(f(x), C1 + C2*x + x**3*(C3 + log(x)/6) + C4*x**2)
    sol1s = constant_renumber(sol1)
    sol2s = constant_renumber(sol2)
    sol3s = constant_renumber(sol3)
    sol4s = constant_renumber(sol4)
    sol5s = constant_renumber(sol5)
    sol6s = constant_renumber(sol6)
    sol7s = constant_renumber(sol7)
    sol8s = constant_renumber(sol8)
    sol9s = constant_renumber(sol9)
    sol10s = constant_renumber(sol10)
    sol11s = constant_renumber(sol11)
    sol12s = constant_renumber(sol12)
    got = dsolve(eq1, hint=hint)
    assert got in (sol1, sol1s), got
    got = dsolve(eq2, hint=hint)
    assert got in (sol2, sol2s), got
    assert dsolve(eq3, hint=hint) in (sol3, sol3s)
    assert dsolve(eq4, hint=hint) in (sol4, sol4s)
    assert dsolve(eq5, hint=hint) in (sol5, sol5s)
    assert dsolve(eq6, hint=hint) in (sol6, sol6s)
    got = dsolve(eq7, hint=hint)
    assert got in (sol7, sol7s), got
    assert dsolve(eq8, hint=hint) in (sol8, sol8s)
    got = dsolve(eq9, hint=hint)
    assert got in (sol9, sol9s), got
    assert dsolve(eq10, hint=hint) in (sol10, sol10s)
    assert dsolve(eq11, hint=hint + '_Integral').doit() in (sol11, sol11s)
    assert dsolve(eq12, hint=hint) in (sol12, sol12s)
    assert checkodesol(eq1, sol1, order=3, solve_for_func=False)[0]
    assert checkodesol(eq2, sol2, order=3, solve_for_func=False)[0]
    assert checkodesol(eq3, sol3, order=1, solve_for_func=False)[0]
    assert checkodesol(eq4, sol4, order=2, solve_for_func=False)[0]
    assert checkodesol(eq5, sol5, order=2, solve_for_func=False)[0]
    assert checkodesol(eq6, sol6, order=2, solve_for_func=False)[0]
    assert checkodesol(eq7, sol7, order=2, solve_for_func=False)[0]
    assert checkodesol(eq8, sol8, order=2, solve_for_func=False)[0]
    assert checkodesol(eq9, sol9, order=3, solve_for_func=False)[0]
    assert checkodesol(eq10, sol10, order=2, solve_for_func=False)[0]
    assert checkodesol(eq12, sol12, order=4, solve_for_func=False)[0]


@slow
def test_nth_linear_constant_coeff_variation_of_parameters_simplify_False():
    # solve_variation_of_parameters shouldn't attempt to simplify the
    # Wronskian if simplify=False.  If wronskian() ever gets good enough
    # to simplify the result itself, this test might fail.
    our_hint = 'nth_linear_constant_coeff_variation_of_parameters_Integral'
    eq = f(x).diff(x, 5) + 2*f(x).diff(x, 3) + f(x).diff(x) - 2*x - exp(I*x)
    sol_simp = dsolve(eq, f(x), hint=our_hint, simplify=True)
    sol_nsimp = dsolve(eq, f(x), hint=our_hint, simplify=False)
    assert sol_simp != sol_nsimp
    assert checkodesol(eq, sol_simp, order=5, solve_for_func=False) == (True, 0)
    assert checkodesol(eq, sol_simp, order=5, solve_for_func=False) == (True, 0)


def test_unexpanded_Liouville_ODE():
    # This is the same as eq1 from test_Liouville_ODE() above.
    eq1 = diff(f(x), x)/x + diff(f(x), x, x)/2 - diff(f(x), x)**2/2
    eq2 = eq1*exp(-f(x))/exp(f(x))
    sol2 = Eq(f(x), C1 + log(x) - log(C2 + x))
    sol2s = constant_renumber(sol2)
    assert dsolve(eq2) in (sol2, sol2s)
    assert checkodesol(eq2, sol2, order=2, solve_for_func=False)[0]

def test_issue_4785():
    from sympy.abc import A
    eq = x + A*(x + diff(f(x), x) + f(x)) + diff(f(x), x) + f(x) + 2
    assert classify_ode(eq, f(x)) == ('1st_linear', 'almost_linear',
        '1st_power_series', 'lie_group',
        'nth_linear_constant_coeff_undetermined_coefficients',
        'nth_linear_constant_coeff_variation_of_parameters',
        '1st_linear_Integral', 'almost_linear_Integral',
        'nth_linear_constant_coeff_variation_of_parameters_Integral')
    # issue 4864
    eq = (x**2 + f(x)**2)*f(x).diff(x) - 2*x*f(x)
    assert classify_ode(eq, f(x)) == ('1st_exact',
        '1st_homogeneous_coeff_best',
        '1st_homogeneous_coeff_subs_indep_div_dep',
        '1st_homogeneous_coeff_subs_dep_div_indep',
        '1st_power_series',
        'lie_group', '1st_exact_Integral',
        '1st_homogeneous_coeff_subs_indep_div_dep_Integral',
        '1st_homogeneous_coeff_subs_dep_div_indep_Integral')

def test_issue_4825():
    raises(ValueError, lambda: dsolve(f(x, y).diff(x) - y*f(x, y), f(x)))
    assert classify_ode(f(x, y).diff(x) - y*f(x, y), f(x), dict=True) == \
        {'order': 0, 'default': None, 'ordered_hints': ()}
    # See also issue 3793, test Z13.
    raises(ValueError, lambda: dsolve(f(x).diff(x), f(y)))
    assert classify_ode(f(x).diff(x), f(y), dict=True) == \
        {'order': 0, 'default': None, 'ordered_hints': ()}


def test_constant_renumber_order_issue_5308():
    from sympy.utilities.iterables import variations

    assert constant_renumber(C1*x + C2*y) == \
        constant_renumber(C1*y + C2*x) == \
        C1*x + C2*y
    e = C1*(C2 + x)*(C3 + y)
    for a, b, c in variations([C1, C2, C3], 3):
        assert constant_renumber(a*(b + x)*(c + y)) == e


def test_issue_5770():
    k = Symbol("k", real=True)
    t = Symbol('t')
    w = Function('w')
    sol = dsolve(w(t).diff(t, 6) - k**6*w(t), w(t))
    assert len([s for s in sol.free_symbols if s.name.startswith('C')]) == 6
    assert constantsimp((C1*cos(x) + C2*cos(x))*exp(x), set([C1, C2])) == \
        C1*cos(x)*exp(x)
    assert constantsimp(C1*cos(x) + C2*cos(x) + C3*sin(x), set([C1, C2, C3])) == \
        C1*cos(x) + C3*sin(x)
    assert constantsimp(exp(C1 + x), set([C1])) == C1*exp(x)
    assert constantsimp(x + C1 + y, set([C1, y])) == C1 + x
    assert constantsimp(x + C1 + Integral(x, (x, 1, 2)), set([C1])) == C1 + x


def test_issue_5112_5430():
    assert homogeneous_order(-log(x) + acosh(x), x) is None
    assert homogeneous_order(y - log(x), x, y) is None


def test_issue_5095():
    f = Function('f')
    raises(ValueError, lambda: dsolve(f(x).diff(x)**2, f(x), 'fdsjf'))


def test_exact_enhancement():
    f = Function('f')(x)
    df = Derivative(f, x)
    eq = f/x**2 + ((f*x - 1)/x)*df
    sol = [Eq(f, (i*sqrt(C1*x**2 + 1) + 1)/x) for i in (-1, 1)]
    assert set(dsolve(eq, f)) == set(sol)
    assert checkodesol(eq, sol, order=1, solve_for_func=False) == [(True, 0), (True, 0)]

    eq = (x*f - 1) + df*(x**2 - x*f)
    sol = [Eq(f, x - sqrt(C1 + x**2 - 2*log(x))),
           Eq(f, x + sqrt(C1 + x**2 - 2*log(x)))]
    assert set(dsolve(eq, f)) == set(sol)
    assert checkodesol(eq, sol, order=1, solve_for_func=False) == [(True, 0), (True, 0)]

    eq = (x + 2)*sin(f) + df*x*cos(f)
    sol = [Eq(f, -asin(C1*exp(-x)/x**2) + pi),
           Eq(f, asin(C1*exp(-x)/x**2))]
    assert set(dsolve(eq, f)) == set(sol)
    assert checkodesol(eq, sol, order=1, solve_for_func=False) == [(True, 0), (True, 0)]


@slow
def test_separable_reduced():
    f = Function('f')
    x = Symbol('x')
    df = f(x).diff(x)
    eq = (x / f(x))*df  + tan(x**2*f(x) / (x**2*f(x) - 1))
    assert classify_ode(eq) == ('separable_reduced', 'lie_group',
        'separable_reduced_Integral')

    eq = x* df  + f(x)* (1 / (x**2*f(x) - 1))
    assert classify_ode(eq) == ('separable_reduced', 'lie_group',
        'separable_reduced_Integral')
    sol = dsolve(eq, hint = 'separable_reduced', simplify=False)
    assert sol.lhs ==  log(x**2*f(x))/3 + log(x**2*f(x) - Rational(3, 2))/6
    assert sol.rhs == C1 + log(x)
    assert checkodesol(eq, sol, order=1, solve_for_func=False)[0]

    eq = f(x).diff(x) + (f(x) / (x**4*f(x) - x))
    assert classify_ode(eq) == ('separable_reduced', 'lie_group',
        'separable_reduced_Integral')
    sol = dsolve(eq, hint = 'separable_reduced')
    # FIXME: This one hangs
    #assert checkodesol(eq, sol, order=1, solve_for_func=False) == [(True, 0)] * 4
    assert len(sol) == 4

    eq = x*df + f(x)*(x**2*f(x))
    sol = dsolve(eq, hint = 'separable_reduced', simplify=False)
    assert sol == Eq(log(x**2*f(x))/2 - log(x**2*f(x) - 2)/2, C1 + log(x))
    assert checkodesol(eq, sol, order=1, solve_for_func=False)[0]

    eq = Eq(f(x).diff(x) + f(x)/x * (1 + (x**(S(2)/3)*f(x))**2), 0)
    sol = dsolve(eq, hint = 'separable_reduced', simplify=False)
    assert sol == Eq(-3*log(x**(S(2)/3)*f(x)) + 3*log(3*x**(S(4)/3)*f(x)**2 + 1)/2, C1 + log(x))
    assert checkodesol(eq, sol, solve_for_func=False) == (True, 0)

    eq = Eq(f(x).diff(x) + f(x)/x * (1 + (x*f(x))**2), 0)
    sol = dsolve(eq, hint = 'separable_reduced')
    assert sol == [Eq(f(x), -sqrt(2)*sqrt(1/(C1 + log(x)))/(2*x)),\
                   Eq(f(x), sqrt(2)*sqrt(1/(C1 + log(x)))/(2*x))]
    assert checkodesol(eq, sol) == [(True, 0)]*2

    eq = Eq(f(x).diff(x) + (x**4*f(x)**2 + x**2*f(x))*f(x)/(x*(x**6*f(x)**3 + x**4*f(x)**2)), 0)
    sol = dsolve(eq, hint = 'separable_reduced')
    assert sol == Eq(f(x), C1 + 1/(2*x**2))
    assert checkodesol(eq, sol) == (True, 0)

    eq = Eq(f(x).diff(x) + (f(x)**2)*f(x)/(x), 0)
    sol = dsolve(eq, hint = 'separable_reduced')
    assert sol == [Eq(f(x), -sqrt(2)*sqrt(1/(C1 + log(x)))/2),\
                  Eq(f(x), sqrt(2)*sqrt(1/(C1 + log(x)))/2)]
    assert checkodesol(eq, sol) == [(True, 0), (True, 0)]

    eq = Eq(f(x).diff(x) + (f(x)+3)*f(x)/(x*(f(x)+2)), 0)
    sol = dsolve(eq, hint = 'separable_reduced', simplify=False)
    assert sol == Eq(-log(f(x) + 3)/3 - 2*log(f(x))/3, C1 + log(x))
    assert checkodesol(eq, sol, solve_for_func=False) == (True, 0)

    eq = Eq(f(x).diff(x) + (f(x)+3)*f(x)/x, 0)
    sol = dsolve(eq, hint = 'separable_reduced')
    assert sol == Eq(f(x), 3/(C1*x**3 - 1))
    assert checkodesol(eq, sol) == (True, 0)

    eq = Eq(f(x).diff(x) + (f(x)**2+f(x))*f(x)/(x), 0)
    sol = dsolve(eq, hint='separable_reduced', simplify=False)
    assert sol == Eq(-log(f(x) + 1) + log(f(x)) + 1/f(x), C1 + log(x))
    assert checkodesol(eq, sol, solve_for_func=False) == (True, 0)


def test_homogeneous_function():
    f = Function('f')
    eq1 = tan(x + f(x))
    eq2 = sin((3*x)/(4*f(x)))
    eq3 = cos(x*f(x)*Rational(3, 4))
    eq4 = log((3*x + 4*f(x))/(5*f(x) + 7*x))
    eq5 = exp((2*x**2)/(3*f(x)**2))
    eq6 = log((3*x + 4*f(x))/(5*f(x) + 7*x) + exp((2*x**2)/(3*f(x)**2)))
    eq7 = sin((3*x)/(5*f(x) + x**2))
    assert homogeneous_order(eq1, x, f(x)) == None
    assert homogeneous_order(eq2, x, f(x)) == 0
    assert homogeneous_order(eq3, x, f(x)) == None
    assert homogeneous_order(eq4, x, f(x)) == 0
    assert homogeneous_order(eq5, x, f(x)) == 0
    assert homogeneous_order(eq6, x, f(x)) == 0
    assert homogeneous_order(eq7, x, f(x)) == None


def test_linear_coeff_match():
    n, d = z*(2*x + 3*f(x) + 5), z*(7*x + 9*f(x) + 11)
    rat = n/d
    eq1 = sin(rat) + cos(rat.expand())
    eq2 = rat
    eq3 = log(sin(rat))
    ans = (4, Rational(-13, 3))
    assert _linear_coeff_match(eq1, f(x)) == ans
    assert _linear_coeff_match(eq2, f(x)) == ans
    assert _linear_coeff_match(eq3, f(x)) == ans

    # no c
    eq4 = (3*x)/f(x)
    # not x and f(x)
    eq5 = (3*x + 2)/x
    # denom will be zero
    eq6 = (3*x + 2*f(x) + 1)/(3*x + 2*f(x) + 5)
    # not rational coefficient
    eq7 = (3*x + 2*f(x) + sqrt(2))/(3*x + 2*f(x) + 5)
    assert _linear_coeff_match(eq4, f(x)) is None
    assert _linear_coeff_match(eq5, f(x)) is None
    assert _linear_coeff_match(eq6, f(x)) is None
    assert _linear_coeff_match(eq7, f(x)) is None


def test_linear_coefficients():
    f = Function('f')
    sol = Eq(f(x), C1/(x**2 + 6*x + 9) - Rational(3, 2))
    eq = f(x).diff(x) + (3 + 2*f(x))/(x + 3)
    assert dsolve(eq, hint='linear_coefficients') == sol
    assert checkodesol(eq, sol, order=1, solve_for_func=False)[0]


def test_constantsimp_take_problem():
    c = exp(C1) + 2
    assert len(Poly(constantsimp(exp(C1) + c + c*x, [C1])).gens) == 2


def test_issue_6879():
    f = Function('f')
    eq = Eq(Derivative(f(x), x, 2) - 2*Derivative(f(x), x) + f(x), sin(x))
    sol = (C1 + C2*x)*exp(x) + cos(x)/2
    assert dsolve(eq).rhs == sol
    assert checkodesol(eq, sol, order=1, solve_for_func=False)[0]


def test_issue_6989():
    f = Function('f')
    k = Symbol('k')

    eq = f(x).diff(x) - x*exp(-k*x)
    csol = Eq(f(x), C1 + Piecewise(
            ((-k*x - 1)*exp(-k*x)/k**2, Ne(k**2, 0)),
            (x**2/2, True)
        ))
    sol = dsolve(eq, f(x))
    assert sol == csol
    assert checkodesol(eq, sol, order=1, solve_for_func=False)[0]

    eq = -f(x).diff(x) + x*exp(-k*x)
    csol = Eq(f(x), C1 + Piecewise(
        ((-k*x - 1)*exp(-k*x)/k**2, Ne(k**2, 0)),
        (x**2/2, True)
    ))
    sol = dsolve(eq, f(x))
    assert sol == csol
    assert checkodesol(eq, sol, order=1, solve_for_func=False)[0]


def test_heuristic1():
    y, a, b, c, a4, a3, a2, a1, a0 = symbols("y a b c a4 a3 a2 a1 a0")
    f = Function('f')
    xi = Function('xi')
    eta = Function('eta')
    df = f(x).diff(x)
    eq = Eq(df, x**2*f(x))
    eq1 = f(x).diff(x) + a*f(x) - c*exp(b*x)
    eq2 = f(x).diff(x) + 2*x*f(x) - x*exp(-x**2)
    eq3 = (1 + 2*x)*df + 2 - 4*exp(-f(x))
    eq4 = f(x).diff(x) - (a4*x**4 + a3*x**3 + a2*x**2 + a1*x + a0)**Rational(-1, 2)
    eq5 = x**2*df - f(x) + x**2*exp(x - (1/x))
    eqlist = [eq, eq1, eq2, eq3, eq4, eq5]

    i = infinitesimals(eq, hint='abaco1_simple')
    assert i == [{eta(x, f(x)): exp(x**3/3), xi(x, f(x)): 0},
        {eta(x, f(x)): f(x), xi(x, f(x)): 0},
        {eta(x, f(x)): 0, xi(x, f(x)): x**(-2)}]
    i1 = infinitesimals(eq1, hint='abaco1_simple')
    assert i1 == [{eta(x, f(x)): exp(-a*x), xi(x, f(x)): 0}]
    i2 = infinitesimals(eq2, hint='abaco1_simple')
    assert i2 == [{eta(x, f(x)): exp(-x**2), xi(x, f(x)): 0}]
    i3 = infinitesimals(eq3, hint='abaco1_simple')
    assert i3 == [{eta(x, f(x)): 0, xi(x, f(x)): 2*x + 1},
        {eta(x, f(x)): 0, xi(x, f(x)): 1/(exp(f(x)) - 2)}]
    i4 = infinitesimals(eq4, hint='abaco1_simple')
    assert i4 == [{eta(x, f(x)): 1, xi(x, f(x)): 0},
        {eta(x, f(x)): 0,
        xi(x, f(x)): sqrt(a0 + a1*x + a2*x**2 + a3*x**3 + a4*x**4)}]
    i5 = infinitesimals(eq5, hint='abaco1_simple')
    assert i5 == [{xi(x, f(x)): 0, eta(x, f(x)): exp(-1/x)}]

    ilist = [i, i1, i2, i3, i4, i5]
    for eq, i in (zip(eqlist, ilist)):
        check = checkinfsol(eq, i)
        assert check[0]


def test_issue_6247():
    eq = x**2*f(x)**2 + x*Derivative(f(x), x)
    sol = Eq(f(x), 2*C1/(C1*x**2 - 1))
    assert dsolve(eq, hint = 'separable_reduced') == sol
    assert checkodesol(eq, sol, order=1)[0]

    eq = f(x).diff(x, x) + 4*f(x)
    sol = Eq(f(x), C1*sin(2*x) + C2*cos(2*x))
    assert dsolve(eq) == sol
    assert checkodesol(eq, sol, order=1)[0]


def test_heuristic2():
    xi = Function('xi')
    eta = Function('eta')
    df = f(x).diff(x)

    # This ODE can be solved by the Lie Group method, when there are
    # better assumptions
    eq = df - (f(x)/x)*(x*log(x**2/f(x)) + 2)
    i = infinitesimals(eq, hint='abaco1_product')
    assert i == [{eta(x, f(x)): f(x)*exp(-x), xi(x, f(x)): 0}]
    assert checkinfsol(eq, i)[0]


@slow
def test_heuristic3():
    xi = Function('xi')
    eta = Function('eta')
    a, b = symbols("a b")
    df = f(x).diff(x)

    eq = x**2*df + x*f(x) + f(x)**2 + x**2
    i = infinitesimals(eq, hint='bivariate')
    assert i == [{eta(x, f(x)): f(x), xi(x, f(x)): x}]
    assert checkinfsol(eq, i)[0]

    eq = x**2*(-f(x)**2 + df)- a*x**2*f(x) + 2 - a*x
    i = infinitesimals(eq, hint='bivariate')
    assert checkinfsol(eq, i)[0]


def test_heuristic_4():
    y, a = symbols("y a")

    eq = x*(f(x).diff(x)) + 1 - f(x)**2
    i = infinitesimals(eq, hint='chi')
    assert checkinfsol(eq, i)[0]


def test_heuristic_function_sum():
    xi = Function('xi')
    eta = Function('eta')
    eq = f(x).diff(x) - (3*(1 + x**2/f(x)**2)*atan(f(x)/x) + (1 - 2*f(x))/x +
       (1 - 3*f(x))*(x/f(x)**2))
    i = infinitesimals(eq, hint='function_sum')
    assert i == [{eta(x, f(x)): f(x)**(-2) + x**(-2), xi(x, f(x)): 0}]
    assert checkinfsol(eq, i)[0]


def test_heuristic_abaco2_similar():
    xi = Function('xi')
    eta = Function('eta')
    F = Function('F')
    a, b = symbols("a b")
    eq = f(x).diff(x) - F(a*x + b*f(x))
    i = infinitesimals(eq, hint='abaco2_similar')
    assert i == [{eta(x, f(x)): -a/b, xi(x, f(x)): 1}]
    assert checkinfsol(eq, i)[0]

    eq = f(x).diff(x) - (f(x)**2 / (sin(f(x) - x) - x**2 + 2*x*f(x)))
    i = infinitesimals(eq, hint='abaco2_similar')
    assert i == [{eta(x, f(x)): f(x)**2, xi(x, f(x)): f(x)**2}]
    assert checkinfsol(eq, i)[0]


def test_heuristic_abaco2_unique_unknown():
    xi = Function('xi')
    eta = Function('eta')
    F = Function('F')
    a, b = symbols("a b")
    x = Symbol("x", positive=True)

    eq = f(x).diff(x) - x**(a - 1)*(f(x)**(1 - b))*F(x**a/a + f(x)**b/b)
    i = infinitesimals(eq, hint='abaco2_unique_unknown')
    assert i == [{eta(x, f(x)): -f(x)*f(x)**(-b), xi(x, f(x)): x*x**(-a)}]
    assert checkinfsol(eq, i)[0]

    eq = f(x).diff(x) + tan(F(x**2 + f(x)**2) + atan(x/f(x)))
    i = infinitesimals(eq, hint='abaco2_unique_unknown')
    assert i == [{eta(x, f(x)): x, xi(x, f(x)): -f(x)}]
    assert checkinfsol(eq, i)[0]

    eq = (x*f(x).diff(x) + f(x) + 2*x)**2 -4*x*f(x) -4*x**2 -4*a
    i = infinitesimals(eq, hint='abaco2_unique_unknown')
    assert checkinfsol(eq, i)[0]


def test_heuristic_linear():
    a, b, m, n = symbols("a b m n")

    eq = x**(n*(m + 1) - m)*(f(x).diff(x)) - a*f(x)**n -b*x**(n*(m + 1))
    i = infinitesimals(eq, hint='linear')
    assert checkinfsol(eq, i)[0]


@XFAIL
def test_kamke():
    a, b, alpha, c = symbols("a b alpha c")
    eq = x**2*(a*f(x)**2+(f(x).diff(x))) + b*x**alpha + c
    i = infinitesimals(eq, hint='sum_function')  # XFAIL
    assert checkinfsol(eq, i)[0]


def test_series():
    C1 = Symbol("C1")
    eq = f(x).diff(x) - f(x)
    sol = Eq(f(x), C1 + C1*x + C1*x**2/2 + C1*x**3/6 + C1*x**4/24 +
            C1*x**5/120 + O(x**6))
    assert dsolve(eq, hint='1st_power_series') == sol
    assert checkodesol(eq, sol, order=1)[0]

    eq = f(x).diff(x) - x*f(x)
    sol = Eq(f(x), C1*x**4/8 + C1*x**2/2 + C1 + O(x**6))
    assert dsolve(eq, hint='1st_power_series') == sol
    assert checkodesol(eq, sol, order=1)[0]

    eq = f(x).diff(x) - sin(x*f(x))
    sol = Eq(f(x), (x - 2)**2*(1+ sin(4))*cos(4) + (x - 2)*sin(4) + 2 + O(x**3))
    assert dsolve(eq, hint='1st_power_series', ics={f(2): 2}, n=3) == sol
    # FIXME: The solution here should be O((x-2)**3) so is incorrect
    #assert checkodesol(eq, sol, order=1)[0]



@XFAIL
@SKIP
def test_lie_group_issue17322_1():
    eq=x*f(x).diff(x)*(f(x)+4) + (f(x)**2) -2*f(x)-2*x
    sol = dsolve(eq, f(x))  # Hangs
    assert checkodesol(eq, sol) == (True, 0)

@XFAIL
@SKIP
def test_lie_group_issue17322_2():
    eq=x*f(x).diff(x)*(f(x)+4) + (f(x)**2) -2*f(x)-2*x
    sol = dsolve(eq)  # Hangs
    assert checkodesol(eq, sol) == (True, 0)

@XFAIL
@SKIP
def test_lie_group_issue17322_3():
    eq=Eq(x**7*Derivative(f(x), x) + 5*x**3*f(x)**2 - (2*x**2 + 2)*f(x)**3, 0)
    sol = dsolve(eq)  # Hangs
    assert checkodesol(eq, sol) == (True, 0)


@XFAIL
def test_lie_group_issue17322_4():
    eq=f(x).diff(x) - (f(x) - x*log(x))**2/x**2 + log(x)
    sol = dsolve(eq)  # NotImplementedError
    assert checkodesol(eq, sol) == (True, 0)


@slow
def test_lie_group():
    C1 = Symbol("C1")
    x = Symbol("x") # assuming x is real generates an error!
    a, b, c = symbols("a b c")
    eq = f(x).diff(x)**2
    sol = dsolve(eq, f(x), hint='lie_group')
    assert checkodesol(eq, sol) == (True, 0)

    eq = Eq(f(x).diff(x), x**2*f(x))
    sol = dsolve(eq, f(x), hint='lie_group')
    assert sol == Eq(f(x), C1*exp(x**3)**Rational(1, 3))
    assert checkodesol(eq, sol) == (True, 0)

    eq = f(x).diff(x) + a*f(x) - c*exp(b*x)
    sol = dsolve(eq, f(x), hint='lie_group')
    assert checkodesol(eq, sol) == (True, 0)

    eq = f(x).diff(x) + 2*x*f(x) - x*exp(-x**2)
    sol = dsolve(eq, f(x), hint='lie_group')
    actual_sol = Eq(f(x), (C1 + x**2/2)*exp(-x**2))
    errstr = str(eq)+' : '+str(sol)+' == '+str(actual_sol)
    assert sol == actual_sol, errstr
    assert checkodesol(eq, sol) == (True, 0)

    eq = (1 + 2*x)*(f(x).diff(x)) + 2 - 4*exp(-f(x))
    sol = dsolve(eq, f(x), hint='lie_group')
    assert sol == Eq(f(x), log(C1/(2*x + 1) + 2))
    assert checkodesol(eq, sol) == (True, 0)

    eq = x**2*(f(x).diff(x)) - f(x) + x**2*exp(x - (1/x))
    sol = dsolve(eq, f(x), hint='lie_group')
    assert checkodesol(eq, sol)[0]

    eq = x**2*f(x)**2 + x*Derivative(f(x), x)
    sol = dsolve(eq, f(x), hint='lie_group')
    assert sol == Eq(f(x), 2/(C1 + x**2))
    assert checkodesol(eq, sol) == (True, 0)

    eq=diff(f(x),x) + 2*x*f(x) - x*exp(-x**2)
    sol = Eq(f(x), exp(-x**2)*(C1 + x**2/2))
    assert sol == dsolve(eq, hint='lie_group')
    assert checkodesol(eq, sol) == (True, 0)

    eq = diff(f(x),x) + f(x)*cos(x) - exp(2*x)
    sol = Eq(f(x), exp(-sin(x))*(C1 + Integral(exp(2*x)*exp(sin(x)), x)))
    assert sol == dsolve(eq, hint='lie_group')
    assert checkodesol(eq, sol) == (True, 0)

    eq = diff(f(x),x) + f(x)*cos(x) - sin(2*x)/2
    sol = Eq(f(x), C1*exp(-sin(x)) + sin(x) - 1)
    assert sol == dsolve(eq, hint='lie_group')
    assert checkodesol(eq, sol) == (True, 0)

    eq = x*diff(f(x),x) + f(x) - x*sin(x)
    sol = Eq(f(x), (C1 - x*cos(x) + sin(x))/x)
    assert sol == dsolve(eq, hint='lie_group')
    assert checkodesol(eq, sol) == (True, 0)

    eq = x*diff(f(x),x) - f(x) - x/log(x)
    sol = Eq(f(x), x*(C1 + log(log(x))))
    assert sol == dsolve(eq, hint='lie_group')
    assert checkodesol(eq, sol) == (True, 0)

    eq = (f(x).diff(x)-f(x)) * (f(x).diff(x)+f(x))
    sol = [Eq(f(x), C1*exp(x)), Eq(f(x), C1*exp(-x))]
    assert set(sol) == set(dsolve(eq, hint='lie_group'))
    assert checkodesol(eq, sol[0]) == (True, 0)
    assert checkodesol(eq, sol[1]) == (True, 0)

    eq = f(x).diff(x) * (f(x).diff(x) - f(x))
    sol = [Eq(f(x), C1*exp(x)), Eq(f(x), C1)]
    assert set(sol) == set(dsolve(eq, hint='lie_group'))
    assert checkodesol(eq, sol[0]) == (True, 0)
    assert checkodesol(eq, sol[1]) == (True, 0)


@XFAIL
def test_lie_group_issue15219():
    eqn = exp(f(x).diff(x)-f(x))
    assert 'lie_group' not in classify_ode(eqn, f(x))


def test_user_infinitesimals():
    x = Symbol("x") # assuming x is real generates an error
    eq = x*(f(x).diff(x)) + 1 - f(x)**2
    sol = Eq(f(x), (C1 + x**2)/(C1 - x**2))
    infinitesimals = {'xi':sqrt(f(x) - 1)/sqrt(f(x) + 1), 'eta':0}
    assert dsolve(eq, hint='lie_group', **infinitesimals) == sol
    assert checkodesol(eq, sol) == (True, 0)


def test_issue_7081():
    eq = x*(f(x).diff(x)) + 1 - f(x)**2
    s = Eq(f(x), -1/(-C1 + x**2)*(C1 + x**2))
    assert dsolve(eq) == s
    assert checkodesol(eq, s) == (True, 0)


@slow
def test_2nd_power_series_ordinary():
    C1, C2 = symbols("C1 C2")

    eq = f(x).diff(x, 2) - x*f(x)
    assert classify_ode(eq) == ('2nd_linear_airy', '2nd_power_series_ordinary')
    sol = Eq(f(x), C2*(x**3/6 + 1) + C1*x*(x**3/12 + 1) + O(x**6))
    assert dsolve(eq, hint='2nd_power_series_ordinary') == sol
    assert checkodesol(eq, sol) == (True, 0)

    sol = Eq(f(x), C2*((x + 2)**4/6 + (x + 2)**3/6 - (x + 2)**2 + 1)
        + C1*(x + (x + 2)**4/12 - (x + 2)**3/3 + S(2))
        + O(x**6))
    assert dsolve(eq, hint='2nd_power_series_ordinary', x0=-2) == sol
    # FIXME: Solution should be O((x+2)**6)
    # assert checkodesol(eq, sol) == (True, 0)

    sol = Eq(f(x), C2*x + C1 + O(x**2))
    assert dsolve(eq, hint='2nd_power_series_ordinary', n=2) == sol
    assert checkodesol(eq, sol) == (True, 0)

    eq = (1 + x**2)*(f(x).diff(x, 2)) + 2*x*(f(x).diff(x)) -2*f(x)
    assert classify_ode(eq) == ('2nd_power_series_ordinary',)
    sol = Eq(f(x), C2*(-x**4/3 + x**2 + 1) + C1*x + O(x**6))
    assert dsolve(eq) == sol
    assert checkodesol(eq, sol) == (True, 0)

    eq = f(x).diff(x, 2) + x*(f(x).diff(x)) + f(x)
    assert classify_ode(eq) == ('2nd_power_series_ordinary',)
    sol = Eq(f(x), C2*(x**4/8 - x**2/2 + 1) + C1*x*(-x**2/3 + 1) + O(x**6))
    assert dsolve(eq) == sol
    # FIXME: checkodesol fails for this solution...
    # assert checkodesol(eq, sol) == (True, 0)

    eq = f(x).diff(x, 2) + f(x).diff(x) - x*f(x)
    assert classify_ode(eq) == ('2nd_power_series_ordinary',)
    sol = Eq(f(x), C2*(-x**4/24 + x**3/6 + 1)
            + C1*x*(x**3/24 + x**2/6 - x/2 + 1) + O(x**6))
    assert dsolve(eq) == sol
    # FIXME: checkodesol fails for this solution...
    # assert checkodesol(eq, sol) == (True, 0)

    eq = f(x).diff(x, 2) + x*f(x)
    assert classify_ode(eq) == ('2nd_linear_airy', '2nd_power_series_ordinary')
    sol = Eq(f(x), C2*(x**6/180 - x**3/6 + 1) + C1*x*(-x**3/12 + 1) + O(x**7))
    assert dsolve(eq, hint='2nd_power_series_ordinary', n=7) == sol
    assert checkodesol(eq, sol) == (True, 0)

def test_Airy_equation():
    eq = f(x).diff(x, 2) - x*f(x)
    sol = Eq(f(x), C1*airyai(x) + C2*airybi(x))
    sols = constant_renumber(sol)
    assert classify_ode(eq) == ("2nd_linear_airy",'2nd_power_series_ordinary')
    assert checkodesol(eq, sol) == (True, 0)
    assert dsolve(eq, f(x)) in (sol, sols)
    assert dsolve(eq, f(x), hint='2nd_linear_airy') in (sol, sols)

    eq = f(x).diff(x, 2) + 2*x*f(x)
    sol = Eq(f(x), C1*airyai(-2**(S(1)/3)*x) + C2*airybi(-2**(S(1)/3)*x))
    sols = constant_renumber(sol)
    assert classify_ode(eq) == ("2nd_linear_airy",'2nd_power_series_ordinary')
    assert checkodesol(eq, sol) == (True, 0)
    assert dsolve(eq, f(x)) in (sol, sols)
    assert dsolve(eq, f(x), hint='2nd_linear_airy') in (sol, sols)

def test_2nd_power_series_regular():
    C1, C2 = symbols("C1 C2")
    eq = x**2*(f(x).diff(x, 2)) - 3*x*(f(x).diff(x)) + (4*x + 4)*f(x)
    sol = Eq(f(x), C1*x**2*(-16*x**3/9 + 4*x**2 - 4*x + 1) + O(x**6))
    assert dsolve(eq, hint='2nd_power_series_regular') == sol
    assert checkodesol(eq, sol) == (True, 0)

    eq = 4*x**2*(f(x).diff(x, 2)) -8*x**2*(f(x).diff(x)) + (4*x**2 +
        1)*f(x)
    sol = Eq(f(x), C1*sqrt(x)*(x**4/24 + x**3/6 + x**2/2 + x + 1) + O(x**6))
    assert dsolve(eq, hint='2nd_power_series_regular') == sol
    assert checkodesol(eq, sol) == (True, 0)

    eq = x**2*(f(x).diff(x, 2)) - x**2*(f(x).diff(x)) + (
        x**2 - 2)*f(x)
    sol = Eq(f(x), C1*(-x**6/720 - 3*x**5/80 - x**4/8 + x**2/2 + x/2 + 1)/x +
            C2*x**2*(-x**3/60 + x**2/20 + x/2 + 1) + O(x**6))
    assert dsolve(eq) == sol
    assert checkodesol(eq, sol) == (True, 0)

    eq = x**2*(f(x).diff(x, 2)) + x*(f(x).diff(x)) + (x**2 - Rational(1, 4))*f(x)
    sol = Eq(f(x), C1*(x**4/24 - x**2/2 + 1)/sqrt(x) +
        C2*sqrt(x)*(x**4/120 - x**2/6 + 1) + O(x**6))
    assert dsolve(eq, hint='2nd_power_series_regular') == sol
    assert checkodesol(eq, sol) == (True, 0)

def test_2nd_linear_bessel_equation():
    eq = x**2*(f(x).diff(x, 2)) + x*(f(x).diff(x)) + (x**2 - 4)*f(x)
    sol = Eq(f(x), C1*besselj(2, x) + C2*bessely(2, x))
    sols = constant_renumber(sol)
    assert dsolve(eq, f(x)) in (sol, sols)
    assert dsolve(eq, f(x), hint='2nd_linear_bessel') in (sol, sols)
    assert checkodesol(eq, sol, order=2, solve_for_func=False) == (True, 0)

    eq = x**2*(f(x).diff(x, 2)) + x*(f(x).diff(x)) + (x**2 +25)*f(x)
    sol = Eq(f(x), C1*besselj(5*I, x) + C2*bessely(5*I, x))
    sols = constant_renumber(sol)
    assert dsolve(eq, f(x)) in (sol, sols)
    assert dsolve(eq, f(x), hint='2nd_linear_bessel') in (sol, sols)
    checkodesol(eq, sol, order=2, solve_for_func=False) == (True, 0)

    eq = x**2*(f(x).diff(x, 2)) + x*(f(x).diff(x)) + (x**2)*f(x)
    sol = Eq(f(x), C1*besselj(0, x) + C2*bessely(0, x))
    sols = constant_renumber(sol)
    assert dsolve(eq, f(x)) in (sol, sols)
    assert dsolve(eq, f(x), hint='2nd_linear_bessel') in (sol, sols)
    assert checkodesol(eq, sol, order=2, solve_for_func=False) == (True, 0)

    eq = x**2*(f(x).diff(x, 2)) + x*(f(x).diff(x)) + (81*x**2 -S(1)/9)*f(x)
    sol = Eq(f(x), C1*besselj(S(1)/3, 9*x) + C2*bessely(S(1)/3, 9*x))
    sols = constant_renumber(sol)
    assert dsolve(eq, f(x)) in (sol, sols)
    assert dsolve(eq, f(x), hint='2nd_linear_bessel') in (sol, sols)
    checkodesol(eq, sol, order=2, solve_for_func=False) == (True, 0)

    eq = x**2*(f(x).diff(x, 2)) + x*(f(x).diff(x)) + (x**4 - 4)*f(x)
    sol = Eq(f(x), C1*besselj(1, x**2/2) + C2*bessely(1, x**2/2))
    sols = constant_renumber(sol)
    assert dsolve(eq, f(x)) in (sol, sols)
    assert dsolve(eq, f(x), hint='2nd_linear_bessel') in (sol, sols)
    assert checkodesol(eq, sol, order=2, solve_for_func=False) == (True, 0)

    eq = x**2*(f(x).diff(x, 2)) + 2*x*(f(x).diff(x)) + (x**4 - 4)*f(x)
    sol = Eq(f(x), (C1*besselj(sqrt(17)/4, x**2/2) + C2*bessely(sqrt(17)/4, x**2/2))/sqrt(x))
    sols = constant_renumber(sol)
    assert dsolve(eq, f(x)) in (sol, sols)
    assert dsolve(eq, f(x), hint='2nd_linear_bessel') in (sol, sols)
    assert checkodesol(eq, sol, order=2, solve_for_func=False) == (True, 0)

    eq = x**2*(f(x).diff(x, 2)) + x*(f(x).diff(x)) + (x**2 - S(1)/4)*f(x)
    sol = Eq(f(x), C1*besselj(S(1)/2, x) + C2*bessely(S(1)/2, x))
    sols = constant_renumber(sol)
    assert dsolve(eq, f(x)) in (sol, sols)
    assert dsolve(eq, f(x), hint='2nd_linear_bessel') in (sol, sols)
    assert checkodesol(eq, sol, order=2, solve_for_func=False) == (True, 0)

    eq = x**2*(f(x).diff(x, 2)) - 3*x*(f(x).diff(x)) + (4*x + 4)*f(x)
    sol = Eq(f(x), x**2*(C1*besselj(0, 4*sqrt(x)) + C2*bessely(0, 4*sqrt(x))))
    sols = constant_renumber(sol)
    assert dsolve(eq, f(x)) in (sol, sols)
    assert dsolve(eq, f(x), hint='2nd_linear_bessel') in (sol, sols)
    assert checkodesol(eq, sol, order=2, solve_for_func=False) == (True, 0)

    eq = x*(f(x).diff(x, 2)) - f(x).diff(x) + 4*x**3*f(x)
    sol = Eq(f(x), x*(C1*besselj(S(1)/2, x**2) + C2*bessely(S(1)/2, x**2)))
    sols = constant_renumber(sol)
    assert dsolve(eq, f(x)) in (sol, sols)
    assert dsolve(eq, f(x), hint='2nd_linear_bessel') in (sol, sols)
    assert checkodesol(eq, sol, order=2, solve_for_func=False) == (True, 0)

    eq = (x-2)**2*(f(x).diff(x, 2)) - (x-2)*f(x).diff(x) + 4*(x-2)**2*f(x)
    sol = Eq(f(x), (x - 2)*(C1*besselj(1, 2*x - 4) + C2*bessely(1, 2*x - 4)))
    sols = constant_renumber(sol)
    assert dsolve(eq, f(x)) in (sol, sols)
    assert dsolve(eq, f(x), hint='2nd_linear_bessel') in (sol, sols)
    assert checkodesol(eq, sol, order=2, solve_for_func=False) == (True, 0)

def test_issue_7093():
    x = Symbol("x") # assuming x is real leads to an error
    sol = [Eq(f(x), C1 - 2*x*sqrt(x**3)/5),
           Eq(f(x), C1 + 2*x*sqrt(x**3)/5)]
    eq = Derivative(f(x), x)**2 - x**3
    assert set(dsolve(eq)) == set(sol)
    assert checkodesol(eq, sol) == [(True, 0)] * 2


def test_dsolve_linsystem_symbol():
    eps = Symbol('epsilon', positive=True)
    eq1 = (Eq(diff(f(x), x), -eps*g(x)), Eq(diff(g(x), x), eps*f(x)))
    sol1 = [Eq(f(x), -C1*eps*cos(eps*x) - C2*eps*sin(eps*x)),
            Eq(g(x), -C1*eps*sin(eps*x) + C2*eps*cos(eps*x))]
    assert checksysodesol(eq1, sol1) == (True, [0, 0])


def test_C1_function_9239():
    t = Symbol('t')
    C1 = Function('C1')
    C2 = Function('C2')
    C3 = Symbol('C3')
    C4 = Symbol('C4')
    eq = (Eq(diff(C1(t), t), 9*C2(t)), Eq(diff(C2(t), t), 12*C1(t)))
    sol = [Eq(C1(t), 9*C3*exp(6*sqrt(3)*t) + 9*C4*exp(-6*sqrt(3)*t)),
           Eq(C2(t), 6*sqrt(3)*C3*exp(6*sqrt(3)*t) - 6*sqrt(3)*C4*exp(-6*sqrt(3)*t))]
    assert checksysodesol(eq, sol) == (True, [0, 0])


def test_issue_15056():
    t = Symbol('t')
    C3 = Symbol('C3')
    assert get_numbered_constants(Symbol('C1') * Function('C2')(t)) == C3


def test_issue_10379():
    t,y = symbols('t,y')
    eq = f(t).diff(t)-(1-51.05*y*f(t))
    sol =  Eq(f(t), (0.019588638589618*exp(y*(C1 - 51.05*t)) + 0.019588638589618)/y)
    dsolve_sol = dsolve(eq, rational=False)
    assert str(dsolve_sol) == str(sol)
    assert checkodesol(eq, dsolve_sol)[0]


def test_issue_10867():
    x = Symbol('x')
    eq = Eq(g(x).diff(x).diff(x), (x-2)**2 + (x-3)**3)
    sol = Eq(g(x), C1 + C2*x + x**5/20 - 2*x**4/3 + 23*x**3/6 - 23*x**2/2)
    assert dsolve(eq, g(x)) == sol
    assert checkodesol(eq, sol, order=2, solve_for_func=False) == (True, 0)


def test_issue_11290():
    eq = cos(f(x)) - (x*sin(f(x)) - f(x)**2)*f(x).diff(x)
    sol_1 = dsolve(eq, f(x), simplify=False, hint='1st_exact_Integral')
    sol_0 = dsolve(eq, f(x), simplify=False, hint='1st_exact')
    assert sol_1.dummy_eq(Eq(Subs(
        Integral(u**2 - x*sin(u) - Integral(-sin(u), x), u) +
        Integral(cos(u), x), u, f(x)), C1))
    assert sol_1.doit() == sol_0
    assert checkodesol(eq, sol_0, order=1, solve_for_func=False)
    assert checkodesol(eq, sol_1, order=1, solve_for_func=False)


def test_issue_4838():
    # Issue #15999
    eq = f(x).diff(x) - C1*f(x)
    sol = Eq(f(x), C2*exp(C1*x))
    assert dsolve(eq, f(x)) == sol
    assert checkodesol(eq, sol, order=1, solve_for_func=False) == (True, 0)

    # Issue #13691
    eq = f(x).diff(x) - C1*g(x).diff(x)
    sol = Eq(f(x), C2 + C1*g(x))
    assert dsolve(eq, f(x)) == sol
    assert checkodesol(eq, sol, f(x), order=1, solve_for_func=False) == (True, 0)

    # Issue #4838
    eq = f(x).diff(x) - 3*C1 - 3*x**2
    sol = Eq(f(x), C2 + 3*C1*x + x**3)
    assert dsolve(eq, f(x)) == sol
    assert checkodesol(eq, sol, order=1, solve_for_func=False) == (True, 0)


@slow
def test_issue_14395():
    eq = Derivative(f(x), x, x) + 9*f(x) - sec(x)
    sol = Eq(f(x), (C1 - x/3 + sin(2*x)/3)*sin(3*x) + (C2 + log(cos(x))
        - 2*log(cos(x)**2)/3 + 2*cos(x)**2/3)*cos(3*x))
    assert dsolve(eq, f(x)) == sol
    # FIXME: assert checkodesol(eq, sol, order=2, solve_for_func=False) == (True, 0)


@slow
def test_nth_order_reducible():
    from sympy.solvers.ode.ode import _nth_order_reducible_match

    eqn = Eq(x*Derivative(f(x), x)**2 + Derivative(f(x), x, 2), 0)
    sol = Eq(f(x),
             C1 - sqrt(-1/C2)*log(-C2*sqrt(-1/C2) + x) + sqrt(-1/C2)*log(C2*sqrt(-1/C2) + x))
    assert checkodesol(eqn, sol, order=2, solve_for_func=False) == (True, 0)
    assert sol == dsolve(eqn, f(x), hint='nth_order_reducible')
    assert sol == dsolve(eqn, f(x))

    F = lambda eq: _nth_order_reducible_match(eq, f(x))
    D = Derivative
    assert F(D(y*f(x), x, y) + D(f(x), x)) is None
    assert F(D(y*f(y), y, y) + D(f(y), y)) is None
    assert F(f(x)*D(f(x), x) + D(f(x), x, 2)) is None
    assert F(D(x*f(y), y, 2) + D(u*y*f(x), x, 3)) is None  # no simplification by design
    assert F(D(f(y), y, 2) + D(f(y), y, 3) + D(f(x), x, 4)) is None
    assert F(D(f(x), x, 2) + D(f(x), x, 3)) == dict(n=2)

    eqn = -exp(x) + (x*Derivative(f(x), (x, 2)) + Derivative(f(x), x))/x
    sol = Eq(f(x), C1 + C2*log(x) + exp(x) - Ei(x))
    assert checkodesol(eqn, sol, order=2, solve_for_func=False) == (True, 0)
    assert sol == dsolve(eqn, f(x))
    assert sol == dsolve(eqn, f(x), hint='nth_order_reducible')

    eqn = Eq(sqrt(2) * f(x).diff(x,x,x) + f(x).diff(x), 0)
    sol = Eq(f(x), C1 + C2*sin(2**Rational(3, 4)*x/2) + C3*cos(2**Rational(3, 4)*x/2))
    assert checkodesol(eqn, sol, order=2, solve_for_func=False) == (True, 0)
    assert sol == dsolve(eqn, f(x))
    assert sol == dsolve(eqn, f(x), hint='nth_order_reducible')

    eqn = f(x).diff(x, 2) + 2*f(x).diff(x)
    sol = Eq(f(x), C1 + C2*exp(-2*x))
    sols = constant_renumber(sol)
    assert checkodesol(eqn, sol, order=2, solve_for_func=False) == (True, 0)
    assert dsolve(eqn, f(x)) in (sol, sols)
    assert dsolve(eqn, f(x), hint='nth_order_reducible') in (sol, sols)

    eqn = f(x).diff(x, 3) + f(x).diff(x, 2) - 6*f(x).diff(x)
    sol = Eq(f(x), C1 + C2*exp(-3*x) + C3*exp(2*x))
    sols = constant_renumber(sol)
    assert checkodesol(eqn, sol, order=2, solve_for_func=False) == (True, 0)
    assert dsolve(eqn, f(x)) in (sol, sols)
    assert dsolve(eqn, f(x), hint='nth_order_reducible') in (sol, sols)

    eqn = f(x).diff(x, 4) - f(x).diff(x, 3) - 4*f(x).diff(x, 2) + \
        4*f(x).diff(x)
    sol = Eq(f(x), C1 + C2*exp(x) + C3*exp(-2*x) + C4*exp(2*x))
    sols = constant_renumber(sol)
    assert checkodesol(eqn, sol, order=2, solve_for_func=False) == (True, 0)
    assert dsolve(eqn, f(x)) in (sol, sols)
    assert dsolve(eqn, f(x), hint='nth_order_reducible') in (sol, sols)

    eqn = f(x).diff(x, 4) + 3*f(x).diff(x, 3)
    sol = Eq(f(x), C1 + C2*x + C3*x**2 + C4*exp(-3*x))
    sols = constant_renumber(sol)
    assert checkodesol(eqn, sol, order=2, solve_for_func=False) == (True, 0)
    assert dsolve(eqn, f(x)) in (sol, sols)
    assert dsolve(eqn, f(x), hint='nth_order_reducible') in (sol, sols)

    eqn = f(x).diff(x, 4) - 2*f(x).diff(x, 2)
    sol = Eq(f(x), C1 + C2*x + C3*exp(x*sqrt(2)) + C4*exp(-x*sqrt(2)))
    sols = constant_renumber(sol)
    assert checkodesol(eqn, sol, order=2, solve_for_func=False) == (True, 0)
    assert dsolve(eqn, f(x)) in (sol, sols)
    assert dsolve(eqn, f(x), hint='nth_order_reducible') in (sol, sols)

    eqn = f(x).diff(x, 4) + 4*f(x).diff(x, 2)
    sol = Eq(f(x), C1 + C2*sin(2*x) + C3*cos(2*x) + C4*x)
    sols = constant_renumber(sol)
    assert checkodesol(eqn, sol, order=2, solve_for_func=False) == (True, 0)
    assert dsolve(eqn, f(x)) in (sol, sols)
    assert dsolve(eqn, f(x), hint='nth_order_reducible') in (sol, sols)

    eqn = f(x).diff(x, 5) + 2*f(x).diff(x, 3) + f(x).diff(x)
    # These are equivalent:
    sol1 = Eq(f(x), C1 + (C2 + C3*x)*sin(x) + (C4 + C5*x)*cos(x))
    sol2 = Eq(f(x), C1 + C2*(x*sin(x) + cos(x)) + C3*(-x*cos(x) + sin(x)) + C4*sin(x) + C5*cos(x))
    sol1s = constant_renumber(sol1)
    sol2s = constant_renumber(sol2)
    assert checkodesol(eqn, sol1, order=2, solve_for_func=False) == (True, 0)
    assert checkodesol(eqn, sol2, order=2, solve_for_func=False) == (True, 0)
    assert dsolve(eqn, f(x)) in (sol1, sol1s)
    assert dsolve(eqn, f(x), hint='nth_order_reducible') in (sol2, sol2s)

    # In this case the reduced ODE has two distinct solutions
    eqn = f(x).diff(x, 2) - f(x).diff(x)**3
    sol = [Eq(f(x), C2 - sqrt(2)*I*(C1 + x)*sqrt(1/(C1 + x))),
           Eq(f(x), C2 + sqrt(2)*I*(C1 + x)*sqrt(1/(C1 + x)))]
    sols = constant_renumber(sol)
    assert checkodesol(eqn, sol, order=2, solve_for_func=False) == [(True, 0), (True, 0)]
    assert dsolve(eqn, f(x)) in (sol, sols)
    assert dsolve(eqn, f(x), hint='nth_order_reducible') in (sol, sols)


# Needs to be a way to know how to combine derivatives in the expression
def test_factoring_ode():
    from sympy import Mul
    eqn = Derivative(x*f(x), x, x, x) + Derivative(f(x), x, x, x)
    # 2-arg Mul!
    soln = Eq(f(x), C1 + C2*x + C3/Mul(2, (x + 1), evaluate=False))
    assert checkodesol(eqn, soln, order=2, solve_for_func=False)[0]
    assert soln == dsolve(eqn, f(x))


def test_issue_11542():
    m = 96
    g = 9.8
    k = .2
    f1 = g * m
    t = Symbol('t')
    v = Function('v')
    v_equation = dsolve(f1 - k * (v(t) ** 2) - m * Derivative(v(t)), 0)
    assert str(v_equation) == \
        'Eq(v(t), -68.585712797929/tanh(C1 - 0.142886901662352*t))'


def test_issue_15913():
    eq = -C1/x - 2*x*f(x) - f(x) + Derivative(f(x), x)
    sol = C2*exp(x**2 + x) + exp(x**2 + x)*Integral(C1*exp(-x**2 - x)/x, x)
    assert checkodesol(eq, sol) == (True, 0)
    sol = C1 + C2*exp(-x*y)
    eq = Derivative(y*f(x), x) + f(x).diff(x, 2)
    assert checkodesol(eq, sol, f(x)) == (True, 0)


def test_issue_16146():
    raises(ValueError, lambda: dsolve([f(x).diff(x), g(x).diff(x)], [f(x), g(x), h(x)]))
    raises(ValueError, lambda: dsolve([f(x).diff(x), g(x).diff(x)], [f(x)]))


def test_dsolve_remove_redundant_solutions():

    eq = (f(x)-2)*f(x).diff(x)
    sol = Eq(f(x), C1)
    assert dsolve(eq) == sol

    eq = (f(x)-sin(x))*(f(x).diff(x, 2))
    sol = {Eq(f(x), C1 + C2*x), Eq(f(x), sin(x))}
    assert set(dsolve(eq)) == sol

    eq = (f(x)**2-2*f(x)+1)*f(x).diff(x, 3)
    sol = Eq(f(x), C1 + C2*x + C3*x**2)
    assert dsolve(eq) == sol


def test_issue_17322():
    eq = (f(x).diff(x)-f(x)) * (f(x).diff(x)+f(x))
    sol = [Eq(f(x), C1*exp(-x)), Eq(f(x), C1*exp(x))]
    assert set(sol) == set(dsolve(eq, hint='lie_group'))
    assert checkodesol(eq, sol) == 2*[(True, 0)]

    eq = f(x).diff(x)*(f(x).diff(x)+f(x))
    sol = [Eq(f(x), C1), Eq(f(x), C1*exp(-x))]
    assert set(sol) == set(dsolve(eq, hint='lie_group'))
    assert checkodesol(eq, sol) == 2*[(True, 0)]

def test_2nd_2F1_hypergeometric():

    eq = x*(x-1)*f(x).diff(x, 2) + (S(3)/2 -2*x)*f(x).diff(x) + 2*f(x)
    sol = Eq(f(x), C1*x**(S(5)/2)*hyper((S(3)/2, S(1)/2), (S(7)/2,), x) + C2*hyper((-1, -2), (-S(3)/2,), x))
    assert sol == dsolve(eq, hint='2nd_hypergeometric')
    assert checkodesol(eq, sol) == (True, 0)

    eq = x*(x-1)*f(x).diff(x, 2) + (S(7)/2*x)*f(x).diff(x) + f(x)
    sol = Eq(f(x), (C1*(1 - x)**(S(5)/2)*hyper((S(1)/2, 2), (S(7)/2,), 1 - x) +
          C2*hyper((-S(1)/2, -2), (-S(3)/2,), 1 - x))/(x - 1)**(S(5)/2))
    assert sol == dsolve(eq, hint='2nd_hypergeometric')
    assert checkodesol(eq, sol) == (True, 0)

    eq = x*(x-1)*f(x).diff(x, 2) + (S(3)+ S(7)/2*x)*f(x).diff(x) + f(x)
    sol = Eq(f(x), (C1*(1 - x)**(S(11)/2)*hyper((S(1)/2, 2), (S(13)/2,), 1 - x) +
          C2*hyper((-S(7)/2, -5), (-S(9)/2,), 1 - x))/(x - 1)**(S(11)/2))
    assert sol == dsolve(eq, hint='2nd_hypergeometric')
    assert checkodesol(eq, sol) == (True, 0)

    eq = x*(x-1)*f(x).diff(x, 2) + (-1+ S(7)/2*x)*f(x).diff(x) + f(x)
    sol = Eq(f(x), (C1 + C2*Integral(exp(Integral((1 - x/2)/(x*(x - 1)), x))/(1 -
          x/2)**2, x))*exp(Integral(1/(x - 1), x)/4)*exp(-Integral(7/(x -
          1), x)/4)*hyper((S(1)/2, -1), (1,), x))
    assert sol == dsolve(eq, hint='2nd_hypergeometric_Integral')
    assert checkodesol(eq, sol) == (True, 0)


    eq = -x**(S(5)/7)*(-416*x**(S(9)/7)/9 - 2385*x**(S(5)/7)/49 + S(298)*x/3)*f(x)/(196*(-x**(S(6)/7) +
         x)**2*(x**(S(6)/7) + x)**2) + Derivative(f(x), (x, 2))
    sol = Eq(f(x), x**(S(45)/98)*(C1*x**(S(4)/49)*hyper((S(1)/3, -S(1)/2), (S(9)/7,), x**(S(2)/7)) +
          C2*hyper((S(1)/21, -S(11)/14), (S(5)/7,), x**(S(2)/7)))/(x**(S(2)/7) - 1)**(S(19)/84))
    assert sol == dsolve(eq, hint='2nd_hypergeometric')
    # assert checkodesol(eq, sol) == (True, 0) #issue-https://github.com/sympy/sympy/issues/17702


def test_issue_5096():
    eq = f(x).diff(x, x) + f(x) - x*sin(x - 2)
    sol = Eq(f(x), C1*sin(x) + C2*cos(x) - x**2*cos(x - 2)/4 + x*sin(x - 2)/4)
    assert sol == dsolve(eq, hint='nth_linear_constant_coeff_undetermined_coefficients')
    assert checkodesol(eq, sol) == (True, 0)

    eq = f(x).diff(x, 2) + f(x) - x**4*sin(x-1)
    sol = Eq(f(x), C1*sin(x) + C2*cos(x) - x**5*cos(x - 1)/10 + x**4*sin(x - 1)/4 + x**3*cos(x - 1)/2 - 3*x**2*sin(x - 1)/4 - 3*x*cos(x - 1)/4)
    assert sol == dsolve(eq, hint='nth_linear_constant_coeff_undetermined_coefficients')
    assert checkodesol(eq, sol) == (True, 0)

    eq = f(x).diff(x, 2) - f(x) - exp(x - 1)
    sol = Eq(f(x), C2*exp(-x) + (C1 + x*exp(-1)/2)*exp(x))
    got = dsolve(eq, hint='nth_linear_constant_coeff_undetermined_coefficients')
    assert sol == got, got
    assert checkodesol(eq, sol) == (True, 0)

    eq = f(x).diff(x, 2)+f(x)-(sin(x-2)+1)
    sol = Eq(f(x), C1*sin(x) + C2*cos(x) - x*cos(x - 2)/2 + 1)
    assert sol == dsolve(eq, hint='nth_linear_constant_coeff_undetermined_coefficients')
    assert checkodesol(eq, sol) == (True, 0)

    eq = 2*x**2*f(x).diff(x, 2) + f(x) + sqrt(2*x)*sin(log(2*x)/2)
    sol = Eq(f(x), sqrt(x)*(C1*sin(log(x)/2) + C2*cos(log(x)/2) + sqrt(2)*log(x)*cos(log(2*x)/2)/2))
    assert sol == dsolve(eq, hint='nth_linear_euler_eq_nonhomogeneous_undetermined_coefficients')
    assert checkodesol(eq, sol) == (True, 0)

    eq = 2*x**2*f(x).diff(x, 2) + f(x) + sin(log(2*x)/2)
    sol = Eq(f(x), C1*sqrt(x)*sin(log(x)/2) + C2*sqrt(x)*cos(log(x)/2) - 2*sin(log(2*x)/2)/5 - 4*cos(log(2*x)/2)/5)
    assert sol == dsolve(eq, hint='nth_linear_euler_eq_nonhomogeneous_undetermined_coefficients')
    assert checkodesol(eq, sol) == (True, 0)


def test_issue_15996():
    eq = f(x).diff(x, 5) + 2*f(x).diff(x, 3) + f(x).diff(x) - 2*x - exp(I*x)
    sol = Eq(f(x), C1 + x**2 + (C2 + x*(C3 - x/8 + 3*exp(I*x)/2 + 3*exp(-I*x)/2) + 5*exp(2*I*x)/16 + 2*I*exp(I*x) - 2*I*exp(-I*x))*sin(x) + (C4 + x*(C5 + I*x/8 + 3*I*exp(I*x)/2 - 3*I*exp(-I*x)/2) + 5*I*exp(2*I*x)/16 - 2*exp(I*x) - 2*exp(-I*x))*cos(x) - I*exp(I*x))
    got = dsolve(eq, hint='nth_linear_constant_coeff_variation_of_parameters')
    assert sol == got, got
    assert checkodesol(eq, sol) == (True, 0)

    eq = f(x).diff(x, 5) + 2*f(x).diff(x, 3) + f(x).diff(x) - exp(I*x)
    sol = Eq(f(x), C1 + (C2 + x*(C3 - x/8) + 5*exp(2*I*x)/16)*sin(x) + (C4 + x*(C5 + I*x/8) + 5*I*exp(2*I*x)/16)*cos(x) - I*exp(I*x))
    got = dsolve(eq, hint='nth_linear_constant_coeff_variation_of_parameters')
    assert sol == got, got
    assert checkodesol(eq, sol) == (True, 0)


def test_issue_18408():
    eq = f(x).diff(x, 3) - f(x).diff(x) - sinh(x)
    sol = Eq(f(x), C1 + C2*exp(-x) + C3*exp(x) + x*sinh(x)/2)
    assert sol == dsolve(eq, hint='nth_linear_constant_coeff_undetermined_coefficients')
    assert checkodesol(eq, sol) == (True, 0)

    eq = f(x).diff(x, 2) - 49*f(x) - sinh(3*x)
    sol = Eq(f(x), C1*exp(-7*x) + C2*exp(7*x) - sinh(3*x)/40)
    assert sol == dsolve(eq, hint='nth_linear_constant_coeff_undetermined_coefficients')
    assert checkodesol(eq, sol) == (True, 0)

    eq = f(x).diff(x, 3) - f(x).diff(x) - sinh(x) - exp(x)
    sol = Eq(f(x), C1 + C3*exp(-x) + x*sinh(x)/2 + (C2 + x/2)*exp(x))
    assert sol == dsolve(eq, hint='nth_linear_constant_coeff_undetermined_coefficients')
    assert checkodesol(eq, sol) == (True, 0)


def test_issue_9446():
    f = Function('f')
    assert dsolve(Eq(f(2 * x), sin(Derivative(f(x)))), f(x)) == \
    [Eq(f(x), C1 + pi*x - Integral(asin(f(2*x)), x)), Eq(f(x), C1 + Integral(asin(f(2*x)), x))]

    assert integrate(-asin(f(2*x)+pi), x) == -Integral(asin(pi + f(2*x)), x)
