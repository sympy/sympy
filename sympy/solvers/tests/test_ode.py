from __future__ import division


from sympy import (acos, acosh, asinh, atan, cos, Derivative, diff, dsolve,
    Dummy, Eq, erf, erfi, exp, Function, I, Integral, LambertW, log, O, pi,
    Rational, RootOf, S, simplify, sin, sqrt, Symbol, tan, asin,
    Piecewise, symbols, Poly)
from sympy.solvers.ode import (_undetermined_coefficients_match, checkodesol,
    classify_ode, constant_renumber, constantsimp,
    homogeneous_order, infinitesimals, checkinfsol)
from sympy.solvers.deutils import ode_order
from sympy.utilities.pytest import XFAIL, skip, raises, slow

C1, C2, C3, C4, C5, C6, C7, C8, C9, C10 = symbols('C1:11')
x, y, z = symbols('x:z', real=True)
f = Function('f')
g = Function('g')

# Note: the tests below may fail (but still be correct) if ODE solver,
# the integral engine, solve(), or even simplify() changes. Also, in
# differently formatted solutions, the arbitrary constants might not be
# equal.  Using specific hints in tests can help to avoid this.

# Tests of order higher than 1 should run the solutions through
# constant_renumber because it will normalize it (constant_renumber causes
# dsolve() to return different results on different machines)


def test_checkodesol():
    # For the most part, checkodesol is well tested in the tests below.
    # These tests only handle cases not checked below.
    raises(ValueError, lambda: checkodesol(f(x, y).diff(x), Eq(f(x, y), x)))
    raises(ValueError, lambda: checkodesol(f(x).diff(x), Eq(f(x, y),
           x), f(x, y)))
    assert checkodesol(f(x).diff(x), Eq(f(x, y), x)) == \
        (False, -f(x).diff(x) + f(x, y).diff(x) - 1)
    assert checkodesol(f(x).diff(x), Eq(f(x), x)) is not True
    assert checkodesol(f(x).diff(x), Eq(f(x), x)) == (False, 1)
    sol1 = Eq(f(x)**5 + 11*f(x) - 2*f(x) + x, 0)
    assert checkodesol(diff(sol1.lhs, x), sol1) == (True, 0)
    assert checkodesol(diff(sol1.lhs, x)*exp(f(x)), sol1) == (True, 0)
    assert checkodesol(diff(sol1.lhs, x, 2), sol1) == (True, 0)
    assert checkodesol(diff(sol1.lhs, x, 2)*exp(f(x)), sol1) == (True, 0)
    assert checkodesol(diff(sol1.lhs, x, 3), sol1) == (True, 0)
    assert checkodesol(diff(sol1.lhs, x, 3)*exp(f(x)), sol1) == (True, 0)
    assert checkodesol(diff(sol1.lhs, x, 3), Eq(f(x), x*log(x))) == \
        (False, 60*x**4*((log(x) + 1)**2 + log(x))*(
        log(x) + 1)*log(x)**2 - 5*x**4*log(x)**4 - 9)
    assert checkodesol(diff(exp(f(x)) + x, x)*x, Eq(exp(f(x)) + x)) == \
        (True, 0)
    assert checkodesol(diff(exp(f(x)) + x, x)*x, Eq(exp(f(x)) + x),
        solve_for_func=False) == (True, 0)
    assert checkodesol(f(x).diff(x, 2), [Eq(f(x), C1 + C2*x),
        Eq(f(x), C2 + C1*x), Eq(f(x), C1*x + C2*x**2)]) == \
        [(True, 0), (True, 0), (False, 2*C2)]
    assert checkodesol(f(x).diff(x, 2), set([Eq(f(x), C1 + C2*x),
        Eq(f(x), C2 + C1*x), Eq(f(x), C1*x + C2*x**2)])) == \
        set([(True, 0), (True, 0), (False, 2*C2)])
    assert checkodesol(f(x).diff(x) - 1/f(x)/2, Eq(f(x)**2, x)) == \
        [(True, 0), (True, 0)]
    assert checkodesol(f(x).diff(x) - f(x), Eq(C1*exp(x), f(x))) == (True, 0)
    # Based on test_1st_homogeneous_coeff_ode2_eq3sol.  Make sure that
    # checkodesol tries back substituting f(x) when it can.
    eq3 = x*exp(f(x)/x) + f(x) - x*f(x).diff(x)
    sol3 = Eq(f(x), log(log(C1/x)**(-x)))
    assert not checkodesol(eq3, sol3)[1].has(f(x))


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
        '1st_linear_Integral', 'almost_linear', 'almost_linear_Integral',
        'best', 'best_hint', 'default', 'lie_group',
        'nth_linear_euler_eq_homogeneous', 'order',
        'separable', 'separable_Integral']
    Integral_keys = ['1st_exact_Integral',
        '1st_homogeneous_coeff_subs_dep_div_indep_Integral',
        '1st_homogeneous_coeff_subs_indep_div_dep_Integral', '1st_linear_Integral',
        'almost_linear_Integral', 'best', 'best_hint', 'default',
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
        ('nth_linear_constant_coeff_homogeneous', 'Liouville',
            '2nd_power_series_ordinary' ,'Liouville_Integral')
    assert classify_ode(f(x), f(x)) == ()
    assert classify_ode(Eq(f(x).diff(x), 0), f(x)) == ('separable',
        '1st_linear', '1st_homogeneous_coeff_best',
        '1st_homogeneous_coeff_subs_indep_div_dep',
        '1st_homogeneous_coeff_subs_dep_div_indep',
        '1st_power_series', 'lie_group',
        'nth_linear_constant_coeff_homogeneous',
        'separable_Integral',
        '1st_linear_Integral',
        '1st_homogeneous_coeff_subs_indep_div_dep_Integral',
        '1st_homogeneous_coeff_subs_dep_div_indep_Integral')
    assert classify_ode(f(x).diff(x)**2, f(x)) == ('lie_group',)
    # 1650: f(x) should be cleared from highest derivative before classifying
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
    assert b == c != ()
    assert classify_ode(
        2*x*f(x)*f(x).diff(x) + (1 + x)*f(x)**2 - exp(x), f(x)
    ) == ('Bernoulli', 'almost_linear', 'lie_group',
        'Bernoulli_Integral', 'almost_linear_Integral')
    assert 'Riccati_special_minus2' in \
        classify_ode(2*f(x).diff(x) + f(x)**2 - f(x)/x + 3*x**(-2), f(x))
    raises(ValueError, lambda: classify_ode(x + f(x, y).diff(x).diff(
        y), f(x, y)))
    # 2077
    k = Symbol('k')
    assert classify_ode(f(x).diff(x)/(k*f(x) + k*x*f(x)) + 2*f(x)/(k*f(x) +
        k*x*f(x)) + x*f(x).diff(x)/(k*f(x) + k*x*f(x)) + z, f(x)) == \
        ('separable', '1st_exact', '1st_power_series', 'lie_group',
         'separable_Integral', '1st_exact_Integral')
    # preprocessing
    ans = ('separable', '1st_exact', '1st_linear', 'Bernoulli',
        '1st_homogeneous_coeff_best',
        '1st_homogeneous_coeff_subs_indep_div_dep',
        '1st_homogeneous_coeff_subs_dep_div_indep',
        'separable_reduced', '1st_power_series', 'lie_group',
        'nth_linear_constant_coeff_undetermined_coefficients',
        'nth_linear_constant_coeff_variation_of_parameters',
        'separable_Integral', '1st_exact_Integral',
        '1st_linear_Integral',
        'Bernoulli_Integral',
        '1st_homogeneous_coeff_subs_indep_div_dep_Integral',
        '1st_homogeneous_coeff_subs_dep_div_indep_Integral',
        'separable_reduced_Integral',
        'nth_linear_constant_coeff_variation_of_parameters_Integral')
    #     w/o f(x) given
    assert classify_ode(diff(f(x) + x, x) + diff(f(x), x)) == ans
    #     w/ f(x) and prep=True
    assert classify_ode(diff(f(x) + x, x) + diff(f(x), x), f(x),
                        prep=True) == ans


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
    # issue 2736: ode_order has to also work for unevaluated derivatives
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
    sol2 = Eq(f(x), C1 + 5*x/3)
    sol3 = Eq(f(x), C1 + 5*x/3)
    sol4 = Eq(f(x), C1*sin(x/3) + C2*cos(x/3))
    sol5 = Eq(f(x), C1*exp(-x/3) + C2*exp(x/3))
    sol6 = Eq(f(x), (C1 - cos(x))/x**3)
    sol7 = Eq(f(x), C1*exp(x) + C2*exp(2*x))
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


def test_1st_linear():
    # Type: first order linear form f'(x)+p(x)f(x)=q(x)
    eq = Eq(f(x).diff(x) + x*f(x), x**2)
    sol = Eq(f(x), (C1 + x*exp(x**2/2)
                    - sqrt(2)*sqrt(pi)*erfi(sqrt(2)*x/2)/2)*exp(-x**2/2))
    assert dsolve(eq, hint='1st_linear') == sol
    assert checkodesol(eq, sol, order=1, solve_for_func=False)[0]


def test_Bernoulli():
    # Type: Bernoulli, f'(x) + p(x)*f(x) == q(x)*f(x)**n
    eq = Eq(x*f(x).diff(x) + f(x) - f(x)**2, 0)
    sol = dsolve(eq, f(x), hint='Bernoulli')
    assert sol == Eq(f(x), 1/(x*(C1 + 1/x)))
    assert checkodesol(eq, sol, order=1, solve_for_func=False)[0]


def test_Riccati_special_minus2():
    # Type: Riccati special alpha = -2, a*dy/dx + b*y**2 + c*y/x +d/x**2
    eq = 2*f(x).diff(x) + f(x)**2 - f(x)/x + 3*x**(-2)
    sol = dsolve(eq, f(x), hint='Riccati_special_minus2')
    assert checkodesol(eq, sol, order=1, solve_for_func=False)[0]


def test_1st_exact1():
    # Type: Exact differential equation, p(x,f) + q(x,f)*f' == 0,
    # where dp/df == dq/dx
    eq1 = sin(x)*cos(f(x)) + cos(x)*sin(f(x))*f(x).diff(x)
    eq2 = (2*x*f(x) + 1)/f(x) + (f(x) - x)/f(x)**2*f(x).diff(x)
    eq3 = 2*x + f(x)*cos(x) + (2*f(x) + sin(x) - sin(f(x)))*f(x).diff(x)
    eq4 = cos(f(x)) - (x*sin(f(x)) - f(x)**2)*f(x).diff(x)
    eq5 = 2*x*f(x) + (x**2 + f(x)**2)*f(x).diff(x)
    sol1 = [Eq(f(x), -acos(C1/cos(x)) + 2*pi), Eq(f(x), acos(C1/cos(x)))]
    sol2 = Eq(f(x), C1*exp(-x**2 + LambertW(C2*x*exp(x**2))))
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
    # issue 1981 needs to be addressed to test these
    # assert checkodesol(eq2, sol2, order=1, solve_for_func=False)[0]
    assert checkodesol(eq2, sol2b, order=1, solve_for_func=False)[0]
    assert checkodesol(eq3, sol3, order=1, solve_for_func=False)[0]
    assert checkodesol(eq4, sol4, order=1, solve_for_func=False)[0]
    assert checkodesol(eq5, sol5, order=1, solve_for_func=False)[0]


@slow
@XFAIL
def test_1st_exact2():
    """
    This is an exact equation that fails under the exact engine. It is caught
    by first order homogeneous albeit with a much contorted solution.  The
    exact engine fails because of a poorly simplified integral of q(0,y)dy,
    where q is the function multiplying f'.  The solutions should be
    Eq(sqrt(x**2+f(x)**2)**3+y**3, C1).  The equation below is
    equivalent, but it is so complex that checkodesol fails, and takes a long
    time to do so.
    """
    eq = (x*sqrt(x**2 + f(x)**2) - (x**2*f(x)/(f(x) -
          sqrt(x**2 + f(x)**2)))*f(x).diff(x))
    sol = dsolve(eq)
    assert sol == Eq(log(x),
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
    sol1 = Eq(f(x), C1*exp(x))
    sol2 = Eq(f(x), C1*x)
    sol3 = Eq(f(x), C1 + cos(x))
    sol4 = Eq(atan(f(x)), C1 + atan(x))
    sol5 = Eq(f(x), -2 + C1*sqrt(1 + tan(x)**2))
    #sol5 = Eq(f(x), C1*(C2 + sqrt(1 + tan(x)**2)))
    #sol5 = Eq(-log(2 + f(x)), C1 - log(1 + tan(x)**2)/2)
    assert dsolve(eq1, hint='separable') == sol1
    assert dsolve(eq2, hint='separable') == sol2
    assert dsolve(eq3, hint='separable') == sol3
    assert dsolve(eq4, hint='separable', simplify=False) == sol4
    assert dsolve(eq5, hint='separable') == simplify(sol5).expand()
    assert checkodesol(eq1, sol1, order=1, solve_for_func=False)[0]
    assert checkodesol(eq2, sol2, order=1, solve_for_func=False)[0]
    assert checkodesol(eq3, sol3, order=1, solve_for_func=False)[0]
    assert checkodesol(eq4, sol4, order=1, solve_for_func=False)[0]
    assert checkodesol(eq5, sol5, order=1, solve_for_func=False)[0]


def test_separable2():
    a = Symbol('a')
    eq6 = f(x)*x**2*f(x).diff(x) - f(x)**3 - 2*x**2*f(x).diff(x)
    eq7 = f(x)**2 - 1 - (2*f(x) + x*f(x))*f(x).diff(x)
    eq8 = x*log(x)*f(x).diff(x) + sqrt(1 + f(x)**2)
    eq9 = exp(x + 1)*tan(f(x)) + cos(f(x))*f(x).diff(x)
    eq10 = (x*cos(f(x)) + x**2*sin(f(x))*f(x).diff(x) -
            a**2*sin(f(x))*f(x).diff(x))
    # solve() messes this one up a little bit, so lets test _Integral here
    # We have to test strings with _Integral because y is a dummy variable.
    sol6str = ("Integral((_y - 2)/_y**3, (_y, f(x))) "
               "== C1 + Integral(x**(-2), x)")
    sol7 = Eq(-log(-1 + f(x)**2)/2, C1 - log(2 + x))
    sol8 = Eq(asinh(f(x)), C1 - log(log(x)))
    # integrate cannot handle the integral on the lhs (cos/tan)
    sol9str = ("Integral(cos(_y)/tan(_y), (_y, f(x)))"
               " == C1 + Integral(-E*exp(x), x)")
    sol10 = Eq(-log(-1 + sin(f(x))**2)/2, C1 - log(x**2 - a**2)/2)
    assert str(dsolve(eq6, hint='separable_Integral')) == sol6str
    assert dsolve(eq7, hint='separable', simplify=False) == sol7
    assert dsolve(eq8, hint='separable', simplify=False) == sol8
    assert str(dsolve(eq9, hint='separable_Integral')) == sol9str
    assert dsolve(eq10, hint='separable', simplify=False) == sol10
    assert checkodesol(eq7, sol7, order=1, solve_for_func=False)[0]
    assert checkodesol(eq8, sol8, order=1, solve_for_func=False)[0]
    assert checkodesol(eq10, sol10, order=1, solve_for_func=False)[0]


def test_separable3():
    eq11 = f(x).diff(x) - f(x)*tan(x)
    eq12 = (x - 1)*cos(f(x))*f(x).diff(x) - 2*x*sin(f(x))
    eq13 = f(x).diff(x) - f(x)*log(f(x))/tan(x)
    sol11 = Eq(f(x), C1*sqrt(1 + tan(x)**2))
    sol12 = Eq(log(-1 + cos(f(x))**2)/2, C1 + 2*x + 2*log(x - 1))
    sol13 = Eq(log(log(f(x))), C1 + log(cos(x)**2 - 1)/2)
    assert dsolve(eq11, hint='separable') == simplify(sol11)
    assert dsolve(eq12, hint='separable', simplify=False) == sol12
    assert dsolve(eq13, hint='separable', simplify=False) == sol13
    assert checkodesol(eq11, sol11, order=1, solve_for_func=False)[0]
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
    sol18 = Eq(-log(-1 + sin(2*f(x))**2)/4, C1 + log(-1 + sin(x)**2)/2)
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
    sol3 = Eq(f(x), C1*LambertW(C2*x))  # Eq(f(x), x*exp(-LambertW(C1*x) + 1))
    sol4 = Eq(log(f(x)), C1 - 2*exp(x/f(x)))
    sol5 = Eq(f(x), C1*exp(LambertW(C2*x**4)/2)/x)
    sol6 = Eq(log(x),
        C1 + exp(-f(x)/x)*sin(f(x)/x)/2 + exp(-f(x)/x)*cos(f(x)/x)/2)
    sol7 = Eq(log(f(x)), C1 - 2*sqrt(-x/f(x) + 1))
    sol8 = Eq(log(x), C1 - log(sqrt(1 + f(x)**2/x**2)) + atan(f(x)/x))
    assert dsolve(eq1, hint='1st_homogeneous_coeff_subs_dep_div_indep') == \
        sol1
    # indep_div_dep actually has a simpler solution for eq2,
    # but it runs too slow
    assert dsolve(eq2, hint='1st_homogeneous_coeff_subs_dep_div_indep',
            simplify=False) == sol2
    assert dsolve(eq3, hint='1st_homogeneous_coeff_best') == sol3
    assert dsolve(eq4, hint='1st_homogeneous_coeff_best') == sol4
    assert dsolve(eq5, hint='1st_homogeneous_coeff_best') == sol5
    assert dsolve(eq6, hint='1st_homogeneous_coeff_subs_dep_div_indep') == \
        sol6
    assert dsolve(eq7, hint='1st_homogeneous_coeff_best') == sol7
    assert dsolve(eq8, hint='1st_homogeneous_coeff_best') == sol8
    # checks are below


@slow
def test_1st_homogeneous_coeff_ode_check134568():
    # These are the checkodesols from test_homogeneous_coeff_ode1.
    eq1 = f(x)/x*cos(f(x)/x) - (x/f(x)*sin(f(x)/x) + cos(f(x)/x))*f(x).diff(x)
    eq3 = f(x) + (x*log(f(x)/x) - 2*x)*diff(f(x), x)
    eq4 = 2*f(x)*exp(x/f(x)) + f(x)*f(x).diff(x) - 2*x*exp(x/f(x))*f(x).diff(x)
    eq5 = 2*x**2*f(x) + f(x)**3 + (x*f(x)**2 - 2*x**3)*f(x).diff(x)
    eq6 = x*exp(f(x)/x) - f(x)*sin(f(x)/x) + x*sin(f(x)/x)*f(x).diff(x)
    eq8 = x + f(x) - (x - f(x))*f(x).diff(x)
    sol1 = Eq(f(x)*sin(f(x)/x), C1)
    sol4 = Eq(log(C1*f(x)) + 2*exp(x/f(x)), 0)
    sol3 = Eq(-f(x)/(1 + log(x/f(x))), C1)
    sol5 = Eq(log(C1*x*sqrt(1/x)*sqrt(f(x))) + x**2/(2*f(x)**2), 0)
    sol6 = Eq(-exp(-f(x)/x)*sin(f(x)/x)/2 + log(C1*x) -
            cos(f(x)/x)*exp(-f(x)/x)/2, 0)
    sol8 = Eq(-atan(f(x)/x) + log(C1*x*sqrt(1 + f(x)**2/x**2)), 0)
    assert checkodesol(eq1, sol1, order=1, solve_for_func=False)[0]
    assert checkodesol(eq3, sol3, order=1, solve_for_func=False)[0]
    assert checkodesol(eq4, sol4, order=1, solve_for_func=False)[0]
    assert checkodesol(eq5, sol5, order=1, solve_for_func=False)[0]
    assert checkodesol(eq6, sol6, order=1, solve_for_func=False)[0]
    assert checkodesol(eq8, sol8, order=1, solve_for_func=False)[0]


def test_1st_homogeneous_coeff_ode_check2():
    eq2 = x*f(x).diff(x) - f(x) - x*sin(f(x)/x)
    sol2 = Eq(x/tan(f(x)/(2*x)), C1)
    assert checkodesol(eq2, sol2, order=1, solve_for_func=False)[0]


@XFAIL
def test_1st_homogeneous_coeff_ode_check3():
    skip('This is a known issue.')
    # checker cannot determine that the following expression is zero:
    # (False,
    #   x*(log(exp(-LambertW(C1*x))) +
    #   LambertW(C1*x))*exp(-LambertW(C1*x) + 1))
    # This is blocked by issue 1981.
    eq3 = f(x) + (x*log(f(x)/x) - 2*x)*diff(f(x), x)
    sol3a = Eq(f(x), x*exp(1 - LambertW(C1*x)))
    assert checkodesol(eq3, sol3a, solve_for_func=True)[0]
    # Checker can't verify this form either
    # (False,
    #   C1*(log(C1*LambertW(C2*x)/x) + LambertW(C2*x) - 1)*LambertW(C2*x))
    # It is because a = W(a)*exp(W(a)), so log(a) == log(W(a)) + W(a) and C2 =
    # -E/C1 (which can be verified by solving with simplify=False).
    sol3b = Eq(f(x), C1*LambertW(C2*x))
    assert checkodesol(eq3, sol3b, solve_for_func=True)[0]


def test_1st_homogeneous_coeff_ode_check7():
    eq7 = (x + sqrt(f(x)**2 - x*f(x)))*f(x).diff(x) - f(x)
    sol7 = Eq(log(C1*f(x)) + 2*sqrt(1 - x/f(x)), 0)
    assert checkodesol(eq7, sol7, order=1, solve_for_func=False)[0]


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
    assert checkodesol(eq1, sol1, order=1, solve_for_func=False)[0]
    assert checkodesol(eq2, sol2, order=1, solve_for_func=False)[0]
    # test for eq3 is in test_1st_homogeneous_coeff_ode2_check3 below


def test_1st_homogeneous_coeff_ode2_check3():
    eq3 = x*exp(f(x)/x) + f(x) - x*f(x).diff(x)
    sol3 = Eq(f(x), log(log(C1/x)**(-x)))
    assert checkodesol(eq3, sol3, order=1, solve_for_func=False)[0]


def test_1st_homogeneous_coeff_ode_check9():
    _u2 = Dummy('u2')
    __a = Dummy('a')
    eq9 = f(x)**2 + (x*sqrt(f(x)**2 - x**2) - x*f(x))*f(x).diff(x)
    sol9 = Eq(-Integral(-1/(-(1 - sqrt(1 - _u2**2))*_u2 + _u2), (_u2, __a,
        x/f(x))) + log(C1*f(x)), 0)
    assert checkodesol(eq9, sol9, order=1, solve_for_func=False)[0]


def test_1st_homogeneous_coeff_ode3():
    # The standard integration engine cannot handle one of the integrals
    # involved (see issue 1452).  meijerg code comes up with an answer, but in
    # unconventional form.
    # checkodesol fails for this equation, so its test is in
    # test_1st_homogeneous_coeff_ode_check9 above. It has to compare string
    # expressions because u2 is a dummy variable.
    eq = f(x)**2 + (x*sqrt(f(x)**2 - x**2) - x*f(x))*f(x).diff(x)
    sol = Eq(log(f(x)), C1 - Piecewise(
            (-acosh(f(x)/x), abs(f(x)**2)/x**2 > 1),
            (I*asin(f(x)/x), True)))
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
    sol1 = Eq(f(x), C1 + C2*exp(-2*x))
    sol2 = Eq(f(x), (C1*exp(x) + C2*exp(2*x)))
    sol3 = Eq(f(x), C1*exp(x) + C2*exp(-x))
    sol4 = Eq(f(x), C1 + C2*exp(-3*x) + C3*exp(2*x))
    sol5 = Eq(f(x), C1*exp(x/2) + C2*exp(4*x/3))
    sol6 = Eq(f(x), C1*exp(x*(-1 + sqrt(2))) + C2*exp(x*(-sqrt(2) - 1)))
    sol7 = Eq(f(x),
        C1*exp(3*x) + C2*exp(x*(-2 - sqrt(2))) + C3*exp(x*(-2 + sqrt(2))))
    sol8 = Eq(f(x), C1 + C2*exp(x) + C3*exp(-2*x) + C4*exp(2*x))
    sol9 = Eq(f(x),
        C1*exp(x) + C2*exp(-x) + C3*exp(x*(-2 + sqrt(2))) +
        C4*exp(x*(-2 - sqrt(2))))
    sol10 = Eq(f(x),
        C1*sin(x*sqrt(a)) + C2*cos(x*sqrt(a)) + C3*exp(x*sqrt(a)) +
        C4*exp(-x*sqrt(a)))
    sol11 = Eq(f(x),
        C1*exp(x*(k - sqrt(k**2 + 2))) + C2*exp(x*(k + sqrt(k**2 + 2))))
    sol12 = Eq(f(x),
        C1*exp(2*x*(-2*abs(k) - k)) + C2*exp(2*x*(2*abs(k) - k)))
    sol13 = Eq(f(x), C1 + C2*x + C3*x**2 + C4*x**3)
    sol14 = Eq(f(x), (C1 + C2*x)*exp(-2*x))
    sol15 = Eq(f(x), (C1 + C2*x)*exp(-x) + C3*exp(x/3))
    sol16 = Eq(f(x), (C1 + C2*x + C3*x**2)*exp(2*x))
    sol17 = Eq(f(x), (C1 + C2*x)*exp(a*x))
    sol18 = Eq(f(x), C1 + C2*x + C3*x**2 + C4*exp(-3*x))
    sol19 = Eq(f(x), C1 + C2*x + C3*exp(x*sqrt(2)) + C4*exp(-x*sqrt(2)))
    sol20 = Eq(f(x), (C1 + C2*x)*exp(-3*x) + (C3 + C4*x)*exp(2*x))
    sol21 = Eq(f(x), C1*exp(x/2) + C2*exp(-x) + C3*exp(-x/3) + C4*exp(5*x/6))
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
    sol1s = constant_renumber(sol1, 'C', 1, 2)
    sol2s = constant_renumber(sol2, 'C', 1, 2)
    sol3s = constant_renumber(sol3, 'C', 1, 2)
    sol4s = constant_renumber(sol4, 'C', 1, 3)
    sol5s = constant_renumber(sol5, 'C', 1, 2)
    sol6s = constant_renumber(sol6, 'C', 1, 2)
    sol7s = constant_renumber(sol7, 'C', 1, 3)
    sol8s = constant_renumber(sol8, 'C', 1, 4)
    sol9s = constant_renumber(sol9, 'C', 1, 4)
    sol10s = constant_renumber(sol10, 'C', 1, 4)
    sol11s = constant_renumber(sol11, 'C', 1, 2)
    sol12s = constant_renumber(sol12, 'C', 1, 2)
    sol13s = constant_renumber(sol13, 'C', 1, 4)
    sol14s = constant_renumber(sol14, 'C', 1, 2)
    sol15s = constant_renumber(sol15, 'C', 1, 3)
    sol16s = constant_renumber(sol16, 'C', 1, 3)
    sol17s = constant_renumber(sol17, 'C', 1, 2)
    sol18s = constant_renumber(sol18, 'C', 1, 4)
    sol19s = constant_renumber(sol19, 'C', 1, 4)
    sol20s = constant_renumber(sol20, 'C', 1, 4)
    sol21s = constant_renumber(sol21, 'C', 1, 4)
    sol22s = constant_renumber(sol22, 'C', 1, 4)
    sol23s = constant_renumber(sol23, 'C', 1, 2)
    sol24s = constant_renumber(sol24, 'C', 1, 2)
    sol25s = constant_renumber(sol25, 'C', 1, 4)
    sol26s = constant_renumber(sol26, 'C', 1, 2)
    sol27s = constant_renumber(sol27, 'C', 1, 4)
    sol28s = constant_renumber(sol28, 'C', 1, 3)
    sol29s = constant_renumber(sol29, 'C', 1, 4)
    sol30s = constant_renumber(sol30, 'C', 1, 5)
    assert dsolve(eq1) in (sol1, sol1s)
    assert dsolve(eq2) in (sol2, sol2s)
    assert dsolve(eq3) in (sol3, sol3s)
    assert dsolve(eq4) in (sol4, sol4s)
    assert dsolve(eq5) in (sol5, sol5s)
    assert dsolve(eq6) in (sol6, sol6s)
    assert dsolve(eq7) in (sol7, sol7s)
    assert dsolve(eq8) in (sol8, sol8s)
    assert dsolve(eq9) in (sol9, sol9s)
    assert dsolve(eq10) in (sol10, sol10s)
    assert dsolve(eq11) in (sol11, sol11s)
    assert dsolve(eq12) in (sol12, sol12s)
    assert dsolve(eq13) in (sol13, sol13s)
    assert dsolve(eq14) in (sol14, sol14s)
    assert dsolve(eq15) in (sol15, sol15s)
    assert dsolve(eq16) in (sol16, sol16s)
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


def test_nth_linear_constant_coeff_homogeneous_RootOf():
    eq = f(x).diff(x, 5) + 11*f(x).diff(x) - 2*f(x)
    sol = Eq(f(x),
        C1*exp(x*RootOf(x**5 + 11*x - 2, 0)) +
        C2*exp(x*RootOf(x**5 + 11*x - 2, 1)) +
        C3*exp(x*RootOf(x**5 + 11*x - 2, 2)) +
        C4*exp(x*RootOf(x**5 + 11*x - 2, 3)) +
        C5*exp(x*RootOf(x**5 + 11*x - 2, 4)))
    assert dsolve(eq) == sol


@XFAIL
def test_nth_linear_constant_coeff_homogeneous_RootOf_sol():
    eq = f(x).diff(x, 5) + 11*f(x).diff(x) - 2*f(x)
    sol = Eq(f(x),
        C1*exp(x*RootOf(x**5 + 11*x - 2, 0)) +
        C2*exp(x*RootOf(x**5 + 11*x - 2, 1)) +
        C3*exp(x*RootOf(x**5 + 11*x - 2, 2)) +
        C4*exp(x*RootOf(x**5 + 11*x - 2, 3)) +
        C5*exp(x*RootOf(x**5 + 11*x - 2, 4)))
    assert checkodesol(eq, sol, order=5, solve_for_func=False)[0]


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
        'test': True, 'trialset': set([x**2*cos(x)*exp(x), x, cos(x), S(1),
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
        {'test': True, 'trialset': set([S(1)])}
    assert _undetermined_coefficients_match(12*exp(x), x) == \
        {'test': True, 'trialset': set([exp(x)])}
    assert _undetermined_coefficients_match(exp(I*x), x) == \
        {'test': True, 'trialset': set([exp(I*x)])}
    assert _undetermined_coefficients_match(sin(x), x) == \
        {'test': True, 'trialset': set([cos(x), sin(x)])}
    assert _undetermined_coefficients_match(cos(x), x) == \
        {'test': True, 'trialset': set([cos(x), sin(x)])}
    assert _undetermined_coefficients_match(8 + 6*exp(x) + 2*sin(x), x) == \
        {'test': True, 'trialset': set([S(1), cos(x), sin(x), exp(x)])}
    assert _undetermined_coefficients_match(x**2, x) == \
        {'test': True, 'trialset': set([S(1), x, x**2])}
    assert _undetermined_coefficients_match(9*x*exp(x) + exp(-x), x) == \
        {'test': True, 'trialset': set([x*exp(x), exp(x), exp(-x)])}
    assert _undetermined_coefficients_match(2*exp(2*x)*sin(x), x) == \
        {'test': True, 'trialset': set([exp(2*x)*sin(x), cos(x)*exp(2*x)])}
    assert _undetermined_coefficients_match(x - sin(x), x) == \
        {'test': True, 'trialset': set([S(1), x, cos(x), sin(x)])}
    assert _undetermined_coefficients_match(x**2 + 2*x, x) == \
        {'test': True, 'trialset': set([S(1), x, x**2])}
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
        {'test': True, 'trialset': set([S(1), x, x**2, exp(-2*x)])}
    assert _undetermined_coefficients_match(x*exp(-x), x) == \
        {'test': True, 'trialset': set([x*exp(-x), exp(-x)])}
    assert _undetermined_coefficients_match(x + exp(2*x), x) == \
        {'test': True, 'trialset': set([S(1), x, exp(2*x)])}
    assert _undetermined_coefficients_match(sin(x) + exp(-x), x) == \
        {'test': True, 'trialset': set([cos(x), sin(x), exp(-x)])}
    assert _undetermined_coefficients_match(exp(x), x) == \
        {'test': True, 'trialset': set([exp(x)])}
    # converted from sin(x)**2
    assert _undetermined_coefficients_match(S(1)/2 - cos(2*x)/2, x) == \
        {'test': True, 'trialset': set([S(1), cos(2*x), sin(2*x)])}
    # converted from exp(2*x)*sin(x)**2
    assert _undetermined_coefficients_match(
        exp(2*x)*(S(1)/2 + cos(2*x)/2), x
    ) == {
        'test': True, 'trialset': set([exp(2*x)*sin(2*x), cos(2*x)*exp(2*x),
        exp(2*x)])}
    assert _undetermined_coefficients_match(2*x + sin(x) + cos(x), x) == \
        {'test': True, 'trialset': set([S(1), x, cos(x), sin(x)])}
    # converted from sin(2*x)*sin(x)
    assert _undetermined_coefficients_match(cos(x)/2 - cos(3*x)/2, x) == \
        {'test': True, 'trialset': set([cos(x), cos(3*x), sin(x), sin(3*x)])}
    assert _undetermined_coefficients_match(cos(x**2), x) == {'test': False}
    assert _undetermined_coefficients_match(2**(x**2), x) == {'test': False}


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
    eq24 = f2 + f(x) - S(1)/2 - cos(2*x)/2
    # exp(2*x)*sin(x)**2
    eq25 = f(x).diff(x, 3) - f(x).diff(x) - exp(2*x)*(S(1)/2 - cos(2*x)/2)
    eq26 = (f(x).diff(x, 5) + 2*f(x).diff(x, 3) + f(x).diff(x) - 2*x -
        sin(x) - cos(x))
    # sin(2*x)*sin(x), skip 3127 for now, match bug
    eq27 = f2 + f(x) - cos(x)/2 + cos(3*x)/2
    eq28 = f(x).diff(x) - 1
    sol1 = Eq(f(x),
        -1 - x + (C1 + C2*x - 3*x**2/32 - x**3/24)*exp(-x) + C3*exp(x/3))
    sol2 = Eq(f(x), -1 - x + (C1 + C2*x - x**2/8)*exp(-x) + C3*exp(x/3))
    sol3 = Eq(f(x), 2 + C1*exp(-x) + C2*exp(-2*x))
    sol4 = Eq(f(x), 2*exp(x) + C1*exp(-x) + C2*exp(-2*x))
    sol5 = Eq(f(x), C1*exp(-x) + C2*exp(-2*x) + (S(1)/10 - 3*I/10)*exp(I*x))
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
    sol17 = Eq(f(x), (C1 + C2*x + x**4/12)*exp(-x))
    sol18 = Eq(f(x), (C1 + C2*x + C3*x**2 - x**5/60 + x**3/3)*exp(-x))
    sol19 = Eq(f(x), S(7)/4 - 3*x/2 + x**2/2 + C1*exp(-x) + (C2 - x)*exp(-2*x))
    sol20 = Eq(f(x), C1*exp(x) + C2*exp(2*x) + (6*x + 5)*exp(-x)/36)
    sol21 = Eq(f(x), -S(1)/36 - x/6 + C1*exp(-3*x) + (C2 + x/5)*exp(2*x))
    sol22 = Eq(f(x), C1*sin(x) + (C2 - x/2)*cos(x) + exp(-x)/2)
    sol23 = Eq(f(x), (C1 + C2*x + C3*x**2 + x**3/6)*exp(x))
    sol24 = Eq(f(x), S(1)/2 - cos(2*x)/6 + C1*sin(x) + C2*cos(x))
    sol25 = Eq(f(x), C1 + C2*exp(-x) + C3*exp(x) +
               (-21*sin(2*x) + 27*cos(2*x) + 130)*exp(2*x)/1560)
    sol26 = Eq(f(x),
        C1 + (C2 + C3*x - x**2/8)*sin(x) + (C4 + C5*x + x**2/8)*cos(x) + x**2)
    sol27 = Eq(f(x), cos(3*x)/16 + C1*cos(x) + (C2 + x/4)*sin(x))
    sol28 = Eq(f(x), C1 + x)
    sol1s = constant_renumber(sol1, 'C', 1, 3)
    sol2s = constant_renumber(sol2, 'C', 1, 3)
    sol3s = constant_renumber(sol3, 'C', 1, 2)
    sol4s = constant_renumber(sol4, 'C', 1, 2)
    sol5s = constant_renumber(sol5, 'C', 1, 2)
    sol6s = constant_renumber(sol6, 'C', 1, 2)
    sol7s = constant_renumber(sol7, 'C', 1, 2)
    sol8s = constant_renumber(sol8, 'C', 1, 2)
    sol9s = constant_renumber(sol9, 'C', 1, 2)
    sol10s = constant_renumber(sol10, 'C', 1, 2)
    sol11s = constant_renumber(sol11, 'C', 1, 2)
    sol12s = constant_renumber(sol12, 'C', 1, 2)
    sol13s = constant_renumber(sol13, 'C', 1, 4)
    sol14s = constant_renumber(sol14, 'C', 1, 2)
    sol15s = constant_renumber(sol15, 'C', 1, 2)
    sol16s = constant_renumber(sol16, 'C', 1, 2)
    sol17s = constant_renumber(sol17, 'C', 1, 2)
    sol18s = constant_renumber(sol18, 'C', 1, 3)
    sol19s = constant_renumber(sol19, 'C', 1, 2)
    sol20s = constant_renumber(sol20, 'C', 1, 2)
    sol21s = constant_renumber(sol21, 'C', 1, 2)
    sol22s = constant_renumber(sol22, 'C', 1, 2)
    sol23s = constant_renumber(sol23, 'C', 1, 3)
    sol24s = constant_renumber(sol24, 'C', 1, 2)
    sol25s = constant_renumber(sol25, 'C', 1, 3)
    sol26s = constant_renumber(sol26, 'C', 1, 5)
    sol27s = constant_renumber(sol27, 'C', 1, 2)
    assert dsolve(eq1, hint=hint) in (sol1, sol1s)
    assert dsolve(eq2, hint=hint) in (sol2, sol2s)
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
    assert dsolve(eq17, hint=hint) in (sol17, sol17s)
    assert dsolve(eq18, hint=hint) in (sol18, sol18s)
    assert dsolve(eq19, hint=hint) in (sol19, sol19s)
    assert dsolve(eq20, hint=hint) in (sol20, sol20s)
    assert dsolve(eq21, hint=hint) in (sol21, sol21s)
    assert dsolve(eq22, hint=hint) in (sol22, sol22s)
    assert dsolve(eq23, hint=hint) in (sol23, sol23s)
    assert dsolve(eq24, hint=hint) in (sol24, sol24s)
    assert dsolve(eq25, hint=hint) in (sol25, sol25s)
    assert dsolve(eq26, hint=hint) in (sol26, sol26s)
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


@XFAIL
def test_nth_linear_constant_coeff_undetermined_coefficients_imaginary_exp():
    # Equivalent to eq26 in
    # test_nth_linear_constant_coeff_undetermined_coefficients above.
    # This fails because the algorithm for undetermined coefficients
    # doesn't know to multiply exp(I*x) by sufficient x because it is linearly
    # dependent on sin(x) and cos(x).
    hint = 'nth_linear_constant_coeff_undetermined_coefficients'
    eq26a = f(x).diff(x, 5) + 2*f(x).diff(x, 3) + f(x).diff(x) - 2*x - exp(I*x)
    sol26 = Eq(f(x),
        C1 + (C2 + C3*x - x**2/8)*sin(x) + (C4 + C5*x + x**2/8)*cos(x) + x**2)
    assert dsolve(eq26a, hint=hint) == sol26
    assert checkodesol(eq26a, sol26, order=5, solve_for_func=False)[0]


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
    sol1 = Eq(f(x),
        -1 - x + (C1 + C2*x - 3*x**2/32 - x**3/24)*exp(-x) + C3*exp(x/3))
    sol2 = Eq(f(x), -1 - x + (C1 + C2*x - x**2/8)*exp(-x) + C3*exp(x/3))
    sol3 = Eq(f(x), C1 + x)
    sol4 = Eq(f(x), 2 + C1*exp(-x) + C2*exp(-2*x))
    sol5 = Eq(f(x), 2*exp(x) + C1*exp(-x) + C2*exp(-2*x))
    sol6 = Eq(f(x), -x*exp(x) - 2*exp(-x) + C1*exp(-2*x) + C2*exp(4*x))
    sol7 = Eq(f(x), (C1 + C2*x + x**4/12)*exp(-x))
    sol8 = Eq(f(x), C1*exp(x) + C2*exp(2*x) + (6*x + 5)*exp(-x)/36)
    sol9 = Eq(f(x), (C1 + C2*x + C3*x**2 + x**3/6)*exp(x))
    sol10 = Eq(f(x), (C1 + x*(C2 + log(x)))*exp(-x))
    sol11 = Eq(f(x), cos(x)*(C2 - Integral(1/cos(x), x)) + sin(x)*(C1 +
        Integral(1/sin(x), x)))
    sol12 = Eq(f(x), C1 + C2*x + x**3*(C3 + log(x)/6) + C4*x**2)
    sol1s = constant_renumber(sol1, 'C', 1, 3)
    sol2s = constant_renumber(sol2, 'C', 1, 3)
    sol3s = constant_renumber(sol3, 'C', 1, 2)
    sol4s = constant_renumber(sol4, 'C', 1, 2)
    sol5s = constant_renumber(sol5, 'C', 1, 2)
    sol6s = constant_renumber(sol6, 'C', 1, 2)
    sol7s = constant_renumber(sol7, 'C', 1, 2)
    sol8s = constant_renumber(sol8, 'C', 1, 2)
    sol9s = constant_renumber(sol9, 'C', 1, 3)
    sol10s = constant_renumber(sol10, 'C', 1, 2)
    sol11s = constant_renumber(sol11, 'C', 1, 2)
    sol12s = constant_renumber(sol12, 'C', 1, 4)
    assert dsolve(eq1, hint=hint) in (sol1, sol1s)
    assert dsolve(eq2, hint=hint) in (sol2, sol2s)
    assert dsolve(eq3, hint=hint) in (sol3, sol3s)
    assert dsolve(eq4, hint=hint) in (sol4, sol4s)
    assert dsolve(eq5, hint=hint) in (sol5, sol5s)
    assert dsolve(eq6, hint=hint) in (sol6, sol6s)
    assert dsolve(eq7, hint=hint) in (sol7, sol7s)
    assert dsolve(eq8, hint=hint) in (sol8, sol8s)
    assert dsolve(eq9, hint=hint) in (sol9, sol9s)
    assert dsolve(eq10, hint=hint) in (sol10, sol10s)
    assert dsolve(eq11, hint=hint + '_Integral') in (sol11, sol11s)
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


def test_nth_linear_constant_coeff_variation_of_parameters_simplify_False():
    # solve_variation_of_parameters shouldn't attempt to simplify the
    # Wronskian if simplify=False.  If wronskian() ever gets good enough
    # to simplify the result itself, this test might fail.
    hint = 'nth_linear_constant_coeff_variation_of_parameters'
    assert dsolve(f(x).diff(x, 5) + 2*f(x).diff(x, 3) + f(x).diff(x) -
        2*x - exp(I*x), f(x), hint + "_Integral", simplify=False) != \
        dsolve(f(x).diff(x, 5) + 2*f(x).diff(x, 3) + f(x).diff(x) -
        2*x - exp(I*x), f(x), hint + "_Integral", simplify=True)


def test_Liouville_ODE():
    hint = 'Liouville'
    # The first part here used to be test_ODE_1() from test_solvers.py
    eq1 = diff(f(x), x)/x + diff(f(x), x, x)/2 - diff(f(x), x)**2/2
    eq1a = diff(x*exp(-f(x)), x, x)
    # compare to test_unexpanded_Liouville_ODE() below
    eq2 = (eq1*exp(-f(x))/exp(f(x))).expand()
    eq3 = diff(f(x), x, x) + 1/f(x)*(diff(f(x), x))**2 + 1/x*diff(f(x), x)
    eq4 = x*diff(f(x), x, x) + x/f(x)*diff(f(x), x)**2 + x*diff(f(x), x)
    eq5 = Eq((x*exp(f(x))).diff(x, x), 0)
    sol1 = Eq(f(x), log(x/(C1 + C2*x)))
    sol1a = Eq(C1 + C2/x - exp(-f(x)), 0)
    sol2 = sol1
    sol3 = set(
        [Eq(f(x), -sqrt(C1 + C2*log(x))), Eq(f(x), sqrt(C1 + C2*log(x)))])
    sol4 = set([Eq(f(x), sqrt(C1 + C2*exp(x))*exp(-x/2)),
                Eq(f(x), -sqrt(C1 + C2*exp(x))*exp(-x/2))])
    sol5 = Eq(f(x), log(C1 + C2/x))
    sol1s = constant_renumber(sol1, 'C', 1, 2)
    sol2s = constant_renumber(sol2, 'C', 1, 2)
    sol3s = constant_renumber(sol3, 'C', 1, 2)
    sol4s = constant_renumber(sol4, 'C', 1, 2)
    sol5s = constant_renumber(sol5, 'C', 1, 2)
    assert dsolve(eq1, hint=hint) in (sol1, sol1s)
    assert dsolve(eq1a, hint=hint) in (sol1, sol1s)
    assert dsolve(eq2, hint=hint) in (sol2, sol2s)
    assert set(dsolve(eq3, hint=hint)) in (sol3, sol3s)
    assert set(dsolve(eq4, hint=hint)) in (sol4, sol4s)
    assert dsolve(eq5, hint=hint) in (sol5, sol5s)
    assert checkodesol(eq1, sol1, order=2, solve_for_func=False)[0]
    assert checkodesol(eq1a, sol1a, order=2, solve_for_func=False)[0]
    assert checkodesol(eq2, sol2, order=2, solve_for_func=False)[0]
    assert all(i[0] for i in checkodesol(eq3, sol3, order=2,
        solve_for_func=False))
    assert all(i[0] for i in checkodesol(eq4, sol4, order=2,
        solve_for_func=False))
    assert checkodesol(eq5, sol5, order=2, solve_for_func=False)[0]
    not_Liouville1 = classify_ode(diff(f(x), x)/x + f(x)*diff(f(x), x, x)/2 -
        diff(f(x), x)**2/2, f(x))
    not_Liouville2 = classify_ode(diff(f(x), x)/x + diff(f(x), x, x)/2 -
        x*diff(f(x), x)**2/2, f(x))
    assert hint not in not_Liouville1
    assert hint not in not_Liouville2
    assert hint + '_Integral' not in not_Liouville1
    assert hint + '_Integral' not in not_Liouville2


def test_unexpanded_Liouville_ODE():
    # This is the same as eq1 from test_Liouville_ODE() above.
    eq1 = diff(f(x), x)/x + diff(f(x), x, x)/2 - diff(f(x), x)**2/2
    eq2 = eq1*exp(-f(x))/exp(f(x))
    sol2 = Eq(f(x), log(x/(C1 + C2*x)))
    sol2s = constant_renumber(sol2, 'C', 1, 2)
    assert dsolve(eq2) in (sol2, sol2s)
    assert checkodesol(eq2, sol2, order=2, solve_for_func=False)[0]


def test_1686():
    from sympy.abc import A
    eq = x + A*(x + diff(f(x), x) + f(x)) + diff(f(x), x) + f(x) + 2
    assert classify_ode(eq, f(x)) == ('1st_linear', 'almost_linear',
        '1st_power_series', 'lie_group',
        'nth_linear_constant_coeff_undetermined_coefficients',
        'nth_linear_constant_coeff_variation_of_parameters',
        '1st_linear_Integral', 'almost_linear_Integral',
        'nth_linear_constant_coeff_variation_of_parameters_Integral')
    # 1765
    eq = (x**2 + f(x)**2)*f(x).diff(x) - 2*x*f(x)
    assert classify_ode(eq, f(x)) == ('1st_exact',
        '1st_homogeneous_coeff_best',
        '1st_homogeneous_coeff_subs_indep_div_dep',
        '1st_homogeneous_coeff_subs_dep_div_indep',
        'separable_reduced', '1st_power_series',
        'lie_group', '1st_exact_Integral',
        '1st_homogeneous_coeff_subs_indep_div_dep_Integral',
        '1st_homogeneous_coeff_subs_dep_div_indep_Integral',
        'separable_reduced_Integral')

def test_1726():
    raises(ValueError, lambda: dsolve(f(x, y).diff(x) - y*f(x, y), f(x)))
    assert classify_ode(f(x, y).diff(x) - y*f(x, y), f(x), dict=True) == \
        {'default': None, 'order': 0}
    # See also issue 694, test Z13.
    raises(ValueError, lambda: dsolve(f(x).diff(x), f(y)))
    assert classify_ode(f(x).diff(x), f(y), dict=True) == \
        {'default': None, 'order': 0}


def test_constant_renumber_order_issue2209():
    from sympy.utilities.iterables import variations

    assert constant_renumber(C1*x + C2*y, "C", 1, 2) == \
        constant_renumber(C1*y + C2*x, "C", 1, 2) == \
        C1*x + C2*y
    e = C1*(C2 + x)*(C3 + y)
    for a, b, c in variations([C1, C2, C3], 3):
        assert constant_renumber(a*(b + x)*(c + y), "C", 1, 3) == e


def test_issue_2671():
    k = Symbol("k", real=True)
    t = Symbol('t')
    w = Function('w')
    sol = dsolve(w(t).diff(t, 6) - k**6*w(t), w(t))
    assert len([s for s in sol.atoms(Symbol) if s.name.startswith('C')]) == 6
    assert constantsimp((C1*cos(x) + C2*cos(x))*exp(x), x, 2) == \
        C1*cos(x)*exp(x)
    assert constantsimp(C1*cos(x) + C2*cos(x) + C3*sin(x), x, 2) == \
        C1*cos(x) + C3*sin(x)
    assert constantsimp(exp(C1 + x), x, 1) == C1*exp(x)
    assert constantsimp(2**(C1 + x), x, 1) == C1*2**x
    assert constantsimp(2**(C1 + x), x, 1) == C1*2**x
    assert constantsimp(x + C1 + y, x, 1) == C1 + x
    assert constantsimp(x + C1 + Integral(x, (x, 1, 2)), x, 1) == C1 + x


def test_issue_2013_2331():
    assert homogeneous_order(-log(x) + acosh(x), x) is None
    assert homogeneous_order(y - log(x), x, y) is None


def test_nth_order_linear_euler_eq_homogeneous():
    x, t, a, b, c = symbols('x t a b c')
    y = Function('y')
    our_hint = "nth_linear_euler_eq_homogeneous"

    eq = diff(f(t), t, 4)*t**4 - 13*diff(f(t), t, 2)*t**2 + 36*f(t)
    assert our_hint in classify_ode(eq)

    eq = a*y(t) + b*t*diff(y(t), t) + c*t**2*diff(y(t), t, 2)
    assert our_hint in classify_ode(eq)

    eq = Eq(-3*diff(f(x), x)*x + 2*x**2*diff(f(x), x, x), 0)
    sol = C1 + C2*x**Rational(5, 2)
    sols = constant_renumber(sol, 'C', 1, 3)
    assert our_hint in classify_ode(eq)
    assert dsolve(eq, f(x), hint=our_hint).rhs in (sol, sols)
    assert checkodesol(eq, sol, order=2, solve_for_func=False)[0]

    eq = Eq(3*f(x) - 5*diff(f(x), x)*x + 2*x**2*diff(f(x), x, x), 0)
    sol = C1*sqrt(x) + C2*x**3
    sols = constant_renumber(sol, 'C', 1, 3)
    assert our_hint in classify_ode(eq)
    assert dsolve(eq, f(x), hint=our_hint).rhs in (sol, sols)
    assert checkodesol(eq, sol, order=2, solve_for_func=False)[0]

    eq = Eq(4*f(x) + 5*diff(f(x), x)*x + x**2*diff(f(x), x, x), 0)
    sol = (C1 + C2*log(x))/x**2
    sols = constant_renumber(sol, 'C', 1, 3)
    assert our_hint in classify_ode(eq)
    assert dsolve(eq, f(x), hint=our_hint).rhs in (sol, sols)
    assert checkodesol(eq, sol, order=2, solve_for_func=False)[0]

    eq = Eq(6*f(x) - 6*diff(f(x), x)*x + 1*x**2*diff(f(x), x, x) + x**3*diff(f(x), x, x, x), 0)
    sol = dsolve(eq, f(x), hint=our_hint)
    sol = C1/x**2 + C2*x + C3*x**3
    sols = constant_renumber(sol, 'C', 1, 4)
    assert our_hint in classify_ode(eq)
    assert dsolve(eq, f(x), hint=our_hint).rhs in (sol, sols)
    assert checkodesol(eq, sol, order=2, solve_for_func=False)[0]

    eq = Eq(-125*f(x) + 61*diff(f(x), x)*x - 12*x**2*diff(f(x), x, x) + x**3*diff(f(x), x, x, x), 0)
    sol = x**5*(C1 + C2*log(x) + C3*log(x)**2)
    sols = constant_renumber(sol, 'C', 1, 4)
    assert our_hint in classify_ode(eq)
    assert dsolve(eq, f(x), hint=our_hint).rhs in (sol, sols)
    assert checkodesol(eq, sol, order=2, solve_for_func=False)[0]

    eq = t**2*diff(y(t), t, 2) + t*diff(y(t), t) - 9*y(t)
    sol = C1*t**3 + C2*t**-3
    sols = constant_renumber(sol, 'C', 1, 3)
    assert our_hint in classify_ode(eq)
    assert dsolve(eq, y(t), hint=our_hint).rhs in (sol, sols)
    assert checkodesol(eq, sol, order=2, solve_for_func=False)[0]


def test_issue_1996():
    f = Function('f')
    raises(ValueError, lambda: dsolve(f(x).diff(x)**2, f(x), 'separable'))
    raises(ValueError, lambda: dsolve(f(x).diff(x)**2, f(x), 'fdsjf'))


def test_almost_linear():
    from sympy import Ei
    A = Symbol('A', positive=True)
    our_hint = 'almost_linear'
    f = Function('f')
    d = f(x).diff(x)
    eq = x**2*f(x)**2*d + f(x)**3 + 1
    sol = dsolve(eq, f(x), hint = 'almost_linear')
    assert sol[0].rhs == (C1*exp(3/x) - 1)**(1/3)
    assert checkodesol(eq, sol, order=1, solve_for_func=False)[0]

    eq = x*f(x)*d + 2*x*f(x)**2 + 1
    sol = dsolve(eq, f(x), hint = 'almost_linear')
    assert sol[0].rhs == -sqrt((C1 - 2*Ei(4*x))*exp(-4*x))
    assert checkodesol(eq, sol, order=1, solve_for_func=False)[0]

    eq = x*d + x*f(x) + 1
    sol = dsolve(eq, f(x), hint = 'almost_linear')
    assert sol.rhs == (C1 - Ei(x))*exp(-x)
    assert checkodesol(eq, sol, order=1, solve_for_func=False)[0]
    assert our_hint in classify_ode(eq, f(x))

    eq = x*exp(f(x))*d + exp(f(x)) + 3*x
    sol = dsolve(eq, f(x), hint = 'almost_linear')
    assert sol.rhs == log(C1/x - 3*x/2)
    assert checkodesol(eq, sol, order=1, solve_for_func=False)[0]

    eq = x + A*(x + diff(f(x), x) + f(x)) + diff(f(x), x) + f(x) + 2
    sol = dsolve(eq, f(x), hint = 'almost_linear')
    assert sol.rhs == (C1 + Piecewise(
        (x, Eq(A + 1, 0)), ((-A*x + A - x - 1)*exp(x)/(A + 1), True)))*exp(-x)
    assert checkodesol(eq, sol, order=1, solve_for_func=False)[0]


def test_exact_enhancement():
    f = Function('f')(x)
    df = Derivative(f, x)
    eq = f/x**2 + ((f*x - 1)/x)*df
    sol = dsolve(eq, f)
    rhs = [eq.rhs for eq in sol]
    assert rhs == [(-sqrt(C1*x**2 + 1) + 1)/x, (sqrt(C1*x**2 + 1) + 1)/x]

    eq = (x*f - 1) + df*(x**2 - x*f)
    rhs = [sol.rhs for sol in dsolve(eq, f)]
    assert rhs[0] == x - sqrt(C1 + x**2 - 2*log(x))
    assert rhs[1] == x + sqrt(C1 + x**2 - 2*log(x))

    eq = (x + 2)*sin(f) + df*x*cos(f)
    rhs = [sol.rhs for sol in dsolve(eq, f)]
    assert rhs == [
        -acos(-sqrt(C1*exp(-2*x) + x**4)/x**2) + 2*pi,
        -acos(sqrt(C1*exp(-2*x) + x**4)/x**2) + 2*pi,
        acos(-sqrt(C1*exp(-2*x) + x**4)/x**2),
        acos(sqrt(C1*exp(-2*x) + x**4)/x**2)]


def test_separable_reduced():
    f = Function('f')
    df = f(x).diff(x)
    eq = (x / f(x))*df  + tan(x**2*f(x) / (x**2*f(x) - 1))
    assert classify_ode(eq) == ('1st_linear', 'separable_reduced', 'lie_group',
        '1st_linear_Integral', 'separable_reduced_Integral')

    eq = x* df  + f(x)* (1 / (x**2*f(x) - 1))
    assert classify_ode(eq) == ('1st_linear', 'separable_reduced', 'lie_group',
        '1st_linear_Integral', 'separable_reduced_Integral')
    sol = dsolve(eq, hint = 'separable_reduced', simplify=False)
    assert sol.lhs ==  log(x**2*f(x))/3 + log(x**2*f(x) - S(3)/2)/6
    assert sol.rhs == C1 + log(x)
    assert checkodesol(eq, sol, order=1, solve_for_func=False)[0]

    eq = df + (f(x) / (x**4*f(x) - x))
    assert classify_ode(eq) == ('1st_linear', 'separable_reduced', 'lie_group',
        '1st_linear_Integral', 'separable_reduced_Integral')
    # generates PolynomialError in solve attempt
    sol = dsolve(eq, hint = 'separable_reduced')
    assert sol.lhs == log(x**3*f(x))/4 + log(x**3*f(x) - S(4)/3)/12
    assert sol.rhs == C1 + log(x)
    assert checkodesol(eq, sol, order=1, solve_for_func=False)[0]

    eq = x*df + f(x)*(x**2*f(x))
    sol = dsolve(eq, hint = 'separable_reduced', simplify=False)
    assert sol == Eq(log(x**2*f(x))/2 - log(x**2*f(x) - 2)/2, C1 + log(x))
    assert checkodesol(eq, sol, order=1, solve_for_func=False)[0]


def test_homogeneous_function():
    f = Function('f')
    eq1 = tan(x + f(x))
    eq2 = sin((3*x)/(4*f(x)))
    eq3 = cos(3*x/4*f(x))
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
    from sympy.solvers.ode import _linear_coeff_match
    n, d = z*(2*x + 3*f(x) + 5), z*(7*x + 9*f(x) + 11)
    rat = n/d
    eq1 = sin(rat) + cos(rat.expand())
    eq2 = rat
    eq3 = log(sin(rat))
    ans = (4, -S(13)/3)
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
    df = f(x).diff(x)
    sol = Eq(f(x), C1/(x**2 + 6*x + 9) - S(3)/2)
    # XXX if force is not used in solve, the following is returned which,
    # for C1 = -81/2, will satisfy the original equation. Should there be
    # another free symbol so a family of solutions can be obtained, e.g.
    # (C2 + C1*x + ... etc)/(...):
    # Eq(f(x), (C1 + C1*x - 3*x**3/2 - 27*x**2/2)/(x**3 + 9*x**2 + 27*x + 27))
    eq = df + (3 + 2*f(x))/(x + 3)
    assert dsolve(eq, hint='linear_coefficients') == sol
    assert checkodesol(eq, sol, order=1, solve_for_func=False)[0]


@XFAIL
def test_constantsimp_take_problem():
    # should this have C1 and C2 or C1 and exp(C1) in addition to x?
    # see note in test_linear_coefficients above
    c = exp(C1) + 2
    assert len(Poly(constantsimp(exp(C1) + c + c*x, x, 2)).gens) > 2


def test_issue_3780():
    f = Function('f')
    eq = Eq(Derivative(f(x), x, 2) - 2*Derivative(f(x), x) + f(x), sin(x))
    sol = (C1 + C2*x)*exp(x) + cos(x)/S(2)
    assert dsolve(eq).rhs == sol
    assert checkodesol(eq, sol, order=1, solve_for_func=False)[0]


def test_issue_3890():
    f = Function('f')
    k = Symbol('k')
    assert dsolve(f(x).diff(x) - x*exp(-k*x), f(x)) == \
        Eq(f(x), C1 + Piecewise(
            (x**2/2, Eq(k**3, 0)),
            ((-k**2*x - k)*exp(-k*x)/k**3, True)
        ))
    assert dsolve(-f(x).diff(x) + x*exp(-k*x), f(x)) == \
        Eq(f(x), Piecewise((C1 + x**2/2, Eq(k**3, 0)),
            (C1 - x*exp(-k*x)/k - exp(-k*x)/k**2, True)
        ))


def test_heuristic1():
    y, a, b, c, a4, a3, a2, a1, a0 = symbols("y a b c a4 a3 a2 a1 a0")
    y = Symbol('y')
    f = Function('f')
    xi = Function('xi')
    eta = Function('eta')
    df = f(x).diff(x)
    eq = Eq(df, x**2*f(x))
    eq1 = f(x).diff(x) + a*f(x) - c*exp(b*x)
    eq2 = f(x).diff(x) + 2*x*f(x) - x*exp(-x**2)
    eq3 = (1 + 2*x)*df + 2 - 4*exp(-f(x))
    eq4 = f(x).diff(x) - (a4*x**4 + a3*x**3 + a2*x**2 + a1*x + a0)**(S(-1)/2)
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
        xi(x, f(x)): sqrt(2*a0 + 2*a1*x + 2*a2*x**2 + 2*a3*x**3 + 2*a4*x**4)}]
    i5 = infinitesimals(eq5, hint='abaco1_simple')
    assert i5 == [{xi(x, f(x)): 0, eta(x, f(x)): exp(-1/x)}]

    ilist = [i, i1, i2, i3, i4, i5]
    for eq, i in (zip(eqlist, ilist)):
        check = checkinfsol(eq, i)
        assert check[0]


@XFAIL
def test_issue_3148():
    eq = x**2*f(x)**2 + x*Derivative(f(x), x)
    sol = dsolve(eq, hint = 'separable_reduced')
    assert checkodesol(eq, sol, order=1)[0]


def test_heuristic2():
    y = Symbol('y')
    xi = Function('xi')
    eta = Function('eta')
    df = f(x).diff(x)

    # This ODE can be solved by the Lie Group method, when there are
    # better assumptions
    eq = df - (f(x)/x)*(x*log(x**2/f(x)) + 2)
    i = infinitesimals(eq, hint='abaco1_product')
    assert i == [{eta(x, f(x)): f(x)*exp(-x), xi(x, f(x)): 0}]
    assert checkinfsol(eq, i)[0]


def test_heuristic3():
    y = Symbol('y')
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
    xi = Function('xi')
    eta = Function('eta')

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
    xi = Function('xi')
    eta = Function('eta')
    F = Function('F')
    a, b, m, n = symbols("a b m n")

    eq = x**(n*(m + 1) - m)*(f(x).diff(x)) - a*f(x)**n -b*x**(n*(m + 1))
    i = infinitesimals(eq, hint='linear')
    assert checkinfsol(eq, i)[0]

@XFAIL
def test_kamke():
    a, b, alpha, c = symbols("a b alpha c")
    eq = x**2*(a*f(x)**2+(f(x).diff(x))) + b*x**alpha + c
    i = infinitesimals(eq, hint='sum_function')
    assert checkinfsol(eq, i)[0]


def test_series():
    C0 = Symbol("C0")
    eq = f(x).diff(x) - f(x)
    assert dsolve(eq, hint='1st_power_series') == Eq(f(x),
        C0 + C0*x + C0*x**2/S(2) + C0*x**3/S(6) + C0*x**4/S(24) +
        C0*x**5/S(120) + O(x**6))
    eq = f(x).diff(x) - x*f(x)
    assert dsolve(eq, hint='1st_power_series') == Eq(f(x),
        C0*x**4/S(8) + C0*x**2/S(2) + C0 + O(x**6))
    eq = f(x).diff(x) - sin(x*f(x))
    assert dsolve(eq, hint='1st_power_series', ics={f(2): 2}, n=3) == Eq(
        f(x), (x - 2)**2*(2*cos(4) + 2*sin(4)*cos(4))/S(2) + (x - 2)*sin(4) +
        2 + O(x**3))


def test_lie_group():
    C1 = Symbol("C1")
    a, b, c = symbols("a b c")
    eq = f(x).diff(x)**2
    sol = dsolve(eq, f(x), hint='lie_group')
    assert checkodesol(eq, sol)[0]

    eq = Eq(f(x).diff(x), x**2*f(x))
    sol = dsolve(eq, f(x), hint='lie_group')
    assert sol == Eq(f(x), C1*exp(x**3/S(3)))
    assert checkodesol(eq, sol)[0]

    eq = f(x).diff(x) + a*f(x) - c*exp(b*x)
    sol = dsolve(eq, f(x), hint='lie_group')
    assert checkodesol(eq, sol)[0]

    eq = f(x).diff(x) + 2*x*f(x) - x*exp(-x**2)
    sol = dsolve(eq, f(x), hint='lie_group')
    assert sol == Eq(f(x), (C1 + x**2/S(2))*exp(-x**2))
    assert checkodesol(eq, sol)[0]

    eq = (1 + 2*x)*(f(x).diff(x)) + 2 - 4*exp(-f(x))
    sol = dsolve(eq, f(x), hint='lie_group')
    assert sol == Eq(f(x), log(C1/(2*x + 1) + 2))
    assert checkodesol(eq, sol)[0]

    eq = x**2*(f(x).diff(x)) - f(x) + x**2*exp(x - (1/x))
    sol = dsolve(eq, f(x), hint='lie_group')
    assert checkodesol(eq, sol)[0]

    eq = x**2*f(x)**2 + x*Derivative(f(x), x)
    sol = dsolve(eq, f(x), hint='lie_group')
    assert sol == Eq(f(x), 2/(C1 + x**2))
    assert checkodesol(eq, sol)[0]


def test_user_infinitesimals():
    C2 = Symbol("C2")
    eq = x*(f(x).diff(x)) + 1 - f(x)**2
    sol = dsolve(eq, hint='lie_group', xi=sqrt(f(x) - 1)/sqrt(f(x) + 1),
        eta=0)
    assert sol == Eq(f(x), (C2 + x**2)/(C1 - x**2))
    raises(ValueError, lambda: dsolve(eq, hint='lie_group', xi=0, eta=f(x)))

@XFAIL
def test_issue_3982():
    eq = x*(f(x).diff(x)) + 1 - f(x)**2
    assert dsolve(eq) == Eq(f(x), (C2 + x**2)/(C1 - x**2))


def test_2nd_power_series_ordinary():
    C0, C1 = symbols("C0 C1")
    eq = f(x).diff(x, 2) - x*f(x)
    assert classify_ode(eq) == ('2nd_power_series_ordinary',)
    assert dsolve(eq) == Eq(f(x),
        C0*(x**3/S(6) + 1) + C1*x*(x**3/S(12) + 1) + O(x**6))
    assert dsolve(eq, x0=-2) == Eq(f(x),
        C0*((x + 2)**4/S(6) + (x + 2)**3/S(6) - (x + 2)**2 + 1)
        + C1*(x + (x + 2)**4/S(12) - (x + 2)**3/3 + S(2))
        + O(x**6))
    assert dsolve(eq, n=2) == Eq(f(x), C1*x + C0 + O(x**2))

    eq = (1 + x**2)*(f(x).diff(x, 2)) + 2*x*(f(x).diff(x)) -2*f(x)
    assert classify_ode(eq) == ('2nd_power_series_ordinary',)
    assert dsolve(eq) == Eq(f(x), C0*(-x**4/S(3) + x**2 + 1) + C1*x
        + O(x**6))

    eq = f(x).diff(x, 2) + x*(f(x).diff(x)) + f(x)
    assert classify_ode(eq) == ('2nd_power_series_ordinary',)
    assert dsolve(eq) == Eq(f(x), C0*(
        x**4/S(8) - x**2/S(2) + 1) + C1*x*(-x**2/S(3) + 1) + O(x**6))

    eq = f(x).diff(x, 2) + f(x).diff(x) - x*f(x)
    assert classify_ode(eq) == ('2nd_power_series_ordinary',)
    assert dsolve(eq) == Eq(f(x), C0*(
        -x**4/S(24) + x**3/S(6) + 1) + C1*x*(x**3/S(24) + x**2/S(6) - x/S(2)
        + 1) + O(x**6))

    eq = f(x).diff(x, 2) + x*f(x)
    assert classify_ode(eq) == ('2nd_power_series_ordinary',)
    assert dsolve(eq, n=7) == Eq(f(x), C0*(
        x**6/S(180) - x**3/S(6) + 1) + C1*x*(-x**3/S(12) + 1) + O(x**7))


def test_2nd_power_series_ordinary():
    C0, C1 = symbols("C0 C1")
    eq = x**2*(f(x).diff(x, 2)) - 3*x*(f(x).diff(x)) + (4*x + 4)*f(x)
    assert dsolve(eq) == Eq(f(x), C0*x**2*(-16*x**3/S(9) +
        4*x**2 - 4*x + 1) + O(x**6))

    eq = 4*x**2*(f(x).diff(x, 2)) -8*x**2*(f(x).diff(x)) + (4*x**2 +
        1)*f(x)
    assert dsolve(eq) == Eq(f(x), C0*sqrt(x)*(
        x**4/S(24) + x**3/S(6) + x**2/S(2) + x + 1) + O(x**6))

    eq = x**2*(f(x).diff(x, 2)) - x**2*(f(x).diff(x)) + (
        x**2 - 2)*f(x)
    assert dsolve(eq) == Eq(f(x), C1*(-x**6/S(720) - 3*x**5/S(80) - x**4/S(8) +
        x**2/S(2) + x/S(2) + 1)/x + C0*x**2*(-x**3/S(60) + x**2/S(20) + x/S(2) + 1)
        + O(x**6))

    eq = x**2*(f(x).diff(x, 2)) + x*(f(x).diff(x)) + (x**2 - 1/S(4))*f(x)
    assert dsolve(eq) == Eq(f(x), C1*(x**4/S(24) - x**2/S(2) + 1)/sqrt(x) +
        C0*sqrt(x)*(x**4/S(120) - x**2/S(6) + 1) + O(x**6))

    eq = x*(f(x).diff(x, 2)) - f(x).diff(x) + 4*x**3*f(x)
    assert dsolve(eq) == Eq(f(x), C1*(-x**4/S(2) + 1) + C0*x**2 + O(x**6))

def test_issue3994():
    sol = Eq(f(x), C1 - 2*x*sqrt(x**3)/5)
    eq = Derivative(f(x), x)**2 - x**3
    assert dsolve(eq) == sol and checkodesol(eq, sol) == (True, 0)
