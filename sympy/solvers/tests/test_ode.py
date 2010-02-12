from sympy import Function, dsolve, Symbol, sin, cos, sinh, acos, tan, cosh, \
        I, exp, log, simplify, normal, together, ratsimp, powsimp, \
        fraction, radsimp, Eq, sqrt, pi, erf,  diff, Rational, asinh, trigsimp, \
        S, RootOf, Poly, Integral, atan, Equality, solve, O, LambertW
from sympy.abc import x, y, z
from sympy.solvers.ode import ode_order, homogeneous_order, \
        _undetermined_coefficients_match, classify_ode, checkodesol, ode_renumber
from sympy.utilities.pytest import XFAIL, skip, raises

C1 = Symbol('C1')
C2 = Symbol('C2')
C3 = Symbol('C3')
C4 = Symbol('C4')
C5 = Symbol('C5')
f = Function('f')
g = Function('g')

# Note that if the ODE solver, the integral engine, solve(), or even simplify(),
# changes, these tests could fail but still be correct, only written differently.
# Also not that in differently formatted solutions, the arbitrary constants
# might not be equal.  Using specific hints in tests can help avoid this.

# Tests order higher than 1 should run the solutions through ode_renumber
# because it will normalize it (ode_renumber causes dsolve() to return different
# results on different machines)
def test_checkodesol():
    # For the most part, checkodesol is well tested in the tests below.
    # These tests only handle cases not checked below.
    raises(ValueError, "checkodesol(f(x).diff(x), f(x), x)")
    raises(ValueError, "checkodesol(f(x).diff(x), f(x, y), Eq(f(x), x))")
    assert checkodesol(f(x).diff(x), f(x), Eq(f(x), x)) is not True
    assert checkodesol(f(x).diff(x), f(x), Eq(f(x), x)) == (False, 1)
    sol1 = Eq(f(x)**5 + 11*f(x) - 2*f(x) + x, 0)
    assert checkodesol(diff(sol1.lhs, x), f(x), sol1) == (True, 0)
    assert checkodesol(diff(sol1.lhs, x)*exp(f(x)), f(x), sol1) == (True, 0)
    assert checkodesol(diff(sol1.lhs, x, 2), f(x), sol1) == (True, 0)
    assert checkodesol(diff(sol1.lhs, x, 2)*exp(f(x)), f(x), sol1) == (True, 0)
    assert checkodesol(diff(sol1.lhs, x, 3), f(x), sol1) == (True, 0)
    assert checkodesol(diff(sol1.lhs, x, 3)*exp(f(x)), f(x), sol1) == (True, 0)
    assert checkodesol(diff(sol1.lhs, x, 3), f(x), Eq(f(x), x*log(x))) == \
        (False, 60*x**4*f(x)**3*log(x) + 180*x**5*f(x)**2*log(x) + \
        60*x**5*f(x)**2 + 180*x**5*f(x)**2*log(x)**2 - 9*x**3 + \
        60*x**4*f(x)**3 + 60*x**5*f(x)**2*log(x)**3 - 5*x**3*f(x)**4)
    assert checkodesol(diff(exp(f(x)) + x, x)*x, f(x), Eq(exp(f(x)) + x)) == (True, 0)
    assert checkodesol(diff(exp(f(x)) + x, x)*x, f(x), Eq(exp(f(x)) + x), \
        solve_for_func=False) == (True, 0)
    assert checkodesol(f(x).diff(x, 2), f(x), [Eq(f(x), C1 + C2*x), \
        Eq(f(x), C2 + C1*x), Eq(f(x), C1*x + C2*x**2)]) == \
            [(True, 0), (True, 0), (False, 2*C2)]
    assert checkodesol(f(x).diff(x) - 1/f(x)/2, f(x), Eq(f(x)**2, x)) == \
        [(True, 0), (True, 0)]
    assert checkodesol(f(x).diff(x) - f(x), f(x), Eq(C1*exp(x), f(x))) == (True, 0)

def test_dsolve_options():
    eq = x*f(x).diff(x) + f(x)
    a = dsolve(eq, f(x), hint='all')
    b = dsolve(eq, f(x), hint='all', simplify=False)
    c = dsolve(eq, f(x), hint='all_Integral')
    keys = ['1st_exact', '1st_exact_Integral', '1st_homogeneous_coeff_best', \
        '1st_homogeneous_coeff_subs_dep_div_indep', \
        '1st_homogeneous_coeff_subs_dep_div_indep_Integral', \
        '1st_homogeneous_coeff_subs_indep_div_dep', \
        '1st_homogeneous_coeff_subs_indep_div_dep_Integral', '1st_linear', \
        '1st_linear_Integral', 'best', 'best_hint', 'default', 'order', \
        'separable', 'separable_Integral']
    Integral_keys = ['1st_exact_Integral',
    '1st_homogeneous_coeff_subs_dep_div_indep_Integral',
    '1st_homogeneous_coeff_subs_indep_div_dep_Integral', '1st_linear_Integral',
    'best', 'best_hint', 'default', 'order', 'separable_Integral']
    assert sorted(a.keys()) == keys
    assert a['order'] == ode_order(eq, f(x))
    assert a['best'] == Eq(f(x), C1/x)
    assert dsolve(eq, f(x), hint='best') == Eq(f(x), C1/x)
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
    assert dsolve(eq, f(x), hint='best', simplify=False) == Eq(f(x), C1/x)
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
    raises(ValueError, "dsolve(eq, f(x), 'notarealhint')")
    raises(ValueError, "dsolve(eq, f(x), 'Liouville')")
    assert dsolve(f(x).diff(x) - 1/f(x)**2, f(x), 'all')['best'] == \
        dsolve(f(x).diff(x) - 1/f(x)**2, f(x), 'best')


def test_classify_ode():
    assert classify_ode(f(x).diff(x, 2), f(x)) == \
        ('nth_linear_constant_coeff_homogeneous', 'Liouville', 'Liouville_Integral')
    assert classify_ode(f(x), f(x)) == ()
    assert classify_ode(Eq(f(x).diff(x), 0), f(x)) == ('separable', '1st_linear', \
        '1st_homogeneous_coeff_best', '1st_homogeneous_coeff_subs_indep_div_dep', \
        '1st_homogeneous_coeff_subs_dep_div_indep', \
        'nth_linear_constant_coeff_homogeneous', 'separable_Integral', \
        '1st_linear_Integral', '1st_homogeneous_coeff_subs_indep_div_dep_Integral', \
        '1st_homogeneous_coeff_subs_dep_div_indep_Integral')
    assert classify_ode(f(x).diff(x)**2, f(x)) == ()
    # 1650: f(x) should be cleared from highest derivative before classifying
    a = classify_ode(Eq(f(x).diff(x) + f(x), x), f(x))
    b = classify_ode(f(x).diff(x)*f(x) + f(x)*f(x) - x*f(x), f(x))
    c = classify_ode(f(x).diff(x)/f(x) + f(x)/f(x) - x/f(x), f(x))
    assert a == b == c != ()
    assert classify_ode(2*x*f(x)*f(x).diff(x) + (1 + x)*f(x)**2 - exp(x), f(x)) ==\
        ('Bernoulli', 'Bernoulli_Integral')
    raises(ValueError, "classify_ode(x + f(x, y).diff(x).diff(y), f(x, y))")

def test_ode_order():
    f = Function('f')
    g = Function('g')
    x = Symbol('x')
    assert ode_order(3*x*exp(f(x)), f(x)) == 0
    assert ode_order(x*diff(f(x),x)+3*x*f(x)-sin(x)/x, f(x)) == 1
    assert ode_order(x**2*f(x).diff(x,x)+x*diff(f(x),x)-f(x),f(x)) == 2
    assert ode_order(diff(x*exp(f(x)),x,x), f(x)) == 2
    assert ode_order(diff(x*diff(x*exp(f(x)), x,x), x), f(x)) == 3
    assert ode_order(diff(f(x), x, x), g(x)) == 0
    assert ode_order(diff(f(x), x, x)*diff(g(x), x), f(x)) == 2
    assert ode_order(diff(f(x), x, x)*diff(g(x), x), g(x)) == 1
    assert ode_order(diff(x*diff(x*exp(f(x)), x,x), x), g(x)) == 0

# In all tests below, checkodesol has the order option set to prevent superfluous
# calls to ode_order(), and the solve_for_func flag set to False because
# dsolve() already tries to solve for the function, unless the simplify=False
# option is set.
def test_old_ode_tests():
    # These are simple tests from the old ode module
    eq1 = Eq(f(x).diff(x), 0)
    eq2 = Eq(3*f(x).diff(x) - 5, 0)
    eq3 = Eq(3*f(x).diff(x), 5)
    eq4 = Eq(9*f(x).diff(x, x) + f(x), 0)
    eq5 = Eq(9*f(x).diff(x, x), f(x))
    # Type: a(x)f'(x)+b(x)*f(x)+c(x)=0
    eq6 = Eq(x**2*f(x).diff(x) + 3*x*f(x) - sin(x)/x, 0)
    eq7 = Eq(f(x).diff(x,x) - 3*diff(f(x),x) + 2*f(x), 0)
    # Type: 2nd order, constant coefficients (two real different roots)
    eq8 = Eq(f(x).diff(x,x) - 4*diff(f(x),x) + 4*f(x), 0)
    # Type: 2nd order, constant coefficients (two real equal roots)
    eq9 = Eq(f(x).diff(x,x)+2*diff(f(x),x)+3*f(x), 0)
    # Type: 2nd order, constant coefficients (two complex roots)
    eq10 = Eq(3*f(x).diff(x) -1,0)
    eq11 = Eq(x*f(x).diff(x) -1,0)
    sol1 = Eq(f(x),C1)
    sol2 = Eq(f(x),C1+5*x/3)
    sol3 = Eq(f(x),C1+5*x/3)
    sol4 = Eq(f(x),C1*sin(x/3) + C2*cos(x/3))
    sol5 = Eq(f(x),C1*exp(-x/3) + C2*exp(x/3))
    sol6 = Eq(f(x),(C1-cos(x))/x**3)
    sol7list = [Eq(f(x),C1*exp(2*x) + C2*exp(x)), Eq(f(x),C1*exp(x) + C2*exp(2*x)),]
    sol8 = Eq(f(x),(C1 + C2*x)*exp(2*x))
    sol9 = Eq(f(x), (C1*sin(x*sqrt(2)) + C2*cos(x*sqrt(2)))*exp(-x))
    sol10 = Eq(f(x), C1 + x/3)
    sol11 = Eq(f(x), C1 + log(x))
    assert dsolve(eq1, f(x)) == sol1
    assert dsolve(eq1.lhs, f(x)) == sol1
    assert dsolve(eq2, f(x)) == sol2
    assert dsolve(eq3, f(x)) == sol3
    assert dsolve(eq4, f(x)) == sol4
    assert dsolve(eq5, f(x)) == sol5
    assert dsolve(eq6, f(x)) == sol6
    assert dsolve(eq7, f(x)) in sol7list
    assert dsolve(eq8, f(x)) == sol8
    assert dsolve(eq9, f(x)) == sol9
    assert dsolve(eq10, f(x)) == sol10
    assert dsolve(eq11, f(x)) == sol11
    assert checkodesol(eq1, f(x), sol1, order=1, solve_for_func=False)[0]
    assert checkodesol(eq2, f(x), sol2, order=1, solve_for_func=False)[0]
    assert checkodesol(eq3, f(x), sol3, order=1, solve_for_func=False)[0]
    assert checkodesol(eq4, f(x), sol4, order=2, solve_for_func=False)[0]
    assert checkodesol(eq5, f(x), sol5, order=2, solve_for_func=False)[0]
    assert checkodesol(eq6, f(x), sol6, order=1, solve_for_func=False)[0]
    assert checkodesol(eq7, f(x), sol7list, order=2, solve_for_func=False) == \
        [(True, 0), (True, 0)]
    assert checkodesol(eq8, f(x), sol8, order=2, solve_for_func=False)[0]
    assert checkodesol(eq9, f(x), sol9, order=2, solve_for_func=False)[0]
    assert checkodesol(eq10, f(x), sol10, order=1, solve_for_func=False)[0]
    assert checkodesol(eq11, f(x), sol11, order=1, solve_for_func=False)[0]

def test_1st_linear():
    # Type: first order linear form f'(x)+p(x)f(x)=q(x)
    eq = Eq(f(x).diff(x) + x*f(x), x**2)
    sol = Eq(f(x),exp(-x**2/2)*(sqrt(2)*sqrt(pi)*I*erf(I*x/sqrt(2))/2 \
    + x*exp(x**2/2) + C1))
    assert dsolve(eq, f(x), hint='1st_linear') == sol
    assert checkodesol(eq, f(x), sol, order=1, solve_for_func=False)[0]


def test_Bernoulli():
    # Type: Bernoulli, f'(x) + p(x)*f(x) == q(x)*f(x)**n
    eq = Eq(x*f(x).diff(x) + f(x) - f(x)**2,0)
    sol = dsolve(eq,f(x), hint='Bernoulli')
    assert sol == Eq(f(x),1/(x*(C1 + 1/x)))
    assert checkodesol(eq, f(x), sol, order=1, solve_for_func=False)[0]

def test_1st_exact1():
    # Type: Exact differential equation, p(x,f) + q(x,f)*f' == 0,
    # where dp/df == dq/dx
    eq1 = sin(x)*cos(f(x)) + cos(x)*sin(f(x))*f(x).diff(x)
    eq2 = (2*x*f(x) + 1)/f(x) + (f(x) - x)/f(x)**2*f(x).diff(x)
    eq3 = 2*x + f(x)*cos(x) + (2*f(x) + sin(x) - sin(f(x)))*f(x).diff(x)
    eq4 = cos(f(x)) - (x*sin(f(x)) - f(x)**2)*f(x).diff(x)
    eq5 = 2*x*f(x) + (x**2 + f(x)**2)*f(x).diff(x)
    sol1 = Eq(f(x),acos((C1)/cos(x)))
    sol2 = Eq(log(f(x))+x/f(x)+x**2,C1)
    sol3 = Eq(f(x)*sin(x)+cos(f(x))+x**2+f(x)**2,C1)
    sol4 = Eq(x*cos(f(x))+f(x)**3/3,C1)
    sol5 = Eq(x**2*f(x) + f(x)**3/3, C1)
    assert dsolve(eq1,f(x), hint='1st_exact') == sol1
    assert dsolve(eq2,f(x), hint='1st_exact') == sol2
    assert dsolve(eq3,f(x), hint='1st_exact') == sol3
    assert dsolve(eq4, f(x), hint='1st_exact') == sol4
    assert dsolve(eq5, f(x), hint='1st_exact', simplify=False) == sol5
    assert checkodesol(eq1, f(x), sol1, order=1, solve_for_func=False)[0]
    assert checkodesol(eq2, f(x), sol2, order=1, solve_for_func=False)[0]
    assert checkodesol(eq3, f(x), sol3, order=1, solve_for_func=False)[0]
    assert checkodesol(eq4, f(x), sol4, order=1, solve_for_func=False)[0]
    assert checkodesol(eq5, f(x), sol5, order=1, solve_for_func=False)[0]

@XFAIL
def test_1st_exact2():
    """
    This is an exact equation that fails under the exact engine. It is caught
    by first order homogeneous albeit with a much contorted solution.  The
    exact engine fails because of a poorly simplified integral of q(0,y)dy,
    where q is the function multiplying f'.  The solutions should be
    Eq((x**2+f(x)**2)**Rational(3,2)+y**3, C1).  The equation below is
    equivalent, but it is so complex that checkodesol fails, and takes a long time
    to do so.
    """
    skip("takes too much time")
    eq = x*sqrt(x**2 + f(x)**2) - (x**2*f(x)/(f(x) - sqrt(x**2 + f(x)**2)))*f(x).diff(x)
    sol = dsolve(eq, f(x))
    assert sol == Eq(log(x),C1 - 9*sqrt(1 + f(x)**2/x**2)*asinh(f(x)/x)/(-27*f(x)/x + \
    27*sqrt(1 + f(x)**2/x**2)) - 9*sqrt(1 + f(x)**2/x**2)*log(1 - sqrt(1 + f(x)**2/x**2)*\
    f(x)/x + 2*f(x)**2/x**2)/(-27*f(x)/x + 27*sqrt(1 + f(x)**2/x**2)) \
    + 9*asinh(f(x)/x)*f(x)/(x*(-27*f(x)/x + 27*sqrt(1 + f(x)**2/x**2))) \
    + 9*f(x)*log(1 - sqrt(1 + f(x)**2/x**2)*f(x)/x + 2*f(x)**2/x**2)/\
    (x*(-27*f(x)/x + 27*sqrt(1 + f(x)**2/x**2))))
    assert checkodesol(eq, f(x), sol, order=1, solve_for_func=False)[0]

def test_separable1():
    # test_separable1-5 are from Ordinary Differential Equations, Tenenbaum and
    # Pollard, pg. 55
    eq1 = f(x).diff(x) - f(x)
    eq2 = x*f(x).diff(x) - f(x)
    eq3 = f(x).diff(x) + sin(x)
    eq4 = f(x)**2 + 1 - (x**2 + 1)*f(x).diff(x)
    eq5 = f(x).diff(x)/tan(x) - f(x) - 2
    sol1 = Eq(f(x), exp(C1 + x))
    sol2 = Eq(f(x), C1*x)
    sol3 = Eq(f(x), C1 + cos(x))
    sol4 = Eq(atan(f(x)), C1 + atan(x))
    sol5 = Eq(f(x), -2 + C1*sqrt(1 + tan(x)**2))
    assert dsolve(eq1, f(x), hint='separable') == sol1
    assert dsolve(eq2, f(x), hint='separable') == sol2
    assert dsolve(eq3, f(x), hint='separable') == sol3
    assert dsolve(eq4, f(x), hint='separable') == sol4
    assert dsolve(eq5, f(x), hint='separable') == sol5
    assert checkodesol(eq1, f(x), sol1, order=1, solve_for_func=False)[0]
    assert checkodesol(eq2, f(x), sol2, order=1, solve_for_func=False)[0]
    assert checkodesol(eq3, f(x), sol3, order=1, solve_for_func=False)[0]
    assert checkodesol(eq4, f(x), sol4, order=1, solve_for_func=False)[0]
    assert checkodesol(eq5, f(x), sol5, order=1, solve_for_func=False)[0]

def test_separable2():
    a = Symbol('a')
    eq6 = f(x)*x**2*f(x).diff(x) - f(x)**3 - 2*x**2*f(x).diff(x)
    eq7 = f(x)**2 - 1 - (2*f(x) + x*f(x))*f(x).diff(x)
    eq8 = x*log(x)*f(x).diff(x) + sqrt(1 + f(x)**2)
    eq9 = exp(x + 1)*tan(f(x)) + cos(f(x))*f(x).diff(x)
    eq10 = x*cos(f(x)) + x**2*sin(f(x))*f(x).diff(x) - a**2*sin(f(x))*f(x).diff(x)
    # solve() messes this one up a little bit, so lets test _Integral here
    # We have to test strings with _Integral because y is a dummy variable.
    sol6str = "Integral(-(2 - _y)/_y**3, (_y, None, f(x))) == C1 + Integral(x**(-2), x)"
    sol7 = Eq(log(1 - f(x)**2)/2, C1 + log(2 + x))
    sol8 = Eq(asinh(f(x)), C1 - log(log(x)))
    # integrate cannot handle the integral on the lhs (cos/tan)
    sol9str = "Integral(cos(_y)/tan(_y), (_y, None, f(x))) == C1 + Integral(-E*exp(x), x)"
    sol10 = Eq(-log(1 - sin(f(x))**2)/2, C1 - log(x**2 - a**2)/2)
    assert str(dsolve(eq6, f(x), hint='separable_Integral')) == sol6str
    assert dsolve(eq7, f(x), hint='separable') == sol7
    assert dsolve(eq8, f(x), hint='separable') == sol8
    assert str(dsolve(eq9, f(x), hint='separable_Integral')) == sol9str
    assert dsolve(eq10, f(x), hint='separable') == sol10
    assert checkodesol(eq7, f(x), sol7, order=1, solve_for_func=False)[0]
    assert checkodesol(eq8, f(x), sol8, order=1, solve_for_func=False)[0]

def test_separable3():
    eq11 = f(x).diff(x) - f(x)*tan(x)
    eq12 = (x - 1)*cos(f(x))*f(x).diff(x) - 2*x*sin(f(x))
    eq13 = f(x).diff(x) - f(x)*log(f(x))/tan(x)
    sol11 = Eq(f(x), C1*sqrt(1 + tan(x)**2))
    sol12 = Eq(-log(1 - cos(f(x))**2)/2, C1 - 2*x - 2*log(1 - x))
    sol13 = Eq(log(log(f(x))), C1 - log(1 + tan(x)**2)/2 + log(tan(x)))
    assert dsolve(eq11, f(x), hint='separable') == sol11
    assert dsolve(eq12, f(x), hint='separable') == sol12
    assert dsolve(eq13, f(x), hint='separable') == sol13
    assert checkodesol(eq11, f(x), sol11, order=1, solve_for_func=False)[0]
    assert checkodesol(eq13, f(x), sol13, order=1, solve_for_func=False)[0]

def test_separable4():
    # This has a slow integral (1/((1 + y**2)*atan(y))), so we isolate it.
    eq14 = x*f(x).diff(x) + (1 + f(x)**2)*atan(f(x))
    sol14 = Eq(log(atan(f(x))), C1 - log(x))
    assert dsolve(eq14, f(x), hint='separable') == sol14
    assert checkodesol(eq14, f(x), sol14, order=1, solve_for_func=False)[0]

def test_separable5():
    eq15 = f(x).diff(x) + x*(f(x) + 1)
    eq16 = exp(f(x)**2)*(x**2 + 2*x + 1) + (x*f(x) + f(x))*f(x).diff(x)
    eq17 = f(x).diff(x) + f(x)
    eq18 = sin(x)*cos(2*f(x)) + cos(x)*sin(2*f(x))*f(x).diff(x)
    eq19 = (1 - x)*f(x).diff(x) - x*(f(x) + 1)
    eq20 = f(x)*diff(f(x), x) + x - 3*x*f(x)**2
    eq21 = f(x).diff(x) - exp(x + f(x))
    sol15 = Eq(f(x), -1 + exp(C1 - x**2/2))
    sol16 = Eq(-exp(-f(x)**2)/2, C1 - x - x**2/2)
    sol17 = Eq(f(x), exp(C1 - x))
    sol18 = Eq(-log(1 - sin(2*f(x))**2)/4, C1 + log(1 - sin(x)**2)/2)
    sol19 = Eq(f(x), -(1 - x - exp(C1 - x))/(1 - x))
    sol20 = Eq(-log(1 - 3*f(x)**2)/6, C1 - x**2/2)
    sol21 = Eq(-exp(-f(x)), C1 + exp(x))
    assert dsolve(eq15, f(x), hint='separable') == sol15
    assert dsolve(eq16, f(x), hint='separable') == sol16
    assert dsolve(eq17, f(x), hint='separable') == sol17
    assert dsolve(eq18, f(x), hint='separable') == sol18
    assert dsolve(eq19, f(x), hint='separable') == sol19
    assert dsolve(eq20, f(x), hint='separable') == sol20
    assert dsolve(eq21, f(x), hint='separable') == sol21
    assert checkodesol(eq15, f(x), sol15, order=1, solve_for_func=False)[0]
    assert checkodesol(eq16, f(x), sol16, order=1, solve_for_func=False)[0]
    assert checkodesol(eq17, f(x), sol17, order=1, solve_for_func=False)[0]
    assert checkodesol(eq18, f(x), sol18, order=1, solve_for_func=False)[0]
    assert checkodesol(eq19, f(x), sol19, order=1, solve_for_func=False)[0]
    assert checkodesol(eq20, f(x), sol20, order=1, solve_for_func=False)[0]
    assert checkodesol(eq21, f(x), sol21, order=1, solve_for_func=False)[0]

@XFAIL
def test_separable_1_5_checkodesol():
    # These fail because trigsimp() cannot reduce the expression to 0 in checkodesol()
    a = Symbol('a')
    eq10 = x*cos(f(x)) + x**2*sin(f(x))*f(x).diff(x) - a**2*sin(f(x))*f(x).diff(x)
    eq12 = (x - 1)*cos(f(x))*f(x).diff(x) - 2*x*sin(f(x))
    sol10 = Eq(-log(1 - sin(f(x))**2)/2, C1 - log(x**2 - a**2)/2)
    sol12 = Eq(-log(1 - cos(f(x))**2)/2, C1 - 2*x - 2*log(1 - x))
    assert checkodesol(eq10, f(x), sol10, order=1, solve_for_func=False)[0]
    assert checkodesol(eq12, f(x), sol12, order=1, solve_for_func=False)[0]

def test_homogeneous_order():
    assert homogeneous_order(exp(y/x) + tan(y/x), x, y) == 0
    assert homogeneous_order(x**2 + sin(x)*cos(y), x, y) == None
    assert homogeneous_order(x - y - x*sin(y/x), x, y) == 1
    assert homogeneous_order((x*y + sqrt(x**4+y**4) + x**2*(log(x) - log(y)))/\
        (pi*x**Rational(2,3)*y**Rational(3,2)), x, y) == Rational(-1,6)
    assert homogeneous_order(y/x*cos(y/x) - x/y*sin(y/x) + cos(y/x), x, y) == 0
    assert homogeneous_order(f(x), x, f(x)) == 1
    assert homogeneous_order(f(x)**2, x, f(x)) == 2
    assert homogeneous_order(x*y*z, x, y) == 2
    assert homogeneous_order(x*y*z, x, y, z) == 3
    assert homogeneous_order(x**2*f(x)/sqrt(x**2 + f(x)**2), f(x)) == None
    assert homogeneous_order(f(x,y)**2, x, f(x,y), y) == 2
    assert homogeneous_order(f(x,y)**2, x, f(x), y) == None
    assert homogeneous_order(f(x,y)**2, x, f(x,y)) == None
    assert homogeneous_order(f(y,x)**2, x, y, f(x, y)) == None
    assert homogeneous_order(f(y), f(x), x) == None
    assert homogeneous_order(-f(x)/x + 1/sin(f(x)/ x), f(x), x) == 0
    assert homogeneous_order(log(1/y) + log(x**2), x, y) == None
    assert homogeneous_order(log(1/y) + log(x), x, y) == 0
    assert homogeneous_order(log(x/y), x, y) == 0
    assert homogeneous_order(2*log(1/y) + 2*log(x), x, y) == 0
    a = Symbol('a')
    assert homogeneous_order(a*log(1/y) + a*log(x), x, y) == 0
    assert homogeneous_order(f(x).diff(x), x, y) == None
    assert homogeneous_order(-f(x).diff(x) + x, x, y) == None
    assert homogeneous_order(O(x), x, y) == None
    assert homogeneous_order(x + O(x**2), x, y) == None


def test_1st_homogeneous_coeff_ode1():
    # Type: First order homogeneous, y'=f(y/x)
    eq1 = f(x)/x*cos(f(x)/x) - (x/f(x)*sin(f(x)/x) + cos(f(x)/x))*f(x).diff(x)
    eq2 = x*f(x).diff(x) - f(x) - x*sin(f(x)/x)
    eq3 = f(x) + (x*log(f(x)/x) - 2*x)*diff(f(x),x)
    eq4 = 2*f(x)*exp(x/f(x)) + f(x)*f(x).diff(x) - 2*x*exp(x/f(x))*f(x).diff(x)
    eq5 = 2*x**2*f(x) + f(x)**3 + (x*f(x)**2 - 2*x**3)*f(x).diff(x)
    eq6 = x*exp(f(x)/x) - f(x)*sin(f(x)/x) + x*sin(f(x)/x)*f(x).diff(x)
    eq7 = (x + sqrt(f(x)**2 - x*f(x)))*f(x).diff(x) - f(x)
    eq8 = x+f(x)-(x-f(x))*f(x).diff(x)
    sol1 = Eq(f(x)*sin(f(x)/x), C1)
    sol2 = Eq(x*sqrt(1 + cos(f(x)/x))/sqrt(-1 + cos(f(x)/x)), C1)
    sol3 = Eq(f(x), x*exp(1 - LambertW(C1*x)))
    sol4 = Eq(log(C1*f(x)) + 2*exp(x/f(x)), 0)
    #sol5 = Eq(log(C1*x*sqrt(1/x)*sqrt(f(x))) + x**2/(2*f(x)**2), 0)
    sol5 = Eq(log(C1*x*sqrt(f(x)/x)) + x**2/(2*f(x)**2), 0)
    sol6 = Eq(-exp(-f(x)/x)*sin(f(x)/x)/2 + log(C1*x) - cos(f(x)/x)*exp(-f(x)/x)/2, 0)
    sol7 = Eq(log(C1*f(x)) + 2*sqrt(1 - x/f(x)), 0)
    sol8 = Eq(-atan(f(x)/x) + log(C1*x*sqrt(1 + f(x)**2/x**2)), 0)
    assert dsolve(eq1, f(x), hint='1st_homogeneous_coeff_subs_dep_div_indep') == sol1
    # indep_div_dep actually has a simpler solution for eq2, but it runs too slow
    assert dsolve(eq2, f(x), hint='1st_homogeneous_coeff_subs_dep_div_indep') == sol2
    assert dsolve(eq3, f(x), hint='1st_homogeneous_coeff_best') == sol3
    assert dsolve(eq4, f(x), hint='1st_homogeneous_coeff_best') == sol4
    assert dsolve(eq5, f(x), hint='1st_homogeneous_coeff_best') == sol5
    assert dsolve(eq6, f(x), hint='1st_homogeneous_coeff_subs_dep_div_indep') == sol6
    assert dsolve(eq7, f(x), hint='1st_homogeneous_coeff_best') == sol7
    assert dsolve(eq8, f(x), hint='1st_homogeneous_coeff_best') == sol8

def test_1st_homogeneous_coeff_ode1_sol():
    skip("This test passes, but it takes too long")
    # These are the checkodesols from test_homogeneous_coeff_ode1.
    eq1 = f(x)/x*cos(f(x)/x) - (x/f(x)*sin(f(x)/x) + cos(f(x)/x))*f(x).diff(x)
    eq3 = f(x) + (x*log(f(x)/x) - 2*x)*diff(f(x),x)
    eq4 = 2*f(x)*exp(x/f(x)) + f(x)*f(x).diff(x) - 2*x*exp(x/f(x))*f(x).diff(x)
    eq5 = 2*x**2*f(x) + f(x)**3 + (x*f(x)**2 - 2*x**3)*f(x).diff(x)
    eq6 = x*exp(f(x)/x) - f(x)*sin(f(x)/x) + x*sin(f(x)/x)*f(x).diff(x)
    eq8 = x+f(x)-(x-f(x))*f(x).diff(x)
    sol1 = Eq(f(x)*sin(f(x)/x), C1)
    sol3 = Eq(-f(x)/(1+log(x/f(x))),C1)
    sol4 = Eq(log(C1*f(x)) + 2*exp(x/f(x)), 0)
    sol5 = Eq(log(C1*x*sqrt(1/x)*sqrt(f(x))) + x**2/(2*f(x)**2), 0)
    sol6 = Eq(-exp(-f(x)/x)*sin(f(x)/x)/2 + log(C1*x) - cos(f(x)/x)*exp(-f(x)/x)/2, 0)
    sol8 = Eq(-atan(f(x)/x) + log(C1*x*sqrt(1 + f(x)**2/x**2)), 0)
    assert checkodesol(eq1, f(x), sol1, order=1, solve_for_func=False)[0]
    assert checkodesol(eq3, f(x), sol3, order=1, solve_for_func=False)[0]
    assert checkodesol(eq4, f(x), sol4, order=1, solve_for_func=False)[0]
    assert checkodesol(eq5, f(x), sol5, order=1, solve_for_func=False)[0]
    assert checkodesol(eq6, f(x), sol6, order=1, solve_for_func=False)[0]
    assert checkodesol(eq8, f(x), sol8, order=1, solve_for_func=False)[0]

@XFAIL
def test_1st_homogeneous_coeff_ode1_sol_fail():
    #skip("Takes too long.")
    _u2 = Symbol('u2', dummy=True)
    __a = Symbol('a', dummy=True)
    eq2 = x*f(x).diff(x) - f(x) - x*sin(f(x)/x)
    eq7 = (x + sqrt(f(x)**2 - x*f(x)))*f(x).diff(x) - f(x)
    # test_1st_homogeneous_coeff_ode3
    eq9 = f(x)**2 + (x*sqrt(f(x)**2 - x**2) - x*f(x))*f(x).diff(x)
    sol2 = Eq(x/tan(f(x)/(2*x)), C1)
    sol7 = Eq(log(C1*f(x)) + 2*sqrt(1 - x/f(x)), 0)
    sol9 = Eq(-Integral(-1/(-(1 - (1 - _u2**2)**(1/2))*_u2 + _u2), (_u2, __a, \
    x/f(x))) + log(C1*f(x)), 0)
    assert checkodesol(eq2, f(x), sol2, order=1, solve_for_func=False)[0]
    assert checkodesol(eq7, f(x), sol7, order=1, solve_for_func=False)[0]
    assert checkodesol(eq9, f(x), sol9, order=1, solve_for_func=False)[0]


def test_1st_homogeneous_coeff_ode2():
    eq1 = f(x).diff(x) - f(x)/x+1/sin(f(x)/x)
    eq2 = x**2 + f(x)**2 - 2*x*f(x)*f(x).diff(x)
    eq3 = x*exp(f(x)/x) + f(x) - x*f(x).diff(x)
    sol1 = Eq(f(x), x*acos(log(C1*x)))
    sol2 = [Eq(f(x), sqrt(C1*x + x**2)), Eq(f(x), -sqrt(C1*x + x**2))]
    sol3 = Eq(f(x), log(log(C1/x)**(-x)))
    # specific hints are applied for speed reasons
    assert dsolve(eq1, f(x), hint='1st_homogeneous_coeff_subs_dep_div_indep') == sol1
    assert dsolve(eq2, f(x), hint='1st_homogeneous_coeff_best') == sol2
    assert dsolve(eq3, f(x), hint='1st_homogeneous_coeff_subs_dep_div_indep') == sol3
    assert checkodesol(eq1, f(x), sol1, order=1, solve_for_func=False)[0]
    assert checkodesol(eq2, f(x), sol2[0], order=1, solve_for_func=False)[0]
    assert checkodesol(eq2, f(x), sol2[1], order=1, solve_for_func=False)[0]

@XFAIL
def test_1st_homogeneous_coeff_ode2_eq3sol():
    # simplify() will need to get way better before it can do this one
    eq3 = x*exp(f(x)/x) + f(x) - x*f(x).diff(x)
    sol3 = Eq(f(x), log(log(C1/x)**(-x)))
    assert checkodesol(eq3, f(x), sol3, order=1, solve_for_func=False)[0]

def test_1st_homogeneous_coeff_ode3():
    # This can be solved explicitly, but the the integration engine cannot handle
    # it (see issue 1452).  The explicit solution is included in an XFAIL test
    # below. checkodesol fails for this equation, so its test is in
    # test_homogeneous_order_ode1_sol above. It has to compare string
    # expressions because u2 is a dummy variable.
    eq = f(x)**2+(x*sqrt(f(x)**2-x**2)-x*f(x))*f(x).diff(x)
    solstr = "f(x) == C1*exp(Integral(-1/((1 - _u2**2)**(1/2)*_u2), (_u2, None, x/f(x))))"
    assert str(dsolve(eq, f(x), hint='1st_homogeneous_coeff_subs_indep_div_dep')) == solstr

@XFAIL
def test_1st_homogeneous_coeff_ode4_explicit():
    eq = f(x)**2+(x*sqrt(f(x)**2-x**2)-x*f(x))*f(x).diff(x)
    sol = Eq(f(x)**2+C1*x, f(x)*sqrt(f(x)**2-x**2))
    # If this XPASSES, it means the integral engine has improved!  Please
    # uncomment the real tests below.
    assert not dsolve(eq, f(x)).has(Integral)
#    assert dsolve(eq, f(x)) == sol
#    assert checkodesol(eq, f(x), order=1, solve_for_func=False)[0]

def test_1st_homogeneous_coeff_corner_case():
    eq1 = f(x).diff(x) - f(x)/x
    eq2 = x*f(x).diff(x) - f(x)
    assert "1st_homogeneous_coeff_subs_dep_div_indep" not in classify_ode(eq1, f(x))
    assert "1st_homogeneous_coeff_subs_indep_div_dep" not in classify_ode(eq1, f(x))
    assert "1st_homogeneous_coeff_subs_dep_div_indep" not in classify_ode(eq2, f(x))
    assert "1st_homogeneous_coeff_subs_indep_div_dep" not in classify_ode(eq2, f(x))

def test_nth_linear_constant_coeff_homogeneous():
    # From Exercise 20, in Ordinary Differential Equations, Tenenbaum and Pollard
    # pg. 220
    a = Symbol('a', positive=True)
    k = Symbol('k', real=True)
    eq1 = f(x).diff(x, 2) + 2*f(x).diff(x)
    eq2 = f(x).diff(x, 2) - 3*f(x).diff(x) + 2*f(x)
    eq3 = f(x).diff(x, 2) - f(x)
    eq4 = f(x).diff(x, 3) + f(x).diff(x, 2) - 6*f(x).diff(x)
    eq5 = 6*f(x).diff(x, 2) - 11*f(x).diff(x) + 4*f(x)
    eq6 = Eq(f(x).diff(x, 2) + 2*f(x).diff(x) - f(x), 0)
    eq7 = diff(f(x), x, 3) + diff(f(x), x, 2) - 10*diff(f(x), x) - 6*f(x)
    eq8 = f(x).diff(x, 4) - f(x).diff(x, 3) - 4*f(x).diff(x, 2) + 4*f(x).diff(x)
    eq9 = f(x).diff(x, 4) + 4*f(x).diff(x, 3) + f(x).diff(x, 2) - \
        4*f(x).diff(x) - 2*f(x)
    eq10 = f(x).diff(x, 4) - a**2*f(x)
    eq11 = f(x).diff(x, 2) - 2*k*f(x).diff(x) - 2*f(x)
    eq12 = f(x).diff(x, 2) + 4*k*f(x).diff(x) - 12*k**2*f(x)
    eq13 = f(x).diff(x, 4)
    eq14 = f(x).diff(x, 2) + 4*f(x).diff(x) + 4*f(x)
    eq15 = 3*f(x).diff(x, 3) + 5*f(x).diff(x, 2) + f(x).diff(x) - f(x)
    eq16 = f(x).diff(x, 3) - 6*f(x).diff(x, 2) +12*f(x).diff(x) - 8*f(x)
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
    sol2 = Eq(f(x), C1*exp(x) + C2*exp(2*x))
    sol3 = Eq(f(x), C1*exp(x) + C2*exp(-x))
    sol4 = Eq(f(x), C1 + C2*exp(-3*x) + C3*exp(2*x))
    sol5 = Eq(f(x), C1*exp(x/2) + C2*exp(4*x/3))
    sol6 = Eq(f(x), C1*exp(-x + x*sqrt(2)) + C2*exp(-x - x*sqrt(2)))
    sol7 = Eq(f(x), C1*exp(3*x) + C2*exp(-2*x + x*sqrt(2)) + C3*exp(-2*x - x*sqrt(2)))
    sol8 = Eq(f(x), C1 + C2*exp(x) + C3*exp(-2*x) + C4*exp(2*x))
    sol9 = Eq(f(x), C1*exp(x) + C2*exp(-x) + C3*exp(-2*x + x*sqrt(2)) + \
        C4*exp(-2*x - x*sqrt(2)))
    sol10 = Eq(f(x), C1*sin(x*sqrt(a)) + C2*cos(x*sqrt(a)) + C3*exp(x*sqrt(a)) + \
        C4*exp(-x*sqrt(a)))
    sol11 = Eq(f(x), C1*exp(k*x + x*sqrt(8 + 4*k**2)/2) + \
        C2*exp(k*x - x*sqrt(8 + 4*k**2)/2))
    sol12 = Eq(f(x), C1*exp(-4*x*abs(k) - 2*k*x) + C2*exp(-2*k*x + 4*x*abs(k)))
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
    sol25 = Eq(f(x), C1*cos(x*sqrt(3)) + C2*sin(x*sqrt(3)) + C3*sin(x*sqrt(2)) + \
    C4*cos(x*sqrt(2)))
    sol26 = Eq(f(x), (C1*sin(4*x) + C2*cos(4*x))*exp(2*x))
    sol27 = Eq(f(x), (C1 + C2*x)*sin(x*sqrt(2)) + (C3 + C4*x)*cos(x*sqrt(2)))
    sol28 = Eq(f(x), (C1*sin(x*sqrt(3)) + C2*cos(x*sqrt(3)))*exp(x) + C3*exp(-2*x))
    sol29 = Eq(f(x), C1 + C2*sin(2*x) + C3*cos(2*x) + C4*x)
    sol30 = Eq(f(x), C1 + (C2 + C3*x)*sin(x) + (C4 + C5*x)*cos(x))
    sol1s = ode_renumber(sol1, 'C', 1, 2)
    sol2s = ode_renumber(sol2, 'C', 1, 2)
    sol3s = ode_renumber(sol3, 'C', 1, 2)
    sol4s = ode_renumber(sol4, 'C', 1, 3)
    sol5s = ode_renumber(sol5, 'C', 1, 2)
    sol6s = ode_renumber(sol6, 'C', 1, 2)
    sol7s = ode_renumber(sol7, 'C', 1, 3)
    sol8s = ode_renumber(sol8, 'C', 1, 4)
    sol9s = ode_renumber(sol9, 'C', 1, 4)
    sol10s = ode_renumber(sol10, 'C', 1, 4)
    sol11s = ode_renumber(sol11, 'C', 1, 2)
    sol12s = ode_renumber(sol12, 'C', 1, 2)
    sol13s = ode_renumber(sol13, 'C', 1, 4)
    sol14s = ode_renumber(sol14, 'C', 1, 2)
    sol15s = ode_renumber(sol15, 'C', 1, 3)
    sol16s = ode_renumber(sol16, 'C', 1, 3)
    sol17s = ode_renumber(sol17, 'C', 1, 2)
    sol18s = ode_renumber(sol18, 'C', 1, 4)
    sol19s = ode_renumber(sol19, 'C', 1, 4)
    sol20s = ode_renumber(sol20, 'C', 1, 4)
    sol21s = ode_renumber(sol21, 'C', 1, 4)
    sol22s = ode_renumber(sol22, 'C', 1, 4)
    sol23s = ode_renumber(sol23, 'C', 1, 2)
    sol24s = ode_renumber(sol24, 'C', 1, 2)
    sol25s = ode_renumber(sol25, 'C', 1, 4)
    sol26s = ode_renumber(sol26, 'C', 1, 2)
    sol27s = ode_renumber(sol27, 'C', 1, 4)
    sol28s = ode_renumber(sol28, 'C', 1, 3)
    sol29s = ode_renumber(sol29, 'C', 1, 4)
    sol30s = ode_renumber(sol30, 'C', 1, 5)
    assert dsolve(eq1, f(x)) in (sol1, sol1s)
    assert dsolve(eq2, f(x)) in (sol2, sol2s)
    assert dsolve(eq3, f(x)) in (sol3, sol3s)
    assert dsolve(eq4, f(x)) in (sol4, sol4s)
    assert dsolve(eq5, f(x)) in (sol5, sol5s)
    assert dsolve(eq6, f(x)) in (sol6, sol6s)
    assert dsolve(eq7, f(x)) in (sol7, sol7s)
    assert dsolve(eq8, f(x)) in (sol8, sol8s)
    assert dsolve(eq9, f(x)) in (sol9, sol9s)
    assert dsolve(eq10, f(x)) in (sol10, sol10s)
    assert dsolve(eq11, f(x)) in (sol11, sol11s)
    assert dsolve(eq12, f(x)) in (sol12, sol12s)
    assert dsolve(eq13, f(x)) in (sol13, sol13s)
    assert dsolve(eq14, f(x)) in (sol14, sol14s)
    assert dsolve(eq15, f(x)) in (sol15, sol15s)
    assert dsolve(eq16, f(x)) in (sol16, sol16s)
    assert dsolve(eq17, f(x)) in (sol17, sol17s)
    assert dsolve(eq18, f(x)) in (sol18, sol18s)
    assert dsolve(eq19, f(x)) in (sol19, sol19s)
    assert dsolve(eq20, f(x)) in (sol20, sol20s)
    assert dsolve(eq21, f(x)) in (sol21, sol21s)
    assert dsolve(eq22, f(x)) in (sol22, sol22s)
    assert dsolve(eq23, f(x)) in (sol23, sol23s)
    assert dsolve(eq24, f(x)) in (sol24, sol24s)
    assert dsolve(eq25, f(x)) in (sol25, sol25s)
    assert dsolve(eq26, f(x)) in (sol26, sol26s)
    assert dsolve(eq27, f(x)) in (sol27, sol27s)
    assert dsolve(eq28, f(x)) in (sol28, sol28s)
    assert dsolve(eq29, f(x)) in (sol29, sol29s)
    assert dsolve(eq30, f(x)) in (sol30, sol30s)
    assert checkodesol(eq1, f(x), sol1, order=2, solve_for_func=False)[0]
    assert checkodesol(eq2, f(x), sol2, order=2, solve_for_func=False)[0]
    assert checkodesol(eq3, f(x), sol3, order=2, solve_for_func=False)[0]
    assert checkodesol(eq4, f(x), sol4, order=3, solve_for_func=False)[0]
    assert checkodesol(eq5, f(x), sol5, order=2, solve_for_func=False)[0]
    assert checkodesol(eq6, f(x), sol6, order=2, solve_for_func=False)[0]
    assert checkodesol(eq7, f(x), sol7, order=3, solve_for_func=False)[0]
    assert checkodesol(eq8, f(x), sol8, order=4, solve_for_func=False)[0]
    assert checkodesol(eq9, f(x), sol9, order=4, solve_for_func=False)[0]
    assert checkodesol(eq10, f(x), sol10, order=4, solve_for_func=False)[0]
    assert checkodesol(eq11, f(x), sol11, order=2, solve_for_func=False)[0]
    assert checkodesol(eq12, f(x), sol12, order=2, solve_for_func=False)[0]
    assert checkodesol(eq13, f(x), sol13, order=4, solve_for_func=False)[0]
    assert checkodesol(eq14, f(x), sol14, order=2, solve_for_func=False)[0]
    assert checkodesol(eq15, f(x), sol15, order=3, solve_for_func=False)[0]
    assert checkodesol(eq16, f(x), sol16, order=3, solve_for_func=False)[0]
    assert checkodesol(eq17, f(x), sol17, order=2, solve_for_func=False)[0]
    assert checkodesol(eq18, f(x), sol18, order=4, solve_for_func=False)[0]
    assert checkodesol(eq19, f(x), sol19, order=4, solve_for_func=False)[0]
    assert checkodesol(eq20, f(x), sol20, order=4, solve_for_func=False)[0]
    assert checkodesol(eq21, f(x), sol21, order=4, solve_for_func=False)[0]
    assert checkodesol(eq22, f(x), sol22, order=4, solve_for_func=False)[0]
    assert checkodesol(eq23, f(x), sol23, order=2, solve_for_func=False)[0]
    assert checkodesol(eq24, f(x), sol24, order=2, solve_for_func=False)[0]
    assert checkodesol(eq25, f(x), sol25, order=4, solve_for_func=False)[0]
    assert checkodesol(eq26, f(x), sol26, order=2, solve_for_func=False)[0]
    assert checkodesol(eq27, f(x), sol27, order=4, solve_for_func=False)[0]
    assert checkodesol(eq28, f(x), sol28, order=3, solve_for_func=False)[0]
    assert checkodesol(eq29, f(x), sol29, order=4, solve_for_func=False)[0]
    assert checkodesol(eq30, f(x), sol30, order=5, solve_for_func=False)[0]

def test_nth_linear_constant_coeff_homogeneous_RootOf():
    # We have to test strings because _m is a dummy variable
    _m = Symbol('_m')
    eq = f(x).diff(x, 5) + 11*f(x).diff(x) - 2*f(x)
    solstr = "f(x) == C1*exp(x*RootOf(_m**5 + 11*_m - 2, _m, domain='ZZ', index=0)) + C2" + \
        "*exp(x*RootOf(_m**5 + 11*_m - 2, _m, domain='ZZ', index=1)) + C3*exp(x*RootOf(_" + \
        "m**5 + 11*_m - 2, _m, domain='ZZ', index=2)) + C4*exp(x*RootOf(_m**5 + 11*_m - " + \
        "2, _m, domain='ZZ', index=3)) + C5*exp(x*RootOf(_m**5 + 11*_m - 2, _m, domain='ZZ', index=4))"
    assert str(dsolve(eq, f(x))) == solstr

@XFAIL
def test_nth_linear_constant_coeff_homogeneous_RootOf_sol():
    # We have to test strings because _m is a dummy variable
    _m = Symbol('_m')
    eq = f(x).diff(x, 5) + 11*f(x).diff(x) - 2*f(x)
    sol = Eq(f(x), C1*exp(x*RootOf(Poly(_m**5 + 11*_m - 2, _m), index=0)) + \
        C2*exp(x*RootOf(Poly(_m**5 + 11*_m - 2, _m), index=1)) + \
        C3*exp(x*RootOf(Poly(_m**5 + 11*_m - 2, _m), index=2)) + \
        C4*exp(x*RootOf(Poly(_m**5 + 11*_m - 2, _m), index=3)) + \
        C5*exp(x*RootOf(Poly(_m**5 + 11*_m - 2, _m), index=4)))
    solstr = "f(x) == C1*exp(x*RootOf(_m**5 + 11*_m - 2, _m, index=0)) + C2" + \
        "*exp(x*RootOf(_m**5 + 11*_m - 2, _m, index=1)) + C3*exp(x*RootOf(_" + \
        "m**5 + 11*_m - 2, _m, index=2)) + C4*exp(x*RootOf(_m**5 + 11*_m - " + \
        "2, _m, index=3)) + C5*exp(x*RootOf(_m**5 + 11*_m - 2, _m, index=4))"
    assert str(sol) == solstr # str(sol) fails
    assert checkodesol(eq, f(x), sol, order=5, solve_for_func=False)[0]

def test_undetermined_coefficients_match():
    assert _undetermined_coefficients_match(g(x), x) == {'test': False}
    assert _undetermined_coefficients_match(sin(2*x + sqrt(5)), x) == \
        {'test': True, 'trialset': set([cos(2*x + sqrt(5)), sin(2*x + sqrt(5))])}
    assert _undetermined_coefficients_match(sin(x)*cos(x), x) == {'test': False}
    assert _undetermined_coefficients_match(sin(x)*(x**2 + x + 1), x) == \
        {'test': True, 'trialset': set([cos(x), x*cos(x), x**2*sin(x), x**2*cos(x),
        x*sin(x), sin(x)])}
    assert _undetermined_coefficients_match(sin(x)*x**2 + sin(x)*x + sin(x), x) == \
        {'test': True, 'trialset': set([x*cos(x), x*sin(x), x**2*cos(x), x**2*sin(x),
        cos(x), sin(x)])}
    assert _undetermined_coefficients_match(exp(2*x)*sin(x)*(x**2 + x + 1), x) == \
        {'test': True, 'trialset': set([exp(2*x)*sin(x), x**2*exp(2*x)*sin(x),
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
    assert _undetermined_coefficients_match(sin(x)*(x + sin(x)), x) == {'test': False}
    assert _undetermined_coefficients_match(sin(x)*(x + sin(2*x)), x) == {'test': False}
    assert _undetermined_coefficients_match(sin(x)*tan(x), x) == {'test': False}
    assert _undetermined_coefficients_match(x**2*sin(x)*exp(x) + x*sin(x) + x, x) == \
        {'test': True, 'trialset': set([x**2*cos(x)*exp(x), x, cos(x), S(1),
        exp(x)*sin(x), sin(x), x*exp(x)*sin(x), x*cos(x), x*cos(x)*exp(x),
        x*sin(x), cos(x)*exp(x), x**2*exp(x)*sin(x)])}
    assert _undetermined_coefficients_match(4*x*sin(x - 2), x) == \
        {'test': True, 'trialset': set([x*cos(x - 2), x*sin(x - 2), cos(x - 2),
        sin(x - 2)])}
    assert _undetermined_coefficients_match(2**x*x, x) == \
        {'test': True, 'trialset': set([2**x, x*2**x])}
    assert _undetermined_coefficients_match(2**x*exp(2*x), x) == \
        {'test': True, 'trialset': set([2**x*exp(2*x)])}
    assert _undetermined_coefficients_match(exp(-x)/x, x) == \
        {'test': False}
    # Below are from Ordinary Differential Equations, Tenenbaum and Pollard, pg. 231
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
        {'test': True, 'trialset': set([x*cos(2*x), x*sin(2*x), cos(2*x), sin(2*x)])}
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
    assert _undetermined_coefficients_match(exp(2*x)*(S(1)/2 + cos(2*x)/2), x) == \
        {'test': True, 'trialset': set([exp(2*x)*sin(2*x), cos(2*x)*exp(2*x),
        exp(2*x)])}
    assert _undetermined_coefficients_match(2*x + sin(x) + cos(x), x) == \
        {'test': True, 'trialset': set([S(1), x, cos(x), sin(x)])}
    # converted from sin(2*x)*sin(x)
    assert _undetermined_coefficients_match(cos(x)/2 - cos(3*x)/2, x) == \
        {'test': True, 'trialset': set([cos(x), cos(3*x), sin(x), sin(3*x)])}

def test_nth_linear_constant_coeff_undetermined_coefficients():
    hint = 'nth_linear_constant_coeff_undetermined_coefficients'
    eq1 = 3*f(x).diff(x, 3) + 5*f(x).diff(x, 2) + f(x).diff(x) - f(x) - x*exp(-x) - x
    eq2 = 3*f(x).diff(x, 3) + 5*f(x).diff(x, 2) + f(x).diff(x) - f(x) - exp(-x) - x
    # 3-27 below are from Ordinary Differential Equations, Tenenbaum and Pollard, pg. 231
    eq3 = f(x).diff(x, 2) + 3*f(x).diff(x) + 2*f(x) - 4
    eq4 = f(x).diff(x, 2) + 3*f(x).diff(x) + 2*f(x) - 12*exp(x)
    eq5 = f(x).diff(x, 2) + 3*f(x).diff(x) + 2*f(x) - exp(I*x)
    eq6 = f(x).diff(x, 2) + 3*f(x).diff(x) + 2*f(x) - sin(x)
    eq7 = f(x).diff(x, 2) + 3*f(x).diff(x) + 2*f(x) - cos(x)
    eq8 = f(x).diff(x, 2) + 3*f(x).diff(x) + 2*f(x) - (8 + 6*exp(x) + 2*sin(x))
    eq9 = f(x).diff(x, 2) + f(x).diff(x) + f(x) - x**2
    eq10 = f(x).diff(x, 2) - 2*f(x).diff(x) - 8*f(x) - 9*x*exp(x) - 10*exp(-x)
    eq11 = f(x).diff(x, 2) - 3*f(x).diff(x) - 2*exp(2*x)*sin(x)
    eq12 = f(x).diff(x, 4) - 2*f(x).diff(x, 2) + f(x) - x + sin(x)
    eq13 = f(x).diff(x, 2) + f(x).diff(x) - x**2 - 2*x
    eq14 = f(x).diff(x, 2) + f(x).diff(x) - x - sin(2*x)
    eq15 = f(x).diff(x, 2) + f(x) - 4*x*sin(x)
    eq16 = f(x).diff(x, 2) + 4*f(x) - x*sin(2*x)
    eq17 = f(x).diff(x, 2) + 2*f(x).diff(x) + f(x) - x**2*exp(-x)
    eq18 = f(x).diff(x, 3) + 3*f(x).diff(x, 2) + 3*f(x).diff(x) + f(x) - 2*exp(-x) + \
        x**2*exp(-x)
    eq19 = f(x).diff(x, 2) + 3*f(x).diff(x) + 2*f(x) - exp(-2*x) - x**2
    eq20 = f(x).diff(x, 2) - 3*f(x).diff(x) + 2*f(x) - x*exp(-x)
    eq21 = f(x).diff(x, 2) + f(x).diff(x) - 6*f(x) - x - exp(2*x)
    eq22 = f(x).diff(x, 2) + f(x) - sin(x) - exp(-x)
    eq23 = f(x).diff(x, 3) - 3*f(x).diff(x, 2) + 3*f(x).diff(x) - f(x) - exp(x)
    # sin(x)**2
    eq24 = f(x).diff(x, 2) + f(x) - S(1)/2 - cos(2*x)/2
    # exp(2*x)*sin(x)**2
    eq25 = f(x).diff(x, 3) - f(x).diff(x) - exp(2*x)*(S(1)/2 - cos(2*x)/2)
    eq26 = f(x).diff(x, 5) + 2*f(x).diff(x, 3) + f(x).diff(x) - 2*x - sin(x) - cos(x)
    # sin(2*x)*sin(x), skip 3127 for now, match bug
    eq27 = f(x).diff(x, 2) + f(x) - cos(x)/2 + cos(3*x)/2
    eq28 = f(x).diff(x) - 1
    sol1 = Eq(f(x), -1 - x + (C1 + C2*x - 3*x**2/32 - x**3/24)*exp(-x) + C3*exp(x/3))
    sol2 = Eq(f(x), -1 - x + (C1 + C2*x - x**2/8)*exp(-x) + C3*exp(x/3))
    sol3 = Eq(f(x), 2 + C1*exp(-x) + C2*exp(-2*x))
    sol4 = Eq(f(x), 2*exp(x) + C1*exp(-x) + C2*exp(-2*x))
    sol5 = Eq(f(x), C1*exp(-x) + C2*exp(-2*x) + exp(I*x)/10 - 3*I*exp(I*x)/10)
    sol6 = Eq(f(x), -3*cos(x)/10 + sin(x)/10 + C1*exp(-x) + C2*exp(-2*x))
    sol7 = Eq(f(x), cos(x)/10 + 3*sin(x)/10 + C1*exp(-x) + C2*exp(-2*x))
    sol8 = Eq(f(x), 4 - 3*cos(x)/5 + sin(x)/5 + exp(x) + C1*exp(-x) + C2*exp(-2*x))
    sol9 = Eq(f(x), -2*x + x**2 + (C1*sin(x*sqrt(3)/2) + C2*cos(x*sqrt(3)/2))*exp(-x/2))
    sol10 = Eq(f(x), -x*exp(x) - 2*exp(-x) + C1*exp(-2*x) + C2*exp(4*x))
    sol11 = Eq(f(x), C1 + (-3*sin(x)/5 - cos(x)/5)*exp(2*x) + C2*exp(3*x))
    sol12 = Eq(f(x), x - sin(x)/4 + (C1 + C2*x)*exp(x) + (C3 + C4*x)*exp(-x))
    sol13 = Eq(f(x), C1 + x**3/3 + C2*exp(-x))
    sol14 = Eq(f(x), C1 - x - sin(2*x)/5 - cos(2*x)/10 + x**2/2 + C2*exp(-x))
    sol15 = Eq(f(x), (C1 + x)*sin(x) + (C2 - x**2)*cos(x))
    sol16 = Eq(f(x), (C1 + x/16)*sin(2*x) + (C2 - x**2/8)*cos(2*x))
    sol17 = Eq(f(x), (C1 + C2*x + x**4/12)*exp(-x))
    sol18 = Eq(f(x), (C1 + C2*x + C3*x**2 + x**3/3 - x**5/60)*exp(-x))
    sol19 = Eq(f(x), S(7)/4 - 3*x/2 + x**2/2 + C1*exp(-x) + (C2 - x)*exp(-2*x))
    sol20 = Eq(f(x), C1*exp(x) + (S(5)/36 + x/6)*exp(-x) + C2*exp(2*x))
    sol21 = Eq(f(x), -S(1)/36 - x/6 + C1*exp(-3*x) + (C2 + x/5)*exp(2*x))
    sol22 = Eq(f(x), C1*sin(x) + (C2 - x/2)*cos(x) + exp(-x)/2)
    sol23 = Eq(f(x), (C1 + C2*x + C3*x**2 + x**3/6)*exp(x))
    sol24 = Eq(f(x), S(1)/2 - cos(2*x)/6 + C1*sin(x) + C2*cos(x))
    sol25 = Eq(f(x), C1 + C2*exp(x) + C3*exp(-x) + (S(1)/12 - 7*sin(2*x)/520 + \
        9*cos(2*x)/520)*exp(2*x))
    sol26 = Eq(f(x), C1 + (C2 + C3*x - x**2/8)*sin(x) + (C4 + C5*x + x**2/8)*cos(x) + x**2)
    sol27 = Eq(f(x), cos(3*x)/16 + C1*cos(x) + (C2 + x/4)*sin(x))
    sol28 = Eq(f(x), C1 + x)
    sol1s = ode_renumber(sol1, 'C', 1, 3)
    sol2s = ode_renumber(sol2, 'C', 1, 3)
    sol3s = ode_renumber(sol3, 'C', 1, 2)
    sol4s = ode_renumber(sol4, 'C', 1, 2)
    sol5s = ode_renumber(sol5, 'C', 1, 2)
    sol6s = ode_renumber(sol6, 'C', 1, 2)
    sol7s = ode_renumber(sol7, 'C', 1, 2)
    sol8s = ode_renumber(sol8, 'C', 1, 2)
    sol9s = ode_renumber(sol9, 'C', 1, 2)
    sol10s = ode_renumber(sol10, 'C', 1, 2)
    sol11s = ode_renumber(sol11, 'C', 1, 2)
    sol12s = ode_renumber(sol12, 'C', 1, 2)
    sol13s = ode_renumber(sol13, 'C', 1, 4)
    sol14s = ode_renumber(sol14, 'C', 1, 2)
    sol15s = ode_renumber(sol15, 'C', 1, 2)
    sol16s = ode_renumber(sol16, 'C', 1, 2)
    sol17s = ode_renumber(sol17, 'C', 1, 2)
    sol18s = ode_renumber(sol18, 'C', 1, 3)
    sol19s = ode_renumber(sol19, 'C', 1, 2)
    sol20s = ode_renumber(sol20, 'C', 1, 2)
    sol21s = ode_renumber(sol21, 'C', 1, 2)
    sol22s = ode_renumber(sol22, 'C', 1, 2)
    sol23s = ode_renumber(sol23, 'C', 1, 3)
    sol24s = ode_renumber(sol24, 'C', 1, 2)
    sol25s = ode_renumber(sol25, 'C', 1, 3)
    sol26s = ode_renumber(sol26, 'C', 1, 5)
    sol27s = ode_renumber(sol27, 'C', 1, 2)
    assert dsolve(eq1, f(x), hint=hint) in (sol1, sol1s)
    assert dsolve(eq2, f(x), hint=hint) in (sol2, sol2s)
    assert dsolve(eq3, f(x), hint=hint) in (sol3, sol3s)
    assert dsolve(eq4, f(x), hint=hint) in (sol4, sol4s)
    assert dsolve(eq5, f(x), hint=hint) in (sol5, sol5s)
    assert dsolve(eq6, f(x), hint=hint) in (sol6, sol6s)
    assert dsolve(eq7, f(x), hint=hint) in (sol7, sol7s)
    assert dsolve(eq8, f(x), hint=hint) in (sol8, sol8s)
    assert dsolve(eq9, f(x), hint=hint) in (sol9, sol9s)
    assert dsolve(eq10, f(x), hint=hint) in (sol10, sol10s)
    assert dsolve(eq11, f(x), hint=hint) in (sol11, sol11s)
    assert dsolve(eq12, f(x), hint=hint) in (sol12, sol12s)
    assert dsolve(eq13, f(x), hint=hint) in (sol13, sol13s)
    assert dsolve(eq14, f(x), hint=hint) in (sol14, sol14s)
    assert dsolve(eq15, f(x), hint=hint) in (sol15, sol15s)
    assert dsolve(eq16, f(x), hint=hint) in (sol16, sol16s)
    assert dsolve(eq17, f(x), hint=hint) in (sol17, sol17s)
    assert dsolve(eq18, f(x), hint=hint) in (sol18, sol18s)
    assert dsolve(eq19, f(x), hint=hint) in (sol19, sol19s)
    assert dsolve(eq20, f(x), hint=hint) in (sol20, sol20s)
    assert dsolve(eq21, f(x), hint=hint) in (sol21, sol21s)
    assert dsolve(eq22, f(x), hint=hint) in (sol22, sol22s)
    assert dsolve(eq23, f(x), hint=hint) in (sol23, sol23s)
    assert dsolve(eq24, f(x), hint=hint) in (sol24, sol24s)
    assert dsolve(eq25, f(x), hint=hint) in (sol25, sol25s)
    assert dsolve(eq26, f(x), hint=hint) in (sol26, sol26s)
    assert dsolve(eq27, f(x), hint=hint) in (sol27, sol27s)
    assert dsolve(eq28, f(x), hint=hint) == sol28
    assert checkodesol(eq1, f(x), sol1, order=3, solve_for_func=False)[0]
    assert checkodesol(eq2, f(x), sol2, order=3, solve_for_func=False)[0]
    assert checkodesol(eq3, f(x), sol3, order=2, solve_for_func=False)[0]
    assert checkodesol(eq4, f(x), sol4, order=2, solve_for_func=False)[0]
    assert checkodesol(eq5, f(x), sol5, order=2, solve_for_func=False)[0]
    assert checkodesol(eq6, f(x), sol6, order=2, solve_for_func=False)[0]
    assert checkodesol(eq7, f(x), sol7, order=2, solve_for_func=False)[0]
    assert checkodesol(eq8, f(x), sol8, order=2, solve_for_func=False)[0]
    assert checkodesol(eq9, f(x), sol9, order=2, solve_for_func=False)[0]
    assert checkodesol(eq10, f(x), sol10, order=2, solve_for_func=False)[0]
    assert checkodesol(eq11, f(x), sol11, order=2, solve_for_func=False)[0]
    assert checkodesol(eq12, f(x), sol12, order=4, solve_for_func=False)[0]
    assert checkodesol(eq13, f(x), sol13, order=2, solve_for_func=False)[0]
    assert checkodesol(eq14, f(x), sol14, order=2, solve_for_func=False)[0]
    assert checkodesol(eq15, f(x), sol15, order=2, solve_for_func=False)[0]
    assert checkodesol(eq16, f(x), sol16, order=2, solve_for_func=False)[0]
    assert checkodesol(eq17, f(x), sol17, order=2, solve_for_func=False)[0]
    assert checkodesol(eq18, f(x), sol18, order=3, solve_for_func=False)[0]
    assert checkodesol(eq19, f(x), sol19, order=2, solve_for_func=False)[0]
    assert checkodesol(eq20, f(x), sol20, order=2, solve_for_func=False)[0]
    assert checkodesol(eq21, f(x), sol21, order=2, solve_for_func=False)[0]
    assert checkodesol(eq22, f(x), sol22, order=2, solve_for_func=False)[0]
    assert checkodesol(eq23, f(x), sol23, order=3, solve_for_func=False)[0]
    assert checkodesol(eq24, f(x), sol24, order=2, solve_for_func=False)[0]
    assert checkodesol(eq25, f(x), sol25, order=3, solve_for_func=False)[0]
    assert checkodesol(eq26, f(x), sol26, order=5, solve_for_func=False)[0]
    assert checkodesol(eq27, f(x), sol27, order=2, solve_for_func=False)[0]
    assert checkodesol(eq28, f(x), sol28, order=1, solve_for_func=False)[0]

@XFAIL
def test_nth_linear_constant_coeff_undetermined_coefficients_imaginary_exp():
    # Equivalent to eq26, in test_nth_linear_constant_coeff_undetermined_coefficients
    # above.  This fails because the algorithm for undetermined coefficients
    # doesn't know to multiply exp(I*x) by sufficient x because it is linearly
    # dependent on sin(x) and cos(x).
    eq26a = f(x).diff(x, 5) + 2*f(x).diff(x, 3) + f(x).diff(x) - 2*x - exp(I*x)
    sol26 = Eq(f(x), C1 + (C2 + C3*x - x**2/8)*sin(x) + (C4 + C5*x + x**2/8)*cos(x) + x**2)
    assert dsolve(eq26a, f(x), hint=hint) == sol26
    assert checkodesol(eq26a, f(x), sol26, order=5, solve_for_func=False)[0]

def test_nth_linear_constant_coeff_variation_of_parameters():
    hint = 'nth_linear_constant_coeff_variation_of_parameters'
    eq1 = 3*f(x).diff(x, 3) + 5*f(x).diff(x, 2) + f(x).diff(x) - f(x) - x*exp(-x) - x
    eq2 = 3*f(x).diff(x, 3) + 5*f(x).diff(x, 2) + f(x).diff(x) - f(x) - exp(-x) - x
    eq3 = f(x).diff(x) - 1
    eq4 = f(x).diff(x, 2) + 3*f(x).diff(x) + 2*f(x) - 4
    eq5 = f(x).diff(x, 2) + 3*f(x).diff(x) + 2*f(x) - 12*exp(x)
    eq6 = f(x).diff(x, 2) - 2*f(x).diff(x) - 8*f(x) - 9*x*exp(x) - 10*exp(-x)
    eq7 = f(x).diff(x, 2) + 2*f(x).diff(x) + f(x) - x**2*exp(-x)
    eq8 = f(x).diff(x, 2) - 3*f(x).diff(x) + 2*f(x) - x*exp(-x)
    eq9 = f(x).diff(x, 3) - 3*f(x).diff(x, 2) + 3*f(x).diff(x) - f(x) - exp(x)
    eq10 = f(x).diff(x, 2) + 2*f(x).diff(x) + f(x) - exp(-x)/x
    eq11 = f(x).diff(x, 2) + f(x) - 1/sin(x)*1/cos(x)
    eq12 = f(x).diff(x, 4)  - 1/x
    sol1 = Eq(f(x), -1 - x - (C1 + C2*x + 3*x**2/32 + x**3/24)*exp(-x) + C3*exp(x/3))
    sol2 = Eq(f(x), -1 - x - (C1 + C2*x + x**2/8)*exp(-x) + C3*exp(x/3))
    sol3 = Eq(f(x), C1 + x)
    sol4 = Eq(f(x), 2 + C1*exp(-x) + C2*exp(-2*x))
    sol5 = Eq(f(x), 2*exp(x) + C1*exp(-x) + C2*exp(-2*x))
    sol6 = Eq(f(x), -x*exp(x) - 2*exp(-x) + C1*exp(-2*x) + C2*exp(4*x))
    sol7 = Eq(f(x), (C1 + C2*x + x**4/12)*exp(-x))
    sol8 = Eq(f(x), C1*exp(x) + (S(5)/36 + x/6)*exp(-x) + C2*exp(2*x))
    sol9 = Eq(f(x), (C1 + C2*x + C3*x**2 + x**3/6)*exp(x))
    sol10 = Eq(f(x), (C1 - x*(C2 - log(x)))*exp(-x))
    sol11 = Eq(f(x), cos(x)*(C1 - Integral(1/cos(x), x)) + sin(x)*(C2 + \
        Integral(1/sin(x), x)))
    sol12 = Eq(f(x), C1 + C2*x - x**3*(C3 - log(x)/6) + C4*x**2)
    sol1s = ode_renumber(sol1, 'C', 1, 3)
    sol2s = ode_renumber(sol2, 'C', 1, 3)
    sol3s = ode_renumber(sol3, 'C', 1, 2)
    sol4s = ode_renumber(sol4, 'C', 1, 2)
    sol5s = ode_renumber(sol5, 'C', 1, 2)
    sol6s = ode_renumber(sol6, 'C', 1, 2)
    sol7s = ode_renumber(sol7, 'C', 1, 2)
    sol8s = ode_renumber(sol8, 'C', 1, 2)
    sol9s = ode_renumber(sol9, 'C', 1, 3)
    sol10s = ode_renumber(sol10, 'C', 1, 2)
    sol11s = ode_renumber(sol11, 'C', 1, 2)
    sol12s = ode_renumber(sol12, 'C', 1, 4)
    assert dsolve(eq1, f(x), hint=hint) in (sol1, sol1s)
    assert dsolve(eq2, f(x), hint=hint) in (sol2, sol2s)
    assert dsolve(eq3, f(x), hint=hint) in (sol3, sol3s)
    assert dsolve(eq4, f(x), hint=hint) in (sol4, sol4s)
    assert dsolve(eq5, f(x), hint=hint) in (sol5, sol5s)
    assert dsolve(eq6, f(x), hint=hint) in (sol6, sol6s)
    assert dsolve(eq7, f(x), hint=hint) in (sol7, sol7s)
    assert dsolve(eq8, f(x), hint=hint) in (sol8, sol8s)
    assert dsolve(eq9, f(x), hint=hint) in (sol9, sol9s)
    assert dsolve(eq10, f(x), hint=hint) in (sol10, sol10s)
    assert dsolve(eq11, f(x), hint=hint+'_Integral') in (sol11, sol11s)
    assert dsolve(eq12, f(x), hint=hint) in (sol12, sol12s)
    assert checkodesol(eq1, f(x), sol1, order=3, solve_for_func=False)[0]
    assert checkodesol(eq2, f(x), sol2, order=3, solve_for_func=False)[0]
    assert checkodesol(eq3, f(x), sol3, order=1, solve_for_func=False)[0]
    assert checkodesol(eq4, f(x), sol4, order=2, solve_for_func=False)[0]
    assert checkodesol(eq5, f(x), sol5, order=2, solve_for_func=False)[0]
    assert checkodesol(eq6, f(x), sol6, order=2, solve_for_func=False)[0]
    assert checkodesol(eq7, f(x), sol7, order=2, solve_for_func=False)[0]
    assert checkodesol(eq8, f(x), sol8, order=2, solve_for_func=False)[0]
    assert checkodesol(eq9, f(x), sol9, order=3, solve_for_func=False)[0]
    assert checkodesol(eq10, f(x), sol10, order=2, solve_for_func=False)[0]
    assert checkodesol(eq12, f(x), sol12, order=4, solve_for_func=False)[0]

@XFAIL # XXX: #@$^&*!
def test_nth_linear_constant_coeff_variation_of_parameters_simplify_False():
    # solve_variation_of_parameters should attempt to simplify the Wronskian
    # if simplify=False.  This test will run considerably slower if this
    # isn't working.
    hint = 'nth_linear_constant_coeff_variation_of_parameters'
    assert len(str(dsolve(f(x).diff(x, 5) + 2*f(x).diff(x, 3) + f(x).diff(x) -
        2*x - exp(I*x), f(x), hint + "_Integral", simplify=False))) == 2522

def test_Liouville_ODE():
    hint = 'Liouville'
    # The first part here used to be test_ODE_1() from test_solvers.py
    eq1 = diff(f(x),x)/x + diff(f(x),x,x)/2 - diff(f(x),x)**2/2
    eq1a = diff(x*exp(-f(x)), x, x)
    eq2 = (eq1*exp(-f(x))/exp(f(x))).expand() # see test_unexpanded_Liouville_ODE() below
    eq3 = diff(f(x), x, x) + 1/f(x)*(diff(f(x), x))**2 + 1/x*diff(f(x), x)
    eq4 = x*diff(f(x), x, x) + x/f(x)*diff(f(x), x)**2 + x*diff(f(x), x)
    eq5 = Eq((x*exp(f(x))).diff(x, x), 0)
    sol1 = Eq(C1 + C2/x - exp(-f(x)), 0)
    # If solve() is ever improved, this is a better solution
    sol1a = Eq(f(x), -log((C1*x+C2)/x))
    sol2 = Eq(C1 + C2/x - exp(-f(x)), 0) # This is equivalent to sol1
    sol3 = [Eq(f(x), -sqrt(C1 + C2*log(x))), Eq(f(x), sqrt(C1 + C2*log(x)))]
    sol4 = [Eq(f(x), sqrt(C1 + C2*exp(-x))), Eq(f(x), -sqrt(C1 + C2*exp(-x)))]
    sol5 = Eq(f(x), -log(x) + log(C1 + C2*x))
    sol1s = ode_renumber(sol1, 'C', 1, 2)
    sol2s = ode_renumber(sol2, 'C', 1, 2)
    sol3s = ode_renumber(sol3, 'C', 1, 2)
    sol4s = ode_renumber(sol4, 'C', 1, 2)
    sol5s = ode_renumber(sol5, 'C', 1, 2)
    assert dsolve(eq1, f(x), hint) in (sol1, sol1s)
    assert dsolve(eq1a, f(x), hint) in (sol1, sol1s)
    assert dsolve(eq2, f(x), hint) in (sol2, sol2s)
    assert dsolve(eq3, f(x), hint) in (sol3, sol3s)
    assert dsolve(eq4, f(x), hint) in (sol4, sol4s)
    assert dsolve(eq5, f(x), hint) in (sol5, sol5s)
    assert checkodesol(sol1, f(x), sol1a, order=2, solve_for_func=False)[0]
    assert checkodesol(eq1, f(x), sol1a, order=2, solve_for_func=False)[0]
    assert checkodesol(eq1a, f(x), sol1a, order=2, solve_for_func=False)[0]
    assert checkodesol(sol2, f(x), sol1a, order=2, solve_for_func=False)[0]
    assert checkodesol(eq2, f(x), sol1a, order=2, solve_for_func=False)[0]
    assert checkodesol(eq3, f(x), sol3[0], order=2, solve_for_func=False)[0]
    assert checkodesol(eq3, f(x), sol3[1], order=2, solve_for_func=False)[0]
    assert checkodesol(eq5, f(x), sol5, order=2, solve_for_func=False)[0]
    not_Liouville1 = classify_ode(diff(f(x),x)/x + f(x)*diff(f(x),x,x)/2 -
        diff(f(x),x)**2/2, f(x))
    not_Liouville2 = classify_ode(diff(f(x),x)/x + diff(f(x),x,x)/2 -
        x*diff(f(x),x)**2/2, f(x))
    assert hint not in not_Liouville1
    assert hint not in not_Liouville2
    assert hint+'_Integral' not in not_Liouville1
    assert hint+'_Integral' not in not_Liouville2


@XFAIL
def test_Liouville_ODE_xfail():
    # This is failling because of how checkodesol is simplifying the
    # differential equation once the solutions is plugged in.
    eq4 = x*diff(f(x), x, x) + x/f(x)*diff(f(x), x)**2 + x*diff(f(x), x)
    sol4 = [Eq(f(x), sqrt(C1 + C2*exp(-x))), Eq(f(x), -sqrt(C1 + C2*exp(-x)))]
    assert checkodesol(eq4, f(x), sol4[0], order=2, solve_for_func=False)[0]
    assert checkodesol(eq4, f(x), sol4[1], order=2, solve_for_func=False)[0]

def test_unexpanded_Liouville_ODE():
    # This is the same as eq1 from test_Liouville_ODE() above.
    eq1 = diff(f(x),x)/x+diff(f(x),x,x)/2- diff(f(x),x)**2/2
    eq2 = eq1*exp(-f(x))/exp(f(x))
    sol2 = Eq(C1 + C2/x - exp(-f(x)), 0)
    sol2s = ode_renumber(sol2, 'C', 1, 2)
    assert dsolve(eq2, f(x)) in (sol2, sol2s)
    assert checkodesol(eq2, f(x), sol2, order=2, solve_for_func=False)[0]

def test_1686():
    from sympy.abc import A
    eq = x + A*(x + diff(f(x), x) + f(x)) + diff(f(x), x) + f(x) + 2
    assert classify_ode(eq, f(x)) == ('1st_linear', \
    'nth_linear_constant_coeff_undetermined_coefficients', \
    'nth_linear_constant_coeff_variation_of_parameters', '1st_linear_Integral', \
    'nth_linear_constant_coeff_variation_of_parameters_Integral')
    # 1765
    eq=(x**2 + f(x)**2)*f(x).diff(x) - 2*x*f(x)
    assert classify_ode(eq, f(x)) == (
        '1st_homogeneous_coeff_best',
        '1st_homogeneous_coeff_subs_indep_div_dep',
        '1st_homogeneous_coeff_subs_dep_div_indep',
        '1st_homogeneous_coeff_subs_indep_div_dep_Integral',
        '1st_homogeneous_coeff_subs_dep_div_indep_Integral')

def test_1726():
    raises(ValueError, "dsolve(f(x, y).diff(x) - y*f(x, y), f(x))")
    assert classify_ode(f(x, y).diff(x) - y*f(x, y), f(x), dict=True) == \
    {'default': None, 'order': 0}
    # See also issue 694, test Z13.
    raises(ValueError, "dsolve(f(x).diff(x), f(y))")
    assert classify_ode(f(x).diff(x), f(y), dict=True) == {'default': None, 'order': 0}

