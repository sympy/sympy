from __future__ import annotations
from sympy.concrete.summations import Sum
from sympy.core.function import (Derivative as D, Function, Subs)
from sympy.core.numbers import I, pi, oo
from sympy.core.relational import Eq
from sympy.core.symbol import (Symbol, symbols)
from sympy.functions.elementary.exponential import (exp, log)
from sympy.functions.elementary.hyperbolic import cosh, sinh
from sympy.functions.elementary.miscellaneous import sqrt
from sympy.functions.elementary.piecewise import Piecewise
from sympy.functions.elementary.trigonometric import (cos, sin)
from sympy.core import S
from sympy.integrals.integrals import Integral
from sympy.series.limits import Limit
from sympy.solvers.pde import (pde_separate, pde_separate_add, pde_separate_mul,
    pdsolve, classify_pde, checkpdesol)
from sympy.testing.pytest import raises


a, b, c, x, y = symbols('a b c x y')

def test_pde_separate_add():
    x, y, z, t = symbols("x,y,z,t")
    F, T, X, Y, Z, u = map(Function, 'FTXYZu')

    eq = Eq(D(u(x, t), x), D(u(x, t), t)*exp(u(x, t)))
    res = pde_separate_add(eq, u(x, t), [X(x), T(t)])
    assert res == [D(X(x), x)*exp(-X(x)), D(T(t), t)*exp(T(t))]


def test_pde_separate():
    x, y, z, t = symbols("x,y,z,t")
    F, T, X, Y, Z, u = map(Function, 'FTXYZu')

    eq = Eq(D(u(x, t), x), D(u(x, t), t)*exp(u(x, t)))
    raises(ValueError, lambda: pde_separate(eq, u(x, t), [X(x), T(t)], 'div'))


def test_pde_separate_mul():
    x, y, z, t = symbols("x,y,z,t")
    c = Symbol("C", real=True)
    Phi = Function('Phi')
    F, R, T, X, Y, Z, u = map(Function, 'FRTXYZu')
    r, theta, z = symbols('r,theta,z')

    # Something simple :)
    eq = Eq(D(F(x, y, z), x) + D(F(x, y, z), y) + D(F(x, y, z), z), 0)

    # Duplicate arguments in functions
    raises(
        ValueError, lambda: pde_separate_mul(eq, F(x, y, z), [X(x), u(z, z)]))
    # Wrong number of arguments
    raises(ValueError, lambda: pde_separate_mul(eq, F(x, y, z), [X(x), Y(y)]))
    # Wrong variables: [x, y] -> [x, z]
    raises(
        ValueError, lambda: pde_separate_mul(eq, F(x, y, z), [X(t), Y(x, y)]))

    assert pde_separate_mul(eq, F(x, y, z), [Y(y), u(x, z)]) == \
        [D(Y(y), y)/Y(y), -D(u(x, z), x)/u(x, z) - D(u(x, z), z)/u(x, z)]
    assert pde_separate_mul(eq, F(x, y, z), [X(x), Y(y), Z(z)]) == \
        [D(X(x), x)/X(x), -D(Z(z), z)/Z(z) - D(Y(y), y)/Y(y)]

    # wave equation
    wave = Eq(D(u(x, t), t, t), c**2*D(u(x, t), x, x))
    res = pde_separate_mul(wave, u(x, t), [X(x), T(t)])
    assert res == [D(X(x), x, x)/X(x), D(T(t), t, t)/(c**2*T(t))]

    # Laplace equation in cylindrical coords
    eq = Eq(1/r * D(Phi(r, theta, z), r) + D(Phi(r, theta, z), r, 2) +
            1/r**2 * D(Phi(r, theta, z), theta, 2) + D(Phi(r, theta, z), z, 2), 0)
    # Separate z
    res = pde_separate_mul(eq, Phi(r, theta, z), [Z(z), u(theta, r)])
    assert res == [D(Z(z), z, z)/Z(z),
            -D(u(theta, r), r, r)/u(theta, r) -
        D(u(theta, r), r)/(r*u(theta, r)) -
        D(u(theta, r), theta, theta)/(r**2*u(theta, r))]
    # Lets use the result to create a new equation...
    eq = Eq(res[1], c)
    # ...and separate theta...
    res = pde_separate_mul(eq, u(theta, r), [T(theta), R(r)])
    assert res == [D(T(theta), theta, theta)/T(theta),
            -r*D(R(r), r)/R(r) - r**2*D(R(r), r, r)/R(r) - c*r**2]
    # ...or r...
    res = pde_separate_mul(eq, u(theta, r), [R(r), T(theta)])
    assert res == [r*D(R(r), r)/R(r) + r**2*D(R(r), r, r)/R(r) + c*r**2,
            -D(T(theta), theta, theta)/T(theta)]


def test_issue_11726():
    x, t = symbols("x t")
    f  = symbols("f", cls=Function)
    X, T = symbols("X T", cls=Function)

    u = f(x, t)
    eq = u.diff(x, 2) - u.diff(t, 2)
    res = pde_separate(eq, u, [T(x), X(t)])
    assert res == [D(T(x), x, x)/T(x),D(X(t), t, t)/X(t)]


def test_pde_classify():
    # When more number of hints are added, add tests for classifying here.
    f = Function('f')
    eq1 = a*f(x,y) + b*f(x,y).diff(x) + c*f(x,y).diff(y)
    eq2 = 3*f(x,y) + 2*f(x,y).diff(x) + f(x,y).diff(y)
    eq3 = a*f(x,y) + b*f(x,y).diff(x) + 2*f(x,y).diff(y)
    eq4 = x*f(x,y) + f(x,y).diff(x) + 3*f(x,y).diff(y)
    eq5 = x**2*f(x,y) + x*f(x,y).diff(x) + x*y*f(x,y).diff(y)
    eq6 = y*x**2*f(x,y) + y*f(x,y).diff(x) + f(x,y).diff(y)
    for eq in [eq1, eq2, eq3]:
        assert classify_pde(eq) == ('1st_linear_constant_coeff_homogeneous',)
    for eq in [eq4, eq5, eq6]:
        assert classify_pde(eq) == ('1st_linear_variable_coeff',)


def test_checkpdesol():
    f, F = map(Function, ['f', 'F'])
    eq1 = a*f(x,y) + b*f(x,y).diff(x) + c*f(x,y).diff(y)
    eq2 = 3*f(x,y) + 2*f(x,y).diff(x) + f(x,y).diff(y)
    eq3 = a*f(x,y) + b*f(x,y).diff(x) + 2*f(x,y).diff(y)
    for eq in [eq1, eq2, eq3]:
        assert checkpdesol(eq, pdsolve(eq))[0]
    eq4 = x*f(x,y) + f(x,y).diff(x) + 3*f(x,y).diff(y)
    eq5 = 2*f(x,y) + 1*f(x,y).diff(x) + 3*f(x,y).diff(y)
    eq6 = f(x,y) + 1*f(x,y).diff(x) + 3*f(x,y).diff(y)
    assert checkpdesol(eq4, [pdsolve(eq5), pdsolve(eq6)]) == [
        (False, (x - 2)*F(3*x - y)*exp(-x/S(5) - 3*y/S(5))),
         (False, (x - 1)*F(3*x - y)*exp(-x/S(10) - 3*y/S(10)))]
    for eq in [eq4, eq5, eq6]:
        assert checkpdesol(eq, pdsolve(eq))[0]
    sol = pdsolve(eq4)
    sol4 = Eq(sol.lhs - sol.rhs, 0)
    raises(NotImplementedError, lambda:
        checkpdesol(eq4, sol4, solve_for_func=False))


def test_solvefun():
    f, F, G, H = map(Function, ['f', 'F', 'G', 'H'])
    eq1 = f(x,y) + f(x,y).diff(x) + f(x,y).diff(y)
    assert pdsolve(eq1) == Eq(f(x, y), F(x - y)*exp(-x/2 - y/2))
    assert pdsolve(eq1, solvefun=G) == Eq(f(x, y), G(x - y)*exp(-x/2 - y/2))
    assert pdsolve(eq1, solvefun=H) == Eq(f(x, y), H(x - y)*exp(-x/2 - y/2))


def test_pde_1st_linear_constant_coeff_homogeneous():
    f, F = map(Function, ['f', 'F'])
    u = f(x, y)
    eq = 2*u + u.diff(x) + u.diff(y)
    assert classify_pde(eq) == ('1st_linear_constant_coeff_homogeneous',)
    sol = pdsolve(eq)
    assert sol == Eq(u, F(x - y)*exp(-x - y))
    assert checkpdesol(eq, sol)[0]

    eq = 4 + (3*u.diff(x)/u) + (2*u.diff(y)/u)
    assert classify_pde(eq) == ('1st_linear_constant_coeff_homogeneous',)
    sol = pdsolve(eq)
    assert sol == Eq(u, F(2*x - 3*y)*exp(-S(12)*x/13 - S(8)*y/13))
    assert checkpdesol(eq, sol)[0]

    eq = u + (6*u.diff(x)) + (7*u.diff(y))
    assert classify_pde(eq) == ('1st_linear_constant_coeff_homogeneous',)
    sol = pdsolve(eq)
    assert sol == Eq(u, F(7*x - 6*y)*exp(-6*x/S(85) - 7*y/S(85)))
    assert checkpdesol(eq, sol)[0]

    eq = a*u + b*u.diff(x) + c*u.diff(y)
    sol = pdsolve(eq)
    assert checkpdesol(eq, sol)[0]


def test_pde_1st_linear_constant_coeff():
    f, F = map(Function, ['f', 'F'])
    u = f(x,y)
    eq = -2*u.diff(x) + 4*u.diff(y) + 5*u - exp(x + 3*y)
    sol = pdsolve(eq)
    assert sol == Eq(f(x,y),
    (F(4*x + 2*y)*exp(x/2) + exp(x + 4*y)/15)*exp(-y))
    assert classify_pde(eq) == ('1st_linear_constant_coeff',
    '1st_linear_constant_coeff_Integral')
    assert checkpdesol(eq, sol)[0]

    eq = (u.diff(x)/u) + (u.diff(y)/u) + 1 - (exp(x + y)/u)
    sol = pdsolve(eq)
    assert sol == Eq(f(x, y), F(x - y)*exp(-x/2 - y/2) + exp(x + y)/3)
    assert classify_pde(eq) == ('1st_linear_constant_coeff',
    '1st_linear_constant_coeff_Integral')
    assert checkpdesol(eq, sol)[0]

    eq = 2*u + -u.diff(x) + 3*u.diff(y) + sin(x)
    sol = pdsolve(eq)
    assert sol == Eq(f(x, y),
         F(3*x + y)*exp(x/5 - 3*y/5) - 2*sin(x)/5 - cos(x)/5)
    assert classify_pde(eq) == ('1st_linear_constant_coeff',
    '1st_linear_constant_coeff_Integral')
    assert checkpdesol(eq, sol)[0]

    eq = u + u.diff(x) + u.diff(y) + x*y
    sol = pdsolve(eq)
    assert sol.expand() == Eq(f(x, y),
        x + y + (x - y)**2/4 - (x + y)**2/4 + F(x - y)*exp(-x/2 - y/2) - 2).expand()
    assert classify_pde(eq) == ('1st_linear_constant_coeff',
    '1st_linear_constant_coeff_Integral')
    assert checkpdesol(eq, sol)[0]
    eq = u + u.diff(x) + u.diff(y) + log(x)
    assert classify_pde(eq) == ('1st_linear_constant_coeff',
    '1st_linear_constant_coeff_Integral')


def test_pdsolve_all():
    f, F = map(Function, ['f', 'F'])
    u = f(x,y)
    eq = u + u.diff(x) + u.diff(y) + x**2*y
    sol = pdsolve(eq, hint = 'all')
    keys = ['1st_linear_constant_coeff',
        '1st_linear_constant_coeff_Integral', 'default', 'order']
    assert sorted(sol.keys()) == keys
    assert sol['order'] == 1
    assert sol['default'] == '1st_linear_constant_coeff'
    assert sol['1st_linear_constant_coeff'].expand() == Eq(f(x, y),
        -x**2*y + x**2 + 2*x*y - 4*x - 2*y + F(x - y)*exp(-x/2 - y/2) + 6).expand()


def test_pdsolve_variable_coeff():
    f, F = map(Function, ['f', 'F'])
    u = f(x, y)
    eq = x*(u.diff(x)) - y*(u.diff(y)) + y**2*u - y**2
    sol = pdsolve(eq, hint="1st_linear_variable_coeff")
    assert sol == Eq(u, F(x*y)*exp(y**2/2) + 1)
    assert checkpdesol(eq, sol)[0]

    eq = x**2*u + x*u.diff(x) + x*y*u.diff(y)
    sol = pdsolve(eq, hint='1st_linear_variable_coeff')
    assert sol == Eq(u, F(y*exp(-x))*exp(-x**2/2))
    assert checkpdesol(eq, sol)[0]

    eq = y*x**2*u + y*u.diff(x) + u.diff(y)
    sol = pdsolve(eq, hint='1st_linear_variable_coeff')
    assert sol == Eq(u, F(-2*x + y**2)*exp(-x**3/3))
    assert checkpdesol(eq, sol)[0]

    eq = exp(x)**2*(u.diff(x)) + y
    sol = pdsolve(eq, hint='1st_linear_variable_coeff')
    assert sol == Eq(u, y*exp(-2*x)/2 + F(y))
    assert checkpdesol(eq, sol)[0]

    eq = exp(2*x)*(u.diff(y)) + y*u - u
    sol = pdsolve(eq, hint='1st_linear_variable_coeff')
    assert sol == Eq(u, exp((-y**2 + 2*y + 2*F(x))*exp(-2*x)/2))


def test_pde_2nd_linear_constant_coeff_separable_exact():
    k, H, L = symbols("k H L", positive=True)
    u = Function("u")

    # dummies
    F1 = Function("F1")
    n = Symbol("n", integer=True, positive=True)
    s = Symbol("s", real=True)

    # 1D Heat, eigenfunction conditions
    eq = u(x, y).diff(x, 2) - k * u(x, y).diff(y)
    conds = {u(0, y): 0, u(L, y): 0, u(x, 0): sin(pi*x/L) - 2*sin(2*pi*x/L)}
    assert classify_pde(eq, ics=conds)[0] == "2nd_linear_constant_coeff_separable"
    expected = Eq(u(x, y), Piecewise((exp(-pi**2*y/(L**2*k))*sin(pi*x/L) - 2*exp(-4*pi**2*y/(L**2*k))*sin(2*pi*x/L), (L >= x) & (x >= 0) & (y >= 0))))
    received = pdsolve(eq, u(x, y), ics=conds)
    assert received.dummy_eq(expected)
    assert checkpdesol(eq, Eq(received.lhs, received.rhs.args[0].expr))[0]

    # Wave with u_y term, eigenfunction conditions
    eq = u(x, y).diff(x, 2) - u(x, y).diff(y, 2) - 2*u(x, y).diff(y)
    conds = {u(0, y): 0, u(pi, y): 0, u(x, 0): sin(x) + sin(2*x), u(x, y).diff(y).subs(y, 0): 0}
    assert classify_pde(eq, ics=conds)[0] == "2nd_linear_constant_coeff_separable"
    expected = Eq(u(x, y), Piecewise(((y + 1)*exp(-y)*sin(x) + (sqrt(3)*sin(sqrt(3)*y)/3 + cos(sqrt(3)*y))*exp(-y)*sin(2*x), (x >= 0) & (x <= pi))))
    received = pdsolve(eq, u(x, y), ics=conds)
    assert received.dummy_eq(expected)
    assert checkpdesol(eq, Eq(received.lhs, received.rhs.args[0].expr))[0]

    # Laplace, symbolic bounds
    eq = u(x, y).diff(x, 2) + u(x, y).diff(y, 2)
    conds = {u(a, y): 0, u(a + L, y): 0, u(x, b): 1, u(x, y).diff(y).subs(y, b + H): 0}
    assert classify_pde(eq, ics=conds)[0] == "2nd_linear_constant_coeff_separable"
    expected = Eq(u(x, y), Piecewise((Sum(-2*sin(n*pi*(-a + x)/L)*((-1)**n - 1)*cosh(n*pi*(-H - b + y)/L)/(n*pi*cosh(n*pi*H/L)), (n, 1, oo)), (a <= x) & (b <= y) & (x <= L + a) & (y <= H + b))))
    received = pdsolve(eq, u(x, y), ics=conds)
    assert received.dummy_eq(expected)
    # checkpdesol stalls

    # Schrodinger's, missing conditions
    eq = u(x, y).diff(y) - I * u(x, y).diff(x, 2)
    conds = {u(x, y).diff(x).subs(x, 0): 0, u(x, y).diff(x).subs(x, pi): 0}
    assert classify_pde(eq, ics=conds)[0] == "2nd_linear_constant_coeff_separable"
    expected = Eq(u(x, y), Piecewise((Integral(F1(s), (s, 0, pi))/pi + Sum(2*exp(-n**2*I*y)*cos(n*x)*Integral(F1(s)*cos(n*s), (s, 0, pi))/pi, (n, 1, oo)), (x >= 0) & (x <= pi))))
    received = pdsolve(eq, u(x, y), ics=conds)
    assert received.dummy_eq(expected)
    assert checkpdesol(eq, Eq(received.lhs, received.rhs.args[0].expr))[0]

    # Wave w/ three-range solution
    eq = u(x, y).diff(x, 2) - u(x, y).diff(y, 2) - 10*u(x, y).diff(y)
    conds = {u(0, y): 0, u(pi, y): 0, u(x, 0): x*(pi-x), u(x, y).diff(y).subs(y, 0): 0}
    assert classify_pde(eq, ics=conds)[0] == "2nd_linear_constant_coeff_separable"
    expected = Eq(u(x, y), Piecewise((8*exp(-5*y)*(5*y + 1)*sin(5*x)/(125*pi) + Sum(-4*exp(-5*y)*((-1)**n - 1)*(cos(y*sqrt(n**2 - 25)) + 5*sin(y*sqrt(n**2 - 25))/sqrt(n**2 - 25))*sin(n*x)/(n**3*pi), (n, 6, oo)) + Sum(-4*exp(-5*y)*((-1)**n - 1)*(cosh(y*sqrt(25 - n**2)) + 5*sinh(y*sqrt(25 - n**2))/sqrt(25 - n**2))*sin(n*x)/(n**3*pi), (n, 1, 4)), (x >= 0) & (x <= pi))))
    received = pdsolve(eq, u(x, y), ics=conds)
    assert received.dummy_eq(expected)
    # checkpdesol stalls

    # Backwards Heat, special first term
    eq = u(x, y).diff(x, 2) + k * u(x, y).diff(y)
    conds = {u(0, y): 0, u(1, y): 0, u(x, 0): x*sin(pi*x)}
    assert classify_pde(eq, ics=conds)[0] == "2nd_linear_constant_coeff_separable"
    expected = Eq(u(x, y), Piecewise((exp(pi**2*y/k)*sin(pi*x)/2 + Sum(-4*n*((-1)**n + 1)*exp(n**2*pi**2*y/k)*sin(n*pi*x)/(pi**2*(n**4 - 2*n**2 + 1)), (n, 2, oo)), (x >= 0) & (x <= 1) & (y <= 0))))
    received = pdsolve(eq, u(x, y), ics=conds)
    assert received.dummy_eq(expected)
    # checkpdesol stalls

def test_pde_2nd_linear_constant_coeff_separable_value_errors():
    k, H, L = symbols("k H L", positive=True)
    m = Symbol("m", integer=True)
    f, g = symbols("f g", cls=Function)
    u = Function("u")
    u_x, u_y = u(x, y).diff(x), u(x, y).diff(y)
    u_xx, u_yy, u_xy = u_x.diff(x), u_y.diff(y), u_x.diff(y)

    BCs = {u(0, y): 0, u(L, y): 0}

    heat = u_xx - k * u_y
    heat_conds = BCs | {u(x, 0): a}
    wave = u_xx - k * u_yy
    wave_conds = heat_conds | {u_y.subs(y, 0): g(x)}
    elliptic = u_xx + u_yy + u(x, y)
    elliptic_conds = heat_conds | {u(x, H): g(x)}

    # confirm no errors
    pdsolve(heat, u(x, y), ics=heat_conds)
    pdsolve(wave, u(x, y), ics=wave_conds)
    pdsolve(elliptic, u(x, y), ics=elliptic_conds)


    # Invalid Individual Condition
    two_fixed_vars = BCs | {u(0, 0): a}
    raises(ValueError, lambda: pdsolve(heat, u(x, y), ics=two_fixed_vars))

    two_fixed_vars_subs = BCs | {Subs(u(x, y), (x, y), (0, 0)): a}
    raises(ValueError, lambda: pdsolve(heat, u(x, y), ics=two_fixed_vars_subs))

    nonlinear = BCs | {Subs(u(x, y)**2, y, 0): a}
    raises(ValueError, lambda: pdsolve(heat, u(x, y), ics=nonlinear))

    wrong_diff = BCs | {u_x.subs(y, 0): a}
    raises(ValueError, lambda: pdsolve(heat, u(x, y), ics=wrong_diff))

    mixed_derivative_numeric = BCs | {u_xy.subs(y, 0): a}
    raises(ValueError, lambda: pdsolve(heat, u(x, y), ics=mixed_derivative_numeric))

    mixed_derivative_symbolic = BCs | {u_xy.subs(y, H): a}
    raises(ValueError, lambda: pdsolve(heat, u(x, y), ics=mixed_derivative_symbolic))

    different_function = BCs | {f(x, 0): a}
    raises(ValueError, lambda: pdsolve(heat, u(x, y), ics=different_function))

    fixed_variable_value = BCs | {u(x, 0): y}
    raises(ValueError, lambda: pdsolve(heat, u(x, y), ics=fixed_variable_value))


    # Invalid Variable/Coefficient
    discrete_variable = heat.subs(x, m)
    discrete_variable_conds = BCs | {u(m, 0): a}
    raises(ValueError, lambda: pdsolve(discrete_variable, u(m, y), ics=discrete_variable_conds))

    ode_in_disguise = u_xx + u(x, y)
    raises(ValueError, lambda: pdsolve(ode_in_disguise, u(x, y), ics=heat_conds))

    ambiguous_coeff_sign = wave.subs(k, a) # `a` can be pos (wave) or neg (laplace)
    raises(ValueError, lambda: pdsolve(ambiguous_coeff_sign, u(x, y), ics=wave_conds))


    # Invalid Condition Arrangement
    raises(ValueError, lambda: pdsolve(heat, u(x, y), ics=wave_conds))

    order_too_high = BCs | {u_y.subs(y, 0): a}
    raises(ValueError, lambda: pdsolve(heat, u(x, y), ics=order_too_high))

    unorderable_points = {u(H, y): 0, u(L, y): 0, u(x, 0): a}
    raises(ValueError, lambda: pdsolve(heat, u(x, y), ics=unorderable_points))

    time_direction_unknown = heat.subs(k, a)
    raises(ValueError, lambda: pdsolve(time_direction_unknown, u(x, y), ics=heat_conds))


    # Invalid Arrangements by PDE Type
    raises(ValueError, lambda: pdsolve(wave, u(x, y), ics=elliptic_conds))
    raises(ValueError, lambda: pdsolve(elliptic, u(x, y), ics=wave_conds))

    laplace = elliptic - u(x, y)
    laplace_all_neumann = {u_x.subs(x, 0): 0, u_x.subs(x, L): 0, u_y.subs(y, 0): 0, u_y.subs(y, H): g(x)}
    raises(ValueError, lambda: pdsolve(laplace, u(x, y), ics=laplace_all_neumann))

    timelike_heat = {u(0, y): 0, u_x.subs(x, 0): 0, u(x, 0): a}
    raises(ValueError, lambda: pdsolve(heat, u(x, y), ics=timelike_heat))

    hyperbolic_needing_two_timelike = wave.subs(k, -1) + 3*u_xy
    raises(ValueError, lambda: pdsolve(hyperbolic_needing_two_timelike, u(x, y), ics=wave_conds))

    timelike_wave = {u(0, y): 0, u_x.subs(x, 0): 0, u(x, 0): a, u_y.subs(y, 0): g(x)}
    raises(ValueError, lambda: pdsolve(wave, u(x, y), ics=timelike_wave))

    insufficient = {u(x, 0): x}
    raises(ValueError, lambda: pdsolve(wave, u(x, y), ics=insufficient))


    # Resonance
    resonant_family = {u(0, y): 0, u(5*pi, y): 0, u(x, 0): a, u(x, 5*pi): 0}
    raises(ValueError, lambda: pdsolve(elliptic, u(x, y), ics=resonant_family))

    resonant_neumann = {u(0, y): 0, u(pi, y): 0, u_y.subs(y, 0): sin(x), u_y.subs(y, 1): 0}
    raises(ValueError, lambda: pdsolve(elliptic, u(x, y), ics=resonant_neumann))

    resonant_eigenvalue = {u_x.subs(x, 0): 0, u_x.subs(x, pi): 0, u(x, 0): a, u(x, pi): g(x)}
    raises(ValueError, lambda: pdsolve(elliptic, u(x, y), ics=resonant_eigenvalue))

    # Should be NotImplementedError
    robin_cond = {u(0, y) + u_x.subs(x, 0): 0, u(L, y): 0, u(x, 0): 1}
    raises(ValueError, lambda: pdsolve(heat, u(x, y), ics=robin_cond))


def test_pde_2nd_linear_constant_coeff_separable_not_implemented_errors():
    k, L = symbols("k L", positive=True)
    u = Function("u")
    u_x, u_y = u(x, y).diff(x), u(x, y).diff(y)
    u_xx, u_yy, u_xy = u_x.diff(x), u_y.diff(y), u_x.diff(y)

    BCs = {u(0, y): 0, u(L, y): 0}
    heat = u_xx - k * u_y

    periodic_cond = BCs | {u(x, 0): u(x, pi)}
    raises(NotImplementedError, lambda: pdsolve(heat, u(x, y), ics=periodic_cond))

    limit_cond = BCs | {Limit(u(x, y), y, oo): 0}
    raises(NotImplementedError, lambda: pdsolve(heat, u(x, y), ics=limit_cond))

    infty_point = BCs | {u(x, oo): 0}
    raises(NotImplementedError, lambda: pdsolve(heat, u(x, y), ics=infty_point))

    infty_value = BCs | {u(x, 0): oo}
    raises(NotImplementedError, lambda: pdsolve(heat, u(x, y), ics=infty_value))

    z = Symbol("z", real=False, complex=True)
    complex_heat = heat.subs(x, z)
    complex_cond = BCs | {u(z, 0): 0}
    raises(NotImplementedError, lambda: pdsolve(complex_heat, u(z, y), ics=complex_cond))

    nonlinear = u_xx*u_x + u_yy
    raises(NotImplementedError, lambda: pdsolve(nonlinear, u(x, y), ics=BCs))

    variable_coeffs = u_xx - u_x*u_y
    raises(NotImplementedError, lambda: pdsolve(variable_coeffs, u(x, y), ics=BCs))

    complex_higher_order = heat * I + u_yy * (1 + I)
    raises(NotImplementedError, lambda: pdsolve(complex_higher_order, u(x, y), ics=BCs))

    first_order = u_xy + u(x, y)
    raises(NotImplementedError, lambda: pdsolve(first_order, u(x, y), ics=BCs))

    parabolic_mixed_term = u_xx + 2*u_xy + u_yy
    raises(NotImplementedError, lambda: pdsolve(parabolic_mixed_term, u(x, y), ics=BCs))

    implicit_eigen_pde = heat + u_x
    implicit_eigen_conds = {u(0, y): 0, u_x.subs(x, L): 0, u(x, 0): a}
    raises(NotImplementedError, lambda: pdsolve(implicit_eigen_pde, u(x, y), ics=implicit_eigen_conds))
