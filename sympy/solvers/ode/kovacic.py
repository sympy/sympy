from .riccati import solve_riccati

from sympy import S, Eq, exp, Integral, Derivative, Function
from sympy.core.symbol import symbols, Wild
from sympy.solvers.ode.ode import get_numbered_constants, _remove_redundant_solutions, constantsimp
from sympy.solvers.ode.subscheck import checkodesol

def match_2nd_order(eq, z, x):
    a = Wild('a', exclude=[z, z.diff(x), z.diff(x)*2])
    b = Wild('b', exclude=[z, z.diff(x), z.diff(x)*2])
    c = Wild('c', exclude=[z, z.diff(x), z.diff(x)*2])
    match = eq.match(a*z.diff(x, 2) + b*z.diff(x) + c*z)

    if a not in match or b not in match or c not in match:
        raise ValueError("Invalid Second Order Linear Homogeneous Equation")
    return match[a], match[b], match[c]

def find_kovacic_simple(x, a, b, c):
    # Find solution for y(x).diff(x, 2) = r(x)*y(x)
    r = (b**2 + 2*a*b.diff(x) - 2*a.diff(x)*b - 4*a*c)/4*a**2

    # Riccati equation to be solved is g(x).diff(x) + g(x)**2 - r
    # Comparing with f(x).diff(x) = b_0(x) + b_1(x)*f(x) + b_2(x)*f(x)**2
    # we can see that b_0 = r, b_1 = 0, b_2 = -1
    b0, b1, b2 = r, S(0), -S(1)
    g = Function('g')
    ric_sol = [sol.rhs for sol in solve_riccati(g(x), x, b0, b1, b2)]

    C1 = symbols('C1')
    return set(map(lambda sol: exp(Integral(sol.subs(C1, 0), x).doit()), ric_sol))

def find_kovacic_sol(eq):
    z = list(eq.atoms(Derivative))[0].args[0]
    x = list(z.free_symbols)[0]
    eq = eq.expand().collect(z)
    a, b, c = match_2nd_order(eq, z, x)

    # Transform the differential equation to a simpler form
    # using z(x) = y(x)*exp(Integral(-b/(2*a))) and find its solution
    ysol = find_kovacic_simple(x, a, b, c)

    zsol = list(map(lambda sol: sol*exp(Integral(-b/(2*a), x).doit()), ysol))
    zsol = _remove_redundant_solutions(Eq(eq, 0), list(map(lambda sol: Eq(z, sol), zsol)), 2, x)

    C1, C2 = symbols('C1 C2')
    print(zsol)
    if len(zsol) == 2:
        return constantsimp(Eq(z, C1*zsol[0].rhs + C2*zsol[1].rhs), [C1, C2])
    if len(zsol) == 1:
        sol1 = zsol[0].rhs
        sol2 = sol1*Integral(exp(Integral(-b, x).doit())/sol1**2, x)
        zsol.append(Eq(z, sol2))
        return constantsimp(Eq(z, C1*zsol[0].rhs + C2*zsol[1].rhs), [C1, C2])
    return zsol

def test_kovacic_sol():
    f = Function('f')
    x = symbols('x')

    eq = f(x).diff(x, 2) - 2*f(x)/x**2
    sol = find_kovacic_sol(eq)
    assert checkodesol(eq, sol)

    eq = f(x).diff(x, 2) + (3/(16*x**2*(x - 1)**2))*f(x)
    sol = find_kovacic_sol(eq)
    assert checkodesol(eq, sol)

