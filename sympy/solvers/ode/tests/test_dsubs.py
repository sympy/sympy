from sympy import (S, symbols, Eq, Function, Derivative, dsolve, dsubs, \
    checkodesol, numer, asin, sin, cos, tan, exp, log, sqrt)
from sympy.abc import a, b, m, n, t, u, w, x, y, E

f, g, h, F, G, H, psi = symbols('f g h F G H psi', cls=Function)
hbar, C1, C2 = symbols("hbar C1 C2")

def test_dsubs():
    # Forward substitution

    # Order 1
    eq = f(x).diff(x) - f(x)
    # Issue 17950
    assert dsubs(eq, {x: exp(y)}) == -f(exp(y)) + exp(-y)*Derivative(f(exp(y)), y)
    assert dsubs(eq, {x: exp(y), f(x): g(y)}) == -g(y) + exp(-y)*g(y).diff(y)
    assert dsubs(eq, {x: exp(y), f(x): tan(g(y))}) == (tan(g(y))**2 + 1)*exp(-y)*g(y).diff(y) - tan(g(y))

    # Order 2
    eq = f(x).diff(x, 2)
    assert dsubs(eq, {x: a*y, f(x): g(y)}, [y]) == g(y).diff(y, 2)/a**2

    assert dsubs(eq, {x: u**3, f(x): g(u)}) == (g(u).diff(u, 2)/(3*u**2) - 2*g(u).diff(u)/(3*u**3))/(3*u**2)

    assert dsubs(eq, {x: t**2, f(x): h(t)**3}) == (3*h(t)**2*h(t).diff(t, 2)/(2*t) + \
        3*h(t)*h(t).diff(t)**2/t - 3*h(t)**2*h(t).diff(t)/(2*t**2))/(2*t)

    assert dsubs(eq, {x: asin(t**3), f(x): g(t)}) == sqrt(1 - t**6)*(-t**3*g(t).diff(t)/sqrt(1 - t**6) \
        + sqrt(1 - t**6)*g(t).diff(t, 2)/(3*t**2) - 2*sqrt(1 - t**6)*g(t).diff(t)/(3*t**3))/(3*t**2)

    # Mailing List https://groups.google.com/g/sympy/c/mj9JRfIafrc/m/S6O_joOgAwAJ
    eq = Eq(n**2*f(x) - x*f(x).diff(x) + (1 - 2*x)*f(x).diff(x, 2), 0)
    assert dsubs(eq, {x: cos(t), f(x): g(t)}) == Eq(n**2*g(t) - (1 - 2*cos(t))*(-g(t).diff(t, 2)/sin(t) + \
        cos(t)*g(t).diff(t)/sin(t)**2)/sin(t) + cos(t)*g(t).diff(t)/sin(t), 0)

    # https://stackoverflow.com/questions/57840957/differential-equation-change-of-variables-with-sympy
    eq = -hbar**2*psi(x).diff(x, 2)/(2*m) + m*w**2*x**2*psi(x)/2 - E*psi(x)
    assert dsubs(eq, {x: u*sqrt(hbar/(m*w)), psi(x): H(u)*exp(-u*u/2)}, [u]).expand() == \
        -E*H(u)*exp(-u**2/2) + hbar*u*w*exp(-u**2/2)*H(u).diff(u) + hbar*w*H(u)*exp(-u**2/2)/2 \
        - hbar*w*exp(-u**2/2)*H(u).diff(u, 2)/2

    # Order 3
    eq = x*f(x).diff(x, 3) + x**(S(5)/2)*f(x).diff(x, 2) + f(x)*f(x).diff(x) + x**2
    assert dsubs(eq, {x: t, f(x): g(t)**2}) == t**(S(5)/2)*(2*g(t)*g(t).diff(t, 2) + 2*g(t).diff(t)**2) \
        + t**2 + t*(2*g(t)*g(t).diff(t, 3) + 6*g(t).diff(t)*g(t).diff(t, 2)) + \
        2*g(t)**3*g(t).diff(t)

    # Order 4
    eq = f(x).diff(x, 4) + f(x)*x*f(x).diff(x, 3) + x**(S(5)/2)*f(x).diff(x, 2)*f(x).diff(x)
    assert dsubs(eq, {x: sqrt(t), f(x): g(t)}).expand() == 8*t**(S(11)/4)*g(t).diff(t)*g(t).diff(t, 2) + \
        4*t**(S(7)/4)*g(t).diff(t)**2 + 8*t**2*g(t)*g(t).diff(t, 3) + 16*t**2*g(t).diff(t, 4) \
        + 12*t*g(t)*g(t).diff(t, 2) + 48*t*g(t).diff(t, 3) + 12*g(t).diff(t, 2)

    # Change multiple variables (Not to be confused with PDE substitution
    eq = f(x).diff(x, 2) + g(y).diff(y) + h(t).diff(t, 3)
    assert dsubs(eq, {t: log(u), x: a**2, y: sqrt(b), h(t): exp(F(u)), f(x): H(a), g(y): G(b)}) == \
        2*sqrt(b)*G(b).diff(b) + u*(u*(u*exp(F(u))*F(u).diff(u)**3 + 3*u*exp(F(u))*F(u).diff(u)*F(u).diff(u, 2) + \
        u*exp(F(u))*F(u).diff(u, 3) + 2*exp(F(u))*F(u).diff(u)**2 + 2*exp(F(u))*F(u).diff(u, 2)) + \
        u*exp(F(u))*F(u).diff(u)**2 + u*exp(F(u))*F(u).diff(u, 2) + exp(F(u))*F(u).diff(u)) + \
        (H(a).diff(a, 2)/(2*a) - H(a).diff(a)/(2*a**2))/(2*a)

    # Reverse Substitution

    eq = f(x).diff(x) - f(x)
    eqt = dsubs(eq, {x: exp(y), f(x): g(y)})
    solt = dsolve(eqt)
    sol = dsubs(solt, {y: log(x), g(y): f(x)})
    assert sol == Eq(f(x), C1*exp(x))
    assert checkodesol(eq, sol)

    eq = f(x).diff(x) - 2*f(x)/x + x**2*f(x)**2
    eqt = numer(dsubs(eq, {f(x): 1/g(t), x: t}).together())
    solt = dsolve(eqt)
    sol = dsubs(solt, {g(t): 1/f(x), t: x})
    assert sol == Eq(1/f(x), C1/x**2 + x**3/5)
    assert checkodesol(eq, sol)

    eq = x**2*f(x).diff(x, 2) + 10*x*f(x).diff(x) + 20*f(x)
    eqt = dsubs(eq, {x: exp(t), f(x): g(t)}).expand()
    solt = dsolve(eqt)
    sol = dsubs(solt, {g(t): f(x), t: log(x)})
    assert sol == Eq(f(x), (C1 + C2/x)/x**4)
    assert checkodesol(eq, sol)
