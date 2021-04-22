from sympy import (S, symbols, Eq, Function, Derivative, asin, sin, cos, tan, exp, log, sqrt)
from sympy.abc import a, b, m, n, t, u, w, x, y, E

f, g, h, F, psi, H = symbols('f g h F psi H', cls=Function)
hbar = symbols("hbar")

def test_dsubs():
    # Order 1
    eq = f(x).diff(x) - f(x)
    # Issue 17950
    assert eq.dsubs({x: exp(y)}) == -f(y) + exp(-y)*Derivative(f(y), y)
    assert eq.dsubs({x: exp(y), f(x): tan(g(y))}) == (tan(g(y))**2 + 1)*exp(-y)*Derivative(g(y), y) - tan(g(y))

    # Order 2
    eq = f(x).diff(x, 2)
    assert eq.dsubs({x: a*y}, [y]) == f(y).diff(y, 2)/a**2

    assert eq.dsubs({x: u**3}) == (f(u).diff(u, 2)/(3*u**2) - 2*f(u).diff(u)/(3*u**3))/(3*u**2)

    assert eq.dsubs({x: t**2, f(x): h(t)**3}) == (3*h(t)**2*h(t).diff(t, 2)/(2*t) + \
        3*h(t)*h(t).diff(t)**2/t - 3*h(t)**2*h(t).diff(t)/(2*t**2))/(2*t)

    assert eq.dsubs({x: asin(t**3)}) == sqrt(1 - t**6)*(-t**3*f(t).diff(t)/sqrt(1 - t**6) \
        + sqrt(1 - t**6)*f(t).diff(t, 2)/(3*t**2) - 2*sqrt(1 - t**6)*f(t).diff(t)/(3*t**3))/(3*t**2)

    # Mailing List https://groups.google.com/g/sympy/c/mj9JRfIafrc/m/S6O_joOgAwAJ
    eq = Eq(n**2*f(x) - x*Derivative(f(x), x) + (1 - 2*x)*Derivative(f(x), (x, 2)), 0)
    assert eq.dsubs({x: cos(t)}) == Eq(n**2*f(t) - (1 - 2*cos(t))*(-f(t).diff(t, 2)/sin(t) + \
        cos(t)*f(t).diff(t)/sin(t)**2)/sin(t) + cos(t)*f(t).diff(t)/sin(t), 0)

    # https://stackoverflow.com/questions/57840957/differential-equation-change-of-variables-with-sympy
    eq = -hbar**2*psi(x).diff(x, 2)/(2*m) + m*w**2*x**2*psi(x)/2 - E*psi(x)
    assert eq.dsubs({x: u*sqrt(hbar/(m*w)), psi(x): H(u)*exp(-u*u/2)}, [u]).expand() == \
        -E*H(u)*exp(-u**2/2) + hbar*u*w*exp(-u**2/2)*Derivative(H(u), u) + hbar*w*H(u)*exp(-u**2/2)/2 \
        - hbar*w*exp(-u**2/2)*Derivative(H(u), (u, 2))/2

    # Order 3
    eq = x*f(x).diff(x, 3) + x**(S(5)/2)*f(x).diff(x, 2) + f(x)*f(x).diff(x) + x**2
    assert eq.dsubs({x: t, f(x): g(t)**2}) == t**(S(5)/2)*(2*g(t)*Derivative(g(t), (t, 2)) + 2*Derivative(g(t), t)**2) \
        + t**2 + t*(2*g(t)*Derivative(g(t), (t, 3)) + 6*Derivative(g(t), t)*Derivative(g(t), (t, 2))) + \
        2*g(t)**3*Derivative(g(t), t)

    # Order 4
    eq = f(x).diff(x, 4) + f(x)*x*f(x).diff(x, 3) + x**(S(5)/2)*f(x).diff(x, 2)*f(x).diff(x)
    assert eq.dsubs({x: sqrt(t)}).expand() == 8*t**(S(11)/4)*Derivative(f(t), t)*Derivative(f(t), (t, 2)) + \
        4*t**(S(7)/4)*Derivative(f(t), t)**2 + 8*t**2*f(t)*Derivative(f(t), (t, 3)) + 16*t**2*Derivative(f(t), (t, 4)) \
        + 12*t*f(t)*Derivative(f(t), (t, 2)) + 48*t*Derivative(f(t), (t, 3)) + 12*Derivative(f(t), (t, 2))

    # Change multiple variables (Not to be confused with PDE substitution, see docstring of dsubs)
    eq = f(x).diff(x, 2) + g(y).diff(y) + h(t).diff(t, 3)
    assert eq.dsubs({t: log(u), x: a**2, y: sqrt(b), h(t): exp(F(u))}) == 2*sqrt(b)*g(b).diff(b) + \
        u*(u*(u*exp(F(u))*F(u).diff(u)**3 + 3*u*exp(F(u))*F(u).diff(u)*F(u).diff(u, 2) + \
        u*exp(F(u))*F(u).diff(u, 3) + 2*exp(F(u))*F(u).diff(u)**2 + 2*exp(F(u))*F(u).diff(u, 2)) + \
        u*exp(F(u))*F(u).diff(u)**2 + u*exp(F(u))*F(u).diff(u, 2) + exp(F(u))*F(u).diff(u)) + \
        (f(a).diff(a, 2)/(2*a) - f(a).diff(a)/(2*a**2))/(2*a)
