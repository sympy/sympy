from sympy import (Symbol, Function, Derivative, Eq, cos, sin)
from sympy.utilities.pytest import raises
from sympy.calculus.euler import euler_equations

def test_euler_interface():
    x = Function('x')
    y = Symbol('y')
    t = Symbol('t')
    raises(TypeError, lambda: euler_equations())
    raises(TypeError, lambda: euler_equations(x(t).diff(t)*y(t), [x(t), y]))
    raises(ValueError, lambda: euler_equations(x(t).diff(t)*x(y), [x(t), x(y)]))
    raises(TypeError, lambda: euler_equations(x(t).diff(t)**2, x(0)))


def test_euler_pendulum():
    x = Function('x')
    t = Symbol('t')
    L = (x(t).diff(t))**2/2 + cos(x(t))
    assert euler_equations(L, x(t), t) == \
        set([Eq(-sin(x(t)) - Derivative(x(t), t, t), 0)])


def test_euler_henonheiles():
    x = Function('x')
    y = Function('y')
    t = Symbol('t')
    L = sum((z(t).diff(t))**2/2 - z(t)**2/2 for z in [x, y])
    L += -x(t)**2*y(t) + y(t)**3/3
    assert euler_equations(L, [x(t), y(t)], t) == \
        set([Eq(-x(t)**2 + y(t)**2 - y(t) - Derivative(y(t), t, t), 0),
             Eq(-2*x(t)*y(t) - x(t) - Derivative(x(t), t, t), 0)])


def test_euler_sineg():
    psi = Function('psi')
    t = Symbol('t')
    x = Symbol('x')
    L = (psi(t, x).diff(t))**2/2 - (psi(t, x).diff(x))**2/2 + cos(psi(t, x))
    assert euler_equations(L, psi(t, x), [t, x]) == \
        set([Eq(-sin(psi(t, x)) - Derivative(psi(t, x), t, t) + \
             Derivative(psi(t, x), x, x), 0)])


def test_euler_high_order():
    # an example from hep-th/0309038
    m = Symbol('m')
    k = Symbol('k')
    x = Function('x')
    y = Function('y')
    t = Symbol('t')
    L = m*Derivative(x(t), t)**2/2 + m*Derivative(y(t), t)**2/2 - \
        k*Derivative(x(t), t)*Derivative(y(t), t, t) + \
        k*Derivative(y(t), t)*Derivative(x(t), t, t)
    assert euler_equations(L, [x(t), y(t)]) == \
        set([Eq(-2*k*Derivative(x(t), t, t, t) - \
        m*Derivative(y(t), t, t), 0), Eq(2*k*Derivative(y(t), t, t, t) - \
        m*Derivative(x(t), t, t), 0)])

    w = Symbol('w')
    L = x(t, w).diff(t, w)**2/2
    assert euler_equations(L) == \
        set([Eq(Derivative(x(t, w), t, t, w, w), 0)])
