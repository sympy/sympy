from sympy import Symbol, Function, Derivative as D, Eq, cos, sin
from sympy.utilities.pytest import raises
from sympy.calculus.hamilton import hamilton_equations as hamilton


def test_euler_shm():
    x = Function('x')
    p = Function('p')
    t = Symbol('t')
    h1 = x(t)**2 / 2 + p(t)**2 / 2
    assert hamilton(h1, x(t), p(t), t) == [
        Eq(x(t) + D(p(t), t), 0), Eq(p(t) - D(x(t), t), 0)
    ]
    y = Function('y')
    q = Function('q')
    h2 = x(t)**2 / 2 + y(t)**2 / 2 + p(t)**2 / 2 + q(t)**2 / 2
    assert hamilton(h2, [x(t), y(t)], [p(t), q(t)], t) == [
        Eq(x(t) + D(p(t), t), 0), Eq(p(t) - D(x(t), t), 0),
        Eq(y(t) + D(q(t), t), 0), Eq(q(t) - D(y(t), t), 0)
    ]
