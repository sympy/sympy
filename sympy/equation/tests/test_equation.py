from sympy import Eqn, Equation, symbols

x, y = symbols('x y')


def test_Eqn():
    assert Eqn(x, y) == Equation(x, y)
    assert Eqn(x, y).lhs == x
    assert Eqn(x, y).rhs == y
