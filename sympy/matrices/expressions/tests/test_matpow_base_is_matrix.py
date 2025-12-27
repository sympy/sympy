from sympy import MatrixSymbol, symbols, Derivative


def test_derivative_subs_matrix_power():
    j = symbols("j")
    A = MatrixSymbol("A", 2, 2)
    d = Derivative(A**j, A)
    d.subs(j, 3)
