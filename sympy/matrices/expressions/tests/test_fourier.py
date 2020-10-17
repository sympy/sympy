from sympy import S, I, ask, Q, Abs, simplify, exp, sqrt
from sympy.core.symbol import symbols
from sympy.matrices.expressions.fourier import DFT, IDFT
from sympy.matrices import det, Matrix, Identity, Determinant
from sympy.testing.pytest import raises


def test_dft_creation():
    assert DFT(2)
    assert DFT(0)
    raises(ValueError, lambda: DFT(-1))
    raises(ValueError, lambda: DFT(2.0))
    raises(ValueError, lambda: DFT(2 + 1j))

    n = symbols('n')
    assert DFT(n)
    n = symbols('n', integer=False)
    raises(ValueError, lambda: DFT(n))
    n = symbols('n', negative=True)
    raises(ValueError, lambda: DFT(n))


def test_dft():
    n, i, j = symbols('n i j')
    assert DFT(4).shape == (4, 4)
    assert ask(Q.unitary(DFT(4)))
    assert Abs(simplify(det(Matrix(DFT(4))))) == 1
    assert DFT(n)*IDFT(n) == Identity(n)
    assert DFT(n)[i, j] == exp(-2*S.Pi*I/n)**(i*j) / sqrt(n)


def test_dft_determinant():
    assert det(DFT(4)) == I
    assert det(DFT(5)) == -1
    assert det(DFT(6)) == 1
    assert det(DFT(7)) == -I
    assert det(DFT(4)) == DFT(4).as_explicit().det()
    assert det(DFT(27)) == I

    assert det(IDFT(4)) == -I
    assert det(IDFT(4)) == IDFT(4).as_explicit().det()

    n = symbols('n')
    assert det(DFT(n)) == Determinant(DFT(n))


def test_dft_eval_eigenvals():
    assert DFT(3)._eval_eigenvals() == {-1: 1, 1: 1, -I: 1}
    assert DFT(44)._eval_eigenvals() == {-1: 11, 1: 12, -I: 11, I: 10}
    assert DFT(37)._eval_eigenvals() == {-1: 9, 1: 10, -I: 9, I: 9}
    assert DFT(26)._eval_eigenvals() == {-1: 7, 1: 7, -I: 6, I: 6}

    DFT(4)._eval_eigenvals() == DFT(4).as_explicit().eigenvals()

    assert IDFT(3)._eval_eigenvals() == {-1: 1, 1: 1, I: 1}
    assert IDFT(44)._eval_eigenvals() == {-1: 11, 1: 12, -I: 10, I: 11}
    assert IDFT(37)._eval_eigenvals() == {-1: 9, 1: 10, -I: 9, I: 9}
    assert IDFT(26)._eval_eigenvals() == {-1: 7, 1: 7, -I: 6, I: 6}

    IDFT(4)._eval_eigenvals() == IDFT(4).as_explicit().eigenvals()

    n = symbols('n')
    assert DFT(n)._eval_eigenvals() is None
