from sympy import Abs, I, Q, S, ask, exp, simplify, sqrt
from sympy.abc import i, j, n
from sympy.matrices import Identity, Matrix, det
from sympy.matrices.expressions.fourier import DFT, IDFT


def test_dft():
    assert DFT(4).shape == (4, 4)
    assert ask(Q.unitary(DFT(4)))
    assert Abs(simplify(det(Matrix(DFT(4))))) == 1
    assert DFT(n)*IDFT(n) == Identity(n)
    assert DFT(n)[i, j] == exp(-2*S.Pi*I/n)**(i*j) / sqrt(n)
