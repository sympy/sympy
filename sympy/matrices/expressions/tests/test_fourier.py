from sympy import S, I, ask, Q, Abs, simplify, exp, sqrt, symbols
from sympy.matrices.expressions.fourier import DFT, IDFT, DFTMatrix
from sympy.matrices import det, Matrix, Identity
from sympy.abc import n, i, j
from sympy.utilities.pytest import warns_deprecated_sympy


def test_deprecated_dft():
    with warns_deprecated_sympy():
        assert DFT(4).shape == (4, 4)
        assert ask(Q.unitary(DFT(4)))
        assert Abs(simplify(det(Matrix(DFT(4))))) == 1
        assert DFT(n)*IDFT(n) == Identity(n)
        assert DFT(n)[i, j] == exp(-2*S.Pi*I/n)**(i*j) / sqrt(n)

def test_dft_matrix():
    n, a, b = symbols('n a b')
    assert ask(Q.unitary(DFTMatrix(4, a=0, b=0))) == False
    assert ask(Q.unitary(DFTMatrix(4, a=0, b=1)))
    assert ask(Q.unitary(DFTMatrix(4, a=0, b=2))) == False
    assert ask(Q.unitary(DFTMatrix(4, a=1, b=0))) == False

    assert ask(Q.unitary(DFTMatrix(n, a=1, b=0))) == False
    assert ask(Q.unitary(DFTMatrix(4, a=a, b=0))) == None
    assert ask(Q.unitary(DFTMatrix(4, a=0, b=b))) == None
    assert ask(Q.unitary(DFTMatrix(n, a=a, b=0))) == None
    assert ask(Q.unitary(DFTMatrix(n, a=1, b=b))) == False
    assert ask(Q.unitary(DFTMatrix(n, a=a, b=b))) == None

    n = symbols('n', integer=True)
    assert ask(Q.unitary(DFTMatrix(n, a=0, b=1))) == True
