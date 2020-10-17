from sympy.core import I
from sympy.core.symbol import symbols
from sympy.matrices.expressions.fourier import DFT
from sympy.matrices.expressions.vandermonde import VandermondeMatrix
from sympy.matrices import det, Matrix
from sympy.testing.pytest import raises


def test_vandermonde_creation():
    assert VandermondeMatrix([1, 2])
    assert VandermondeMatrix([1, 2], 4)
    assert VandermondeMatrix((1, 2), 4)
    assert VandermondeMatrix(range(4))
    assert VandermondeMatrix(range(4), 4)
    assert VandermondeMatrix(range(4), range(4))
    assert VandermondeMatrix(range(4), range(1, 5))
    raises(ValueError, lambda: VandermondeMatrix(1, 1))
    raises(ValueError, lambda: VandermondeMatrix([1, 2, 3], 2.0))
    raises(ValueError, lambda: VandermondeMatrix([1, 2, 3], -1))

    k, l, m, n = symbols('k l m n')
    assert VandermondeMatrix([k, l, m])
    assert VandermondeMatrix([k, l, m], n)
    assert VandermondeMatrix([k, l, m], 4)
    assert VandermondeMatrix([k, l, m], range(1, 2, 9))

    assert VandermondeMatrix([1, 2, 3], 4).shape == (3, 4)
    assert VandermondeMatrix([k, l, m], n).shape == (3, n)

    n = symbols('n', integer=False)
    raises(ValueError, lambda: VandermondeMatrix([k, l, m], n))
    n = symbols('n', negative=True)
    raises(ValueError, lambda: VandermondeMatrix([k, l, m], n))


def test_vandermonde():
    x = symbols('x(1:4)')
    x1, x2, x3 = x
    assert VandermondeMatrix(x).as_explicit() == Matrix([[1, x1, x1**2],
                                                         [1, x2, x2**2],
                                                         [1, x3, x3**2]])
    assert VandermondeMatrix(x, range(3)).equals(VandermondeMatrix(x))
    assert VandermondeMatrix(x, [0, 1, 2]).equals(VandermondeMatrix(x))
    assert VandermondeMatrix(x, 2).as_explicit() == Matrix([[1, x1],
                                                            [1, x2],
                                                            [1, x3]])
    assert VandermondeMatrix(range(3)).as_explicit() == Matrix([[1, 0, 0],
                                                                [1, 1, 1],
                                                                [1, 2, 4]])
    assert VandermondeMatrix(range(3), range(4)).as_explicit() == Matrix([
                                                                [1, 0, 0, 0],
                                                                [1, 1, 1, 1],
                                                                [1, 2, 4, 8]])


def test_vandermonde_dft():
    assert VandermondeMatrix([1, -I, -1, I], 4).equals(2*DFT(4))


def test_vandermonde_determinant():
    k, l, m = symbols('k l m')
    assert det(VandermondeMatrix([k, l, m], 3)).equals(VandermondeMatrix([k, l, m], 3).as_explicit().det())
    assert VandermondeMatrix([k, l, m], 3)._eval_determinant().equals(VandermondeMatrix([k, l, m], 3).as_explicit().det())
