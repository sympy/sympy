
from sympy.core.symbol import symbols
from sympy.physics.control.routh_table import RouthHurwitz
from sympy.matrices.dense import Matrix
from sympy.polys import Poly

s = symbols('s')

def test_table():
    b0, b1, b2, b3, b4 = symbols('b_0 b_1 b_2 b_3 b_4')
    p1 = b4 * s**4 + b3 * s**3 + b2 * s**2 + b1 * s + b0

    # generic polynomial tests
    t1 = RouthHurwitz(p1, s)
    expected1 = Matrix([
        [b4, b2, b0], [b3, b1, 0], [-b1*b4/b3 + b2, b0, 0],
        [(b0*b3**2 + b1*(b1*b4 - b2*b3))/(b1*b4 - b2*b3), 0, 0],
        [b0, 0, 0]])

    assert t1.equals(expected1)

    expected1_1 = Matrix([
        [b4], [b3], [-b1*b4/b3 + b2],
        [(b0*b3**2 + b1*(b1*b4 - b2*b3))/(b1*b4 - b2*b3)],
        [b0]])

    assert t1[:, 0].equals(expected1_1)
    assert t1.zero_row_case is False
    assert t1.auxiliary_polynomials is None

    # zero in the first column test case
    p2 = s**4 + s**3 + 3*s**2 + 3*s + 3
    t2 = RouthHurwitz(p2, s)

    expected2 = Matrix([
        [1, 3, 3], [1, 3, 0],[-3, 3, 0],
        [4, 0, 0], [3, 0, 0]])

    assert t2.equals(expected2)
    assert t2.zero_row_case is False

    # zero row test case
    p3 = s**6 + 2*s**5 + 8*s**4 + 12*s**3 + 20*s**2 + 16*s + 16
    t3 = RouthHurwitz(p3, s)

    expected3 = Matrix([
        [1, 8, 20, 16], [2, 12, 16, 0], [2, 12, 16, 0],
        [8, 24, 0, 0], [6, 16, 0, 0], [8/3, 0, 0, 0], [16, 0, 0, 0]])

    assert t3.equals(expected3)
    assert t3.zero_row_case is True
    assert t3.auxiliary_polynomials == [Poly(2*s**4 + 12*s**2 + 16, s)]
