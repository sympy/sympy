from sympy import QQ
from sympy.abc import theta, x
from sympy.matrices import Matrix
from sympy.polys import cyclotomic_poly


def test_round_two():
    cases = (
        # A couple of cyclotomic fields:
        (cyclotomic_poly(5), Matrix.eye(4), 125),
        (cyclotomic_poly(7), Matrix.eye(6), -16807),
        # A couple of quadratic fields (one 1 mod 4, one 3 mod 4):
        (x ** 2 - 5, Matrix([[1, 1 / 2], [0, 1 / 2]]), 5),
        (x ** 2 - 7, Matrix([[1, 0], [0, 1]]), 28),
        # Dedekind's example of a field with 2 as essential disc divisor:
        (x ** 3 + x ** 2 - 2 * x + 8, Matrix([[1, 0, 0], [0, 1, 0], [0, 1/2, 1/2]]).T, -503),
        # A bunch of cubics with various forms for F -- all of these require
        # second or third enlargements. (Five of them require a third, while the rest require just a second.)
        # F = 2^2
        (x**3 + 3 * x**2 - 4 * x + 4, Matrix([(1 / 2, 1 / 4, 1 / 4), (0, 1 / 2, 1 / 2), (0, 0, 1)]).T, -83),
        # F = 2^2 * 3
        (x**3 + 3 * x**2 + 3 * x - 3, Matrix([(1 / 2, 0, 1 / 2), (0, 1, 0), (0, 0, 1)]).T, -108),
        # F = 2^3
        (x**3 + 5 * x**2 - x + 3, Matrix([(1 / 4, 0, 3 / 4), (0, 1 / 2, 1 / 2), (0, 0, 1)]).T, -31),
        # F = 2^2 * 5
        (x**3 + 5 * x**2 - 5 * x - 5, Matrix([(1 / 2, 0, 1 / 2), (0, 1, 0), (0, 0, 1)]).T, 1300),
        # F = 3^2
        (x**3 + 3 * x**2 + 5, Matrix([(1 / 3, 1 / 3, 1 / 3), (0, 1, 0), (0, 0, 1)]).T, -135),
        # F = 3^3
        (x**3 + 6 * x**2 + 3 * x - 1, Matrix([(1 / 3, 1 / 3, 1 / 3), (0, 1, 0), (0, 0, 1)]).T, 81),
        # F = 2^2 * 3^2
        (x**3 + 6 * x**2 + 4, Matrix([(1 / 3, 2 / 3, 1 / 3), (0, 1, 0), (0, 0, 1 / 2)]).T, -108),
        # F = 2^3 * 7
        (x**3 + 7 * x**2 + 7 * x - 7, Matrix([(1 / 4, 0, 3 / 4), (0, 1 / 2, 1 / 2), (0, 0, 1)]).T, 49),
        # F = 2^2 * 13
        (x**3 + 7 * x**2 - x + 5, Matrix([(1 / 2, 0, 1 / 2), (0, 1, 0), (0, 0, 1)]).T, -2028),
        # F = 2^4
        (x**3 + 7 * x**2 - 5 * x + 5, Matrix([(1 / 4, 0, 3 / 4), (0, 1 / 2, 1 / 2), (0, 0, 1)]).T, -140),
        # F = 5^2
        (x**3 + 4 * x**2 - 3 * x + 7, Matrix([(1 / 5, 4 / 5, 4 / 5), (0, 1, 0), (0, 0, 1)]).T, -175),
        # F = 7^2
        (x**3 + 8 * x**2 + 5 * x - 1, Matrix([(1 / 7, 6 / 7, 2 / 7), (0, 1, 0), (0, 0, 1)]).T, 49),
        # F = 2 * 5 * 7
        (x**3 + 8 * x**2 - 2 * x + 6, Matrix([(1, 0, 0), (0, 1, 0), (0, 0, 1)]).T, -14700),
        # F = 2^2 * 3 * 5
        (x**3 + 6 * x**2 - 3 * x + 8, Matrix([(1, 0, 0), (0, 1 / 4, 1 / 4), (0, 0, 1)]).T, -675),
        # F = 2 * 3^2 * 7
        (x**3 + 9 * x**2 + 6 * x - 8, Matrix([(1, 0, 0), (0, 1 / 2, 1 / 2), (0, 0, 1)]).T, 3969),
        # F = 2^2 * 3^2 * 7
        (x**3 + 15 * x**2 - 9 * x + 13, Matrix([(1 / 6, 1 / 3, 1 / 6), (0, 1, 0), (0, 0, 1)]).T, -5292),
    )
    for f, B_exp, d_exp in cases:
        K = QQ.algebraic_field((f, theta))
        B = K.integral_basis
        d = K.discriminant
        assert d == d_exp
        # The computed basis need not equal the expected one, but their quotient
        # must be unimodular:
        assert (B.inv()*B_exp).det()**2 == 1
