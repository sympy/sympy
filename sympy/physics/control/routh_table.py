from sympy.matrices.dense import MutableDenseMatrix
from sympy import Poly
from sympy import simplify
from sympy import symbols, Symbol

__all__ = ['RouthHurwitz']


class RouthHurwitz(MutableDenseMatrix):
    r"""
        A class for creating a Routh-Hurwitz table from a given polynomial.
        It handle special cases with methods discussed in [1].

        Note: When a row of the table is zero, the property ``zero_rows_case`` is set to True.

        Explanation
        ============

        TODO: Explain the special cases and why the ``zero_rows_case`` property is important.

        Parameters
        ==========

        polynomial : :py:class:`~.Expr`, :py:class:`~.Number`
            The polynomial whose Routh-Hurwitz table is to be created.
        var : :py:class:`~.Symbol`
            The symbol representing the variable in the polynomial.

        Examples
        ========

        >>> from sympy import symbols
        >>> from sympy.physics.control import RouthHurwitz
        >>> b1, b2, b3, b4 = symbols('b_1 b_2 b_3 b_4')
        >>> s = symbols("s")

        Here is a generic example of how to use the ``RouthHurwitz`` class:

        >>> p = b1*s**3 + b2*s**2 + b3*s + b4
        >>> RouthHurwitz(p, s)
        Matrix([
        [               b_1, b_3],
        [               b_2, b_4],
        [-b_1*b_4/b_2 + b_3,   0],
        [               b_4,   0]])
        >>> RouthHurwitz(p, s).first_column
        Matrix([
        [               b_1],
        [               b_2],
        [-b_1*b_4/b_2 + b_3],
        [               b_4]])

        Here you can see how the table appears in the first column zero case:

        >>> p1 = s**4 + s**3 + 3*s**2 + 3*s + 3
        >>> RouthHurwitz(p1, s)
        Matrix([
        [            1, 3, 3],
        [            1, 3, 0],
        [      epsilon, 3, 0],
        [3 - 3/epsilon, 0, 0],
        [            3, 0, 0]])
        >>> RouthHurwitz(p1, s).zero_rows_case
        False

        Here you can see how the table appears in the full row zero case (poles with only imaginary part):

        >>> p2 = s**6 + 2*s**5 + 8*s**4 + 12*s**3 + 20*s**2 + 16*s + 16
        >>> RouthHurwitz(p2, s)
        Matrix([
        [  1,  8, 20, 16],
        [  2, 12, 16,  0],
        [  2, 12, 16,  0],
        [  8, 24,  0,  0],
        [  6, 16,  0,  0],
        [8/3,  0,  0,  0],
        [ 16,  0,  0,  0]])
        >>> RouthHurwitz(p2, s).zero_rows_case
        True

        References
        ==========
        .. [1] https://en.wikipedia.org/wiki/Routh-Hurwitz_stability_criterion

    """
    def __new__(cls, polynomial, var):
        if not isinstance(var, Symbol):
            raise ValueError("var must be a Symbol")
        n = Poly(polynomial, var).degree()

        return super().__new__(cls, n + 1, n//2 + 1, [0]*(n + 1)*(n//2 + 1))

    def __init__(self, polynomial, var):
        self._polynomial = Poly(polynomial, var)
        self._poly_degree = self._polynomial.degree()
        self._coeffs = self._polynomial.all_coeffs()

        self._zero_row_case = False
        self._aux_poly_degree = 0

        self._inf_element = symbols("epsilon", dummy=True)

        self._build_table()

    def _build_table(self):
        self._initialize()
        self._calculate()

    def _initialize(self):
        row, col = 0, 0
        for coeff in self._coeffs:
            self[row, col] = coeff

            row = (row+1) % 2
            col = col + 1 - row

        if self[1, 0] != 0:
            return

        self._handle_special_cases(1)

    def _calculate(self):
        for i in range(2, self.rows):
            self._calculate_row(i)
            self._handle_special_cases(i)

    def _calculate_row(self, i):
        trailing_columns = i//2
        for j in range(self.cols - trailing_columns):
            num = (self[i-1, 0] * self[i-2, j+1]
                   - self[i-2, 0] * self[i-1, j+1])
            den = self[i-1, 0]

            self[i, j] = simplify(num / den)

    def _handle_special_cases(self, i):
        if all(self[i, j] == 0 for j in range(self.cols)):
            self._zero_row_case = True
            self._aux_poly_degree = self._poly_degree - i + 1
            # calculate the row using the auxiliary polynomial coefficients degrees
            for j in range(self.cols):
                aux_coeff_deg = self._aux_poly_degree - 2*j

                if aux_coeff_deg < 0:
                    continue

                self[i, j] = self[i - 1, j] * aux_coeff_deg

            return

        if self[i, 0] == 0:
            self[i, 0] = self._inf_element

    @property
    def first_column(self):
        """Return the first column of the table."""
        return self[:, 0]

    @property
    def zero_rows_case(self):
        """Return True if a row of the table was zero, else False."""
        return self._zero_row_case

    @property
    def infinitesimal_element(self):
        """Return the infinitesimal element used in the table."""
        return self._inf_element

    @property
    def auxiliary_polynomial(self):
        """If zero_rows_case is True, return the auxiliary polynomial, else None """
        # TODO
        pass
