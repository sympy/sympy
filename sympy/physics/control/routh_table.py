from sympy.matrices.dense import MutableDenseMatrix
from sympy import Poly
from sympy import simplify
from sympy import symbols, Symbol

__all__ = ['RouthHurwitz']


class RouthHurwitz(MutableDenseMatrix):

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

        for j in range(self.cols):
            if self[1, j] != 0:
                self[1, 0] = self[1, j] * (-1)**j
                return

    def _calculate(self):
        for i in range(2, self.rows):
            self._calculate_row(i)
            self._handle_special_cases(i)

    def _calculate_row(self, i):
        for j in range(self.cols):
            num = (self[i-1, 0] * self._get_safe_element(i-2, j+1)
                   - self[i-2, 0] * self._get_safe_element(i-1, j+1))
            den = self[i-1, 0]

            self[i, j] = simplify(num / den)

    def _handle_special_cases(self, i):
        if all(self[i, j] == 0 for j in range(self.cols)):
            self._zero_row_case = True
            # calculate the row using the auxiliary polynomial coefficients degrees
            for j in range(self.cols):
                aux_coeff_deg = self._poly_degree-i+1 - 2*j

                if aux_coeff_deg < 0:
                    continue

                self[i, j] = self[i - 1, j] * aux_coeff_deg

            return

        if self[i, 0] == 0:
            self[i, 0] = self._inf_element

    def _get_safe_element(self, i, j):
        if j >= self.cols:
            return 0

        return self[i, j]

    @property
    def first_column(self):
        return self[:, 0]

    @property
    def zero_rows_case(self):
        return self._zero_row_case

    @property
    def infinitesimal_element(self):
        return self._inf_element