from sympy.matrices.dense import MutableDenseMatrix
from sympy.polys import Poly
from sympy.polys.polytools  import cancel
from sympy import symbols, Symbol
from sympy.logic.boolalg import false, true, Or
from sympy.series import limit
from sympy.core.relational import Unequality

__all__ = ['RouthHurwitz', 'negative_real_root_conditions']

def negative_real_root_conditions(polynomial: Poly):
    """
    Returns the conditions for a polynomial to have all its roots with
    negative real parts.

    Examples
    ========

    >>> from sympy.physics.control import negative_real_root_conditions
    >>> from sympy import symbols, Poly
    >>> s = symbols('s')
    >>> p1 = Poly((s+1)*(s+2)*(s+2.5), s)
    >>> negative_real_root_conditions(p1)
    [[True], [False]]
    >>> k = symbols('k')
    >>> p2 = Poly((s+1)*(s+k), s)
    >>> negative_real_root_conditions(p2)
    [[k + 1 > 0, k > 0], [False]]
    >>> a = symbols('a', negative = True)
    >>> p3 = Poly((s+1)*(s+a), s)
    >>> negative_real_root_conditions(p3)
    [[a + 1 > 0, False], [False]]

    """
    table = RouthHurwitz(polynomial)

    if table.zero_row_case:
        # There is a pole with a real part equal to zero or
        # a pole with a positive real part.
        return [[False], [False]]

    num_conditions = [true, true]
    var_conditions: list[list] = [[], []]

    _add_zero_col_conds(var_conditions, table)

    for eq in table[:, 0]:
        limit_p = limit(eq, table.infinitesimal_element, 0)
        if limit_p == 0:
            continue

        limit_n = limit(eq, table.infinitesimal_element, 0, dir = "-")

        if limit_p.is_number:
            num_conditions[0] &= limit_p > 0
            num_conditions[1] &= limit_n < 0

        else:
            var_conditions[0].append(limit_p > 0)
            var_conditions[1].append(limit_n < 0)

        if Or(num_conditions[0], num_conditions[1]) is false:
            return [[False], [False]]

    # if we are here, num_conditions != [[False], [False]]
    if num_conditions[0] is false:
        var_conditions[0] = [False]
        if len(var_conditions[1]) == 0:
            var_conditions[1] = [num_conditions[1]]

    elif num_conditions[1] is false:
        var_conditions[1] = [False]
        if len(var_conditions[0]) == 0:
            var_conditions[0] = [num_conditions[0]]

    return var_conditions

def _add_zero_col_conds(var_conditions, table):
    """
    Add checks for rows where the first-column entry is zero,
    ensuring the row isn't entirely zeros and preventing it from
    being classified as a Full Row Zero Case.

    """
    for i in table.zero_col_rows:
        conds_i = Or(*[Unequality(table[i, j], 0)
                       for j in range(1, table.shape[1])])
        var_conditions[0].append(conds_i)
        var_conditions[1].append(conds_i)

class RouthHurwitz(MutableDenseMatrix):
    r"""
    A class for creating a Routh-Hurwitz table from a given polynomial.
    It handle special cases with methods discussed in [1].

    Note: When a row of the table is zero, the property ``zero_row_case``
    is set to True.

    Explanation
    ============

    In mathematics, the Routh-Hurwitz table is used to determine the number of
    roots of a polynomial that have positive or negative real parts.

    It's crucial in the control system theory because it can be used to
    retrieve necessary and sufficient conditions for the stability of a linear
    time-invariant control system.

    Once the table is constructed, the stability of the system can be assessed
    by counting the number of sign changes in the first column.
    Each sign change corresponds to a root with a positive real part, whereas
    each preservation of sign corresponds to a root with a negative real part.

    There are two special cases to consider during the construction of the
    table:

    1.  First Column Zero Case:
        If a zero appears in the first column of a row (while the row is not
        entirely zero), the zero is replaced with a small symbolic infinitesimal
        (:math:`\epsilon`, by default).
        When checking for sign changes, the limit :math:`\epsilon \to 0`
        is taken.
        Specifically, if the other elements in the column are positive, then
        :math:`\epsilon \to 0^+`; if they are negative, then
        :math:`\epsilon \to 0^-`; otherwise, :math:`\epsilon \to 0^\pm`.

    2.  Full Row Zero Case:
        If an entire row becomes zero, we can substitute the row with the
        coefficients of the derivative of an auxiliary polynomial.
        The auxiliary polynomial is constructed using the row immediately
        above the zero row.
        For instance, consider the following example:

        .. math::
            \begin{matrix}3\\2\\1\end{matrix}\begin{bmatrix}b_3&b_1\\
            b_2&b_0\\ 0&0\end{bmatrix}

        The auxiliary polynomial will be: :math:`a(s) = b_2 s^2 + b_0`
        The characteristic is that, if :math:`p(s)` is the polynomial we are
        analyzing, :math:`p(s)=a(s)\cdot other(s)` and
        the roots of :math:`a(s)` are symmetric about the origin in the
        s-plane, so when we
        fall in this case, we should note that there could be poles with only
        imaginary part or poles with negative and positive real parts.


    The table is constructed for a polynomial of the form
    :math:`p(s) = b_n s^n + b_{n-1} s^{n-1} + \ldots + b_1 s + b_0`
    and the table has :math:`n+1` rows and the following structure:

    .. math::
        \begin{bmatrix}b_n&b_{n-2}&b_{n-4}&\cdots\\
        b_{n-1}&b_{n-3}&b_{n-5}&\cdots\\ c_{1}&c_2&c_3&\cdots\\
            d_{1}&d_2&d_3&\cdots\\ \vdots&\vdots&\vdots&\ddots\end{bmatrix}

    In this table, the elements in the subsequent rows are computed using the
    formulas:
    :math:`c_i = \frac{b_{n-1}\cdot b_{n-2i}-b_n\cdot b_{n-(2i+1)}}{b_{n-1}}`

    :math:`d_i = \frac{c_1 \cdot b_{n-(2i+1)}-b_{n-1}\cdot c_{i+1}}{c_1}`

    Parameters
    ==========

    polynomial : :py:class:`~.Expr`, :py:class:`~.Number`
        The polynomial whose Routh-Hurwitz table is to be created.
    var : :py:class:`~.Symbol`
        The symbol representing the variable in the polynomial.
    infinitesimal_element : None, :py:class:`~.Symbol`, optional
        The symbol representing the infinitesimal element for the first column
        zero case.
        If not provided, a default symbol ``epsilon`` will be used.

    Examples
    ========

    >>> from sympy import symbols, Poly
    >>> from sympy.physics.control import RouthHurwitz
    >>> b1, b2, b3, b4 = symbols('b_1 b_2 b_3 b_4')
    >>> s = symbols("s")

    Here is a generic example of how to use the ``RouthHurwitz`` class:

    >>> p = Poly(b1*s**3 + b2*s**2 + b3*s + b4, s)
    >>> RouthHurwitz(p)
    Matrix([
    [                     b_1, b_3],
    [                     b_2, b_4],
    [(-b_1*b_4 + b_2*b_3)/b_2,   0],
    [                     b_4,   0]])
    >>> RouthHurwitz(p)[:, 0]
    Matrix([
    [                     b_1],
    [                     b_2],
    [(-b_1*b_4 + b_2*b_3)/b_2],
    [                     b_4]])

    Here you can see how the table appears in the first column zero case:

    >>> p1 = Poly(s**4 + s**3 + 3*s**2 + 3*s + 3, s)
    >>> RouthHurwitz(p1)
    Matrix([
    [                      1, 3, 3],
    [                      1, 3, 0],
    [                epsilon, 3, 0],
    [(3*epsilon - 3)/epsilon, 0, 0],
    [                      3, 0, 0]])
    >>> RouthHurwitz(p1).zero_row_case
    False

    Here you can see how the table appears in the full row zero case
    (poles with only imaginary part):

    >>> p2 = Poly(s**6 + 2*s**5 + 8*s**4 + 12*s**3 + 20*s**2 + 16*s + 16, s)
    >>> RouthHurwitz(p2)
    Matrix([
    [  1,  8, 20, 16],
    [  2, 12, 16,  0],
    [  2, 12, 16,  0],
    [  8, 24,  0,  0],
    [  6, 16,  0,  0],
    [8/3,  0,  0,  0],
    [ 16,  0,  0,  0]])
    >>> RouthHurwitz(p2).zero_row_case
    True
    >>> RouthHurwitz(p2).auxiliary_polynomial
    Poly(2*s**4 + 12*s**2 + 16, s, domain='ZZ')

    References
    ==========
    .. [1] https://en.wikipedia.org/wiki/Routh-Hurwitz_stability_criterion
    .. [2] https://www.circuitbread.com/tutorials/routh-hurwitz-criterion-part-2-3-3

    """
    def __new__(cls, polynomial, infinitesimal_element = None):
        n = polynomial.degree()

        return super().__new__(cls, n + 1, n//2 + 1, [0]*(n + 1)*(n//2 + 1))

    def __init__(self, polynomial, infinitesimal_element = None):
        [var, ] = polynomial.gens
        if not isinstance(var, Symbol):
            raise ValueError("var must be a Symbol")

        self._var = var
        self._polynomial = Poly(polynomial, var)
        self._poly_degree = self._polynomial.degree()
        self._coeffs = self._polynomial.all_coeffs()

        self._zero_row_case = False
        self._zero_col_rows = []
        self._aux_poly_degree = 0

        self._inf_element = infinitesimal_element
        if infinitesimal_element is None:
            self._inf_element = symbols("epsilon", dummy=True)

        if self._poly_degree < 1:
            self[0, 0] = self._coeffs[0]
            return

        self._build_table()

    def _build_table(self):
        """Build the Routh-Hurwitz table."""
        self._initialize()
        self._calculate()

    def _initialize(self):
        """"Initialize the table with the coefficients of the polynomial."""
        row, col = 0, 0
        for coeff in self._coeffs:
            self[row, col] = coeff

            row = (row+1) % 2
            col = col + 1 - row

        if self[1, 0] != 0:
            return

        self._handle_special_cases(1)

    def _calculate(self):
        """Calculate the table using the first 2 rows."""
        for i in range(2, self.rows):
            self._calculate_row(i)
            self._handle_special_cases(i)

    def _calculate_row(self, i):
        trailing_columns = i//2
        for j in range(self.cols - trailing_columns):
            num = (self[i-1, 0] * self[i-2, j+1]
                   - self[i-2, 0] * self[i-1, j+1])
            den = self[i-1, 0]

            self[i, j] = cancel(num / den)

    def _handle_special_cases(self, i):
        """Handle the first column zero case and the full row zero case."""
        if all(self[i, j] == 0 for j in range(self.cols)):
            self._zero_row_case = True
            self._aux_poly_degree = self._poly_degree - i + 1
            # calculate the row using the auxiliary polynomial coefficients
            # degrees
            for j in range(self.cols):
                aux_coeff_deg = self._aux_poly_degree - 2*j

                if aux_coeff_deg < 0:
                    continue

                self[i, j] = self[i - 1, j] * aux_coeff_deg

            return

        if self[i, 0] == 0:
            self[i, 0] = self._inf_element
            self._zero_col_rows.append(i)

    @property
    def zero_col_rows(self):
        """
        Return a list of rows where the First Column Zero Case occurs.

        """
        return self._zero_col_rows

    @property
    def zero_row_case(self):
        """
        Return True if during the building of the table the Full Row Zero Case
        (see the explanation section) has been encountered, else False.

        """
        return self._zero_row_case

    @property
    def infinitesimal_element(self):
        """Return the infinitesimal element used in the table."""
        return self._inf_element

    @property
    def auxiliary_polynomial(self):
        """
        If the zero_row_case is True, return the auxiliary polynomial
        associated with the Full Row Zero Case.
        Otherwise, return None.

        It is used to handle the Full Row Zero Case during the
        construction of the Routh-Hurwitz table.

        """
        if self.zero_row_case is False:
            return None

        aux_poly = 0
        aux_poly_row = self._poly_degree - self._aux_poly_degree

        for j, exp in enumerate(range(self._aux_poly_degree, -1, -2)):
            aux_poly += self[aux_poly_row, j] * self._var**exp

        return Poly(aux_poly, self._var)
