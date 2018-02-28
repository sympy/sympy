"""
This module contains functions for the  of two multivariate resultants. These
are:

- Dixon's resultant.
- Maucalay's resultant.

Multivariate resultants are used to identify weather a multivariate system has
common roots. That is when the resultant is equal to zero.
"""

import sympy as sym
import functools

from sympy.polys.monomials import itermonomials
from sympy.polys.orderings import monomial_key
from sympy.functions.combinatorial.factorials import binomial

class DixonResultant():
    """
    A class for retrieving the Dixon's resultant of a multivariate system.

    Examples
    ========
    >>> from sympy.core import symbols
    >>> from sympy.utilities.lambdify import lambdify

    >>> from sympy.polys.multivariate_resultants import DixonResultant
    >>> x, y = symbols('x, y')

    >>> p = lambdify((x, y), x + y)
    >>> q = lambdify((x, y), x ** 2 + y **3)
    >>> h = lambdify((x, y), x ** 2 + y)

    >>> dixon = DixonResultant(variables=[x, y], polynomials=[p, q, h])
    >>> poly = dixon.get_dixon_polynomial()
    >>> matrix = dixon.get_dixon_matrix(polynomial=poly)
    >>> matrix
    Matrix([
    [ 0,  0, -1,  0, -1],
    [ 0, -1,  0, -1,  0],
    [-1,  0,  1,  0,  0],
    [ 0, -1,  0,  0,  1],
    [-1,  0,  0,  1,  0]])
    >>> matrix.det()
    0

    Reference
    ==========
    1. [Kapur1994]_
    2. [Palancz08]_

    See Also
    ========
    Notebook in examples:
    """

    def __init__(self, polynomials, variables):
        """
        A class that takes two lists, a list of polynomials and list of
        variables. Returns the Dixon matrix of the multivariate system.

        Parameters
        ----------
        variables: list
            A list of all n variables
        polynomials : list of sympy polynomials
            A  list of m n-degree polynomials
        """
        self.polynomials = polynomials
        self.variables = variables

        self.n = len(self.variables)
        self.m = len(self.polynomials)

        self.dummy_variables = self.get_dummy_variables()
        self.max_degrees = self.get_max_degrees()

    def get_dummy_variables(self):
        """
        Returns
        -------

        dummy_variables: list
            A list of n alpha variables. These are the replacing variables
        """
        a = sym.IndexedBase("alpha")
        dummy_variables = [a[i] for i in range(self.n)]

        return dummy_variables

    def get_max_degrees(self):
        r"""
        Returns
        -------

        degrees: list
            A list of the d_max of each variable. The max degree is the
            max(degree(p_1, x_i), ..., degree(p_m, x_i))
        """
        return [sym.Poly(f(*self.variables)).degree() for f in self.polynomials]

    def get_dixon_polynomial(self):
        r"""
        Returns
        -------

        dixon_polynomial: sympy polynomial
            A sympy polynomial which formulated by Dixon's formulation.
            Dixon's polynomial is calculated as:

            delta = Delta(A) / ((x_1 - a_1) ... (x_n - a_n)) where,

            A =  |p_1(x_1,... x_n), ..., p_n(x_1,... x_n)|
                 |p_1(a_1,... x_n), ..., p_n(a_1,... x_n)|
                 |...             , ...,              ...|
                 |p_1(a_1,... a_n), ..., p_n(a_1,... a_n)|
        """
        if self.m != (self.n + 1):
            raise ValueError('Method invalid for given combination.')

        # first row
        rows = [[poly(*self.variables) for poly in self.polynomials]]

        temp = [i for i in self.variables]
        iterator = iter(self.dummy_variables)

        for idx in (idx for idx in range(self.n)):
            temp[idx] = next(iterator)
            rows.append([poly(*temp) for poly in self.polynomials])

        A = sym.Matrix(rows)
        product_of_differences = functools.reduce(lambda x, y: x * y,
                                                  [a - b for a, b in zip(self.variables, self.dummy_variables)])
        dixon_polynomial = (A.det() / product_of_differences).factor()
        return sym.Poly(dixon_polynomial, *self.dummy_variables)

    def get_coefficients_of_alpha(self, polynomial):
        r"""
        Returns
        --------
        coefficients: list
            A list of coefficients (in x_i, ..., x_n terms) of the power products
            a_1, ..., a_n in Dixon's polynomial
        """
        coefficients = []
        for powers in polynomial.monoms():
            monomial = functools.reduce(lambda i, j: i * j,
                                    [a ** b for a, b in zip(self.dummy_variables, powers)])
            coefficients.append(polynomial.coeff_monomial(monomial))

        return coefficients

    def get_upper_degree(self):
        list_of_products = [self.variables[i] ** ((i + 1) * self.max_degrees[i] -1)
                            for i in range(self.n)]
        product = sym.prod(list_of_products)
        product = sym.Poly(product).monoms()

        return sym.polys.monomials.monomial_deg(*product)

    def get_dixon_matrix(self, polynomial):
        r"""
        Construct the Dixon matrix from the coefficients of polynomial \alpha.
        Each coefficient is viewed as a polynomial of x_1,..., x_n.
        """
        coefficients = self.get_coefficients_of_alpha(polynomial)
        monomials = list(itermonomials(self.variables, self.get_upper_degree()))
        monomials = sorted(monomials, key=monomial_key('lex', self.variables))[::-1]

        dixon_matrix = sym.Matrix([[sym.Poly(c, *self.variables).coeff_monomial(m) for m in monomials]
                                   for c in coefficients])

        keep = [column for column in range(dixon_matrix.shape[-1])
                if all([elements == 0 for elements in dixon_matrix[:, column]]) is False]
        return dixon_matrix[:, keep]

class MacaulayResultant():
    """
    A class for calculating the Macaulay resultant. Note that the coefficients
    of the polynomials must be given as symbols.

    Examples
    ========
    >>> from sympy.core import symbols
    >>> from sympy.utilities.lambdify import lambdify

    >>> from sympy.polys.multivariate_resultants import MacaulayResultant
    >>> x, y, z = symbols('x, y, z')

    >>> a_0, a_1, a_2 = symbols('a_0, a_1, a_2')
    >>> b_0, b_1, b_2 = symbols('b_0, b_1, b_2')
    >>> c_0, c_1, c_2,c_3, c_4 = symbols('c_0, c_1, c_2, c_3, c_4')

    >>> f = lambdify((x, y, z), a_0 * y -  a_1 * x + a_2 * z)
    >>> g = lambdify((x, y, z), b_1 * x ** 2 + b_0 * y ** 2 - b_2 * z ** 2)
    >>> h = lambdify((x, y, z), c_0 * y - c_1 * x ** 3 + c_2 * x ** 2 * z - c_3 * x * z ** 2 + c_4 * z ** 3)

    >>> mac = MacaulayResultant(polynomials=[f, g, h], variables=[x, y, z])
    >>> mac.get_monomials_set()
    >>> matrix = mac.get_matrix()
    >>> submatrix = mac.get_submatrix(matrix)
    >>> submatrix
    Matrix([
    [-a_1,  a_0,  a_2,    0],
    [   0, -a_1,    0,    0],
    [   0,    0, -a_1,    0],
    [   0,    0,    0, -a_1]])

    Reference
    ==========
    1. [Bruce97]_
    2. [Stiller96]_

    See Also
    ========
    Notebook in examples:
    """
    def __init__(self, polynomials, variables):
        """
        Parameters
        ----------
        variables: list
            A list of all n variables
        polynomials : list of sympy polynomials
            A  list of m n-degree polynomials
        """
        self.polynomials = polynomials
        self.variables = variables
        self.n = len(variables)

        self.degrees = self.get_max_degrees()
        self.degree_m = self.get_degree_m()
        self.monomials_size = self.get_size()

    def get_max_degrees(self):
        r"""
        Returns
        -------
        degrees: list
            A list of the d_max of each polynomial
        """
        degrees = [self.get_polynomial_degree(poly) for poly in self.polynomials]
        return degrees

    def get_polynomial_degree(self, poly):
        """
        Returns
        -------
        degree: int
            The degree of a polynomial
        """
        return sym.Poly(poly(*self.variables)).degree()

    def get_degree_m(self):
        r"""
        Returns
        -------
        degree_m: int
            The degree_m is calculated as  1 + \sum_1 ^ n (d_i - 1), where
            d_i is the degree of the i polynomial
        """
        return 1 + sum([d - 1 for d in self.degrees])

    def get_size(self):
        r"""
        Returns
        -------
        size: int
            The size of set T. Set T is the set of all possible monomials of
            the n variables for degree equal to the degree_m
        """
        return binomial(self.degree_m + self.n - 1, self.n - 1)

    def get_monomials_of_certain_degree(self, degree):
        """
        Returns
        -------
        monomials: list
            A list of monomials of a certain degree. Sympy returns up to a degree
        """
        monomials = list(itermonomials(self.variables, degree) -
                         itermonomials(self.variables, degree - 1))

        return sorted(monomials, key=monomial_key('lex', self.variables))[::-1]

    def get_monomials_set(self):
        r"""
        Returns
        -------
        self.monomial_set: set
            The set T. Set of all possible monomials of degree degree_m
        """
        monomial_set = self.get_monomials_of_certain_degree(self.degree_m)
        self.monomial_set = monomial_set

    def get_row_coefficients(self):
        """
        Returns
        -------
        row_coefficients: list
            The row coefficients of Macaulay's matrix
        """
        row_coefficients = []
        divisible = []
        for i in range(self.n):
            if i == 0:
                row_coefficients.append(self.get_monomials_of_certain_degree(self.degree_m - self.degrees[i]))
            else:
                divisible.append(self.variables[i - 1] ** self.degrees[i - 1])
                poss_rows = self.get_monomials_of_certain_degree(self.degree_m - self.degrees[i])
                for div in divisible:
                    for p in poss_rows:
                        if sym.fraction((p / div).expand())[1] == 1:
                            poss_rows = [item for item in poss_rows if item != p]
                row_coefficients.append(poss_rows)
        return row_coefficients

    def get_matrix(self):
        """
        Returns
        -------
        macaulay_matrix: sym Matrix
            The Macaulay's matrix
        """
        rows = []
        row_coefficients = self.get_row_coefficients()
        for i in range(self.n):
            for multiplier in row_coefficients[i]:
                coefficients = []
                poly = sym.Poly(self.polynomials[i](*self.variables) * multiplier, *self.variables)

                for mono in self.monomial_set:
                    coefficients.append(poly.coeff_monomial(mono))
                rows.append(coefficients)

        macaulay_matrix = sym.Matrix(rows)
        return macaulay_matrix

    def get_reduced_nonreduced(self):
        r"""
        Returns
        -------
        reduced: list
            A list of the reduced monomials
        non_reduced: list
            A list of the monomials that are not reduced

        Defition.
        ---------
        A polynomial is said to be reduced in x_i, if its degree (the maximum
        degree of its monomials) in x_i is less than d_i. A polynomial that is
        reduced in all variables but one is said simply to be reduced.
        """
        divisible = []
        for m in self.monomial_set:
            temp = []
            for i, v in enumerate(self.variables):
                temp.append(sym.Poly(m, v).degree() >= self.degrees[i])
            divisible.append(temp)

        reduced = [i for i, r in enumerate(divisible) if sum(r) < self.n - 1]
        non_reduced = [i for i, r in enumerate(divisible) if sum(r) >= self.n -1]

        return reduced, non_reduced

    def get_submatrix(self, matrix):
        r"""
        Returns
        -------
        macaulay_submatrix: sym Matrix
            The Macaulay's matrix. Columns that are non reduced are kept. The row
            which contain one if the a_{i}s is dropoed. a_{i}s are the coefficients
            of x_i ^ {d_i}.
        """
        reduced, non_reduced = self.get_reduced_nonreduced()

        reduction_set = [v ** self.degrees[i] for i, v in enumerate(self.variables)]

        ais = list([self.polynomials[i](*self.variables).coeff(reduction_set[i])
                    for i in range(self.n)])

        reduced_matrix = matrix[:, reduced]
        keep = []
        for row in range(reduced_matrix.rows):
            check = [ai in reduced_matrix[row, :] for ai in ais]
            if True not in check:
                keep.append(row)

        return matrix[keep, non_reduced]
