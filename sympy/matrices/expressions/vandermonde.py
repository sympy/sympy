from sympy.core.compatibility import iterable
from sympy.core.sympify import sympify
from sympy.matrices.expressions import MatrixExpr, Determinant
from sympy import S, Pow, Range

class VandermondeMatrix(MatrixExpr):
    r""" Generates a Vandermonde matrix.

        A Vandermonde matrix is a matrix with the terms of a geometric
        progression in each row. It typically arises in a polynomial
        least squares fitting, where each row corresponds to a polynomial point
        and each column corresponds to a power of that point.

        A square Vandermonde matrix for the points `x_1` to `x_n` is
         .. math ::
             \begin{bmatrix}
             1 & x_1 & x_1^2 & \hdots % x_1^{n-1} \\
             1 & x_2 & x_2^2 & \hdots % x_2^{n-1} \\
            \vdots & \vdots & \vdots & \ddots & \vdots \\
             1 & x_n & x_n^2 & \hdots % x_n^{n-1} \\
             \end{bmatrix}

        Parameters
        ==========

        points : iterable
            The points to generate the VandermondeMatrix for.

        n : integer, iterable or None
            Powers to use. For an integer, the powers range from 0 to n-1.
            For an iterable, use those powers. If ``None`` (default), the
            powers range from 0 to ``len(points)``-1.

        Examples
        ========

        >>> from sympy.matrices.expressions.vandermonde import VandermondeMatrix
        >>> v = VandermondeMatrix([0, 1, -1])
        >>> v.as_explicit()
        Matrix([
        [1,  0, 0],
        [1,  1, 1],
        [1, -1, 1]])

        >>> from sympy import symbols
        >>> x = symbols('x(1:4)')
        >>> VandermondeMatrix(x, range(4)).as_explicit()
        Matrix([
        [1, x1, x1**2, x1**3],
        [1, x2, x2**2, x2**3],
        [1, x3, x3**2, x3**3]])

        Sometimes, a Vandermonde matrix is defined without the power of zero
        column:

        >>> VandermondeMatrix(x, range(1,4)).as_explicit()
        Matrix([
        [x1, x1**2, x1**3],
        [x2, x2**2, x2**3],
        [x3, x3**2, x3**3]])

        As it is possible to specify an iterable for the powers, generalized
        versions can be rapildy created:

        >>> VandermondeMatrix([-1, 2, 3], [1, 3, 6]).as_explicit()
        Matrix([
        [-1, -1,   1],
        [ 2,  8,  64],
        [ 3, 27, 729]])

        It is also possible to use a symbolic column count:
        >>> n = symbols('n')
        >>> VandermondeMatrix([-1, 2, 3], n)
        VandermondeMatrix((-1, 2, 3), n)

        or symbolic powers

        >>> VandermondeMatrix(x, [0, 1, n]).as_explicit()
        Matrix([
        [1, x1, x1**n],
        [1, x2, x2**n],
        [1, x3, x3**n]])

        References
        ==========

        .. [1] https://en.wikipedia.org/wiki/Vandermonde_matrix

    """
    def __new__(cls, points, n=None):
        if not iterable(points):
            raise ValueError("First argument must be an iterable")
        points = sympify(tuple(points))
        n = sympify(n)
        if n:
            if iterable(n):
                n = tuple(n)
                cls.colpows = n
            else:
                cls._check_dim(n)
                cls.colpows = Range(n)
        else:
            cls.colpows = Range(len(points))

        if n is None:
            obj = super().__new__(cls, points)
        else:
            obj = super().__new__(cls, points, n)
        return obj

    @property
    def shape(self):
        rows = len(self.args[0])
        if len(self.args) == 1 or self.args[1] is None:
            return (rows, rows)
        elif iterable(self.args[1]):
            return (rows, len(self.args[1]))
        else:
            return (rows, self.args[1])

    def _entry(self, i, j, **kwargs):
        return Pow(self.args[0][i], self.colpows[j])

    def _eval_determinant(self):
        if self.is_square and self.colpows == Range(self.shape[1]):
            tmp = S.One
            for i in range(self.shape[1]-1):
                for j in range(i+1, self.shape[1]):
                    tmp *= (self.args[0][j] - self.args[0][i])
            return tmp
        return Determinant(self)
