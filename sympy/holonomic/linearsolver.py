""" Linear Solver for Holonomic Functions"""

from __future__ import print_function, division

from sympy.matrices.dense import MutableDenseMatrix
from sympy.utilities.iterables import flatten, numbered_symbols
from sympy.core.symbol import Symbol, Dummy, symbols
from sympy import S


class NewMatrix(MutableDenseMatrix):
    """
    Supports elements which can't be Sympified.
    """
    @staticmethod
    def _sympify(a):
        return a

    def row_join(self, rhs):
        """Concatenates two matrices along self's last and rhs's first column

        Examples
        ========

        >>> from sympy import zeros, ones
        >>> M = zeros(3)
        >>> V = ones(3, 1)
        >>> M.row_join(V)
        Matrix([
        [0, 0, 0, 1],
        [0, 0, 0, 1],
        [0, 0, 0, 1]])

        See Also
        ========

        row
        col_join
        """
        from sympy.matrices import MutableMatrix
        # Allows you to build a matrix even if it is null matrix
        if not self:
            return type(self)(rhs)

        if self.rows != rhs.rows:
            raise ShapeError(
                "`self` and `rhs` must have the same number of rows.")
        newmat = NewMatrix.zeros(self.rows, self.cols + rhs.cols)
        newmat[:, :self.cols] = self
        newmat[:, self.cols:] = rhs
        return type(self)(newmat)

    def col_join(self, bott):
        """Concatenates two matrices along self's last and bott's first row

        Examples
        ========

        >>> from sympy import zeros, ones
        >>> M = zeros(3)
        >>> V = ones(1, 3)
        >>> M.col_join(V)
        Matrix([
        [0, 0, 0],
        [0, 0, 0],
        [0, 0, 0],
        [1, 1, 1]])

        See Also
        ========

        col
        row_join
        """
        from sympy.matrices import MutableMatrix
        # Allows you to build a matrix even if it is null matrix
        if not self:
            return type(self)(bott)

        if self.cols != bott.cols:
            raise ShapeError(
                "`self` and `bott` must have the same number of columns.")
        newmat = NewMatrix.zeros(self.rows + bott.rows, self.cols)
        newmat[:self.rows, :] = self
        newmat[self.rows:, :] = bott
        return type(self)(newmat)

    def gauss_jordan_solve(self, b, freevar=False):
        """
        Solves Ax = b using Gauss Jordan elimination.

        There may be zero, one, or infinite solutions.  If one solution
        exists, it will be returned. If infinite solutions exist, it will
        be returned parametrically. If no solutions exist, It will throw
        ValueError.

        Parameters
        ==========

        b : Matrix
            The right hand side of the equation to be solved for.  Must have
            the same number of rows as matrix A.

        freevar : List
            If the system is underdetermined (e.g. A has more columns than
            rows), infinite solutions are possible, in terms of an arbitrary
            values of free variables. Then the index of the free variables
            in the solutions (column Matrix) will be returned by freevar, if
            the flag `freevar` is set to `True`.

        Returns
        =======

        x : Matrix
            The matrix that will satisfy Ax = B.  Will have as many rows as
            matrix A has columns, and as many columns as matrix B.

        params : Matrix
            If the system is underdetermined (e.g. A has more columns than
            rows), infinite solutions are possible, in terms of an arbitrary
            parameters. These arbitrary parameters are returned as params
            Matrix.

        Examples
        ========

        >>> from sympy import Matrix
        >>> A = Matrix([[1, 2, 1, 1], [1, 2, 2, -1], [2, 4, 0, 6]])
        >>> b = Matrix([7, 12, 4])
        >>> sol, params = A.gauss_jordan_solve(b)
        >>> sol
        Matrix([
        [-2*_tau0 - 3*_tau1 + 2],
        [                 _tau0],
        [           2*_tau1 + 5],
        [                 _tau1]])
        >>> params
        Matrix([
        [_tau0],
        [_tau1]])

        >>> A = Matrix([[1, 2, 3], [4, 5, 6], [7, 8, 10]])
        >>> b = Matrix([3, 6, 9])
        >>> sol, params = A.gauss_jordan_solve(b)
        >>> sol
        Matrix([
        [-1],
        [ 2],
        [ 0]])
        >>> params
        Matrix(0, 1, [])

        See Also
        ========

        lower_triangular_solve
        upper_triangular_solve
        cholesky_solve
        diagonal_solve
        LDLsolve
        LUsolve
        QRsolve
        pinv

        References
        ==========

        .. [1] http://en.wikipedia.org/wiki/Gaussian_elimination

        """
        from sympy.matrices import Matrix, zeros

        aug = self.hstack(self.copy(), b.copy())
        row, col = aug[:, :-1].shape

        # solve by reduced row echelon form
        A, pivots = aug.rref()
        A, v = A[:, :-1], A[:, -1]
        pivots = list(filter(lambda p: p < col, pivots))
        rank = len(pivots)

        # Bring to block form
        permutation = Matrix(range(col)).T
        A = A.vstack(A, permutation)

        for i, c in enumerate(pivots):
            A.col_swap(i, c)

        A, permutation = A[:-1, :], A[-1, :]

        # check for existence of solutions
        # rank of aug Matrix should be equal to rank of coefficient matrix
        if not v[rank:, 0].is_zero:
            raise ValueError("Linear system has no solution")

        # Get index of free symbols (free parameters)
        free_var_index = permutation[len(pivots):]  # non-pivots columns are free variables

        # Free parameters
        tau = NewMatrix([S(1) for k in range(col - rank)]).reshape(col - rank, 1)

        # Full parametric solution
        V = A[:rank, rank:]
        vt = v[:rank, 0]
        free_sol = tau.vstack(vt - V*tau, tau)

        # Undo permutation
        sol = NewMatrix.zeros(col, 1)
        for k, v in enumerate(free_sol):
            sol[permutation[k], 0] = v

        if freevar:
            return sol, tau, free_var_index
        else:
            return sol, tau
