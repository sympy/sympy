""" Linear Solver for Holonomic Functions"""


from sympy.core import S
from sympy.matrices.common import ShapeError
from sympy.matrices.dense import MutableDenseMatrix


class NewMatrix(MutableDenseMatrix):
    """
    Supports elements which can't be Sympified.
    See docstrings in sympy/matrices/matrices.py
    """

    @staticmethod
    def _sympify(a):
        return a

    @classmethod
    def _new(cls, *args, **kwargs):
        # if the `copy` flag is set to False, the input
        # was rows, cols, [list].  It should be used directly
        # without creating a copy.
        if kwargs.get('copy', True) is False:
            if len(args) != 3:
                raise TypeError("'copy=False' requires a matrix be initialized as rows,cols,[list]")
            rows, cols, flat_list = args
        else:
            rows, cols, flat_list = cls._handle_creation_inputs(*args, **kwargs)
            flat_list = list(flat_list) # create a shallow copy
        self = object.__new__(cls)
        self.rows = rows
        self.cols = cols
        self._mat = flat_list
        return self

    def row_join(self, rhs):
        # Allows you to build a matrix even if it is null matrix
        if not self:
            return type(self)(rhs)

        if self.rows != rhs.rows:
            raise ShapeError(
                "`self` and `rhs` must have the same number of rows.")
        self_list = self.tolist()
        rhs_list = rhs.tolist()
        newmat_list = [srow + rrow for srow, rrow in zip(self_list, rhs_list)]
        return self._new(newmat_list)

    def col_join(self, bott):
        # Allows you to build a matrix even if it is null matrix
        if not self:
            return type(self)(bott)

        if self.cols != bott.cols:
            raise ShapeError(
                "`self` and `bott` must have the same number of columns.")
        newmat_list = self.tolist() + bott.tolist()
        return type(self)(newmat_list)

    def gauss_jordan_solve(self, b, freevar=False):
        from sympy.matrices import Matrix

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
        if not v[rank:, 0].is_zero_matrix:
            raise ValueError("Linear system has no solution")

        # Get index of free symbols (free parameters)
        free_var_index = permutation[len(pivots):]  # non-pivots columns are free variables

        # Free parameters
        tau = NewMatrix([S.One for k in range(col - rank)]).reshape(col - rank, 1)

        # Full parametric solution
        V = A[:rank, rank:]
        vt = v[:rank, 0]
        assert type(vt - V*tau) == NewMatrix
        free_sol = tau.vstack(vt - V*tau, tau)

        # Undo permutation
        sol = [0] * col
        for k, v in enumerate(free_sol):
            sol[permutation[k]] = v
        sol = NewMatrix(sol)

        if freevar:
            return sol, tau, free_var_index
        else:
            return sol, tau
