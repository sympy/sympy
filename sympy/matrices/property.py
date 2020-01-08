from sympy.core.basic import Atom
from sympy.core.logic import fuzzy_and
from sympy.core.symbol import Symbol

from .common import MatrixRequired


class MatrixProperties(MatrixRequired):
    """Provides basic properties of a matrix.

    Notes
    =====

    Every matrix properties named like ``is_foo``, ``is_bar`` should be
    taking no optional keyword arguments and should be decorated with
    ``@property`` or any equivalent one (because they are incompatible).
    """

    def _eval_atoms(self, *types):
        result = set()
        for i in self:
            result.update(i.atoms(*types))
        return result

    def _eval_free_symbols(self):
        return set().union(*(i.free_symbols for i in self))

    def _eval_has(self, *patterns):
        return any(a.has(*patterns) for a in self)

    def _eval_is_anti_symmetric(self):
        def pred():
            for i in range(self.rows):
                for j in range(self.cols):
                    if i == j:
                        yield self[i, j].is_zero
                    else:
                        yield (self[i, j] + self[j, i]).is_zero
        return fuzzy_and(pred())

    def _eval_is_diagonal(self):
        def pred():
            for i in range(self.rows):
                for j in range(self.cols):
                    if i != j:
                        yield self[i, j].is_zero
        return fuzzy_and(pred())

    def _eval_is_hermitian_matrix(self):
        def pred():
            for i in range(self.rows):
                for j in range(self.cols):
                    if i == j:
                        yield self[i, j].is_real
                    else:
                        yield (self[i, j] - self[j, i].conjugate()).is_zero
        return fuzzy_and(pred())

    def _eval_is_Identity(self):
        def pred():
            for i in range(self.rows):
                for j in range(self.cols):
                    if i == j:
                        yield (self[i, j] - 1).is_zero
                    else:
                        yield self[i, j].is_zero
        return fuzzy_and(pred())

    def _eval_is_lower_hessenberg(self):
        def pred():
            for i in range(self.rows):
                for j in range(i+2, self.cols):
                    yield self[i, j].is_zero
        return fuzzy_and(pred())

    def _eval_is_lower(self):
        def pred():
            for i in range(self.rows):
                for j in range(i+1, self.cols):
                    yield self[i, j].is_zero
        return fuzzy_and(pred())

    def _eval_is_symbolic(self):
        return self.has(Symbol)

    def _eval_is_symmetric(self):
        def pred():
            for i in range(self.rows):
                for j in range(self.cols):
                    if i != j:
                        yield (self[i, j] - self[j, i]).is_zero
        return fuzzy_and(pred())

    def _eval_is_zero_matrix(self):
        def pred():
            for i in range(self.rows):
                for j in range(self.cols):
                    yield self[i, j].is_zero
        return fuzzy_and(pred())

    def _eval_is_upper_hessenberg(self):
        def pred():
            for i in range(2, self.rows):
                for j in range(min(self.cols, (i - 1))):
                    yield self[i, j].is_zero
        return fuzzy_and(pred())

    def _eval_is_upper(self):
        def pred():
            for i in range(1, self.rows):
                for j in range(min(self.cols, i)):
                    yield self[i, j].is_zero
        return fuzzy_and(pred())

    def _eval_values(self):
        return [i for i in self if not i.is_zero]

    def atoms(self, *types):
        """Returns the atoms that form the current object.

        Examples
        ========

        >>> from sympy.abc import x, y
        >>> from sympy.matrices import Matrix
        >>> Matrix([[x]])
        Matrix([[x]])
        >>> _.atoms()
        {x}
        """

        types = tuple(t if isinstance(t, type) else type(t) for t in types)
        if not types:
            types = (Atom,)
        return self._eval_atoms(*types)

    @property
    def free_symbols(self):
        """Returns the free symbols within the matrix.

        Examples
        ========

        >>> from sympy.abc import x
        >>> from sympy.matrices import Matrix
        >>> Matrix([[x], [1]]).free_symbols
        {x}
        """
        return self._eval_free_symbols()

    def has(self, *patterns):
        """Test whether any subexpression matches any of the patterns.

        Examples
        ========

        >>> from sympy import Matrix, SparseMatrix, Float
        >>> from sympy.abc import x, y
        >>> A = Matrix(((1, x), (0.2, 3)))
        >>> B = SparseMatrix(((1, x), (0.2, 3)))
        >>> A.has(x)
        True
        >>> A.has(y)
        False
        >>> A.has(Float)
        True
        >>> B.has(x)
        True
        >>> B.has(y)
        False
        >>> B.has(Float)
        True
        """
        return self._eval_has(*patterns)

    @property
    def is_anti_symmetric(self):
        """Check if the matrix is an anti-symmetric matrix

        Examples
        ========

        An example of an anti-symmetric matrix

        >>> from sympy import Matrix, symbols
        >>> m = Matrix([[0, 1], [-1, 0]])
        >>> m.is_anti_symmetric
        True

        An example of a logically undecidable matrix:

        >>> x, y = symbols('x y')
        >>> m = Matrix([[0, x], [y, 0]])
        >>> m.is_anti_symmetric

        An example of a matrix that can be logically decided, but
        should be simplified.

        >>> m = Matrix([
        ...     [0, x**2 + 2*x + 1, y],
        ...     [-(x + 1)**2 , 0, x*y],
        ...     [-y, -x*y, 0]])
        >>> m.is_anti_symmetric
        >>> m.expand().is_anti_symmetric
        True

        Notes
        =====

        A matrix $M$ is a anti-symmetric matrix if all
        $M_{i, j} = -M_{j, i}$. It may also imply that the all diagonal
        entries are zero.

        The return value can be a fuzzy value ``None`` if some pairs are
        logically undetermined. You can try out using ``simplify`` or
        ``expand`` to the matrix if you face such issues,
        or otherwise, try to build up the full logic programmatically
        and simplify some logically undecidable pairs.
        """
        if not self.is_square:
            return False

        return self._eval_is_anti_symmetric()

    @property
    def is_diagonal(self):
        """Check if the matrix is diagonal.

        Examples
        ========

        An example of a square diagonal matrix:

        >>> from sympy import Matrix
        >>> m = Matrix([[1, 0], [0, 2]])
        >>> m.is_diagonal
        True

        An example of a rectangular diagonal matrix:

        >>> m = Matrix([[1, 0, 0], [0, 2, 0]])
        >>> m.is_diagonal
        True

        An example of a non-diagonal matrix:

        >>> m = Matrix([[1, 1], [0, 2]])
        >>> m.is_diagonal
        False

        Notes
        =====

        This does not check whether a matrix is square.

        See Also
        ========

        is_lower
        is_upper
        sympy.matrices.matrices.MatrixEigen.is_diagonalizable
        diagonalize
        """
        return self._eval_is_diagonal()

    @property
    def is_hermitian(self):
        """Checks if the matrix is Hermitian.

        Examples
        ========

        >>> from sympy.matrices import Matrix
        >>> from sympy import I
        >>> from sympy.abc import x

        An example of a Hermitian matrix:

        >>> a = Matrix([[1, I], [-I, 1]])
        >>> a.is_hermitian
        True

        An example of a non-Hermitian matrix:

        >>> a = Matrix([[2*I, I], [-I, 1]])
        >>> a.is_hermitian
        False

        An example of a logically undecidable matrix:

        >>> a = Matrix([[x, I], [-I, 1]])
        >>> a.is_hermitian

        An example of a matrix with a logically decidable non-hermitian
        pairs:

        >>> a = Matrix([[x, 1], [-I, 1]])
        >>> a.is_hermitian
        False

        Notes
        =====

        In a Hermitian matrix $M$, the element $M_{i, j}$ is the complex
        conjugate of element $M_{j, i}$. This also implies that all
        diagonal entries are real.
        """
        if not self.is_square:
            return False

        return self._eval_is_hermitian_matrix()

    @property
    def is_Identity(self):
        if not self.is_square:
            return False
        return self._eval_is_Identity()

    @property
    def is_lower_hessenberg(self):
        r"""Checks if the matrix is in the lower-Hessenberg form.

        Examples
        ========

        >>> from sympy.matrices import Matrix
        >>> a = Matrix([
        ...     [1, 2, 0, 0],
        ...     [5, 2, 3, 0],
        ...     [3, 4, 3, 7],
        ...     [5, 6, 1, 1]])
        >>> a.is_lower_hessenberg
        True

        Notes
        =====

        The lower hessenberg matrix has zero entries
        above the first superdiagonal.

        See Also
        ========

        is_upper_hessenberg
        is_lower
        """
        return self._eval_is_lower_hessenberg()

    @property
    def is_lower(self):
        """Check if matrix is a lower triangular matrix.

        Examples
        ========

        >>> from sympy import Matrix
        >>> m = Matrix([[1, 0], [0, 1]])
        >>> m.is_lower
        True

        >>> m = Matrix([[0, 0, 0], [2, 0, 0], [1, 4, 0], [6, 6, 5]])
        >>> m.is_lower
        True

        An example of a logically undecidable matrix:

        >>> from sympy.abc import x, y
        >>> m = Matrix([[x, y], [0, y]])
        >>> m.is_lower

        Notes
        =====

        ``True`` can be returned even if the matrix is not square.

        See Also
        ========

        is_upper
        is_diagonal
        is_lower_hessenberg
        """
        return self._eval_is_lower()

    @property
    def is_square(self):
        """Checks if a matrix is square.

        Examples
        ========

        >>> from sympy import Matrix
        >>> a = Matrix([[1, 2, 3], [4, 5, 6]])
        >>> b = Matrix([[1, 2, 3], [4, 5, 6], [7, 8, 9]])
        >>> c = Matrix([])
        >>> a.is_square
        False
        >>> b.is_square
        True
        >>> c.is_square
        True

        Notes
        =====

        A matrix is square if the number of rows equals the number of
        columns. The empty matrix is square by definition, since the
        number of rows and the number of columns are both zero.
        """
        return self.rows == self.cols

    @property
    def is_symbolic(self):
        """Checks if any matrix elements contain Symbols.

        Examples
        ========

        >>> from sympy.matrices import Matrix
        >>> from sympy.abc import x, y
        >>> M = Matrix([[x, y], [1, 0]])
        >>> M.is_symbolic
        True

        """
        return self._eval_is_symbolic()

    @property
    def is_symmetric(self):
        """Check if matrix is a symmetric matrix.

        Examples
        ========

        >>> from sympy import Matrix
        >>> m = Matrix([[0, 1], [1, 2]])
        >>> m.is_symmetric
        True

        >>> m = Matrix([[0, 1], [2, 0]])
        >>> m.is_symmetric
        False

        >>> m = Matrix([[0, 0, 0], [0, 0, 0]])
        >>> m.is_symmetric
        False

        >>> from sympy.abc import x, y
        >>> m = Matrix([
        ...     [1, x**2 + 2*x + 1, y],
        ...     [(x + 1)**2 , 2, 0],
        ...     [y, 0, 3]])
        >>> m.is_symmetric
        >>> m.expand().is_symmetric
        True

        Notes
        =====

        In a symmetric matrix $M$, the element $M_{i, j}$ is the element
        $M_{j, i}$.
        """
        if not self.is_square:
            return False

        return self._eval_is_symmetric()

    @property
    def is_upper_hessenberg(self):
        """Checks if the matrix is the upper-Hessenberg form.

        Examples
        ========

        >>> from sympy.matrices import Matrix
        >>> a = Matrix([
        ...     [1, 4, 2, 3],
        ...     [3, 4, 1, 7],
        ...     [0, 2, 3, 4],
        ...     [0, 0, 1, 3]])
        >>> a.is_upper_hessenberg
        True

        Notes
        =====

        The upper hessenberg matrix has zero entries below the first
        subdiagonal.

        See Also
        ========

        is_lower_hessenberg
        is_upper
        """
        return self._eval_is_upper_hessenberg()

    @property
    def is_upper(self):
        """Check if matrix is an upper triangular matrix.

        Examples
        ========

        >>> from sympy import Matrix
        >>> m = Matrix([[1, 0], [0, 1]])
        >>> m.is_upper
        True

        >>> m = Matrix([[5, 1, 9], [0, 4, 6], [0, 0, 5], [0, 0, 0]])
        >>> m.is_upper
        True

        >>> m = Matrix([[4, 2, 5], [6, 1, 1]])
        >>> m.is_upper
        False

        Notes
        =====

        ``True`` can be returned even if the matrix is not square.

        See Also
        ========

        is_lower
        is_diagonal
        is_upper_hessenberg
        """
        return self._eval_is_upper()

    @property
    def is_zero_matrix(self):
        """Checks if a matrix is a zero matrix.

        Examples
        ========

        >>> from sympy import Matrix, zeros
        >>> from sympy.abc import x
        >>> a = Matrix([[0, 0], [0, 0]])
        >>> b = zeros(3, 4)
        >>> c = Matrix([[0, 1], [0, 0]])
        >>> d = Matrix([])
        >>> e = Matrix([[x, 0], [0, 0]])
        >>> a.is_zero_matrix
        True
        >>> b.is_zero_matrix
        True
        >>> c.is_zero_matrix
        False
        >>> d.is_zero_matrix
        True
        >>> e.is_zero_matrix

        Notes
        =====

        A matrix is zero if every element is zero.  A matrix need not be
        square to be considered zero.
        The empty matrix is zero by the principle of vacuous truth.
        For a matrix that may or may not be zero
        (e.g. contains a symbol), this will be ``None``.
        """
        return self._eval_is_zero_matrix()

    def values(self):
        """Return non-zero values of self."""
        return self._eval_values()