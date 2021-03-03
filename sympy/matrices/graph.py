from sympy.utilities.iterables import \
    flatten, connected_components

from .common import NonSquareMatrixError


def _connected_components(M):
    """Returns the list of connected vertices of the graph when
    a square matrix is viewed as a weighted graph.

    Examples
    ========

    >>> from sympy import symbols, Matrix
    >>> a, b, c, d, e, f, g, h = symbols('a:h')
    >>> A = Matrix([
    ...     [a, 0, b, 0],
    ...     [0, e, 0, f],
    ...     [c, 0, d, 0],
    ...     [0, g, 0, h]])
    >>> A.connected_components()
    [[0, 2], [1, 3]]

    Notes
    =====

    Even if any symbolic elements of the matrix can be indeterminate
    to be zero mathematically, this only takes the account of the
    structural aspect of the matrix, so they will considered to be
    nonzero.
    """
    if not M.is_square:
        raise NonSquareMatrixError

    V = range(M.rows)
    E = sorted(M.todok().keys())
    return connected_components((V, E))


def _connected_components_decomposition(M):
    """Decomposes a square matrix into block diagonal form only
    using the permutations.

    Explanation
    ===========

    The decomposition is in a form of $A = P B P^{-1}$ where $P$ is a
    permutation matrix and $B$ is a block diagonal matrix.

    Returns
    =======

    P, B : PermutationMatrix, BlockDiagMatrix
        *P* is a permutation matrix for the similarity transform
        as in the explanation. And *B* is the block diagonal matrix of
        the result of the permutation.

        If you would like to get the diagonal blocks from the
        BlockDiagMatrix, see
        :meth:`~sympy.matrices.expressions.blockmatrix.BlockDiagMatrix.get_diag_blocks`.

    Examples
    ========

    >>> from sympy import symbols, Matrix
    >>> a, b, c, d, e, f, g, h = symbols('a:h')
    >>> A = Matrix([
    ...     [a, 0, b, 0],
    ...     [0, e, 0, f],
    ...     [c, 0, d, 0],
    ...     [0, g, 0, h]])

    >>> P, B = A.connected_components_decomposition()
    >>> P = P.as_explicit()
    >>> P_inv = P.inv().as_explicit()
    >>> B = B.as_explicit()

    >>> P
    Matrix([
    [1, 0, 0, 0],
    [0, 0, 1, 0],
    [0, 1, 0, 0],
    [0, 0, 0, 1]])
    >>> B
    Matrix([
    [a, b, 0, 0],
    [c, d, 0, 0],
    [0, 0, e, f],
    [0, 0, g, h]])
    >>> P * B * P_inv
    Matrix([
    [a, 0, b, 0],
    [0, e, 0, f],
    [c, 0, d, 0],
    [0, g, 0, h]])

    Notes
    =====

    This problem corresponds to the finding of the connected components
    of a graph, when a matrix is viewed as a weighted graph.
    """
    from sympy.combinatorics.permutations import Permutation
    from sympy.matrices.expressions.blockmatrix import BlockDiagMatrix
    from sympy.matrices.expressions.permutation import PermutationMatrix

    iblocks = M.connected_components()

    p = Permutation(flatten(iblocks))
    P = PermutationMatrix(p)

    blocks = []
    for b in iblocks:
        blocks.append(M[b, b])
    B = BlockDiagMatrix(*blocks)
    return P, B
