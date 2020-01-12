from __future__ import division, print_function

from sympy.core.compatibility import reduce
from sympy.core.sympify import sympify

from .utilities import _iszero


def _columnspace(M, simplify=False, dotprodsimp=None):
    """Returns a list of vectors (Matrix objects) that span columnspace of ``M``

    Parameters
    ==========

    dotprodsimp : bool, optional
        Specifies whether intermediate term algebraic simplification is used
        during matrix multiplications to control expression blowup and thus
        speed up calculation.

    Examples
    ========

    >>> from sympy.matrices import Matrix, columnspace
    >>> M = Matrix(3, 3, [1, 3, 0, -2, -6, 0, 3, 9, 6])
    >>> M
    Matrix([
    [ 1,  3, 0],
    [-2, -6, 0],
    [ 3,  9, 6]])
    >>> columnspace(M)
    [Matrix([
    [ 1],
    [-2],
    [ 3]]), Matrix([
    [0],
    [0],
    [6]])]
    >>> columnspace(M) == M.columnspace()
    True

    See Also
    ========

    nullspace
    rowspace
    """

    reduced, pivots = M.echelon_form(simplify=simplify, with_pivots=True,
            dotprodsimp=dotprodsimp)

    return [M.col(i) for i in pivots]


def _nullspace(M, simplify=False, iszerofunc=_iszero, dotprodsimp=None):
    """Returns list of vectors (Matrix objects) that span nullspace of ``M``

    Parameters
    ==========

    dotprodsimp : bool, optional
        Specifies whether intermediate term algebraic simplification is used
        during matrix multiplications to control expression blowup and thus
        speed up calculation.

    Examples
    ========

    >>> from sympy.matrices import Matrix, nullspace
    >>> M = Matrix(3, 3, [1, 3, 0, -2, -6, 0, 3, 9, 6])
    >>> M
    Matrix([
    [ 1,  3, 0],
    [-2, -6, 0],
    [ 3,  9, 6]])
    >>> nullspace(M)
    [Matrix([
    [-3],
    [ 1],
    [ 0]])]
    >>> nullspace(M) == M.nullspace()
    True

    See Also
    ========

    columnspace
    rowspace
    """

    reduced, pivots = M.rref(iszerofunc=iszerofunc, simplify=simplify,
            dotprodsimp=dotprodsimp)

    free_vars = [i for i in range(M.cols) if i not in pivots]
    basis     = []

    for free_var in free_vars:
        # for each free variable, we will set it to 1 and all others
        # to 0.  Then, we will use back substitution to solve the system
        vec           = [M.zero] * M.cols
        vec[free_var] = M.one

        for piv_row, piv_col in enumerate(pivots):
            vec[piv_col] -= reduced[piv_row, free_var]

        basis.append(vec)

    return [M._new(M.cols, 1, b) for b in basis]


def _rowspace(M, simplify=False, dotprodsimp=None):
    """Returns a list of vectors that span the row space of ``M``.

    Parameters
    ==========

    dotprodsimp : bool, optional
        Specifies whether intermediate term algebraic simplification is used
        during matrix multiplications to control expression blowup and thus
        speed up calculation.

    Examples
    ========

    >>> from sympy.matrices import Matrix, rowspace
    >>> M = Matrix(3, 3, [1, 3, 0, -2, -6, 0, 3, 9, 6])
    >>> M
    Matrix([
    [ 1,  3, 0],
    [-2, -6, 0],
    [ 3,  9, 6]])
    >>> rowspace(M)
    [Matrix([[1, 3, 0]]), Matrix([[0, 0, 6]])]
    >>> rowspace(M) == M.rowspace()
    True
    """

    reduced, pivots = M.echelon_form(simplify=simplify, with_pivots=True,
            dotprodsimp=dotprodsimp)

    return [reduced.row(i) for i in range(len(pivots))]


def _orthogonalize(cls, *vecs, **kwargs):
    """Apply the Gram-Schmidt orthogonalization procedure
    to vectors supplied in ``vecs``.

    Parameters
    ==========

    vecs
        vectors to be made orthogonal

    normalize : bool
        If ``True``, return an orthonormal basis.

    rankcheck : bool
        If ``True``, the computation does not stop when encountering
        linearly dependent vectors.

        If ``False``, it will raise ``ValueError`` when any zero
        or linearly dependent vectors are found.

    Returns
    =======

    list
        List of orthogonal (or orthonormal) basis vectors.

    Examples
    ========

    >>> from sympy import I
    >>> from sympy.matrices import Matrix, orthogonalize
    >>> v = [Matrix([1, I]), Matrix([1, -I])]
    >>> orthogonalize(Matrix, *v)
    [Matrix([
    [1],
    [I]]), Matrix([
    [ 1],
    [-I]])]
    >>> orthogonalize(Matrix, *v) == Matrix.orthogonalize(*v)
    True

    See Also
    ========

    MatrixBase.QRdecomposition

    References
    ==========

    .. [1] https://en.wikipedia.org/wiki/Gram%E2%80%93Schmidt_process
    """

    normalize = kwargs.get('normalize', False)
    rankcheck = kwargs.get('rankcheck', False)

    def project(a, b):
        return b * (a.dot(b, hermitian=True) / b.dot(b, hermitian=True))

    def perp_to_subspace(vec, basis):
        """projects vec onto the subspace given
        by the orthogonal basis ``basis``"""

        components = [project(vec, b) for b in basis]

        if len(basis) == 0:
            return vec

        return vec - reduce(lambda a, b: a + b, components)

    ret  = []
    vecs = list(vecs) # make sure we start with a non-zero vector

    while len(vecs) > 0 and vecs[0].is_zero:
        if rankcheck is False:
            del vecs[0]
        else:
            raise ValueError("GramSchmidt: vector set not linearly independent")

    for vec in vecs:
        perp = perp_to_subspace(vec, ret)

        if not perp.is_zero:
            ret.append(cls(perp))
        elif rankcheck is True:
            raise ValueError("GramSchmidt: vector set not linearly independent")

    if normalize:
        ret = [vec / vec.norm() for vec in ret]

    return ret


# The following are top level stand-alone interface functions which sympify the
# matrix where needed (so it becomes immutable and returns immutable), otherwise
# the implementations above do not sympify due to problems in other parts of the
# codebase using matrices of unsympifiable objects.

def columnspace(self, simplify=False, dotprodsimp=None):
    return [sympify(v) for v in _columnspace(self, simplify=simplify,
            dotprodsimp=dotprodsimp)]

def nullspace(self, simplify=False, iszerofunc=_iszero, dotprodsimp=None):
    return [sympify(v) for v in _nullspace(self, simplify=simplify,
            iszerofunc=iszerofunc, dotprodsimp=dotprodsimp)]

def rowspace(self, simplify=False, dotprodsimp=None):
    return [sympify(v) for v in _rowspace(self, simplify=simplify,
            dotprodsimp=dotprodsimp)]

def orthogonalize(cls, *vecs, **kwargs):
    return [sympify(v) for v in _orthogonalize(cls, *vecs, **kwargs)]


columnspace.__doc__   = _columnspace.__doc__
nullspace.__doc__     = _nullspace.__doc__
rowspace.__doc__      = _rowspace.__doc__
orthogonalize.__doc__ = _orthogonalize.__doc__
