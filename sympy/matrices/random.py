from ..core.mul import Mul
from ..core.random import sample
from ..tensor import shape
from .dense import eye, diag

__all__ = 'regular_to_singular', 'elementary', 'triangular', 'square'


# === random number generator functions ===


def _sample(scalars, k=None):
    _k = k or 1
    if hasattr(scalars, 'sample'):
        smpl = scalars.sample(_k)
    else:
        smpl = sample(scalars, _k)
    return smpl if k else smpl[0]


# === default sets ===


_elementary_scalars = -1, 1
_elementary_units = _elementary_scalars


# === fundamental matrix functions ===


def regular_to_singular(mat, rank=None):
    r""" build matrix with linear combination of matrix columns to meet rank

    Explanation
    ===========
    Let $\mathbf{I}$ be the $n \times n$ identity matrix
    and $\mathbf{D}$ an diagonal matrix with only zero and one entries
    of rank $r$ (i.e. all entries but $n-r$ diagonal entries are zero).

    If a full rank matrix $\mathbf{A}$ is given, then

    .. math::

        \mathbf{B}
        = \mathbf{D} \cdot
        (\mathbf{I} + \mathbf{A} \cdot (\mathbf{I}-\mathbf{D}))
        = \mathbf{D} + \mathbf{D} \cdot \mathbf{A} \cdot(\mathbf{I}-\mathbf{D})

    is of rank $r$, too. Note, $\mathbf{A}\cdot\mathbf{B}$
    contains $r$ columns of $\mathbf{A}$.

    If $\mathbf{A}$ is an upper triangular matrix $\mathbf{B}$
    as well $\mathbf{A} \cdot \mathbf{B}$ is one, too.

    Here, $\mathbf{A} \cdot \mathbf{B}$ is returned.
    Note, if $n=r$ the matrix $\mathbf{B}$ is the identity.

    Examples
    ========

    .. ..testsetup::

       >>> from sympy.core.random import rng, seed
       >>> _rng_state = rng.getstate()
       >>> seed(1)

    >>> from sympy import Matrix
    >>> from sympy.matrices.random import regular_to_singular

    >>> m = Matrix([[1, 2, 3],[4, 5, 6], [7, 8 ,0]])
    >>> m.rank()
    3
    >>> n = regular_to_singular(m, 2)
    >>> n
    Matrix([
    [1, 26, 3],
    [4, 56, 6],
    [7, 14, 0]])
    >>> n.rank()
    2

    .. ..testcleanup::

       >>> assert not rng.getstate() == _rng_state
       >>> rng.setstate(_rng_state)
       >>> assert rng.getstate() == _rng_state

    Parameters
    ==========
    mat : Matrix
        a complex matrix
    rank : integer
        rank of matrix

    """
    if rank is None:
        return mat
    dim = shape(mat)[0]
    if rank == dim:
        return mat
    if dim < rank:
        raise ValueError(
            'rank of a matrix has to be less than or equal to dimension')

    i = eye(dim)
    # permutation
    p = eye(dim).permute(_sample(range(dim), dim))

    # d = projection(dim, index=(0, rank))
    spec = [1] * rank + [0] * (dim - rank)
    d = diag(*spec)

    d = p.inv() * d * p
    return mat * d * (i + mat * (i - d))


# === base matrices ===


def elementary(dim,
               *,
               index=None,
               scalar=None):
    r"""an elementary matrix n x n for Gauss elimination

    Explanantion
    ============
    Elementary matrices are matrices for Gauss elimination operation.

    In two dimensions any matrix of the following types are elemenarty

    .. math::

       \mathbf{T}
       = \left[\begin{array}{cc}0 & 1 \\ 1 & 0\end{array}\right] \\

       \mathbf{M}
       = \left[\begin{array}{cc}\lambda & 0 \\ 0 & 1\end{array}\right] \\

       \mathbf{A}
       = \left[\begin{array}{cc}1 & \mu \\ 0 & 1\end{array}\right] \\

    In higher dimensions, any matrix $\mathbf{A}$
    looking like identity matrix but with entries

    .. math::

        A_{ii} = a, \quad A_{ij} = b, \quad A_{ji} = -c, \quad A_{jj} = d

    and $ad - bc \neq 0$.
    So multiplication with $\mathbf{A}$ gives for

    * $a=0=d$ and $b=1=-c$ a transposition,
      i.e. swapping rows $i$ and $j$
    * $a=1$, $d=\lambda$ and $b=0=-c$ this matrix describes scaling the
      row $j$ by $\lambda$
    * $a=1=d$, $b=\mu$, $-c=0$ adding the $\mu$ multiple
      of the row $i$ to row $j$.

    Examples
    ========

    .. ..testsetup::

       >>> from sympy.core.random import rng, seed
       >>> _rng_state = rng.getstate()
       >>> seed(1)

    >>> from sympy.matrices.random import elementary

    >>> elementary(3)
    Matrix([
    [0, 0, 1],
    [0, 1, 0],
    [1, 0, 0]])

    >>> elementary(3, index=(0,2), scalar=5)  # no more random
    Matrix([
    [1, 0, 5],
    [0, 1, 0],
    [0, 0, 1]])

    >>> elementary(3, index=(2,2), scalar=5)  # no more random
    Matrix([
    [1, 0, 0],
    [0, 1, 0],
    [0, 0, 5]])

    >>> elementary(3, index=(2,1), scalar=4)  # no more random
    Matrix([
    [1, 0, 0],
    [0, 1, 0],
    [0, 4, 1]])

    >>> elementary(3, index=(0,1), scalar=None)  # no more random
    Matrix([
    [0, 1, 0],
    [1, 0, 0],
    [0, 0, 1]])

    .. ..testcleanup::

       >>> assert not rng.getstate() == _rng_state
       >>> rng.setstate(_rng_state)
       >>> assert rng.getstate() == _rng_state

    Parameters
    ==========
    dim : integer
        dimension of matrix
    index : tuple(integer, integer) (optional)
        coordinates ``(i,j)`` of matrix operation
        if ``i`` equals ``j`` the *elementary matrix*
        will be of type $\mathbf{M}$,
        else of type $\mathbf{A}$ or $\mathbf{T}$

        if not given, **index** is drawn randomly
    scalar : symbol (optional)
        value of elementary entry,
        defaults to values -1 or 1
        if **scalar** is **None** the *elementary matrix*
        will be a transposition $\mathbf{T}$ or identity

    See Also
    ========
    triangular
    square

    """

    row, col = index or _sample(range(dim), 2)

    obj = eye(dim)
    if row == col:
        # identity or multiply by scalar
        obj[row, col] = scalar or 1 # if scalar is not None else 1
    elif scalar is None:
        # transposition
        obj[row, row] = obj[col, col] = 0
        obj[row, col] = obj[col, row] = 1
    else:
        # add scalar multiple to another
        obj[row, col] = scalar

    return obj


# === compound matrices, i.e. product of base matrices ===


def triangular(dim,
               *,
               rank=None,
               scalars=_elementary_scalars,
               units=_elementary_units,
               length=None):
    r"""n x n upper triangular matrix with random entries.

    Explanation
    ===========
    Matrix with values placed above and on the diagonal. Constructed as
    product of random upper *triangular* :func:`elementary` matrices.

    It is constructed from an invertible triangular matrix $\mathbf{S}$
    of which $r$ basis columns are randomly choosen
    and the $n-r$ columns are build as linear combinations
    of the basis colums.

    Examples
    ========

    .. ..testsetup::

       >>> from sympy.core.random import rng, seed
       >>> _rng_state = rng.getstate()
       >>> seed(1)

    >>> from sympy.matrices.random import triangular

    >>> triangular(3, scalars=(2, -2, 0))
    Matrix([
    [-1, -2,  2],
    [ 0,  1, -2],
    [ 0,  0, -1]])

    >>> triangular(3, scalars=(2, -2, 0), units=(1, 2, 3))
    Matrix([
    [1, 4, 0],
    [0, 4, 0],
    [0, 0, 9]])

    >>> triangular(3, rank=2, scalars=(2, -2, 0), units=(1, 2, 3))
    Matrix([
    [0, 2, 0],
    [0, 2, 0],
    [0, 0, 9]])

    >>> triangular(3, length=99)
    Matrix([
    [1, 5, -3],
    [0, 1, -1],
    [0, 0,  1]])

    .. ..testcleanup::

       >>> assert not rng.getstate() == _rng_state
       >>> rng.setstate(_rng_state)
       >>> assert rng.getstate() == _rng_state

    Parameters
    ==========
    dim : integer
        dimension of matrix
    rank : integer (optional with default dim)
        rank of matrix
    scalars : tuple or list of symbols (optional)
        default values for random choosen scalars of
        elementary matrices to build the invertible matrix
    units : tuple or list of symbols (optional)
        default values for random choosen scalar of
        elementary matrices to build the invertible matrix
        and non-zero diagonal entries
        which are assumed to have a multiplicative inverse
    length : integer (optional with default 2 * **dim**)
        number of invertible matrices to build the resulting matrix.

    See Also
    ========
    elementary
    square

    """

    if length == 0:
        return regular_to_singular(eye(dim), rank)
    scalars = scalars or (0,)
    units = units or (1,)
    length = length or 2 * dim
    row_indicies = [_sample(range(dim)) for _ in range(length)]
    indicies = [(i, _sample(range(i, dim))) for i in row_indicies]
    scalars = [_sample(units) if i == j
               else _sample(scalars) for i, j in indicies]
    items = [elementary(dim, index=ix, scalar=s)
             for ix, s in zip(indicies, scalars)]
    return regular_to_singular(Mul(*items), rank)


def square(dim,
           *,
           rank=None,
           scalars=_elementary_scalars,
           units=_elementary_units,
           length=None):
    r"""a square matrix n x n with a given rank

    Explanation
    ===========
    A full rank square matrix has **rank=dim**.
    Such a matrix is constructed as
    product of :func:`elementary` matrices.

    If **rank** is less than **dim**,
    it is constructed from a full rank matrix $\mathbf{S}$
    by choosing randomly $r$ basis columns
    and build $n-r$ columes as a linear combinations these.

    Such a square matrix $\mathbf{S}$ represents an endomorphism
    $$f_S:V \rightarrow V, v \mapsto \mathbf{S}\cdot v$$
    of a $n$ dimensional vector space $V$.

    Parameters
    ==========
    dim : integer
        dimension of matrix
    rank : integer (optional with default **dim**)
        rank of matrix
    scalars : tuple or list of symbols (optional)
        default values for random choosen scalar of
        elementary matrices to build the matrix
    units : tuple or list of symbols (optional)
        default values for random choosen scalar of
        elementary matrices to build the matrix
        which are assumed to have a multiplicative inverse
    length : integer (optional with default 2 * **dim**)
        number of invertible matrices
        to build the resulting matrix.

    Examples
    ========

    .. ..testsetup::

       >>> from sympy.core.random import rng, seed
       >>> _rng_state = rng.getstate()
       >>> seed(1)

    >>> from sympy.matrices.random import square

    >>> square(3)
    Matrix([
    [-1, -1, -1],
    [-2, -1, -1],
    [ 0,  0,  1]])

    >>> m = square(4, rank=2)
    >>> m.rank()
    2
    >>> m
    Matrix([
    [ 1, 0,  0,  0],
    [ 0, 0,  1, -1],
    [ 0, 0, -1,  1],
    [-1, 0,  0,  0]])

    .. ..testcleanup::

       >>> assert not rng.getstate() == _rng_state
       >>> rng.setstate(_rng_state)
       >>> assert rng.getstate() == _rng_state

    See Also
    ========
    triangular

    """

    if length == 0:
        return regular_to_singular(eye(dim), rank)

    length = length or 2 * dim
    lwr = triangular(
        dim, rank=None, scalars=scalars, units=units, length=length // 2)
    upr = triangular(
        dim,
        rank=rank, scalars=scalars, units=units, length=length - length // 2)
    return lwr.T * upr
