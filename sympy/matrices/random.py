import random

from ..core import I as _i
from ..core.mul import Mul as _multiply
from ..functions import sqrt as _sqrt, re as _re, im as _im, \
    transpose as _t, adjoint as _c
from ..tensor import shape as _shape
from .dense import eye as _eye

__all__ = 'identity', 'projection', 'jordan', 'transposition', 'permutation', \
          'elementary', 'rotation', 'reflection', \
          'diagonal_normal', 'jordan_normal', 'isometry_normal', \
          'triangular', 'invertible', 'singular', \
          'idempotent', 'nilpotent', 'diagonalizable', 'trigonalizable', \
          'orthogonal', 'unitary', 'normal', \
          'symmetric', 'hermite', \
          'regular_to_singular', 'complex_to_real'

_TEST = False

_eps = 1e-13


def _get_value(o, i, j):
    return o[i, j]


def _set_value(o, i, j, v):
    o[i, j] = v


def _inv(o):
    return o.inv()


def _cs(scalar):
    sgn = random.choice((-1, 1))
    if isinstance(scalar, (tuple, list)) and len(scalar) == 2:
        return scalar
    if isinstance(scalar, complex) or getattr(scalar, 'is_real', True) is False:
        return _re(scalar), _im(scalar)
    c = scalar
    s = sgn * _sqrt(1 - c * c)
    return c, s

    # msg = "rotation scalar argument must have norm equal to 1"
    # msg += " or norm less than 1 and real"
    # msg += " not abs(%s)=%s" % (str(scalar), str(abs(scalar)))
    # raise ValueError(msg)


_elementary_scalar_set = -1, 1
_rotation_scalar_set = (_sqrt(2) / 2, _sqrt(2) / 2), (0, -1),
_unitary_scalar_set = tuple(c * 1 + s * _i for c, s in _rotation_scalar_set)


# === fundamental constructor ===

def super_elementary_matrix(dim,
                            index=None,
                            value=None,
                            *scalar_set):
    r"""super elementary matrix n x n, i.e. identity with on 2 x 2 block

    Explanation
    ===========
    The super elementary matrix $A$ is a gerealization of elementary matrices
    as well as of a rotation matrices that rotate only a single plane.

    In two dimensions any invertible matrix

    .. math::

       A = \left[\begin{array}{cc}a & b \\ -c & d\end{array}\right]

    is super elementary, i.e. $ad+bc \neq 0$. In higher dimensions,
    any matrix $A$ looking like identity matrix but with entries

    .. math::

        A[i,i] = a, \quad A[i,j] = b, \quad A[j,i] = -c, \quad A[j,j] = d

    and $ad - bc \neq 0$.

    This inculdes elementary matrices of matrix operation for Gauss elimination.
    So multiplication with $A$ gives for

    * $a=0=d$ and $b=1=-c$ a *transposition*, i.e. swapping rows $i$ and $j$
    * $a=1$, $d=\lambda$ and $b=0=-c$ this matrix describes scaling the
      row $j$ by $\lambda$
    * $a=1=d$, $b=\lambda$, $-c=0$ adding the $\lambda$ multiple
      of the row $i$ to row $j$.

    Moreover, a simple rotation by $\phi$ in $(i,j)$ plane is super elementary,
    given by $a=\cos \phi = d$ and $b = \sin \phi = c$. Hence,

    .. math::

       \left[\begin{array}{cc}a & b \\ -c & d\end{array}\right]
       =
       \left[\begin{array}{cc}
        \cos \phi & \sin \phi \\ -\sin \phi & \cos \phi \end{array}\right]

    In this module the super elementary matrix serves as the base class to
    create futher matrices of given type.

    Examples
    ========

    >>> from sympy.matrices.random import super_elementary_matrix
    >>> import random
    >>> random.seed(1)

    >>> super_elementary_matrix(3, (0,1))
    Matrix([
    [0, 1, 0],
    [1, 0, 0],
    [0, 0, 1]])

    >>> super_elementary_matrix(3, (0,2), value=5)
    Matrix([
    [1, 0, 5],
    [0, 1, 0],
    [0, 0, 1]])

    >>> super_elementary_matrix(3, (2,2), value=5)
    Matrix([
    [1, 0, 0],
    [0, 1, 0],
    [0, 0, 5]])

    >>> super_elementary_matrix(3, (1,2), 1,2,3,4)
    Matrix([
    [1, 0, 0],
    [0, 1, 2],
    [0, -3, 4]])

    Parameters
    ==========
    dim : integer
        dimension of matrix
    index : tuple(int, int)
        coordinates (i,j) of super elementary square
    value : symbol (optional)
        value as elementary entry
    *scalar_set : symbols (optional)
        up to three additional values as additional super elementry entries

    See Also
    ========
    identity
    diagonal
    elementary
    transposition
    rotation

    """
    obj = _eye(dim)
    if index is None:
        return obj

    if not isinstance(index, (list, tuple, set)) or not len(index) == 2:
        raise ValueError(
            "index argument must be tuple of two matrix index integer.")
    row, col = index
    if not isinstance(row, int) or not isinstance(col, int):
        raise ValueError(
            "index argument must be tuple of two matrix index integer.")

    if row == col:
        # identity or _multiply by scalar
        _set_value(obj, row, col, value or 1)
        if scalar_set:
            raise ValueError(
                "if row==col no further args than value may be supplied.")
    elif not scalar_set:
        # elementary
        if value is None:
            # transposition
            _set_value(obj, row, row, 0)
            _set_value(obj, row, col, 1)
            _set_value(obj, col, row, 1)
            _set_value(obj, col, col, 0)
        else:
            # add scalar multiple to another
            _set_value(obj, row, col, value)
    elif len(scalar_set) == 3:
        _set_value(obj, row, row, value)
        y, v, u = scalar_set
        _set_value(obj, row, col, y)
        _set_value(obj, col, row, -v)
        _set_value(obj, col, col, u)
    else:
        msg = "either no, one or three additional scalar_set arguments " \
              "may be supplied, but not %i"
        raise ValueError(msg % len(scalar_set))
    return obj


# === fundamental matrix functions ===


def complex_to_real(mat=None):
    r""" returns real version of a complex matrix

    Explanation
    ===========
    Any complex vector space $V$ with basis $(b_1, \dots, b_n)$
    can be seen as a real vector space $V^{prime}$ with basis
    $(b_1, b_1 \ i, \dots, b_n, b_n \ i)$
    of dimension $2 \ n$.

    Since multiplication with complex number is linear on $V$
    resp. $V^{\prime}$,
    any complex homomorphism (i.e. complex matrix $A$) on $V$
    becomes a real homomorphism (i.e. matrix $A^{\prime}$) on $V^{\prime}$.

    This transforms a complex matrix $A$
    to the corresponding real matrix $A^{\prime}$.

    Examples
    ========
    >>> from sympy import Matrix, I
    >>> from sympy.abc import a, b
    >>> from sympy.matrices.random import complex_to_real

    >>> z = a + b * I
    >>> A = Matrix([[z]])
    >>> A
    Matrix([[a + I*b]])

    >>> complex_to_real(A)
    Matrix([
    [ re(a) - im(b), re(b) + im(a)],
    [-re(b) - im(a), re(a) - im(b)]])

    >>> A = Matrix([[z, 0],[1, I]])
    >>> A
    Matrix([
    [a + I*b, 0],
    [      1, I]])
    >>> complex_to_real(A)
    Matrix([
    [ re(a) - im(b), re(b) + im(a),  0, 0],
    [-re(b) - im(a), re(a) - im(b),  0, 0],
    [             1,             0,  0, 1],
    [             0,             1, -1, 0]])

    Parameters
    ==========
    mat : Matrix
        a complex matrix

    """
    dim = _shape(mat)[0]
    mat = mat or identity(dim)
    obj = super_elementary_matrix(2 * dim)
    for i in range(dim):
        for j in range(dim):
            z_value = _get_value(mat, i, j)
            a, b = _re(z_value), _im(z_value)
            _set_value(obj, 2 * i, 2 * j, a)
            _set_value(obj, 2 * i, 2 * j + 1, b)
            _set_value(obj, 2 * i + 1, 2 * j, -b)
            _set_value(obj, 2 * i + 1, 2 * j + 1, a)
    return obj


def regular_to_singular(mat, rank=None):
    r""" build matrix with linear combination of matrix columns to meet rank

    Explanation
    ===========
    Let $I$ be the $n \times n$ identity matrix
    and $D$ an diagonal matrix with only zero and one entries
    of rank $r$ (i.e. all entries but $n-r$ diagonal entries are zero).

    If a full rank matrix $A$ is given, then

    .. math::
        B = D * (I + A * (I-D)) = D + D * A * (I-D)

    is of rank $r$, too.  Some for $A*B$ which contains $r$ columns with $A$.

    If $A$ is an upper triangular matrix $B$ as well $A*B$ is one, too.

    Here, $A*B$ is returned. Note, if $n=r$ the matrix $B$ is the identity.

    Examples
    ========
    >>> from sympy import Matrix
    >>> from sympy.matrices.random import regular_to_singular
    >>> m = Matrix([[1, 2, 3],[4, 5, 6], [7, 8 ,0]])
    >>> m.rank()
    3
    >>> n = regular_to_singular(m, 2)
    >>> n
    Matrix([
    [1, 2, 15],
    [4, 5, 42],
    [7, 8, 69]])
    >>> n.rank()
    2

    Parameters
    ==========
    mat : Matrix
        a complex matrix
    rank : integer
        rank of matrix

    """
    if rank is None:
        return mat
    dim = _shape(mat)[0]
    if rank == dim:
        return mat
    i = identity(dim)
    p = permutation(dim)
    d = projection(dim, (0, rank))
    d = _multiply(p.inv(), d, p)
    md = _multiply(mat, d)
    return md + _multiply(md, mat, i - d)


# === base matrices ===


def identity(dim):
    r"""identity matrix n x n

    Explanation
    ===========

    Creates identity matrix with only ones on the diagonal
    and zeros anywhere else.

    Examples
    ========

    >>> from sympy.matrices.random import identity

    >>> identity(3)
    Matrix([
    [1, 0, 0],
    [0, 1, 0],
    [0, 0, 1]])

    Parameters
    ==========
    dim : integer
        dimension n of matrix

    See Also
    ========
    diagonal
    sympy.matrices.dense.eye

    """

    return super_elementary_matrix(dim)


def projection(dim,
               index=None):
    r"""a projection matrix n x n

    Explanation
    ===========
    A projection is a identity like matrix
    but with zero diagonal entiries off **index**.

    Examples
    ========
    >>> from sympy.matrices.random import projection
    >>> import random
    >>> random.seed(1)

    >>> projection(3)
    Matrix([
    [0, 0, 0],
    [0, 1, 0],
    [0, 0, 0]])

    >>> projection(3, index=(1,3))
    Matrix([
    [0, 0, 0],
    [0, 1, 0],
    [0, 0, 1]])

    Parameters
    ==========
    dim : integer
        dimension n of matrix
    index : tuple(integer, integer) (optional)
        coordinates ``(i,j)`` for start (included) and end (excluded)
        to project onto

    See Also
    ========
    identity
    diagonal

    """

    index = index or random.sample(range(dim + 1), 2)
    obj = identity(dim)
    start, end = sorted(index)
    for i in range(dim):
        v = 1 if start <= i < end else 0
        _set_value(obj, i, i, v)
    return obj


def jordan(dim,
           index=None,
           scalar=None):
    r"""n x n matrix with a single Jordan block of given eigenvalue

    Explanation
    ===========
    A matrix with a single Jordan block matrix
    with only non-zero blocks on the diagonal
    which have the form of an Jordan block $J$ which is

    .. math::

        J = \left[\begin{array}{cccccc}
            \lambda  & 1       & 0         & \dots     &         & 0     \\
             0       & \lambda & 1         & 0         & \dots   & 0     \\
             0       & \ddots  & \ddots    & \ddots    & \ddots  & 0     \\
             0       & \dots   & 0         & 0         & \lambda & 1     \\
             0       & \dots   &           &           & 0       & \lambda
             \end{array}\right]

    Finally, each $\lambda$ of each Jordan block will be an eigenvalue.

    >>> from sympy.matrices.random import jordan
    >>> import random
    >>> random.seed(1)

    >>> jordan(3)
    Matrix([
    [1, 1, 0],
    [0, 1, 0],
    [0, 0, 1]])
    >>> jordan(3, scalar=2)
    Matrix([
    [2, 0, 0],
    [0, 1, 0],
    [0, 0, 1]])
    >>> jordan(4, index=(1,4), scalar=2)
    Matrix([
    [1, 0, 0, 0],
    [0, 2, 1, 0],
    [0, 0, 2, 1],
    [0, 0, 0, 2]])

    Parameters
    ==========
    dim : integer
        dimension of matrix
    index : tuple(integer, integer) (optional)
        coordinates ``(i,j)`` for start (included) and end (excluded)
        of the Jordan square
    scalar : symbol (optional)
        eigenvalue used for the Jordan block
        defaults to values -1 or 1

    See Also
    ========
    jordan_normal

    """
    index = index or random.sample(range(dim), 2)
    scalar = random.choice(_elementary_scalar_set) if scalar is None else scalar
    obj = identity(dim)
    start, end = sorted(index)
    _set_value(obj, start, start, scalar)
    for i in range(start + 1, end):
        _set_value(obj, i - 1, i, 1)
        _set_value(obj, i, i, scalar)
    return obj


def transposition(dim,
                  index=None):
    r"""n x n transposition matrix

    Explamation
    ===========
    A transposition matrix is an identity matrix where two columns (or rows)
    are swapped. It is a special permutation matrix.
    Moreover, any permutaion matrix is a product of transposition matrices.

    In two dimensions it is

    .. math::

       T = \left[\begin{array}{cc}0 & 1 \\ 1 & 0\end{array}\right] \\

    If index is given, it sets which column or row will be swapped.

    Examples
    ========

    >>> from sympy.matrices.random import transposition
    >>> import random
    >>> random.seed(1)

    >>> transposition(4)
    Matrix([
    [1, 0, 0, 0],
    [0, 0, 1, 0],
    [0, 1, 0, 0],
    [0, 0, 0, 1]])

    >>> transposition(3, (0,1))
    Matrix([
    [0, 1, 0],
    [1, 0, 0],
    [0, 0, 1]])

    Parameters
    ==========
    dim : integer
        dimension n of matrix
    index : tuple(integer, integer) (optional)
        coordinates (i,j) of transposition

    See Also
    ========
    identity
    elementary
    permutation

    """
    index = index or random.sample(range(dim), 2)
    return super_elementary_matrix(dim, index)


def permutation(dim,
                perm=None):
    r"""permutation matrix n x n

    Explamation
    ===========

    A permutation matrix is a matrix with only 0 or 1 entries
    and each column and each row consists of only one non-zero entry.

    Such a matrix can be obtained by shuffeling the rows
    of an identiy matrix .i.e ``permutation(dim, perm)`` is eqivalent to
    ``identity(dim).permute(perm)``.

    Examples
    ========
    >>> import random
    >>> random.seed(1)
    >>> from sympy.matrices.random import permutation

    >>> permutation(3)
    Matrix([
    [1, 0, 0],
    [0, 0, 1],
    [0, 1, 0]])

    >>> permutation(3, (2,0,1))
    Matrix([
    [0, 0, 1],
    [1, 0, 0],
    [0, 1, 0]])

    Parameters
    ==========
    dim : integer
        dimension n of matrix
    perm : tuple or list of integer (optional)
        permutation as list of n integers, e.g. ``[2,1,0,3]`` as a
        permutation of ``[0,1,2,3]``

    See Also
    ========
    identity
    transposition

    """
    perm = perm or random.sample(range(dim), dim)
    obj = identity(dim)
    for i, j in enumerate(perm):
        _set_value(obj, i, i, 0)  # aka zeros(dim)
        _set_value(obj, i, j, 1)
    return obj  # aka identity(dim).permute(perm)


def elementary(dim,
               index=None,
               scalar=None):
    r"""an elementary matrix n x n for Gauss elimination

    Explanantion
    ============
    Elementary matrices are matrices for Gauss elimination operation.

    In two dimensions any matrix of the following types are elemenarty

    .. math::

       T = \left[\begin{array}{cc}0 & 1 \\ 1 & 0\end{array}\right] \\

       M = \left[\begin{array}{cc}\lambda & 0 \\ 0 & 1\end{array}\right] \\

       A = \left[\begin{array}{cc}1 & \mu \\ 0 & 1\end{array}\right] \\

    In higher dimensions, any matrix $A$ looking like identity matrix
    but with entries

    .. math::

        A[i,i] = a, \quad A[i,j] = b, \quad A[j,i] = -c, \quad A[j,j] = d

    and $ad - bc \neq 0$.
    So multiplication with $A$ gives for

    * $a=0=d$ and $b=1=-c$ a ``transposition``, i.e. swapping rows $i$ and $j$
    * $a=1$, $d=\lambda$ and $b=0=-c$ this matrix describes scaling the
        row $j$ by $\lambda$
    * $a=1=d$, $b=\mu, $-c=0$ adding the $\mu multiple
        of the row $i$ to row $j$.

    Examples
    ========
    >>> from sympy.matrices.random import elementary
    >>> import random
    >>> random.seed(1)

    >>> elementary(3)
    Matrix([
    [0, 0, 1],
    [0, 1, 0],
    [1, 0, 0]])

    >>> elementary(3, (0,2), scalar=5)
    Matrix([
    [1, 0, 5],
    [0, 1, 0],
    [0, 0, 1]])

    >>> elementary(3, (2,2), scalar=5)
    Matrix([
    [1, 0, 0],
    [0, 1, 0],
    [0, 0, 5]])

    >>> elementary(3, (2,1), 4)
    Matrix([
    [1, 0, 0],
    [0, 1, 0],
    [0, 4, 1]])

    >>> elementary(3, (0,1), None)
    Matrix([
    [0, 1, 0],
    [1, 0, 0],
    [0, 0, 1]])

    Parameters
    ==========
    dim : integer
        dimension of matrix
    index : tuple(integer, integer) (optional)
        coordinates ``(i,j)`` of matrix operation
        if ``i`` equals ``j`` the *elementary matrix*
        will be of type $M$,
        else of type $A$ or $T$
    scalar : symbol (optional)
        value of elementary entry,
        defaults to values -1 or 1
        if **scalar** is **None** the *elementary matrix*
        will be a transposition $T$ or identity

    See Also
    ========
    diagonal
    transposition

    """

    index = index or random.sample(range(dim), 2)
    # scalar = scalar or random.choice(_elementary_scalar_set + (None,))
    return super_elementary_matrix(dim, index, scalar)


def rotation(dim,
             index=None,
             scalar=None):
    r"""a square matrix n x n of a plane rotation

    Explanation
    ===========
    The matrix decribes in n dimensional Euclidian space
    a rotation of $(i,j)$ plane given by a single rotation square

    .. math::

        \left[\begin{array}{cc} c & s \\ -s & c \end{array}\right]


    Examples
    ========
    >>> from sympy import sqrt, cos, symbols, eye
    >>> from sympy.matrices.random import rotation
    >>> import random
    >>> random.seed(1)

    >>> rotation(3)
    Matrix([
    [0, 0, -1],
    [0, 1,  0],
    [1, 0,  0]])

    >>> cos_a = sin_a = sqrt(2)/2
    >>> rotation(3, scalar=cos_a)
    Matrix([
    [1,          0,         0],
    [0,  sqrt(2)/2, sqrt(2)/2],
    [0, -sqrt(2)/2, sqrt(2)/2]])

    >>> r = rotation(3, scalar=(cos_a, sin_a))
    >>> r
    Matrix([
    [1,         0,          0],
    [0, sqrt(2)/2, -sqrt(2)/2],
    [0, sqrt(2)/2,  sqrt(2)/2]])

    >>> r.T * r == eye(3)
    True

    works with symbols too

    >>> cos_phi = cos(symbols('phi'))
    >>> random.seed(1)
    >>> r =rotation(3, scalar=cos_phi)
    >>> r
    Matrix([
    [              cos(phi), 0, sqrt(1 - cos(phi)**2)],
    [                     0, 1,                     0],
    [-sqrt(1 - cos(phi)**2), 0,              cos(phi)]])

    >>> r.det()
    1

    >>> r.T * r == eye(3)
    True

    Parameters
    ==========
    dim : integer
        dimension of matrix
    index : tuple(integer, integer) (optional)
        coordinates ``(i,j)`` of rotation plane
    scalar : tuple(symbol, symbol) or symbol (optional)
        either a tuple of cosine value $c$ and sine value $s$ of rotation square
        or just cosine value $c$ of rotation square.

        If given, the **scalar** $c$ is be between -1 and 1.
        The resulting rotation square

        .. math::

            \left[\begin{array}{cc} c & \pm s \\ \mp s & c \end{array}\right]

        takes **scalar** for $c$ and $\pm \sqrt{1-c^2}$ for $s$.

    See Also
    ========
    reflection
    isometry_normal
    orthogonal

    """
    index = index or random.sample(range(dim), 2)
    scalar = scalar or random.choice(_rotation_scalar_set)
    c, s = _cs(scalar)
    return super_elementary_matrix(dim, index, c, s, s, c)


def reflection(dim,
               index=None,
               scalar=None):
    r"""a square matrix n x n of a hyperplane reflection

    Explanation
    ===========
    The matrix describes a reflaction in n dimensional Euclidian space.
    Constructed as ``reflection =  rotation * transposition``, i.e.
    there will be the reflection square

    .. math::

        \left[\begin{array}{cc} s & c \\ c & -s \end{array}\right]

    Examples
    ========
    >>> from sympy import sqrt, eye
    >>> from sympy.matrices.random import reflection
    >>> import random
    >>> random.seed(1)

    >>> reflection(3)
    Matrix([
    [-1, 0, 0],
    [ 0, 1, 0],
    [ 0, 0, 1]])

    >>> c = s = sqrt(2)/2
    >>> reflection(3, scalar=c)
    Matrix([
    [1,         0,          0],
    [0, sqrt(2)/2,  sqrt(2)/2],
    [0, sqrt(2)/2, -sqrt(2)/2]])

    >>> r = reflection(3, scalar=(c, s))
    >>> r
    Matrix([
    [1,          0,         0],
    [0, -sqrt(2)/2, sqrt(2)/2],
    [0,  sqrt(2)/2, sqrt(2)/2]])

    >>> r.det()
    -1

    >>> r.T * r == eye(3)
    True

    Parameters
    ==========
    dim : integer
        dimension of matrix
    index : tuple(integer, integer) (optional)
        coordinates ``(i,j)`` of transposition and rotation
    scalar : tuple(symbol, symbol) or symbol (optional)
        either a tuple of cosine value $c$ and sine value $s$ of rotation square
        or just cosine value $c$ of rotation square.

        If given, the **scalar** $c$ as be between -1 and 1.
        The resulting rotation square

        .. math::

            \left[\begin{array}{cc} c & \pm s \\ \mp s & c \end{array}\right]

        takes **scalar** for $c$ and $\pm \sqrt{1-c^2}$ for $s$.

    See Also
    ========
    rotation

    """

    index = index or random.sample(range(dim), 2)
    return rotation(dim, index, scalar) * transposition(dim, index)


# === normal form matrices, i.e. defined by eigenvalues ===


def diagonal_normal(dim,
                    spec=None):
    r"""n x n matrix with random values placed on the diagonal.

    .. _matrices-random-diagonal:

    Explanation
    ===========
    Creates a square diagonal matrix, i.e. matrix with only zero non-diagonal
    entries.
    The  diagonal values will be choose from **spec**,
    the *spectrum* argument which is the set of eigenvalues.

    Examples
    ========

    >>> from sympy.matrices.random import diagonal_normal
    >>> import random
    >>> random.seed(1)

    >>> diagonal_normal(3)
    Matrix([
    [-1,  0, 0],
    [ 0, -1, 0],
    [ 0,  0, 1]])

    >>> diagonal_normal(3, (4,-3,2))
    Matrix([
    [4,  0, 0],
    [0, -3, 0],
    [0,  0, 2]])

    >>> diagonal_normal(3, (-2,2))
    Matrix([
    [-2, 0, 0],
    [ 0, 2, 0],
    [ 0, 0, 2]])

    Parameters
    ==========
    dim : integer
        dimension n of matrix
    spec : tuple or list of symbols (optional)
        set of values of which scalars (diagonal entries) are choosen.

        If **dim** meets the length of **spec**,
        spec will be the eigenvalues (diagonal entries) as it is.

        If **dim** and the length of **spec** differs,
        the diagonal entries will be choosen randomly from **spec**.

        If not given **spec** defaults to $\{ -1, 1 \}$.

    See Also
    ========
    identity
    sympy.matrices.dense.diag

    """
    spec = spec or _elementary_scalar_set

    # choose spec randomly if dim and len(spec) does not meet
    if not dim == len(spec):
        spec = tuple(random.choice(spec) for _ in range(dim))

    # set diagonal entries
    if _TEST:
        # multiplicative matrix construction (only for testing)
        items = list()
        for start, scalar in enumerate(spec):
            items.append(elementary(dim, index=(start, start), scalar=scalar))
        return _multiply(*items)

    else:
        obj = identity(dim)
        for start, scalar in enumerate(spec):
            _set_value(obj, start, start, scalar)
        return obj  # eq. to gauss_jordan(dim, rank, eigenvalue_set, dim)


def jordan_normal(dim,
                  spec=None):
    r"""n x n matrix in Jordan normal form

    Explanation
    ===========
    A matrix in Jordan normal form is a block matrix
    with only non-zero blocks on the diagonal
    which have the form of an Jordan block $J$ which is

    .. math::

        J = \left[\begin{array}{cccccc}
            \lambda  & 1       & 0         & \dots     &         & 0     \\
             0       & \lambda & 1         & 0         & \dots   & 0     \\
             0       & \ddots  & \ddots    & \ddots    & \ddots  & 0     \\
             0       & \dots   & 0         & 0         & \lambda & 1     \\
             0       & \dots   &           &           & 0       & \lambda
             \end{array}\right]

    Finally, each $\lambda$ of each Jordan block will be an eigenvalue.

    To set the eigenvalues and block sizes use **spec**.

    As long as an eigenvalue is repeated in **spec**
    it mandates the same *Jordan* block.
    Such a *Jordan* block is interrupted by **None** entries.

    E.g. **spec** = (1,1,2, None, 2) would lead to

    .. math::

        J = \left[\begin{array}{cccc}
            1 & 1 & 0 & 0 \\
            0 & 1 & 0 & 0 \\
            0 & 0 & 2 & 0 \\
            0 & 0 & 0 & 2
        \end{array}\right]

    Alternativly, a sequence of pairs *(eigenvalue, block size)*,
    i.e. tuple of length 2, can be provided.

    So the above example is equivalent to
    **spec** = [(1,2),(2,1),(2,1)].

    .. math::

    Note, if the length of **spec** (**None** entries excluded)
    does not meet **dim**, blocks will be choosen randomly.

    To meet **dim** the final block might be truncated.

    Examples
    ========

    >>> from sympy.matrices.random import jordan_normal
    >>> import random
    >>> random.seed(1)

    >>> jordan_normal(3, spec=(2,2,2))
    Matrix([
    [2, 1, 0],
    [0, 2, 1],
    [0, 0, 2]])
    >>> jordan_normal(6, spec=(2,None,2,2,2,2,0))
    Matrix([
    [2, 0, 0, 0, 0, 0],
    [0, 2, 1, 0, 0, 0],
    [0, 0, 2, 1, 0, 0],
    [0, 0, 0, 2, 1, 0],
    [0, 0, 0, 0, 2, 0],
    [0, 0, 0, 0, 0, 0]])
    >>> # equivalent to
    >>> jordan_normal(6, spec=((2,1),(2,4),(0,1)))
    Matrix([
    [2, 0, 0, 0, 0, 0],
    [0, 2, 1, 0, 0, 0],
    [0, 0, 2, 1, 0, 0],
    [0, 0, 0, 2, 1, 0],
    [0, 0, 0, 0, 2, 0],
    [0, 0, 0, 0, 0, 0]])

    Parameters
    ==========
    dim : integer
        dimension n of matrix
    spec : tuple or list of symbols (optional)
        set of values of which scalars (diagonal entries) are choosen.

        If **dim** meets the length of **spec**,
        spec will be the eigenvalues (diagonal entries) as it is.

        If **dim** and the length of **spec** differs,
        the diagonal entries will be choosen randomly from **spec**.

        If not given **spec** defaults to $\{ -1, 1 \}$.

    See Also
    ========
    diagonal_normal

    """
    spec = spec or _elementary_scalar_set

    # make spec list of (scalar, size) tuples
    if spec and spec[0] is None:
        raise ValueError("spec argument must not start with None")

    block_list = list()
    val, cnt = None, 0
    for s in spec:
        if s == val:
            cnt += 1
            continue
        elif s is None:
            block_list.append((val, cnt))
            val, cnt = None, 0
            continue
        elif cnt:
            block_list.append((val, cnt))
        if isinstance(s, (tuple, list)) and len(s) == 2:
            block_list.append(s)
            val, cnt = None, 0
            continue
        val, cnt = s, 1
    if cnt:
        block_list.append((val, cnt))
    spec = block_list
    del block_list

    # choose block randomly if dim and len(spec) does not meet
    if not dim == sum(cnt for val, cnt in spec):
        spec = tuple(random.choice(spec) for _ in range(dim))

    # set entries of jordan blocks
    if _TEST:
        # multiplicative matrix construction (only for testing)
        items = list()
        start = 0
        for val, cnt in spec:
            end = min(dim, start + cnt)
            items.append(jordan(dim, index=(start, end), scalar=val))
            if end == dim:
                break
            start = end
        return _multiply(*items)

    else:
        obj = identity(dim)
        start = 0
        for scalar, size in spec:
            end = min(dim, start + size)
            _set_value(obj, start, start, scalar)
            for i in range(start + 1, end):
                _set_value(obj, i - 1, i, 1)
                _set_value(obj, i, i, scalar)
            if end == dim:
                break
            start = end
        return obj


def isometry_normal(dim,
                    spec=None):
    r""" isometry matrix n x n in normal form

    Explanation
    ===========

    Isometries preserve the geometric structure of vectorspaces.
    Let $<-,->$ be either the standard scalar product (over the reals)
    or the standard *Herminte* form (for complex vectorspaces).

    For vectors $x,y$ and a *isometry matrix* $A$ we have

    .. math::

        <Ax,Ay> = <x,y>

    Hence, $A^tA = I$ for real $A$ which is called *orthogonal*
    and $\bar{A}^tA = I$ for complex $A$ (called *unitary*).


    In normal from complex isometries have only diagonal entries of
    norm 1, i.e. for a diagonal element $z$ the following holds.

    .. math::

        |z| = z * \bar{z} = 1

    As - in normal form - real isometry (othogonal) matrices
    may have non diagonal entries but $2 \times 2$ rotation blocks

    .. math::

        \left[\begin{array}{cc} c & s \\ -s & c \end{array}\right]

    with $c^2 + s^2 = 1$.

    .. math::

    So, **spec** is assumed to be a list of entries which are either
    scalars $z$ with $|z| = z * \bar{z} = 1$
    or
    pairs $(c,s)$ with $c^2 + s^2 = 1$
    to form the above rotaion square.

    If a **spec** entry $c$ is neither such a pair nor $|c| = 1$
    a corresponding $s = \pm \sqrt{1-c^2}$ is choosen randomly.

    .. math::

    Note: Testing $|c| = 1$ might be not exact.
    To avoid such testing, provide $c$ as a 1-tuple in spec.

    In stead of

        **spec** = ( 1, -1, xi, (1/sqrt(2), 1/sqrt(2)) )

    checking abs(xi) == 1,

        **spec** = ((1,), (-1,), (xi,), (1/sqrt(2), 1/sqrt(2)))

    will not check abs(x) == 1.

    Examples
    ========

    >>> from sympy.matrices.random import isometry_normal
    >>> import random
    >>> random.seed(1)

    >>> isometry_normal(3)
    Matrix([
    [ sqrt(2)/2, sqrt(2)/2, 0],
    [-sqrt(2)/2, sqrt(2)/2, 0],
    [         0,         0, 1]])
    >>> isometry_normal(3, spec=(0.5, -0.5))
    Matrix([
    [              -0.5, 0.866025403784439, 0],
    [-0.866025403784439,              -0.5, 0],
    [                 0,                 0, 1]])
    >>> isometry_normal(3, spec=(-1, 1j, 1))
    Matrix([
    [-1,     0, 0],
    [ 0, 1.0*I, 0],
    [ 0,     0, 1]])
    >>> z = 1+1j
    >>> isometry_normal(3, spec=((-1,), (z/abs(z),), (1,)))
    Matrix([
    [-1,                                       0, 0],
    [ 0, 0.707106781186547 + 0.707106781186547*I, 0],
    [ 0,                                       0, 1]])

    Parameters
    ==========
    dim : integer
        dimension n of matrix
    spec : tuple or list of symbols (optional)
        set of values of which scalars (diagonal entries) are choosen.

        If **dim** meets the length of **spec**,
        spec will be the eigenvalues (diagonal entries) as it is.

        If **dim** and the length of **spec** differs,
        the diagonal entries will be choosen randomly from **spec**.

    See Also
    ========
    rotation
    diagonal_normal

    """
    spec = spec or _rotation_scalar_set

    # make spec list of (scalar,) or (scalar, +/-sqrt(1-scalar**2)) tuples
    block_list = list()
    for c in spec:
        if isinstance(c, (tuple, list)):
            block_list.append(c)
        elif abs(c) == 1:
            block_list.append((c,))
        else:
            block_list.append(_cs(c))
    spec = block_list
    del block_list

    # choose spec randomly if dim and len(spec) does not meet
    if not dim == sum(len(s) for s in spec):
        spec = tuple(random.choice(spec) for _ in range(dim))

    # set entries of rotation blocks

    if _TEST:
        # multiplicative matrix construction (only for testing)
        items = list()
        start = 0
        for s in spec:
            end = start + len(s)
            if dim < end:
                # leave the last diagonal entry one
                break
            if len(s) == 1:
                scalar, = s
                items.append(
                    elementary(dim, index=(start, start), scalar=scalar))
            else:
                c, s = s
                items.append(rotation(dim, (start, start + 1), (c, s)))
            start = end
        return _multiply(*items)

    else:
        obj = identity(dim)
        start = 0
        for s in spec:
            end = start + len(s)
            if dim < end:
                # leave the last diagonal entry one
                break
            if len(s) == 1:
                scalar, = s
                _set_value(obj, start, start, scalar)
            else:
                c, s = s
                _set_value(obj, start, start, c)
                _set_value(obj, start, start + 1, s)
                _set_value(obj, start + 1, start, -1 * s)
                _set_value(obj, start + 1, start + 1, c)
            start = end
        return obj


# === compound matrices, i.e. product of base matrices ===

def triangular(dim,
               rank=None,
               scalar_set=_elementary_scalar_set,
               unit_set=_elementary_scalar_set,
               length=None):
    r"""n x n upper triangular matrix with random entries.

    Explanation
    ===========
    Matrix with values placed above and on the diagonal. Constructed as
    product of random upper *triangular* :func:`elementary` matrices.

    It is constructed from an invertible triangular matrix $S$
    of which $r$ basis columns are randomly choosen
    and the $n-r$ columes are build as linear combinations
    of the basis colums.

    Examples
    ========
    >>> from sympy.matrices.random import triangular
    >>> import random
    >>> random.seed(1)

    >>> triangular(3, scalar_set=(2, -2, 0))
    Matrix([
    [-1, -2,  2],
    [ 0,  1, -2],
    [ 0,  0, -1]])

    Parameters
    ==========
    dim : integer
        dimension of matrix
    rank : integer (optional with default dim)
        rank of matrix
    scalar_set : tuple or list of symbols (optional)
        default values for random choosen scalars of
        elementary matrices to build the invertible matrix
    unit_set : tuple or list of symbols (optional)
        default values for random choosen scalar of
        elementary matrices to build the invertible matrix
        and non-zero diagonal entries
        which are assumed to have a multiplicative inverse
    length : integer (optional with default 2 * **dim**)
        number of invertible matrices to build the resulting matrix.

    See Also
    ========
    identity
    elementary
    square

    """
    if length == 0:
        return regular_to_singular(identity(dim), rank)
    scalar_set = scalar_set or (0,)
    unit_set = unit_set or (1,)
    length = length or 2 * dim
    row_indicies = [random.choice(range(dim)) for _ in range(length)]
    indicies = [(i, random.choice(range(i, dim))) for i in row_indicies]
    scalars = [random.choice(unit_set) if i == j else random.choice(scalar_set)
               for i, j in indicies]
    items = [elementary(dim, ix, s) for ix, s in zip(indicies, scalars)]
    return regular_to_singular(_multiply(*items), rank)


def square(dim,
           rank=None,
           scalar_set=_elementary_scalar_set,
           unit_set=_elementary_scalar_set,
           length=None):
    r"""a square matrix n x n

    Explanation
    ===========
    A square matrix represents an endomorphism
    in n dimensional vector space. It can be constructed as
    product of :func:`elementary` matrices.

    It is constructed from an :func:`invertible` matrix $S$
    of which $r$ basis columns are randomly choosen
    and the $n-r$ columes are build as linear combinations
    of the basis colums.

    Parameters
    ==========
    dim : integer
        dimension of matrix
    rank : integer (optional with default dim)
        rank of matrix
    scalar_set : tuple or list of symbols (optional)
        default values for random choosen scalar of
        elementary matrices to build the matrix
    unit_set : tuple or list of symbols (optional)
        default values for random choosen scalar of
        elementary matrices to build the matrix
        which are assumed to have a multiplicative inverse
    length : integer (optional with default 2 * **dim**)
        number of invertible matrices
        to build the resulting matrix.

    Examples
    ========
    >>> from sympy.matrices.random import square
    >>> import random
    >>> random.seed(1)

    >>> square(3)
    Matrix([
    [-1, -1, -1],
    [-2, -1, -1],
    [ 0,  0,  1]])

    See Also
    ========
    triangular
    invertible
    singular

    """
    if length == 0:
        return regular_to_singular(identity(dim), rank)

    length = length or 2 * dim
    lwr = triangular(dim, None, scalar_set, unit_set, int(length / 2))
    upr = triangular(dim, rank, scalar_set, unit_set, length - int(length / 2))
    return _multiply(_t(lwr), upr)


def invertible(dim,
               scalar_set=_elementary_scalar_set,
               unit_set=_elementary_scalar_set,
               length=None):
    r"""an invertible matrix n x n

    Explanation
    ===========
    An invertible matrix represents an isomorphism
    in n dimensional vector space. It can be constructed as
    product of :func:`elementary` matrices.

    Parameters
    ==========
    dim : integer
        dimension of matrix
    scalar_set : tuple or list of symbols (optional)
        default values for random choosen scalar of
        elementary matrices to build the matrix
    unit_set : tuple or list of symbols (optional)
        default values for random choosen scalar of
        elementary matrices to build the matrix
        which are assumed to have a multiplicative inverse
    length : integer (optional with default 2 * **dim**)
        number of invertible matrices
        to build the resulting matrix.

    Examples
    ========
    >>> from sympy import symbols, pprint
    >>> from sympy.matrices.random import invertible
    >>> import random
    >>> random.seed(1)

    >>> invertible(3)
    Matrix([
    [-1, -1, -1],
    [-2, -1, -1],
    [ 0,  0,  1]])

    >>> phi = symbols('phi')
    >>> m = invertible(3, scalar_set=(1, -1, phi, -phi))
    >>> m
    Matrix([
    [1, 0,     -phi],
    [1, 1, -phi - 1],
    [0, 0,        1]])
    >>> m.det()
    1
    >>> m.inv()
    Matrix([
    [ 1, 0, phi],
    [-1, 1,   1],
    [ 0, 0,   1]])

    """
    return square(dim, None, scalar_set, unit_set, length)


def singular(dim,
             rank=None,
             scalar_set=_elementary_scalar_set,
             unit_set=_elementary_scalar_set,
             length=None):
    r"""a singular square matrix n x n of a given rank.

    Explanation
    ===========
    A singular matrix is non-invertible and has determinant of zero.

    A matrix of rank $r$ has $r$ linear independend rows
    (i.e. a row space of dimension $r$).

    It is constructed from an :func:`invertible` matrix $S$
    of which $r$ basis columns are randomly choosen
    and the $n-r$ columes are build as linear combinations
    of the basis colums.

    Examples
    ========
    >>> from sympy.matrices.random import singular
    >>> import random
    >>> random.seed(1)

    >>> m = singular(3, 2)
    >>> m
    Matrix([
    [-1, 1, -1],
    [-2, 2, -1],
    [ 0, 0,  1]])
    >>> m.rank()
    2

    Parameters
    ==========
    dim : integer
        dimension of matrix
    rank : integer (optional with default dim-1)
        rank of matrix
    scalar_set : tuple or list of symbols (optional)
        default values for random choosen scalars of
        elementary matrices to build the matrix
    unit_set : tuple or list of symbols (optional)
        default values for random choosen scalar of
        elementary matrices to build the matrix
        which are assumed to have a multiplicative inverse
    length : integer (optional with default 2 * **dim**)
        number of invertible matrices to build the resulting matrix.

    See Also
    ========
    invertible
    triangular

    """
    rank = dim - 1 if rank is None else rank
    return square(dim, rank, scalar_set, unit_set, length)


# === conjugate matrices, i.e. defined by similarity to a normal from ==


def idempotent(dim,
               rank,
               scalar_set=_elementary_scalar_set,
               unit_set=_elementary_scalar_set,
               length=None):
    r"""an idempotent square matrix of a given rank s.th. $A*A=A$

    Explanation
    ===========
    An idempotent matrix is a matrix $A$ such that $A^2 = A$.
    It is constructed as product

    .. math::

        S \cdot P \cdot S^{-1}

    of an :func:`invertible` matrices $S$
    and a :func:`projection` matrix $P$ of given rank .

    Examples
    ========
    >>> from sympy.matrices.random import idempotent
    >>> import random
    >>> random.seed(1)

    >>> A = idempotent(3, 2)
    >>> A
    Matrix([
    [1, 0, 0],
    [0, 1, 1],
    [0, 0, 0]])
    >>> A.rank()
    2
    >>> A*A == A
    True

    Parameters
    ==========
    dim : integer
        dimension of matrix
    rank : integer (optional with default ``dim``)
        see :func:`projection`
    scalar_set : tuple or list of symbols (optional)
        see :func:`invertible`
    unit_set : tuple or list of symbols (optional)
        see :func:`invertible`
    length : integer (optional with default 2 * **dim**)
        see :func:`invertible`

    See Also
    ========
    elementary
    projection
    invertible
    nilpotent

    References
    ==========
    .. [1] https://en.wikipedia.org/wiki/Idempotence#Idempotent_functions

    """
    normal_form = projection(dim, (0, rank))
    s = invertible(dim, scalar_set, unit_set, length)
    return _multiply(_inv(s), normal_form, s)


def nilpotent(dim,
              rank,
              scalar_set=_elementary_scalar_set,
              unit_set=_elementary_scalar_set,
              length=None):
    r"""a nilpotent matrix of dinemsion n x n.

    Explanation
    ===========
    A *nilpotent* matrix is a matrix $A$ such that there is
    an integer $n$ with $A^n = 0$.
    It is build as the product

    .. math::

        S \cdot T \cdot S^{-1}

    where $S$ is an :func:`invertible`,
    $T$ a :func:`jordan_normal` matrix
    with zero diagonal entries.

    Examples
    ========
    >>> from sympy import zeros
    >>> from sympy.matrices.random import nilpotent
    >>> import random
    >>> random.seed(1)

    >>> A = nilpotent(3, 2)
    >>> A
    Matrix([
    [-2, -1, -2],
    [ 4,  2,  3],
    [ 0,  0,  0]])
    >>> A.rank()
    2
    >>> A*A*A == zeros(3)
    True

    Parameters
    ==========
    dim : integer
        dimension of matrix
    rank : integer (optional with default dim)
        rank of the matrix
    scalar_set : tuple or list of symbols (optional)
        see :func:`invertible`
    unit_set : tuple or list of symbols (optional)
        see :func:`invertible`
    length : integer (optional with default 2 * **dim**)
        see :func:`invertible`

    See Also
    ========
    elementary
    jordan_normal
    invertible
    idempotent

    References
    ==========
    .. [1] Lang, S., Algebra, DOI 10.1007/978-1-4613-0041-0, (2002)
    .. [2] https://en.wikipedia.org/wiki/nilpotent

    """

    def numbers_with_sum(n, k):
        """n numbers with sum dim"""
        if n == 1:
            return [k]
        num = random.choice(range(k - n)) + 1
        return [num] + numbers_with_sum(n - 1, k - num)

    index = numbers_with_sum(dim - rank, dim)
    spec = tuple((0, i) for i in index)
    normal_form = jordan_normal(dim, spec=spec)
    s = invertible(dim, scalar_set, unit_set, length)
    return _multiply(_inv(s), normal_form, s)


def diagonalizable(dim,
                   spec=_elementary_scalar_set,
                   scalar_set=_elementary_scalar_set,
                   unit_set=_elementary_scalar_set,
                   length=None):
    r"""a square matrix n x n of a given rank

    Explanation
    ===========
    Creates a square matrix of a given rank
    which is diagonalizable.

    The matrix will be a product

    .. math::

        S^{-1}*D*S

    of an :func:`invertible` matrix $S$ and
    a :func:`diagonal_normal` matrix of given rank
    and eigenvalues (diagonal entries) from **spec**.

    To obtain a matrix with *complex* *eigenvalues*
    but non-complex entries,
    e.g. $\lambda_{1/2}=a \pm b \ i$,
    you better may use a product of an invertible matrix $S$
    with non-complex entries (**scalar_set**)
    and where $D$ is a non-complex block diagonal matrix with a block $D_z$

    .. math::

        D_z = \left[\begin{array}{cc}a & b \\ -b & a\end{array}\right]

    Examples
    ========
    >>> from sympy.matrices.random import diagonalizable
    >>> import random
    >>> random.seed(1)

    >>> diagonalizable(3)
    Matrix([
    [-1,  0, 0],
    [-2, -3, 4],
    [-2, -2, 3]])

    >>> m = diagonalizable(3, spec=(1,2,2), scalar_set=(1,2,3), length=10)
    >>> m
    Matrix([
    [ 4,  2, -6],
    [-6, -4, 18],
    [-1, -1,  5]])
    >>> m.eigenvals(multiple=True)
    [1, 2, 2]
    >>> m.jordan_form(calc_transform=False)
    Matrix([
    [1, 0, 0],
    [0, 2, 0],
    [0, 0, 2]])

    Parameters
    ==========
    dim : integer
        dimension of matrix
    spec : tuple or list of symbols (optional)
        set of values of which scalars (diagonal entries) are choosen
        see :func:`diagonal_normal`
    scalar_set : tuple or list of symbols (optional)
        see :func:`invertible`
    unit_set : tuple or list of symbols (optional)
        see :func:`invertible`
    length : integer (optional with default 2 * **dim**)
        see :func:`invertible`

    """
    normal_form = diagonal_normal(dim, spec)
    s = invertible(dim, scalar_set, unit_set, length)
    return _multiply(_inv(s), normal_form, s)


def trigonalizable(dim,
                   spec=_elementary_scalar_set,
                   scalar_set=_elementary_scalar_set,
                   unit_set=_elementary_scalar_set,
                   length=None):
    r"""Creates a square matrix n x n of a given rank

    Explanation
    ===========
    Creates a square matrix n x n of a given rank
    which is trigonalizable (has echelon form)
    and has eigenvalues from a given set.

    The matrix will be a product

    .. math::

        S^{-1}*T*S

    of an invertible matrix $S$ and a *Jordan* matrix $T$
    and eigenvalues (diagonal entries) from **spec**.

    Examples
    ========
    >>> from sympy import sqrt
    >>> from sympy.matrices.random import trigonalizable
    >>> import random
    >>> random.seed(1)

    >>> trigonalizable(3)
    Matrix([
    [-1,  0, 0],
    [-2, -3, 4],
    [-2, -2, 3]])

    >>> m = trigonalizable(3, spec=(1, 1, 3), scalar_set=(1, sqrt(2), 2))
    >>> m
    Matrix([
    [            1,           1, 0],
    [            0,           1, 0],
    [2*sqrt(2) + 4, sqrt(2) + 2, 3]])
    >>> m.eigenvals(multiple=True)
    [1, 1, 3]
    >>> m.jordan_form(calc_transform=False)
    Matrix([
    [1, 1, 0],
    [0, 1, 0],
    [0, 0, 3]])

    Parameters
    ==========
    dim : integer
        dimension of matrix
    spec : tuple or list of symbols (optional)
        set of values of which scalars (diagonal entries) are choosen
        see :func:`jordan_normal`
    scalar_set : tuple or list of symbols (optional)
        see :func:`invertible`
    unit_set : tuple or list of symbols (optional)
        see :func:`invertible`
    length : integer (optional with default 2 * **dim**)
        see :func:`invertible`

    See Also
    ========
    jordan_normal
    invertible

    """
    normal_form = jordan_normal(dim, spec)
    s = invertible(dim, scalar_set, unit_set, length)
    return _multiply(_inv(s), normal_form, s)


# === matrices conjugate by isometries ==

def orthogonal(dim,
               spec=None,
               scalar_set=_rotation_scalar_set,
               length=None):
    r"""an orthogonal matrix n x n

    Explanation
    ===========
    An orthogonal matrix is an isometry
    in n dimensional Euclidian space. Constructed as
    product of random :func:`rotation`
    and a :func:`reflection` matrices.

    The values for random choosen scalar of rotations
    are from **scalar_set**.

    If **spec** is given,

    Examples
    ========
    >>> from sympy import sqrt, eye, expand
    >>> from sympy.matrices.random import orthogonal
    >>> import random
    >>> random.seed(1)

    >>> orthogonal(3, scalar_set=(-1,1))
    Matrix([
    [-1, 0,  0],
    [ 0, 1,  0],
    [ 0, 0, -1]])

    >>> s = sqrt(2)/2
    >>> random.seed(1)
    >>> orthogonal(3, scalar_set=(s,-s), length=2)
    Matrix([
    [sqrt(2)/2,       -1/2,      -1/2],
    [sqrt(2)/2,        1/2,       1/2],
    [        0, -sqrt(2)/2, sqrt(2)/2]])

    >>> random.seed(1)
    >>> o = orthogonal(3, spec=(-1, (s, s)), scalar_set=(s,), length=2)
    >>> expand(o * o.T) == eye(3)
    True

    Parameters
    ==========
    dim : integer
        dimension of matrix
    spec : tuple or list of symbols (optional),
        set of values of which scalars (diagonal entries)
        of :func:`isometry_normal` matrix
        see :func:`isometry_normal` for details
    scalar_set : tuple or list of symbols (optional),
        default values for random choosen scalar of
        rotations to build the orthogonal matrix
        :func:`rotation`
    length : integer (optional with default **dim**)
        number rotations to build the matrix.

    See Also
    ========
    rotation
    reflection
    isometry_normal
    unitary

    """
    if spec is None:
        if length == 0:
            return identity(dim)
        length = length or dim
        scalars = [random.choice(scalar_set) for _ in range(length)]
        items = [rotation(dim, scalar=s) for s in scalars]
        return _multiply(*items)
    normal_form = isometry_normal(dim, spec)
    s = orthogonal(dim, scalar_set=scalar_set, length=length)
    return _multiply(_t(s), normal_form, s)


def unitary(dim,
            spec=None,
            scalar_set=_unitary_scalar_set,
            length=None):
    r"""an unitary matrix n x n

    Explanation
    ===========
    An unitary matrix is an isometry
    in n dimensional unitary space (complex vectorspace with Hermite form).
    Constructed as product $A=UDV$ of two orthogonal matrices $U$ an $V$
    and a diagonal unitary martix $D$ (see :func:`isometry_normal`.
    All entries are taken from **scalar_set**.

    But if **spec** is given, another diagonal unitary matrix $D'$ is
    build with entries from **spec**.
    Note, to give an *unitary* matrix,
    **spec** must consist of roots of unity only.
    A root of unity is a complex number $z$,
    s.th. $|z| = z * \bar{z} = 1$.

    Finally, the resulting unitar matrix $B$ is given as

    .. math::

        B = \bar{A}^t D' A = \bar{V}^t \bar{D}^t \bar{U}^t D' UDV

    (see Fuehr/Rzeszotnik, "A note on factoring unitary matrices", 2018
    https://doi.org/10.1016%2Fj.laa.2018.02.017)

    Examples
    ========
    >>> from sympy import sqrt, I, simplify, eye, expand, re
    >>> from sympy.matrices.random import unitary
    >>> import random
    >>> random.seed(1)
    >>> u = unitary(3)

    >>> random.seed(1)
    >>> u = unitary(3)
    >>> simplify(u)
    Matrix([
    [ 0, -sqrt(2)*I/2, -sqrt(2)*I/2],
    [-I,            0,            0],
    [ 0,    1/2 + I/2,   -1/2 - I/2]])

    >>> u = unitary(3, scalar_set=(-1,))
    >>> simplify(u) # doctest: +SKIP
    Matrix([
    [-1, 0, 0],
    [ 0, 1, 0],
    [ 0, 0, 1]])

    >>> s = sqrt(2)/2
    >>> z = simplify(s + s * I)
    >>> z
    sqrt(2)*(1 + I)/2
    >>> abs(z)
    1
    >>> spec = [-1, z, 1]
    >>> unitary(3, spec=spec, length=0)
    Matrix([
    [-1,                 0, 0],
    [ 0, sqrt(2)*(1 + I)/2, 0],
    [ 0,                 0, 1]])

    >>> random.seed(1)
    >>> u = unitary(3, spec=spec)
    >>> simplify(u)
    Matrix([
    [sqrt(2)*(1 + I)/2,  0,  0],
    [                0,  0, -1],
    [                0, -1,  0]])
    >>> simplify(u * u.C)
    Matrix([
    [1, 0, 0],
    [0, 1, 0],
    [0, 0, 1]])

    >>> simplify(u.det()) == -z
    True
    >>> ev = [simplify(val) for val in u.eigenvals(multiple=True)]
    >>> ev = sorted(ev, key=(lambda x: re(x.evalf())))
    >>> ev
    [-1, sqrt(2)*(1 + I)/2, 1]
    >>> spec == ev
    True

    Parameters
    ==========
    dim : integer
        dimension of matrix
    spec : tuple or list of symbols (optional),
        set of values of which scalars (diagonal entries)
        of :func:`isometry_normal` matrix.
        see :func:`isometry_normal` for details
    scalar_set : tuple or list of symbols (optional),
        default values for random choosen scalar of
        rotations to build the orthogonal matrix
        see :func:`orthogonal`
    length : integer (optional with default 2 * **dim**)
        number rotations to build the matrix.
        see :func:`orthogonal`

    See Also
    ========
    isometry_normal
    orthogonal

    """
    if spec is None:
        if length == 0:
            return identity(dim)
        length = length or dim

        real_scalar_set = tuple(_cs(x) for x in scalar_set)
        complex_scalar_set = tuple(c * 1 + s * _i for c, s in real_scalar_set)
        scalars = [random.choice(complex_scalar_set) for _ in range(dim)]

        half = int(length/2)
        u = orthogonal(dim, scalar_set=real_scalar_set, length=half)
        d = diagonal_normal(dim, scalars)
        v = orthogonal(dim, scalar_set=real_scalar_set, length=length-half)
        return _multiply(u, d, v)

    normal_form = isometry_normal(dim, spec)
    s = unitary(dim, scalar_set=scalar_set, length=length)
    return _multiply(_c(s), normal_form, s)


def normal(dim,
           spec=_elementary_scalar_set,
           scalar_set=None,
           length=None):
    r""" normal n x n matrix

    Explanation
    ===========
    By definition for a *normal* matrix $A$

    .. math::

        \bar{A}^t A = A \bar{A}^t

    holds. Note, this extends the notion of an :func:`unitary` matrix $U$
    since $\bar{U}^t U = I = U \bar{U}^t$.

    A matrix $A$ is *normal* if and only if it is diagonalizable
    by an *unitary* matrix $U$,
    i.e. $\bar{U}^t A U = D$ with diagonal matrix $D$.

    Hence, since $\bar{U}^t = U^{-1}$,

    .. math::

        A = U D \bar{U}^t

    The entries of $D$ are taken from **spec**
    and the entries to build $U$ are taken from **scalar_set**.

    Since **orthogonal** matrices are **unitary**,
    $U$ will be **orthogonal**
    if **scalar_set** has only eal entries.

    So the final *normal* matrix will be real
    if **spec** consists only of real entries, too.

    Examples
    ========
    >>> from sympy import sqrt, I
    >>> from sympy.matrices.random import normal
    >>> import random
    >>> random.seed(1)

    >>> normal(3)
    Matrix([
    [ 0, -1,  0],
    [-1,  0,  0],
    [ 0,  0, -1]])

    >>> z = complex(1, 1)
    >>> random.seed(1)
    >>> n = normal(2, spec=(2, 3),  scalar_set=(sqrt(2)/2,))
    >>> n.T * n == n * n.T
    True

    Parameters
    ==========
    dim : integer
        dimension of matrix
    spec : tuple or list of symbols (optional)
        set of values of which scalars (diagonal entries)
        see :func:`isometry_normal` for details
    scalar_set : tuple or list of symbols (optional)
        see :func:`orthogonal` and :func:`unitary`
    length : integer (optional with default **dim**)
        see :func:`orthogonal` and :func:`unitary`

    See Also
    ========
    isometry_normal
    orthogonal
    unitary

    """
    normal_form = diagonal_normal(dim, spec)
    if any(isinstance(s, (list, tuple)) for s in spec):
        scalar_set = scalar_set or _rotation_scalar_set
        s = orthogonal(dim, None, scalar_set, length)
    elif all(s == _c(s) for s in spec):
        scalar_set = scalar_set or _rotation_scalar_set
        s = orthogonal(dim, None, scalar_set, length)
    else:
        scalar_set = scalar_set or _unitary_scalar_set
        s = unitary(dim, None, scalar_set, length)
    return _multiply(_c(s), normal_form, s)


# === symmetric or complex adjoined matrices ===

def symmetric(dim,
              scalar_set=_elementary_scalar_set,
              unit_set=_elementary_scalar_set,
              length=None):
    r"""Creates a symmetric square matrix n x n of a given rank.

    Explanation
    ===========
    A *symmetric* matrix is a matrix such that $A^t = A$
    The matrix is simply constructed as

    .. math::

        S*D*S^t

    by an ``invertible`` matrix $S$ and a ``diagonal`` matrix $D$ of given rank.

    To obtain a matrix with given *eigenvalues*,
    you better may use :func:`normal` matrix, which is a product

    .. math::

        O*D*O^t

    by an *orthogonal* matrix $O$ and
    a *diagonal* matrix $D$ with given eigenvalues (diagonal entries).

    Examples
    ========
    >>> from sympy.matrices.random import symmetric
    >>> import random
    >>> random.seed(1)

    >>> symmetric(3)
    Matrix([
    [5, 3, 3],
    [3, 2, 2],
    [3, 2, 3]])

    For a symmetric matrix with given eigenvalues
    :func:`normal` works

    >>> from sympy.matrices.random import normal
    >>> import random
    >>> random.seed(1)

    >>> n = normal(3, spec=(1,2,3))
    >>> sorted(n.eigenvals(multiple=True))
    [1, 2, 3]

    Parameters
    ==========
    dim : integer
        dimension of matrix
    scalar_set : tuple or list of symbols (optional)
        see :func:`invertible`
    unit_set : tuple or list of symbols (optional)
        see :func:`invertible`
    length : integer (optional with default 2 * **dim**)
        see :func:`invertible`

    See Also
    ========
    invertible
    normal
    hermite

    """
    s = invertible(dim, scalar_set, unit_set, length)
    return _multiply(_t(s), s)


def hermite(dim,
            scalar_set=_elementary_scalar_set,
            unit_set=_elementary_scalar_set,
            length=None):
    r"""

    Explanation
    ===========

    Examples
    ========
    >>> from sympy.matrices.random import hermite
    >>> import random
    >>> random.seed(1)

    >>> z = complex(1,-1)
    >>> hermite(3, (z, ), length=1)
    Matrix([
    [          1, 0,                     1.0 - 1.0*I],
    [          0, 1,                               0],
    [1.0 + 1.0*I, 0, 1 + (1.0 - 1.0*I)*(1.0 + 1.0*I)]])

    Parameters
    ==========
    dim : integer
        dimension of matrix
    scalar_set : tuple or list of symbols (optional)
        see :func:`invertible`
    unit_set : tuple or list of symbols (optional)
        see :func:`invertible`
    length : integer (optional with default 2 * **dim**)
        see :func:`invertible`

    See Also
    ========
    invertible
    normal
    hermite

    """
    s = invertible(dim, scalar_set, unit_set, length)
    return _multiply(_c(s), s)
