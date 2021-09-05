"""Known matrices related to physics"""

from sympy import Matrix, I, pi, sqrt
from sympy.functions import exp
from sympy.core.decorators import deprecated


def msigma(i):
    r"""Returns a Pauli matrix `\sigma_i` with ``i=1,2,3``.

    References
    ==========

    .. [1] https://en.wikipedia.org/wiki/Pauli_matrices

    Examples
    ========

    >>> from sympy.physics.matrices import msigma
    >>> msigma(1)
    Matrix([
    [0, 1],
    [1, 0]])
    """
    if i == 1:
        mat = ( (
            (0, 1),
            (1, 0)
        ) )
    elif i == 2:
        mat = ( (
            (0, -I),
            (I, 0)
        ) )
    elif i == 3:
        mat = ( (
            (1, 0),
            (0, -1)
        ) )
    else:
        raise IndexError("Invalid Pauli index")
    return Matrix(mat)


def pat_matrix(m, dx, dy, dz):
    """Returns the Parallel Axis Theorem matrix to translate the inertia
    matrix a distance of `(dx, dy, dz)` for a body of mass m.

    Examples
    ========

    To translate a body having a mass of 2 units a distance of 1 unit along
    the `x`-axis we get:

    >>> from sympy.physics.matrices import pat_matrix
    >>> pat_matrix(2, 1, 0, 0)
    Matrix([
    [0, 0, 0],
    [0, 2, 0],
    [0, 0, 2]])

    """
    dxdy = -dx*dy
    dydz = -dy*dz
    dzdx = -dz*dx
    dxdx = dx**2
    dydy = dy**2
    dzdz = dz**2
    mat = ((dydy + dzdz, dxdy, dzdx),
           (dxdy, dxdx + dzdz, dydz),
           (dzdx, dydz, dydy + dxdx))
    return m*Matrix(mat)


def mgamma(mu, lower=False):
    r"""Returns a Dirac gamma matrix `\gamma^\mu` in the standard
    (Dirac) representation.

    Explanation
    ===========

    If you want `\gamma_\mu`, use ``gamma(mu, True)``.

    We use a convention:

    `\gamma^5 = i \cdot \gamma^0 \cdot \gamma^1 \cdot \gamma^2 \cdot \gamma^3`

    `\gamma_5 = i \cdot \gamma_0 \cdot \gamma_1 \cdot \gamma_2 \cdot \gamma_3 = - \gamma^5`

    References
    ==========

    .. [1] https://en.wikipedia.org/wiki/Gamma_matrices

    Examples
    ========

    >>> from sympy.physics.matrices import mgamma
    >>> mgamma(1)
    Matrix([
    [ 0,  0, 0, 1],
    [ 0,  0, 1, 0],
    [ 0, -1, 0, 0],
    [-1,  0, 0, 0]])
    """
    if not mu in [0, 1, 2, 3, 5]:
        raise IndexError("Invalid Dirac index")
    if mu == 0:
        mat = (
            (1, 0, 0, 0),
            (0, 1, 0, 0),
            (0, 0, -1, 0),
            (0, 0, 0, -1)
        )
    elif mu == 1:
        mat = (
            (0, 0, 0, 1),
            (0, 0, 1, 0),
            (0, -1, 0, 0),
            (-1, 0, 0, 0)
        )
    elif mu == 2:
        mat = (
            (0, 0, 0, -I),
            (0, 0, I, 0),
            (0, I, 0, 0),
            (-I, 0, 0, 0)
        )
    elif mu == 3:
        mat = (
            (0, 0, 1, 0),
            (0, 0, 0, -1),
            (-1, 0, 0, 0),
            (0, 1, 0, 0)
        )
    elif mu == 5:
        mat = (
            (0, 0, 1, 0),
            (0, 0, 0, 1),
            (1, 0, 0, 0),
            (0, 1, 0, 0)
        )
    m = Matrix(mat)
    if lower:
        if mu in [1, 2, 3, 5]:
            m = -m
    return m

#Minkowski tensor using the convention (+,-,-,-) used in the Quantum Field
#Theory
minkowski_tensor = Matrix( (
    (1, 0, 0, 0),
    (0, -1, 0, 0),
    (0, 0, -1, 0),
    (0, 0, 0, -1)
))


@deprecated(issue=20246, useinstead="DFT(n).as_mutable(), DFT(n), DFT(n).as_explicit()",
            deprecated_since_version="1.9")
def mdft(n):
    r"""
    Deprecated. Use DFT from sympy.matrices.expressions.fourier instead.

    To get identical behavior to ``mdft(n)``, use ``DFT(n).as_mutable()``.
    """
    mat = [[None for x in range(n)] for y in range(n)]
    base = exp(-2*pi*I/n)
    mat[0] = [1]*n
    for i in range(n):
        mat[i][0] = 1
    for i in range(1, n):
        for j in range(i, n):
            mat[i][j] = mat[j][i] = base**(i*j)
    return (1/sqrt(n))*Matrix(mat)
