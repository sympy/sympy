"""Known matrices related to physics"""

from sympy import Matrix, I

def msigma(i):
    """Returns a Pauli matrix sigma_i. i=1,2,3

    See Also
    ========

    http://en.wikipedia.org/wiki/Pauli_matrices

    Examples
    ========

    >>> from sympy.physics.matrices import msigma
    >>> msigma(1)
    [0, 1]
    [1, 0]
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
    matrix a distance of (dx, dy, dz) for a body of mass m.

    Examples
    --------
    If the point we want the inertia about is a distance of 2 units of
    length and 1 unit along the x-axis we get:
    >>> from sympy.physics.matrices import pat_matrix
    >>> pat_matrix(2,1,0,0)
    [0, 0, 0]
    [0, 2, 0]
    [0, 0, 2]

    In case we want to find the inertia along a vector of (1,1,1):
    >>> pat_matrix(2,1,1,1)
    [ 4, -2, -2]
    [-2,  4, -2]
    [-2, -2,  4]
    """
    dxdy = -dx*dy ; dydz = -dy*dz ; dzdx = -dz*dx
    dxdx =  dx**2 ; dydy =  dy**2 ; dzdz =  dz**2
    mat = ((dydy + dzdz, dxdy, dzdx),
           (dxdy, dxdx + dzdz, dydz),
           (dzdx, dydz, dydy + dxdx))
    return m*Matrix(mat)

def mgamma(mu,lower=False):
    """Returns a Dirac gamma matrix gamma^mu in the standard
    (Dirac) representation.

    If you want gamma_mu, use gamma(mu, True).

    We use a convention:

    gamma^5 = I * gamma^0 * gamma^1 * gamma^2 * gamma^3
    gamma_5 = I * gamma_0 * gamma_1 * gamma_2 * gamma_3 = - gamma^5

    See Also
    ========

    http://en.wikipedia.org/wiki/Gamma_matrices

    Examples
    ========

    >>> from sympy.physics.matrices import mgamma
    >>> mgamma(1)
    [ 0,  0, 0, 1]
    [ 0,  0, 1, 0]
    [ 0, -1, 0, 0]
    [-1,  0, 0, 0]
    """
    if not mu in [0,1,2,3,5]:
        raise IndexError("Invalid Dirac index")
    if mu == 0:
        mat = (
                (1,0,0,0),
                (0,1,0,0),
                (0,0,-1,0),
                (0,0,0,-1)
                )
    elif mu == 1:
        mat = (
                (0,0,0,1),
                (0,0,1,0),
                (0,-1,0,0),
                (-1,0,0,0)
                )
    elif mu == 2:
        mat = (
                (0,0,0,-I),
                (0,0,I,0),
                (0,I,0,0),
                (-I,0,0,0)
                )
    elif mu == 3:
        mat = (
                (0,0,1,0),
                (0,0,0,-1),
                (-1,0,0,0),
                (0,1,0,0)
                )
    elif mu == 5:
        mat = (
                (0,0,1,0),
                (0,0,0,1),
                (1,0,0,0),
                (0,1,0,0)
                )
    m= Matrix(mat)
    if lower:
        if mu in [1,2,3,5]:
            m = - m
    return m

#Minkowski tensor using the convention (+,-,-,-) used in the Quantum Field
#Theory
minkowski_tensor = Matrix( (
    (1,0,0,0),
    (0,-1,0,0),
    (0,0,-1,0),
    (0,0,0,-1)
    ))
