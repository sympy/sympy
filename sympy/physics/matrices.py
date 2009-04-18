"""Known matrices related to physics"""

from sympy import Matrix, I

def msigma(i):
    """Returns a Pauli matrix sigma_i. i=1,2,3
    See also:
    http://en.wikipedia.org/wiki/Pauli_matrices
    """
    if i==1:
        mat=( (
            (0, 1),
            (1, 0)
            ) )
    elif i==2:
        mat=( (
            (0, -I),
            (I, 0)
            ) )
    elif i==3:
        mat=( (
            (1, 0),
            (0, -1)
            ) )
    else:
        raise IndexError("Invalid Pauli index")
    return Matrix(mat)

def mgamma(mu,lower=False):
    """Returns a Dirac gamma matrix gamma^mu in the standard
    (Dirac) representation.

    If you want gamma_mu, use gamma(mu, True).

    We use a convention:

    gamma^5 = I * gamma^0 * gamma^1 * gamma^2 * gamma^3
    gamma_5 = I * gamma_0 * gamma_1 * gamma_2 * gamma_3 = - gamma^5

    See also:

    http://en.wikipedia.org/wiki/Gamma_matrices

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

