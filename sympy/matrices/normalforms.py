'''Functions returning normal forms of matrices'''

from sympy.core.numbers import igcdex
from sympy.polys.polytools import Poly
from sympy.polys.matrices import DomainMatrix
from sympy.polys.matrices.normalforms import (
        smith_normal_form as _snf,
        invariant_factors as _invf,
    )


def _to_domain(m, domain=None):
    """Convert Matrix to DomainMatrix"""
    # XXX: deprecated support for RawMatrix:
    ring = getattr(m, "ring", None)
    m = m.applyfunc(lambda e: e.as_expr() if isinstance(e, Poly) else e)

    dM = DomainMatrix.from_Matrix(m)

    domain = domain or ring
    if domain is not None:
        dM = dM.convert_to(domain)
    return dM


def smith_normal_form(m, domain=None):
    '''
    Return the Smith Normal Form of a matrix `m` over the ring `domain`.
    This will only work if the ring is a principal ideal domain.

    Examples
    ========

    >>> from sympy import Matrix, ZZ
    >>> from sympy.matrices.normalforms import smith_normal_form
    >>> m = Matrix([[12, 6, 4], [3, 9, 6], [2, 16, 14]])
    >>> print(smith_normal_form(m, domain=ZZ))
    Matrix([[1, 0, 0], [0, 10, 0], [0, 0, -30]])

    '''
    dM = _to_domain(m, domain)
    return _snf(dM).to_Matrix()


def invariant_factors(m, domain=None):
    '''
    Return the tuple of abelian invariants for a matrix `m`
    (as in the Smith-Normal form)

    References
    ==========

    [1] https://en.wikipedia.org/wiki/Smith_normal_form#Algorithm
    [2] http://sierra.nmsu.edu/morandi/notes/SmithNormalForm.pdf

    '''
    dM = _to_domain(m, domain)
    factors = _invf(dM)
    factors = tuple(dM.domain.to_sympy(f) for f in factors)
    # XXX: deprecated.
    if hasattr(m, "ring"):
        if m.ring.is_PolynomialRing:
            K = m.ring
            to_poly = lambda f: Poly(f, K.symbols, domain=K.domain)
            factors = tuple(to_poly(f) for f in factors)
    return factors


def _gcdex(a, b):
    """
    This supports the functions that compute Hermite Normal Form.

    Let x, y be the coefficients returned by the extended Euclidean
    Algorithm, so that x*a + y*b = g. In the algorithms for computing HNF,
    it is critical that x, y not only satisfy the condition of being small
    in magnitude -- namely that |x| <= |b|/g, |y| <- |a|/g -- but also that
    y == 0 when a | b.
    """
    x, y, g = igcdex(a, b)
    if a != 0 and b % a == 0:
        y = 0
        x = -1 if a < 0 else 1
    return x, y, g


def _hermite_normal_form(A):
    '''
    Compute the Hermite Normal Form of a matrix `A`.
    '''
    m = A.rows
    n = A.cols
    i = m - 1
    j = n - 1
    k = j
    ell = 0 if m <= n else m - n

    state = 2
    while True:
        if state == 2:
            if j == 0:
                state = 4
            else:
                j -= 1
                if A[i, j] != 0:
                    state = 3
        elif state == 3:
            u, v, d = _gcdex(A[i, k], A[i, j])
            B = u * A.col(k) + v * A.col(j)
            r, s = A[i, k] // d, A[i, j] // d
            for ii in range(m):
                A[ii, j] = r * A[ii, j] - s * A[ii, k]
                A[ii, k] = B[ii, 0]
            state = 2
        else:
            assert state == 4
            b = A[i, k]
            if b < 0:
                for ii in range(m):
                    A[ii, k] = -A[ii, k]
                b = -b
            if b == 0:
                k += 1
            else:
                for jj in range(k+1, n):
                    q = A[i, jj] // b
                    for ii in range(m):
                        A[ii, jj] -= q*A[ii, k]
            if i == ell:
                W = A.zeros(m, n - k)
                for jj in range(n - k):
                    for ii in range(m):
                        W[ii, jj] = A[ii, jj + k]
                return W
            else:
                i -= 1
                k -= 1
                j = k
                state = 2


# FIXME:
#  The modulo D algorithm is broken.
#  E.g. let it be used by `sympy.polys.numberfields.basis.round_two`,
#  and you will see bad results!
def _hermite_normal_form_modulo_D(A, D):
    '''
    Compute the Hermite Normal Form `W` of an m x n matrix `A` having integer
    coefficients, and rank m, given a positive integer `D` known in advance to
    be a multiple of `det(W)`.
    '''
    m = A.rows
    n = A.cols
    i = m - 1
    j = n - 1
    k = j
    R = D

    # FIXME:
    #  What if m == 1?
    #  Alg 2.4.8 in Cohen does not handle that case.
    #  Does our modification to Step 4 fix this case too?
    W = A.zeros(m, m)

    state = 2
    while True:
        if state == 2:
            # We will want residues mod R in (-R/2, R/2], not [0, R):
            # FIXME:
            #  Use `sympy.ntheory.modular.symmetric_residue` instead?
            res_adjust = (R - 2 + (R % 2)) // 2
            if j == 0:
                state = 4
            else:
                j -= 1
                if A[i, j] != 0:
                    state = 3
        elif state == 3:
            u, v, d = _gcdex(A[i, k], A[i, j])
            B = u*A.col(k) + v*A.col(j)
            r, s = A[i, k] // d, A[i, j] // d
            for ii in range(m):
                if R > 1:
                    A[ii, j] = (r*A[ii, j] - s*A[ii, k]) % R - res_adjust
                    A[ii, k] = B[ii, 0] % R - res_adjust
                else:
                    A[ii, j] = r * A[ii, j] - s * A[ii, k]
                    A[ii, k] = B[ii, 0]
            state = 2
        else:
            assert state == 4
            u, v, d = _gcdex(A[i, k], R)
            for ii in range(m):
                if R > 1:
                    # This time we want residue in [0, R) except on the diagonal,
                    # where 0 must be replaced with R.
                    r = u * A[ii, k] % R
                    if r > 0 or ii != i:
                        W[ii, i] = r
                    else:
                        W[ii, i] = R
                else:
                    W[ii, i] = u * A[ii, k]
            for jj in range(i + 1, m):
                q = W[i, jj] // W[i, i]
                for ii in range(m):
                    W[ii, jj] -= q * W[ii, i]
            if i == 0:
                break
            else:
                R //= d
                i -= 1
                k -= 1
                j = k
                if A[i, k] == 0:
                    A[i, k] = R
                state = 2

    return W


def hermite_normal_form(A, D=None):
    '''
    Return the Hermite Normal Form `W` of a matrix `A` all of whose entries
    are integers.

    If known in advance, a positive integer `D` being any multiple of `det(W)`
    may be provided. In this case, if `A` also has rank equal to its number of
    rows, then we may use an alternative algorithm that prevents coefficient
    explosion.

    References
    ==========

    [1] Cohen, H. A Course in Computational Algebraic Number Theory.

    '''
    if D is not None and A.rank() == A.rows:
        # FIXME:
        #  Shutting off the mod D option for now, since it is still broken...
        #return _hermite_normal_form_modulo_D(A, D)
        return _hermite_normal_form(A)
    else:
        return _hermite_normal_form(A)
