"""

Module for the ddm_* routines for operating on a matrix in list of lists
matrix representation.

These routines are used internally by the DDM class which also provides a
friendlier interface for them. The idea here is to implement core matrix
routines in a way that can be applied to any simple list representation
without the need to use any particular matrix class. For example we can
compute the RREF of a matrix like:

    >>> from sympy.polys.matrices.dense import ddm_irref
    >>> M = [[1, 2, 3], [4, 5, 6]]
    >>> pivots = ddm_irref(M)
    >>> M
    [[1.0, 0.0, -1.0], [0, 1.0, 2.0]]

These are lower-level routines that work mostly in place.The routines at this
level should not need to know what the domain of the elements is but should
ideally document what operations they will use and what functions they need to
be provided with.

The next-level up is the DDM class which uses these routines but wraps them up
with an interface that handles copying etc and keeps track of the Domain of
the elements of the matrix:

    >>> from sympy.polys.domains import QQ
    >>> from sympy.polys.matrices.ddm import DDM
    >>> M = DDM([[QQ(1), QQ(2), QQ(3)], [QQ(4), QQ(5), QQ(6)]], (2, 3), QQ)
    >>> M
    [[1, 2, 3], [4, 5, 6]]
    >>> Mrref, pivots = M.rref()
    >>> Mrref
    [[1, 0, -1], [0, 1, 2]]

"""
from __future__ import annotations
from operator import mul
from .exceptions import (
    DMShapeError,
    DMNonInvertibleMatrixError,
    DMNonSquareMatrixError,
)
from typing import Sequence, TypeVar
from sympy.polys.matrices._typing import RingElement


T = TypeVar('T')
R = TypeVar('R', bound=RingElement)


def ddm_transpose(matrix: Sequence[Sequence[T]]) -> list[list[T]]:
    """matrix transpose"""
    return list(map(list, zip(*matrix)))


def ddm_iadd(a: list[list[R]], b: Sequence[Sequence[R]]) -> None:
    """a += b"""
    for ai, bi in zip(a, b):
        for j, bij in enumerate(bi):
            ai[j] += bij


def ddm_isub(a: list[list[R]], b: Sequence[Sequence[R]]) -> None:
    """a -= b"""
    for ai, bi in zip(a, b):
        for j, bij in enumerate(bi):
            ai[j] -= bij


def ddm_ineg(a: list[list[R]]) -> None:
    """a  <--  -a"""
    for ai in a:
        for j, aij in enumerate(ai):
            ai[j] = -aij


def ddm_imul(a: list[list[R]], b: R) -> None:
    for ai in a:
        for j, aij in enumerate(ai):
            ai[j] = aij * b


def ddm_irmul(a: list[list[R]], b: R) -> None:
    for ai in a:
        for j, aij in enumerate(ai):
            ai[j] = b * aij


def ddm_imatmul(
    a: list[list[R]], b: Sequence[Sequence[R]], c: Sequence[Sequence[R]]
) -> None:
    """a += b @ c"""
    cT = list(zip(*c))

    for bi, ai in zip(b, a):
        for j, cTj in enumerate(cT):
            ai[j] = sum(map(mul, bi, cTj), ai[j])


def ddm_irref(a, _partial_pivot=False):
    """a  <--  rref(a)"""
    # a is (m x n)
    m = len(a)
    if not m:
        return []
    n = len(a[0])

    i = 0
    pivots = []

    for j in range(n):
        # Proper pivoting should be used for all domains for performance
        # reasons but it is only strictly needed for RR and CC (and possibly
        # other domains like RR(x)). This path is used by DDM.rref() if the
        # domain is RR or CC. It uses partial (row) pivoting based on the
        # absolute value of the pivot candidates.
        if _partial_pivot:
            ip = max(range(i, m), key=lambda ip: abs(a[ip][j]))
            a[i], a[ip] = a[ip], a[i]

        # pivot
        aij = a[i][j]

        # zero-pivot
        if not aij:
            for ip in range(i+1, m):
                aij = a[ip][j]
                # row-swap
                if aij:
                    a[i], a[ip] = a[ip], a[i]
                    break
            else:
                # next column
                continue

        # normalise row
        ai = a[i]
        aijinv = aij**-1
        for l in range(j, n):
            ai[l] *= aijinv # ai[j] = one

        # eliminate above and below to the right
        for k, ak in enumerate(a):
            if k == i or not ak[j]:
                continue
            akj = ak[j]
            ak[j] -= akj # ak[j] = zero
            for l in range(j+1, n):
                ak[l] -= akj * ai[l]

        # next row
        pivots.append(j)
        i += 1

        # no more rows?
        if i >= m:
            break

    return pivots


def ddm_idet(a, K):
    """a  <--  echelon(a); return det"""
    # Bareiss algorithm
    # https://www.math.usm.edu/perry/Research/Thesis_DRL.pdf

    # a is (m x n)
    m = len(a)
    if not m:
        return K.one
    n = len(a[0])

    exquo = K.exquo
    # uf keeps track of the sign change from row swaps
    uf = K.one

    for k in range(n-1):
        if not a[k][k]:
            for i in range(k+1, n):
                if a[i][k]:
                    a[k], a[i] = a[i], a[k]
                    uf = -uf
                    break
            else:
                return K.zero

        akkm1 = a[k-1][k-1] if k else K.one

        for i in range(k+1, n):
            for j in range(k+1, n):
                a[i][j] = exquo(a[i][j]*a[k][k] - a[i][k]*a[k][j], akkm1)

    return uf * a[-1][-1]


def ddm_iinv(ainv, a, K):
    if not K.is_Field:
        raise ValueError('Not a field')

    # a is (m x n)
    m = len(a)
    if not m:
        return
    n = len(a[0])
    if m != n:
        raise DMNonSquareMatrixError

    eye = [[K.one if i==j else K.zero for j in range(n)] for i in range(n)]
    Aaug = [row + eyerow for row, eyerow in zip(a, eye)]
    pivots = ddm_irref(Aaug)
    if pivots != list(range(n)):
        raise DMNonInvertibleMatrixError('Matrix det == 0; not invertible.')
    ainv[:] = [row[n:] for row in Aaug]


def ddm_ilu_split(L, U, K):
    """L, U  <--  LU(U)"""
    m = len(U)
    if not m:
        return []
    n = len(U[0])

    swaps = ddm_ilu(U)

    zeros = [K.zero] * min(m, n)
    for i in range(1, m):
        j = min(i, n)
        L[i][:j] = U[i][:j]
        U[i][:j] = zeros[:j]

    return swaps


def ddm_ilu(a):
    """a  <--  LU(a)"""
    m = len(a)
    if not m:
        return []
    n = len(a[0])

    swaps = []

    for i in range(min(m, n)):
        if not a[i][i]:
            for ip in range(i+1, m):
                if a[ip][i]:
                    swaps.append((i, ip))
                    a[i], a[ip] = a[ip], a[i]
                    break
            else:
                # M = Matrix([[1, 0, 0, 0], [0, 0, 0, 0], [0, 0, 1, 1], [0, 0, 1, 2]])
                continue
        for j in range(i+1, m):
            l_ji = a[j][i] / a[i][i]
            a[j][i] = l_ji
            for k in range(i+1, n):
                a[j][k] -= l_ji * a[i][k]

    return swaps


def ddm_ilu_solve(x, L, U, swaps, b):
    """x  <--  solve(L*U*x = swaps(b))"""
    m = len(U)
    if not m:
        return
    n = len(U[0])

    m2 = len(b)
    if not m2:
        raise DMShapeError("Shape mismtch")
    o = len(b[0])

    if m != m2:
        raise DMShapeError("Shape mismtch")
    if m < n:
        raise NotImplementedError("Underdetermined")

    if swaps:
        b = [row[:] for row in b]
        for i1, i2 in swaps:
            b[i1], b[i2] = b[i2], b[i1]

    # solve Ly = b
    y = [[None] * o for _ in range(m)]
    for k in range(o):
        for i in range(m):
            rhs = b[i][k]
            for j in range(i):
                rhs -= L[i][j] * y[j][k]
            y[i][k] = rhs

    if m > n:
        for i in range(n, m):
            for j in range(o):
                if y[i][j]:
                    raise DMNonInvertibleMatrixError

    # Solve Ux = y
    for k in range(o):
        for i in reversed(range(n)):
            if not U[i][i]:
                raise DMNonInvertibleMatrixError
            rhs = y[i][k]
            for j in range(i+1, n):
                rhs -= U[i][j] * x[j][k]
            x[i][k] = rhs / U[i][i]


def ddm_berk(M, K):
    m = len(M)
    if not m:
        return [[K.one]]
    n = len(M[0])

    if m != n:
        raise DMShapeError("Not square")

    if n == 1:
        return [[K.one], [-M[0][0]]]

    a = M[0][0]
    R = [M[0][1:]]
    C = [[row[0]] for row in M[1:]]
    A = [row[1:] for row in M[1:]]

    q = ddm_berk(A, K)

    T = [[K.zero] * n for _ in range(n+1)]
    for i in range(n):
        T[i][i] = K.one
        T[i+1][i] = -a
    for i in range(2, n+1):
        if i == 2:
            AnC = C
        else:
            C = AnC
            AnC = [[K.zero] for row in C]
            ddm_imatmul(AnC, A, C)
        RAnC = [[K.zero]]
        ddm_imatmul(RAnC, R, AnC)
        for j in range(0, n+1-i):
            T[i+j][j] = -RAnC[0][0]

    qout = [[K.zero] for _ in range(n+1)]
    ddm_imatmul(qout, T, q)
    return qout
