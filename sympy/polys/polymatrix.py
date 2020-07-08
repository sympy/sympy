from __future__ import print_function

from sympy.core.sympify import _sympify

from sympy.matrices.dense import MutableDenseMatrix
from sympy.matrices import NonSquareMatrixError
from sympy.matrices.common import NonInvertibleMatrixError, ShapeError

from sympy.polys.constructor import construct_domain
from sympy.polys.polytools import Poly
from sympy.polys.domains import EX


class MutablePolyDenseMatrix(MutableDenseMatrix):
    """
    A mutable matrix of objects from poly module or to operate with them.

    Examples
    ========

    >>> from sympy.polys.polymatrix import PolyMatrix
    >>> from sympy import Symbol, Poly, ZZ
    >>> x = Symbol('x')
    >>> pm1 = PolyMatrix([[Poly(x**2, x), Poly(-x, x)], [Poly(x**3, x), Poly(-1 + x, x)]])
    >>> v1 = PolyMatrix([[1, 0], [-1, 0]])
    >>> pm1*v1
    Matrix([
    [    Poly(x**2 + x, x, domain='ZZ'), Poly(0, x, domain='ZZ')],
    [Poly(x**3 - x + 1, x, domain='ZZ'), Poly(0, x, domain='ZZ')]])

    >>> pm1.ring
    ZZ[x]

    >>> v1*pm1
    Matrix([
    [ Poly(x**2, x, domain='ZZ'), Poly(-x, x, domain='ZZ')],
    [Poly(-x**2, x, domain='ZZ'),  Poly(x, x, domain='ZZ')]])

    >>> pm2 = PolyMatrix([[Poly(x**2, x, domain='QQ'), Poly(0, x, domain='QQ'), Poly(1, x, domain='QQ'), \
            Poly(x**3, x, domain='QQ'), Poly(0, x, domain='QQ'), Poly(-x**3, x, domain='QQ')]])
    >>> v2 = PolyMatrix([1, 0, 0, 0, 0, 0], ring=ZZ)
    >>> v2.ring
    ZZ
    >>> pm2*v2
    Matrix([[Poly(x**2, x, domain='QQ')]])

    """
    _class_priority = 10
    # we don't want to sympify the elements of PolyMatrix
    _sympify = staticmethod(lambda x: x)

    def __init__(self, *args, **kwargs):
        # if any non-Poly element is given as input then
        # 'ring' defaults 'EX'
        ring = kwargs.get('ring', EX)
        if all(isinstance(p, Poly) for p in self._mat) and self._mat:
            domain = tuple([p.domain[p.gens] for p in self._mat])
            ring = domain[0]
            for i in range(1, len(domain)):
                ring = ring.unify(domain[i])
        self.ring = ring

    def _eval_matrix_mul(self, other):
        self_cols = self.cols
        other_rows, other_cols = other.rows, other.cols
        other_len = other_rows*other_cols
        new_mat_rows = self.rows
        new_mat_cols = other.cols

        new_mat = [0]*new_mat_rows*new_mat_cols

        if self.cols != 0 and other.rows != 0:
            mat = self._mat
            other_mat = getattr(other, "_mat", None)
            if other_mat is None:
                other_mat = other._flat()
            for i in range(len(new_mat)):
                row, col = i // new_mat_cols, i % new_mat_cols
                row_indices = range(self_cols*row, self_cols*(row+1))
                col_indices = range(col, other_len, other_cols)
                vec = (mat[a]*other_mat[b] for a,b in zip(row_indices, col_indices))
                # 'Add' shouldn't be used here
                new_mat[i] = sum(vec)

        return self.__class__(new_mat_rows, new_mat_cols, new_mat, copy=False)

    def _eval_scalar_mul(self, other):
        mat = [Poly(a.as_expr()*other, *a.gens) if isinstance(a, Poly) else a*other for a in self._mat]
        return self.__class__(self.rows, self.cols, mat, copy=False)

    def _eval_scalar_rmul(self, other):
        mat = [Poly(other*a.as_expr(), *a.gens) if isinstance(a, Poly) else other*a for a in self._mat]
        return self.__class__(self.rows, self.cols, mat, copy=False)


MutablePolyMatrix = PolyMatrix = MutablePolyDenseMatrix


class DomainMatrix:

    def __init__(self, rows, shape, domain):
        self.shape = shape
        self.rows = [[item for item in row] for row in rows]
        self.domain = domain

    @classmethod
    def from_list_sympy(cls, rows):
        nrows = len(rows)
        ncols = len(rows[0])
        assert len(rows) == nrows
        assert all(len(row) == ncols for row in rows)

        items_sympy = [_sympify(item) for row in rows for item in row]

        domain, items_domain = cls.get_domain(items_sympy, field=True, extension=True)

        domain_rows = [[items_domain[ncols*r + c] for c in range(ncols)] for r in range(nrows)]

        return DomainMatrix(domain_rows, (nrows, ncols), domain)


    @classmethod
    def from_list_sympy_2(cls, nrows, ncols, rows):
        assert len(rows) == nrows
        assert all(len(row) == ncols for row in rows)

        items_sympy = [_sympify(item) for row in rows for item in row]

        domain, items_domain = cls.get_domain(items_sympy)

        domain_rows = [[items_domain[ncols*r + c] for c in range(ncols)] for r in range(nrows)]

        return DomainMatrix(domain_rows, (nrows, ncols), domain)

    @classmethod
    def get_domain(cls, items_sympy, **kwargs):
        K, items_K = construct_domain(items_sympy, **kwargs)
        return K, items_K

    def convert_to(self, K):
        Kold = self.domain
        new_rows = [[K.convert_from(e, Kold) for e in row] for row in self.rows]
        return DomainMatrix(new_rows, self.shape, K)

    def to_field(self):
        K = self.domain.get_field()
        return self.convert_to(K)

    def unify(self, other):
        K1 = self.domain
        K2 = other.domain
        if K1 == K2:
            return self, other
        K = K1.unify(K2)
        if K1 != K:
            self = self.convert_to(K)
        if K2 != K:
            other = other.convert_to(K)
        return self, other

    def to_Matrix(self):
        rows_sympy = [[self.domain.to_sympy(e) for e in row] for row in self.rows]
        return MutableDenseMatrix(rows_sympy)

    def __repr__(self):
        rows_str = ['[%s]' % (', '.join(map(str, row))) for row in self.rows]
        rowstr = '[%s]' % ', '.join(rows_str)
        return 'DomainMatrix(%s, %r, %r)' % (rowstr, self.shape, self.domain)

    def __add__(A, B):
        if not isinstance(B, DomainMatrix):
            return NotImplemented
        if A.shape != B.shape:
            raise ShapeError("shape")
        if A.domain != B.domain:
            raise ValueError("domain")
        rows = [[a+b for a, b in zip(row1, row2)]
                    for row1, row2 in zip(A.rows, B.rows)]
        return type(A)(rows, A.shape, A.domain)

    def __mul__(A, B):
        """A * B"""
        if not isinstance(B, DomainMatrix):
            return NotImplemented
        rows, shape = matrix_mul(A.rows, A.shape, B.rows, B.shape)
        domain = A.domain.unify(B.domain)
        return type(A)(rows, shape, domain)

    def __pow__(A, n):
        """A ** n"""
        if not isinstance(n, int):
            return NotImplemented
        if n < 0:
            raise NotImplementedError('Negative powers')
        elif n == 0:
            m, n = A.shape
            rows = [[A.domain.zero] * m for _ in range(m)]
            for i in range(m):
                rows[i][i] = A.domain.one
            return type(A)(rows, A.shape, A.domain)
        elif n == 1:
            return A
        elif n % 2 == 1:
            return A * A**(n - 1)
        else:
            sqrtAn = A ** (n // 2)
            return sqrtAn * sqrtAn

    def rref(self):
        if not self.domain.is_Field:
            raise ValueError('Not a field')
        rref_rows, pivots = rref(self.rows, self.shape)
        rref_matrix = type(self)(rref_rows, self.shape, self.domain)
        pivots = tuple(pivots)
        return rref_matrix, pivots

    def inv(self):
        if not self.domain.is_Field:
            raise ValueError('Not a field')
        m, n = self.shape
        if m != n:
            raise NonSquareMatrixError
        dom = self.domain
        rows = self.rows
        eye = [[dom.one if i==j else dom.zero for j in range(n)] for i in range(n)]
        rows = [row + eyerow for row, eyerow in zip(rows, eye)]
        Aaug = DomainMatrix(rows, (n, 2*n), dom)
        Aaug_rref, pivots = Aaug.rref()
        if pivots != tuple(range(n)):
            raise NonInvertibleMatrixError('Matrix det == 0; not invertible.')
        Ainv_rows = [row[n:] for row in Aaug_rref.rows]
        Ainv = DomainMatrix(Ainv_rows, (n, n), dom)
        return Ainv

    def det(self):
        m, n = self.shape
        if m != n:
            raise NonSquareMatrixError
        # Fraction-free Gaussian elimination
        # https://www.math.usm.edu/perry/Research/Thesis_DRL.pdf
        a = [row[:] for row in self.rows]
        dom = self.domain
        is_field = dom.is_Field
        # uf keeps track of the effect of row swaps and multiplies
        uf = dom.one
        for j in range(n-1):
            # if zero on the diagonal need to swap
            if not a[j][j]:
                for l in range(j+1, n):
                    if a[l][j]:
                        a[j], a[l] = a[l], a[j]
                        uf = -uf
                        break
                else:
                    # unable to swap: det = 0
                    return dom.zero
            for i in range(j+1, n):
                if a[i][j]:
                    if not is_field:
                        d = dom.gcd(a[j][j], a[i][j])
                        b = a[j][j] // d
                        c = a[i][j] // d
                    else:
                        b = a[j][j]
                        c = a[i][j]
                    # account for multiplying row i by b
                    uf = b * uf
                    for k in range(j+1, n):
                        a[i][k] = b*a[i][k] - c*a[j][k]
        # triangular det is product of diagonal
        prod = dom.one
        for i in range(n):
            prod = prod * a[i][i]
        # incorporate swaps and multiplies
        if not is_field:
            D = prod // uf
        else:
            D = prod / uf
        return D

    def __eq__(A, B):
        """A == B"""
        if not isinstance(B, DomainMatrix):
            return NotImplemented
        return A.rows == B.rows


def matrix_mul(items1, shape1, items2, shape2):
    m, n1 = shape1
    n2, o = shape2
    if n1 != n2:
        raise NonSquareMatrixError
    n = n1
    shape3 = (m, o)
    items3 = [[None] * o for _ in range(m)]
    for i in range(m):
        for j in range(o):
            items3[i][j] = sum(items1[i][k] * items2[k][j] for k in range(n))
    return items3, shape3


def rref(rows, shape):
    nrows, ncols = shape
    rows = [[item for item in row] for row in rows]
    pivots = []
    ri = 0
    for ci in range(ncols):
        for rj in range(ri, nrows):
            if rows[rj][ci]:
                # Row swap for pivot
                if rj != ri:
                    rows[rj], rows[ri] = rows[ri], rows[rj]
                # Record pivot
                pivots.append(ci)
                break
        else:
            # No pivot
            continue
        # Normalise row
        pivoti = rows[ri][ci]
        for ck in range(ci, ncols):
            rows[ri][ck] = rows[ri][ck] / pivoti
        # Eliminate above and below from col to the right
        for rk in range(nrows):
            pivotk = rows[rk][ci]
            if rk != ri and pivotk:
                for ck in range(ci, ncols):
                    rows[rk][ck] = rows[rk][ck] - pivotk * rows[ri][ck]
        ri += 1
    return rows, pivots


def lu_decomp(M):

    dom = M.domain

    if not dom.is_Field:
        raise ValueError("Domain")

    rows = M.rows
    N, O = M.shape

    lu = [row[:] for row in M.rows]
    swaps = []

    for n in range(min(N, O)):
        if not lu[n][n]:
            for m in range(n+1, N):
                if lu[m][n]:
                    swaps.append((n, m))
                    lu[n], lu[m] = lu[m], lu[n]
                    break
            else:
                break
        for i in range(n+1, N):
            l_in = lu[i][n] / lu[n][n]
            lu[i][n] = l_in
            for j in range(n+1, O):
                lu[i][j] -= l_in * lu[n][j]

    L = [[dom.one] + [dom.zero] * (N-1)]
    U = [lu[0]]
    for i in range(1, N):
        L.append(lu[i][:i] + [dom.one] + [dom.zero] * (N-i))
        U.append([dom.zero] * i + lu[i][i:])
    L = DomainMatrix(L, (N, N), dom)
    U = DomainMatrix(U, (N, O), dom)
    return L, U, swaps
