from __future__ import print_function

from sympy.core.add import Add
from sympy.core.numbers import I
from sympy.core.power import Pow
from sympy.core.sympify import _sympify

from sympy.matrices.dense import MutableDenseMatrix

from sympy.polys.polytools import Poly
from sympy.polys.domains import EX, QQ, ZZ, AlgebraicField, PolynomialRing


class DomainMatrixDomainError(NotImplementedError):
    """Raised when DomainMatrix is unable to construct a domain"""
    pass


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
            other_mat = other._mat
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

        rows_sympy = [[_sympify(item) for item in row] for row in rows]
        domain = cls.get_domain([item for row in rows_sympy for item in row])
        domain_rows = [[domain.from_sympy(item) for item in row] for row in rows_sympy]

        return DomainMatrix(domain_rows, (nrows, ncols), domain)

    def to_Matrix(self):
        rows_sympy = [[self.domain.to_sympy(e) for e in row] for row in self.rows]
        return MutableDenseMatrix(rows_sympy)

    @classmethod
    def get_domain(cls, items):
        ext = set()
        syms = set()
        dom = ZZ
        items = items[:]
        for item in items:
            if item.has(I):
                items.extend(Poly(item, I).coeffs())
                ext.add(I)
                dom = dom.unify(QQ)
                continue
            elif item.is_polynomial():
                itemsyms = item.free_symbols
                if itemsyms:
                    items.extend(Poly(item, *itemsyms).coeffs())
                    syms |= itemsyms
                    continue
            if item.has(Pow):
                for pow in item.atoms(Pow):
                    b, e = pow.as_base_exp()
                    if b.is_Number and e.is_Rational and e.is_positive:
                        items.extend(Poly(item, pow).coeffs())
                        ext.add(pow)
                        dom = dom.unify(QQ)
                        break
                    else:
                        raise DomainMatrixDomainError
                continue
            if item.is_Integer:
                pass
            elif item.is_Rational:
                dom = dom.unify(QQ)
            else:
                raise DomainMatrixDomainError
        if ext:
            dom = AlgebraicField(dom, *ext)
        if syms:
            dom = PolynomialRing(dom, list(syms))
        return dom

    def __repr__(self):
        return str(self.to_Matrix())

    def __mul__(A, B):
        """A * B"""
        if not isinstance(B, DomainMatrix):
            return NotImplemented
        rows, shape = matrix_mul(A.rows, A.shape, B.rows, B.shape)
        domain = A.domain.unify(B.domain)
        return type(A)(rows, shape, domain)

    def __pow__(A, n):
        """A ** n"""
        if n == 1:
            return A
        elif n % 2 == 1:
            return A * A**(n - 1)
        else:
            sqrtAn = A ** (n // 2)
            return sqrtAn * sqrtAn

    def rref(self):
        rref_domain = self.domain.get_field()
        convert = lambda item: rref_domain.convert(item)
        converted_rows = [[convert(item) for item in row] for row in self.rows]
        rref_rows, pivots = rref(converted_rows, self.shape)
        rref_matrix = type(self)(rref_rows, self.shape, rref_domain)
        pivots = tuple(pivots)
        return rref_matrix, pivots

    def __eq__(A, B):
        """A == B"""
        if not isinstance(B, DomainMatrix):
            return NotImplemented
        return A.rows == B.rows


def matrix_mul(items1, shape1, items2, shape2):
    m, n1 = shape1
    n2, o = shape2
    assert n1 == n2
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
            if rk != ri:
                for ck in range(ci, ncols):
                    rows[rk][ck] = rows[rk][ck] - pivotk * rows[ri][ck]
        ri += 1
    return rows, pivots


def linsolve_domain(system, symbols):
    # Possibly raises DomainMatrixDomainError
    Aaugdm = DomainMatrix.from_list_sympy(system.tolist())

    Aaugdm_rref, pivots = Aaugdm.rref()
    Aaug = Aaugdm_rref.to_Matrix()

    if Aaug.is_zero_matrix:
        return {sym:sym for sym in symbols}

    rows, cols = Aaug.shape
    while all(Aaug[rows-1, c] == 0 for c in range(cols)):
        rows -= 1
    if all(Aaug[rows-1, c] == 0 for c in range(cols-1)):
        return None
    else:
        sol = {}
        for col, sym in zip(reversed(range(cols-1)), reversed(symbols)):
            if col not in pivots:
                sol[sym] = sym
            else:
                terms = [Aaug[rows-1, cols-1]] # rhs
                for c in range(col+1, cols-1):
                    terms.append( - Aaug[rows-1, c] * sol[symbols[c]])
                sol[sym] = Add(*terms)
                rows -= 1
        return sol
