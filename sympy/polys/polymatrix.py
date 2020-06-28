from __future__ import print_function

from sympy.core.singleton import S
from sympy.core.symbol import Symbol
from sympy.core.sympify import _sympify
from sympy.matrices.dense import MutableDenseMatrix
from sympy.polys.fields import sfield
from sympy.polys.polytools import Poly
from sympy.polys.domains import EX, PolynomialRing

import sympy.polys.polyoptions as polyoptions


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

    @classmethod
    def _new(cls, *args, ring=None, **kwargs):
        if kwargs.get('copy', True) is False:
            if len(args) != 3:
                raise TypeError(
                    "copy=False requires a matrix be initialized as "
                    "rows, cols, [list].")
            rows, cols, flat_list = args
        else:
            rows, cols, flat_list = cls._handle_creation_inputs(*args, **kwargs)
            flat_list = list(flat_list) # create a shallow copy

        self = object.__new__(cls)
        self.rows = rows
        self.cols = cols
        self._mat = flat_list

        if ring:
            ring = polyoptions.Domain.preprocess(ring)
        # XXX PolyMatrix must always have Poly instances or sympifiable
        # objects that can be comprehended as Poly.
        if all(isinstance(p, Poly) for p in self._mat) and self._mat:
            domains = set([p.domain[p.gens] for p in self._mat])
            if not ring:
                ring = domains.pop()
            for dom in domains:
                ring = ring.unify(dom)
        # XXX PolyMatrix should not have EX domain
        if not ring:
            ring = EX
        self.ring = ring

        if isinstance(ring, PolynomialRing):
            self.zero = Poly(0, *ring.symbols, domain=ring.domain)
            self.one = Poly(1, *ring.symbols, domain=ring.domain)
        # XXX PolyMatrix should not have EX domain
        elif ring == EX:
            self.zero = S.Zero
            self.one = S.One
        else:
            self.zero = ring.zero
            self.one = ring.one
        return self

    def _eval_Abs(self):
        raise NotImplementedError

    def _eval_Mod(self, other):
        raise NotImplementedError

    def _eval_adjoint(self):
        raise NotImplementedError

    def _eval_conjugate(self):
        raise NotImplementedError

    @classmethod
    def _eval_companion(cls, size, poly):
        coeffs = poly.all_coeffs()

        ring = poly.domain[poly.gens]
        zero = Poly(0, *poly.gens, domain=poly.domain)
        one = Poly(1, *poly.gens, domain=poly.domain)

        def entry(i, j):
            if j == size - 1:
                return Poly(-coeffs[-1 - i], *poly.gens, domain=poly.domain)
            elif i == j + 1:
                return one
            return zero

        return cls._new(size, size, entry, ring=ring)

    def _eval_as_real_imag(self):
        raise NotImplementedError

    def _eval_is_anti_symmetric(self, simpfunc):
        for i in range(self.rows):
            if self[i, i] != self.zero:
                return False
            for j in range(i+1, self.cols):
                if self[i, j] + self[j, i] != self.zero:
                    return False
        return True

    def _eval_is_matrix_hermitian(self, simpfunc):
        raise NotImplementedError

    def _eval_is_symmetric(self, simpfunc):
        for i in range(self.rows):
            for j in range(i+1, self.cols):
                if self[i, j] != self[j, i]:
                    return False
        return True

    @classmethod
    def _eval_jordan_block(cls, rows, cols, eigenvalue, band='upper'):
        if isinstance(eigenvalue, Poly):
            ring = eigenvalue.domain[eigenvalue.gens]
            zero = Poly(0, *eigenvalue.gens, domain=eigenvalue.domain)
            one = Poly(1, *eigenvalue.gens, domain=eigenvalue.domain)
        else:
            ring = EX
            zero, one = S.Zero, S.One

        if band == 'lower':
            def entry(i, j):
                if i == j:
                    return eigenvalue
                elif j + 1 == i:
                    return one
                return zero
        else:
            def entry(i, j):
                if i == j:
                    return eigenvalue
                elif i + 1 == j:
                    return one
                return zero

        return cls._new(rows, cols, entry, ring=ring)

    def _eval_matrix_mul(self, other):
        if not isinstance(other, PolyMatrix):
            other = PolyMatrix(other)
        ring = self.ring.unify(other.ring)
        self_cols = self.cols
        other_rows, other_cols = other.rows, other.cols
        other_len = other_rows*other_cols
        new_mat_rows = self.rows
        new_mat_cols = other.cols
        new_mat = [ring.zero]*new_mat_rows*new_mat_cols
        if self.cols and other.rows:
            mat = self._mat
            other_mat = other._mat
            for i in range(len(new_mat)):
                row, col = i // new_mat_cols, i % new_mat_cols
                row_indices = range(self_cols*row, self_cols*(row+1))
                col_indices = range(col, other_len, other_cols)
                vec = (mat[a] * other_mat[b] for a, b in
                       zip(row_indices, col_indices))
                # 'Add' shouldn't be used here
                new_mat[i] = sum(vec)

        return self.__class__(
            new_mat_rows, new_mat_cols, new_mat, ring=ring, copy=False)

    def _eval_scalar_mul(self, other):
        ring = self.ring
        if isinstance(other, Poly):
            ring = ring.unify(other.domain[other.gens])
        elif other not in self.ring:
            ring = ring.unify(EX)

        mat = [Poly(a.as_expr()*other, *a.gens) if isinstance(a, Poly) else a*other for a in self._mat]
        return self.__class__(self.rows, self.cols, mat, ring=ring, copy=False)

    def _eval_scalar_rmul(self, other):
        ring = self.ring
        if isinstance(other, Poly):
            ring = ring.unify(other.domain[other.gens])
        elif other not in self.ring:
            ring = ring.unify(EX)

        mat = [Poly(other*a.as_expr(), *a.gens) if isinstance(a, Poly) else other*a for a in self._mat]
        return self.__class__(self.rows, self.cols, mat, ring=ring, copy=False)

    def _eval_simplify(self, *args, **kwargs):
        return self

    def _eval_trigsimp(self, *args, **kwargs):
        return self

    def __pow__(self, exp):
        return self.pow(exp, method='multiply')

    def _eval_det_bareiss(self, *args, **kwargs):
        from sympy.matrices.determinant import _poly_det_bareiss
        return _poly_det_bareiss(self)

    def _eval_det_berkowitz(self, *args, **kwargs):
        from sympy.matrices.determinant import _poly_det_berkowitz
        return _poly_det_berkowitz(self)

    def _eval_det_lu(self, *args, **kwargs):
        raise NotImplementedError

    def _eval_pow_by_cayley(self, exp):
        raise NotImplementedError

    def _eval_pow_by_recursion_dotprodsimp(self, exp):
        raise NotImplementedError

    def _eval_pow_by_jordan_blocks(self, exp):
        raise NotImplementedError

    @property
    def D(self):
        raise NotImplementedError

    def LDLdecomposition(self, *args, **kwargs):
        raise NotImplementedError

    def LDLsolve(self, *args, **kwargs):
        raise NotImplementedError

    def LUdecomposition(self, *args, **kwargs):
        raise NotImplementedError

    def LUdecomposition_Simple(self, *args, **kwargs):
        raise NotImplementedError

    def LUdecompositionFF(self, *args, **kwargs):
        raise NotImplementedError

    def LUsolve(self, *args, **kwargs):
        raise NotImplementedError

    def QRdecomposition(self, *args, **kwargs):
        raise NotImplementedError

    def QRsolve(self, *args, **kwargs):
        raise NotImplementedError

    def analytic_func(self, *args, **kwargs):
        raise NotImplementedError

    def bidiagonalize(self, *args, **kwargs):
        raise NotImplementedError

    def bidiagonal_decomposition(self, *args, **kwargs):
        raise NotImplementedError

    def charpoly(self, x='lambda', **kwargs):
        from sympy.matrices.determinant import _poly_charpoly
        return _poly_charpoly(self, x)

    def cholesky(self, *args, **kwargs):
        raise NotImplementedError

    def cholesky_solve(self, *args, **kwargs):
        raise NotImplementedError

    def condition_number(self):
        raise NotImplementedError

    def diagonal_solve(self, *args, **kwargs):
        raise NotImplementedError

    def diff(self, *args, **kwargs):
        if not all(isinstance(x, Poly) for x in self._mat):
            raise NotImplementedError
        if not all(isinstance(x, Symbol) for x in args):
            raise NotImplementedError
        return self.applyfunc(lambda p: p.diff(*args, **kwargs))

    def echelon_form(self, with_pivots=False, **kwargs):
        from sympy.matrices.reductions import _poly_echelon_form
        return _poly_echelon_form(self, with_pivots=with_pivots)

    def eigenvals(self, *args, **kwargs):
        raise NotImplementedError

    def eigenvects(self, *args, **kwargs):
        raise NotImplementedError

    def evalf(self, *args, **kwargs):
        raise NotImplementedError

    def exp(self, *args, **kwargs):
        raise NotImplementedError

    def expand(self, *args, **kwargs):
        raise NotImplementedError

    def gauss_jordan_solve(self, *args, **kwargs):
        raise NotImplementedError

    def integrate(self, *args, **kwargs):
        if not all(isinstance(x, Poly) for x in self._mat):
            raise NotImplementedError
        return self.applyfunc(lambda p: p.integrate(*args, **kwargs))

    def inv(self, *args, **kwargs):
        raise NotImplementedError

    def inv_mod(self, *args, **kwargs):
        raise NotImplementedError

    def inverse_ADJ(self, *args, **kwargs):
        raise NotImplementedError

    def inverse_BLOCK(self, *args, **kwargs):
        raise NotImplementedError

    def inverse_CH(self, *args, **kwargs):
        raise NotImplementedError

    def inverse_GE(self, *args, **kwargs):
        raise NotImplementedError

    def inverse_LDL(self, *args, **kwargs):
        raise NotImplementedError

    def inverse_LU(self, *args, **kwargs):
        raise NotImplementedError

    def inverse_QR(self, *args, **kwargs):
        raise NotImplementedError

    @property
    def is_indefinite(self):
        raise NotImplementedError

    @property
    def is_negative_definite(self):
        raise NotImplementedError

    @property
    def is_negative_semidefinite(self):
        raise NotImplementedError

    @property
    def is_nilpotent(self):
        raise NotImplementedError

    @property
    def is_positive_definite(self):
        raise NotImplementedError

    @property
    def is_positive_semidefinite(self):
        raise NotImplementedError

    @property
    def is_strongly_diagonally_dominant(self):
        raise NotImplementedError

    @property
    def is_weakly_diagonally_dominant(self):
        raise NotImplementedError

    def jacobian(self, other):
        raise NotImplementedError

    def limit(self, *args, **kwargs):
        raise NotImplementedError

    def log(self, *args, **kwargs):
        raise NotImplementedError

    def lower_triangular_solve(self, *args, **kwargs):
        raise NotImplementedError

    def normalized(self, *args, **kwargs):
        raise NotImplementedError

    @classmethod
    def orthogonalize(self, *args, **kwargs):
        raise NotImplementedError

    def pinv(self, *args, **kwargs):
        raise NotImplementedError

    def pinv_solve(self, *args, **kwargs):
        raise NotImplementedError

    def rref(self, pivots=True, normalize_last=True, **kwargs):
        if self.ring == EX:
            return super().rref(
                pivots=pivots, normalize_last=normalize_last, **kwargs)

        M = self._to_domain_matrix()
        rref, pivots = M.rref()
        rref = self._from_domain_matrix(rref)
        if pivots:
            return rref, pivots
        return rref

    def simplify(self, *args, **kwargs):
        return

    def solve(self, *args, **kwargs):
        raise NotImplementedError

    def solve_least_squares(self, *args, **kwargs):
        raise NotImplementedError

    def upper_triangular_solve(self, *args, **kwargs):
        raise NotImplementedError

    def _to_field(self, dummy='X'):
        """Helper function to convert all polynomials of the matrix to
        fraction field elements such that the division is possible."""
        from sympy.core.symbol import uniquely_named_symbol

        ring = self.ring
        if ring.is_Field:
            return self

        if not all(isinstance(x, Poly) for x in self._mat):
            return self

        if all(x.degree() <= 0 for x in self._mat):
            ring, gens = ring.domain, ring.symbols
            ring = ring.get_field()
            new_mat = []
            for x in self._mat:
                x = Poly(x.as_expr(), *gens, domain=ring)
                new_mat.append(x)
            ring = ring[gens]
            return self._new(self.rows, self.cols, new_mat, ring=ring)
        elif all(x.degree() > 0 for x in self._mat):
            ring = ring.get_field()
            dummy = uniquely_named_symbol(dummy, self)
            new_mat = []
            for x in self._mat:
                x = Poly(x.as_expr(), dummy, domain=ring)
                new_mat.append(x)
            ring = ring[dummy]
            return self._new(self.rows, self.cols, new_mat, ring=ring)
        else:
            # XXX FractionField has some issues with
            # ZZ[a, b][x].get_field()(a + b)
            raise NotImplementedError

    def _to_domain_matrix(self):
        rows = []

        ring = self.ring
        if all(isinstance(x, Poly) for x in self._mat) and \
            all(x.degree() <= 0 for x in self._mat):
            ring = ring.domain

        for i in range(self.rows):
            row = []
            for j in range(self.cols):
                if isinstance(self[i, j], Poly):
                    # XXX Poly to Domain conversion should be more
                    # straightforward than Expr to Domain.
                    x = ring(self[i, j].as_expr())
                else:
                    x = ring(self[i, j])
                row.append(x)
            rows.append(row)
        return DomainMatrix(rows, self.shape, ring)

    def _from_domain_matrix(self, M):
        flat = []
        ring = self.ring
        for row in M.rows:
            for x in row:
                x = M.domain.to_sympy(x)
                if isinstance(ring, PolynomialRing):
                    x = Poly(x, *ring.symbols, domain=ring.domain)
                flat.append(x)

        return self._new(self.rows, self.cols, flat, ring=self.ring)


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

        domain, items_domain = cls.get_domain(items_sympy)

        domain_rows = [[items_domain[ncols*r + c] for c in range(ncols)] for r in range(nrows)]

        return DomainMatrix(domain_rows, (nrows, ncols), domain)

    @classmethod
    def get_domain(cls, items_sympy):
        K, items_K = sfield(items_sympy, field=True, extension=True)

        if K.gens:
            domain = K.to_domain()
        else:
            domain = K.domain

            def convert(item):
                if not item:
                    return domain.zero
                else:
                    return item.numer[()] / item.denom[()]

            items_K = [convert(item) for item in items_K]

        return domain, items_K

    def to_Matrix(self):
        rows_sympy = [[self.domain.to_sympy(e) for e in row] for row in self.rows]
        return MutableDenseMatrix(rows_sympy)

    def __repr__(self):
        rows_str = ['[%s]' % (', '.join(map(str, row))) for row in self.rows]
        rowstr = '[%s]' % ', '.join(rows_str)
        return 'DomainMatrix(%s, %r, %r)' % (rowstr, self.shape, self.domain)

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
        rref_rows, pivots = rref(self.rows, self.shape)
        rref_matrix = type(self)(rref_rows, self.shape, self.domain)
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
            if rk != ri and pivotk:
                for ck in range(ci, ncols):
                    rows[rk][ck] = rows[rk][ck] - pivotk * rows[ri][ck]
        ri += 1
    return rows, pivots
