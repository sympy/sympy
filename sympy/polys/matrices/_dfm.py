#
# sympy.polys.matrices.dfm
#
# This modules defines the DFM class which is a wrapper for dense flint
# matrices as found in python-flint.
#
# As of python-flint 0.4.1 matrices over the following domains can be supported
# by python-flint:
#
#   ZZ: flint.fmpz_mat
#   QQ: flint.fmpq_mat
#   GF(p): flint.nmod_mat (p prime and p < ~2**62)
#
# The underlying flint library has many more domains, but these are not yet
# supported by python-flint.
#
# The DFM class is a wrapper for the flint matrices and provides a common
# interface for all supported domains that is interchangeable with the DDM
# and SDM classes so that DomainMatrix can be used with any as its internal
# matrix representation.
#

from sympy.external.importtools import import_module
from sympy.utilities.decorator import doctest_depends_on

from sympy.polys.domains import ZZ, QQ
from .exceptions import DMBadInputError

flint = import_module('flint')


__all__ = ['DFM']


@doctest_depends_on(ground_types=['flint'])
class DFM:
    """
    Dense FLINT matrix. This class is a wrapper for matrices from python-flint.

    >>> from sympy.polys.domains import ZZ
    >>> from sympy.polys.matrices.dfm import DFM
    >>> dfm = DFM([[ZZ(1), ZZ(2)], [ZZ(3), ZZ(4)]], (2, 2), ZZ)
    >>> dfm
    [[1, 2], [3, 4]]
    >>> dfm.rep
    [1, 2]
    [3, 4]
    >>> type(dfm.rep)
    <class 'flint._flint.fmpz_mat'>

    Usually, the DFM class is not instantiated directly, but is created as the
    internal representation of :class:`DomainMatrix`. When `SYMPY_GROUND_TYPES`
    is set to `flint` and `python-flint` is installed, the DFM class is used
    automatically as the internal representation of :class:`DomainMatrix`in
    dense format if the domain is supported by python-flint.

    >>> from sympy.polys.matrices.domainmatrix import DM
    >>> dM = DM([[1, 2], [3, 4]], ZZ)
    >>> dM.rep
    [[1, 2], [3, 4]]

    A :class:`DomainMatrix` can be converted to a DFM by calling the
    :meth:`to_dfm` method::

    >>> dM.to_dfm()
    [[1, 2], [3, 4]]

    """

    fmt = 'dense'
    is_DFM = True
    is_DDM = False

    def __new__(cls, rowslist, shape, domain):
        flint_mat = cls._get_flint_func(domain)

        if 0 not in shape:
            rep = flint_mat(rowslist)
        else:
            rep = flint_mat(*shape)

        return cls._new(rep, shape, domain)

    @classmethod
    def _new(cls, rep, shape, domain):
        cls._check(rep, shape, domain)
        obj = object.__new__(cls)
        obj.rep = rep
        obj.shape = obj.rows, obj.cols = shape
        obj.domain = domain
        return obj

    @classmethod
    def _check(cls, rep, shape, domain):
        repshape = (rep.nrows(), rep.ncols())
        if repshape != shape:
            raise DMBadInputError("Shape of rep does not match shape of DFM")
        if domain == ZZ and not isinstance(rep, flint.fmpz_mat):
            raise RuntimeError("Rep is not a flint.fmpz_mat")
        elif domain == QQ and not isinstance(rep, flint.fmpq_mat):
            raise RuntimeError("Rep is not a flint.fmpq_mat")
        elif domain not in (ZZ, QQ):
            raise NotImplementedError("Only ZZ and QQ are supported by DFM")

    @classmethod
    def _supports_domain(cls, domain):
        return domain in (ZZ, QQ)

    @classmethod
    def _get_flint_func(cls, domain):
        if domain == ZZ:
            return flint.fmpz_mat
        elif domain == QQ:
            return flint.fmpq_mat
        else:
            raise NotImplementedError("Only ZZ and QQ are supported by DFM")

    def __str__(self):
        """Return ``str(self)``."""
        return str(self.to_ddm())

    def __repr__(self):
        """Return ``repr(self)``."""
        return f'DFM{repr(self.to_ddm())[3:]}'

    def __eq__(self, other):
        """Return ``self == other``."""
        if not isinstance(other, DFM):
            return NotImplemented
        # Compare domains first because we do *not* want matrices with
        # different domains to be equal but e.g. a flint fmpz_mat and fmpq_mat
        # with the same entries will compare equal.
        return self.domain == other.domain and self.rep == other.rep

    def copy(self):
        """Return a copy of self."""
        # XXX: fmpz_mat and fmpq_mat do not have a copy method
        return self.to_ddm().to_dfm()

    def to_ddm(self):
        """Convert to a DDM."""
        return DDM(self.rep.tolist(), self.shape, self.domain)

    def to_sdm(self):
        """Convert to a SDM."""
        return self.to_ddm().to_sdm()

    def to_dfm(self):
        """Return self."""
        return self

    def to_dfm_or_ddm(self):
        """Return self."""
        return self

    def to_list(self):
        """Convert to a nested list."""
        return self.rep.tolist()

    def to_list_flat(self):
        """Convert to a flat list."""
        return self.to_ddm().to_list_flat()

    def to_flat_nz(self):
        """Convert to a flat list of non-zeros."""
        return self.to_ddm().to_flat_nz()

    @classmethod
    def from_flat_nz(cls, elements, data, domain):
        """Inverse of :meth:`to_flat_nz`."""
        return DDM.from_flat_nz(elements, data, domain).to_dfm()

    def to_dok(self):
        """Convert to a DOK."""
        return self.to_ddm().to_dok()

    def convert_to(self, domain):
        """Convert to a new domain."""
        # XXX: fmpz_mat and fmpq_mat do not have conversion methods. You can do
        # fmpq_mat(fmpz_mat) but not fmpz_mat(fmpq_mat). The domain will be
        # checked in __new__.
        #
        # It is the responsibility of the caller to ensure that the conversion
        # is possible (i.e. DomainMatrix should have checked already).
        return self.to_ddm().convert_to(domain).to_dfm()

    def getitem(self, i, j):
        """Get the ``(i, j)``-th entry."""
        # XXX: flint matrices do not support negative indices
        # XXX: They also raise ValueError instead of IndexError
        m, n = self.shape
        if i < 0:
            i += m
        if j < 0:
            j += n
        try:
            return self.rep[i, j]
        except ValueError:
            raise IndexError(f"Invalid indices ({i}, {j}) for Matrix of shape {self.shape}")

    def setitem(self, i, j, value):
        """Set the ``(i, j)``-th entry."""
        # XXX: flint matrices do not support negative indices
        # XXX: They also raise ValueError instead of IndexError
        m, n = self.shape
        if i < 0:
            i += m
        if j < 0:
            j += n
        try:
            self.rep[i, j] = value
        except ValueError:
            raise IndexError(f"Invalid indices ({i}, {j}) for Matrix of shape {self.shape}")

    def extract(self, rowslist, colslist):
        """Extract a submatrix."""
        # XXX: flint matrices do not support fancy indexing
        return self.to_ddm().extract(rowslist, colslist).to_dfm()

    def extract_slice(self, rowslice, colslice):
        """Extract a submatrix."""
        # XXX: flint matrices do not support slicing
        return self.to_ddm().extract_slice(rowslice, colslice).to_dfm()

    def neg(self):
        """Negate a DFM matrix."""
        return self._new(-self.rep, self.shape, self.domain)

    def add(self, other):
        """Add two DFM matrices."""
        return self._new(self.rep + other.rep, self.shape, self.domain)

    def sub(self, other):
        """Subtract two DFM matrices."""
        return self._new(self.rep - other.rep, self.shape, self.domain)

    def mul(self, other):
        """Multiply a DFM matrix from the right by a scalar."""
        return self._new(self.rep * other, self.shape, self.domain)

    def rmul(self, other):
        """Multiply a DFM matrix from the left by a scalar."""
        return self._new(other * self.rep, self.shape, self.domain)

    def mul_elementwise(self, other):
        """Elementwise multiplication of two DFM matrices."""
        return self.to_ddm().mul_elementwise(other.to_ddm()).to_dfm()

    def matmul(self, other):
        """Multiply two DFM matrices."""
        shape = (self.rows, other.cols)
        return self._new(self.rep * other.rep, shape, self.domain)

    def __neg__(self):
        """Negate a DFM matrix."""
        return self.neg()

    @classmethod
    def zeros(cls, shape, domain):
        """Return a zero DFM matrix."""
        func = cls._get_flint_func(domain)
        return cls(func(*shape), shape, domain)

    @classmethod
    def eye(cls, n, domain):
        """Return the identity matrix of size n."""
        # XXX: DDM and SDM have inconsistent signatures for eye because SDM
        # assumes a shape argument while DDM expects an integer size.
        return DDM.eye(n, domain).to_dfm()

    def applyfunc(self, func, domain):
        """Apply a function to each entry of a DFM matrix."""
        return self.to_ddm().applyfunc(func, domain).to_dfm()

    def transpose(self):
        """Transpose a DFM matrix."""
        m, n = self.shape
        return self._new(self.rep.transpose(), (n, m), self.domain)

    def hstack(self, *others):
        """Horizontally stack matrices."""
        return self.to_ddm().hstack(*[o.to_ddm() for o in others]).to_dfm()

    def vstack(self, *others):
        """Vertically stack matrices."""
        return self.to_ddm().vstack(*[o.to_ddm() for o in others]).to_dfm()

    def diagonal(self):
        """Return the diagonal of a DFM matrix."""
        return self.to_ddm().diagonal()

    def is_upper(self):
        """Return ``True`` if the matrix is upper triangular."""
        return self.to_ddm().is_upper()

    def is_lower(self):
        """Return ``True`` if the matrix is lower triangular."""
        return self.to_ddm().is_lower()

    def is_diagonal(self):
        """Return ``True`` if the matrix is diagonal."""
        return self.to_ddm().is_diagonal()

    def is_zero_matrix(self):
        """Return ``True`` if the matrix is the zero matrix."""
        return self.to_ddm().is_zero_matrix()

    def scc(self):
        """Return the strongly connected components of the matrix."""
        return self.to_ddm().scc()

    def det(self):
        """Return the determinant of the matrix."""
        # XXX: Use the flint det method!!!
        return self.to_ddm().det()

    def inv(self):
        """Return the inverse of the matrix."""
        # XXX: Use the flint inv method!!!
        return self.to_ddm().inv().to_dfm()

    def charpoly(self):
        """Return the characteristic polynomial of the matrix."""
        # XXX: Use the flint charpoly method!!!
        return self.to_ddm().charpoly()

    def lu(self):
        """Return the LU decomposition of the matrix."""
        L, U, swaps = self.to_ddm().lu()
        return L.to_dfm(), U.to_dfm(), swaps

    def lu_solve(self, rhs):
        """Solve a linear system using LU decomposition."""
        return self.to_ddm().lu_solve(rhs.to_ddm()).to_dfm()

    def nullspace(self):
        """Return a basis for the nullspace of the matrix."""
        ddm, nonpivots = self.to_ddm().nullspace()
        return ddm.to_dfm(), nonpivots

    def nullspace_from_rref(self, pivots=None):
        """Return a basis for the nullspace of the matrix."""
        sdm, nonpivots = self.to_sdm().nullspace_from_rref(pivots=pivots)
        return sdm.to_dfm(), nonpivots

    def particular(self):
        """Return a particular solution to the system."""
        return self.to_ddm().particular().to_dfm()

    def lll(self, delta=QQ(3, 4)):
        """Return an LLL-reduced basis for the matrix."""
        # XXX: Use the flint lll method!!!
        return self.to_ddm().lll(delta=delta).to_dfm()

    def lll_transform(self, delta=QQ(3, 4)):
        """Return the LLL-reduced basis and transform matrix."""
        # XXX: Use the flint lll_transform method!!!
        y, T = self.to_ddm().lll_transform(delta=delta)
        return y.to_dfm(), T.to_dfm()


# Avoid circular imports
from sympy.polys.matrices.ddm import DDM
