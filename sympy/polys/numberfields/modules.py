r"""Modules in number fields.

The classes defined here allow us to work with finitely generated, free
modules, whose generators are algebraic numbers.

There is an abstract base class called ``Module``, which has two concrete
subclasses, ``PowerBasis`` and ``Submodule``.

Every module is defined by its basis, or set of generators:

* For a ``PowerBasis``, the generators are the first $n$ powers (starting with
  the zeroth) of an algebraic integer $\theta$ of degree $n$. The
  ``PowerBasis`` is constructed by passing the minimal polynomial of $\theta$.

* For a ``Submodule``, the generators are a set of $\mathbb{Q}$-linear
  combinations of the generators of another module. That other module is then
  the "parent" of the ``Submodule``. The coefficients of the
  $\mathbb{Q}$-linear combinations may be given by an integer matrix, and a
  positive integer denominator. Each column of the matrix defines a generator.

>>> from sympy.polys import Poly, cyclotomic_poly, ZZ
>>> from sympy.abc import x
>>> from sympy.polys.matrices import DomainMatrix, DM
>>> from sympy.polys.numberfields.modules import PowerBasis
>>> T = Poly(cyclotomic_poly(5, x))
>>> A = PowerBasis(T)
>>> print(A)
PowerBasis(x**4 + x**3 + x**2 + x + 1)
>>> B = A.submodule_from_matrix(2 * DomainMatrix.eye(4, ZZ), denom=3)
>>> print(B)
Submodule[[2, 0, 0, 0], [0, 2, 0, 0], [0, 0, 2, 0], [0, 0, 0, 2]]/3
>>> print(B.parent)
PowerBasis(x**4 + x**3 + x**2 + x + 1)

Thus, every module is either a ``PowerBasis``, or a ``Submodule``, some
ancestor of which is a ``PowerBasis``. (If ``S`` is a ``Submodule``, then its
ancestors are ``S.parent``, ``S.parent.parent``, and so on).

The ``ModuleElement`` class represents a linear combination of the generators
of any module. Critically, the coefficients of this linear combination are not
restricted to be integers, but may be any rational numbers. This is necessary
so that any and all algebraic integers be representable, starting from the
power basis in a primitive element $\theta$ for the number field in question.
For example, in a quadratic field $\mathbb{Q}(\sqrt{d})$ where
$d \equiv 1 \mod{4}$, a denominator of $2$ is needed.

A ``ModuleElement`` can be constructed from an integer column vector and a
denominator:

>>> U = Poly(x**2 - 5)
>>> M = PowerBasis(U)
>>> e = M(DM([[1], [1]], ZZ), denom=2)
>>> print(e)
[1, 1]/2
>>> print(e.module)
PowerBasis(x**2 - 5)

The ``PowerBasisElement`` class is a subclass of ``ModuleElement`` that
represents elements of a ``PowerBasis``, and adds functionality pertinent to
elements represented over powers of the primitive element $\theta$.


Arithmetic with module elements
===============================

While a ``ModuleElement`` represents a linear combination over the generators
of a particular module, recall that every module is either a ``PowerBasis``
or a descendant (along a chain of ``Submodule`` objects) thereof, so that in
fact every ``ModuleElement`` represents an algebraic number in some field
$\mathbb{Q}(\theta)$, where $\theta$ is the defining element of some
``PowerBasis``. It thus makes sense to talk about the number field to which a
given ``ModuleElement`` belongs.

This means that any two ``ModuleElement`` instances can be added, subtracted,
multiplied, or divided, provided they belong to the same number field.
Similarly, since $\mathbb{Q}$ is a subfield of every number field, any
``ModuleElement`` may be added, multiplied, etc. by any rational number.

>>> from sympy import QQ
>>> from sympy.polys.numberfields.modules import to_col
>>> T = Poly(cyclotomic_poly(5))
>>> A = PowerBasis(T)
>>> C = A.submodule_from_matrix(3 * DomainMatrix.eye(4, ZZ))
>>> e = A(to_col([1, 2, 3, 4]), denom=6)
>>> f = A(to_col([5, 6, 7, 8]), denom=10)
>>> g = C(to_col([1, 1, 1, 1]), denom=2)
>>> print(e + f)
[10, 14, 18, 22]/15
>>> print(e - f)
[-5, -4, -3, -2]/15
>>> print(e + g)
[10, 11, 12, 13]/6
>>> print(e + QQ(7, 10))
[26, 10, 15, 20]/30
>>> print(g + QQ(7, 10))
[22, 15, 15, 15]/10

However, care must be taken with arithmetic operations on ``ModuleElemnt``,
because the module $C$ to which the result will belong will be the nearest
common ancestor (NCA) of the modules $A$, $B$ to which the two operands belong,
and $C$ may be different from either or both of $A$ and $B$.

>>> A = PowerBasis(T)
>>> B = A.submodule_from_matrix(2 * DomainMatrix.eye(4, ZZ))
>>> C = A.submodule_from_matrix(3 * DomainMatrix.eye(4, ZZ))
>>> print((B(0) * C(0)).module == A)
True

Before the arithmetic operation is performed, copies of the two operands are
automatically converted into elements of the NCA (the operands themselves are
not modified). This upward conversion along an ancestor chain is easy: it just
requires the successive multiplication by the defining matrix of each
``Submodule``.

Conversely, downward conversion, i.e. representing a given ``ModuleElement`` in
a submodule, is also supported -- namely by the
:py:meth:`~sympy.polys.numberfields.modules.Submodule.represent`
method -- but is not guaranteed to succeed in general, since the given element
may not belong to the submodule. The main circumstance in which this issue
tends to arise is with multiplication, since modules, while closed under
addition, need not be closed under multiplication.


Multiplication
--------------

Generally speaking, a module need not be closed under multiplication, i.e. need
not form a ring. However, many of the modules we work with in the context of
number fields are in fact rings, and our classes do support multiplication.

Specifically, any ``Module`` can attempt to compute its own multiplication
table, but this does not happen unless an attempt is made to multiply two
``ModuleElement`` instances belonging to it.

>>> A = PowerBasis(T)
>>> print(A._mult_tab is None)
True
>>> a = A(0)*A(1)
>>> print(A._mult_tab is None)
False

Every ``PowerBasis`` is, by its nature, closed under multiplication, so
instances of ``PowerBasis`` can always successfully compute their
multiplication table.

When a ``Submodule`` attempts to compute its multiplication table, it converts
each of its own generators into elements of its parent module, multiplies
them there, in every possible pairing, and then tries to represent the results
in itself, i.e. as $\mathbb{Z}$-linear combinations over its own generators.
This will succeed if and only if the submodule is in fact closed under
multiplication.


Module Homomorphisms
====================

Many important number theoretic algorithms require the calculation of the
kernel of one or more module homomorphisms. Accordingly we have several
lightweight classes, :py:class:`~.ModuleHomomorphism`,
:py:class:`~.ModuleEndomorphism`, :py:class:`~.InnerEndomorphism`, and
:py:class:`~.EndomorphismRing`, which provide the minimal necessary machinery
to support this.

"""

from sympy.core.numbers import igcd, ilcm
from sympy.core.symbol import Dummy
from sympy.polys.polytools import Poly
from sympy.polys.densetools import dup_clear_denoms
from sympy.polys.domains.finitefield import FF
from sympy.polys.domains.rationalfield import QQ
from sympy.polys.domains.integerring import ZZ
from sympy.polys.matrices.domainmatrix import DomainMatrix
from sympy.polys.matrices.exceptions import DMBadInputError
from sympy.polys.matrices.normalforms import hermite_normal_form
from sympy.polys.polyerrors import CoercionFailed
from sympy.polys.polyutils import IntegerPowerable
from .exceptions import ClosureFailure, MissingUnityError
from .utilities import AlgIntPowers, is_int, is_rat, get_num_denom


def to_col(coeffs):
    """Transform a list of integer coefficients into a column vector."""
    return DomainMatrix([[ZZ(c) for c in coeffs]], (1, len(coeffs)), ZZ).transpose()


class Module:
    """
    Generic finitely-generated module.

    This is an abstract base class, and should not be instantiated directly.
    The two concrete subclasses are :py:class:`~.PowerBasis` and
    :py:class:`~.Submodule`.

    Every ``Submodule`` is derived from another module,
    referenced by its ``parent`` attribute. If ``S`` is a submodule, then
    we refer to ``S.parent``, ``S.parent.parent``, and so on, as the
    "ancestors" of ``S``. Thus, every ``Module`` is either a ``PowerBasis``
    or a ``Submodule``, some ancestor of which is a ``PowerBasis``.
    """

    @property
    def n(self):
        """The number of generators of this module."""
        raise NotImplementedError

    def mult_tab(self):
        """
        Get the multiplication table for this module (if closed under mult).

        Examples
        ========

        >>> from sympy.polys import Poly, cyclotomic_poly
        >>> from sympy.polys.numberfields.modules import PowerBasis
        >>> T = Poly(cyclotomic_poly(5))
        >>> A = PowerBasis(T)
        >>> print(A.mult_tab())  # doctest: +SKIP
        {0: {0: [1, 0, 0, 0], 1: [0, 1, 0, 0], 2: [0, 0, 1, 0],     3: [0, 0, 0, 1]},
                          1: {1: [0, 0, 1, 0], 2: [0, 0, 0, 1],     3: [-1, -1, -1, -1]},
                                           2: {2: [-1, -1, -1, -1], 3: [1, 0, 0, 0]},
                                                                3: {3: [0, 1, 0, 0]}}

        Returns
        =======

        Dictionary of dictionaries of lists ``M``, representing the upper
        triangular half of the multiplication table. In other words, if
        ``0 <= i <= j < self.n``, then ``M[i][j]`` is the list ``c`` of
        coefficients such that
        ``g[i] * g[j] == sum(c[k]*g[k], k in range(self.n))``,
        where ``g`` is the list of generators of this module.
        If ``j < i`` then ``M[i][j]`` is undefined.

        """
        raise NotImplementedError

    @property
    def parent(self):
        """
        The parent ``Module``, if any, for this ``Module``.

        For a ``Submodule`` this is its ``parent`` attribute; for a
        ``PowerBasis`` this is ``None``.

        See Also
        ========

        Module

        """
        return None

    def represent(self, elt):
        """
        Represent an element of an ancestor module as an integer-linear
        combination over the generators of this module.

        See Also
        ========

        Module

        """
        raise NotImplementedError

    def ancestors(self, include_self=False):
        """
        Return the list of ancestor ``Module``s of this ``Module``, from the
        foundational ``PowerBasis`` downward, optionally including ``self``.

        See Also
        ========

        Module

        """
        c = self.parent
        a = [] if c is None else c.ancestors(include_self=True)
        if include_self:
            a.append(self)
        return a

    def power_basis_ancestor(self):
        """
        Return the ``PowerBasis`` that is an ancestor of this module.

        See Also
        ========

        Module

        """
        if isinstance(self, PowerBasis):
            return self
        c = self.parent
        if c is not None:
            return c.power_basis_ancestor()
        return None

    def nearest_common_ancestor(self, other):
        """
        Locate the nearest common ancestor of this Module and another.

        Returns
        =======

        Module or None

        See Also
        ========

        Module

        """
        sA = self.ancestors(include_self=True)
        oA = other.ancestors(include_self=True)
        nca = None
        for sa, oa in zip(sA, oA):
            if sa == oa:
                nca = sa
            else:
                break
        return nca

    def is_compat_col(self, col):
        """Say whether *col* is a suitable column vector for this module."""
        return isinstance(col, DomainMatrix) and col.shape == (self.n, 1) and col.domain.is_ZZ

    def __call__(self, spec, denom=1):
        """
        Generate a ModuleElement belonging to this module.

        Parameters
        ==========

        spec: either a compatible column vector, or else an integer ``j``,
          ``0 <= j < self.n``, which we interpret as the jth elementary basis vector.
          In either case, we return the ModuleElement defined by this vector.
        denom: optional denominator for the ModuleElement.

        """
        if isinstance(spec, int) and 0 <= spec < self.n:
            spec = DomainMatrix.eye(self.n, ZZ)[:, spec].to_dense()
        if not self.is_compat_col(spec):
            raise ValueError('Expect compatible column vector specification.')
        return make_mod_elt(self, spec, denom=denom)

    def starts_with_unity(self):
        """Say whether the module's first generator equals unity."""
        raise NotImplementedError

    def basis_elements(self):
        """
        Get list of ``ModuleElement`` being the generators of this module.
        """
        return [self(j) for j in range(self.n)]

    def zero(self):
        """Return a ModuleElement representing zero."""
        return self(0) * 0

    def one(self):
        """
        Return a ModuleElement representing unity, and belonging to the first
        ancestor of this module (including itself) that starts with unity.
        """
        return self.element_from_rational(1)

    def element_from_rational(self, a):
        """
        Return a ModuleElement representing a rational number.

        Explanation
        ===========

        The returned ModuleElement will belong to the first module on this
        module's ancestor chain (including this module itself) that starts
        with unity.

        Parameters
        ==========

        a : int, ZZ, QQ

        Returns
        =======

        ModuleElement

        """
        raise NotImplementedError

    def submodule_from_gens(self, gens, hnf=True, hnf_modulus=None):
        """
        Form the submodule generated by a list of :py:class:`~.ModuleElement`
        belonging to this Module.

        Parameters
        ==========

        gens : list of :py:class:`~.ModuleElement` belonging to this Module.
        hnf : boolean, optional (default=True)
            If True, we will reduce the matrix into Hermite normal form before
            forming the Submodule.
        hnf_modulus : int, None, optional (default=None)
            Modulus for use in the HNF reduction algorithm.

        Returns
        =======

        Submodule

        """
        if not all(g.module == self for g in gens):
            raise ValueError('Generators must belong to this module.')
        n = len(gens)
        if n == 0:
            raise ValueError('Need at least one generator.')
        m = gens[0].n
        d = gens[0].denom if n == 1 else ilcm(*[g.denom for g in gens])
        B = DomainMatrix.zeros((m, 0), ZZ).hstack(*[(d // g.denom) * g.col for g in gens])
        if hnf:
            B = hermite_normal_form(B, D=hnf_modulus)
        return self.submodule_from_matrix(B, denom=d)

    def submodule_from_matrix(self, B, denom=1):
        """
        Form the submodule generated by elements of this Module indicated by
        the columns of a matrix, with an optional denominator.

        Parameters
        ==========

        B: DomainMatrix over :ref:`ZZ` whose columns define, along with the denom, the
          elements of this module that generate the submodule. Thus, the number
          of rows of *B* must equal the number of generators of the present module.
        denom: see *B*.

        Returns
        =======

        Submodule

        Raises
        ======

        ValueError if the given matrix *B* is not over :ref:`ZZ` or its dimension is
        mismatched with this module.

        """
        m, n = B.shape
        if not B.domain.is_ZZ:
            raise ValueError('Matrix must be over ZZ.')
        if not m == self.n:
            raise ValueError('Matrix row count must match base module.')
        return Submodule(self, B, denom=denom)

    def whole_submodule(self):
        """Return a submodule equal to this entire module."""
        B = DomainMatrix.eye(self.n, ZZ)
        return self.submodule_from_matrix(B)

    def endomorphism_ring(self):
        """Form the endomorphism ring for this module."""
        return EndomorphismRing(self)


class PowerBasis(Module):
    """The module generated by the powers of an algebraic integer."""

    def __init__(self, T):
        """
        Parameters
        ==========

        T: monic, irreducible, univariate polynomial over :ref:`ZZ`

        """
        self.T = T
        self._n = T.degree()
        self._mult_tab = None

    def __repr__(self):
        return f'PowerBasis({self.T.as_expr()})'

    def __eq__(self, other):
        if isinstance(other, PowerBasis):
            return self.T == other.T
        return NotImplemented

    @property
    def n(self):
        return self._n

    def mult_tab(self):
        if self._mult_tab is None:
            self.compute_mult_tab()
        return self._mult_tab

    def compute_mult_tab(self):
        theta_pow = AlgIntPowers(self.T)
        M = {}
        n = self.n
        for u in range(n):
            M[u] = {}
            for v in range(u, n):
                M[u][v] = theta_pow[u + v]
        self._mult_tab = M

    def represent(self, elt):
        """
        In our system, to "represent" always means to write a ModuleElement as
        a ``ZZ``-linear combination over the generators of the present Module.
        Furthermore, the incoming ModuleElement must belong to an ancestor of
        the present Module.

        Therefore, when the present Module is a PowerBasis, it is an odd case,
        and one which tends not to arise in practice, except maybe when using a
        ModuleEndomorphism on a PowerBasis.

        In such a case, (1) the incoming ModuleElement must belong to the
        PowerBasis itself (since the latter has no proper ancestors) and
        (2) it is "representable" iff it belongs to ``ZZ[theta]`` (although
        generally PowerBasisElements may represent any element of ``QQ[theta]``,
        i.e. any algebraic number).
        """
        if elt.module == self and elt.denom == 1:
            return elt.column()
        else:
            raise ClosureFailure('Element not representable in ZZ[theta].')

    def starts_with_unity(self):
        return True

    def element_from_rational(self, a):
        return self(0) * a

    def element_from_poly(self, f):
        """
        Produce an element of this module, representing *f* after reduction mod
        our defining minimal polynomial.

        Parameters
        ==========

        f: Poly over :ref:`ZZ` in same var as our defining poly.

        Returns
        =======

        PowerBasisElement

        """
        n, k = self.n, f.degree()
        if k >= n:
            f = f % self.T
        if f == 0:
            return self.zero()
        d, c = dup_clear_denoms(f.rep.rep, QQ, convert=True)
        c = list(reversed(c))
        ell = len(c)
        z = [ZZ(0)] * (n - ell)
        col = to_col(c + z)
        return self(col, denom=d)


class Submodule(Module, IntegerPowerable):
    """A submodule of another module."""

    def __init__(self, parent, matrix, denom=1, mult_tab=None):
        """
        Parameters
        ==========

        parent: the Module from which this one is derived
        matrix: the DomainMatrix over :ref:`ZZ` whose columns define this submodule's
          generators as linear combinations over the parent's generators.
        denom: optional denominator.
        mult_tab: optional multiplication table for this module.

        """
        self._parent = parent
        self._matrix = matrix
        self._denom = denom
        self._mult_tab = mult_tab
        self._n = matrix.shape[1]
        self._QQ_matrix = None
        self._starts_with_unity = None
        self._is_sq_maxrank_HNF = None

    def __repr__(self):
        r = 'Submodule' + repr(self.matrix.transpose().to_Matrix().tolist())
        if self.denom > 1:
            r += f'/{self.denom}'
        return r

    def reduced(self):
        """Produce a reduced version of this Submodule."""
        if self.denom == 1:
            return self
        g = igcd(self.denom, *self.coeffs)
        if g == 1:
            return self
        return type(self)(self.parent, (self.matrix / g).convert_to(ZZ), denom=self.denom // g, mult_tab=self._mult_tab)

    def discard_before(self, r):
        """
        Produce a new module by discarding all generators before a given index *r*.
        """
        W = self.matrix[:, r:]
        s = self.n - r
        M = None
        mt = self._mult_tab
        if mt is not None:
            M = {}
            for u in range(s):
                M[u] = {}
                for v in range(u, s):
                    M[u][v] = mt[r + u][r + v][r:]
        return Submodule(self.parent, W, denom=self.denom, mult_tab=M)

    @property
    def n(self):
        return self._n

    def mult_tab(self):
        if self._mult_tab is None:
            self.compute_mult_tab()
        return self._mult_tab

    def compute_mult_tab(self):
        gens = self.basis_element_pullbacks()
        M = {}
        n = self.n
        for u in range(n):
            M[u] = {}
            for v in range(u, n):
                M[u][v] = self.represent(gens[u] * gens[v]).flat()
        self._mult_tab = M

    @property
    def parent(self):
        return self._parent

    @property
    def matrix(self):
        return self._matrix

    @property
    def coeffs(self):
        return self.matrix.flat()

    @property
    def denom(self):
        return self._denom

    @property
    def QQ_matrix(self):
        """Matrix over :ref:`QQ`, equal to ``self.matrix / self.denom``, and always dense."""
        if self._QQ_matrix is None:
            self._QQ_matrix = (self.matrix / self.denom).to_dense()
        return self._QQ_matrix

    def starts_with_unity(self):
        if self._starts_with_unity is None:
            self._starts_with_unity = self(0).equiv(1)
        return self._starts_with_unity

    def is_sq_maxrank_HNF(self):
        if self._is_sq_maxrank_HNF is None:
            self._is_sq_maxrank_HNF = is_sq_maxrank_HNF(self._matrix)
        return self._is_sq_maxrank_HNF

    def is_power_basis_submodule(self):
        return isinstance(self.parent, PowerBasis)

    def element_from_rational(self, a):
        if self.starts_with_unity():
            return self(0) * a
        else:
            return self.parent.element_from_rational(a)

    def basis_element_pullbacks(self):
        """Return list of this Submodule's basis elements as elements of the parent."""
        return [e.to_parent() for e in self.basis_elements()]

    def represent(self, elt):
        """
        Represent a ModuleElement belonging to an ancestor module, as a
        :ref:`ZZ`-linear combination over the generators of this Submodule.

        Parameters
        ==========

        elt: ModuleElement to be represented. Must belong to some ancestor
          module of this Submodule (including this Submodule itself).

        Returns
        =======

        DomainMatrix
            This will be a column vector, representing the coefficients of a
            linear combination.

        Raises
        ======

        ClosureFailure if the given element cannot be represented as an element
        of this Submodule.

        """
        if elt.module == self:
            return elt.column()
        elif elt.module == self.parent:
            try:
                # The given element should be a ZZ-linear combination over our
                # basis vectors; however, due to the presence of denominators,
                # we need to solve over QQ.
                A = self.QQ_matrix
                b = elt.QQ_col
                x = A._solve(b)[0].transpose()
                x = x.convert_to(ZZ)
            except DMBadInputError:
                raise ClosureFailure('Element outside QQ-span of this basis.')
            except CoercionFailed:
                raise ClosureFailure('Element in QQ-span but not ZZ-span of this basis.')
            return x
        elif isinstance(self.parent, Submodule):
            coeffs_in_parent = self.parent.represent(elt)
            parent_element = self.parent(coeffs_in_parent)
            return self.represent(parent_element)
        else:
            raise ClosureFailure('Element outside ancestor chain of this module.')

    def is_compat_submodule(self, other):
        return isinstance(other, Submodule) and other.parent == self.parent

    def __eq__(self, other):
        if self.is_compat_submodule(other):
            return other.QQ_matrix == self.QQ_matrix
        return NotImplemented

    def add(self, other, hnf=True, hnf_modulus=None):
        d, e = self.denom, other.denom
        m = ilcm(d, e)
        a, b = m // d, m // e
        B = (a * self.matrix).hstack(b * other.matrix)
        if hnf:
            B = hermite_normal_form(B, D=hnf_modulus)
        return self.parent.submodule_from_matrix(B, denom=m)

    def __add__(self, other):
        if self.is_compat_submodule(other):
            return self.add(other)
        return NotImplemented

    __radd__ = __add__

    def mul(self, other, hnf=True, hnf_modulus=None):
        if is_rat(other):
            a, b = get_num_denom(other)
            if a == b == 1:
                return self
            else:
                return Submodule(self.parent,
                             self.matrix * a, denom=self.denom * b,
                             mult_tab=None).reduced()
        elif isinstance(other, ModuleElement) and other.module == self.parent:
            # The submodule is multiplied by an element of the parent module.
            # We presume this means we want a new submodule of the parent module.
            gens = [other * e for e in self.basis_element_pullbacks()]
            return self.parent.submodule_from_gens(gens, hnf=hnf, hnf_modulus=hnf_modulus)
        elif self.is_compat_submodule(other):
            # This case usually means you're multiplying ideals, and want another
            # ideal, i.e. another submodule of the same parent module.
            alphas, betas = self.basis_element_pullbacks(), other.basis_element_pullbacks()
            gens = [a * b for a in alphas for b in betas]
            return self.parent.submodule_from_gens(gens, hnf=hnf, hnf_modulus=hnf_modulus)
        return NotImplemented

    def __mul__(self, other):
        return self.mul(other)

    __rmul__ = __mul__

    def _first_power(self):
        return self


def is_sq_maxrank_HNF(dm):
    """
    Say whether a DomainMatrix is in that special case of Hermite Normal Form,
    in which the matrix is also square and of maximal rank.
    """
    if dm.domain.is_ZZ and dm.is_square and dm.is_upper:
        n = dm.shape[0]
        for i in range(n):
            d = dm[i, i].element
            if d <= 0:
                return False
            for j in range(i + 1, n):
                if not (0 <= dm[i, j].element < d):
                    return False
        return True
    return False


def make_mod_elt(module, col, denom=1):
    """
    Factory function which builds a ModuleElement, but ensures that it is a
    PowerBasisElement if the module is a PowerBasis.
    """
    if isinstance(module, PowerBasis):
        return PowerBasisElement(module, col, denom=denom)
    else:
        return ModuleElement(module, col, denom=denom)


class ModuleElement(IntegerPowerable):
    """
    Element of a :py:class:`~.Module`.

    NOTE: Should not be constructed directly. Use :py:func:`make_mod_elt()`
    factory function instead.
    """

    def __init__(self, module, col, denom=1):
        """
        Parameters
        ==========

        module: Module
        col: column vector over :ref:`ZZ`

        """
        self.module = module
        self.col = col
        self.denom = denom
        self._QQ_col = None

    def __repr__(self):
        r = str([int(c) for c in self.col.flat()])
        if self.denom > 1:
            r += f'/{self.denom}'
        return r

    def reduced(self):
        """
        Produce a reduced version of this ModuleElement, i.e. one in which the
        gcd of the denominator together with all numerator coefficients is 1.
        """
        if self.denom == 1:
            return self
        g = igcd(self.denom, *self.coeffs)
        if g == 1:
            return self
        return type(self)(self.module,
                            (self.col / g).convert_to(ZZ),
                            denom=self.denom // g)

    def reduced_mod_p(self, p):
        """
        Produce a version of this ModuleElement in which all numerator coefficients
        have been reduced mod *p*.
        """
        return make_mod_elt(self.module,
                            self.col.convert_to(FF(p)).convert_to(ZZ),
                            denom=self.denom)

    @classmethod
    def from_int_list(cls, module, coeffs, denom=1):
        """Make a module element from a list of ints (instead of a column vect)."""
        col = to_col(coeffs)
        return cls(module, col, denom=denom)

    @property
    def n(self):
        """The length of this element's column."""
        return self.module.n

    def __len__(self):
        return self.n

    def column(self, domain=None):
        """Get a copy of this element's column, optionally converting to a domain."""
        return self.col.convert_to(domain)

    @property
    def coeffs(self):
        return self.col.flat()

    @property
    def QQ_col(self):
        """Column vector over :ref:`QQ`, equal to ``self.col / self.denom``, and always dense."""
        if self._QQ_col is None:
            self._QQ_col = (self.col / self.denom).to_dense()
        return self._QQ_col

    def to_parent(self):
        """Transform into a ModuleElement belonging to the parent of our module."""
        if not isinstance(self.module, Submodule):
            raise ValueError('Not an element of a Submodule.')
        return make_mod_elt(
            self.module.parent, self.module.matrix * self.col,
            denom=self.module.denom * self.denom)

    def to_ancestor(self, anc):
        """Transform into a ModuleElement belonging to a given ancestor of our module."""
        if anc == self.module:
            return self
        else:
            return self.to_parent().to_ancestor(anc)

    def over_power_basis(self):
        """Transform into a PowerBasisElement over our PowerBasis ancestor."""
        e = self
        while not isinstance(e.module, PowerBasis):
            e = e.to_parent()
        return e

    def is_compat(self, other):
        """Test whether other is another ModuleElement with same module."""
        return isinstance(other, ModuleElement) and other.module == self.module

    def make_compat(self, other):
        """
        Try to make a compatible pair of ``ModuleElement``, one equivalent to
        this one, and one equivalent to the other. This means finding the
        nearest common ancestor Module for the pair of elements, and
        representing each one there.

        Returns
        =======

        Pair of ``ModuleElement`` belonging to a common Module, or else the
        pair ``(None, None)`` if no common ancestor could be found.

        """
        if self.module == other.module:
            return self, other
        nca = self.module.nearest_common_ancestor(other.module)
        if nca is not None:
            return self.to_ancestor(nca), other.to_ancestor(nca)
        return None, None

    def __eq__(self, other):
        if self.is_compat(other):
            return self.QQ_col == other.QQ_col
        return NotImplemented

    def equiv(self, other):
        if self == other:
            return True
        elif isinstance(other, ModuleElement):
            a, b = self.make_compat(other)
            if a is None:
                return False
            return a == b
        elif is_rat(other):
            if isinstance(self, PowerBasisElement):
                return self == self.module(0) * other
            else:
                return self.over_power_basis().equiv(other)
        return False

    def __add__(self, other):
        """
        A ModuleElement can be added to a rational number, or to another
        ModuleElement.

        When the other summand is a rational number, it will be converted into
        a ModuleElement (belonging to the first ancestor of this module that
        starts with unity).

        In all cases, the sum belongs to the nearest common ancestor (NCA) of
        the modules of the two summands. If the NCA does not exist, we return
        ``NotImplemented``.
        """
        if self.is_compat(other):
            d, e = self.denom, other.denom
            m = ilcm(d, e)
            u, v = m // d, m // e
            col = to_col([u * a + v * b for a, b in zip(self.coeffs, other.coeffs)])
            return type(self)(self.module, col, denom=m).reduced()
        elif isinstance(other, ModuleElement):
            a, b = self.make_compat(other)
            if a is None:
                return NotImplemented
            return a + b
        elif is_rat(other):
            return self + self.module.element_from_rational(other)
        return NotImplemented

    __radd__ = __add__

    def __neg__(self):
        return self * -1

    def __sub__(self, other):
        return self + (-other)

    def __rsub__(self, other):
        return -self + other

    def __mul__(self, other):
        """
        A ModuleElement can be multiplied by a rational number, or by another
        ModuleElement.

        When the multiplier is a rational number, the product is computed by
        operating directly on the coefficients of this ModuleElement.

        When the multiplier is another ModuleElement, the product will belong
        to the nearest common ancestor (NCA) of the modules of the two operands,
        and that NCA must have a multiplication table. If the NCA does not
        exist, we return ``NotImplemented``. If the NCA does not have a mult.
        table, ``ClosureFailure`` will be raised.
        """
        if self.is_compat(other):
            M = self.module.mult_tab()
            A, B = self.col.flat(), other.col.flat()
            n = self.n
            C = [0] * n
            for u in range(n):
                for v in range(u, n):
                    c = A[u] * B[v]
                    if v > u:
                        c += A[v] * B[u]
                    if c != 0:
                        R = M[u][v]
                        for k in range(n):
                            C[k] += c * R[k]
            d = self.denom * other.denom
            return self.from_int_list(self.module, C, denom=d)
        elif isinstance(other, ModuleElement):
            a, b = self.make_compat(other)
            if a is None:
                return NotImplemented
            return a * b
        elif is_rat(other):
            a, b = get_num_denom(other)
            if a == b == 1:
                return self
            else:
                return make_mod_elt(self.module,
                                 self.col * a, denom=self.denom * b).reduced()
        return NotImplemented

    __rmul__ = __mul__

    def _zeroth_power(self):
        return self.module.one()

    def _first_power(self):
        return self

    def __floordiv__(self, a):
        if is_rat(a):
            a = QQ(a)
            return self * (1/a)
        return NotImplemented

    def __rfloordiv__(self, a):
        return a // self.over_power_basis()

    def __mod__(self, a):
        r"""
        Reducing mod an integer *a* reduces all numerator coeffs mod $d a$, where
        $d$ is our denominator.

        For example, if we represent

        $$ \frac{15 a_1 + a_0}{2} = \frac{a_1 + a_0}{2} + 7 a_1 $$

        then reducing mod 7 should mean throwing away that part that is a poly
        in the basis elements $a_i$, with content divisible by 7. But this is
        achieved by reducing our coeffs mod $2 \times 7$.
        """
        if is_int(a):
            m = a * self.denom
            col = to_col([c % m for c in self.coeffs])
            return type(self)(self.module, col, denom=self.denom)
        return NotImplemented


class PowerBasisElement(ModuleElement):
    """
    Subclass for ModuleElement instances whose module is a PowerBasis.
    """

    @property
    def T(self):
        return self.module.T

    def numerator(self, x=None):
        """Obtain the numerator as a polynomial over :ref:`ZZ`."""
        x = x or self.T.gen
        return Poly(reversed(self.coeffs), x, domain=ZZ)

    def poly(self, x=None):
        """Obtain the number as a polynomial over :ref:`QQ`."""
        return self.numerator(x=x) // self.denom

    def norm(self, T=None):
        """Compute the norm of this number."""
        T = T or self.T
        x = T.gen
        A = self.numerator(x=x)
        return T.resultant(A) // self.denom ** self.n

    def inverse(self):
        f = self.poly()
        f_inv = f.invert(self.T)
        return self.module.element_from_poly(f_inv)

    def __rfloordiv__(self, a):
        return self.inverse() * a

    def _negative_power(self, e, modulo=None):
        return self.inverse() ** abs(e)


class ModuleHomomorphism:
    """A homomorphism from one module to another."""

    def __init__(self, domain, codomain, mapping):
        """
        Parameters
        ==========

        domain: :py:class:`~.Module` being the domain of the mapping.

        codomain: :py:class:`~.Module` being the codomain of the mapping.

        mapping: an arbitrary callable is accepted, but should be chosen so as
            to represent an actual module homomorphism. In particular, should
            accept elements of ``domain`` and return elements of ``codomain``.

        Examples
        ========

        >>> from sympy import Poly, cyclotomic_poly
        >>> from sympy.polys.numberfields.modules import PowerBasis, ModuleHomomorphism
        >>> T = Poly(cyclotomic_poly(5))
        >>> A = PowerBasis(T)
        >>> B = A.submodule_from_gens([2*A(j) for j in range(4)])
        >>> phi = ModuleHomomorphism(A, B, lambda x: 6*x)
        >>> print(phi.matrix())  # doctest: +SKIP
        DomainMatrix([[3, 0, 0, 0], [0, 3, 0, 0], [0, 0, 3, 0], [0, 0, 0, 3]], (4, 4), ZZ)

        """
        self.domain = domain
        self.codomain = codomain
        self.mapping = mapping

    def matrix(self, modulus=None):
        """
        Compute the matrix of this homomorphism.

        Parameters
        ==========

        modulus: (optional) a positive prime number $p$ if the matrix should
            be reduced mod $p$.

        Returns
        =======

        :py:class:`~.DomainMatrix` over :ref:`ZZ`, or over :ref:`GF(p)` if a
            modulus was given.

        """
        basis = self.domain.basis_elements()
        cols = [self.codomain.represent(self.mapping(elt)) for elt in basis]
        if not cols:
            return DomainMatrix.zeros((self.codomain.n, 0), ZZ).to_dense()
        M = cols[0].hstack(*cols[1:])
        if modulus:
            M = M.convert_to(FF(modulus))
        return M

    def kernel(self, modulus=None):
        """
        Compute a Submodule representing the kernel of this homomorphism.

        Parameters
        ==========

        modulus: (optional) a positive prime number $p$ if the kernel should
            be computed mod $p$.

        Returns
        =======

        :py:class:`~.Submodule` whose generators span the kernel of this
            homomorphism over :ref:`ZZ`, or over :ref:`GF(p)` if a
            modulus was given.

        """
        M = self.matrix(modulus=modulus)
        if modulus is None:
            M = M.convert_to(QQ)
        # Note: Even when working over a finite field, what we want here is
        # the pullback into the integers, so in this case the conversion to ZZ
        # below is appropriate. When working over ZZ, the kernel should be a
        # ZZ-submodule, so, while the conversion to QQ above was required in
        # order for the nullspace calculation to work, conversion back to ZZ
        # afterward should always work.
        # TODO:
        #  Watch <https://github.com/sympy/sympy/issues/21834>, which calls
        #  for fraction-free algorithms. If this is implemented, we can skip
        #  the conversion to `QQ` above.
        K = M.nullspace().convert_to(ZZ).transpose()
        return self.domain.submodule_from_matrix(K)


class ModuleEndomorphism(ModuleHomomorphism):
    """A homomorphism from one module to itself."""

    def __init__(self, domain, mapping):
        """
        Parameters
        ==========

        domain: :py:class:`~.Module` being the common domain and codomain of
            the mapping.

        mapping: an arbitrary callable is accepted, but should be chosen so as
            to represent an actual module endomorphism. In particular, should
            accept and return elements of ``domain``.

        """
        super().__init__(domain, domain, mapping)


class InnerEndomorphism(ModuleEndomorphism):
    """
    An inner endomorphism on a module, i.e. the endomorphism corresponding to
    multiplication by a fixed element.
    """

    def __init__(self, domain, multiplier):
        r"""
        Parameters
        ==========

        domain: :py:class:`~.Module` being the domain and codomain of the
            endomorphism.

        multiplier: :py:class:`~.ModuleElement` $a$ defining the mapping
            as $x \mapsto a x$.

        """
        super().__init__(domain, lambda x: multiplier * x)
        self.multiplier = multiplier


class EndomorphismRing:
    """The ring of endomorphisms on a module."""

    def __init__(self, domain):
        """
        Parameters
        ==========

        domain: :py:class:`~.Module` being the domain and codomain of the
            endomorphisms.

        """
        self.domain = domain

    def inner_endomorphism(self, multiplier):
        """
        Form an inner endomorphism belonging to this endomorphism ring.

        Parameters
        ==========

        multiplier: :py:class:`~.ModuleElement` $a$ defining the inner
            endomorphism $x \mapsto a x$.

        Returns
        =======

        :py:class:`~.InnerEndomorphism`

        """
        return InnerEndomorphism(self.domain, multiplier)

    def represent(self, element):
        r"""
        Represent an element of this endomorphism ring, as a single column
        vector.

        Explanation
        ===========

        Let $M$ be a module, and $E$ its ring of endomorphisms. Let $N$ be
        another module, and consider a homomorphism $\varphi: N \rightarrow E$.
        In the event that $\varphi$ is to be represented by a matrix $A$, each
        column of $A$ must represent an element of $E$. This is possible when
        the elements of $E$ are themselves representable as matrices, by
        stacking the columns of such a matrix into a single column.

        This method supports calculating such matrices $A$, by representing
        an element of this endomorphism ring first as a matrix, and then
        stacking that matrix's columns into a single column.

        Parameters
        ==========

        element: :py:class:`~.ModuleEndomorphism` belonging to this ring.

        Returns
        =======

        DomainMatrix, in the shape of a column vector.

        """
        if isinstance(element, ModuleEndomorphism) and element.domain == self.domain:
            M = element.matrix()
            # Transform the matrix into a single column, which should reproduce
            # the original columns, one after another.
            m, n = M.shape
            if n == 0:
                return M
            return M[:, 0].vstack(*[M[:, j] for j in range(1, n)])
        raise NotImplementedError


def find_min_poly(alpha, domain, x=None, powers=None):
    """
    Find a polynomial of least degree (not necessarily irreducible) satisfied
    by an element of a finitely-generated ring with unity.

    Parameters
    ==========

    alpha: :py:class:`~.ModuleElement` whose min poly is to be found, and whose
        module has multiplication and starts with unity.

    domain: The desired :py:class:`~.Domain` of the polynomial.

    x: (optional) :py:class:`~.Symbol`, the desired variable for the
        polynomial.

    powers: (optional) If desired, pass an empty list. The powers of alpha
        (as ModuleElement instances) from the zeroth up to the degree of the
        min poly will be recorded here, as we compute them.

    Returns
    =======

    :py:class:`~.Poly`, being the minimal polynomial for alpha, or ``None`` if
        no polynomial could be found over the desired domain.

    Raises
    ======

    MissingUnityError
        If the module to which alpha belongs does not start with unity.
    ClosureFailure
        If the module to which alpha belongs is not closed under multiplication.

    """
    R = alpha.module
    if not R.starts_with_unity():
        raise MissingUnityError("alpha must belong to finitely generated ring with unity.")
    if powers is None:
        powers = []
    one = R(0)
    powers.append(one)
    powers_matrix = one.column(domain=domain)
    ak = alpha
    m = None
    for k in range(1, R.n + 1):
        powers.append(ak)
        ak_col = ak.column(domain=domain)
        try:
            X = powers_matrix._solve(ak_col)[0]
        except DMBadInputError:
            # This means alpha^k still isn't in the domain-span of the lower powers.
            powers_matrix = powers_matrix.hstack(ak_col)
            ak *= alpha
        else:
            # alpha^k is in the domain-span of the lower powers, so we have found a
            # minimal-degree poly for alpha.
            coeffs = [1] + [-c for c in reversed(X.to_list_flat())]
            x = x or Dummy('x')
            if domain.is_FF:
                m = Poly(coeffs, x, modulus=domain.mod)
            else:
                m = Poly(coeffs, x, domain=domain)
            break
    return m
