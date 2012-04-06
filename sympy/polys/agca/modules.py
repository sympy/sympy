"""
Computations with modules over polynomial rings.

This module implements various classes that encapsulate groebner basis
computations for modules. Most of them should not be instantiated by hand.
Instead, use the constructing routines on objects you already have.

For example, to construct a free module over ``QQ[x, y]``, call
``QQ[x, y].free_module(rank)`` instead of the ``FreeModule`` constructor.
In fact ``FreeModule`` is an abstract base class that should not be
instantiated, the ``free_module`` method instead returns the implementing class
``FreeModulePolyRing``.

In general, the abstract base classes implement most functionality in terms of
a few non-implemented methods. The concrete base classes supply only these
non-implemented methods. They may also supply new implementations of the
convenience methods, for example if there are faster algorithms available.
"""

from copy import copy

from sympy.polys.polyerrors import CoercionFailed
from sympy.polys.monomialtools import ProductOrder, monomial_key
from sympy.polys.domains import Field

from sympy.core.compatibility import iterable

# TODO
# - intersection
# - quotient modules
# - modules over quotient rings
# - module quotient / saturation
# - free resoltutions / syzygies
# - homomorphisms
# - finding small/minimal generating sets
# - ...

class Module(object):
    """
    Abstract base class for modules.

    Do not instantiate - use ring explicit constructors instead:

    >>> from sympy import QQ
    >>> from sympy.abc import x
    >>> QQ[x].free_module(2)
    QQ[x]**2

    Attributes:

    - dtype - type of elements
    - ring - containing ring

    Non-implemented methods:

    - submodule
    - quotient_module
    - is_zero
    - is_submodule

    The method convert likely needs to be changed in subclasses.
    """

    def __init__(self, ring):
        self.ring = ring

    def convert(self, elem, M=None):
        """
        Convert ``elem`` into internal representation of this module.

        If ``M`` is not None, it should be a module containing it.
        """
        if not isinstance(elem, self.dtype):
            raise CoercionFailed
        return elem

    def submodule(self, *gens):
        """Generate a submodule."""
        raise NotImplementedError

    def quotient_module(self, other):
        """Generate a quotient module."""
        raise NotImplementedError

    def contains(self, elem):
        """Return True if ``elem`` is an element of this module."""
        try:
            self.convert(elem)
            return True
        except CoercionFailed:
            return False

    def __contains__(self, elem):
        return self.contains(elem)

    def subset(self, other):
        """
        Returns True if ``other`` is is a subset of ``self``.

        >>> from sympy.abc import x
        >>> from sympy import QQ
        >>> F = QQ[x].free_module(2)
        >>> F.subset([(1, x), (x, 2)])
        True
        >>> F.subset([(1/x, x), (x, 2)])
        False
        """
        return all(self.contains(x) for x in other)

    def __eq__(self, other):
        return self.is_submodule(other) and other.is_submodule(self)

    def __ne__(self, other):
        return not (self == other)

    def is_zero(self):
        """Returns True if ``self`` is a zero module."""
        raise NotImplementedError

    def is_submodule(self, other):
        """Returns True if ``other`` is a submodule of ``self``."""
        raise NotImplementedError

class ModuleElement(object):
    """
    Base class for module element wrappers.

    Use this class to wrap primitive data types as module elements. It stores
    a reference to the containing module, and implements all the arithmetic
    operators.

    Attributes:

    - module - containing module
    - data - internal data

    Methods that likely need change in subclasses:

    - add
    - mul
    - div
    - eq
    """

    def __init__(self, module, data):
        self.module = module
        self.data = data

    def add(self, d1, d2):
        """Add data ``d1`` and ``d2``."""
        return d1 + d2

    def mul(self, m, d):
        """Multiply module data ``m`` by coefficient d."""
        return m * d

    def div(self, m, d):
        """Divide module data ``m`` by coefficient d."""
        return m / d

    def eq(self, d1, d2):
        """Return true if d1 and d2 represent the same element."""
        return d1 == d2

    def __add__(self, om):
        if not isinstance(om, self.__class__) or om.module != self.module:
            om = self.module.convert(om)
        return self.__class__(self.module, self.add(self.data, om.data))

    __radd__ = __add__

    def __neg__(self):
        return self.__class__(self.module, self.mul(self.data,
                       self.module.ring.convert(-1)))

    def __sub__(self, om):
        return self.__add__(-om)

    def __rsub__(self, om):
        return (-self).__add__(om)

    def __mul__(self, o):
        if not isinstance(o, self.module.ring.dtype):
            o = self.module.ring.convert(o)
        return self.__class__(self.module, self.mul(self.data, o))

    __rmul__ = __mul__

    def __div__(self, o):
        if not isinstance(o, self.module.ring.dtype):
            o = self.module.ring.convert(o)
        return self.__class__(self.module, self.div(self.data, o))

    __truediv__ = __div__

    def __eq__(self, om):
        if not isinstance(om, self.__class__) or om.module != self.module:
            om = self.module.convert(om)
        return self.eq(self.data, om.data)

    def __ne__(self, om):
        return not self.__eq__(om)

class FreeModuleElement(ModuleElement):
    """Element of a free module. Data stored as a tuple."""

    def add(self, d1, d2):
        return tuple(x + y for x, y in zip(d1, d2))

    def mul(self, d, p):
        return tuple(x * p for x in d)

    def div(self, d, p):
        return tuple(x / p for x in d)

    def __repr__(self):
        from sympy import sstr
        return '[' + ', '.join(sstr(x) for x in self.data) + ']'

    def __iter__(self):
        return self.data.__iter__()

    def __getitem__(self, idx):
        return self.data[idx]

class FreeModule(Module):
    """
    Abstract base class for free modules.

    Additional attributes:

    - rank - rank of the free module

    Non-implemented methods:

    - quotient_module
    - submodule
    """

    dtype = FreeModuleElement

    def __init__(self, ring, rank):
        Module.__init__(self, ring)
        self.rank = rank

    def __repr__(self):
        return repr(self.ring) + "**" + repr(self.rank)

    def is_submodule(self, other):
        """
        Returns True if ``other`` is a submodule of ``self``.

        >>> from sympy.abc import x
        >>> from sympy import QQ
        >>> F = QQ[x].free_module(2)
        >>> M = F.submodule([2, x])
        >>> F.is_submodule(F)
        True
        >>> F.is_submodule(M)
        True
        >>> M.is_submodule(F)
        False
        """
        if isinstance(other, SubModule):
            return other.container == self
        if isinstance(other, FreeModule):
            return other.ring == self.ring and other.rank == self.rank
        return False

    def convert(self, elem, M=None):
        """
        Convert ``elem`` into the internal representation.

        This method is called implicitly whenever computations involve elements
        not in the internal representation.

        >>> from sympy.abc import x
        >>> from sympy import QQ
        >>> F = QQ[x].free_module(2)
        >>> F.convert([1, 0])
        [1, 0]
        """
        if isinstance(elem, FreeModuleElement):
            if elem.module == self:
                return elem
            if elem.module.rank != self.rank:
                raise CoercionFailed
            return FreeModuleElement(self,
                     tuple(self.ring.convert(x, elem.module.ring) for x in elem.data))
        elif iterable(elem):
            tpl = tuple(self.ring.convert(x) for x in elem)
            if len(tpl) != self.rank:
                raise CoercionFailed
            return FreeModuleElement(self, tpl)
        elif elem == 0:
            return FreeModuleElement(self, (self.ring.convert(0),)*self.rank)
        else:
            raise CoercionFailed

    def is_zero(self):
        """
        Returns True if ``self`` is a zero module.

        (If, as this implementation assumes, the coefficient ring is not the
        zero ring, then this is equivalent to the rank being zero.)

        >>> from sympy.abc import x
        >>> from sympy import QQ
        >>> QQ[x].free_module(0).is_zero()
        True
        >>> QQ[x].free_module(1).is_zero()
        False
        """
        return self.rank == 0

    def basis(self):
        """
        Return a set of basis elements.

        >>> from sympy.abc import x
        >>> from sympy import QQ
        >>> QQ[x].free_module(3).basis()
        ([1, 0, 0], [0, 1, 0], [0, 0, 1])
        """
        from sympy.matrices import eye
        M = eye(self.rank)
        return tuple(self.convert(M.row(i)) for i in range(self.rank))

class FreeModulePolyRing(FreeModule):
    """
    Free module over a generalized polynomial ring.

    Do not instantiate this, use the constructor method of the ring instead:

    >>> from sympy.abc import x
    >>> from sympy import QQ
    >>> F = QQ[x].free_module(3)
    >>> F
    QQ[x]**3
    >>> F.contains([x, 1, 0])
    True
    >>> F.contains([1/x, 0, 1])
    False
    """

    def __init__(self, ring, rank):
        from sympy.polys.domains.polynomialring import PolynomialRingBase
        FreeModule.__init__(self, ring, rank)
        if not isinstance(ring, PolynomialRingBase):
            raise NotImplementedError('This implementation only works over '
                                      + 'polynomial rings, got %s' % ring)
        if not isinstance(ring.dom, Field):
            raise NotImplementedError('Ground domain must be a field, '
                                      + 'got %s' % ring.dom)

    def submodule(self, *gens, **opts):
        """
        Generate a submodule.

        >>> from sympy.abc import x, y
        >>> from sympy import QQ
        >>> M = QQ[x, y].free_module(2).submodule([x, x + y])
        >>> M
        <[x, x + y]>
        >>> M.contains([2*x, 2*x + 2*y])
        True
        >>> M.contains([x, y])
        False
        """
        return SubModulePolyRing(gens, self, **opts)

class SubModule(Module):
    """
    Base class for submodules.

    Attributes:

    - container - containing module
    - gens - generators (subset of containing module)
    - rank - rank of containing module

    Non-implemented methods:

    - _contains
    - _syzygies
    - _in_terms_of_generators
    """

    def __init__(self, gens, container):
        Module.__init__(self, container.ring)
        self.gens = tuple(container.convert(x) for x in gens)
        self.container = container
        self.rank = container.rank
        self.ring = container.ring
        self.dtype = container.dtype

    def __repr__(self):
        return "<" + ", ".join(repr(x) for x in self.gens) + ">"

    def _contains(self, other):
        """Implementation of containment.
           Other is guaranteed to be FreeModuleElement."""
        raise NotImplementedError

    def _syzygies(self):
        """Implementation of syzygy computation wrt self generators."""
        raise NotImplementedError

    def _in_terms_of_generators(self, e):
        """Implementation of expression in terms of generators."""
        raise NotImplementedError

    def convert(self, elem, M=None):
        """
        Convert ``elem`` into the internal represantition.

        Mostly called implicitly.

        >>> from sympy.abc import x
        >>> from sympy import QQ
        >>> M = QQ[x].free_module(2).submodule([1, x])
        >>> M.convert([2, 2*x])
        [2, 2*x]
        """
        r = copy(self.container.convert(elem, M))
        r.module = self
        if not self._contains(r):
            raise CoercionFailed
        return r

    def _intersect(self, other):
        """Implementation of intersection.
           Other is guaranteed to be a submodule of same free module."""
        raise NotImplementedError

    def intersect(self, other):
        """Returns the intersection of ``self`` with submodule ``other``."""
        if not isinstance(other, SubModule):
            raise TypeError('%s is not a SubModule' % other)
        if other.container != self.container:
            raise ValueError('%s is contained in a different free module' % other)
        return self._intersect(other)

    def union(self, other):
        """
        Returns the module generated by the union of ``self`` and ``other``.

        >>> from sympy.abc import x
        >>> from sympy import QQ
        >>> F = QQ[x].free_module(1)
        >>> M = F.submodule([x**2 + x]) # <x(x+1)>
        >>> N = F.submodule([x**2 - 1]) # <(x-1)(x+1)>
        >>> M.union(N) == F.submodule([x+1])
        True
        """
        if not isinstance(other, SubModule):
            raise TypeError('%s is not a SubModule' % other)
        if other.container != self.container:
            raise ValueError('%s is contained in a different free module' % other)
        return self.__class__(self.gens + other.gens, self.container)

    def is_zero(self):
        """
        Return True if ``self`` is a zero module.

        >>> from sympy.abc import x
        >>> from sympy import QQ
        >>> F = QQ[x].free_module(2)
        >>> F.submodule([x, 1]).is_zero()
        False
        >>> F.submodule([0, 0]).is_zero()
        True
        """
        return all(x == 0 for x in self.gens)

    def submodule(self, *gens):
        """
        Generate a submodule.

        >>> from sympy.abc import x
        >>> from sympy import QQ
        >>> M = QQ[x].free_module(2).submodule([x, 1])
        >>> M.submodule([x**2, x])
        <[x**2, x]>
        """
        if not self.subset(gens):
            raise ValueError('%s not a subset of %s' % (gens, self))
        return self.__class__(gens, self.container)

    def is_full_module(self):
        """
        Return True if ``self`` is the entire free module.

        >>> from sympy.abc import x
        >>> from sympy import QQ
        >>> F = QQ[x].free_module(2)
        >>> F.submodule([x, 1]).is_full_module()
        False
        >>> F.submodule([1, 1], [1, 2]).is_full_module()
        True
        """
        return all(self.contains(x) for x in self.container.basis())

    def is_submodule(self, other):
        """
        Returns True if ``other`` is a submodule of ``self``.

        >>> from sympy.abc import x
        >>> from sympy import QQ
        >>> F = QQ[x].free_module(2)
        >>> M = F.submodule([2, x])
        >>> N = M.submodule([2*x, x**2])
        >>> M.is_submodule(M)
        True
        >>> M.is_submodule(N)
        True
        >>> N.is_submodule(M)
        False
        """
        if isinstance(other, SubModule):
            return self.container == other.container and \
                   all(self.contains(x) for x in other.gens)
        if isinstance(other, FreeModule):
            return self.container == other and self.is_full_module()
        return False

    def syzygy_module(self, **opts):
        r"""
        Compute the syzygy module of the generators of ``self``.

        Suppose `M` is generated by `f_1, \dots, f_n` over the ring
        `R`. Consider the homomorphism `\phi: R^n \to M`, given by
        sending `(r_1, \dots, r_n) \to r_1 f_1 + \dots + r_n f_n`.
        The syzygy module is defined to be the kernel of `\phi`.

        The syzygy module is zero iff the generators generate freely a free
        submodule:

        >>> from sympy.abc import x, y
        >>> from sympy import QQ
        >>> QQ[x].free_module(2).submodule([1, 0], [1, 1]).syzygy_module().is_zero()
        True

        A slightly more interesting example:

        >>> M = QQ[x, y].free_module(2).submodule([x, 2*x], [y, 2*y])
        >>> S = QQ[x, y].free_module(2).submodule([y, -x])
        >>> M.syzygy_module() == S
        True
        """
        F = self.ring.free_module(len(self.gens))
        # NOTE we filter out zero syzygies. This is for convenience of the
        # _syzygies function and not meant to replace any real "generating set
        # reduction" algorithm
        return F.submodule(*[x for x in self._syzygies() if F.convert(x) != 0],
                           **opts)

    def in_terms_of_generators(self, e):
        """
        Express element ``e`` of ``self`` in terms of the generators.

        >>> from sympy.abc import x
        >>> from sympy import QQ
        >>> F = QQ[x].free_module(2)
        >>> M = F.submodule([1, 0], [1, 1])
        >>> M.in_terms_of_generators([x, x**2])
        [-x**2 + x, x**2]
        """
        try:
            e = self.convert(e)
        except CoercionFailed:
            raise ValueError('%s is not an element of %s' % (e, self))
        return self._in_terms_of_generators(e)

_subs0 = lambda x: x[0]
_subs1 = lambda x: x[1:]
class ModuleOrder(ProductOrder):
    """A product monomial order with a zeroth term as module index."""

    def __init__(self, o1, o2):
        ProductOrder.__init__(self, (o1, _subs0), (o2, _subs1))

class SubModulePolyRing(SubModule):
    """
    Submodule of a free module over a generalized polynomial ring.

    Do not instantiate this, use the constructor method of FreeModule instead:

    >>> from sympy.abc import x, y
    >>> from sympy import QQ
    >>> F = QQ[x, y].free_module(2)
    >>> F.submodule([x, y], [1, 0])
    <[x, y], [1, 0]>

    Attributes:

    - order - monomial order used
    """

    #self._gb - cached groebner basis

    def __init__(self, gens, container, order="lex"):
        SubModule.__init__(self, gens, container)
        if not isinstance(container, FreeModulePolyRing):
            raise NotImplementedError('This implementation is for submodules of '
                             + 'FreeModulePolyRing, got %s' % container)
        self.order = ModuleOrder(monomial_key(order), self.ring.order)
        self._gb = None

    def __eq__(self, other):
        if isinstance(other, SubModulePolyRing) and self.order != other.order:
            return False
        return SubModule.__eq__(self, other)

    def _groebner(self):
        """Returns a standard basis in sdm form."""
        from sympy.polys.distributedmodules import sdm_groebner, sdm_nf_mora
        if self._gb is None:
            self._gb = tuple(sdm_groebner(
               [self.ring._vector_to_sdm(x, self.order) for x in self.gens],
               sdm_nf_mora, self.order, self.ring.dom))
        return self._gb

    def _groebner_vec(self):
        """Returns a standard basis in element form."""
        return [self.convert(self.ring._sdm_to_vector(x, self.rank)) \
                for x in self._groebner()]

    def _contains(self, x):
        from sympy.polys.distributedmodules import sdm_zero, sdm_nf_mora
        return sdm_nf_mora(self.ring._vector_to_sdm(x, self.order),
                           self._groebner(), self.order, self.ring.dom) == \
               sdm_zero()

    def _syzygies(self):
        """Compute syzygies. See [SCA, algorithm 2.5.4]."""
        # NOTE if self.gens is a standard basis, this can be done more
        #      efficiently using Schreyer's theorem
        from sympy.matrices import eye

        # First bullet point
        k = len(self.gens)
        r = self.rank
        im = eye(k)
        Rkr = self.ring.free_module(r + k)
        newgens = []
        for j, f in enumerate(self.gens):
            m = [0]*(r+k)
            for i, v in enumerate(f):
                m[i] = f[i]
            for i in range(k):
                m[r + i] = im[j, i]
            newgens.append(Rkr.convert(m))
        # Note: we need *descending* order on module index
        F = Rkr.submodule(*newgens, **{'order': 'ilex'})

        # Second bullet point: standard basis of F
        G = F._groebner_vec()

        # Third bullet point: G0 = G intersect the new k components
        G0 = [x[r:] for x in G if all(y == self.ring.convert(0) for y in x[:r])]

        # Fourth and fifth bullet points: we are done
        return G0

    def _in_terms_of_generators(self, e):
        """Expression in terms of generators. See [SCA, 2.8.1]."""
        # NOTE: if gens is a standard basis, this can be done more efficiently
        M = self.ring.free_module(self.rank).submodule(*((e,) + self.gens))
        S = M.syzygy_module(order="ilex") # We want decreasing order!
        G = S._groebner_vec()
        # This list cannot not be empty since e is an element
        e = list(filter(lambda x: self.ring.is_unit(x[0]), G))[0]
        return [-x/e[0] for x in e[1:]]
