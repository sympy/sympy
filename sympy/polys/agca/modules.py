"""Computations with modules over polynomial rings."""

from sympy.polys.polyerrors import CoercionFailed

from sympy.core.compatibility import iterable

class Module(object):
    """Base class for modules."""

    #self.dtype - type representing elements
    #self.ring - containing ring

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

    def contains(self, elem):
        """Return True if ``elem`` is an element of this module."""
        try:
            self.convert(elem)
            return True
        except CoercionFailed:
            return False

    def __contains__(self, elem):
        return self.contains(elem)

class ModuleElement(object):
    """Base class for module element wrappers."""

    #self.module - containing module
    #self.data - data

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
        raise NotImplementedError

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

class FreeModule(Module):
    """Base class for free modules."""

    #self.rank - rank of the free module
    dtype = FreeModuleElement

    def __init__(self, ring, rank):
        Module.__init__(self, ring)
        self.rank = rank

    def __repr__(self):
        return repr(self.ring) + "**" + repr(self.rank)

    def __eq__(self, other):
        if not isinstance(other, FreeModule):
            return False
        return self.ring == other.ring and self.rank == other.rank

    def convert(self, elem, M=None):
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
        else:
            raise CoercionFailed
