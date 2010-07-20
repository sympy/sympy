"""Implementation of :class:`ModularInteger` class. """

class ModularInteger(object):
    """A class representing an modular integers. """

    mod, dom, sym = None, None, None

    __slots__ = ['val']

    def __init__(self, val):
        self.val = val % self.mod

    def __repr__(self):
        return "%s(%s)" % (self.__class__.__name__, self.val)

    def __str__(self):
        return "%s mod %s" % (self.val, self.mod)

    def __int__(self):
        return int(self.to_int())

    def to_int(self):
        if self.sym:
            if self.val <= self.mod // 2:
                return self.val
            else:
                return self.val - self.mod
        else:
            return self.val

    def __pos__(self):
        return self

    def __neg__(self):
        return self.__class__(-self.val)

    def __add__(self, other):
        return self.__class__(self.val + other.val)

    def __sub__(self, other):
        return self.__class__(self.val - other.val)

    def __mul__(self, other):
        return self.__class__(self.val * other.val)

    def __div__(self, other):
        return self.__class__(self.val * self.dom.invert(other.val, self.mod))

    __truediv__ = __div__

    def __mod__(self, other):
        return self.__class__(self.val % other.val)

    def __pow__(self, exp):
        if not exp:
            return self.__class__(self.dom.one)

        if exp < 0:
            val, exp = self.dom.invert(self.val, self.mod), -exp
        else:
            val = self.val

        return self.__class__(val**exp)

    def __eq__(self, other):
        return isinstance(other, ModularInteger) and self.val == other.val

    def __ne__(self, other):
        return not isinstance(other, ModularInteger) or self.val != other.val

    def __nonzero__(self):
        return bool(self.val)

def ModularIntegerFactory(_mod, _dom, _sym):
    """Create custom class for specific integer modulus."""

    class cls(ModularInteger):
        mod, dom, sym = _dom.convert(_mod), _dom, _sym

    if _sym:
        cls.__name__ = "SymmetricModularInteger%s" % _mod
    else:
        cls.__name__ = "ModularInteger%s" % _mod

    return cls

