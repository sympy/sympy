"""Implementation of :class:`ModularInteger` class. """

class ModularInteger(object):
    """A class representing an modular integers. """

    mod, dom = None, None

    __slots__ = ['val']

    def __init__(self, val):
        self.val = val % self.mod

    def __repr__(self):
        return "%s(%s)" % (self.__class__.__name__, self.val)

    def __str__(self):
        return "%s mod %s" % (self.val, self.mod)

    def __int__(self):
        """Return the unique integer `-m/2 < i <= m/2`. """
        if self.val <= self.mod // 2:
            return self.val
        else:
            return self.val - self.mod

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

    def __pow__(self, exp):
        if not exp:
            return self.__class__(self.dom.one)

        if exp < 0:
            val, exp = self.dom.invert(self.val, self.mod), -exp
        else:
            val = self.val

        return self.__class__(val**exp)

    def __eq__(self, other):
        return self.val == other.val

    def __ne__(self, other):
        return self.val != other.val

    def __nonzero__(self):
        return bool(self.val)

def ModularIntegerFactory(_mod, _dom):
    """Create custom class for specific integer modulus."""

    class cls(ModularInteger):
        mod, dom = _dom.convert(_mod), _dom

    cls.__name__ = "ModularInteger%s" % _mod

    return cls
