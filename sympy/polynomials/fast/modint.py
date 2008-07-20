"""Class for modular integer arithmetic."""

from sympy.core.numbers import igcdex

class ModularInteger(object):
    """Modular integer arithmetic, based on Python integers."""

    modulus = 0

    def __init__(self, value):
        self.value = value % self.modulus

    def __repr__(self):
        return "%s(%s)" % (self.__class__.__name__, self.value)

    def __str__(self):
        return "%s mod %s" % (int(self), self.modulus)

    def __int__(self):
        """Return the unique integer -m/2 < i <= m/2."""
        if self.value <= self.modulus/2:
            return self.value
        else:
            return self.value - self.modulus

    def __pos__(self):
        return self

    def __neg__(self):
        return self.__class__(-self.value)

    def __add__(self, other):
        return self.__class__(self.value + other.value)

    def __sub__(self, other):
        return self.__class__(self.value - other.value)

    def __mul__(self, other):
        return self.__class__(self.value * other.value)

    def __div__(self, other):
        x, y, g = igcdex(other.value, self.modulus)
        assert g == 1, "Zero division!"
        return self.__class__(self.value * x)

    def __pow__(self, exponent):
        """Repeated squaring."""
        exponent = int(exponent)
        if exponent < 0:
            x, y, g = igcdex(self.value, self.exponent)
            assert g == 1, "Zero division!"
            value = x
            exponent *= -1
        if exponent == 0:
            return self.__class__(1)
        else:
            value = self.value
        value = pow(value, exponent, self.modulus)
        return self.__class__(value)

    def __eq__(self, other):
        return self.value == other.value

    def __ne__(self, other):
        return self.value != other.value

    def __nonzero__(self):
        return bool(self.value)

def ModularIntegerFactory(m):
    """Create custom class for specific integer modulus."""

    class newClass(ModularInteger):
        modulus = m

    newClass.__name__ = "IntMod%s" % m
    return newClass
