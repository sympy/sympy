"""Implementation of :class:`Field` class. """

from sympy.polys.domains.ring import Ring

class Field(Ring):
    """Represents a field domain. """

    has_Field = True

    def get_field(self):
        """Returns a field associated with `self`. """
        return self

    def exquo(self, a, b):
        """Exact quotient of `a` and `b`, implies `__div__`.  """
        return a / b

    def quo(self, a, b):
        """Quotient of `a` and `b`, implies `__div__`. """
        return a / b

    def rem(self, a, b):
        """Remainder of `a` and `b`, implies nothing.  """
        return self.zero

    def div(self, a, b):
        """Division of `a` and `b`, implies `__div__`. """
        return a / b, self.zero

    def gcd(self, a, b):
        """Returns GCD of `a` and `b`. """
        return self.one

    def lcm(self, a, b):
        """Returns LCM of `a` and `b`. """
        return a*b
