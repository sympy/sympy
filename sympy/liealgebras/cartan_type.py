from sympy.core import Basic, Symbol, Dict, Tuple


class CartanType_generator(Basic):
    """
    Constructor for actually creating things
    """
    def __call__(self, *args):
        c = args[0]
        c = list(c)

        letter, n = c[0], int(c[1])
        if letter == "A":
            if n >= 0:
                from . import type_a
                return type_a.TypeA(n)
        if letter == "B":
            if n >= 0:
                from . import type_b
                return type_b.TypeB(n)

        if letter == "C":
            if n >= 0:
                from . import type_C
                return type_C.CartanType(n)

        if letter == "D":
            if n >= 0:
                from . import type_D
                return type_D.CartanType(n)


        if letter == "E":
            if n >= 6 and n <= 8:
                from . import type_E
                return type_E.CartanType(n)

        if letter == "F":
            if n == 4:
                from . import type_F
                return type_F.CartanType(n)

        if letter == "G":
            if n == 2:
                from . import type_G
                return type_G.CartanType(n)

CartanType = CartanType_generator()


class Standard_Cartan(Basic):
    """
    Concrete base class for Cartan types such as A4, etc
    """

    def __init__(self, series, n):
        self.series = series
        self.n = n

    def rank(self):
        """
        Returns the rank of the Lie algebra
        """
        return self.n

    def type(self):
        """
        Returns the type of the Lie algebra
        """
        return self.series
