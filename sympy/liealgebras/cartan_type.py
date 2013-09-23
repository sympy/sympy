from __future__ import print_function, division

from sympy.core import Basic, Symbol, Dict, Tuple


class CartanType_generator(Basic):
    """
    Constructor for actually creating things
    """
    def __call__(self, *args):
        c = args[0]
        c = list(c)

        letter, n = c[0], int(c[1])
        if n < 0:
            raise ValueError("Lie algebra rank cannot be negative")
        if letter == "A":
            if not (n >= 1):
                raise ValueError("The rank of the A series must be greater than or equal to 1.")
            from . import type_a
            return type_a.TypeA(n)
        if letter == "B":
            if not (n >= 2):
                raise ValueError("The rank of the B series must be greater than or equal to 2.")
            from . import type_b
            return type_b.TypeB(n)

        if letter == "C":
            if not (n >= 3):
                raise ValueError("The rank of the C series must be greater than or equal to 3.")
            from . import type_c
            return type_c.TypeC(n)

        if letter == "D":
            if not (n >= 4):
                raise ValueError("The rank of the D series must be greater than or equal to 4.")
            from . import type_d
            return type_d.TypeD(n)

        if letter == "E":
            if not (6 <= n <= 8):
                raise ValueError("The rank of a Lie algebra of type E must be 6, 7, or 8.")
            from . import type_e
            return type_e.TypeE(n)

        if letter == "F":
            if not (n == 4):
                raise ValueError("F4 is the only Lie algebra of type F.")
            from . import type_f
            return type_f.TypeF(n)

        if letter == "G":
            if not (n == 2):
                raise ValueError("G2 is the only Lie algebra of type G.")
            from . import type_g
            return type_g.TypeG(n)

CartanType = CartanType_generator()


class Standard_Cartan(Basic):
    """
    Concrete base class for Cartan types such as A4, etc
    """

    def __new__(cls, series, n):
        obj = Basic.__new__(cls, series, n)
        obj.n = n
        obj.series = series
        return obj

    def rank(self):
        """
        Returns the rank of the Lie algebra
        """
        return self.n

    def series(self):
        """
        Returns the type of the Lie algebra
        """
        return self.series
