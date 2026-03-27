from __future__ import annotations
from sympy.core import Atom, Basic


class CartanType_generator():
    """
    Constructor for actually creating Cartan types
    """

    def __call__(self, *args):
        if len(args) != 1:
            raise TypeError("CartanType takes exactly one argument")

        c = args[0]

        # -------- Handle list input --------
        if isinstance(c, list):
            if len(c) != 2:
                raise ValueError(f"Invalid Cartan type '{c}'")

            letter, n = c

            if not isinstance(letter, str):
                raise TypeError("Cartan type letter must be a string")

            letter = letter.strip().upper()

            if not isinstance(n, int):
                raise TypeError("Rank must be an integer")

        # -------- Handle string input --------
        elif isinstance(c, str):
            c = c.strip()

            if len(c) < 2:
                raise ValueError(f"Invalid Cartan type '{c}'")

            letter = c[0].upper()
            rank_str = c[1:]

            if not rank_str.isdecimal():
                raise ValueError(f"Invalid Cartan type '{c}'")

            n = int(rank_str)

        else:
            raise TypeError(
                "Argument must be a string (e.g. 'A3') or a list (e.g. ['A', 3])"
            )

        # -------- Rank validation --------
        if n < 0:
            raise ValueError(
                f"Invalid Cartan type '{c}' (rank must be non-negative)"
            )

        # -------- Dispatch --------
        if letter == "A":
            if n < 1:
                raise ValueError(f"Invalid Cartan type '{c}'")
            from .type_a import TypeA
            return TypeA(n)

        if letter == "B":
            if n < 2:
                raise ValueError(f"Invalid Cartan type '{c}'")
            from .type_b import TypeB
            return TypeB(n)

        if letter == "C":
            if n < 2:
                raise ValueError(f"Invalid Cartan type '{c}'")
            from .type_c import TypeC
            return TypeC(n)

        if letter == "D":
            if n < 3:
                raise ValueError(f"Invalid Cartan type '{c}'")
            from .type_d import TypeD
            return TypeD(n)

        if letter == "E" and 6 <= n <= 8:
            from .type_e import TypeE
            return TypeE(n)

        if letter == "F" and n == 4:
            from .type_f import TypeF
            return TypeF(n)

        if letter == "G" and n == 2:
            from .type_g import TypeG
            return TypeG(n)

        raise ValueError(f"Invalid Cartan type '{c}'")


CartanType = CartanType_generator()


class Standard_Cartan(Atom):
    """
    Concrete base class for Cartan types such as A4, etc
    """

    def __new__(cls, series, n):
        if not isinstance(series, str):
            raise TypeError("Series must be a string")
        if not isinstance(n, int):
            raise TypeError("Rank must be an integer")

        obj = Basic.__new__(cls)
        obj.n = n
        obj._series = series.upper()
        return obj

    def rank(self):
        return self.n

    def series(self):
        return self._series

    # important for tests & debugging
    def __repr__(self):
        return f"{self._series}{self.n}"

    __str__ = __repr__
