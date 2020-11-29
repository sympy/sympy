from sympy.core import Basic
from .type_a import TypeA
from .type_b import TypeB
from .type_c import TypeC
from .type_d import TypeD
from .type_e import TypeE
from .type_f import TypeF
from .type_g import TypeG


class CartanType_generator(Basic):
    """
    Constructor for initializing any of the lie algebra classes by way of passing
    string args.

    Examples
    ========

    >>> from sympy.liealgebras import CartanType # CartanType = CartanType_generator()
    >>> CartanType("A4")
    TypeA('A', 4)
    """
    def __call__(self, *args):
        c = args[0]
        if type(c) == list:
            letter, n = c[0], int(c[1])
        elif type(c) == str:
            letter, n = c[0], int(c[1:])

        else:
            raise TypeError("Argument must be a string (e.g. 'A3') or a list (e.g. ['A', 3])")

        if n < 0:
            raise ValueError("Lie algebra rank cannot be negative")
        if letter == "A":
            return TypeA(n)
        if letter == "B":
            return TypeB(n)
        if letter == "C":
            return TypeC(n)
        if letter == "D":
            return TypeD(n)
        if letter == "E":
            if n >= 6 and n <= 8:
                return TypeE(n)
        if letter == "F":
            if n == 4:
                return TypeF(n)
        if letter == "G":
            if n == 2:
                return TypeG(n)

        raise TypeError("Undefined Lie algebra")

CartanType = CartanType_generator()
