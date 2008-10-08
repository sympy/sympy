
from sympy.core.basic import Basic
from sympy.core.function import Function, diff


class Piecewise(Function):
    """
    Represents a piecewise function.

    Usage
    =====
      Piecewise(x, (-1, 0, f(x)), (0, oo, g(x))) -> Returns piecewise function
        - The first argument is the variable of the intervals.
        - The subsequent arguments are tuples defining each piece
          (begin, end, function)

    Examples
    ========
      >>> from sympy import *
      >>> x = Symbol('x')
      >>> f = x**2
      >>> g = log(x)
      >>> p = Piecewise(x, (-1,0,f), (0,oo,g))
      >>> p.diff(x)
      Piecewise(x, (-1, 0, 2*x), (0, oo, 1/x))
      >>> f*p
      x**2*Piecewise(x, (-1, 0, x**2), (0, oo, log(x)))
    """

    nargs=1

    @classmethod
    def canonize(cls, *args):
        if not args[0].is_Symbol:
            raise TypeError, "First argument must be symbol"
        for piece in args[1:]:
            if not isinstance(piece,tuple) or len(piece) != 3:
                raise TypeError, "Must use 3-tuples for intervals"
        return None

    def _eval_derivative(self, s):
        new_pieces = []
        for start, end, f in self.args[1:]:
            t = (start, end, diff(f, s))
            new_pieces.append( t )
        return Piecewise(self.args[0], *new_pieces)

    def _eval_subs(self, old, new):
        if self == old:
            return new
        new_pieces = []
        for start, end, f in self.args[1:]:
            new_pieces.append( (start, end, f.subs(old,new)) )
        return Piecewise(self.args[0], *new_pieces)
