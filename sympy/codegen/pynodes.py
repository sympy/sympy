from sympy.core import Tuple


class List(Tuple):
    """Represents a (frozen) (python) list (for code printing purposes)."""
    def __eq__(self, other):
        if isinstance(other, list):
            return self == List(*other)
        else:
            return self.args == other
