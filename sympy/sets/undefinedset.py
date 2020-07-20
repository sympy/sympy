from sympy.core import Tuple
from sympy.core.symbol import Str
from sympy.core.sympify import _sympify
from .sets import Set

class UndefinedSet(Set):
    """
    An abstract set with no condition defined, which is commonly
    encountered in set theory and abstract algebra.

    Parameters
    ==========

    name : str
        Name of the set.

    supersets : tuple of Sets, optional
        List of supersets of the set.

    Examples
    ========

    >>> from sympy import Set
    >>> A = Set('A')
    >>> B = Set('B', (A,))

    >>> A
    A
    >>> B
    B

    >>> A.is_superset(B)
    True
    >>> B.is_subset(A)
    True

    """
    def __new__(cls, name, supersets=()):
        if not isinstance(name, Str):
            name = Str(name)
        supersets = Tuple(*[_sympify(s) for s in supersets])

        return super().__new__(cls, name, supersets)

    @property
    def name(self):
        return self.args[0]

    @property
    def supersets(self):
        return self.args[1]

    def _eval_is_subset(self, other):
        for s in self.supersets:
            if s.is_subset(other):
                return True
        return False
