from sympy.core import Basic, Tuple
from sympy.core.symbol import Str
from sympy.core.sympify import _sympify
from .sets import Set

__all__ = [
    'UndefinedSet', 'SetElement'
]

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

    def _contains(self, other):
        if hasattr(other, 'set'):
            s = other.set
            if s.is_subset(self):
                return True
        return False

class SetElement(Basic):
    """
    An arbitrary element of a set.

    .. note::
       SetElement of certain sets (e.g. `S.Reals`) returns `Symbol`,
       since `Symbol('x', real=True)` is indeed an arbitrary element of `S.Reals`.
       Operation between SetElements is not defined. To define it,
       construct algebraic structure using sympy.algebras module.

    Parameters
    ==========

    name : str
        Name of the element.

    s : Set
        Set of which the element belongs to.

    Examples
    ========

    >>> from sympy import Set
    >>> A = Set('A')
    >>> B = Set('B', (A,))
    >>> b = B.element('b')

    >>> b
    b
    >>> b in B
    True
    >>> b in A
    True

    """

    def __new__(cls, name, s, **kwargs):
        if not isinstance(name, Str):
            name = Str(name)
        s = _sympify(s)

        v = s._element(name.name)
        if v is not None:
            return v

        return super().__new__(cls, name, s)

    @property
    def name(self):
        return self.args[0]

    @property
    def set(self):
        return self.args[1]
