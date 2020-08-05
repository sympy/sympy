from sympy.core import Tuple
from sympy.core.symbol import Str
from sympy.core.sympify import _sympify
from sympy.sets import Set, Union, ProductSet

__all__ = [
    'AlgebraicStructure',
]

class AlgebraicStructure(Set):
    r"""
    A base class for algebriac structure, i.e. group, ring, etc.

    Explanation
    ===========

    Algebraic structure is a mathematical structure which consists of one or more
    nonempty sets, and a collection of operations on these sets [1]. The operations
    must be closed to the set, i.e. $ f: \mathbb{S}^{n} \rightarrow \mathbb{S} $, where
    $f$ is the operator and $\mathbb{S}$ is the set.
    If every set of structure $B$ is subset of any set of structure $A$, and every operator
    of $B$ is restricted operator of $A$, $B$ is substructure of $A$.

    Examples
    ========

    >>> from sympy import Set, Map, AlgebraicStructure
    >>> A = Set('A')
    >>> a = A.element('a')
    >>> f = Map('f', domain=A*A, codomain=A)
    >>> S = AlgebraicStructure('S', (A,), (f,))

    >>> a in S
    True
    >>> f(a, a) in S
    True

    >>> B = Set('B', (A,))
    >>> f2 = Map('f', domain=B*B, codomain=B)
    >>> S2 = AlgebraicStructure('S2', (B,), (f2,))
    >>> S2.is_substructure(S)
    True

    Parameters
    ==========

    name : str
        Name of the structure used for printing.

    sets : tuple of Sets
        See sympy.sets module.

    operators : tuple of Maps
        See sympy.map module.

    References
    ==========

    https://en.wikipedia.org/wiki/Algebraic_structure

    """
    def __new__(cls, name, sets, operators, **kwargs):

        if not isinstance(name, Str):
            name = Str(name)

        sets = Tuple(*[_sympify(a) for a in sets])
        operators = Tuple(*[_sympify(a) for a in operators])

        if not sets:
            raise TypeError("At least one set must be provided.")

        if not operators:
            # Set is a degenerate algebraic structure with no operations.
            # If no operator is given, return the set.
            return Union(*sets)

        obj = super().__new__(cls, name, sets, operators)
        obj._check_closure()
        return obj

    @property
    def name(self):
        return self.args[0]

    @property
    def sets(self):
        return self.args[1]

    @property
    def operators(self):
        return self.args[2]

    @property
    def domain(self):
        sets = []
        for s in self.sets:
            if isinstance(s, AlgebraicStructure):
                sets.append(s.domain)
            else:
                sets.append(s)
        return Union(*sets)

    def _check_closure(self):
        for o in self.operators:
            n = o.arity

            # check domain
            if n == 1 and not isinstance(o.domain, ProductSet):
                domains = (o.domain,)
            else:
                domains = o.domain.args
            for d in domains:
                if not d.is_subset(self.domain):
                    raise TypeError(
                "%s is not closed on the structure." % o
                )

            if not o.codomain.is_subset(self.domain):
                raise TypeError(
                "%s is not closed on the structure." % o
                )

    def _contains(self, other):
        return self.domain.contains(other)

    def element(self, name):
        return self.domain.element(name)

    def is_substructure(self, other):

        if self == other:
            return True

        val = self._eval_is_substructure(other)
        if val is not None:
            return val

        if not issubclass(self.func, other.func):
            return False
        if not self.domain.is_subset(other.domain):
            return False

        if len(other.operators) < len(self.operators):
            return False
        zip_ops = zip(self.operators, other.operators[:len(self.operators)])
        for self_o, other_o in zip_ops:
            if not self_o.is_restriction(other_o):
                return False

        return True
    is_subset = is_substructure

    def _eval_is_substructure(self, other):
        return

    def is_superstructure(self, other):

        if self == other:
            return True

        val = self._eval_is_superstructure(other)
        if val is not None:
            return val

        if not issubclass(other.func, self.func):
            return False
        if not self.domain.is_superset(other.domain):
            return False

        if len(self.operators) < len(other.operators):
            return False
        zip_ops = zip(self.operators[:len(other.operators)], other.operators)
        for self_o, other_o in zip_ops:
            if not other_o.is_restriction(self_o):
                return False

        return True
    is_superset = is_superstructure

    def _eval_is_superstructure(self, other):
        return
