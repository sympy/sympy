from __future__ import print_function, division

from sympy.core.sympify import _sympify

from .sets import Set, FiniteSet


class UnderlyingSetOf(Set):
    """A symbolic notation for the underlying set of an algebraic
    structure.

    Parameters
    ==========

    parent : Group
        A symbolic group object referring to
    """
    def __new__(cls, parent):
        parent = _sympify(parent)
        obj = super(UnderlyingSetOf, cls).__new__(cls, parent)
        return obj

    @property
    def parent(self):
        return self.args[0]

    def _eval_rewrite_as_FiniteSet(self, *args, **kwargs):
        from sympy.combinatorics.perm_groups import PermutationGroup
        parent = self.parent

        if isinstance(parent, PermutationGroup):
            return FiniteSet(*self.parent.generate())

    def _contains(self, other, **kwargs):
        from sympy.combinatorics.perm_groups import PermutationGroup
        parent = self.parent

        if isinstance(parent, PermutationGroup):
            return parent.contains(other, **kwargs)
