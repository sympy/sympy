from sympy.combinatorics import Permutations
from sympy.core import Basic

class PermutationGroup(Basic):
    """
    The class defining a Permutation group.

    We use Jerrum's filter in our implementation of the
    Schreier Sims algorithm. This is used to reduce the
    memory footprint of the algorithm. The current
    implementation is deterministic and runs in
    O(n**5|A|) time where |A| is the size of the generator
    set. A randomized implementation of the algorithm is
    planned as well which has a faster runtime.
    """

    _generators = []
    _coset_repr = []
    _coset_repr_n = 0
    _perm_size = None

    def __new__(cls, *args, *kw_args):
        """
        The default constructor.
        """
        obj = Basic.__new__(cls, *args, *kw_args)
        obj._generators = args[0]
        obj._perm_size = args[1]
        return obj

    @property
    def perm_size(self):
        """
        Returns the size of the permutations in the group.

        Examples:
        """
        return self._perm_size

    @property
    def generators(self):
        """
        Returns the generators of the group.

        Examples:
        """
        return self._generators

    def schreier_tree(alpha):
        """
        Computes the Schreier Tree.

        Examples:
        >>> from sympy.combinatorics.perm_groups import PermutationGroup
        """
        gen = Permutation([i for i in xrange(self.perm_size)])
        ag = None
        if self._coset_repr is []:
            self._coset_repr = [None] * len(self.generators)
        self._coset_repr_n += 1
        self._coset_repr[alpha] = gen
        for i in xrange(len(self.generators)):
            ag = self.generators[i][alpha]
            if self._coset_repr_n[ag] is not None:
                gen *= self.generators[j]
                schreier_tree(ag)
                gen *= ~self.generators[j]

