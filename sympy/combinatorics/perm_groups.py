from sympy.combinatorics import Permutation
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
    _is_abelian = None

    def __new__(cls, *args, **kw_args):
        """
        The default constructor.
        """
        obj = Basic.__new__(cls, *args, **kw_args)
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

    @property
    def is_abelian(self):
        """
        Checks if the generators are Abelian or not.

        Examples:
        >>> from sympy.combinatorics.permutations import Permutation
        >>> from sympy.combinatorics.perm_groups import PermutationGroup
        >>> a = Permutation([0,1])
        >>> b = Permutation([1,0])
        >>> c = PermutationGroup([a,b], 2)
        >>> c.is_abelian
        True
        """
        if self._is_abelian is not None:
            return self._is_abelian

        self._is_abelian = True
        for i in self.generators:
            for j in self.generators:
                if i == j:
                    continue
                if i * j != j * i:
                    self._is_abelian = False
                    break
        return self._is_abelian

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

    def generate(self, n = 1024):
        """
        An implementation of Dimino's algorithm to generate all the elements
        of the group given the generators.

        Reference:
        [1] The implementation of various algorithms for Permutation Groups in
        the Computer Algebra System: AXIOM, NJ Doye, M.Sc. Thesis

        Examples:
        >>> from sympy.combinatorics.permutations import Permutation
        >>> from sympy.combinatorics.perm_groups import PermutationGroup
        >>> a = Permutation([0,2,1,3])
        >>> b = Permutation([0,2,3,1])
        >>> c = PermutationGroup([a,b], 4)
        >>> list(c.generate())
        [Permutation([0, 2, 1, 3]), Permutation([0, 2, 3, 1]), Permutation([0,\
        3, 2, 1]), Permutation([0, 1, 3, 2]), Permutation([0, 3, 1, 2])]
        """
        id = Permutation(list(xrange(0, self.perm_size)))
        order = 0
        element_list = [id]

        for i in xrange(len(self.generators)):
            D = element_list[:]
            N = [Permutation(list(xrange(0, self.perm_size)))]
            while len(N) > 0:
                A = N[:]
                N = []
                for  a in A:
                    for g in self.generators[:i+1]:
                        ag   = a * g
                        if ag not in element_list:
                            for d in D:
                                order += 1
                                yield d * ag
                                element_list.append(d * ag)
                                if order >= n:
                                    yield d * ag
                                    return
                                N.append(ag)

    def pointwise_stabilizers(self, points):
        """
        Returns the point wise stabilizers under action of the group.

        Examples:
        >>> from sympy.combinatorics.permutations import Permutation
        >>> from sympy.combinatorics.perm_groups import PermutationGroup
        >>> a = Permutation([0,2,1,3])
        >>> b = Permutation([2,0,1,3])
        >>> c = PermutationGroup([a,b], 4)
        >>> c.pointwise_stabilizers([2])
        [Permutation([1, 0, 2, 3])]
        """
        stabs = []
        for g in self.generate():
            stabilizes = True
            for point in points:
                if g.array_form[point] != point:
                    stabilizes = False
                    break
            if stabilizes:
                stabs.append(g)
        return stabs

    def setwise_stabilizers(self, point_set, npoints):
        """
        Returns setwise stabilizer elements under action of group elements.

        Examples:
        >>> from sympy.combinatorics.permutations import Permutation
        >>> from sympy.combinatorics.perm_groups import PermutationGroup
        >>> a = Permutation([0,2,1,3])
        >>> b = Permutation([2,0,1,3])
        >>> c = PermutationGroup([a,b], 4)
        >>> c.setwise_stabilizers([2], 1)
        [Permutation([2, 0, 1, 3]), Permutation([2, 1, 0, 3])]
        """
        stabs = []
        for g in self.generate():
            stabilizes = True
            for point in xrange(0, npoints):
                if g.array_form[point] not in point_set:
                    stabilizes = False
                    break
            if stabilizes:
                stabs.append(g)
        return stabs
