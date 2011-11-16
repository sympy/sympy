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
    _center = []
    _commutators = []

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
        >>> from sympy.combinatorics.permutations import Permutation
        >>> from sympy.combinatorics.perm_groups import PermutationGroup
        >>> a = Permutation([0,1])
        >>> b = Permutation([1,0])
        >>> c = PermutationGroup([a,b], 2)
        >>> c.perm_size
        2
        """
        return self._perm_size

    @property
    def generators(self):
        """
        Returns the generators of the group.

        Examples:
        >>> from sympy.combinatorics.permutations import Permutation
        >>> from sympy.combinatorics.perm_groups import PermutationGroup
        >>> a = Permutation([0,1])
        >>> b = Permutation([1,0])
        >>> c = PermutationGroup([a,b], 2)
        >>> c.generators
        [Permutation([0, 1]), Permutation([1, 0])]
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

    def coset_repr(self, alpha=None):
        """
        Computes the Schreier Tree for a given alpha and returns the
        corresponding coset representation.

        Examples:
        >>> from sympy.combinatorics.perm_groups import PermutationGroup
        >>> from sympy.combinatorics.permutations import Permutation
        >>> a = Permutation([0,2,1,3])
        >>> b = Permutation([0,2,3,1])
        >>> c = Permutation([2,0,3,1])
        >>> d = PermutationGroup([a,b,c],4)
        >>> d.schreier_tree(3)
        >>> d.coset_repr
        set[Permutation([2, 3, 1, 0]), Permutation([0, 2, 3, 1]), \
        Permutation([0, 1, 3, 2]), Permutation([0, 1, 2, 3])]
        """
        if (alpha != None):
            self.schreier_tree(alpha)
        return set(self._coset_repr)

    def schreier_tree(self, alpha, gen=None, ag=None):
        """
        Computes the Schreier Tree.

        Examples:
        >>> from sympy.combinatorics.perm_groups import PermutationGroup
        >>> from sympy.combinatorics.permutations import Permutation
        >>> a = Permutation([0,2,1,3,4])
        >>> b = Permutation([0,2,3,1,4])
        >>> c = Permutation([2,4,3,1,0])
        >>> d = Permutation([0,2,4,3,1])
        >>> e = PermutationGroup([a,b,c,d],5)
        >>> e.schreier_tree(3)
        >>> e._coset_repr
        [Permutation([2, 3, 1, 0, 4]), Permutation([0, 2, 3, 1, 4]), \
        Permutation([0, 1, 3, 2, 4]), Permutation([0, 1, 2, 3, 4]), \
        Permutation([0, 2, 3, 4, 1])]
        """
        if gen == None:
            gen = Permutation([i for i in xrange(self.perm_size)])
        if self._coset_repr == []:
            self._coset_repr = [None] * self.perm_size
        self._coset_repr_n += 1
        self._coset_repr[alpha] = gen
        for i in xrange(len(self.generators)):
            ag = self.generators[i](alpha)
            if self._coset_repr[ag] == None:
                gen *= self.generators[i]
                self.schreier_tree(ag, gen, ag)
                gen *= ~self.generators[i]

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

    @property
    def center(self):
        """
        Returns the central subgroup of the group.

        TODO: This is a very naive implementation and is quadratic
        in order. Superior algorithms must exist.

        Examples:
        >>> from sympy.combinatorics.perm_groups import PermutationGroup
        >>> from sympy.combinatorics.permutations import Permutation
        >>> a = Permutation([0,1,2,3])
        >>> b = Permutation([2,1,3,0])
        >>> c = PermutationGroup([a,b], 4)
        >>> c.center
        [Permutation([2, 1, 3, 0]), Permutation([3, 1, 0, 2])]
        >>> d = Permutation([2, 1, 0, 3])
        >>> d * a == a * d
        True
        """
        if self._center != []:
            return self._center
        group_elems = list(self.generate())
        for x in group_elems:
            commute = True
            for y in group_elems:
                xy = x * y
                yx = y * x
                if xy != yx:
                    commute = False
                    break
            if commute:
                self._center.append(x)
        return self._center

    @property
    def commutators(self):
        """
        Gets the commutators of the group.

        Examples:
        >>> from sympy.combinatorics.perm_groups import PermutationGroup
        >>> from sympy.combinatorics.permutations import Permutation
        >>> a = Permutation([1,0,2,3])
        >>> b = Permutation([1,0,3,2])
        >>> c = PermutationGroup([a,b], 4)
        >>> c.commutators
        [Permutation([0, 1, 2, 3])]
        """
        if self._commutators != []:
            return self._commutators
        group_elems = list(self.generate())
        for x in group_elems:
            for y in group_elems:
                c = x * y * ~x * ~y
                if not c in self._commutators:
                    self._commutators.append(c)
        return self._commutators
