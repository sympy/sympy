from sympy.combinatorics import Permutation
from sympy.core import Basic
from sympy.combinatorics.permutations import perm_af_mul, \
 _new_from_array_form, perm_af_commutes_with, perm_af_invert


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

    def __new__(cls, *args, **kw_args):
        """
        The default constructor.
        """
        obj = Basic.__new__(cls, *args, **kw_args)
        obj._generators = args[0]
        obj._order = None
        obj._center = []
        obj._is_abelian = None
        obj._commutators = []
        obj._coset_repr = []
        obj._coset_repr_n = 0
        size = len(args[0][0].array_form)
        if not all(len(args[0][i].array_form)==size for i in xrange(1, len(args[0]))):
                raise ValueError("Permutation group size is not correct")
        obj._perm_size = size
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
        gens = [p.array_form for p in self.generators]
        for x in gens:
            for y in gens:
                if y <= x:
                    continue
                if not perm_af_commutes_with(x, y):
                    self._is_abelian = False
                    return self._is_abelian
        return self._is_abelian

    def order(self):
        if self._order:
            return self._order
        self._order = len(list(self.generate(af=True)))
        return self._order

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
        >>> d.coset_repr()
        set([Permutation([0, 1, 2, 3]), Permutation([0, 1, 3, 2]), Permutation([0, 2, 3, 1]), Permutation([2, 3, 1, 0])])
        """
        if (alpha != None):
            self.schreier_tree(alpha)
        return set(self._coset_repr)

    def schreier_tree(self, alpha, gen=None, ag=None):
        """
        Computes the Schreier Tree.

        A Schreier tree with root a for S is a representation of the orbit of a.
        Its rooted at a with the elements of a^{G} as its vertices, and its edges
        describe the elements of S needed to get from a to each vertex, each edge
        {i, j} in the tree with i closer to the root than j is labeled by a
        generator s \{\in} S moving from i to j.

        In the current implementation, the entire tree is not stored. We store only
        a coset representation for each point in the orbit.
        
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
        if not self._coset_repr:
            self._coset_repr = [None]*self.perm_size
        self._coset_repr_n += 1
        self._coset_repr[alpha] = gen
        for i in xrange(len(self.generators)):
            ag = self.generators[i](alpha)
            if self._coset_repr[ag] == None:
                gen *= self.generators[i]
                self.schreier_tree(ag, gen, ag)
                gen *= ~self.generators[i]

    def generate(self, af=False):
        """
        An implementation of Dimino's algorithm to generate all the elements
        of the group given the generators.

        If af = True it yields the array form of the permutations

        Reference:
        [1] The implementation of various algorithms for Permutation Groups in
        the Computer Algebra System: AXIOM, N.J. Doye, M.Sc. Thesis

        Examples:
        >>> from sympy.combinatorics.permutations import Permutation
        >>> from sympy.combinatorics.perm_groups import PermutationGroup
        >>> a = Permutation([0, 2, 1, 3])
        >>> b = Permutation([0, 2, 3, 1])
        >>> g = PermutationGroup([a, b], 4)
        >>> list(g.generate())
        [Permutation([0, 2, 1, 3]), Permutation([0, 2, 3, 1]), Permutation([0, 1, 3, 2]), Permutation([0, 3, 2, 1]), Permutation([0, 3, 1, 2])]
        """
        idn = range(self.perm_size)
        order = 0
        element_list = [idn]
        set_element_list = set([tuple(idn)])
        if af:
            yield idn
        else:
            yield _new_from_array_form(idn)
        gens = [p.array_form for p in self.generators]

        for i in xrange(len(gens)):
            # D elements of the subgroup G_i generated by gens[:i]
            D = element_list[:]
            N = [idn]
            while N:
                A = N
                N = []
                for a in A:
                    for g in gens[:i+1]:
                        ag = perm_af_mul(a, g)
                        if tuple(ag) not in set_element_list:
                            # produce G_i*g
                            for d in D:
                                order += 1
                                ap = perm_af_mul(d, ag)
                                if af:
                                    yield ap
                                else:
                                    p = _new_from_array_form(ap)
                                    yield p
                                element_list.append(ap)
                                set_element_list.add(tuple(ap))
                                N.append(ap)
        self._order = len(element_list)


    def pointwise_stabilizers(self, points, af=False):
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
        for g in self.generate(af=True):
            stabilizes = True
            for point in points:
                if g[point] != point:
                    stabilizes = False
                    break
            if stabilizes:
                if af:
                    stabs.append(g)
                else:
                    stabs.append(_new_from_array_form(g))
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
            for point in xrange(npoints):
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
        >>> a = Permutation([1, 0, 2, 4, 3])
        >>> b = Permutation([0, 1, 3, 2, 4])
        >>> g = PermutationGroup([a, b], 5)
        >>> g.center
        [Permutation([1, 0, 2, 3, 4]), Permutation([0, 1, 2, 3, 4])]
        >>> d = Permutation([1, 0, 2, 3, 4])
        >>> d * b == b * d
        True
        """
        if self._center:
            return self._center
        gens = list(self.generate(af=True))
        for x in gens:
            for y in gens:
                if not perm_af_commutes_with(x, y):
                    break
            else:
                self._center.append(_new_from_array_form(list(x)))
        return self._center

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
        if self._commutators:
            return self._commutators
        group_elems = list(self.generate())
        gv = [p.array_form for p in group_elems]
        n = len(gv)
        gv_inv = [(~p).array_form for p in group_elems]
        set_commutators = set()
        for i in range(n):
            for j in range(n):
                x = gv[i]
                xinv = gv_inv[i]
                y = gv[j]
                yinv = gv_inv[j]
                c = [x[y[xinv[k]]] for k in yinv]

                ct = tuple(c)
                if not ct in set_commutators:
                    set_commutators.add(ct)
                    self._commutators.append(_new_from_array_form(c))
        return self._commutators

