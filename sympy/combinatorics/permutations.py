from sympy.core import Basic
from sympy.utilities.iterables import rotate_left

import itertools

class Permutation(Basic):
    is_Permutation = True

    def __mul__(self, other):
        """
        Routine for multiplication of permutations

        Examples:
        >>> from sympy.combinatorics.permutations import Permutation
        >>> p = Permutation([0,1,2,3])
        >>> q = Permutation([3,2,1,0])
        >>> p*q
        Permutation([3, 2, 1, 0])

        If one of the permutations is in a cyclic form then it is first
        converted to an array form and then multiplied.
        >>> from sympy.combinatorics.permutations import Permutation
        >>> q = Permutation([[1,3,2],[0]])
        >>> p = Permutation([0,3,1,2])
        >>> p*q
        Permutation([1, 0, 3, 2])
        """
        mul1, mul2 = self, other
        if not mul1.is_ArrayForm:
            mul1 = mul1.to_array()
        if not mul2.is_ArrayForm:
            mul2 = mul2.to_array()
        if len(mul1.args[0]) != len(mul2.args[0]):
            raise ValueError("The permutations must have equal \
            number of elements")
        return_val = [None] * len(mul1.args[0])
        mul1_form = mul1.args[0]
        mul2_form = mul2.args[0]
        for i in range(len(mul1.args[0])):
             return_val[i] = mul2_form[mul1_form[i]]
        return Permutation(return_val)


    def __pow__(self, n):
        """
        Routine for finding powers of a permutation.

        Examples:
        >>> from sympy.combinatorics.permutations import Permutation
        >>> p = Permutation([2,0,3,1])
        >>> p**4
        Permutation([0, 1, 2, 3])
        """
        val = self
        for i in range(n-1):
            val = self*val
        if self.is_CyclicForm:
            return val.to_cycles()
        return val

    def __invert__(self):
        inv_form = [0] * len(self.args[0])
        self_form = self.args[0]
        for i in range(len(self.args[0])):
            inv_form[self_form[i]] = i
        return Permutation(inv_form)

    def to_array(self):
        """
        This is used to convert from cyclic notation to the
        canonical notation.
        Currently singleton cycles need to be written
        explicitly.

        Examples:
        >>> from sympy.combinatorics.permutations import Permutation
        >>> p = Permutation([[2,0],[3,1]])
        >>> p.to_array()
        Permutation([1, 3, 0, 2])
        """
        if self.is_ArrayForm:
            return
        cycles = self.args[0]
        linear_form = []
        for cycle in cycles:
            min_element = min(cycle)
            while cycle[0] != min_element:
                cycle = rotate_left(cycle, 1)
            linear_form.append(cycle)
        linear_form.sort(key=lambda t: -t[0])
        return Permutation(list(itertools.chain(*linear_form)))

    def to_cycles(self):
        """
        This is used to convert to the cyclic notation
        from the canonical notation.

        Examples:
        >>> from sympy.combinatorics.permutations import Permutation
        >>> p = Permutation([0,3,1,2])
        >>> p.to_cycles()
        Permutation([[1, 3, 2], [0]])
        """
        if self.is_CyclicForm:
            return
        linear_rep = self.args[0]
        P = [True] * len(linear_rep)
        cyclic_form = []
        for i in xrange(len(linear_rep)):
            if P[i]:
                cycle = []
                cycle.append(i)
                P[i] = False
                j = i
                while P[linear_rep[j]]:
                    j = linear_rep[j]
                    cycle.append(j)
                    P[j] = False
                cyclic_form.append(cycle)
        cyclic_form.sort(key=lambda t: -t[0])
        return Permutation(cyclic_form)

    @property
    def is_ArrayForm(self):
        return not isinstance(self.args[0][0], list)

    @property
    def is_CyclicForm(self):
        return isinstance(self.args[0][0], list)

    def atoms(self):
        """
        Returns all the elements of a permutation
        """
        if self.is_ArrayForm:
            return set(self.args[0])
        return self.to_array().atoms()

    def unrank_nonlex(self, r):
        """
        This is a linear time unranking algorithm that does not
        respect lexicographic order [1].

        [1] Wendy Myrvold and Frank Ruskey. 2001. Ranking and unranking \
            permutations in linear time. Inf. Process. Lett. 79, 6 (September 2001),\
            281-284. DOI=10.1016/S0020-0190(01)00141-7

        Examples:
        >>> from sympy.combinatorics.permutations import Permutation
        >>> p = Permutation([0,1,2,3])
        >>> p.unrank_nonlex(5)
        Permutation([2, 0, 3, 1])

        Consider in the cyclic form
        >>> from sympy.combinatorics.permutations import Permutation
        >>> p = Permutation([[0,1],[2,3]])
        >>> p.unrank_nonlex(5)
        Permutation([2, 0, 3, 1])
        """
        temp = self
        if temp.is_CyclicForm:
            temp = self.to_array()
        n = len(temp.args[0])
        id_perm = [i for i in range(n)]
        while n > 1:
            id_perm[n-1],id_perm[r % n] = id_perm[r % n], id_perm[n-1]
            n -= 1
            r = r/n
        return Permutation(id_perm)

    def rank_nonlex(self, inv_perm = None, n = 0):
        """
        This is a linear time ranking algorithm that does not
        enforce lexicographic order [1].

        Examples
        >>> from sympy.combinatorics.permutations import Permutation
        >>> p = Permutation([2,0,3,1])
        >>> p.rank_nonlex()
        5
        """
        temp = self
        if temp.is_CyclicForm:
            temp = temp.to_array()
        temp_inv = inv_perm
        if temp_inv is None:
            temp_inv = ~temp
            if temp_inv.is_CyclicForm:
                temp_inv = temp_inv.to_array()
        if n == 0:
            n = len(temp.args[0])
        if n == 1:
            return 0
        perm_form = temp.args[0]
        temp_inv_form = temp_inv.args[0]
        s = perm_form[n-1]
        perm_form[n-1], perm_form[temp_inv_form[n-1]] = \
                        perm_form[temp_inv_form[n-1]], perm_form[n-1]
        temp_inv_form[s], temp_inv_form[n-1] = \
                          temp_inv_form[n-1], temp_inv_form[s]
        return s + n*temp.rank_nonlex(temp_inv, n - 1)

    @property
    def is_Singleton(self):
        return len(self.args[0]) == 1

    @property
    def is_Empty(self):
        return len(self.args[0]) == 0

    @property
    def is_Identity(self):
        if self.is_CyclicForm:
            return len(self.atoms()) == \
                   len(self.args[0])
        return len(self.atoms()) == \
               len(self.to_cycles().args[0])

    @property
    def ascents(self):
        """
        Returns the positions of ascents in a permutation, ie, the location
        where p_{i} < p_{i+1}
        """
        pos = []
        temp = self
        if temp.is_CyclicForm:
            temp = temp.to_array()
        for i in xrange(len(temp.args[0])-1):
            if temp.args[0][i] < temp.args[0][i+1]:
                pos.append(i)
        return pos

    @property
    def descents(self):
        """
        Returns the positions of descents in a permutation, ie, the location
        where p_{i} > p_{i+1}
        """
        pos = []
        temp = self
        if temp.is_CyclicForm:
            temp = temp.to_array()
        for i in xrange(len(temp.args[0])-1):
            if temp.args[0][i] > temp.args[0][i+1]:
                pos.append(i)
        return pos

    @property
    def max(self):
        """
        The maximum element moved by the permutation.
        """
        temp = self
        if temp.is_CyclicForm:
            temp = temp.to_array()
        temp_form = temp.args[0]
        max = 0
        for i in xrange(len(temp_form)):
            if temp_form[i] != i and temp_form[i] > max:
                max = temp_form[i]
        return max

    @property
    def min(self):
        """
        The minimum element moved by the permutation
        """
        temp = self
        if temp.is_CyclicForm:
            temp = temp.to_array()
        temp_form = temp.args[0]
        min = len(temp_form)
        for i in xrange(len(temp_form)):
            if temp_form[i] != i and temp_form[i] < min:
                min = temp_form[i]
        return min

    @property
    def inversions(self):
        """
        Computes the inversions of a permutation.

        An inversion is where i > j but p_{i} < p_{j}.

        Examples:
        >>> from sympy.combinatorics.permutations import Permutation
        >>> p = Permutation([0,1,2,3,4,5])
        >>> p.inversions
        0

        In the cyclic form
        >>> from sympy.combinatorics.permutations import Permutation
        >>> p = Permutation([[4,0,2],[3,5],[1]])
        >>> p.inversions
        8
        """
        inversions = 0
        temp = self
        if temp.is_CyclicForm:
            temp = temp.to_array()
        temp_form = temp.args[0]
        for i in xrange(len(temp_form)):
            for j in xrange(len(temp_form)):
                if i==j:
                    continue
                if i < j and temp_form[i] > temp_form[j]:
                    inversions += 1
        return inversions

