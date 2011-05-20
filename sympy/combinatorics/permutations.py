from sympy.core import Basic
from sympy.utilities.iterables import rotate_left
from sympy.polys.polytools import lcm

import itertools

class Permutation(Basic):
    """
    A permutation, alternatively known as an 'arrangement number'
    or 'ordering' is an arrangement of the elements of an ordered list
    into a one-to-one mapping with  itself. The number of permutations
    on a set of n elements is given by n!.

    A representation of a permutation as a product of permutation cycles
    is unique (up to the ordering of the cycles). An example of a cyclic
    decomposition is the permutation [4, 2, 1, 3] of the set [1, 2, 3, 4].
    This is denoted as [[2], [1, 4, 3]], corresponding to the disjoint
    permutation cycles [2] and [1, 4, 3]. We can choose the cyclic form as
    we want since the cycles are disjoint and can therefore be specified
    in any order and a rotation of a given cycle specifies the same cycle [1]
    Therefore, (431)(2), (314)(2), (143)(2), (2)(431), (2)(314), and (2)(143)
    all describe the same permutation.

    Another notation that explicitly identifies the positions occupied by
    elements before and after application of a permutation on n elements uses a
    2xn matrix, where the first row is the identity permutation  and the second
    row is the new arrangement [2].

    Any permutation is also a product of transpositions.

    Permutations are commonly denoted in lexicographic or transposition order.

    [1] Skiena, S. 'Permutations.' 1.1 in Implementing Discrete Mathematics
        Combinatorics and Graph Theory with Mathematica.
        Reading, MA: Addison-Wesley, pp. 3-16, 1990.
    [2] Knuth, D. E. The Art of Computer Programming, Vol. 4: Combinatorial Algorithms,
        1st ed. Reading, MA: Addison-Wesley, 2011.
    """
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
        for i in xrange(len(mul1.args[0])):
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
        if val.is_CyclicForm:
            val = val.to_array()
        ref = val
        for i in xrange(n-1):
            val = val * ref
        if self.is_CyclicForm:
            return val.to_cycles()
        return val

    def __invert__(self):
        inv_form = [0] * len(self.args[0])
        self_form = self.args[0]
        for i in xrange(len(self.args[0])):
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
        id_perm = [i for i in xrange(n)]
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
        Computes the number of inversions of a permutation.

        An inversion is where i > j but p[i] < p[j].

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

    @property
    def signature(self):
        """
        Gives the signature of the permutation needed to place the
        elements of the permutation in canonical order.

        Example:
        >>> from sympy.combinatorics.permutations import Permutation
        >>> p = Permutation([0,1,2])
        >>> p.signature
        1
        >>> q = Permutation([0,2,1])
        >>> q.signature
        -1
        """
        return (-1)**self.inversions

    @property
    def order(self):
        """
        Computes the order of a permutation.

        When the permutation is raised to the power of its
        order it equals the identity permutation.

        Example:
        >>> from sympy.combinatorics.permutations import Permutation
        >>> p = Permutation([3,1,5,2,4,0])
        >>> p.order
        4
        >>> (p**(p.order)).is_Identity
        True
        """
        temp = self
        if temp.is_ArrayForm:
            temp = temp.to_cycles()
        order = 1
        for cycle in temp.args[0]:
            order = lcm(order, len(cycle))
        return order

    @property
    def length(self):
        """
        Returns the number of integers moved by a permutation.
        """
        length = 0
        temp = self
        if temp.is_CyclicForm:
            temp = temp.to_array()
        temp_form = temp.args[0]
        for i in xrange(len(temp_form)):
            if temp_form[i] != i:
                length += 1
        return length

    @property
    def is_Positive(self):
        return self.signature > 0

    @property
    def is_Negative(self):
        return self.signature < 0

    @property
    def cycles(self):
        temp = self
        if temp.is_ArrayForm:
            temp = temp.to_cycles()
        return len(temp.args[0])

    def runs(self):
        """
        Returns the run of a permutation.

        A set of ascending sequences in a permutation is called a run [1]

        [1] Graham, R. L.; Knuth, D. E.; and Patashnik, O.
            Concrete Mathematics: A Foundation for Computer Science, 2nd ed.
            Reading, MA: Addison-Wesley, 1994.

        Examples:
        >>> from sympy.combinatorics.permutations import Permutation
        >>> p = Permutation([2,5,7,3,6,0,1,4,8])
        >>> p.runs()
        [[2, 5, 7], [3, 6], [0, 1, 4, 8]]
        >>> q = Permutation([1,3,2,0])
        >>> q.runs()
        [[1, 3], [2], [0]]
        """
        temp = self
        if temp.is_CyclicForm:
            temp = temp.to_array()
        temp_form = temp.args[0]
        cycles = []
        temp_cycle = []
        for i in xrange(len(temp_form) - 1):
            current_elem = temp_form[i]
            next_elem    = temp_form[i+1]

            if current_elem < next_elem:
                temp_cycle.append(current_elem)
                continue

            if current_elem > next_elem:
                if temp_cycle != [] and \
                       temp_cycle[-1] < current_elem:
                    temp_cycle.append(current_elem)
                    cycles.append(temp_cycle)
                    temp_cycle = []
                    continue
                else:
                    if temp_cycle != []:
                        cycles.append(temp_cycle)
                    cycles.append([current_elem])
                    temp_cycle = []
                    continue

        if current_elem < next_elem:
            temp_cycle.append(next_elem)
            cycles.append(temp_cycle)
        else:
            if temp_cycle != []:
                cycles.append(temp_cycle)
            cycles.append([next_elem])
        return cycles
