from sympy.core import Basic
from sympy.utilities.iterables import rotate_left
from sympy.polys.polytools import lcm
from sympy.matrices import Matrix, zeros

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

    _array_form = None
    _cyclic_form = None

    @property
    def array_form(self):
        """
        This is used to convert from cyclic notation to the
        canonical notation.
        Currently singleton cycles need to be written
        explicitly.

        Examples:
        >>> from sympy.combinatorics.permutations import Permutation
        >>> p = Permutation([[2,0],[3,1]])
        >>> p.array_form
        [1, 3, 0, 2]
        """
        if self._array_form is not None:
            return self._array_form
        if not isinstance(self.args[0][0], list):
            self._array_form = self.args[0]
            return self._array_form
        cycles = self.args[0]
        linear_form = []
        for cycle in cycles:
            min_element = min(cycle)
            while cycle[0] != min_element:
                cycle = rotate_left(cycle, 1)
            linear_form.append(cycle)
        linear_form.sort(key=lambda t: -t[0])
        self._array_form = list(itertools.chain(*linear_form))
        return self._array_form

    @property
    def cyclic_form(self):
        """
        This is used to convert to the cyclic notation
        from the canonical notation.

        Examples:
        >>> from sympy.combinatorics.permutations import Permutation
        >>> p = Permutation([0,3,1,2])
        >>> p.cyclic_form
        [[1, 3, 2], [0]]
        """
        if self._cyclic_form is not None:
            return self._cyclic_form
        if isinstance(self.args[0][0], list):
            self._cyclic_form = self.args[0]
            return self._cyclic_form
        linear_rep = self.args[0]
        unchecked = [True] * len(linear_rep)
        cyclic_form = []
        for i in xrange(len(linear_rep)):
            if unchecked[i]:
                cycle = []
                cycle.append(i)
                unchecked[i] = False
                j = i
                while unchecked[linear_rep[j]]:
                    j = linear_rep[j]
                    cycle.append(j)
                    unchecked[j] = False
                cyclic_form.append(cycle)
        cyclic_form.sort(key=lambda t: -t[0])
        self._cyclic_form = cyclic_form
        return self.cyclic_form

    @property
    def size(self):
        return len(self.array_form)

    def __mul__(self, other):
        """
        Routine for multiplication of permutations.

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
        if self.size != other.size:
            raise ValueError("The number of elements in the permutations \
            dont match")
        return Permutation([other.array_form[self.array_form[i]]
                            for i in range(self.size)])

    def __pow__(self, n):
        """
        Routine for finding powers of a permutation.

        Examples:
        >>> from sympy.combinatorics.permutations import Permutation
        >>> p = Permutation([2,0,3,1])
        >>> p**4
        Permutation([0, 1, 2, 3])
        """
        return reduce(lambda x, y: x*y, [self]*n)

    def __invert__(self):
        """
        Finds the invert of a permutation.

        An invert of a permutation when multiplied by it
        results in the identity permutation.

        Examples:
        >>> from sympy.combinatorics.permutations import Permutation
        >>> p = Permutation([[2,0],[3,1]])
        >>> ~p
        Permutation([2, 0, 3, 1])
        >>> p*(~p) == Permutation([0,1,2,3])
        True
        """
        self_form = self.array_form
        inv_form = [0] * self.size
        for i in xrange(self.size):
            inv_form[self_form[i]] = i
        return Permutation(inv_form)

    def atoms(self):
        """
        Returns all the elements of a permutation
        """
        return set(self.array_form)

    def unrank_nonlex(self, r):
        """
        This is a linear time unranking algorithm that does not
        respect lexicographic order [3].

        [3] See below

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
        temp = self.array_form
        n = self.size
        id_perm = [i for i in xrange(n)]
        while n > 1:
            id_perm[n-1],id_perm[r % n] = id_perm[r % n], id_perm[n-1]
            n -= 1
            r = r/n
        return Permutation(id_perm)

    def rank_nonlex(self, inv_perm = None, n = 0):
        """
        This is a linear time ranking algorithm that does not
        enforce lexicographic order [3].

        [3] Wendy Myrvold and Frank Ruskey. 2001. Ranking and unranking
            permutations in linear time. Inf. Process. Lett. 79, 6 (September 2001),
            281-284. DOI=10.1016/S0020-0190(01)00141-7

        Examples
        >>> from sympy.combinatorics.permutations import Permutation
        >>> p = Permutation([0,1,2,3])
        >>> p.rank_nonlex()
        23
        """
        temp = self.array_form
        temp_inv = inv_perm
        if temp_inv is None:
            temp_inv = (~self).array_form
        if n == 0:
            n = self.size
        if n == 1:
            return 0
        perm_form = temp[:]
        temp_inv_form = temp_inv[:]
        s = perm_form[n-1]
        perm_form[n-1], perm_form[temp_inv_form[n-1]] = \
            perm_form[temp_inv_form[n-1]], perm_form[n-1]
        temp_inv_form[s], temp_inv_form[n-1] = \
            temp_inv_form[n-1], temp_inv_form[s]
        return s + n*self.rank_nonlex(temp_inv, n - 1)

    @property
    def is_Singleton(self):
        return self.size == 1

    @property
    def is_Empty(self):
        return self.size == 0

    @property
    def is_Identity(self):
        return self.size == len(self.cyclic_form)

    @property
    def ascents(self):
        """
        Returns the positions of ascents in a permutation, ie, the location
        where p[i] < p[i+1]

        Example:
        >>> from sympy.combinatorics.permutations import Permutation
        >>> p = Permutation([4,0,1,3,2])
        >>> p.ascents
        [1, 2]
        """
        pos = []
        temp = self.array_form
        for i in xrange(self.size-1):
            if temp[i] < temp[i+1]:
                pos.append(i)
        return pos

    @property
    def descents(self):
        """
        Returns the positions of descents in a permutation, ie, the location
        where p[i] > p[i+1]

        Example:
        >>> from sympy.combinatorics.permutations import Permutation
        >>> p = Permutation([4,0,1,3,2])
        >>> p.descents
        [0, 3]
        """
        pos = []
        temp = self.array_form
        for i in xrange(len(temp)-1):
            if temp[i] > temp[i+1]:
                pos.append(i)
        return pos

    @property
    def max(self):
        """
        The maximum element moved by the permutation.

        Example:
        >>> from sympy.combinatorics.permutations import Permutation
        >>> p = Permutation([1,0,2,3,4])
        >>> p.max
        1
        """
        temp = self.array_form
        max = 0
        for i in xrange(self.size):
            if temp[i] != i and temp[i] > max:
                max = temp[i]
        return max

    @property
    def min(self):
        """
        The minimum element moved by the permutation

        Example:
        >>> from sympy.combinatorics.permutations import Permutation
        >>> p = Permutation([0,1,4,3,2])
        >>> p.min
        2
        """
        temp = self.array_form
        min = self.size
        for i in xrange(self.size):
            if temp[i] != i and temp[i] < min:
                min = temp[i]
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
        temp = self.array_form
        for i in xrange(self.size - 1):
            for j in xrange(i + 1, self.size):
                if temp[i] > temp[j]:
                    inversions += 1
        return inversions

    @property
    def signature(self):
        """
        Gives the signature of the permutation needed to place the
        elements of the permutation in canonical order.

        The signature is calculated as (-1)^<# no. of inversions>

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
        >>> (p**(p.order))
        Permutation([0, 1, 2, 3, 4, 5])
        """
        return reduce(lcm,[1]+[len(cycle) for cycle in self.cyclic_form])

    @property
    def length(self):
        """
        Returns the number of integers moved by a permutation.
        """
        length = 0
        temp = self.array_form
        for i in xrange(self.size):
            if temp[i] != i:
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
        return len(self.cyclic_form)

    def runs(self):
        """
        Returns the runs of a permutation.

        An ascending sequence in a permutation is called a run [1]

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
        temp = self.array_form
        temp_form = temp
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

    def get_precedence_matrix(self):
        """
        Gets the precedence matrix. This is used for computing the
        distance between two permutations.

        Examples:
        >>> from sympy.combinatorics.permutations import josephus
        >>> p = josephus(3,6,1)
        >>> p.get_precedence_matrix()
        [0, 0, 0, 0, 0, 0]
        [1, 0, 0, 0, 1, 0]
        [1, 1, 0, 1, 1, 1]
        [1, 1, 0, 0, 1, 0]
        [1, 0, 0, 0, 0, 0]
        [1, 1, 0, 1, 1, 0]
        """
        m = zeros(self.size)
        perm = self.array_form
        for i in xrange(m.rows):
            for j in xrange(i + 1, m.cols):
                m[perm[i], perm[j]] = 1
        return m

    def get_precedence_distance(self, other):
        """
        Computes the precedence distance between two permutations.

        Suppose p and p' represent n jobs. The precedence metric
        counts the number of times a job j is prededed by job i
        in both p and p'. This metric is commutative.

        Examples:
        >>> from sympy.combinatorics.permutations import Permutation
        >>> p = Permutation([2, 0, 4, 3, 1])
        >>> q = Permutation([3, 1, 2, 4, 0])
        >>> p.get_precedence_distance(q)
        7
        >>> q.get_precedence_distance(p)
        7
        """
        if self.size != other.size:
            raise ValueError("The permutations must be of the same size")
        self_prec_mat = self.get_precedence_matrix()
        other_prec_mat = other.get_precedence_matrix()
        n_prec = 0
        for i in xrange(self.size):
            for j in xrange(self.size):
                if i == j:
                    continue
                if self_prec_mat[i, j] * other_prec_mat[i, j] == 1:
                    n_prec += 1
        d = self.size * (self.size - 1)/2 - n_prec
        return d

    def get_adjacency_matrix(self):
        """
        Computes the adjacency matrix of a permutation.

        If job i is adjacent to job j in a permutation p
        then we set m[i, j] = 1 where m is the adjacency
        matrix of p.

        Examples:
        >>> from sympy.combinatorics.permutations import josephus
        >>> p = josephus(3,6,1)
        >>> p.get_adjacency_matrix()
        [0, 0, 0, 0, 0, 0]
        [0, 0, 0, 0, 1, 0]
        [0, 0, 0, 0, 0, 1]
        [0, 1, 0, 0, 0, 0]
        [1, 0, 0, 0, 0, 0]
        [0, 0, 0, 1, 0, 0]
        >>> from sympy.combinatorics.permutations import Permutation
        >>> q = Permutation([0, 1, 2, 3])
        >>> q.get_adjacency_matrix()
        [0, 1, 0, 0]
        [0, 0, 1, 0]
        [0, 0, 0, 1]
        [0, 0, 0, 0]
        """
        m = zeros(self.size)
        perm = self.array_form
        for i in xrange(self.size - 1):
            m[perm[i], perm[i + 1]] = 1
        return m

    def get_adjacency_distance(self, other):
        """
        Computes the adjacency distance between two permutations.

        This metric counts the number of times a pair i,j of jobs is
        adjacent in both p and p'. If n_adj is this quantity then
        the adjacency distance is n - n_adj - 1 [1]

        [1] Reeves, Colin R. Landscapes, Operators and Heuristic search, Annals
        of Operational Research, 86, pp 473-490. (1999)

        Examples:
        >>> from sympy.combinatorics.permutations import Permutation, josephus
        >>> p = Permutation([0, 3, 1, 2, 4])
        >>> q = josephus(4, 5, 2)
        >>> p.get_adjacency_distance(q)
        3
        """
        if self.size != other.size:
            raise ValueError("The permutations must be of the same size")
        self_adj_mat = self.get_adjacency_matrix()
        other_adj_mat = other.get_adjacency_matrix()
        n_adj = 0
        for i in xrange(self.size):
            for j in xrange(self.size):
                if i == j:
                    continue
                if self_adj_mat[i, j] * other_adj_mat[i, j] == 1:
                    n_adj += 1
        d = self.size - n_adj - 1
        return d

    def get_positional_distance(self, other):
        """
        Computes the positional distance between two permutations.

        Examples:
        >>> from sympy.combinatorics.permutations import Permutation, josephus
        >>> p = Permutation([0, 3, 1, 2, 4])
        >>> q = josephus(4, 5, 2)
        >>> p.get_positional_distance(q)
        12
        """
        if self.size != other.size:
            raise ValueError("The permutations must be of the same size")
        perm_array = self.array_form
        other_array = other.array_form
        return sum([abs(perm_array[i] - other_array[i]) for i in xrange(self.size)])

def josephus(m, n, s = 1):
    """
    Computes the Josephus permutation for a given number of
    prisoners, frequency of removal and desired number of
    survivors.

    Examples:
    >>> from sympy.combinatorics.permutations import josephus
    >>> josephus(3, 40, 1)
    Permutation([2, 5, 8, 11, 14, 17, 20, 23, 26, 29, 32, 35, \
    38, 1, 6, 10, 15, 19, 24, 28, 33, 37, 3, 9, 16, 22, 30, 36, \
    4, 13, 25, 34, 7, 21, 39, 18, 0, 31, 12, 27])
    """
    from collections import deque
    m -= 1
    if s <= 0:
        s = 1
    Q = deque([i for i in xrange(0, n)])
    perm = []
    while len(Q) > s:
        for dp in xrange(0, m):
            Q.append(Q.popleft())
        perm.append(Q.popleft())
    perm.extend(list(Q))
    return Permutation(perm)
