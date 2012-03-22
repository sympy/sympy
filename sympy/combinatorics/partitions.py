from sympy.core import Basic, C
from sympy.matrices import zeros
from sympy.functions import floor

import random

class Partition(C.FiniteSet):
    """
    This class represents an abstract partition.

    A partition is a set of disjoint sets whose union equals a given set.
    """

    _rank = None
    _super = None
    _super_list = None

    def __new__(cls, *args, **kw_args):
        """
        Generates a new partition object.

        This method also verifies if the arguments passed are
        valid and if it is found that they are not then an exception is raised.

        Examples:
        >>> from sympy.combinatorics.partitions import Partition
        >>> a = Partition([[1,2],[3]])
        >>> str(a)
        '{{3}, {1, 2}}'
        """
        partition = args[0]
        super_set = C.FiniteSet(sum(partition, []))

        check = []
        for part in partition:
            if not isinstance(part, list):
                raise ValueError("The input has been provided incorrectly")
            check.extend(part)

        check.sort()
        
        if C.FiniteSet(check) != super_set:
            raise ValueError("The partition provided is not valid.")
        obj = C.FiniteSet.__new__(cls, map(lambda x: C.FiniteSet(x), partition, **kw_args))
        obj.partition_list_form = sorted([sorted(i) for i in partition])
        return obj

    def next(self):
        """
        Generates the next partition.

        Examples:
        >>> from sympy.combinatorics.partitions import Partition
        >>> a = Partition([[1,2],[3,4,5]])
        >>> a.next()
        {{1, 2}, {3, 4}, {5}}
        """
        return self + 1

    def previous(self):
        """
        Generates the previous partition.

        Examples:
        >>> from sympy.combinatorics.partitions import Partition
        >>> a = Partition([[1,2],[3,4],[5]])
        >>> a.previous()
        {{1, 2}, {3, 4, 5}}
        """
        return self - 1

    @property
    def size(self):
        """
        Gets the size of the partition.

        Examples:
        >>> from sympy.combinatorics.partitions import Partition
        >>> a = Partition([[1,2],[3]])
        >>> a.size
        2
        """
        return len(self)

    @property
    def partition(self):
        """
        Gets the partition itself.

        Examples:
        >>> from sympy.combinatorics.partitions import Partition
        >>> a = Partition([[1,2],[3]])
        >>> a.partition
        {{3}, {1, 2}}
        """
        return self

    @property
    def partition_set(self):
        """
        Gets the set of the partition.

        Examples:
        >>> from sympy.combinatorics.partitions import Partition
        >>> a = Partition([[1,2],[3]])
        >>> a.partition_set
        {1, 2, 3}
        """
        if self._super != None:
            return self._super
        e = C.FiniteSet()
        for elem in self:
            e = e.union(elem)
        self._super = e
        return self._super

    @property
    def partition_set_size(self):
        """
        Gets the total number of elements in the partition set.

        Examples:
        >>> from sympy.combinatorics.partitions import Partition
        >>> a = Partition([[1,2],[3]])
        >>> a.partition_set_size
        3
        """
        return len(self.partition_set)

    @property
    def partition_set_list(self):
        """
        Gets the partition set as a list.

        Examples:
        >>> from sympy.combinatorics.partitions import Partition
        >>> a = Partition([[1,2],[3]])
        >>> a.partition_set_list
        [1, 2, 3]
        """
        if self._super_list != None:
            return self._super_list
        self._super_list = sorted(list(self.partition_set))
        return self._super_list

    def __repr__(self):
        return str(self.elements)

    def __add__(self, other):
        """
        Routine to add partitions.

        Examples:
        >>> from sympy.combinatorics.partitions import Partition
        >>> a = Partition([[1,2],[3]])
        >>> a.rank
        1
        >>> a = a + 3
        >>> a
        {{1}, {2}, {3}}
        >>> a.rank
        4
        >>> a = a + a
        >>> a
        {{2, 3}, {1}}
        """
        return self._partition_op(other)

    def __sub__(self, other):
        """
        Routine to add partitions.

        Examples:
        >>> from sympy.combinatorics.partitions import Partition
        >>> a = Partition([[1,2],[3]])
        >>> a.rank
        1
        >>> a = a - 1
        >>> a
        {{1, 2, 3}}
        >>> a.rank
        0
        >>> a = a + Partition([[1,2],[3]])
        >>> str(a)
        '{{1, 2}, {3}}'
        """
        return self._partition_op(other, 1)

    def _partition_op(self, other, op=0):
        """
        Helper method for the __add__ and __sub__ methods.
        """
        if isinstance(other, Partition):
            if other.partition_set != self.partition_set:
                raise ValueError("Partition sets are not equal.")
            if op == 0:
                offset = self.rank + other.rank
            else:
                offset = self.rank - other.rank
            result = RGS_unrank((offset) %
                                RGS_enum(self.partition_set_size),
                                self.partition_set_size)
        elif isinstance(other, int):
            if op == 0:
                offset = self.rank + other
            else:
                offset = self.rank - other
            result = RGS_unrank((offset) %
                                RGS_enum(self.partition_set_size),
                                self.partition_set_size)
        return Partition.partition_from_rgs(result, self.partition_set)


    def _compare(self, other):
        """
        Compares two partitions.

        The basis for comparison of two partitions is rank.
        Partitions are sorted in an ascending order wrt their ranks.

        Examples:
        >>> from sympy.combinatorics.partitions import Partition
        >>> a = Partition([[1,2],[3,4,5]])
        >>> b = Partition([[1],[2,3],[4],[5]])
        >>> a.compare(b)
        -1
        >>> a.compare(a)
        0
        >>> b.compare(a)
        1
        """
        if self.rank > other.rank:
            return -1
        elif self.rank < other.rank:
            return 1
        return 0

    def __eq__(self, other):
        """
        Checks for equality of two partitions.

        Examples:
        >>> from sympy.combinatorics.partitions import Partition
        >>> a = Partition([[1,2],[3,4,5]])
        >>> b = Partition([[1],[2,3],[4],[5]])
        >>> a == b
        False
        >>> a == a
        True
        """
        return self._compare(other) == 0

    def __ne__(self, other):
        """
        Checks for inequality of two partitions.

        Examples:
        >>> from sympy.combinatorics.partitions import Partition
        >>> a = Partition([[1,2],[3,4,5]])
        >>> b = Partition([[1],[2,3],[4],[5]])
        >>> a != b
        True
        >>> a != a
        False
        """
        return self._compare(other) != 0

    def __gt__(self, other):
        """
        Checks if a partition is greater than the other.

        Examples:
        >>> from sympy.combinatorics.partitions import Partition
        >>> a = Partition([[1,2],[3,4,5]])
        >>> b = Partition([[1],[2,3],[4],[5]])
        >>> a > b
        True
        """
        return self._compare(other) > 0

    def __lt__(self, other):
        """
        Checks if a partition is less than the other.

        Examples:
        >>> from sympy.combinatorics.partitions import Partition
        >>> a = Partition([[1,2],[3,4,5]])
        >>> b = Partition([[1],[2,3],[4],[5]])
        >>> a < b
        False
        """
        return self._compare(other) < 0

    def __ge__(self, other):
        """
        Checks if a partition is greater than or equal to
        the other partition.

        Examples:
        >>> from sympy.combinatorics.partitions import Partition
        >>> a = Partition([[1,2],[3,4,5]])
        >>> b = Partition([[1],[2,3],[4],[5]])
        >>> a >= a
        True
        >>> a >= b
        True
        """
        return self == other or self > other

    def __le__(self, other):
        """
        Checks if a partition is less than or equal to
        the other partition.

        Examples:
        >>> from sympy.combinatorics.partitions import Partition
        >>> a = Partition([[1,2],[3,4,5]])
        >>> b = Partition([[1],[2,3],[4],[5]])
        >>> a <= a
        True
        >>> a <= b
        False
        """
        return self == other or self < other

    @property
    def rank(self):
        """
        Gets the rank of a partition.

        Examples:
        >>> from sympy.combinatorics.partitions import Partition
        >>> a = Partition([[1,2],[3],[4,5]])
        >>> a.rank
        13
        """
        if self._rank != None:
            return self._rank
        self._rank = RGS_rank(self.RGS)
        return self._rank

    @property
    def RGS(self):
        """
        Returns the restricted growth string of the partition.

        For a set partition of consisting of n elements, the n-character string
        a_1, a_2,..., a_n, in which each character gives the set block (B_0, B_1, ...)
        in which the corresponding element belongs is called the restricted growth
        string or the restricted growth function. 

        Examples:
        >>> from sympy.combinatorics.partitions import Partition
        >>> a = Partition([[1,2],[3],[4,5]])
        >>> a.RGS
        [0, 0, 1, 2, 2]
        >>> a = Partition([[1,4],[2],[3,5]])
        >>> a.RGS
        [0, 1, 2, 0, 2]
        """
        rgs = [0] * self.partition_set_size
        a = 0
        for part in self.partition_list_form:
            for i in part:
                rgs[self.partition_set_list.index(i)] = a
            a += 1
        return rgs

    @classmethod
    def partition_from_rgs(self, rgs, superset):
        """
        Creates a set partition from a restricted growth string.

        Examples:
        >>> from sympy.combinatorics.partitions import Partition
        >>> Partition.partition_from_rgs([0,1,2,0,1],['a','b','c','d','e'])
        {{a, d}, {b, e}, {c}}
        >>> a = Partition([[1,4],[2],[3,5]])
        >>> Partition.partition_from_rgs(a.RGS, a.partition_set_list)
        {{1, 4}, {2}, {3, 5}}
        """
        superset = list(superset)
        max_elem = max(rgs) + 1
        partition = [[] for i in xrange(max_elem)]
        j = 0
        for i in rgs:
            partition[i].append(superset[j])
            j += 1
        return Partition(partition)


class IntegerPartition(Basic):
    """
    This class represents an integer partition.

    In number theory and combinatorics, a partition of a positive integer n,
    also called an integer partition, is a way of writing n as a sum of positive
    integers. Two sums that differ only in the order of their summands are
    considered to be the same partition; if order matters then the sum becomes a
    composition. For example, 4 can be partitioned in five distinct ways:
    [[4],[3,1],[2,2],[2,1,1],[1,1,1,1]].
    
    The order-dependent composition [1, 3] is the same partition as [3, 1],
    while [1, 2, 1] and [1, 1, 2] are the same partition as [2, 1, 1].

    Reference: http://en.wikipedia.org/wiki/Partition_(number_theory)
    """

    def __new__(cls, *args, **kw_args):
        """
        Generates a new IntegerPartition object.
        
        It also verifies if the arguments passed are valid and if it is
        found that they are not then an exception is raised. The arguments
        taken are the integer partition and the integer itself. We simply
        check if the argument types are what we expect and if the partition
        is valid.

        Examples:
        >>> from sympy.combinatorics.partitions import IntegerPartition
        >>> a = IntegerPartition([5,4,3,1,1,1], 15)
        >>> a
        IntegerPartition([5, 4, 3, 1, 1, 1], 15)
        >>> b = IntegerPartition([5,4,3,1], 13)
        >>> print b
        [5, 4, 3, 1]
        """
        partition = args[0]
        integer_rep = args[1]

        if not isinstance(partition, list) or sum(partition) != integer_rep:
            raise ValueError("The partition is not valid")

        list.sort(args[0], key=lambda x: -x)

        obj = Basic.__new__(cls, *args, **kw_args)
        obj.partition = partition
        return obj

    def next(self):
        """
        Generates the next partition.

        Examples:
        """
        raise NotImplementedError("The method to generate the next integer partition \
        is not implemented yet")

    def previous(self):
        """
        Generates the previous partition.

        Examples:
        """
        raise NotImplementedError("The method to generate the previous integer \
        partition is not implemented yet")

    @property
    def size(self):
        """
        Gets the size of the partition.

        Examples:
        >>> from sympy.combinatorics.partitions import IntegerPartition
        >>> a = IntegerPartition([1,3,4], 8)
        >>> a.size
        3
        """
        return len(self.args[0])

    @property
    def partition_set(self):
        """
        Gets the set of the partition.

        Examples:
        >>> from sympy.combinatorics.partitions import IntegerPartition
        >>> a = IntegerPartition([1,3,4], 8)
        >>> a.partition_set
        8
        """
        return self.args[1]

    @property
    def partition_array(self):
        """
        Gets the array of partitions from the
        partition object

        Examples:
        >>> from sympy.combinatorics.partitions import IntegerPartition
        >>> a = IntegerPartition([1,3,4], 8)
        >>> a.partition_array
        [4, 3, 1]
        """
        return self.args[0]

    @property
    def is_self_conjugate(self):
        """
        Checks if the conjugate of a partition equals itself.

        Examples:
        >>> from sympy.combinatorics.partitions import IntegerPartition
        >>> a = IntegerPartition([6,3,3,2,1], 15)
        >>> a.is_self_conjugate
        False
        >>> a = IntegerPartition([3,2,1], 6)
        >>> a.is_self_conjugate
        True
        """
        return self.conjugate_partition == self.partition_array

    @property
    def conjugate_partition(self):
        """
        Computes the conjugate partition of itself.

        Examples:
        >>> from sympy.combinatorics.partitions import \
        IntegerPartition
        >>> a = IntegerPartition([6,3,3,2,1], 15)
        >>> a.conjugate_partition
        [5, 4, 3, 1, 1, 1]
        >>> b = IntegerPartition(a.conjugate_partition, 15)
        >>> print b.ferrers_representation()
        #####
        ####
        ###
        #
        #
        #
        """
        j = 1
        temp_arr = self.partition_array[:] + [0]
        k = temp_arr[0]
        b = [0] * (k)
        while k > 0:
            while k > temp_arr[j]:
                b[k - 1] = j
                k -= 1
            j += 1
        return b

    def _compare(self, other):
        """
        Compares two partitions.

        The basis for comparison of two integer partitions is the
        majorizing concept. A majorization is a partial order on
        vectors of real numbers.

        References:

        [1] Inequalities: Theory of Majorization and Its Applications (1980)
        Albert W. Marshall, Ingram Olkin, Academic Press 

        Examples:
        >>> from sympy.combinatorics.partitions import IntegerPartition
        >>> a = IntegerPartition([1,1,1,1,1], 5)
        >>> b = IntegerPartition([1,1,1,2], 5)
        >>> a.compare(b)
        -1
        >>> a < b
        True
        """
        k = min(self.size, other.size)
        val_self = 0
        val_other = 0
        for i in xrange(k):
            val_self += self.partition_array[i]
            val_other += other.partition_array[i]
            if val_self > val_other:
                return 1
            elif val_self < val_other:
                return -1
        return 1

    @property
    def rank(self):
        """
        Gets the rank of a partition.

        Examples:
        """
        raise NotImplementedError()

    def ferrers_representation(self):
        """
        Prints the ferrer diagram of a partition.

        Examples:
        >>> from sympy.combinatorics.partitions import IntegerPartition
        >>> a = IntegerPartition([3,2,1], 6)
        >>> b = IntegerPartition([5,1,1], 7)
        >>> print a.ferrers_representation()
        ###
        ##
        #
        >>> print b.ferrers_representation()
        #####
        #
        #
        """
        return "\n".join(['#'*i for i in self.partition_array])

    def __str__(self):
        return str(self.partition)

def random_integer_partition(n):
    """
    Generates a random integer partition.
    """
    partition = []
    while(n > 0):
        k = random.randint(1, n)
        partition.append(k)
        n -= k
    list.sort(partition)
    return partition

def RGS_generalized(m):
    """
    Computes the generalized unrestricted growth strings.

    Examples:
    >>> from sympy.combinatorics.partitions import \
    RGS_generalized
    >>> RGS_generalized(6)
    [  1,   1,   1,  1,  1, 1, 1]
    [  1,   2,   3,  4,  5, 6, 0]
    [  2,   5,  10, 17, 26, 0, 0]
    [  5,  15,  37, 77,  0, 0, 0]
    [ 15,  52, 151,  0,  0, 0, 0]
    [ 52, 203,   0,  0,  0, 0, 0]
    [203,   0,   0,  0,  0, 0, 0]
    """
    d = zeros(m + 1)
    for i in xrange(0, m + 1):
        d[0, i] = 1

    for i in xrange(1, m + 1):
        for j in xrange(m):
            if j <= m - i:
                d[i, j] = j * d[i - 1, j] \
                          + d[i - 1, j + 1]
            else:
                d[i, j] = 0
    return d

def RGS_enum(m):
    """
    RGS_enum computes the total number of restricted growth strings
    possible for a superset of size m.

    Examples:
    >>> from sympy.combinatorics.partitions import RGS_enum
    >>> RGS_enum(4)
    15
    >>> RGS_enum(5)
    52
    >>> RGS_enum(6)
    203
    """
    m += 1
    if (m < 0):
        return 0
    elif (m == 0):
        return 1
    else:
        b = [1] * (m)
        for j in xrange(1, m):
            for i in xrange(1, j):
                b[j] += C.binomial(j - 1, i) * b[i]

        nrgf = b[m - 1]
    return nrgf

def RGS_unrank(rank, m):
    """
    Gives the unranked restricted growth string for a given
    superset size.

    Examples:
    >>> from sympy.combinatorics.partitions import RGS_unrank
    >>> RGS_unrank(14, 4)
    [0, 1, 2, 3]
    >>> RGS_unrank(0, 4)
    [0, 0, 0, 0]
    """
    if m < 1:
        raise ValueError("The superset size must be >= 1")
    if rank < 0 or RGS_enum(m) <= rank:
        raise ValueError("Invalid arguments")

    L = [1] * (m + 1)
    j = 1
    D = RGS_generalized(m)
    for i in xrange(2,m+1):
        v = D[m - i,j]
        cr = j*v
        if cr <= rank:
            L[i] = j + 1
            rank -= cr
            j += 1
        else:
            L[i] = int(rank / v + 1)
            rank %= v
    return map(lambda x: x - 1, L[1:])

def RGS_rank(rgs):
    """
    Computes the rank of a restricted growth string.

    Examples:
    >>> from sympy.combinatorics.partitions import RGS_rank, RGS_unrank
    >>> RGS_rank([0, 1, 2, 1, 3])
    42
    >>> RGS_rank(RGS_unrank(4,7))
    4
    """
    rgs_size = len(rgs)
    rank = 0
    D = RGS_generalized(rgs_size)
    for i in xrange(1, rgs_size):
        n = len(rgs[(i + 1):])
        m = max(rgs[0:i])
        rank += D[n, m + 1] * rgs[i]
    return rank
