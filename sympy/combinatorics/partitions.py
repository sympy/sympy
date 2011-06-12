from sympy.core import Basic

import random

class Partition(Basic):
    """
    This class represents an abstract partition.
    A partition is a set of disjoint sets whose
    union equals a given set.
    """

    def next(self):
        """
        Generates the next partition.

        Examples:
        """
        raise NotImplementedError()

    def previous(self):
        """
        Generates the previous partition.

        Examples:
        """
        raise NotImplementedError()

    @property
    def size(self):
        """
        Gets the size of the partition.

        Examples:
        >>> from sympy.combinatorics.partitions import Partition
        >>> a = Partition([[1,2],[3]],[1,2,3])
        >>> a.size
        2
        """
        return len(self.args[0])

    @property
    def partition(self):
        """
        Gets the partition itself.

        Examples:
        >>> from sympy.combinatorics.partitions import Partition
        >>> a = Partition([[1,2],[3]],[1,2,3])
        >>> a.partition
        [[1, 2], [3]]
        """
        return self.args[0]

    @property
    def partition_set(self):
        """
        Gets the set of the partition.

        Examples:
        >>> from sympy.combinatorics.partitions import Partition
        >>> a = Partition([[1,2],[3]],[1,2,3])
        >>> a.partition_set
        [1, 2, 3]
        """
        return self.args[1]

    @property
    def partition_set_size(self):
        """
        Gets the set of the partition.

        Examples:
        >>> from sympy.combinatorics.partitions import Partition
        >>> a = Partition([[1,2],[3]],[1,2,3])
        >>> a.partition_set_size
        3
        """
        return len(self.args[1])

    def __str__(self):
        return str(self.partition)

    def __repr__(self):
        return str(self)

    def __new__(cls, *args, **kw_args):
        """
        Generates a new partition object.
        It also verifies if the arguments passed are
        valid and if it is found that they are not then
        an exception is raised.

        Examples:
        >>> from sympy.combinatorics.partitions import Partition
        >>> a = Partition([[1,2],[3]],[1,2,3])
        >>> str(a)
        '[[1, 2], [3]]'
        """
        partition = args[0]
        super_set = args[1][:]

        check = []
        for part in partition:
            if not isinstance(part, list):
                raise ValueError("The input has been provided incorrectly")
            check.extend(part)

        check.sort()
        super_set.sort()

        if check != super_set:
            raise ValueError("The partition provided is not valid.")
        obj = Basic.__new__(cls, *args, **kw_args)
        return obj

    def _compare(self, other):
        """
        Compares two partitions.

        The basis for comparison of two partitions is rank.
        A partition with a lesser rank is greater than a
        partition with a greater rank.

        Examples:
        """
        raise NotImplementedError()

    def __eq__(self, other):
        """
        Checks for equality of two partitions.

        Examples:
        """
        return self._compare(other) == 0

    def __ne__(self, other):
        """
        Checks for inequality of two partitions.

        Examples:
        """
        return self._compare(other) != 0

    def __gt__(self, other):
        """
        Checks if a partition is greater than the other.

        Examples:
        """
        return self._compare(other) > 0

    def __lt__(self, other):
        """
        Checks if a partition is less than the other.

        Examples:
        """
        return self._compare(other) < 0

    def __ge__(self, other):
        """
        Checks if a partition is greater than or equal to
        the other partition.

        Examples:
        """
        return self == other or self > other

    def __le__(self, other):
        """
        Checks if a partition is less than or equal to
        the other partition.

        Examples:
        """
        return self == other or self < other

    @property
    def rank(self):
        """
        Gets the rank of a partition.

        Examples:
        """
        raise NotImplementedError()

    @property
    def RGS(self):
        """
        Returns the restricted growth string of the partition.

        Examples:
        >>> from sympy.combinatorics.partitions import Partition
        >>> a = Partition([[1,2],[3],[4,5]], [1,2,3,4,5])
        >>> a.RGS
        [0, 0, 1, 2, 2]
        >>> a = Partition([[1,4],[2],[3,5]], [1,2,3,4,5])
        >>> a.RGS
        [0, 1, 2, 0, 2]
        """
        rgs = [0] * self.partition_set_size
        a = 0
        for part in self.partition:
            for i in part:
                rgs[self.partition_set.index(i)] = a
            a += 1
        return rgs

def from_RGS(rgs, superset):
    """
    Creates a set partition from a restricted growth string.

    Examples:
    >>> from sympy.combinatorics.partitions import *
    >>> from_RGS([0,1,2,0,1],['a','b','c','d','e'])
    Partition([['a', 'd'], ['b', 'e'], ['c']], \
    ['a', 'b', 'c', 'd', 'e'])
    >>> a = Partition([[1,4],[2],[3,5]], [1,2,3,4,5])
    >>> from_RGS(a.RGS, a.partition_set)
    Partition([[1, 4], [2], [3, 5]], [1, 2, 3, 4, 5])
    """
    max_elem = max(rgs) + 1
    partition = [[] for i in xrange(max_elem)]
    j = 0
    for i in rgs:
        partition[i].append(superset[j])
        j += 1
    return Partition(partition, superset)

class IntegerPartition(Partition):
    """
    This class represents an abstract partition.
    A partition is a set of disjoint sets whose
    union equals a given set.
    """

    def next(self):
        """
        Generates the next partition.

        Examples:
        """
        raise NotImplementedError()

    def previous(self):
        """
        Generates the previous partition.

        Examples:
        """
        raise NotImplementedError()

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
    def conjugate(self):
        """
        Find the conjugate of a partition.
        This is the vector that satisfies
        len(p) = max(conjugate(p)) and vice versa.

        Examples:
        >>> from sympy.combinatorics.partitions import IntegerPartition
        >>> a = IntegerPartition([1,3,4], 8)
        >>> a.conjugate
        [3, 2, 2, 1]
        """
        result = []
        j = len(self.partition_array)
        if j <= 0:
            return result
        while True:
            result.append(j)
            while len(result) >= self.partition_array[j-1]:
                j -= 1
                if j == 0:
                    return result

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
        >>> a = IntegerPartition([5,4,3,1,1,1], 15)
        >>> a.conjugate_partition
        [6, 3, 3, 2, 1]
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

    def __new__(cls, *args, **kw_args):
        """
        Generates a new partition object.
        It also verifies if the arguments passed are
        valid and if it is found that they are not then
        an exception is raised.

        Examples:
        """
        partition = args[0]
        integer_rep = args[1]

        if not isinstance(partition, list) or sum(partition) != integer_rep:
            raise ValueError("The partition is not valid")

        list.sort(args[0], key = lambda x: -x)

        obj = Basic.__new__(cls, *args, **kw_args)
        return obj

    def _compare(self, other):
        """
        Compares two partitions.

        The basis for comparison of two integer partitions is the
        majorizing concept.

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
