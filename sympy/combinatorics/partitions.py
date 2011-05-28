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
        [1, 3, 4]
        """
        return self.args[0]

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
