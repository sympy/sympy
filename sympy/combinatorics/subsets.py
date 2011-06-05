from sympy.core import Basic
from sympy.combinatorics.graycode import GrayCode

import itertools

class Subset(Basic):
    """
    Represents a basic subset object.

    We generate subsets using essentially two techniques,
    binary enumeration and lexicographic enumeration.
    """

    _rank_binary = None
    _rank_lex = None
    _rank_graycode = None

    def next_binary(self):
        """
        Generates the next binary ordered subset.

        Examples:
        >>> from sympy.combinatorics.subsets import Subset
        >>> a = Subset(['c','d'],['a','b','c','d'])
        >>> a.next_binary()
        ['b']
        >>> a = Subset(['a','b','c','d'],['a','b','c','d'])
        >>> a.next_binary()
        []
        """
        bin_list = get_bitlist_from_subset(self.subset, self.superset)
        next_bin_list = list(bin((int(reduce(lambda x, y:
                                             x + y, bin_list), 2) + 1)
                                 % 2**self.superset_size))[2:]
        next_bin_list = [0] * (self.superset_size - len(next_bin_list)) + \
                        next_bin_list
        return get_subset_from_bitlist(self.superset, next_bin_list)

    def prev_binary(self):
        """
        Generates the previous binary ordered subset.

        Examples:
        >>> from sympy.combinatorics.subsets import Subset
        >>> a = Subset([],['a','b','c','d'])
        >>> a.prev_binary()
        ['a', 'b', 'c', 'd']
        >>> a = Subset(['c','d'],['a','b','c','d'])
        >>> a.prev_binary()
        ['c']
        """
        bin_list = get_bitlist_from_subset(self.subset, self.superset)
        next_bin_list = list(bin((int(reduce(lambda x, y:
                                             x + y, bin_list), 2) - 1)
                                 % 2**self.superset_size))[2:]
        next_bin_list = [0] * (self.superset_size - len(next_bin_list)) + \
                        next_bin_list
        return get_subset_from_bitlist(self.superset, next_bin_list)

    def next_lexicographic(self):
        """
        Generates the next lexicographically ordered
        subset.

        Examples:
        """
        raise NotImplementedError()

    def prev_lexicographic(self):
        """
        Generates the previous lexicographically ordered
        subset.

        Examples:
        """
        raise NotImplementedError()

    def next_graycode(self):
        """
        Generates the next gray code ordered subset.

        Examples:
        >>> from sympy.combinatorics.subsets import *
        >>> a = Subset([1,2],[1,2,3,4])
        >>> a.next_graycode()
        [1, 2, 4]
        """
        bin_list = get_bitlist_from_subset(self.subset, self.superset)
        gc_iter = GrayCode(self.superset_size).\
                  generate_bitlist(bin_list)
        gc_iter = gc_iter.next()
        return get_subset_from_bitlist(self.superset, gc_iter.next())

    def prev_graycode(self):
        """
        Generates the previous gray code ordered subset.

        Examples:
        """
        raise NotImplementedError()

    @property
    def rank_binary(self):
        """
        Computes the binary ordered rank.

        Examples:
        >>> from sympy.combinatorics.subsets import Subset
        >>> a = Subset([],['a','b','c','d'])
        >>> a.rank_binary
        0
        >>> a = Subset(['c','d'],['a','b','c','d'])
        >>> a.rank_binary
        3
        """
        if self._rank_binary == None:
            self._rank_binary = int("".join(
                get_bitlist_from_subset(self.subset,
                                        self.superset)), 2)
        return self._rank_binary

    @property
    def rank_lexicographic(self):
        """
        Computes the lexicographic ranking of the subset.

        Examples:
        >>> from sympy.combinatorics.subsets import Subset
        >>> a = Subset(['c','d'],['a','b','c','d'])
        >>> a.rank_lexicographic
        14
        >>> a = Subset([2,4,5],[1,2,3,4,5,6])
        >>> a.rank_lexicographic
        43
        """
        if self._rank_lex == None:
            superset_index = list(xrange(self.superset_size))
            subset_index = [self.superset.index(j) for j in self.subset]
            def ranklex(self, subset_index, i, n):
                if subset_index == [] or i > n:
                    return 0
                if i in subset_index:
                    subset_index.remove(i)
                    return 1 + ranklex(self, subset_index, i + 1, n)
                return 2**(n - i - 1) + ranklex(self, subset_index, i + 1, n)
            self._rank_lex = ranklex(self, subset_index, 0, self.superset_size)
        return self._rank_lex

    @property
    def rank_graycode(self):
        """
        Computes the gray code ranking of the subset.

        Examples:
        >>> from sympy.combinatorics.subsets import Subset
        >>> a = Subset(['c','d'],['a','b','c','d'])
        >>> a.rank_graycode
        8
        >>> a = Subset([2,4,5],[1,2,3,4,5,6])
        >>> a.rank_graycode
        19
        """
        if self._rank_graycode == None:
            bits = get_bitlist_from_subset(self.subset, self.superset)
            self._rank_graycode = GrayCode(start = bits).rank
        return self._rank_graycode

    @property
    def subset(self):
        """
        Gets the subset represented by the current instance.
        """
        return self.args[0]

    @property
    def size(self):
        """
        Gets the size of the subset.
        """
        return len(self.subset)

    @property
    def superset(self):
        """
        Gets the superset of the subset.
        """
        return self.args[1]

    @property
    def superset_size(self):
        """
        Returns the size of the superset.
        """
        return len(self.superset)

def get_subset_from_bitlist(super_set, bitlist):
    """
    Gets the subset defined by the bitlist.

    Examples:
    >>> from sympy.combinatorics.subsets import get_subset_from_bitlist
    >>> get_subset_from_bitlist(['a','b','c','d'],['0','0','1','1'])
    ['c', 'd']
    """
    if len(super_set) != len(bitlist):
        raise ValueError("The sizes of the lists are not equal")
    ret_set = []
    for i in xrange(len(bitlist)):
        if bitlist[i] == '1':
            ret_set.append(super_set[i])
    return ret_set

def get_bitlist_from_subset(subset, superset):
    """
    Gets the bitlist corresponding to a subset.

    Examples:
    >>> from sympy.combinatorics.subsets import get_bitlist_from_subset
    >>> get_bitlist_from_subset(['c','d'],['a','b','c','d'])
    ['0', '0', '1', '1']
    """
    bitlist = ['0'] * len(superset)
    for i in subset:
        bitlist[superset.index(i)] = '1'
    return bitlist

def unrank_binary(rank, superset):
    """
    Gets the binary ordered subset of the
    specified rank.

    Examples:
    >>> from sympy.combinatorics.subsets import unrank_binary
    >>> unrank_binary(4,['a','b','c','d'])
    ['b']
    """
    bin_list = list(bin(rank))[2:]
    bin_list = [0] * (len(superset) - len(bin_list)) + bin_list
    return get_subset_from_bitlist(superset, bin_list)

def ksubsets(superset, k):
    """
    Finds the subsets of size k in lexicographic order.

    This uses the itertools generator.

    Examples:
    >>> from sympy.combinatorics.subsets import ksubsets
    >>> list(ksubsets([1,2,3], 2))
    [(1, 2), (1, 3), (2, 3)]
    >>> list(ksubsets([1,2,3,4,5], 2))
    [(1, 2), (1, 3), (4, 5), (1, 4), (1, 5), (2, 3), \
    (2, 5), (3, 4), (2, 4), (3, 5)]
    """
    return set(itertools.combinations(superset, k))
