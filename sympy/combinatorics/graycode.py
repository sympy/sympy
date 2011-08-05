from sympy.core import Basic
from sympy.functions import ceiling, log

import random

class GrayCode(Basic):
    """
    A Gray code is essentially a Hamiltonian walk on
    a n-dimensional cube with edge length of one.
    The vertices of the cube are represented by vectors
    whose values are binary. The Hamilton walk visits
    each vertex exactly once. The Gray code for a 3d
    cube is [[0,0,0],[1,0,0],[1,1,0],[0,1,0],[0,1,1],
    [1,1,1],[1,0,1],[0,0,1]].

    A Gray code solves the problem of sequentially
    generating all possible subsets of n objects in such
    a way that each subset is obtained from the previous
    one by either deleting or adding a single object.
    In the above example, 1 indicates that the object is
    present, and 0 indicates that its absent.

    Gray codes have applications in statistics as well when
    we want to compute various statistics related to subsets
    in an efficient manner.

    References:
    [1] Nijenhuis,A. and Wilf,H.S.(1978).
    Combinatorial Algorithms. Academic Press.
    [2] Knuth, D. (2011). The Art of Computer Programming, Vol 4
    Addison Wesley

    Examples:
    >>> from sympy.combinatorics.graycode import GrayCode
    >>> a = GrayCode(3)
    >>> list(a.generate_bitlist())
    [['0', '0', '0'], ['0', '0', '1'], ['0', '1', '1'], \
    ['0', '1', '0'], ['1', '1', '0'], ['1', '1', '1'], \
    ['1', '0', '1'], ['1', '0', '0']]
    >>> a = GrayCode(4)
    >>> list(a.generate_bitstring())
    ['0000', '0001', '0011', '0010', '0110', '0111', '0101', '0100', \
    '1100', '1101', '1111', '1110', '1010', '1011', '1001', '1000']
    """

    _reset = False
    _current = 0
    _rank = None

    def __new__(cls, *args, **kw_args):
        """
        Default constructor
        >>> from sympy.combinatorics.graycode import GrayCode
        >>> a = GrayCode(3)
        >>> a
        GrayCode(3)
        >>> a = GrayCode(3, start=['1','0','0'])
        >>> a.selections
        8
        """
        obj = Basic.__new__(cls, *args, **kw_args)
        if kw_args.has_key("start"):
            obj._current = kw_args["start"]
        return obj

    @property
    def selections(self):
        """
        Returns the number of bit vectors in the
        Gray code.

        Examples:
        >>> from sympy.combinatorics.graycode import GrayCode
        >>> a = GrayCode(3)
        >>> a.selections
        8
        """
        return 2**self.n

    @property
    def n(self):
        """
        Returns the dimension of the gray code.

        Examples:
        >>> from sympy.combinatorics.graycode import GrayCode
        >>> a = GrayCode(5)
        >>> a.n
        5
        """
        return self.args[0]

    def generate_bitlist(self, start=None):
        """
        Generates the sequence of bit vectors of a Gray Code.

        [1] Knuth, D. (2011). The Art of Computer Programming,
        Vol 4, Addison Wesley

        Examples:
        >>> from sympy.combinatorics.graycode import GrayCode
        >>> a = GrayCode(3)
        >>> list(a.generate_bitlist())
        [['0', '0', '0'], ['0', '0', '1'], ['0', '1', '1'], \
        ['0', '1', '0'], ['1', '1', '0'], ['1', '1', '1'], \
        ['1', '0', '1'], ['1', '0', '0']]
        >>> list(a.generate_bitlist(['0', '1', '1']))
        [['0', '1', '1'], ['0', '1', '0'], ['1', '1', '0'], \
        ['1', '1', '1'], ['1', '0', '1'], ['1', '0', '0']]
        """
        bits = self.args[0]
        if start != None:
            self._current = start
        if isinstance(self._current, list):
            graycode_int = gray_to_bin(self._current)
            self._current = int(''.join(self._current), 2)
        else:
            graycode_int = gray_to_bin(list(bin(self._current))[2:])
        graycode_int = int(''.join(graycode_int), 2)
        for i in xrange(graycode_int, 1 << bits):
            if self._reset:
                self._current = 0
                self._reset = False
                break
            retlist = list(bin(self._current))[2:]
            yield ['0'] * (self.n - len(retlist)) + retlist
            bbtc = (i ^ (i + 1))
            gbtc = (bbtc ^ (bbtc >> 1))
            self._current = (self._current ^ gbtc)

    def generate_bitstring(self, start=None):
        """
        Generates bitstrings instead of bitlists.

        Examples:
        >>> from sympy.combinatorics.graycode import GrayCode
        >>> a = GrayCode(3)
        >>> list(a.generate_bitstring())
        ['000', '001', '011', '010', '110', '111', '101', '100']
        """
        return [''.join(r) for r in self.generate_bitlist(start)]

    def reset(self):
        """
        Resets the bit generation.

        Examples:
        >>> from sympy.combinatorics.graycode import GrayCode
        >>> a = GrayCode(3)
        >>> for i in a.generate_bitlist():
        ...     if i == ['0', '1', '0']:
        ...         a.reset()
        ...     print i
        ...
        ['0', '0', '0']
        ['0', '0', '1']
        ['0', '1', '1']
        ['0', '1', '0']
        """
        self._reset = True

    @property
    def rank(self):
        """
        Ranks the gray code.

        The rank of a gray code represents the number of gray
        codes that precede it in its family. Thus a gray code
        may be common to both the binary-reflected group as
        well as the Frank-Ruskey group buts its rank in these
        two families will differ.

        Examples:
        >>> from sympy.combinatorics.graycode import GrayCode
        >>> a = GrayCode(3, start=['1','0','0'])
        >>> a.rank
        1
        >>> a = GrayCode(5, start=['1','0','1','0','0'])
        >>> a.rank
        6
        """
        if self._rank is None:
            if isinstance(self._current, int):
                self._current = list(bin(self._current))[2:]
            if len(self._current)==0:
                return 0
            elif self._current[-1]=='0':
                self._rank = GrayCode(start = self._current[:-1]).rank
            else:
                self._rank =  2**len(self._current) - \
                             GrayCode(start = self._current[:-1]).rank - 1
        return self._rank

    @property
    def current(self):
        """
        Returns the currently referenced gray code.
        """
        return self._current

    @classmethod
    def unrank_gray_code(self, k, n):
        """
        Unranks an n-bit sized gray code of rank k.

        We generate in reverse order to allow for tail-call
        optimization.

        Examples:
        >>> from sympy.combinatorics.graycode import GrayCode
        >>> GrayCode.unrank_gray_code(3, 5).current
        ['0', '1', '0', '0', '0']
        >>> GrayCode.unrank_gray_code(7, 10).rank
        7
        """
        def unrank(k, n):
            if n == 1:
                return [str(k % 2)]
            m = 2**(n - 1)
            if k < m:
                return ["0"] + unrank(k, n - 1)
            return ["1"] + unrank(m - (k % m) - 1, n - 1)
        ret_list = unrank(k, n)
        list.reverse(ret_list)
        return GrayCode(start = ret_list)


def random_bitlist(n):
    """
    Generates a random bitlist of length n.
    """
    return [random.choice('01') for i in xrange(n)]

def gray_to_bin(bin_list):
    """
    Convert from Gray coding to binary coding.

    We assume big endian encoding.

    Examples:
    >>> from sympy.combinatorics.graycode import gray_to_bin
    >>> gray_to_bin(['1','0','0'])
    ['1', '1', '1']
    """
    b = bin_list[0]
    for i in xrange(1, len(bin_list)):
        b += str(int(b[i-1] != bin_list[i]))

    return list(b)

def bin_to_gray(bin_list):
    """
    Convert from binary coding to gray coding.

    We assume big endian encoding.

    Examples:
    >>> from sympy.combinatorics.graycode import bin_to_gray
    >>> bin_to_gray(['1','1','1'])
    ['1', '0', '0']
    """
    b = bin_list[0]
    for i in xrange(0, len(bin_list) - 1):
        b += str(int(bin_list[i]) ^ int(b[i - 1]))

    return list(b)

def get_subset_from_bitlist(super_set, bitlist):
    """
    Gets the subset defined by the bitlist.

    Examples:
    >>> from sympy.combinatorics.graycode import get_subset_from_bitlist
    >>> get_subset_from_bitlist(['a','b','c','d'],['0','0','1','1'])
    ['c', 'd']
    >>> get_subset_from_bitlist(['c','a','c','c'],['1','1','0','0'])
    ['c', 'a']
    """
    if len(super_set) != len(bitlist):
        raise ValueError("The sizes of the lists are not equal")
    return [super_set[i] for i, j in enumerate(bitlist) if bitlist[i] == '1']

def gray_code_subsets(gray_code_set):
    """
    Generates the subsets as enumerated by a Gray code.

     Examples:
    >>> from sympy.combinatorics.graycode import gray_code_subsets
    >>> list(gray_code_subsets(['a','b','c']))
    [[], ['c'], ['b', 'c'], ['b'], ['a', 'b'], ['a', 'b', 'c'], \
    ['a', 'c'], ['a']]
    >>> list(gray_code_subsets(['a','b','c','c']))
    [[], ['c'], ['c', 'c'], ['c'], ['b', 'c'], ['b', 'c', 'c'], \
    ['b', 'c'], ['b'], ['a', 'b'], ['a', 'b', 'c'], ['a', 'b', 'c', 'c'], \
    ['a', 'b', 'c'], ['a', 'c'], ['a', 'c', 'c'], ['a', 'c'], ['a']]
    """
    return [get_subset_from_bitlist(gray_code_set, bitlist) for \
            bitlist in list(GrayCode(len(gray_code_set)).generate_bitlist())]
