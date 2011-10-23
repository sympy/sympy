from sympy.core import Basic

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

    The Gray code solves the problem of sequentially
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
    """

    _reset = False
    _current = 0

    def __new__(cls, *args, **kw_args):
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
        Returns the number of objects.

        Examples:
        >>> from sympy.combinatorics.graycode import GrayCode
        >>> a = GrayCode(5)
        >>> a.n
        5
        """
        return self.args[0]

    def generate_bitlist(self, start = None):
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
            graycode_bin = gray_to_bin(self._current)
            self._current = int(''.join(self._current), 2)
        else:
            graycode_bin = gray_to_bin(list(bin(self._current))[2:])
        graycode_bin = int(''.join(graycode_bin), 2)
        for i in xrange(graycode_bin, 1 << bits):
            if self._reset:
                self._current = 0
                self._reset = False
                break
            retlist = list(bin(self._current))[2:]
            yield ['0'] * (self.n - len(retlist)) + retlist
            bbtc = (i ^ (i + 1))
            gbtc = (bbtc ^ (bbtc >> 1))
            self._current = (self._current ^ gbtc)

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

def random_bitlist(n):
    """
    Generates a random bitlist of length n.
    """
    return [random.randint(0, 1) for i in xrange(n)]

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
