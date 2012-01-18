from sympy.core import Basic, S
from sympy.core.compatibility import is_sequence
from sympy.utilities.iterables import flatten
from sympy.polys.polytools import lcm
from sympy.matrices import zeros
from sympy.mpmath.libmp.libintmath import ifac

import random

def perm_af_parity(pi):
    """
    Computes the parity of a permutation in array form.

    The parity of a permutation reflects the parity of the
    number of inversions in the permutation, i.e., the
    number of pairs of x and y such that x > y but p[x] < p[y].

    Examples
    ========

    >>> from sympy.combinatorics.permutations import perm_af_parity
    >>> perm_af_parity([0,1,2,3])
    0
    >>> perm_af_parity([3,2,0,1])
    1

    See Also
    ========
    Permutation
    """
    n = len(pi)
    a = [0] * n
    c = 0
    for j in xrange(n):
        if a[j] == 0:
            c += 1
            a[j] = 1
            i = j
            while pi[i] != j:
                i = pi[i]
                a[i] = 1
    return (n - c) % 2

def perm_af_mul(a, b):
    """
    Product of two permutations in array form

    Examples
    ========

    >>> from sympy.combinatorics.permutations import perm_af_mul
    >>> perm_af_mul([1,2,3,0], [3,2,0,1])
    [0, 3, 1, 2]

    See Also
    ========
    Permutation
    """
    if len(a) != len(b):
        raise ValueError("The number of elements in the permutations \
don\'t match.")

    return [a[i] for i in b]

class Permutation(Basic):
    """
    A permutation, alternatively known as an 'arrangement number'
    or 'ordering' is an arrangement of the elements of an ordered list
    into a one-to-one mapping with itself. The number of permutations
    on a set of n elements is given by n!.

    A representation of a permutation as a product of permutation cycles
    is unique (up to the ordering of the cycles). An example of a cyclic
    decomposition is the permutation [3, 1, 0, 2] of the set [0, 1, 2, 3].
    This is denoted as [[1], [0, 3, 2]], corresponding to the disjoint
    permutation cycles [1] and [0, 3, 2]. We can choose the cyclic form as
    we want since the cycles are disjoint and can therefore be specified
    in any order and a rotation of a given cycle specifies the same cycle [see 1.]
    Therefore, (320)(1), (203)(1), (032)(1), (1)(320), (1)(203), and (1)(032)
    all describe the same permutation.

    Another notation that explicitly identifies the positions occupied by
    elements before and after application of a permutation on n elements uses a
    2xn matrix, where the first row is the identity permutation and the second
    row is the new arrangement [see 2.].

    Any permutation is also a product of transpositions.

    Permutations are commonly denoted in lexicographic or transposition order.

    The product of two permutations a and q is defined as their composition as
    functions, (p*q)(i) = p(q(i)) [see 6.].

    (1) Skiena, S. 'Permutations.' 1.1 in Implementing Discrete Mathematics
        Combinatorics and Graph Theory with Mathematica.
        Reading, MA: Addison-Wesley, pp. 3-16, 1990.
    (2) Knuth, D. E. The Art of Computer Programming, Vol. 4: Combinatorial
        Algorithms, 1st ed. Reading, MA: Addison-Wesley, 2011.
    (3) Wendy Myrvold and Frank Ruskey. 2001. Ranking and unranking
        permutations in linear time. Inf. Process. Lett. 79, 6 (September 2001),
        281-284. DOI=10.1016/S0020-0190(01)00141-7
    (4) D. L. Kreher, D. R. Stinson 'Combinatorial Algorithms'
        CRC Press, 1999
    (5) Graham, R. L.; Knuth, D. E.; and Patashnik, O.
        Concrete Mathematics: A Foundation for Computer Science, 2nd ed.
        Reading, MA: Addison-Wesley, 1994.
    (6) http://en.wikipedia.org/wiki/Permutation#Product_and_inverse

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

        Examples
        ========

        >>> from sympy.combinatorics.permutations import Permutation
        >>> p = Permutation([[2,0],[3,1]])
        >>> p.array_form
        [2, 3, 0, 1]
        >>> Permutation([[2,0,3,1]]).array_form
        [3, 2, 0, 1]
        >>> Permutation([2,0,3,1]).array_form
        [2, 0, 3, 1]

        See Also
        ========
        cyclic_form
        """
        if self._array_form is not None:
            return self._array_form
        if not isinstance(self.args[0][0], list):
            self._array_form = self.args[0]
            return self._array_form
        size = 0
        cycles = self.args[0]
        for c in cycles:
            size += len(c)
        perm = [None]*size
        for c in cycles:
            for i in range(len(c)-1):
                perm[c[i]] = c[i+1]
            perm[c[-1]] = c[0]
        self._array_form = perm
        return perm

    @property
    def cyclic_form(self):
        """
        This is used to convert to the cyclic notation
        from the canonical notation.

        Examples
        ========

        >>> from sympy.combinatorics.permutations import Permutation
        >>> p = Permutation([0,3,1,2])
        >>> p.cyclic_form
        [[1, 3, 2], [0]]

        See Also
        ========
        array_form
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
        """
        Returns the number of numbers in the permutation

        Examples
        ========

        >>> from sympy.combinatorics import Permutation
        >>> Permutation([[3, 2], [0, 1]]).size
        4
        """
        return len(self.array_form)

    def __new__(cls, *args, **kw_args):
        """
        Constructor for the Permutation object.

        Examples
        ========

        >>> from sympy.combinatorics.permutations import Permutation
        >>> p = Permutation([0,1,2])
        >>> p
        Permutation([0, 1, 2])
        >>> q = Permutation([[0,1],[2]])
        >>> q
        Permutation([[0, 1], [2]])
        """
        if not args or not is_sequence(args[0]) or len(args) > 1 or \
           len(set(is_sequence(a) for a in args[0])) > 1:
            raise ValueError('Permutation argument must be a list of ints or a list of lists.')

        # 0, 1, ..., n-1 should all be present
        temp = [int(i) for i in flatten(args[0])]
        if set(range(len(temp))) != set(temp):
            raise ValueError("Integers 0 through %s must be present." % len(temp))

        cform = aform = None
        if args[0] and is_sequence(args[0][0]):
            cform = [list(a) for a in args[0]]
        else:
            aform = list(args[0])

        ret_obj = Basic.__new__(cls, (cform or aform), **kw_args)
        ret_obj._cyclic_form, ret_obj._array_form = cform, aform
        return ret_obj


    def __add__(self, other):
        """
        Routine for addition of permutations.

        This is defined in terms of the Lehmer code of a
        permutation. The Lehmer code is nothing but the
        inversion vector of a permutation. In this scheme
        the identity permutation is like a zero element.
        # TODO add a reference

        Examples
        ========

        >>> from sympy.combinatorics.permutations import Permutation
        >>> p = Permutation([0,1,2,3])
        >>> q = Permutation([2,1,3,0])
        >>> p+q == q
        True
        >>> from sympy.combinatorics.permutations import Permutation
        >>> a = Permutation([0, 3, 1, 2])
        >>> b = Permutation([2, 1, 0, 3])
        >>> a+b
        Permutation([2, 0, 1, 3])
        """
        n = self.size
        if n != other.size:
            raise ValueError("The permutations must be of equal size.")
        a = self.inversion_vector()
        b = other.inversion_vector()
        result_inv = [(a[i] + b[i]) % (n - i) for i in xrange(n - 1)]
        return Permutation.from_inversion_vector(result_inv)

    def __sub__(self, other):
        """
        Routine for subtraction of permutations.

        The idea behind this is the same as the above.
        With subtraction however,

        Examples
        ========

        >>> from sympy.combinatorics.permutations import Permutation
        >>> p = Permutation([0,1,2,3])
        >>> q = Permutation([2,1,3,0])
        >>> q-p==q
        True
        """
        n = self.size
        if n != other.size:
            raise ValueError("The permutations must be of equal size.")
        a = self.inversion_vector()
        b = other.inversion_vector()
        result_inv = [(a[i] - b[i]) % (n - i) for i in xrange(n - 1)]
        return Permutation.from_inversion_vector(result_inv)

    def __mul__(self, other):
        """
        Routine for multiplication of permutations.

        Examples
        ========

        >>> from sympy.combinatorics.permutations import Permutation
        >>> p = Permutation([1,2,3,0])
        >>> q = Permutation([3,2,0,1])
        >>> p*q
        Permutation([0, 3, 1, 2])

        If one of the permutations is in a cyclic form then it is first
        converted to an array form and then multiplied. ::

            >>> q = Permutation([[1,3,2],[0]])
            >>> p*q
            Permutation([1, 0, 2, 3])

        """
        a = self.array_form
        b = other.array_form
        if len(a) != len(b):
            raise ValueError("The number of elements in the permutations \
don\'t match.")

        perm = [a[i] for i in b]
        return _new_from_array_form(perm)

    def commutes_with(self, other):
        """
        Checks if the elements are commuting.

        Examples
        ========

        >>> from sympy.combinatorics.permutations import Permutation
        >>> a = Permutation([1,4,3,0,2,5])
        >>> b = Permutation([0,1,2,3,4,5])
        >>> a.commutes_with(b)
        True
        >>> b = Permutation([2,3,5,4,1,0])
        >>> a.commutes_with(b)
        False

        """
        a = self.array_form
        b = other.array_form
        if len(a) != len(b):
            raise ValueError("The number of elements in the permutations \
don\'t match.")
        for i in range(len(a)-1):
            if a[b[i]] != b[a[i]]:
                return False
        return True

    def __pow__(self, n):
        """
        Routine for finding powers of a permutation.

        Examples
        ========

        >>> from sympy.combinatorics.permutations import Permutation
        >>> p = Permutation([2,0,3,1])
        >>> p**4
        Permutation([0, 1, 2, 3])
        """
        n = int(n)
        if n == 0:
            return Permutation(range(self.size))
        if n < 0:
            return pow(~self, -n)
        a = self.array_form
        if n == 2:
            b = [a[i] for i in a]
        elif n == 3:
            b = [a[a[i]] for i in a]
        elif n == 4:
            b = [a[a[a[i]]] for i in a]
        else:
            b = range(len(a))
            while 1:
                if n&1:
                    b = [b[i] for i in a]
                    n -= 1
                    if not n:
                        break
                if n%4 == 0:
                    a = [a[a[a[i]]] for i in a]
                    n = n // 4
                elif n%2 == 0:
                    a = [a[i] for i in a]
                    n = n // 2
        return _new_from_array_form(b)

    def transpositions(self):
        """a list of transpositions representing the permutation

        Examples
        ========

        >>> from sympy.combinatorics.permutations import Permutation
        >>> p = Permutation([[1,2,3],[0,4,5,6,7]])
        >>> p.transpositions()
        [(1, 3), (1, 2), (0, 7), (0, 6), (0, 5), (0, 4)]
        """
        a = self.cyclic_form
        res = []
        for x in a:
            nx = len(x)
            if nx == 2:
                res.append(tuple(x))
            elif nx > 2:
                first = x[0]
                for y in x[nx-1:0:-1]:
                    res.append((first,y))
        return res

    def __invert__(self):
        """
        Finds the invert of a permutation.

        A permutation multiplied by its invert equals
        the identity permutation.

        Examples
        ========

        >>> from sympy.combinatorics.permutations import Permutation
        >>> p = Permutation([[2,0],[3,1]])
        >>> ~p
        Permutation([2, 3, 0, 1])
        >>> p*(~p) == Permutation([0,1,2,3])
        True
        """
        a = self.array_form
        n = len(a)
        inv_form = [0] * n
        for i in xrange(n):
            inv_form[a[i]] = i
        return _new_from_array_form(inv_form)

    def __call__(self, arg):
        """
        Allows applying a permutation instance as a bijective function.

        Examples
        ========

        >>> from sympy.combinatorics.permutations import Permutation
        >>> p = Permutation([[2,0],[3,1]])
        >>> p(3)
        1
        """
        if not isinstance(arg, int):
            raise ValueError("Arguments must be integers.")
        return self.array_form[arg]

    def atoms(self):
        """
        Returns all the elements of a permutation

        Examples
        ========

        >>> from sympy.combinatorics import Permutation
        >>> Permutation([0, 1, 2, 3, 4, 5]).atoms()
        set([0, 1, 2, 3, 4, 5])
        >>> Permutation([[0, 1], [2, 3], [4, 5]]).atoms()
        set([0, 1, 2, 3, 4, 5])
        """
        return set(self.array_form)

    def next_lex(self):
        """
        Returns the next permutation in lexicographical order.
        If self is the last permutation in lexicographical order
        it returns None.
        See [4] section 2.4.


        Examples
        ========

        >>> from sympy.combinatorics.permutations import Permutation
        >>> p = Permutation([2, 3, 1, 0])
        >>> p = Permutation([2, 3, 1, 0]); p.rank()
        17
        >>> p = p.next_lex(); p.rank()
        18

        See Also
        ========
        rank, unrank_lex
        """
        perm = self.array_form[:]
        n = len(perm)
        i = n - 2
        while perm[i+1] < perm[i]:
            i -= 1
        if i == -1:
            return None
        else:
            j = n - 1
            while perm[j] < perm[i]:
                j -= 1
            perm[j], perm[i] = perm[i], perm[j]
            i += 1
            j = n - 1
            while i < j:
                perm[j], perm[i] = perm[i], perm[j]
                i += 1
                j -= 1
        return _new_from_array_form(perm)

    @classmethod
    def unrank_nonlex(self, n, r):
        """
        This is a linear time unranking algorithm that does not
        respect lexicographic order [3].

        Examples
        ========

        >>> from sympy.combinatorics.permutations import Permutation
        >>> Permutation.unrank_nonlex(4, 5)
        Permutation([2, 0, 3, 1])
        >>> Permutation.unrank_nonlex(4, -1)
        Permutation([0, 1, 2, 3])

        See Also
        ========
        next_nonlex, rank_nonlex
        """
        def _unrank1(n, r, a):
            if n > 0:
                a[n - 1], a[r % n] = a[r % n], a[n - 1]
                _unrank1(n - 1, r//n, a)

        id_perm = range(n)
        n = int(n)
        r = r % ifac(n)
        _unrank1(n, r, id_perm)
        return _new_from_array_form(id_perm)

    def rank_nonlex(self, inv_perm = None):
        """
        This is a linear time ranking algorithm that does not
        enforce lexicographic order [3].


        Examples
        ========

        >>> from sympy.combinatorics.permutations import Permutation
        >>> p = Permutation([0,1,2,3])
        >>> p.rank_nonlex()
        23

        See Also
        ========
        next_nonlex, unrank_nonlex
        """
        def _rank1(n, perm, inv_perm):
            if n == 1:
                return 0
            s = perm[n - 1]
            t = inv_perm[n - 1]
            perm[n - 1], perm[t] = perm[t], s
            inv_perm[n - 1], inv_perm[s] = inv_perm[s], t
            return s + n*_rank1(n - 1, perm, inv_perm)

        if inv_perm is None:
            inv_perm = (~self).array_form
        if not inv_perm:
            return 0
        perm = self.array_form[:]
        r = _rank1(len(perm), perm, inv_perm)
        return r

    def next_nonlex(self):
        """
        Returns the next permutation in nonlex order [3].
        If self is the last permutation in this order it returns None.

        Examples
        ========

        >>> from sympy.combinatorics.permutations import Permutation
        >>> p = Permutation([2, 0, 3, 1]); p.rank_nonlex()
        5
        >>> p = p.next_nonlex(); p
        Permutation([3, 0, 1, 2])
        >>> p.rank_nonlex()
        6

        See Also
        ========
        rank_nonlex, unrank_nonlex
        """
        r = self.rank_nonlex()
        if r == ifac(self.size) - 1:
            return None
        return Permutation.unrank_nonlex(self.size, r+1)

    def rank(self):
        """
        Returns the lexicographic rank of the permutation.

        Examples
        ========

        >>> from sympy.combinatorics.permutations import Permutation
        >>> p = Permutation([0,1,2,3])
        >>> p.rank()
        0
        >>> p = Permutation([3,2,1,0])
        >>> p.rank()
        23

        See Also
        ========
        next_lex, unrank_lex
        """
        rank = 0
        rho = self.array_form[:]
        n = self.size - 1
        size = n + 1
        psize = int(ifac(n))
        for j in xrange(size - 1):
            rank += rho[j]*psize
            for i in xrange(j + 1, size):
                if rho[i] > rho[j]:
                    rho[i] -= 1
            psize //= n
            n -= 1
        return rank

    @property
    def cardinality(self):
        """
        Returns the number of all possible permutations.

        Examples
        ========

        >>> from sympy.combinatorics.permutations import Permutation
        >>> p = Permutation([0,1,2,3])
        >>> p.cardinality
        24
        """
        return int(ifac(self.size))

    def parity(self):
        """
        Computes the parity of a permutation.

        The parity of a permutation reflects the parity of the
        number of inversions in the permutation, i.e., the
        number of pairs of x and y such that x > y but p[x] < p[y].

        Examples
        ========

        >>> from sympy.combinatorics.permutations import Permutation
        >>> p = Permutation([0,1,2,3])
        >>> p.parity()
        0
        >>> p = Permutation([3,2,0,1])
        >>> p.parity()
        1

        See Also
        ========
        perm_af_parity
        """
        if self._cyclic_form is not None:
            return (self.size - len(self._cyclic_form)) % 2

        return perm_af_parity(self.array_form)

    @property
    def is_even(self):
        """
        Checks if a permutation is even.

        Examples
        ========

        >>> from sympy.combinatorics.permutations import Permutation
        >>> p = Permutation([0,1,2,3])
        >>> p.is_even
        True
        >>> p = Permutation([3,2,1,0])
        >>> p.is_even
        True

        See Also
        ========
        is_odd
        """
        return S(self.parity()).is_even
    @property
    def is_odd(self):
        """
        Checks if a permutation is odd.

        Examples
        ========

        >>> from sympy.combinatorics.permutations import Permutation
        >>> p = Permutation([0,1,2,3])
        >>> p.is_odd
        False
        >>> p = Permutation([3,2,0,1])
        >>> p.is_odd
        True

        See Also
        ========
        is_even
        """
        return S(self.parity()).is_odd

    @property
    def is_Singleton(self):
        """
        Checks to see if the permutation contains only one number
        Therefore there is only one possible permutation of this set of numbers

        Examples
        ========

        >>> from sympy.combinatorics import Permutation
        >>> Permutation([0]).is_Singleton
        True
        >>> Permutation([0, 1]).is_Singleton
        False

        See Also
        ========
        is_Empty
        """
        return self.size == 1

    @property
    def is_Empty(self):
        """
        Checks to see if the permutation is a set with zero elements

        Examples
        ========

        >>> from sympy.combinatorics import Permutation
        >>> Permutation([]).is_Empty
        True
        >>> Permutation([0]).is_Empty
        False

        See Also
        ========
        is_Singleton
        """
        return self.size == 0

    @property
    def is_Identity(self):
        """
        Returns True if the Permutation is an identity permutation.

        Examples
        ========
        >>> from sympy.combinatorics.permutations import Permutation
        >>> p = Permutation([[0],[1],[2]])
        >>> p.is_Identity
        True
        >>> p = Permutation([0,1,2])
        >>> p.is_Identity
        True
        >>> p = Permutation([0,2,1])
        >>> p.is_Identity
        False

        See Also
        ========
        order
        """
        if self._cyclic_form:
            return self.size == len(self._cyclic_form)
        a = self.array_form
        return a == range(len(a))

    def ascents(self):
        """
        Returns the positions of ascents in a permutation, ie, the location
        where p[i] < p[i+1]

        Examples
        ========

        >>> from sympy.combinatorics.permutations import Permutation
        >>> p = Permutation([4,0,1,3,2])
        >>> p.ascents()
        [1, 2]

        See Also
        ========
        descents, inversions, min, max
        """
        a = self.array_form
        pos = [i for i in xrange(len(a)-1) if a[i] < a[i+1]]
        return pos

    def descents(self):
        """
        Returns the positions of descents in a permutation, ie, the location
        where p[i] > p[i+1]

        Examples
        ========

        >>> from sympy.combinatorics.permutations import Permutation
        >>> p = Permutation([4,0,1,3,2])
        >>> p.descents()
        [0, 3]

        See Also
        ========
        ascents, inversions, min, max
        """
        a = self.array_form
        pos = [i for i in xrange(len(a)-1) if a[i] > a[i+1]]
        return pos

    def max(self):
        """
        The maximum element moved by the permutation.

        Examples
        ========

        >>> from sympy.combinatorics.permutations import Permutation
        >>> p = Permutation([1,0,2,3,4])
        >>> p.max()
        1

        See Also
        ========
        min, descents, ascents, inversions
        """
        max = 0
        a = self.array_form
        for i in xrange(len(a)):
            if a[i] != i and a[i] > max:
                max = a[i]
        return max

    def min(self):
        """
        The minimum element moved by the permutation.

        Examples
        ========

        >>> from sympy.combinatorics.permutations import Permutation
        >>> p = Permutation([0,1,4,3,2])
        >>> p.min()
        2

        See Also
        ========
        max, descents, ascents, inversions
        """
        a = self.array_form
        min = len(a)
        for i in xrange(len(a)):
            if a[i] != i and a[i] < min:
                min = a[i]
        return min

    def inversions(self):
        """
        Computes the number of inversions of a permutation.

        An inversion is where i > j but p[i] < p[j].

        Examples
        ========

        >>> from sympy.combinatorics.permutations import Permutation
        >>> p = Permutation([0,1,2,3,4,5])
        >>> p.inversions()
        0

        See Also
        ========
        descents, ascents, min, max
        """
        inversions = 0
        a = self.array_form
        n = len(a)
        for i in xrange(n - 1):
            b = a[i]
            for j in xrange(i + 1, n):
                if b > a[j]:
                    inversions += 1
        return inversions

    def conjugate(self, x):
        """
        Computes the conjugate Permutation `~x*p*x'

        Examples
        ========

        >>> from sympy.combinatorics.permutations import Permutation
        >>> a = Permutation([0,2,1,3])
        >>> b = Permutation([0,2,3,1])
        >>> a.conjugate(b)
        Permutation([0, 3, 2, 1])
        >>> ~b*a*b
        Permutation([0, 3, 2, 1])
        """

        a = self.array_form
        b = x.array_form
        n = len(a)
        if len(b) != n:
            raise ValueError("The number of elements in the permutations \
don\'t match.")
        invb = [None]*n
        for i in xrange(n):
            invb[b[i]] = i
        return _new_from_array_form([invb[a[i]] for i in b])

    def signature(self):
        """
        Gives the signature of the permutation needed to place the
        elements of the permutation in canonical order.

        The signature is calculated as (-1)^<# no. of inversions>

        Examples
        ========

        >>> from sympy.combinatorics.permutations import Permutation
        >>> p = Permutation([0,1,2])
        >>> p.signature()
        1
        >>> q = Permutation([0,2,1])
        >>> q.signature()
        -1

        See Also
        ========
        inversions
        """
        return (-1)**self.inversions()

    def order(self):
        """
        Computes the order of a permutation.

        When the permutation is raised to the power of its
        order it equals the identity permutation.

        Examples
        ========

        >>> from sympy.combinatorics.permutations import Permutation
        >>> p = Permutation([3,1,5,2,4,0])
        >>> p.order()
        4
        >>> (p**(p.order()))
        Permutation([0, 1, 2, 3, 4, 5])

        See Also
        ========
        identity
        """
        return reduce(lcm,[1]+[len(cycle) for cycle in self.cyclic_form])

    def length(self):
        """
        Returns the number of integers moved by a permutation.

        Examples
        ========

        >>> from sympy.combinatorics import Permutation
        >>> Permutation([0, 3, 2, 1]).length()
        2
        >>> Permutation([[0, 1], [2, 3]]).length()
        4

        See Also
        ========
        min, max
        """
        length = 0
        a = self.array_form
        for i in xrange(len(a)):
            if a[i] != i:
                length += 1
        return length

    @property
    def is_Positive(self):
        """
        Checks if the permutation is positive

        Examples
        ========

        >>> from sympy.combinatorics import Permutation
        >>> Permutation([0, 1, 2]).is_Positive
        True

        See Also
        ========
        is_Negative
        """
        return self.signature() > 0

    @property
    def is_Negative(self):
        """
        Checks if the permutation is negative

        Examples
        ========

        >>> from sympy.combinatorics import Permutation
        >>> Permutation([0, 1, 2]).is_Negative
        False

        See Also
        ========
        is_Positive
        """
        return self.signature() < 0

    @property
    def cycles(self):
        """
        Returns the number of cycles that the permutation
        has been decomposed into.

        Examples
        ========

        >>> from sympy.combinatorics import Permutation
        >>> Permutation([0, 1, 2]).cycles
        3
        >>> Permutation([[0, 1], [2, 3]]).cycles
        2
        """
        return len(self.cyclic_form)

    def index(self):
        """
        Returns the index of a permutation.

        The index of a permutation is the sum of all
        subscripts j such that p[j] is greater than
        p[j+1].

        Examples
        ========

        >>> from sympy.combinatorics.permutations import Permutation
        >>> p = Permutation([3, 0, 2, 1, 4])
        >>> p.index()
        2
        """
        a = self.array_form

        return sum([j for j in xrange(len(a) - 1) if a[j] > a[j+1]])

    def runs(self):
        """
        Returns the runs of a permutation.

        An ascending sequence in a permutation is called a run [5]


        Examples
        ========

        >>> from sympy.combinatorics.permutations import Permutation
        >>> p = Permutation([2,5,7,3,6,0,1,4,8])
        >>> p.runs()
        [[2, 5, 7], [3, 6], [0, 1, 4, 8]]
        >>> q = Permutation([1,3,2,0])
        >>> q.runs()
        [[1, 3], [2], [0]]
        """
        cycles = []
        temp_cycle = []
        perm = self.array_form
        for i in xrange(len(perm) - 1):
            current_elem = perm[i]
            next_elem    = perm[i+1]

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

    def inversion_vector(self):
        """
        Gets the inversion vector of the permutation.

        The inversion vector consists of elements whose value
        indicates the number of elements in the permutation
        that are lesser than it and lie on its right hand side.

        Examples
        ========

        >>> from sympy.combinatorics.permutations import Permutation
        >>> p = Permutation([4,8,0,7,1,5,3,6,2])
        >>> p.inversion_vector()
        [4, 7, 0, 5, 0, 2, 1, 1]
        >>> p = Permutation([3,2,1,0])
        >>> p.inversion_vector()
        [3, 2, 1]
        """
        self_array_form = self.array_form
        n = len(self_array_form)
        inversion_vector = [0] * (n - 1)

        for i in xrange(n - 1):
            val = 0
            for j in xrange(i+1, n):
                if self_array_form[j] < self_array_form[i]:
                    val += 1
            inversion_vector[i] = val
        return inversion_vector

    def rank_trotterjohnson(self):
        """
        Returns the Trotter Johnson rank, which we get from the minimal
        change algorithm. See [4] section 2.4.

        Examples
        ========

        >>> from sympy.combinatorics.permutations import Permutation
        >>> p = Permutation([0,1,2,3])
        >>> p.rank_trotterjohnson()
        0
        >>> p = Permutation([0,2,1,3])
        >>> p.rank_trotterjohnson()
        7

        See Also
        ========
        unrank_trotterjohnson, next_trotterjohnson
        """
        if self.array_form == [] or self.is_Identity:
            return 0
        if self.array_form == [1, 0]:
            return 1
        perm = self.array_form
        n = self.size
        rank = 0
        for j in range(1, n):
            k = 1
            i = 0
            while perm[i] != j:
                if perm[i] < j:
                    k += 1
                i += 1
            j1 = j + 1
            if rank % 2 == 0:
                rank = j1*rank + j1 - k
            else:
                rank = j1*rank + k - 1
        return rank

    @classmethod
    def unrank_trotterjohnson(self, size, rank):
        """
        Trotter Johnson permutation unranking. See [4] section 2.4.

        Examples
        ========

        >>> from sympy.combinatorics.permutations import Permutation
        >>> Permutation.unrank_trotterjohnson(5,10)
        Permutation([0, 3, 1, 2, 4])

        See Also
        ========
        rank_trotterjohnson, next_trotterjohnson
        """
        perm = [0]*size
        r2 = 0
        n = ifac(size)
        pj = 1
        for j in range(2, size+1):
            pj *= j
            r1 = (rank * pj) // n
            k = r1 - j*r2
            if r2 % 2 == 0:
                for i in range(j-1, j-k-1, -1):
                    perm[i] = perm[i-1]
                perm[j-k-1] = j-1
            else:
                for i in range(j-1, k, -1):
                    perm[i] = perm[i-1]
                perm[k] = j-1
            r2 = r1
        return _new_from_array_form(perm)

    def next_trotterjohnson(self):
        """
        Returns the next permutation in Trotter-Johnson order.
        If self is the last permutation it returns None.
        See [4] section 2.4.

        Examples
        ========

        >>> from sympy.combinatorics.permutations import Permutation
        >>> p = Permutation([3, 0, 2, 1])
        >>> p.rank_trotterjohnson()
        4
        >>> p = p.next_trotterjohnson(); p
        Permutation([0, 3, 2, 1])
        >>> p.rank_trotterjohnson()
        5

        See Also
        ========
        rank_trotterjohnson, unrank_trotterjohnson
        """
        pi = self.array_form[:]
        n = len(pi)
        st = 0
        rho = pi[:]
        done = False
        m = n-1
        while m > 0 and not done:
            d = rho.index(m)
            for i in range(d, m):
                rho[i] = rho[i+1]
            par = perm_af_parity(rho[:m])
            if par == 1:
                if d == m:
                    m -= 1
                else:
                    pi[st+d], pi[st+d+1] = pi[st+d+1], pi[st+d]
                    done = True
            else:
                if d == 0:
                    m -= 1
                    st += 1
                else:
                    pi[st+d], pi[st+d-1] = pi[st+d-1], pi[st+d]
                    done = True
        if m == 0:
            return None
        return _new_from_array_form(pi)

    def get_precedence_matrix(self):
        """
        Gets the precedence matrix. This is used for computing the
        distance between two permutations.

        Examples
        ========

        >>> from sympy.combinatorics.permutations import Permutation
        >>> p = Permutation.josephus(3,6,1)
        >>> p
        Permutation([2, 5, 3, 1, 4, 0])
        >>> p.get_precedence_matrix()
        [0, 0, 0, 0, 0, 0]
        [1, 0, 0, 0, 1, 0]
        [1, 1, 0, 1, 1, 1]
        [1, 1, 0, 0, 1, 0]
        [1, 0, 0, 0, 0, 0]
        [1, 1, 0, 1, 1, 0]

        See Also
        ========
        get_precedence_distance, get_adjacency_matrix, get_adjacency_distance
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

        Examples
        ========

        >>> from sympy.combinatorics.permutations import Permutation
        >>> p = Permutation([2, 0, 4, 3, 1])
        >>> q = Permutation([3, 1, 2, 4, 0])
        >>> p.get_precedence_distance(q)
        7
        >>> q.get_precedence_distance(p)
        7

        See Also
        ========
        get_precedence_matrix, get_adjacency_matrix, get_adjacency_distance
        """
        if self.size != other.size:
            raise ValueError("The permutations must be of the same size.")
        self_prec_mat = self.get_precedence_matrix()
        other_prec_mat = other.get_precedence_matrix()
        n_prec = 0
        for i in xrange(self.size):
            for j in xrange(self.size):
                if i == j:
                    continue
                if self_prec_mat[i, j] * other_prec_mat[i, j] == 1:
                    n_prec += 1
        d = self.size * (self.size - 1)//2 - n_prec
        return d

    def get_adjacency_matrix(self):
        """
        Computes the adjacency matrix of a permutation.

        If job i is adjacent to job j in a permutation p
        then we set m[i, j] = 1 where m is the adjacency
        matrix of p.

        Examples
        ========

        >>> from sympy.combinatorics.permutations import Permutation
        >>> p = Permutation.josephus(3,6,1)
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

        See Also
        ========
        get_precedence_matrix, get_precedence_distance, get_adjacency_distance
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


        Examples
        ========

        >>> from sympy.combinatorics.permutations import Permutation
        >>> p = Permutation([0, 3, 1, 2, 4])
        >>> q = Permutation.josephus(4, 5, 2)
        >>> p.get_adjacency_distance(q)
        3
        >>> r = Permutation([0, 2, 1, 4, 3])
        >>> p.get_adjacency_distance(r)
        4

        See Also
        ========
        get_precedence_matrix, get_precedence_distance, get_adjacency_matrix
        """
        if self.size != other.size:
            raise ValueError("The permutations must be of the same size.")
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

        Examples
        ========

        >>> from sympy.combinatorics.permutations import Permutation
        >>> p = Permutation([0, 3, 1, 2, 4])
        >>> q = Permutation.josephus(4, 5, 2)
        >>> r = Permutation([3, 1, 4, 0, 2])
        >>> p.get_positional_distance(q)
        12
        >>> p.get_positional_distance(r)
        12

        See Also
        ========
        get_precedence_distance, get_adjacency_distance
        """
        a = self.array_form
        b = other.array_form
        if len(a) != len(b):
            raise ValueError("The permutations must be of the same size.")
        return sum([abs(a[i] - b[i]) for i in xrange(len(a))])

    @classmethod
    def josephus(self, m, n, s = 1):
        """
        Computes the Josephus permutation for a given number of
        prisoners, frequency of removal and desired number of
        survivors.

        There are people standing in a circle waiting to be executed.
        After the first person is executed, certain number of people
        are skipped and another person is executed. Then again, people
        are skipped and a person is executed. The elimination proceeds
        around the circle (which is becoming smaller and smaller as the
        executed people are removed), until only the last person
        remains, who is given freedom.

        References:
        [1] http://en.wikipedia.org/wiki/Flavius_Josephus
        [2] http://en.wikipedia.org/wiki/Josephus_problem

        Examples
        ========

        >>> from sympy.combinatorics.permutations import Permutation
        >>> Permutation.josephus(3, 40, 1)
        Permutation([2, 5, 8, 11, 14, 17, 20, 23, 26, 29, 32, 35, \
        38, 1, 6, 10, 15, 19, 24, 28, 33, 37, 3, 9, 16, 22, 30, 36, \
        4, 13, 25, 34, 7, 21, 39, 18, 0, 31, 12, 27])
        """
        from collections import deque
        m -= 1
        if s <= 0:
            s = 1
        Q = deque(range(n))
        perm = []
        while len(Q) > s:
            for dp in xrange(m):
                Q.append(Q.popleft())
            perm.append(Q.popleft())
        perm.extend(list(Q))
        return Permutation(perm)

    @classmethod
    def from_inversion_vector(self, inversion):
        """
        Calculates the permutation from the inversion
        vector.

        Examples
        ========

        >>> from sympy.combinatorics.permutations import Permutation
        >>> Permutation.from_inversion_vector([3,2,1,0,0])
        Permutation([3, 2, 1, 0, 4, 5])

        """
        size = len(inversion) + 1
        N = [i for i in xrange(size)]
        perm = []
        try:
            for k in xrange(size - 1):
                val = N[inversion[k]]
                perm.append(val)
                N.remove(val)
        except IndexError:
            raise ValueError("The inversion vector is not valid.")
        perm.extend(N)
        return _new_from_array_form(perm)

    @classmethod
    def random_permutation(self, n):
        """
        Generates a random permutation.

        Uses the underlying Python psuedo-random
        number generator.

        Examples
        ========

        >>> from sympy.combinatorics.permutations import Permutation
        >>> a = Permutation.random_permutation(5)
        >>> (a*(~a)).is_Identity
        True

        """
        perm_array = range(n)
        random.shuffle(perm_array)
        return _new_from_array_form(perm_array)

    @classmethod
    def unrank_lex(self, size, rank):
        """
        Lexicographic permutation unranking.

        Examples
        ========

        >>> from sympy.combinatorics.permutations import Permutation
        >>> a = Permutation.unrank_lex(5,10)
        >>> a.rank()
        10
        >>> a
        Permutation([0, 2, 4, 1, 3])

        See Also
        ========
        rank, next_lex
        """
        perm_array = [0] * size
        psize = 1
        for i in xrange(size):
            new_psize = psize*(i + 1)
            d = (rank % new_psize) // psize
            rank -= d*psize
            perm_array[size - i - 1] = d
            for j in xrange(size - i, size):
                if perm_array[j] > d-1:
                    perm_array[j] += 1
            psize = new_psize
        return _new_from_array_form(perm_array)

def _new_from_array_form(perm):
    """
    factory function to produce a Permutation object from a list;
    the list is bound to the _array_form attribute, so it must
    not be modified; this method is meant for internal use only;
    the list a is supposed to be generated as a temporary value in a method,
    so p = _new_from_array_form(a) is the only object to hold a
    reference to a::

        >>> from sympy.combinatorics.permutations import _new_from_array_form
        >>> a = [2,1,3,0]
        >>> p = _new_from_array_form(a)
        >>> p
        Permutation([2, 1, 3, 0])

    """
    p = Basic.__new__(Permutation, perm)
    p._array_form = perm
    return p
