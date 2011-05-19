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
