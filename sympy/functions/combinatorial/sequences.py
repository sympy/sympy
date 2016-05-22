"""
This module implements some special sequences that commonly appear in
combinatorial contexts. The functions return the initial segment of
given size for these sequences.

The need for such a module arises from the fact that sequences often
can be computed much more efficiently (for example via recurrences)
than by list comprehension.

A convincing example is the computation of the Euler numbers.
To compute the first 500 nonzero Euler numbers with 'euler_sequence'
is more than 1000 times faster than with list comprehension of 'euler'.
(Benchmark with SymPy 1.0 and Python 3.5.)

Moreover the typical use case often requires the initial segment of
of sequences. For example Bernoulli numbers rarely come alone;
their typical use is in the Taylor-MacLaurin expansion of functions
which needs B0, B2, B4,..., B2n up to some n.

Contrary to the 'bernoulli' and 'euler' functions the sequences here
are not depended on the mpmath module.

Currently the following sequences are implemented:
    * andre_sequence
    * bell_sequence
    * bernoulli_sequence
    * euler_sequence
    * tangent_sequence

"""

from __future__ import print_function, division

from sympy import sympify
from sympy.core import Rational, Integer
from sympy.core.function import Function


class bernoulli_sequence(Function):
    """
        Input
        =====

        - size -- positive integer

        Output
        ======

        If 0 < size then return the list of the even Bernoulli numbers

            [B0, B2, B4, ..., Bn] where n = 2*size-2,

        else return [].

        Is equivalent to [bernoulli(n) for n in range(0, size, 2)].
        Generally 'bernoulli_number_list' is (much) more efficient than
        computing such a list by list comprehension and 'bernoulli'.

        The algorithm used is Ludwig Seidel's 'boustrophedon transform'.

        Examples
        ========

        >>> from sympy.functions.combinatorial.sequences import bernoulli_sequence
        >>> bernoulli_sequence(0)
        []

        >>> bernoulli_sequence(7)
        [1, 1/6, -1/30, 1/42, -1/30, 5/66, -691/2730]

        References
        ==========

        .. [1] Ludwig Seidel, Ueber eine einfache Entstehungsweise der
               Bernoulli'schen Zahlen und einiger verwandten Reihen.
        .. [2] Peter Luschny, Computation and asymptotics of
               Bernoulli numbers, http://goo.gl/AwskX9
        .. [3] OEIS Bernoulli numbers B_2n
               [ http://oeis.org/A000367 / http://oeis.org/A002445 ]
    """

    @classmethod
    def eval(cls, size):
        s = sympify(size)
        if not s.is_Integer or s <= Integer(0):
            return []

        h, p, q, s, b = 1, 1, 1, -2, True
        D = [0] * (size + 2)
        D[1] = 1
        R = [Rational(1, 1)]

        for i in range(2 * size - 2):
            if b:
                h += 1
                p *= 4
                s = -s
                q = s * (p - 1)
                for k in range(h, 0, -1):
                    D[k] += D[k + 1]
            else:
                for k in range(1, h, 1):
                    D[k] += D[k - 1]
                R.append(Rational(D[h - 1], q))
            b = not b
        return R


class andre_sequence(Function):
    """
        Input
        =====

        - size -- positive integer

        Output
        ======

        If 0 < size then return the list of Andre numbers

            [A0, A1, A2, ..., A{n-1}]

        else return the empty list [].

        The Andre numbers are defined by OEIS A000111. At even indices
        these are the secant numbers OEIS A000364 and at odd indices
        the tangent numbers OEIS A000182.

        Andre numbers count the number of alternating permutations.
        The Andre numbers are named after Desire Andre who studied
        2-alternating permutations in 1879 ("Developpements de sec x
        et tg x") and in 1881 ("Memoire sur les permutations alternees").
        They are also called Euler or up/down numbers.

        The Andre numbers can be expressed as values of the polygamma
        function for n>0:

            A(n) = 2 * I**(n + 1) * polylog(-n, -I).

        Analytically the A(n) are the coefficients in the expansion
        of sec(x) + tan(x).

        The Andre numbers can also be seen as the values of the
        Andre polynomials A(n,x) at x = 1.

        The algorithm used is Ludwig Seidel's 'boustrophedon transform'.

        Examples
        ========

        >>> from sympy.functions.combinatorial.sequences import andre_sequence
        >>> andre_sequence(0)
        []

        >>> a = andre_sequence(12)
        >>> print(a)  # numbers of alternating permutations
        [1, 1, 1, 2, 5, 16, 61, 272, 1385, 7936, 50521, 353792]

        >>> print(a[::2])  # secant numbers
        [1, 1, 5, 61, 1385, 50521]

        >>> print(a[1::2])  # tangent numbers
        [1, 2, 16, 272, 7936, 353792]

        References
        ==========

        .. [1] Ludwig Seidel, Ueber eine einfache Entstehungsweise der
               Bernoulli'schen Zahlen und einiger verwandten Reihen.
        .. [2] Peter Luschny, The Seidel transform. http://goo.gl/IxYZl0
        .. [3] OEIS A000111, Euler up/down numbers. http://oeis.org/A000111
        .. [4] OEIS A094503, Andre polynomials. http://oeis.org/A094503
    """

    @classmethod
    def eval(cls, size):
        s = sympify(size)
        if not s.is_Integer or s <= Integer(0):
            return []

        k, e, R = 0, 1, []
        A = {-1: 0, 0: 1}
        for i in range(size):
            ak = 0
            A[k + e] = 0
            e = -e
            for _ in range(i + 1):
                ak += A[k]
                A[k] = ak
                k += e
            R.append(ak)
        return R


class euler_sequence(Function):
    """
        Input
        =====

        - size -- positive integer

        Output
        ======

        If 0 < size then return the list of nonzero secant numbers
        which is also the list of the even unsigned Euler numbers

            [E0, E2, E4, ..., E{2size-2}]

        else return the empty list [].

        The list is defined as the initial segment of OEIS A000364.

        The Euler numbers can be expressed with Bernoulli polynomials:

        E(n) = (2**(1 + 2*n)/(n + 1)) * (B(n + 1, 3/4) - B(n + 1, 1/4))

        Is equivalent to [abs(euler(n)) for n in range(0,2*size,2)].
        Is equivalent to andre_sequence(2*size)[::2].

        Examples
        ========

        >>> from sympy.functions.combinatorial.sequences import euler_sequence
        >>> euler_sequence(0)
        []

        >>> euler_sequence(9)
        [1, 1, 5, 61, 1385, 50521, 2702765, 199360981, 19391512145]

        References
        ==========

        .. [1] OEIS A000364, Euler secant numbers.
               https://oeis.org/A000364

        See Also
        ========

        - Andre numbers, OEIS A000111
        - Tangent numbers, OEIS A000182
    """

    @classmethod
    def eval(cls, size):
        s = sympify(size)
        if not s.is_Integer or s <= Integer(0):
            return []

        S = [0] * size
        S[0] = 1
        for k in range(1, size):
            S[k] = k * S[k - 1]
        for k in range(1, size - 1):
            for j in range(k, size):
                S[j] = (j - k) * S[j - 1] + (j - k + 1) * S[j]
        return S


class tangent_sequence(Function):
    """
        Input
        =====

        - size -- positive integer

        Output
        ======

        If 0 < size then return the list of nonzero tangent numbers

            [T1, T3, T5, ..., T{2*size-1}]

        else return the empty list [].

        This list is defined as the initial segment of OEIS A000182.

        Is equivalent to andre_sequence(2*size)[1::2].

        The tangent numbers can be expressed with the zeta function:

            T(n) = ((-16)**n - (-4)**n) * zeta(1 - 2 * n)

        Analytically the T(n) are the nonzero coefficients in the
        expansion of tan(x).

        Examples
        ========

        >>> from sympy.functions.combinatorial.sequences import tangent_sequence
        >>> tangent_sequence(0)
        []

        >>> tangent_sequence(9)
        [1, 2, 16, 272, 7936, 353792, 22368256, 1903757312, 209865342976]

         References
        ==========

        .. [1] OEIS A000182, tangent numbers.
               https://oeis.org/A000182

        See Also
        ========

        - Andre numbers, OEIS A000111
        - Euler numbers, OEIS A000364
    """

    @classmethod
    def eval(cls, size):
        s = sympify(size)
        if not s.is_Integer or s <= Integer(0):
            return []

        T = [0] * (size + 3)
        T[1] = 1
        for k in range(2, size + 3):
            T[k] = (k - 1) * T[k - 1]
        for k in range(2, size + 1):
            for j in range(k, size + 2):
                T[j] = (j - k) * T[j - 1] + (j - k + 2) * T[j]
        return T[1:-2]


class bell_sequence(Function):
    """
        Input
        =====

        - size -- positive integer

        Output
        ======

        If 0 < size then return the list of Bell numbers

            [B0, B1, B2, ..., B{n-1}]

        else return the empty list [].

        Is equivalent to [bell(n) for n in range(size)].

        Examples
        ========

        >>> from sympy.functions.combinatorial.sequences import bell_sequence
        >>> bell_sequence(0)
        []

        >>> bell_sequence(7)
        [1, 1, 2, 5, 15, 52, 203]

        References
        ==========

        .. [1] OEIS A000110: Bell or exponential numbers.
            https://oeis.org/A000110
    """

    @classmethod
    def eval(cls, size):
        s = sympify(size)
        if not s.is_Integer or s <= Integer(0):
            return []

        A = [0] * size
        A[0] = 1
        R = [1]
        for n in range(size - 1):
            A[n] = A[0]
            for k in range(n, 0, -1):
                A[k - 1] += A[k]
            R.append(A[0])
        return R
