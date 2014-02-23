from __future__ import print_function, division

from sympy.core import S, C, sympify
from sympy.core.function import Function, ArgumentIndexError
from sympy.ntheory import sieve
from math import sqrt as _sqrt

from sympy.core.compatibility import reduce, as_int, xrange
from sympy.core.cache import cacheit



class CombinatorialFunction(Function):
    """Base class for combinatorial functions. """

    def _eval_simplify(self, ratio, measure):
        from sympy.simplify.simplify import combsimp
        expr = combsimp(self)
        if measure(expr) <= ratio*measure(self):
            return expr
        return self

###############################################################################
######################## FACTORIAL and MULTI-FACTORIAL ########################
###############################################################################


class factorial(CombinatorialFunction):
    """Implementation of factorial function over nonnegative integers.
       By convention (consistent with the gamma function and the binomial
       coefficients), factorial of a negative integer is complex infinity.

       The factorial is very important in combinatorics where it gives
       the number of ways in which `n` objects can be permuted. It also
       arises in calculus, probability, number theory, etc.

       There is strict relation of factorial with gamma function. In
       fact n! = gamma(n+1) for nonnegative integers. Rewrite of this
       kind is very useful in case of combinatorial simplification.

       Computation of the factorial is done using two algorithms. For
       small arguments naive product is evaluated. However for bigger
       input algorithm Prime-Swing is used. It is the fastest algorithm
       known and computes n! via prime factorization of special class
       of numbers, called here the 'Swing Numbers'.

       Examples
       ========

       >>> from sympy import Symbol, factorial, S
       >>> n = Symbol('n', integer=True)

       >>> factorial(0)
       1

       >>> factorial(7)
       5040

       >>> factorial(-2)
       zoo

       >>> factorial(n)
       factorial(n)

       >>> factorial(2*n)
       factorial(2*n)

       >>> factorial(S(1)/2)
       factorial(1/2)

       See Also
       ========

       factorial2, RisingFactorial, FallingFactorial
    """

    def fdiff(self, argindex=1):
        if argindex == 1:
            return C.gamma(self.args[0] + 1)*C.polygamma(0, self.args[0] + 1)
        else:
            raise ArgumentIndexError(self, argindex)

    _small_swing = [
        1, 1, 1, 3, 3, 15, 5, 35, 35, 315, 63, 693, 231, 3003, 429, 6435, 6435, 109395,
        12155, 230945, 46189, 969969, 88179, 2028117, 676039, 16900975, 1300075,
        35102025, 5014575, 145422675, 9694845, 300540195, 300540195
    ]

    @classmethod
    def _swing(cls, n):
        if n < 33:
            return cls._small_swing[n]
        else:
            N, primes = int(_sqrt(n)), []

            for prime in sieve.primerange(3, N + 1):
                p, q = 1, n

                while True:
                    q //= prime

                    if q > 0:
                        if q & 1 == 1:
                            p *= prime
                    else:
                        break

                if p > 1:
                    primes.append(p)

            for prime in sieve.primerange(N + 1, n//3 + 1):
                if (n // prime) & 1 == 1:
                    primes.append(prime)

            L_product = R_product = 1

            for prime in sieve.primerange(n//2 + 1, n + 1):
                L_product *= prime

            for prime in primes:
                R_product *= prime

            return L_product*R_product

    @classmethod
    def _recursive(cls, n):
        if n < 2:
            return 1
        else:
            return (cls._recursive(n//2)**2)*cls._swing(n)

    @classmethod
    def eval(cls, n):
        n = sympify(n)

        if n.is_Number:
            if n is S.Zero:
                return S.One
            elif n is S.Infinity:
                return S.Infinity
            elif n.is_Integer:
                if n.is_negative:
                    return S.ComplexInfinity
                else:
                    n, result = n.p, 1

                    if n < 20:
                        for i in range(2, n + 1):
                            result *= i
                    else:
                        N, bits = n, 0

                        while N != 0:
                            if N & 1 == 1:
                                bits += 1

                            N = N >> 1

                        result = cls._recursive(n)*2**(n - bits)

                    return C.Integer(result)

    def _eval_rewrite_as_gamma(self, n):
        return C.gamma(n + 1)

    def _eval_is_integer(self):
        return self.args[0].is_integer

    def _eval_is_positive(self):
        if self.args[0].is_integer and self.args[0].is_positive:
            return True


class MultiFactorial(CombinatorialFunction):
    pass


class subfactorial(CombinatorialFunction):
    """The subfactorial counts the derangements of n items and is
    defined for non-negative integers as::

              ,
             |  1                             for n = 0
        !n = {  0                             for n = 1
             |  (n - 1)*(!(n - 1) + !(n - 2)) for n > 1
              `

    It can also be written as int(round(n!/exp(1))) but the recursive
    definition with caching is implemented for this function.

    References
    ==========
    .. [1] http://en.wikipedia.org/wiki/Subfactorial

    Examples
    ========

    >>> from sympy import subfactorial
    >>> from sympy.abc import n
    >>> subfactorial(n + 1)
    subfactorial(n + 1)
    >>> subfactorial(5)
    44

    See Also
    ========
    factorial, sympy.utilities.iterables.generate_derangements
    """

    @classmethod
    @cacheit
    def _eval(self, n):
        if not n:
            return 1
        elif n == 1:
            return 0
        return (n - 1)*(self._eval(n - 1) + self._eval(n - 2))

    @classmethod
    def eval(cls, arg):
        try:
            arg = as_int(arg)
            if arg < 0:
                raise ValueError
            return C.Integer(cls._eval(arg))
        except ValueError:
            if sympify(arg).is_Number:
                raise ValueError("argument must be a nonnegative integer")


class factorial2(CombinatorialFunction):
    """The double factorial n!!, not to be confused with (n!)!

    The double factorial is defined for integers >= -1 as::

               ,
              |  n*(n - 2)*(n - 4)* ... * 1    for n odd
        n!! = {  n*(n - 2)*(n - 4)* ... * 2    for n even
              |  1                             for n = 0, -1
               `

    Examples
    ========

    >>> from sympy import factorial2, var
    >>> var('n')
    n
    >>> factorial2(n + 1)
    factorial2(n + 1)
    >>> factorial2(5)
    15
    >>> factorial2(-1)
    1

    See Also
    ========

    factorial, RisingFactorial, FallingFactorial
    """

    @classmethod
    def eval(cls, arg):
        if arg.is_Number:
            if arg == S.Zero or arg == S.NegativeOne:
                return S.One
            return factorial2(arg - 2)*arg


###############################################################################
######################## RISING and FALLING FACTORIALS ########################
###############################################################################


class RisingFactorial(CombinatorialFunction):
    """Rising factorial (also called Pochhammer symbol) is a double valued
       function arising in concrete mathematics, hypergeometric functions
       and series expansions. It is defined by:

                   rf(x, k) = x * (x+1) * ... * (x + k-1)

       where 'x' can be arbitrary expression and 'k' is an integer. For
       more information check "Concrete mathematics" by Graham, pp. 66
       or visit http://mathworld.wolfram.com/RisingFactorial.html page.

       Examples
       ========

       >>> from sympy import rf
       >>> from sympy.abc import x

       >>> rf(x, 0)
       1

       >>> rf(1, 5)
       120

       >>> rf(x, 5) == x*(1 + x)*(2 + x)*(3 + x)*(4 + x)
       True

       See Also
       ========

       factorial, factorial2, FallingFactorial
    """

    @classmethod
    def eval(cls, x, k):
        x = sympify(x)
        k = sympify(k)

        if x is S.NaN:
            return S.NaN
        elif x is S.One:
            return factorial(k)
        elif k.is_Integer:
            if k is S.NaN:
                return S.NaN
            elif k is S.Zero:
                return S.One
            else:
                if k.is_positive:
                    if x is S.Infinity:
                        return S.Infinity
                    elif x is S.NegativeInfinity:
                        if k.is_odd:
                            return S.NegativeInfinity
                        else:
                            return S.Infinity
                    else:
                        return reduce(lambda r, i: r*(x + i), xrange(0, int(k)), 1)
                else:
                    if x is S.Infinity:
                        return S.Infinity
                    elif x is S.NegativeInfinity:
                        return S.Infinity
                    else:
                        return 1/reduce(lambda r, i: r*(x - i), xrange(1, abs(int(k)) + 1), 1)

    def _eval_rewrite_as_gamma(self, x, k):
        return C.gamma(x + k) / C.gamma(x)


class FallingFactorial(CombinatorialFunction):
    """Falling factorial (related to rising factorial) is a double valued
       function arising in concrete mathematics, hypergeometric functions
       and series expansions. It is defined by

                   ff(x, k) = x * (x-1) * ... * (x - k+1)

       where 'x' can be arbitrary expression and 'k' is an integer. For
       more information check "Concrete mathematics" by Graham, pp. 66
       or visit http://mathworld.wolfram.com/FallingFactorial.html page.

       >>> from sympy import ff
       >>> from sympy.abc import x

       >>> ff(x, 0)
       1

       >>> ff(5, 5)
       120

       >>> ff(x, 5) == x*(x-1)*(x-2)*(x-3)*(x-4)
       True

       See Also
       ========

       factorial, factorial2, RisingFactorial
    """

    @classmethod
    def eval(cls, x, k):
        x = sympify(x)
        k = sympify(k)

        if x is S.NaN:
            return S.NaN
        elif k.is_Integer:
            if k is S.NaN:
                return S.NaN
            elif k is S.Zero:
                return S.One
            else:
                if k.is_positive:
                    if x is S.Infinity:
                        return S.Infinity
                    elif x is S.NegativeInfinity:
                        if k.is_odd:
                            return S.NegativeInfinity
                        else:
                            return S.Infinity
                    else:
                        return reduce(lambda r, i: r*(x - i), xrange(0, int(k)), 1)
                else:
                    if x is S.Infinity:
                        return S.Infinity
                    elif x is S.NegativeInfinity:
                        return S.Infinity
                    else:
                        return 1/reduce(lambda r, i: r*(x + i), xrange(1, abs(int(k)) + 1), 1)

    def _eval_rewrite_as_gamma(self, x, k):
        return (-1)**k * C.gamma(-x + k) / C.gamma(-x)

rf = RisingFactorial
ff = FallingFactorial

###############################################################################
########################### BINOMIAL COEFFICIENTS #############################
###############################################################################


class binomial(CombinatorialFunction):
    """Implementation of the binomial coefficient. It can be defined
       in two ways depending on its desired interpretation:

           C(n,k) = n!/(k!(n-k)!)   or   C(n, k) = ff(n, k)/k!

       First, in a strict combinatorial sense it defines the
       number of ways we can choose 'k' elements from a set of
       'n' elements. In this case both arguments are nonnegative
       integers and binomial is computed using an efficient
       algorithm based on prime factorization.

       The other definition is generalization for arbitrary 'n',
       however 'k' must also be nonnegative. This case is very
       useful when evaluating summations.

       For the sake of convenience for negative 'k' this function
       will return zero no matter what valued is the other argument.

       To expand the binomial when n is a symbol, use either
       expand_func() or expand(func=True). The former will keep the
       polynomial in factored form while the latter will expand the
       polynomial itself. See examples for details.

       Examples
       ========

       >>> from sympy import Symbol, Rational, binomial, expand_func
       >>> n = Symbol('n', integer=True)

       >>> binomial(15, 8)
       6435

       >>> binomial(n, -1)
       0

       >>> [ binomial(0, i) for i in range(1)]
       [1]
       >>> [ binomial(1, i) for i in range(2)]
       [1, 1]
       >>> [ binomial(2, i) for i in range(3)]
       [1, 2, 1]
       >>> [ binomial(3, i) for i in range(4)]
       [1, 3, 3, 1]
       >>> [ binomial(4, i) for i in range(5)]
       [1, 4, 6, 4, 1]

       >>> binomial(Rational(5,4), 3)
       -5/128

       >>> binomial(n, 3)
       binomial(n, 3)

       >>> binomial(n, 3).expand(func=True)
       n**3/6 - n**2/2 + n/3

       >>> expand_func(binomial(n, 3))
       n*(n - 2)*(n - 1)/6

    """

    def fdiff(self, argindex=1):
        if argindex == 1:
            # http://functions.wolfram.com/GammaBetaErf/Binomial/20/01/01/
            n, k = self.args
            return binomial(n, k)*(C.polygamma(0, n + 1) - C.polygamma(0, n - k + 1))
        elif argindex == 2:
            # http://functions.wolfram.com/GammaBetaErf/Binomial/20/01/02/
            n, k = self.args
            return binomial(n, k)*(C.polygamma(0, n - k + 1) - C.polygamma(0, k + 1))
        else:
            raise ArgumentIndexError(self, argindex)

    @classmethod
    def eval(cls, n, k):
        n, k = map(sympify, (n, k))

        if k.is_Number:
            if k.is_Integer:
                if k < 0:
                    return S.Zero
                elif k == 0 or n == k:
                    return S.One
                elif n.is_Integer and n >= 0:
                    n, k = int(n), int(k)

                    if k > n:
                        return S.Zero
                    elif k > n // 2:
                        k = n - k

                    M, result = int(_sqrt(n)), 1

                    for prime in sieve.primerange(2, n + 1):
                        if prime > n - k:
                            result *= prime
                        elif prime > n // 2:
                            continue
                        elif prime > M:
                            if n % prime < k % prime:
                                result *= prime
                        else:
                            N, K = n, k
                            exp = a = 0

                            while N > 0:
                                a = int((N % prime) < (K % prime + a))
                                N, K = N // prime, K // prime
                                exp = a + exp

                            if exp > 0:
                                result *= prime**exp

                    return C.Integer(result)
                elif n.is_Number:
                    result = n - k + 1
                    for i in xrange(2, k + 1):
                        result *= n - k + i
                        result /= i
                    return result

        elif k.is_negative:
            return S.Zero
        elif (n - k).simplify().is_negative:
            return S.Zero
        else:
            d = n - k

            if d.is_Integer:
                return cls.eval(n, d)

    def _eval_expand_func(self, **hints):
        """
        Function to expand binomial(n,k) when m is positive integer
        Also,
        n is self.args[0] and k is self.args[1] while using binomial(n, k)
        """
        n = self.args[0]
        if n.is_Number:
            return binomial(*self.args)

        k = self.args[1]
        if k.is_Add and n in k.args:
            k = n - k

        if k.is_Integer:
            if k == S.Zero:
                return S.One
            elif k < 0:
                return S.Zero
            else:
                n = self.args[0]
                result = n - k + 1
                for i in xrange(2, k + 1):
                    result *= n - k + i
                    result /= i
                return result
        else:
            return binomial(*self.args)

    def _eval_rewrite_as_factorial(self, n, k):
        return C.factorial(n)/(C.factorial(k)*C.factorial(n - k))

    def _eval_rewrite_as_gamma(self, n, k):
        return C.gamma(n + 1)/(C.gamma(k + 1)*C.gamma(n - k + 1))

    def _eval_is_integer(self):
        return self.args[0].is_integer and self.args[1].is_integer
