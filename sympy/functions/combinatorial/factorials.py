from __future__ import print_function, division

from sympy.core import S, sympify, Dummy
from sympy.core.function import Function, ArgumentIndexError
from sympy.core.logic import fuzzy_and
from sympy.core.numbers import Integer
from sympy.ntheory import sieve
from math import sqrt as _sqrt

from sympy.core.compatibility import reduce, range
from sympy.core.cache import cacheit

from sympy.polys.polytools import poly_from_expr
from sympy.polys.polyerrors import PolificationFailed


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
        from sympy import gamma, polygamma
        if argindex == 1:
            return gamma(self.args[0] + 1)*polygamma(0, self.args[0] + 1)
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

                    return Integer(result)

    def _eval_rewrite_as_gamma(self, n):
        from sympy import gamma
        return gamma(n + 1)

    def _eval_rewrite_as_Product(self, n):
        from sympy import Product
        if n.is_nonnegative and n.is_integer:
            i = Dummy('i', integer=True)
            return Product(i, (i, 1, n))

    def _eval_is_integer(self):
        if self.args[0].is_integer and self.args[0].is_nonnegative:
            return True

    def _eval_is_positive(self):
        if self.args[0].is_integer and self.args[0].is_nonnegative:
            return True

    def _eval_is_composite(self):
        x = self.args[0]
        if x.is_integer:
            return (x - 3).is_nonnegative

    def _eval_is_real(self):
        x = self.args[0]
        if x.is_nonnegative or x.is_noninteger:
            return True


class MultiFactorial(CombinatorialFunction):
    pass


class subfactorial(CombinatorialFunction):
    r"""The subfactorial counts the derangements of n items and is
    defined for non-negative integers as::

              ,
             |  1                             for n = 0
        !n = {  0                             for n = 1
             |  (n - 1)*(!(n - 1) + !(n - 2)) for n > 1
              `

    It can also be written as int(round(n!/exp(1))) but the recursive
    definition with caching is implemented for this function.

    An interesting analytic expression is the following [2]_

    .. math:: !x = \Gamma(x + 1, -1)/e

    which is valid for non-negative integers x. The above formula
    is not very useful incase of non-integers. :math:`\Gamma(x + 1, -1)` is
    single-valued only for integral arguments x, elsewhere on the positive real
    axis it has an infinite number of branches none of which are real.

    References
    ==========

    .. [1] http://en.wikipedia.org/wiki/Subfactorial
    .. [2] http://mathworld.wolfram.com/Subfactorial.html

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

    sympy.functions.combinatorial.factorials.factorial,
    sympy.utilities.iterables.generate_derangements,
    sympy.functions.special.gamma_functions.uppergamma
    """

    @classmethod
    @cacheit
    def _eval(self, n):
        if not n:
            return S.One
        elif n == 1:
            return S.Zero
        return (n - 1)*(self._eval(n - 1) + self._eval(n - 2))

    @classmethod
    def eval(cls, arg):
        if arg.is_Number:
            if arg.is_Integer and arg.is_nonnegative:
                return cls._eval(arg)
            elif arg is S.NaN:
                return S.NaN
            elif arg is S.Infinity:
                return S.Infinity

    def _eval_is_even(self):
        if self.args[0].is_odd and self.args[0].is_nonnegative:
            return True

    def _eval_is_integer(self):
        if self.args[0].is_integer and self.args[0].is_nonnegative:
            return True

    def _eval_rewrite_as_uppergamma(self, arg):
        from sympy import uppergamma
        return uppergamma(arg + 1, -1)/S.Exp1

    def _eval_is_nonnegative(self):
        if self.args[0].is_integer and self.args[0].is_nonnegative:
            return True

    def _eval_is_odd(self):
        if self.args[0].is_even and self.args[0].is_nonnegative:
            return True


class factorial2(CombinatorialFunction):
    """The double factorial n!!, not to be confused with (n!)!

    The double factorial is defined for nonnegative integers and for odd
    negative integers as::

               ,
              |  n*(n - 2)*(n - 4)* ... * 1    for n positive odd
        n!! = {  n*(n - 2)*(n - 4)* ... * 2    for n positive even
              |  1                             for n = 0
              |  (n+2)!! / (n+2)               for n negative odd
               `

    References
    ==========
    .. [1] https://en.wikipedia.org/wiki/Double_factorial

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
    >>> factorial2(-5)
    1/3

    See Also
    ========

    factorial, RisingFactorial, FallingFactorial
    """

    @classmethod
    def eval(cls, arg):
        # TODO: extend this to complex numbers?
        if arg.is_integer:
            if arg.is_infinite:
                return

            # This implementation is faster than the recursive one
            # It also avoids "maximum recursion depth exceeded" runtime error
            if arg.is_nonnegative:
                if arg.is_even:
                    k = arg / 2
                    return 2 ** k * factorial(k)

                return factorial(arg) / factorial2(arg - 1)

            if arg.is_even:
                raise ValueError("argument must be nonnegative or odd")

            return arg * (S.NegativeOne) ** ((1 - arg) / 2) / factorial2(-arg)

    def _eval_is_even(self):
        # Double factorial is even for every positive even input
        n = self.args[0]
        if n.is_integer:
            if n.is_odd:
                return False
            if n.is_even:
                if n.is_positive:
                    return True
                if n.is_zero:
                    return False

    def _eval_is_integer(self):
        # Double factorial is an integer for every nonnegative input, and for
        # -1 and -3
        n = self.args[0]
        if n.is_integer:
            if (n + 1).is_nonnegative:
                return True
            if n.is_odd:
                return (n + 3).is_nonnegative

    def _eval_is_odd(self):
        # Double factorial is odd for every odd input not smaller than -3, and
        # for 0
        n = self.args[0]
        if n.is_odd:
            return (n + 3).is_nonnegative
        if n.is_even:
            if n.is_positive:
                return False
            if n.is_zero:
                return True

    def _eval_is_positive(self):
        # Double factorial is positive for every nonnegative input, and for
        # every odd negative input which is of the form -1-4k for an
        # nonnegative integer k
        n = self.args[0]
        if n.is_integer:
            if (n + 1).is_nonnegative:
                return True
            if n.is_odd:
                return ((n + 1) / 2).is_even


###############################################################################
######################## RISING and FALLING FACTORIALS ########################
###############################################################################


class RisingFactorial(CombinatorialFunction):
    """Rising factorial (also called Pochhammer symbol) is a double valued
    function arising in concrete mathematics, hypergeometric functions
    and series expansions. It is defined by:

                rf(x, k) = x * (x + 1) * ... * (x + k - 1)

    where 'x' can be arbitrary expression and 'k' is an integer. For
    more information check "Concrete mathematics" by Graham, pp. 66
    or visit http://mathworld.wolfram.com/RisingFactorial.html page.

    When x is a polynomial f of a single variable y of order >= 1,
    rf(x,k) = f(y) * f(y+1) * ... * f(x+k-1) as described in
    Peter Paule, "Greatest Factorial Factorization and Symbolic Summation",
    Journal of Symbolic Computation, vol. 20, pp. 235-268, 1995.

    Examples
    ========

    >>> from sympy import rf, symbols, factorial, ff, binomial
    >>> from sympy.abc import x
    >>> n, k = symbols('n k', integer=True)
    >>> rf(x, 0)
    1
    >>> rf(1, 5)
    120
    >>> rf(x, 5) == x*(1 + x)*(2 + x)*(3 + x)*(4 + x)
    True
    >>> rf(x**3, 2)
    Poly(x**6 + 3*x**5 + 3*x**4 + x**3, x, domain='ZZ')

    Rewrite

    >>> rf(x, k).rewrite(ff)
    FallingFactorial(k + x - 1, k)
    >>> rf(x, k).rewrite(binomial)
    binomial(k + x - 1, k)*factorial(k)
    >>> rf(n, k).rewrite(factorial)
    factorial(k + n - 1)/factorial(n - 1)

    See Also
    ========

    factorial, factorial2, FallingFactorial

    References
    ==========

    .. [1] https://en.wikipedia.org/wiki/Pochhammer_symbol
    """

    @classmethod
    def eval(cls, x, k):
        x = sympify(x)
        k = sympify(k)

        if x is S.NaN or k is S.NaN:
            return S.NaN
        elif x is S.One:
            return factorial(k)
        elif k.is_Integer:
            if k is S.Zero:
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
                        try:
                            F, opt = poly_from_expr(x)
                        except PolificationFailed:
                            return reduce(lambda r, i: r*(x + i), range(0, int(k)), 1)
                        if len(opt.gens) > 1 or F.degree() <= 1:
                            return reduce(lambda r, i: r*(x + i), range(0, int(k)), 1)
                        else:
                            v = opt.gens[0]
                            return reduce(lambda r, i:
                                          r*(F.subs(v, v + i).expand()),
                                          range(0, int(k)), 1)
                else:
                    if x is S.Infinity:
                        return S.Infinity
                    elif x is S.NegativeInfinity:
                        return S.Infinity
                    else:
                        try:
                            F, opt = poly_from_expr(x)
                        except PolificationFailed:
                            return 1/reduce(lambda r, i:
                                            r*(x - i),
                                            range(1, abs(int(k)) + 1), 1)
                        if len(opt.gens) > 1 or F.degree() <= 1:
                            return 1/reduce(lambda r, i:
                                            r*(x - i),
                                            range(1, abs(int(k)) + 1), 1)
                        else:
                            v = opt.gens[0]
                            return 1/reduce(lambda r, i:
                                            r*(F.subs(v, v - i).expand()),
                                            range(1, abs(int(k)) + 1), 1)

    def _eval_rewrite_as_gamma(self, x, k):
        from sympy import gamma
        return gamma(x + k) / gamma(x)

    def _eval_rewrite_as_FallingFactorial(self, x, k):
        return FallingFactorial(x + k - 1, k)

    def _eval_rewrite_as_factorial(self, x, k):
        if x.is_integer and k.is_integer:
            return factorial(k + x - 1) / factorial(x - 1)

    def _eval_rewrite_as_binomial(self, x, k):
        if k.is_integer:
            return factorial(k) * binomial(x + k - 1, k)

    def _eval_is_integer(self):
        return fuzzy_and((self.args[0].is_integer, self.args[1].is_integer,
                          self.args[1].is_nonnegative))

    def _sage_(self):
        import sage.all as sage
        return sage.rising_factorial(self.args[0]._sage_(), self.args[1]._sage_())


class FallingFactorial(CombinatorialFunction):
    """Falling factorial (related to rising factorial) is a double valued
    function arising in concrete mathematics, hypergeometric functions
    and series expansions. It is defined by

                ff(x, k) = x * (x-1) * ... * (x - k+1)

    where 'x' can be arbitrary expression and 'k' is an integer. For
    more information check "Concrete mathematics" by Graham, pp. 66
    or visit http://mathworld.wolfram.com/FallingFactorial.html page.

    When x is a polynomial f of a single variable y of order >= 1,
    ff(x,k) = f(y) * f(y-1) * ... * f(x-k+1) as described in
    Peter Paule, "Greatest Factorial Factorization and Symbolic Summation",
    Journal of Symbolic Computation, vol. 20, pp. 235-268, 1995.

    >>> from sympy import ff, factorial, rf, gamma, polygamma, binomial, symbols
    >>> from sympy.abc import x, k
    >>> n, m = symbols('n m', integer=True)
    >>> ff(x, 0)
    1
    >>> ff(5, 5)
    120
    >>> ff(x, 5) == x*(x-1)*(x-2)*(x-3)*(x-4)
    True
    >>> ff(x**2, 2)
    Poly(x**4 - 2*x**3 + x**2, x, domain='ZZ')
    >>> ff(n, n)
    factorial(n)

    Rewrite

    >>> ff(x, k).rewrite(gamma)
    (-1)**k*gamma(k - x)/gamma(-x)
    >>> ff(x, k).rewrite(rf)
    RisingFactorial(-k + x + 1, k)
    >>> ff(x, m).rewrite(binomial)
    binomial(x, m)*factorial(m)
    >>> ff(n, m).rewrite(factorial)
    factorial(n)/factorial(-m + n)

    See Also
    ========

    factorial, factorial2, RisingFactorial

    References
    ==========

    .. [1] http://mathworld.wolfram.com/FallingFactorial.html
    """

    @classmethod
    def eval(cls, x, k):
        x = sympify(x)
        k = sympify(k)

        if x is S.NaN or k is S.NaN:
            return S.NaN
        elif k.is_integer and x == k:
            return factorial(x)
        elif k.is_Integer:
            if k is S.Zero:
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
                        try:
                            F, opt = poly_from_expr(x)
                        except PolificationFailed:
                            return reduce(lambda r, i: r*(x - i),
                                          range(0, int(k)), 1)
                        if len(opt.gens) > 1 or F.degree() <= 1:
                            return reduce(lambda r, i: r*(x - i),
                                          range(0, int(k)), 1)
                        else:
                            v = opt.gens[0]
                            return reduce(lambda r, i:
                                          r*(F.subs(v, v - i).expand()),
                                          range(0, int(k)), 1)
                else:
                    if x is S.Infinity:
                        return S.Infinity
                    elif x is S.NegativeInfinity:
                        return S.Infinity
                    else:
                        try:
                            F, opt = poly_from_expr(x)
                        except PolificationFailed:
                            return 1/reduce(lambda r, i: r*(x + i),
                                            range(1, abs(int(k)) + 1), 1)
                        if len(opt.gens) > 1 or F.degree() <= 1:
                            return 1/reduce(lambda r, i: r*(x + i),
                                            range(1, abs(int(k)) + 1), 1)
                        else:
                            v = opt.gens[0]
                            return 1/reduce(lambda r, i:
                                            r*(F.subs(v, v + i).expand()),
                                            range(1, abs(int(k)) + 1), 1)

    def _eval_rewrite_as_gamma(self, x, k):
        from sympy import gamma
        return (-1)**k*gamma(k - x) / gamma(-x)

    def _eval_rewrite_as_RisingFactorial(self, x, k):
        return rf(x - k + 1, k)

    def _eval_rewrite_as_binomial(self, x, k):
        if k.is_integer:
            return factorial(k) * binomial(x, k)

    def _eval_rewrite_as_factorial(self, x, k):
        if x.is_integer and k.is_integer:
            return factorial(x) / factorial(x - k)

    def _eval_is_integer(self):
        return fuzzy_and((self.args[0].is_integer, self.args[1].is_integer,
                          self.args[1].is_nonnegative))

    def _sage_(self):
        import sage.all as sage
        return sage.falling_factorial(self.args[0]._sage_(),
                                      self.args[1]._sage_())


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
    >>> n = Symbol('n', integer=True, positive=True)

    >>> binomial(15, 8)
    6435

    >>> binomial(n, -1)
    0

    Rows of Pascal's triangle can be generated with the binomial function:

    >>> for N in range(8):
    ...     print([ binomial(N, i) for i in range(N + 1)])
    ...
    [1]
    [1, 1]
    [1, 2, 1]
    [1, 3, 3, 1]
    [1, 4, 6, 4, 1]
    [1, 5, 10, 10, 5, 1]
    [1, 6, 15, 20, 15, 6, 1]
    [1, 7, 21, 35, 35, 21, 7, 1]

    As can a given diagonal, e.g. the 4th diagonal:

    >>> N = -4
    >>> [ binomial(N, i) for i in range(1 - N)]
    [1, -4, 10, -20, 35]

    >>> binomial(Rational(5, 4), 3)
    -5/128
    >>> binomial(Rational(-5, 4), 3)
    -195/128

    >>> binomial(n, 3)
    binomial(n, 3)

    >>> binomial(n, 3).expand(func=True)
    n**3/6 - n**2/2 + n/3

    >>> expand_func(binomial(n, 3))
    n*(n - 2)*(n - 1)/6

    """

    def fdiff(self, argindex=1):
        from sympy import polygamma
        if argindex == 1:
            # http://functions.wolfram.com/GammaBetaErf/Binomial/20/01/01/
            n, k = self.args
            return binomial(n, k)*(polygamma(0, n + 1) - \
                polygamma(0, n - k + 1))
        elif argindex == 2:
            # http://functions.wolfram.com/GammaBetaErf/Binomial/20/01/02/
            n, k = self.args
            return binomial(n, k)*(polygamma(0, n - k + 1) - \
                polygamma(0, k + 1))
        else:
            raise ArgumentIndexError(self, argindex)

    @classmethod
    def _eval(self, n, k):
        # n.is_Number and k.is_Integer and k != 1 and n != k
        if k.is_Integer:
            if n.is_Integer and n >= 0:
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
                return Integer(result)
            else:
                d = result = n - k + 1
                for i in range(2, k + 1):
                    d += 1
                    result *= d
                    result /= i
                return result

    @classmethod
    def eval(cls, n, k):
        n, k = map(sympify, (n, k))
        d = n - k
        if d.is_zero or k.is_zero:
            return S.One
        elif d.is_zero is False:
            if (k - 1).is_zero:
                return n
            elif k.is_negative:
                return S.Zero
            elif n.is_integer and n.is_nonnegative and d.is_negative:
                return S.Zero
        if k.is_Integer and k > 0 and n.is_Number:
            return cls._eval(n, k)

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
                for i in range(2, k + 1):
                    result *= n - k + i
                    result /= i
                return result
        else:
            return binomial(*self.args)

    def _eval_rewrite_as_factorial(self, n, k):
        return factorial(n)/(factorial(k)*factorial(n - k))

    def _eval_rewrite_as_gamma(self, n, k):
        from sympy import gamma
        return gamma(n + 1)/(gamma(k + 1)*gamma(n - k + 1))

    def _eval_rewrite_as_FallingFactorial(self, n, k):
        if k.is_integer:
            return ff(n, k) / factorial(k)

    def _eval_is_integer(self):
        n, k = self.args
        if n.is_integer and k.is_integer:
            return True
        elif k.is_integer is False:
            return False
