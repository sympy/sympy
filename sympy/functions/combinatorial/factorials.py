from typing import List

from sympy.core import S, sympify, Mod
from sympy.core.cache import cacheit
from sympy.core.compatibility import reduce, HAS_GMPY, as_int
from sympy.core.function import Function, ArgumentIndexError, _mexpand
from sympy.core.logic import fuzzy_and
from sympy.core.mul import Mul
from sympy.core.numbers import Float, Integer, pi
from sympy.core.relational import Eq
from sympy.core.symbol import Symbol, Dummy
from sympy.functions.elementary.piecewise import Piecewise
from sympy.functions.elementary.integers import floor
from sympy.logic.boolalg import Or, And
from sympy.ntheory import sieve
from sympy.polys.polytools import Poly
from sympy.utilities.iterables import numbered_symbols

from math import sqrt as _sqrt


def int_like(x):
    if x.is_integer is False:
        if x.is_Float and int(x) == x:
            return True
        return False
    return x.is_integer


def igen(i='i'):
    """return unnumbered symbol and then numbered versions

    Examples
    ========

    >>> from sympy.functions.combinatorial.factorials import igen
    >>> [i for j,i in zip((1, 2, 3), igen())]
    [i, i_0, i_1]
    >>> [i for j,i in zip((1, 2, 3), igen('x1'))]
    [x1, x1_0, x1_1]
    """
    yield Symbol(i)
    i += '_'
    for ii in numbered_symbols(i):
        yield ii

class CombinatorialFunction(Function):
    """Base class for combinatorial functions. """

    def _eval_simplify(self, **kwargs):
        from sympy.simplify.combsimp import combsimp
        # combinatorial function with non-integer arguments is
        # automatically passed to gammasimp
        expr = combsimp(self)
        measure = kwargs['measure']
        if measure(expr) <= kwargs['ratio']*measure(self):
            return expr
        return self


###############################################################################
######################## FACTORIAL and MULTI-FACTORIAL ########################
###############################################################################


class factorial(CombinatorialFunction):
    r"""Implementation of factorial function over nonnegative integers.
       By convention (consistent with the gamma function and the binomial
       coefficients), factorial of a negative integer is complex infinity.

       The factorial is very important in combinatorics where it gives
       the number of ways in which `n` objects can be permuted. It also
       arises in calculus, probability, number theory, etc.

       There is strict relation of factorial with gamma function. In
       fact `n! = gamma(n+1)` for nonnegative integers. Rewrite of this
       kind is very useful in case of combinatorial simplification.

       Computation of the factorial is done using two algorithms. For
       small arguments a precomputed look up table is used. However for bigger
       input algorithm Prime-Swing is used. It is the fastest algorithm
       known and computes `n!` via prime factorization of special class
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

    _small_factorials = []  # type: List[int]

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
            if n.is_zero:
                return S.One
            elif n is S.Infinity:
                return S.Infinity
            elif n.is_Integer:
                if n.is_negative:
                    return S.ComplexInfinity
                else:
                    n = n.p

                    if n < 20:
                        if not cls._small_factorials:
                            result = 1
                            for i in range(1, 20):
                                result *= i
                                cls._small_factorials.append(result)
                        result = cls._small_factorials[n-1]

                    # GMPY factorial is faster, use it when available
                    elif HAS_GMPY:
                        from sympy.core.compatibility import gmpy
                        result = gmpy.fac(n)

                    else:
                        bits = bin(n).count('1')
                        result = cls._recursive(n)*2**(n - bits)

                    return Integer(result)

    def _facmod(self, n, q):
        res, N = 1, int(_sqrt(n))

        # Exponent of prime p in n! is e_p(n) = [n/p] + [n/p**2] + ...
        # for p > sqrt(n), e_p(n) < sqrt(n), the primes with [n/p] = m,
        # occur consecutively and are grouped together in pw[m] for
        # simultaneous exponentiation at a later stage
        pw = [1]*N

        m = 2 # to initialize the if condition below
        for prime in sieve.primerange(2, n + 1):
            if m > 1:
                m, y = 0, n // prime
                while y:
                    m += y
                    y //= prime
            if m < N:
                pw[m] = pw[m]*prime % q
            else:
                res = res*pow(prime, m, q) % q

        for ex, bs in enumerate(pw):
            if ex == 0 or bs == 1:
                continue
            if bs == 0:
                return 0
            res = res*pow(bs, ex, q) % q

        return res

    def _eval_Mod(self, q):
        n = self.args[0]
        if n.is_integer and n.is_nonnegative and q.is_integer:
            aq = abs(q)
            d = aq - n
            if d.is_nonpositive:
                return S.Zero
            else:
                isprime = aq.is_prime
                if d == 1:
                    # Apply Wilson's theorem (if a natural number n > 1
                    # is a prime number, then (n-1)! = -1 mod n) and
                    # its inverse (if n > 4 is a composite number, then
                    # (n-1)! = 0 mod n)
                    if isprime:
                        return S(-1 % q)
                    elif isprime is False and (aq - 6).is_nonnegative:
                        return S.Zero
                elif n.is_Integer and q.is_Integer:
                    n, d, aq = map(int, (n, d, aq))
                    if isprime and (d - 1 < n):
                        fc = self._facmod(d - 1, aq)
                        fc = pow(fc, aq - 2, aq)
                        if d % 2:
                            fc = -fc
                    else:
                        fc = self._facmod(n, aq)

                    return S(fc % q)

    def _eval_rewrite_as_gamma(self, n, **kwargs):
        from sympy import gamma
        return gamma(n + 1)

    def _eval_rewrite_as_Product(self, n, **kwargs):
        from sympy.concrete.products import Product
        if n.is_nonnegative and int_like(n):
            for i in igen():
                if not n.has(i):
                    break
            return Product(i, (i, 1, n))

    def _eval_is_integer(self):
        x = self.args[0]
        if x.is_nonnegative:
            if x.is_integer:
                return True
            if int_like(x):
                # this handles is_even and is_composite, too
                # because neither of those is T if self is
                # not an integer
                return False

    def _eval_is_positive(self):
        if int_like(self.args[0]) and self.args[0].is_nonnegative:
            return True

    def _eval_is_even(self):
        x = self.args[0]
        if x.is_nonnegative:
            if x.is_integer:
                return (x - 2).is_nonnegative

    def _eval_is_composite(self):
        x = self.args[0]
        if x.is_nonnegative:
            if x.is_integer:
                return (x - 3).is_nonnegative

    def _eval_is_real(self):
        x = self.args[0]
        if x.is_nonnegative or x.is_real and int_like(x) is False:
            return True

    def _eval_as_leading_term(self, x):
        from sympy import Order
        arg = self.args[0]
        arg_1 = arg.as_leading_term(x)
        if Order(x, x).contains(arg_1):
            return S.One
        if Order(1, x).contains(arg_1):
            return self.func(arg_1)
        ####################################################
        # The correct result here should be 'None'.        #
        # Indeed arg is not bounded as x tends to 0.       #
        # Consequently the series expansion does not admit #
        # the leading term.                                #
        # For compatibility reasons, the return value here #
        # is the original function, i.e. factorial(arg),   #
        # instead of None.                                 #
        ####################################################
        return self.func(arg)

class MultiFactorial(CombinatorialFunction):
    pass


class subfactorial(CombinatorialFunction):
    r"""The subfactorial counts the derangements of n items and is
    defined for non-negative integers as:

    .. math:: !n = \begin{cases} 1 & n = 0 \\ 0 & n = 1 \\
                    (n-1)(!(n-1) + !(n-2)) & n > 1 \end{cases}

    It can also be written as ``int(round(n!/exp(1)))`` but the
    recursive definition with caching is implemented for this function.

    An interesting analytic expression is the following [2]_

    .. math:: !x = \Gamma(x + 1, -1)/e

    which is valid for non-negative integers `x`. The above formula
    is not very useful incase of non-integers. :math:`\Gamma(x + 1, -1)` is
    single-valued only for integral arguments `x`, elsewhere on the positive
    real axis it has an infinite number of branches none of which are real.

    References
    ==========

    .. [1] https://en.wikipedia.org/wiki/Subfactorial
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
    def _eval(cls, n):
        if not n:
            return S.One
        elif n == 1:
            return S.Zero
        else:
            z1, z2 = 1, 0
            for i in range(2, n + 1):
                z1, z2 = z2, (i - 1)*(z2 + z1)
            return z2

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

    def _eval_rewrite_as_uppergamma(self, arg, **kwargs):
        from sympy import uppergamma
        return uppergamma(arg + 1, -1)/S.Exp1

    def _eval_is_nonnegative(self):
        if self.args[0].is_integer and self.args[0].is_nonnegative:
            return True

    def _eval_is_odd(self):
        if self.args[0].is_even and self.args[0].is_nonnegative:
            return True


class factorial2(CombinatorialFunction):
    r"""The double factorial `n!!`, not to be confused with `(n!)!`

    The double factorial is defined for nonnegative integers and for odd
    negative integers as:

    .. math:: n!! = \begin{cases} 1 & n = 0 \\
                    n(n-2)(n-4) \cdots 1 & n\ \text{positive odd} \\
                    n(n-2)(n-4) \cdots 2 & n\ \text{positive even} \\
                    (n+2)!!/(n+2) & n\ \text{negative odd} \end{cases}

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

        if arg.is_Number:
            if not arg.is_Integer:
                raise ValueError("argument must be nonnegative integer "
                                    "or negative odd integer")

            # This implementation is faster than the recursive one
            # It also avoids "maximum recursion depth exceeded" runtime error
            if arg.is_nonnegative:
                if arg.is_even:
                    k = arg / 2
                    return 2**k * factorial(k)
                return factorial(arg) / factorial2(arg - 1)

            if arg.is_odd:
                return arg*(S.NegativeOne)**((1 - arg)/2) / factorial2(-arg)
            raise ValueError("argument must be nonnegative integer "
                                "or negative odd integer")

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
                # -1, 0, 1, ...
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

    def _eval_rewrite_as_gamma(self, n, **kwargs):
        from sympy import gamma, Piecewise, sqrt
        return 2**(n/2)*gamma(n/2 + 1) * Piecewise((1, Eq(Mod(n, 2), 0)),
                (sqrt(2/pi), Eq(Mod(n, 2), 1)))


###############################################################################
######################## RISING and FALLING FACTORIALS ########################
###############################################################################


class RisingFactorial(CombinatorialFunction):
    r"""
    Rising factorial (also called Pochhammer symbol) is a double valued
    function arising in concrete mathematics, hypergeometric functions
    and series expansions. It is defined by:

    .. math:: rf(x,k) = x \cdot (x+1) \cdots (x+k-1)

    where `x` can be arbitrary expression and `k` is an integer. For
    more information check "Concrete mathematics" by Graham, pp. 66
    or visit http://mathworld.wolfram.com/RisingFactorial.html page.

    When `x` is a Poly instance of degree >= 1 with a single variable,
    `rf(x,k) = x(y) \cdot x(y+1) \cdots x(y+k-1)`, where `y` is the
    variable of `x`. This is as described in Peter Paule, "Greatest
    Factorial Factorization and Symbolic Summation", Journal of
    Symbolic Computation, vol. 20, pp. 235-268, 1995.

    Examples
    ========

    >>> from sympy import rf, symbols, factorial, ff, binomial, Poly
    >>> from sympy.abc import x
    >>> n, k = symbols('n k', integer=True)
    >>> rf(x, 0)
    1
    >>> rf(1, 5)
    120
    >>> rf(x, 5) == x*(1 + x)*(2 + x)*(3 + x)*(4 + x)
    True
    >>> rf(Poly(x**3, x), 2)
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
            if k.is_zero:
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
                        if isinstance(x, Poly):
                            gens = x.gens
                            if len(gens)!= 1:
                                raise ValueError("rf only defined for "
                                            "polynomials on one generator")
                            else:
                                return reduce(lambda r, i:
                                              r*(x.shift(i)),
                                              range(0, int(k)), 1)
                        else:
                            return reduce(lambda r, i: r*(x + i),
                                        range(0, int(k)), 1)

                else:
                    if x is S.Infinity:
                        return S.Infinity
                    elif x is S.NegativeInfinity:
                        return S.Infinity
                    else:
                        if isinstance(x, Poly):
                            gens = x.gens
                            if len(gens)!= 1:
                                raise ValueError("rf only defined for "
                                            "polynomials on one generator")
                            else:
                                return 1/reduce(lambda r, i:
                                                r*(x.shift(-i)),
                                                range(1, abs(int(k)) + 1), 1)
                        else:
                            return 1/reduce(lambda r, i:
                                            r*(x - i),
                                            range(1, abs(int(k)) + 1), 1)

        if k.is_integer == False:
            if x.is_integer and x.is_negative:
                return S.Zero

    def _eval_rewrite_as_gamma(self, x, k, **kwargs):
        from sympy import gamma
        return gamma(x + k) / gamma(x)

    def _eval_rewrite_as_FallingFactorial(self, x, k, **kwargs):
        return FallingFactorial(x + k - 1, k)

    def _eval_rewrite_as_factorial(self, x, k, **kwargs):
        if x.is_integer and k.is_integer:
            return factorial(k + x - 1) / factorial(x - 1)

    def _eval_rewrite_as_binomial(self, x, k, **kwargs):
        if k.is_integer:
            return factorial(k) * binomial(x + k - 1, k)

    def _eval_is_integer(self):
        return fuzzy_and((self.args[0].is_integer, self.args[1].is_integer,
                          self.args[1].is_nonnegative))

    def _sage_(self):
        import sage.all as sage
        return sage.rising_factorial(self.args[0]._sage_(),
                                     self.args[1]._sage_())


class FallingFactorial(CombinatorialFunction):
    r"""
    Falling factorial (related to rising factorial) is a double valued
    function arising in concrete mathematics, hypergeometric functions
    and series expansions. It is defined by

    .. math:: ff(x,k) = x \cdot (x-1) \cdots (x-k+1)

    where `x` can be arbitrary expression and `k` is an integer. For
    more information check "Concrete mathematics" by Graham, pp. 66
    or visit http://mathworld.wolfram.com/FallingFactorial.html page.

    When `x` is a Poly instance of degree >= 1 with single variable,
    `ff(x,k) = x(y) \cdot x(y-1) \cdots x(y-k+1)`, where `y` is the
    variable of `x`. This is as described in Peter Paule, "Greatest
    Factorial Factorization and Symbolic Summation", Journal of
    Symbolic Computation, vol. 20, pp. 235-268, 1995.

    >>> from sympy import ff, factorial, rf, gamma, polygamma, binomial, symbols, Poly
    >>> from sympy.abc import x, k
    >>> n, m = symbols('n m', integer=True)
    >>> ff(x, 0)
    1
    >>> ff(5, 5)
    120
    >>> ff(x, 5) == x*(x-1)*(x-2)*(x-3)*(x-4)
    True
    >>> ff(Poly(x**2, x), 2)
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
            if k.is_zero:
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
                        if isinstance(x, Poly):
                            gens = x.gens
                            if len(gens)!= 1:
                                raise ValueError("ff only defined for "
                                            "polynomials on one generator")
                            else:
                                return reduce(lambda r, i:
                                              r*(x.shift(-i)),
                                              range(0, int(k)), 1)
                        else:
                            return reduce(lambda r, i: r*(x - i),
                                          range(0, int(k)), 1)
                else:
                    if x is S.Infinity:
                        return S.Infinity
                    elif x is S.NegativeInfinity:
                        return S.Infinity
                    else:
                        if isinstance(x, Poly):
                            gens = x.gens
                            if len(gens)!= 1:
                                raise ValueError("rf only defined for "
                                            "polynomials on one generator")
                            else:
                                return 1/reduce(lambda r, i:
                                                r*(x.shift(i)),
                                                range(1, abs(int(k)) + 1), 1)
                        else:
                            return 1/reduce(lambda r, i: r*(x + i),
                                            range(1, abs(int(k)) + 1), 1)

    def _eval_rewrite_as_gamma(self, x, k, **kwargs):
        from sympy import gamma
        return (-1)**k*gamma(k - x) / gamma(-x)

    def _eval_rewrite_as_RisingFactorial(self, x, k, **kwargs):
        return rf(x - k + 1, k)

    def _eval_rewrite_as_binomial(self, x, k, **kwargs):
        if k.is_integer:
            return factorial(k) * binomial(x, k)

    def _eval_rewrite_as_factorial(self, x, k, **kwargs):
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
    r"""Implementation of the binomial coefficient.

    In a strict combinatorial sense, the binomial coefficient gives
    the number of ways we can choose ``k`` elements from a set of
    ``n`` elements. In this case both arguments are nonnegative
    integers and the value of the binomial can be computed in terms
    of factorials::

    .. math:: \binom{n}{k} = \frac{n!}{k!(n - k)!}

    For arbitrary ``n`` and integer ``k``, we use Newton's Generalized
    Binomial Theorem. See [4]_ ::

    .. math:: \binom{n}{k} = \frac{ff(n, k)}{k!}

    .. math:: (x+y)^n & =\sum_{k=0}^\infty \binom{n}{k} x^{n-k} y^k

    Using the aforementioned theorem, one can interpret the binomial
    coefficient :math:`\binom{n}{k}` as the coefficient
    of :math:`x^k, y^k, x^(n - k), y^(n - k)` in the Series expansion
    of :math:`(x + y)^n`.

    When ``k`` is negative and ``n - k`` is a non-negative integer
    then the following identity is used before the coefficient is
    calculated::

    .. math:: \binom{n}{k} = \binom{n}{n - k}

    .. math:: \binom{n}{k} = \frac{ff(n, n - k)}{(n - k)!}

    The binomial coefficient is zero in the following cases, provided
    both ``n`` and ``k`` are integers. See [3]_ ::

    .. math:: k > n \geq 0

    .. math:: n \geq 0 > k

    .. math:: 0 > k > n

    If neither ``k`` nor ``n - k`` is a nonnegative integer then
    the binomial is expressed in terms of the gamma function. See [2]_ ::

    .. math:: \binom{n}{k} = \frac{\gamma(n + 1)}{\gamma(k + 1)\gamma(n - k + 1)}


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
    ...     print([binomial(N, i) for i in range(N + 1)])
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
    >>> [binomial(N, i) for i in range(1 - N)]
    [1, -4, 10, -20, 35]

    >>> binomial(Rational(5, 4), 3)
    -5/128
    >>> binomial(Rational(-5, 4), 3)
    -195/128

    >>> b = binomial(n, 3); b
    binomial(n, 3)

    When ``k`` is positive the polynomial expression
    can be viewed in either expanded or factored form:

    >>> b.expand(func=True)
    n**3/6 - n**2/2 + n/3

    >>> expand_func(b)
    n*(n - 2)*(n - 1)/6


    References
    ==========

    .. [1] https://www.johndcook.com/blog/binomial_coefficients/

    .. [2] http://functions.wolfram.com/GammaBetaErf/Binomial/02/

    .. [3] https://core.ac.uk/download/pdf/82732883.pdf#page=9

    .. [4] https://en.wikipedia.org/wiki/Binomial_theorem#Newton%27s_generalized_binomial_theorem

    See Also
    ========

    sympy.ntheory.multinomial.binomial_coefficients

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
    def _eval(cls, n, k):
        # n.is_Number and k.is_Integer and k != 1 and n != k
        # and k > 0 (and, if n is positive, k < n - k)
        Ntype = Float if k.has(Float) or n.has(Float) else Integer
        k = as_int(k, strict=False)
        direct = n.is_nonnegative
        if direct:
            try:
                n = as_int(n, strict=False)
            except ValueError:
                direct = False
        if direct:

            if HAS_GMPY:
                from sympy.core.compatibility import gmpy
                return Ntype(gmpy.bincoef(n, k))

            d, result = n - k, 1
            for i in range(1, k + 1):
                d += 1
                result = result * d // i
            return Ntype(result)
        else:
            return Mul(*[n + i for i in range(1 - k, 1)])/factorial(k)

    @classmethod
    def eval(cls, n, k):
        from sympy import Min
        n, k = map(sympify, (n, k))
        if k == 1:
            return n
        n_k = n - k
        if n_k == 1:
            return n
        if k.is_zero:
            return S.One
        if n_k.is_zero:
            return S.One
        if -1 in (k, n_k):
            if (n + 1).is_zero is False:
                # includes n.is_integer=False, too
                # since adding 1 to it cannot give 0
                return S.Zero
        if all(int_like(x) for x in (n, k)):
            if (n*k).is_positive and n_k.is_negative: # 0 > k > n
                return S.Zero
            if n.is_nonnegative and k.is_negative: # n >= 0 > k
                return S.Zero
            if n.is_nonnegative and n_k.is_negative: # k > n >= 0
                return S.Zero
        elif n.is_negative:
            if int_like(n) and int_like(k) is False:
                return S.ComplexInfinity
            elif int_like(n) is False:
                if k.is_negative and int_like(k):
                    return S.Zero
                if n_k.is_negative and int_like(n_k):
                    return S.Zero

        if n.is_positive and k.is_positive:
            kk = Min(k, n_k)
            if kk.is_positive:
                k = kk
        elif n.is_nonzero and k.is_negative:
            k = n_k

        if k.is_number:
            if int_like(k) and n.is_number:
                res = cls._eval(n, k)
                return _mexpand(res) if res else res
            elif int_like(k) is False:
                from sympy import gamma
                rv = gamma(n + 1)/(gamma(k + 1)*gamma(n - k + 1))
                if not rv.has(gamma):
                     return rv

    def _eval_Mod(self, q):
        n, k = self.args

        if any(int_like(x) is False for x in (n, k, q)):
            raise ValueError("Integers expected for binomial Mod")

        Ntype = Integer if not any(i.atoms(Float) for i in (n,k,q)) else Float

        if all(int_like(x) and x.is_number for x in (n, k)):
            n, k = map(int, (n, k))

            # handle negative integers k or n
            if k < 0:
                if n > -1 or k > n:
                    return S.Zero
                k = n - k  # nonneg since n <= -1 and k <= n

            if n < 0:
                neg = True
                n = -n + k - 1
                res = -1 if k%2 else 1
            else:
                neg = False
                res = 1

            # nonnegative n and k
            if k > n:
                return S.Zero
            if n == k:
                return Ntype(-1 if neg else 1) % q
            if k == 1:
                return Ntype(-1 if neg else 1)*n % q

            # nonnegative n and k and res
            if int_like(q) and q.is_number:
                aq = abs(q)
                isprime = aq.is_prime
                aq = int(aq)
                if isprime:
                    if aq < n:
                        # use Lucas Theorem
                        N, K = n, k
                        while N or K:
                            res = res*binomial(N % aq, K % aq) % aq
                            N, K = N // aq, K // aq
                    else:
                        # use Factorial Modulo
                        d = n - k
                        kf = 1
                        for i in range(2, k + 1):
                            kf = kf*i % aq
                        df = kf
                        for i in range(k + 1, d + 1):
                            df = df*i % aq
                        res *= df
                        for i in range(d + 1, n + 1):
                            res = res*i % aq

                        res *= pow(kf*df % aq, aq - 2, aq)
                        res %= aq

                else:
                    # Binomial Factorization is performed by calculating the
                    # exponents of primes <= n in `n! /(k! (n - k)!)`,
                    # for non-negative integers n and k. As the exponent of
                    # prime in n! is e_p(n) = [n/p] + [n/p**2] + ...
                    # the exponent of prime in binomial(n, k) would be
                    # e_p(n) - e_p(k) - e_p(n - k)
                    M = int(_sqrt(n))
                    if n - k < k:
                        k = n - k
                    for prime in sieve.primerange(2, n + 1):
                        if prime > n - k:
                            res = res*prime % aq
                        elif prime > n // 2:
                            continue
                        elif prime > M:
                            if n % prime < k % prime:
                                res = res*prime % aq
                        else:
                            N, K = n, k
                            exp = a = 0

                            while N > 0:
                                a = int((N % prime) < (K % prime + a))
                                N, K = N // prime, K // prime
                                exp += a

                            if exp > 0:
                                res *= pow(prime, exp, aq)
                                res %= aq

            return Ntype(res % q)  # not sure how to track precision

    def _eval_expand_func(self, **hints):
        """return binomial(n, k) as a number or a product
        if the summation limit is a nonnegative integer.
        """
        from sympy import Product
        n, k = self.args
        # see if re-evaluation evaluates
        e = self.func(n, k)
        if e.func != self.func:
            return e

        # see if there is a non-negative k or n - k
        if k.is_number and int_like(k) and k.is_nonnegative:
            pass
        elif (n - k).is_number and int_like(n - k) and (n - k).is_nonnegative:
            k = n - k
        else:
            return self
        i = Dummy('never_seen')
        return (Product(n - k + i, (i, 1, k))/factorial(k)).doit()

    def _eval_rewrite_as_Piecewise(self, n, k, **kwargs):
        from sympy import gamma, Abs, Lt
        isint = lambda x: Eq(x, floor(x))
        nneg = lambda x: Or(Eq(x, 0), Eq(x/Abs(x), 1))
        for i in igen():
            if not n.has(i) and not k.has(i):
                break
        return Piecewise(
            (factorial(n)/factorial(k)/factorial(n - k),
                And(isint(n), isint(k),
                nneg(n), nneg(k), nneg(n - k))),
            (ff(n, k)/factorial(k),
                And(isint(k), nneg(k))),
            (ff(n, n - k)/factorial(n - k),
                And(isint(n - k), nneg(n - k))),
            (0, And(*[i for x in (n, k, n - k)
                for i in [isint(x), Lt(x, 0)]])),
            (gamma(n + 1)/gamma(k + 1)/gamma(n - k + 1),
                True))

    def _eval_rewrite_as_factorial(self, n, k, **kwargs):
        rv = self._eval_rewrite_as_gamma(n, k)
        if rv is not None:
            return rv.rewrite(factorial)

    def _eval_rewrite_as_gamma(self, n, k, **kwargs):
        from sympy import gamma, Pow
        ki = k.is_integer
        nki = (n - k).is_integer
        if ki is None:
            return
        if nki is None:
            return
        if n.is_integer and k.is_integer:
            if n.is_negative:
                if k.is_nonnegative:
                    return(Pow(-1, k)*gamma(k - n)/(gamma(-n)*gamma(k + 1)))
                if k.is_negative:
                    return (Pow(-1, n - k)*gamma(-k)/(gamma(-n)*gamma(n - k + 1)))
            elif n.is_nonnegative:
                if k.is_negative:
                    return S.Zero
                if k.is_nonnegative:
                    return (gamma(n + 1)/(gamma(k + 1)*gamma(n - k + 1)))
        elif ki is False and nki is False:
            return (gamma(n + 1)/(gamma(k + 1)*gamma(n - k + 1)))
        elif ki and (k + 1).is_positive:
            return (gamma(n + 1)/(gamma(k + 1)*gamma(n - k + 1)))
        elif nki and (n - k + 1).is_positive:
            return (gamma(n + 1)/(gamma(k + 1)*gamma(n - k + 1)))

    def _eval_rewrite_as_tractable(self, n, k, **kwargs):
        g = self._eval_rewrite_as_gamma(n, k)
        if g is not None:
            return g.rewrite('tractable')

    def _eval_rewrite_as_FallingFactorial(self, n, k, **kwargs):
        if int_like(k) and k.is_nonnegative:
            return ff(n, k) / factorial(k)
        if int_like(n - k) and (n - k).is_nonnegative:
            return ff(n, n - k) / factorial(n - k)

    def _eval_is_integer(self):
        n, k = self.args
        if k.is_integer:
            if n.is_integer or k.is_negative:
                return True
                # and if not, it *may* be true but it might
                # only be obvious after expansion, e.g.
                # (3 + sqrt(11)*I, 3) = -10
        elif (n - k).is_integer and (n - k).is_negative:
            # (x, -int) == 0
            return True
        elif k.is_infinite:
            if k.is_extended_positive:
                if (n + 1).is_zero:
                    # (-1, oo) -> zoo
                    return False
                elif (n + 1).is_zero is False:
                    if n.is_infinite:
                        if n.is_extended_negative:
                            return True
                        elif n.is_extended_positive:
                            return False
                    elif n.is_finite:
                        return True
        elif k.is_real and k.is_integer is False and (
                n - k).is_integer is False:
            return False

    def _eval_is_nonnegative(self):
        n, k = self.args
        if (n - k).is_zero or k.is_zero or n.is_nonnegative:
            return True
        if n.is_negative:
            if k.is_negative:
                if (n - k).is_positive:
                    if (k.is_even == n.is_even) and k.is_even is not None:
                        return True
                    elif (k.is_even != n.is_even) and (
                        k.is_even is not None) and (n.is_even is not None):
                        return False
                else:
                    return True
            elif k.is_positive:
                if k.is_odd:
                    return False
                elif k.is_even:
                    return True
