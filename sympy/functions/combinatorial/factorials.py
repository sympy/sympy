
from sympy.core.basic import Basic, S, cache_it, cache_it_immutable
from sympy.core.function import DefinedFunction, Apply, Lambda
from sympy.ntheory import sieve
from math import sqrt

###############################################################################
######################## FACTORIAL and MULTI-FACTORIAL ########################
###############################################################################

class Factorial(DefinedFunction):
    """Implementation of factorial function over nonnegative integers.
       For the sake of convenience and simplicity of procedures using
       this function it is defined for negative integers and returns
       zero in this case.

       The factorial is very important in combinatorics where it gives
       the number of ways in which 'n' objects can be permuted. It also
       arises in calculus, probability, number theory etc.

       There is strict relation of factorial with gamma function. In
       fact n! = gamma(n+1) for nonnegarive integers. Rewrite of this
       kind is very useful in case of combinatorial simplification.

       Computation of the factorial is done using two algorithms. For
       small arguments naive product is evaluated. However for bigger
       input algorithm Prime-Swing is used. It is the fastest algorithm
       known and computes n! via prime factorization of special class
       of numbers, called here the 'Swing Numbers'.

       >>> from sympy import *
       >>> n = Symbol('n', integer=True)

       >>> factorial(-2)
       0

       >>> factorial(0)
       1

       >>> factorial(7)
       5040

       >>> factorial(n)
       n!

       >>> factorial(2*n)
       (2*n)!

    """

    nofargs = 1

    _small_swing = [
        1,1,1,3,3,15,5,35,35,315,63,693,231,3003,429,6435,6435,109395,
        12155,230945,46189,969969,88179,2028117,676039,16900975,1300075,
        35102025,5014575,145422675,9694845,300540195,300540195
    ]

    def _swing(self, n):
        if n < 33:
            return self._small_swing[n]
        else:
            N, primes = int(sqrt(n)), []

            for prime in sieve.primerange(3, N+1):
                p, q = 1, n

                while True:
                    q /= prime

                    if q > 0:
                        if q & 1 == 1:
                            p *= prime
                    else:
                        break

                if p > 1:
                    primes.append(p)

            for prime in sieve.primerange(N+1, n/3 + 1):
                if (n / prime) & 1 == 1:
                    primes.append(prime)

            L_product = R_product = 1

            for prime in sieve.primerange(n/2 + 1, n+1):
                L_product *= prime

            for prime in primes:
                R_product *= prime

            return L_product*R_product

    def _recursive(self, n):
        if n < 2:
            return 1
        else:
            return (self._recursive(n/2)**2)*self._swing(n)

    def _eval_apply(self, n):
        n = Basic.sympify(n)

        if isinstance(n, Basic.Number):
            if isinstance(n, Basic.Zero):
                return S.One
            elif isinstance(n, Basic.Integer):
                if n.is_negative:
                    return S.Zero
                else:
                    n, result = n.p, 1

                    if n < 20:
                        for i in range(2, n+1):
                            result *= i
                    else:
                        N, bits = n, 0

                        while N != 0:
                            if N & 1 == 1:
                                bits += 1

                            N = N >> 1

                        result = self._recursive(n)*2**(n-bits)

                    return Basic.Integer(result)

        if n.is_integer:
            if n.is_negative:
                return S.Zero
        else:
            return S.Gamma(n+1)

class ApplyFactorial(Apply):

    def _eval_expand_func(self, *args):
        return S.Gamma(1 + self.args[0])

    def tostr(self, level=0):
        return '%s!' % self.args[0].tostr(self.precedence)

    def _eval_is_integer(self):
        return self.args[0].is_integer

Basic.singleton['factorial'] = Factorial

###############################################################################
######################## RISING and FALLING FACTORIALS ########################
###############################################################################

class RisingFactorial(DefinedFunction):
    """Rising factorial (also called Pochhammer symbol) is a double valued
       function arising in concrete mathematics, hypergeometric functions
       and series expanansions. It is defined by

                   rf(x, k) = x * (x+1) * ... * (x + k-1)

       where 'x' can be arbitrary expression and 'k' is an integer. For
       more information check "Concrete mathematics" by Graham, pp. 66
       or visit http://mathworld.wolfram.com/RisingFactorial.html page.

       >>> from sympy import *
       >>> x = Symbol('x')

       >>> rf(x, 0)
       1

       >>> rf(1, 5)
       120

       >>> rf(x, 5)
       x*(1 + x)*(2 + x)*(3 + x)*(4 + x)

    """

    nofargs = 2

    def _eval_apply(self, x, k):
        x = Basic.sympify(x)
        k = Basic.sympify(k)

        if isinstance(x, Basic.NaN):
            return S.NaN
        elif isinstance(x, Basic.One):
            return S.Factorial(k)
        elif isinstance(k, Basic.Integer):
            if isinstance(k, Basic.NaN):
                return S.NaN
            if isinstance(k, Basic.Zero):
                return S.One
            else:
                if k.is_positive:
                    if isinstance(x, Basic.Infinity):
                        return S.Infinity
                    elif isinstance(x, Basic.NegativeInfinity):
                        if k.is_odd:
                            return S.NegativeInfinity
                        else:
                            return S.Infinity
                    else:
                        return reduce(lambda r, i: r*(x+i), xrange(0, int(k)), 1)
                else:
                    if isinstance(x, Basic.Infinity):
                        return S.Infinity
                    elif isinstance(x, Basic.NegativeInfinity):
                        return S.Infinity
                    else:
                        return 1/reduce(lambda r, i: r*(x-i), xrange(1, abs(int(k))+1), 1)

class ApplyRisingFactorial(Apply):

    def _eval_expand_func(self, *args):
        x, k = self.args[0], self.args[1]
        return S.Gamma(x+k)/S.Gamma(x)

    def tostr(self, level=0):
        r = 'rf(%s)' % ', '.join([a.tostr() for a in self.args])

        if self.precedence <= level:
            return '(%s)' % (r)
        else:
            return r

class FallingFactorial(DefinedFunction):
    """Falling factorial (related to rising factorial) is a double valued
       function arising in concrete mathematics, hypergeometric functions
       and series expanansions. It is defined by

                   ff(x, k) = x * (x-1) * ... * (x - k+1)

       where 'x' can be arbitrary expression and 'k' is an integer. For
       more information check "Concrete mathematics" by Graham, pp. 66
       or visit http://mathworld.wolfram.com/FallingFactorial.html page.

       >>> from sympy import *
       >>> x = Symbol('x')

       >>> ff(x, 0)
       1

       >>> ff(5, 5)
       120

       >>> ff(x, 5)
       x*(1 - x)*(2 - x)*(3 - x)*(4 - x)

    """

    nofargs = 2

    def _eval_apply(self, x, k):
        x = Basic.sympify(x)
        k = Basic.sympify(k)

        if isinstance(x, Basic.NaN):
            return S.NaN
        elif isinstance(k, Basic.Integer):
            if isinstance(k, Basic.NaN):
                return S.NaN
            if isinstance(k, Basic.Zero):
                return S.One
            else:
                result = S.One

                if k.is_positive:
                    if isinstance(x, Basic.Infinity):
                        return S.Infinity
                    elif isinstance(x, Basic.NegativeInfinity):
                        if k.is_odd:
                            return S.NegativeInfinity
                        else:
                            return S.Infinity
                    else:
                        return reduce(lambda r, i: r*(x-i), xrange(0, int(k)), 1)
                else:
                    if isinstance(x, Basic.Infinity):
                        return S.Infinity
                    elif isinstance(x, Basic.NegativeInfinity):
                        return S.Infinity
                    else:
                        return 1/reduce(lambda r, i: r*(x+i), xrange(1, abs(int(k))+1), 1)

class ApplyFallingFactorial(Apply):

    def _eval_expand_func(self, *args):
        x, k = self.args[0], self.args[1]
        return (-1)**k*S.Gamma(-x+k)/S.Gamma(-x)

    def tostr(self, level=0):
        r = 'ff(%s)' % ', '.join([a.tostr() for a in self.args])

        if self.precedence <= level:
            return '(%s)' % (r)
        else:
            return r

Basic.singleton['rf'] = RisingFactorial
Basic.singleton['ff'] = FallingFactorial

###############################################################################
########################### BINOMIAL COEFFICIENTS #############################
###############################################################################

class Binomial(DefinedFunction):
    """Implementation of the binomial coefficient. It can be defined
       in two ways depending on its desired interpretation:

           C(n,k) = n!/(k!(n-k)!)   or   C(n, k) = ff(n, k)/k!

       First formula has strict combinatorial meaning, definig the
       number of ways we can choose 'k' elements from 'n' element
       set. In this case both arguments are nonnegative integers
       and binomial is computed using efficient algorithm based
       on prime factorisation.

       The other definition is generalisation for arbitaty 'n',
       however 'k' must be also nonnegative. This case is very
       useful in case of evaluating summations.

       For the sake of convenience for negative 'k' this function
       will return zero no matter what valued is the other argument.

       >>> from sympy import *
       >>> n = symbols('n', integer=True)

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

       >>> binomial(Rational(5,4), 3)
       -5/128

       >>> binomial(n, 3)
       (1/6)*n*(1 - n)*(2 - n)

    """

    nofargs = 2

    def _eval_apply(self, r, k):
        r, k = map(Basic.sympify, (r, k))

        if isinstance(k, Basic.Number):
            if isinstance(k, Basic.Zero):
                return S.One
            elif isinstance(k, Basic.Integer):
                if k.is_negative:
                    return S.Zero
                else:
                    if isinstance(r, Basic.Integer) and r.is_nonnegative:
                        r, k = int(r), int(k)

                        if k > r:
                            return S.Zero
                        elif k > r / 2:
                            k = r - k

                        M, result = int(sqrt(r)), 1

                        for prime in sieve.primerange(2, r+1):
                            if prime > r - k:
                                result *= prime
                            elif prime > r / 2:
                                continue
                            elif prime > M:
                                if r % prime < k % prime:
                                    result *= prime
                            else:
                                R, K = r, k
                                exp = a = 0

                                while R > 0:
                                    a = int((R % prime) < (K % prime + a))
                                    R, K, exp = R / prime, K / prime, exp + a

                                if exp > 0:
                                    result *= prime**exp

                        return Basic.Integer(result)
                    else:
                        result = r - k + 1

                        for i in xrange(2, k+1):
                            result *= r-k+i
                            result /= i

                        return result

        if k.is_integer:
            if k.is_negative:
                return S.Zero
        else:
            return S.Gamma(r+1)/(S.Gamma(r-k+1)*S.Gamma(k+1))

class ApplyBinomial(Apply):

    def _eval_expand_func(self, *args):
        r, k = self.args[0], self.args[1]
        return S.Gamma(r+1)/(S.Gamma(r-k+1)*S.Gamma(k+1))

    def _eval_is_integer(self):
        return self.args[0].is_integer and self.args[1].is_integer

Basic.singleton['binomial'] = Binomial
