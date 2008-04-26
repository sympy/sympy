from sympy.core import Function, Add, Pow, Mul, Rational, Integer, pi, \
        oo, Real, Symbol, sympify
from sympy.core.basic import S, C

from sympy.functions.elementary.miscellaneous import sqrt
from sympy.utilities.memoization import recurrence_memo

# Lanczos approximation for low-precision numerical factorial
_lanczos_coef = [0.99999999999980993, 676.5203681218851, -1259.1392167224028,
  771.32342877765313, -176.61502916214059, 12.507343278686905,
    -0.13857109526572012, 9.9843695780195716e-6, 1.5056327351493116e-7]

def _lanczos(z):
    from cmath import pi, sin, log, exp
    if z.real < 0:
        return pi*z / (sin(pi*z) * _lanczos(-z))
    else:
        x = _lanczos_coef[0]
        for i in range(1, 9):
            x += _lanczos_coef[i]/(z+i)
        logw = 0.91893853320467267+(z+0.5)*log(z+7.5)+log(x)-z-7.5
        return exp(logw)

class _Factorial(Function):
    """
    Factorials and multiple factorials

    Usage
    =====
        factorial(x) gives the factorial of x, defined as
        x! = 1*2*3*...*x if x is a positive integer. If x is not a
        positive integer, the value of x! is defined in terms of the
        gamma function.

        factorial(x, m) returns the m-th order multifactorial of x;
        i.e., 1 * ... * (x-2*m) * (x-m) * x. In particular,
        factorial(x, 2) returns the double factorial of x,

          * x!! = 2*4*...*(x-2)*x if x is a positive even integer
          * x!! = 1*3*...*(x-2)*x if x is a positive odd integer

        The argument m must be a positive integer.

    Examples
    ========
        >>> factorial(5)
        120
        >>> factorial(0)
        1
        >>> factorial(Rational(5,2))
        (15/8)*pi**(1/2)

        >>> factorial(5, 2)
        15
        >>> factorial(6, 2)
        48

    """

    @staticmethod
    @recurrence_memo([1])
    def _fac1(n, prev):
        return n * prev[-1]

    # generators for multifactorials
    _generators = {}

    @classmethod
    def canonize(cls, x, m=1):
        m = sympify(m)

        # the usual case
        if m == 1 and x.is_integer:
            # handle factorial poles, also for symbols with assumptions
            if x.is_negative:
                return oo
            if x.is_Integer:
                return Integer(cls._fac1(int(x)))

        # half-integer case of the ordinary factorial
        if m == 1 and x.is_Rational and x.q == 2:
            n = (x.p + 1) / 2
            if n < 0:
                return (-1)**(-n+1) * pi * x / factorial(-x)
            return sqrt(pi) * Rational(1, 2**n) * factorial(2*n-1, 2)

        # multifactorials are only defined for integers
        if (not m.is_Integer and m > 0) or not x.is_Integer:
            return

        x = int(x)
        m = int(m)

        if x == 0:
            return S.One

        # e.g., for double factorials, the start value is 1 for odd x
        # and 2 for even x
        start_value = x % m
        if start_value == 0:
            start_value = m

        # we can define the m-multifactorials for negative numbers
        # to satisfy f(x,m) = f(x+m,m) / (x+m)
        if x < 0:
            if start_value == m:
                return oo
            return factorial(x+m, m) / (x+m)

        # positive case
        if not (m, start_value) in cls._generators:
            @recurrence_memo([start_value])
            def _f(k, prev):
                return (k*m + start_value) * prev[-1]
            cls._generators[m, start_value] = _f

        return Integer(cls._generators[m, start_value]((x-start_value) // m))

    # This should give a series expansion around x = oo. Needs fixing
    # def series(self, x, n):
    #    return sqrt(2*pi*x) * x**x * exp(-x) * (1 + O(1/x))

    # Derivatives are given in terms of polygamma functions
    # (XXX: only for order m = 1)
    def fdiff(self, argindex=1):
        if argindex == 1:
            from sympy.functions import gamma, polygamma
            return gamma(self.args[0]+1)*polygamma(0,self.args[0]+1)
        else:
            raise ArgumentIndexError(self, argindex)

    @classmethod
    def _eval_apply_evalf(cls, x, m=1):
        """Return a low-precision numerical approximation."""
        assert m == 1
        a, b = x.as_real_imag()
        y = _lanczos(complex(a, b))
        if b == 0:
            return Real(y.real)
        else:
            Real(y.real) + I*Real(y.imag)

factorial = _Factorial

class UnevaluatedFactorial(_Factorial):

    @classmethod
    def canonize(cls, x, m=1):
        return None

unfac = UnevaluatedFactorial



# factorial_simplify helpers; could use refactoring

def _isfactorial(expr):
    return isinstance(expr, _Factorial)

def _collect_factors(expr):
    assert isinstance(expr, Mul)
    numer_args = []
    denom_args = []
    other = []
    for x in expr.args:
        if x.is_Mul:
            n, d, o = _collect_factors(x)
            numer_args += n
            denom_args += d
            other += o
        elif x.is_Pow:
            base, exp = x.args
            if _isfactorial(base) and \
                isinstance(exp, Rational) and exp.is_integer:
                if exp > 0:
                    for i in xrange(exp.p): numer_args.append(base[0])
                else:
                    for i in xrange(-exp.p): denom_args.append(base.args[0])
            else:
                other.append(x)
        elif _isfactorial(x):
            numer_args.append(x.args[0])
        else:
            other.append(x)
    return numer_args, denom_args, other

# handle x!/(x+n)!
def _simplify_quotient(na, da, other):
    while 1:
        candidates = []
        for i, y in enumerate(na):
            for j, x in enumerate(da):
                #delta = simplify(y - x)
                delta = y - x
                if delta.is_Integer:
                    candidates.append((delta, i, j))
        if candidates:
            # There may be multiple candidates. Choose the quotient pair
            # that minimizes the work.
            candidates.sort(key=lambda x: abs(x[0]))
            delta, i, j = candidates[0]
            p = na[i]
            q = da[j]
            t = S.One
            if delta > 0:
                for k in xrange(1, int(delta)+1):
                    t *= (q+k)
            else:
                for k in xrange(1, -int(delta)+1):
                    t /= (p+k)
            other.append(t)
            del na[i], da[j]
        else:
            return

# handle x!*(x+1) and x!/x
def _simplify_recurrence(facs, other, reciprocal=False):
    # this should probably be rewritten more elegantly
    i = 0
    while i < len(facs):
        j = 0
        while j < len(other):
            othr = other[j]
            fact = facs[i]
            if reciprocal:
                othr = 1/othr
            if   othr - fact == 1: facs[i] += 1; del other[j]; j -= 1
            elif -othr - fact == 1: facs[i] += 1; del other[j]; other.append(-1); j -= 1
            elif 1/othr - fact == 0: facs[i] -= 1; del other[j]; j -= 1
            elif -1/othr - fact == 0: facs[i] -= 1; del other[j]; other.append(-1); j -= 1
            j += 1
        i += 1

def factorial_simplify(expr):
    """
    This function takes an expression containing factorials
    and removes as many of them as possible by combining
    products and quotients of factorials into polynomials
    or other simpler expressions.

    TODO: handle reflection formula, duplication formula
    double factorials
    """
    expr = sympify(expr)

    if expr.is_Add:
        return Add(*(factorial_simplify(x) for x in expr))

    if expr.is_Pow:
        return Pow(factorial_simplify(expr.args[0]), expr.args[1])

    if expr.is_Mul:
        na, da, other = _collect_factors(expr)
        _simplify_quotient(na, da, other)
        _simplify_recurrence(na, other)
        _simplify_recurrence(da, other, reciprocal=True)

        result = S.One
        for n in na: result *= factorial(n)
        for d in da: result /= factorial(d)
        for o in other: result *= o
        return result

    expr = expr.subs(unfac, factorial)

    return expr

class Rising_factorial(Function):
    """
    Usage
    =====
        Calculate the rising factorial (x)^(n) = x(x+1)...(x+n-1)
        as a quotient of factorials

    Examples
    ========
        >>> rising_factorial(3, 2)
        12

    """
    nargs = 2

    @classmethod
    def canonize(cls, x, n):
        return factorial_simplify(unfac(x+n-1) / unfac(x-1))


rising_factorial = Rising_factorial


class Falling_factorial(Function):
    """
    Usage
    =====
        Calculate the falling factorial (x)_(n) = x(x-1)...(x-n+1)
        as a quotient of factorials

    Examples
    ========
        >>> falling_factorial(5, 3)
        60

    """
    nargs = 2

    @classmethod
    def canonize(cls, x, n):
        return factorial_simplify(unfac(x) / unfac(x-n))


falling_factorial = Falling_factorial


class Binomial2(Function):
    """
    Usage
    =====
        Calculate the binomial coefficient C(n,k) = n!/(k!(n-k)!)

    Notes
    =====
        When n and k are positive integers, the result is always
        a positive integer

        If k is a positive integer, the result is a polynomial in n
        that is evaluated explicitly.

    Examples
    ========
        >>> binomial2(15,8)
        6435
        >>> # Building Pascal's triangle
        >>> [binomial2(0,k) for k in range(1)]
        [1]
        >>> [binomial2(1,k) for k in range(2)]
        [1, 1]
        >>> [binomial2(2,k) for k in range(3)]
        [1, 2, 1]
        >>> [binomial2(3,k) for k in range(4)]
        [1, 3, 3, 1]
        >>> # n can be arbitrary if k is a positive integer
        >>> binomial2(Rational(5,4), 3)
        -5/128
        >>> x = Symbol('x')
        >>> binomial2(x, 3)
        (1/6)*x*(1 - x)*(2 - x)

    """
    nargs = 2

    @classmethod
    def canonize(cls, n, k):

        # TODO: move these two cases to factorial_simplify as well
        if n == 0 and k != 0:
            return C.sin(pi*k)/(pi*k)

        return factorial_simplify(unfac(n) / unfac(k) / unfac(n-k))

binomial2 = Binomial2
