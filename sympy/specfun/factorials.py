
from sympy.core import *

# Factorial and gamma related functions

def sqrt(arg):
    return arg**(Rational(1,2))


# Lanczos approximation for low-precision numerical factorial
# This implementation is not particularly numerically stable
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


class Factorial_(DefinedFunction):
    """
    Usage
    =====
        factorial(x) -> Returns the factorial of x, defined as
        x! = 1*2*3*...*x if x is a positive integer

    Notes
    =====
        factorial(x) is evaluated explicitly if x is an integer or
        half an integer. If x is a negative integer, the value is
        infinite.

    Examples
    ========
        >>> from sympy import *
        >>> from sympy.specfun.factorials import *
        >>> factorial(5)
        120
        >>> factorial(0)
        1
        >>> factorial(Rational(5,2))
        15/8*pi**(1/2)

    """
    nofargs = 1

    def _eval_apply(self, x):
        if isinstance(x, Rational):
            if x.is_integer:
                if x < 0:
                    return oo
                y = 1
                for m in xrange(1, x.p+1):
                    y *= m
                return Rational(y)
            if x.q == 2:
                n = (x.p + 1) / 2
                if n < 0:
                    return (-1)**(-n+1) * pi * x / factorial_(-x)
                return sqrt(pi) * Rational(1, 2**n) * factorial2(2*n-1)

    def diff(self, sym):
        return gamma(self._args+1).diff(sym)

    def evalf(self):
        """Return a low-precision approximation of self."""
        a, b = self._args.get_re_im()
        y = _lanczos(complex(a, b))
        return Real(y.real) + I*Real(y.imag)

    # This should give a series expansion around x = oo. Needs fixing
    def _series(self, x, n):
        return sqrt(2*pi*x) * x**x * exp(-x) * (1 + O(1/x))

    def __str__(self):
        x = self._args
        if (isinstance(x, Rational) and x.is_integer and x >= 0) or \
            isinstance(x, Symbol):
            s = str(x)
        else:
            s = "(" + str(x) + ")"
        return s + "!"

    def __latex__(self):
        x = self._args
        if (isinstance(x, Rational) and x.is_integer and x >= 0) or \
            isinstance(x, Symbol):
            s = x.__latex__()
        else:
            s = "(" + x.__latex__() + ")"
        return s + "!"


class UnevalatedFactorial(Factorial_):
    def _eval_apply(self, x):
        return None

unfac = UnevalatedFactorial()


class Factorial2(DefinedFunction):
    """
    Usage
    =====
        factorial2(x) -> Returns the double factorial of x, defined as
        x!! = 2*4*6*...*x if x is a positive even integer and as
        x!! = 1*3*5*...*x if x is a positive odd integer.

    Notes
    =====
        Also defined for negative odd integers, but infinite for
        negative even integers.

    Examples
    ========
        >>> from sympy import *
        >>> from sympy.specfun.factorials import *
        >>> factorial2(5)
        15
        >>> factorial2(6)
        48

    """
    nofargs = 1

    def _eval_apply(self, x):
        if isinstance(x, Rational) and x.is_integer:
            if int(x) % 2 == 0:
                if x < 0:
                    return oo
                else:
                    return 2**(x/2) * factorial(x/2)
            else:
                if x < 0:
                    return factorial2(x+2) / (x+2)
                else:
                    return factorial(x) / 2**((x-1)/2) / factorial((x-1)/2)

    def __latex__(self):
        x = self._args
        if (isinstance(x, Rational) and x.is_integer and x >= 0) or \
            isinstance(x, Symbol):
            s = x.__latex__()
        else:
            s = "(" + x.__latex__() + ")"
        return s + "!!"


# factorial_simplify helpers; could use refactoring

def _isfactorial(expr):
    return isinstance(expr, Apply) and isinstance(expr[0], Factorial_)

def _collect_factors(expr):
    assert isinstance(expr, Mul)
    numer_args = []
    denom_args = []
    other = []
    for x in expr:
        if isinstance(x, Mul):
            n, d, o = _collect_factors(x)
            numer_args += n
            denom_args += d
            other += o
        elif isinstance(x, Pow):
            base, exp = x[:]
            if _isfactorial(base) and \
                isinstance(exp, Rational) and exp.is_integer:
                if exp > 0:
                    for i in xrange(exp.p): numer_args.append(base.args[0])
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
                if isinstance(delta, Rational) and delta.is_integer:
                    candidates.append((delta, i, j))
        if candidates:
            # There may be multiple candidates. Choose the quotient pair
            # that minimizes the work.
            candidates.sort(key=lambda x: abs(x[0]))
            delta, i, j = candidates[0]
            p = na[i]
            q = da[j]
            t = Rational(1)
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
    
    if isinstance(expr, Add):
        return Add(*(factorial_simplify(x) for x in expr))

    if isinstance(expr, Factorial_):
        #return expr.eval()
        return expr

    if isinstance(expr, Pow):
        return Pow(factorial_simplify(expr[0]), expr[1])

    if isinstance(expr, Mul):
        na, da, other = _collect_factors(expr)
        _simplify_quotient(na, da, other)
        _simplify_recurrence(na, other)
        _simplify_recurrence(da, other, reciprocal=True)

        result = Rational(1)
        for n in na: result *= factorial_(n)
        for d in da: result /= factorial_(d)
        for o in other: result *= o
        return result

    expr = expr.subs(unfac, factorial_)

    return expr

class Rising_factorial(DefinedFunction):
    """
    Usage
    =====
        Calculate the rising factorial (x)^(n) = x(x+1)...(x+n-1)
        as a quotient of factorials

    Examples
    ========
        >>> from sympy.specfun.factorials import *
        >>> rising_factorial(3, 2)
        12

    """
    nofargs = 2

    def _eval_apply(self, x, n):
        return factorial_simplify(unfac(x+n-1) / unfac(x-1))

    def __latex__(self):
        x, n = self._args
        return "{(%s)}^{(%s)}" % (x.__latex__(), n.__latex__())


class Falling_factorial(DefinedFunction):
    """
    Usage
    =====
        Calculate the falling factorial (x)_(n) = x(x-1)...(x-n+1)
        as a quotient of factorials

    Examples
    ========
        >>> from sympy.specfun.factorials import *
        >>> falling_factorial(5, 3)
        60

    """
    nofargs = 2

    def _eval_apply(self, x, n):
        return factorial_simplify(unfac(x) / unfac(x-n))

    def __latex__(self):
        x, n = self._args
        return "{(%s)}_{(%s)}" % (x.__latex__(), n.__latex__())


class Binomial2(DefinedFunction):
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
        >>> from sympy import *
        >>> from sympy.specfun.factorials import *
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
        1/6*x*(-2+x)*(-1+x)

    """
    nofargs = 2

    def _eval_apply(self, n, k):

        # TODO: move these two cases to factorial_simplify as well
        if n == 0 and k != 0:
            return sin(pi*k)/(pi*k)

        return factorial_simplify(unfac(n) / unfac(k) / unfac(n-k))

    def __latex__(self):
        n, k = self._args
        return r"{{%s}\choose{%s}}" % (n.__latex__(), k.__latex__())


class Gamma(DefinedFunction):
    """
    Usage
    =====
        gamma(x) -> calculate the gamma function of x

    Notes
    =====
        gamma(x) = (x-1)!

        When x is an integer or half-integer, the result is automatically
        simplified to the corresponding factorial

    Examples
    ========
        >>> from sympy import *
        >>> from sympy.specfun.factorials import *
        >>> gamma(3)
        2
        >>> gamma(Rational(1,2))
        pi**(1/2)

    """
    nofargs = 1

    def _eval_apply(self, x):
        y = factorial_(x-1)
        try:
            if not isinstance(y.func, Factorial_):
                return y
        except:
            return y

    def diff(self, sym):
        from zeta_functions import polygamma
        x = self._args
        return gamma(x)*polygamma(0,x)*x.diff(sym)

    def __latex__(self):
        return "\Gamma(" + self._args.__latex__() + ")"


class LowerGamma(DefinedFunction):
    """
    Lower incomplete gamma function

    gamma(a, x)
    """
    nofargs = 2

    def _eval_apply(self, a, x):
        if a == 1:
            return 1 - exp(-x)
        if a.is_integer and a > 1:
            b = a-1
            return b*lower_gamma(b, x) - x**b * exp(-x)
        #return self


class UpperGamma(DefinedFunction):
    """
    Upper incomplete gamma function

    Gamma(a, x)
    """
    nofargs = 2

    def _eval_apply(self, a, x):
        if x == 0:
            return gamma(a)
        if a == 1:
            return exp(-x)
        if a.is_integer and a > 1:
            b = a-1
            return b*upper_gamma(b, x) + x**b * exp(-x)
        #return self


factorial_ = Factorial_()
factorial2 = Factorial2()
rising_factorial = Rising_factorial()
falling_factorial = Falling_factorial()
binomial2 = Binomial2()
upper_gamma = UpperGamma()
lower_gamma = LowerGamma()
gamma = Gamma()