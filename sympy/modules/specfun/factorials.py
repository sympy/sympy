from sympy.core.function import DefinedFunction
from sympy.core.numbers import Number, Real, Rational, pi, I, oo
from sympy import Symbol, Add, Mul, Pow, Basic
#from sympy.modules.simplify import simplify
#from sympy import O
#from sympy.modules.trigonometric import sin

# Factorial and gamma related functions


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


class Factorial(DefinedFunction):
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
        >>> from sympy.modules.specfun.factorials import *
        >>> factorial(5)
        120
        >>> factorial(0)
        1
        >>> factorial(Rational(5,2))
        15/8*pi**(1/2)

    """
    def _eval_apply(self, args):
        x = args
        if isinstance(x, Rational):
            if x.is_integer:
                if x < 0:
                    return oo
                y = Rational(1)
                for m in range(1, x+1):
                    y *= m
                return y
            if x.q == 2:
                n = (x.p + 1) / 2
                if n < 0:
                    return (-1)**(-n+1) * pi * x / factorial(-x)
                return sqrt(pi) * Rational(1, 2**n) * factorial2(2*n-1)
        return self

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

factorial = Factorial()

def _fac(x):
    return factorial(x, evaluate=False)


class factorial2(DefinedFunction):
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
        >>> from sympy.modules.specfun.factorials import *
        >>> factorial2(5)
        15
        >>> factorial2(6)
        48

    """
    def eval(self):
        x = self._args
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
        return self

    def __latex__(self):
        x = self._args
        if (isinstance(x, Rational) and x.is_integer and x >= 0) or \
            isinstance(x, Symbol):
            s = x.__latex__()
        else:
            s = "(" + x.__latex__() + ")"
        return s + "!!"


# factorial_simplify helpers; could use refactoring

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
            if isinstance(base, factorial) and \
                isinstance(exp, Rational) and exp.is_integer:
                if exp > 0:
                    for i in range(exp): numer_args.append(base._args)
                else:
                    for i in range(-exp): denom_args.append(base._args)
            else:
                other.append(x)
        elif isinstance(x, factorial):
            numer_args.append(x._args)
        else:
            other.append(x)
    return numer_args, denom_args, other

# handle x!/(x+n)!
def _simplify_quotient(na, da, other):
    while 1:
        candidates = []
        for i, y in enumerate(na):
            for j, x in enumerate(da):
                delta = simplify(y - x)
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
                for k in range(1, delta+1):
                    t *= (q+k)
            else:
                for k in range(1, -delta+1):
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
            if reciprocal:
                if simplify(other[j] - facs[i]) == 0:
                    facs[i] -= 1; del other[j]; j = -1
                elif simplify(1/other[j] - facs[i]) == 1:
                    facs[i] += 1; del other[j]; j = -1
            else:
                if simplify(other[j] - facs[i]) == 1:
                    facs[i] += 1; del other[j]; j = -1
                elif simplify(1/other[j] - facs[i]) == 0:
                    facs[i] -= 1; del other[j]; j = -1
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

    if isinstance(expr, factorial):
        return expr.eval()

    if isinstance(expr, Pow):
        return Pow(factorial_simplify(expr[0]), expr[1])

    if isinstance(expr, Mul):
        na, da, other = _collect_factors(expr)

        _simplify_quotient(na, da, other)
        _simplify_recurrence(na, other)
        _simplify_recurrence(da, other, reciprocal=True)

        result = Rational(1)
        for n in na: result *= factorial(n).eval()
        for d in da: result /= factorial(d).eval()
        for o in other: result *= o
        return result

    return expr


# This class is a temporary solution
class Function2(DefinedFunction):

    def __init__(self, x, y):
        Basic.__init__(self, is_commutative=True)
        self._args = self.sympify(x), self.sympify(y)

    def atoms(self, s=[], type=None):
        x, y = self._args

        s_temp = list(set(x.atoms()) | set(y.atoms()))

        if type is not None:
            return filter(lambda x : isinstance(x, type), s_temp)

        return s_temp

    def subs(self, old, new):
        x, y = self._args
        x = x.subs(old, new)
        y = y.subs(old, new)
        return self.__class__(x, y)

class rising_factorial(Function2):
    """
    Usage
    =====
        Calculate the rising factorial (x)^(n) = x(x+1)...(x+n-1)
        as a quotient of factorials

    Examples
    ========
        >>> from sympy.modules.specfun.factorials import *
        >>> rising_factorial(3, 2)
        12

    """

    def eval(self):
        x, n = self._args
        return factorial_simplify(_fac(x+n-1) / _fac(x-1))

    def __latex__(self):
        x, n = self._args
        return "{(%s)}^{(%s)}" % (x.__latex__(), n.__latex__())


class falling_factorial(Function2):
    """
    Usage
    =====
        Calculate the falling factorial (x)_(n) = x(x-1)...(x-n+1)
        as a quotient of factorials

    Examples
    ========
        >>> from sympy.modules.specfun.factorials import *
        >>> falling_factorial(5, 3)
        60

    """
    def eval(self):
        x, n = self._args
        return factorial_simplify(_fac(x) / _fac(x-n))

    def __latex__(self):
        x, n = self._args
        return "{(%s)}_{(%s)}" % (x.__latex__(), n.__latex__())


class binomial(Function2):
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
        >>> from sympy.modules.specfun.factorials import *
        >>> binomial(15,8)
        6435
        >>> # Building Pascal's triangle
        >>> [binomial(0,k) for k in range(1)]
        [1]
        >>> [binomial(1,k) for k in range(2)]
        [1, 1]
        >>> [binomial(2,k) for k in range(3)]
        [1, 2, 1]
        >>> [binomial(3,k) for k in range(4)]
        [1, 3, 3, 1]
        >>> # n can be arbitrary if k is a positive integer
        >>> binomial(Rational(5,4), 3)
        -5/128
        >>> x = Symbol('x')
        >>> binomial(x, 3)
        1/6*x*(-2+x)*(-1+x)

    """
    def eval(self):
        n, k = self._args

        # TODO: move these two cases to factorial_simplify as well
        if n == 0 and k != 0:
            return sin(pi*k)/(pi*k)

        return factorial_simplify(_fac(n) / _fac(k) / _fac(n-k))

    def __latex__(self):
        n, k = self._args
        return r"{{%s}\choose{%s}}" % (n.__latex__(), k.__latex__())

class gamma(DefinedFunction):
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
        >>> from sympy.modules.specfun.factorials import *
        >>> gamma(3)
        2
        >>> gamma(Rational(1,2))
        pi**(1/2)

    """

    def eval(self):
        x = self._args
        y = factorial(x-1)
        if isinstance(y, factorial):
            return self
        else:
            return y

    def diff(self, sym):
        from zeta_functions import polygamma
        x = self._args
        return gamma(x)*polygamma(0,x)*x.diff(sym)

    def __latex__(self):
        return "\Gamma(" + self._args.__latex__() + ")"


class lower_gamma(Function2):
    """
    Lower incomplete gamma function
    
    gamma(a, x)
    """
    def eval(self):
        a, x = self._args
        if a == 1:
            return 1 - exp(-x)
        if a.is_integer and a > 1:
            b = a-1
            return b*lower_gamma(b, x) - x**b * exp(-x)
        return self


class upper_gamma(Function2):
    """
    Upper incomplete gamma function

    Gamma(a, x)
    """
    def eval(self):
        a, x = self._args
        if x == 0:
            return gamma(a)
        if a == 1:
            return exp(-x)
        if a.is_integer and a > 1:
            b = a-1
            return b*upper_gamma(b, x) + x**b * exp(-x)
        return self
