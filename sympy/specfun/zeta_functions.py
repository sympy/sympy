from factorials import factorial_, binomial2

from sympy.core import *

# first some utilities for calculating Bernoulli numbers

class _memo:
    def __init__(self, f):
        self.func = f
        self.cache = {}
        self.highest = 0
    def __call__(self, m):
        if m < self.highest:
            return self.cache[m]
        else:
            x = self.func(m)
            self.cache[m] = x
            self.highest = m
            return x

def _product(low, high):
    p = 1
    for k in xrange(low, high+1):
        p *= k
    return p

def _bernoulli_sum(m, mtop):
    s = 0
    a = int(binomial2(m+3, m-6))
    for j in xrange(1, mtop+1):
        s += a*bernoulli2(m-6*j)
        a *= _product(m-6 - 6*j + 1, m-6*j)
        a //= _product(6*j+4, 6*j+9)
    return s

@_memo
def _b0mod6(m):
    return (Rational(m+3, 3) - _bernoulli_sum(m, m//6)) / binomial2(m+3, m)

@_memo
def _b2mod6(m):
    return (Rational(m+3, 3) - _bernoulli_sum(m, (m-2)//6)) / binomial2(m+3, m)

@_memo
def _b4mod6(m):
    return (-Rational(m+3, 6) - _bernoulli_sum(m, (m-4)//6)) / binomial2(m+3, m)

class Bernoulli2(DefinedFunction):
    """
    Usage
    =====
        bernoulli2(n) -> nth Bernoulli number, B_n

    Notes
    =====
        * Bernoulli numbers are rational numbers
        * For odd integers n > 1, bernoulli(n) = 0

    Examples
    ========
        >>> from sympy.specfun.zeta_functions import *
        >>> [bernoulli2(n) for n in range(11)]
        [1, -1/2, 1/6, 0, -1/30, 0, 1/42, 0, -1/30, 0, 5/66]
        >>> bernoulli2(1000001)
        0

    """
    nofargs = 1

    def _eval_apply(self, m):
        if m.is_integer and m >= 0:
            m = int(m)
            if m == 0: return 1
            if m == 1: return -Rational(1,2)
            if m % 6 == 0: return _b0mod6(m)
            if m % 6 == 2: return _b2mod6(m)
            if m % 6 == 4: return _b4mod6(m)
            return 0


# TODO: speed up
class BernoulliPoly(DefinedFunction):
    """
    bernoulli_poly(n, x) - nth Bernoulli polynomial of x
    """
    nofargs = 2

    def _eval_apply(self, n, x):
        if isinstance(n, Rational) and n.is_integer:
            s = 0
            for k in xrange(n.p+1):
                s += binomial2(n,k)*bernoulli2(k)*x**(n-k)
            return s


class Zeta(DefinedFunction):
    """
    Usage
    =====
        zeta(s) -> Riemann zeta function of s

    Notes
    =====
        * The zeta function has a pole at s = 1
        * For even positive integers n, zeta(n) is a rational number
          times a power of pi
        * For nonpositive integers n, zeta(n) is a rational number

    Examples
    ========
        >>> from sympy.specfun.zeta_functions import *
        >>> zeta(4)
        1/90*pi**4
        >>> zeta(5)  # no closed form known
        zeta(5)

    """
    nofargs = 1

    def _eval_apply(self, s):
        if s.is_integer:
            if s == 0:
                return Rational(-1,2)
            if s == 1:
                return oo
            if s > 1 and int(s) % 2 == 0:
                return abs(bernoulli2(s)) * 2**(s-1) / factorial(s) * pi**s
            if s < 1:
                return -bernoulli2(-s+1)/(-s+1)


class DirichletEta(DefinedFunction):
    """
    Dirichlet eta function
    """
    nofargs = 1

    def _eval_apply(self, s):
        if s == 1:
            return log(2)
        else:
            return (1-2**(1-s)) * zeta(s)


class Harmonic(DefinedFunction):
    """
    harmonic(n, m=1) -- nth harmonic number (of order m)
    """
    nofargs = 2

    #def __repr__(self):
    #    return "harmonic(%r, %r)" % self._args

    #__str__ = __repr__

    def _eval_apply(self, n, m):
        if n == oo:
            return zeta(m)
        if isinstance(n, Rational) and n >= 0 and \
           isinstance(m, Rational) and m >= 0:
            print "yeska", n, m
            if n == 0:
                return 0
            s = 0
            for i in xrange(1, n.p+1):
                s += Rational(1)/i**m
            return s


class PolyGamma(DefinedFunction):
    """
    polygamma(m, z) -- mth order polygamma function of z
    """
    nofargs = 2

    #def __repr__(self):
    #    return "polygamma(%r, %r)" % self._args

    #__str__ = __repr__

    def _eval_apply(self, m, z):
        # TODO: rational arguments, reflection formula
        if m.is_integer and m >= 0 and z == 0:
            return oo
        if m == 0:
            if isinstance(z, Rational):
                #if z < 0:
                #    return polygamma(0, 1-z) - pi/tan(pi*z)
                if z.is_integer and z > 0:
                    return -EulerGamma + harmonic(z-1, 1)
        # if m == 1:
        #    if isinstance(z, Rational) and z < 0:
        #        return -polygamma(1, 1-z) + pi**2 / sin(pi*z)**2
        if m.is_integer and m > 0 and z.is_integer and z > 0:
            return (-1)**(m+1)*factorial_(m)*(zeta(m+1)-harmonic(z-1, m+1))
        if m.is_integer and z == Rational(1,2):
            return (-1)**(m+1)*factorial_(m)*(2**(m+1)-1)*zeta(m+1)

    def fdiff(self, argindex=2):
        if argindex == 2:
            m = Basic.Symbol('m', dummy=True)
            z = Basic.Symbol('z', dummy=True)
            return Lambda(polygamma(m+1, z), m, z)
        else:
            raise NotImplementedError


def digamma(z):
    return polygamma(0, z)

def trigamma(z):
    return polygamma(1, z)

def tetragamma(z):
    return polygamma(2, z)


bernoulli2 = Bernoulli2()
bernoulli_poly = BernoulliPoly()
zeta = Zeta()
dirichlet_eta = DirichletEta()
harmonic = Harmonic()
polygamma = PolyGamma()
