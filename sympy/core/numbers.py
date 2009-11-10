from basic import Atom, SingletonMeta, S, Basic
from decorators import _sympifyit
from cache import Memoizer, cacheit, clear_cache
import sympy.mpmath as mpmath
import sympy.mpmath.libmpf as mlib
import sympy.mpmath.libmpc as mlibc
from sympy.mpmath.libelefun import mpf_pow, mpf_pi, mpf_e, phi_fixed
import decimal

rnd = mlib.round_nearest


# TODO: we should use the warnings module
_errdict = {"divide": False}
def seterr(divide=False):
    """
    Should sympy raise an exception on 0/0 or return a nan?

    divide == True .... raise an exception
    divide == False ... return nan
    """
    if _errdict["divide"] != divide:
        clear_cache()
        _errdict["divide"] = divide

# (a,b) -> gcd(a,b)
_gcdcache = {}

# TODO caching with decorator, but not to degrade performance
def igcd(a, b):
    """Computes integer greatest common divisor of two numbers.

       The algorithm is based on the well known Euclid's algorithm. To
       improve speed, igcd() has its own caching mechanism implemented.
    """
    try:
        return _gcdcache[(a,b)]
    except KeyError:
        if a and b:
            if b < 0:
                b = -b

            while b:
                a, b = b, a % b
        else:
            a = abs(a or b)

        _gcdcache[(a,b)] = a
        return a

def ilcm(a, b):
    """Computes integer least common multiple of two numbers. """
    if a == 0 and b == 0:
        return 0
    else:
        return a * b // igcd(a, b)

def igcdex(a, b):
    """Returns x, y, g such that g = x*a + y*b = gcd(a, b).

       >>> from sympy.core.numbers import igcdex
       >>> igcdex(2, 3)
       (-1, 1, 1)
       >>> igcdex(10, 12)
       (-1, 1, 2)

       >>> x, y, g = igcdex(100, 2004)
       >>> x, y, g
       (-20, 1, 4)
       >>> x*100 + y*2004
       4

    """
    if (not a) and (not b):
        return (0, 1, 0)

    if not a:
        return (0, b//abs(b), abs(b))
    if not b:
        return (a//abs(a), 0, abs(a))

    if a < 0:
        a, x_sign = -a, -1
    else:
        x_sign = 1

    if b < 0:
        b, y_sign = -b, -1
    else:
        y_sign = 1

    x, y, r, s = 1, 0, 0, 1

    while b:
        (c, q) = (a % b, a // b)
        (a, b, r, s, x, y) = (b, c, x-q*r, y-q*s, r, s)

    return (x*x_sign, y*y_sign, a)

@Memoizer((int, long), return_value_converter = lambda d: d.copy())
def factor_trial_division(n):
    """
    Factor any integer into a product of primes, 0, 1, and -1.
    Returns a dictionary {<prime: exponent>}.
    """
    if not n:
        return {0:1}
    factors = {}
    if n < 0:
        factors[-1] = 1
        n = -n
    if n==1:
        factors[1] = 1
        return factors
    d = 2
    while n % d == 0:
        try:
            factors[d] += 1
        except KeyError:
            factors[d] = 1
        n //= d
    d = 3
    while n > 1 and d*d <= n:
        if n % d:
            d += 2
        else:
            try:
                factors[d] += 1
            except KeyError:
                factors[d] = 1
            n //= d
    if n>1:
        try:
            factors[n] += 1
        except KeyError:
            factors[n] = 1
    return factors


class Number(Atom):
    """
    Represents any kind of number in sympy.

    Floating point numbers are represented by the Real class.
    Integer numbers (of any size), together with rational numbers (again, there
    is no limit on their size) are represented by the Rational class.

    If you want to represent for example 1+sqrt(2), then you need to do:

    Rational(1) + sqrt(Rational(2))
    """
    is_commutative = True
    is_comparable = True
    is_bounded = True
    is_finite = True

    __slots__ = []

    # Used to make max(x._prec, y._prec) return x._prec when only x is a float
    _prec = -1

    is_Number = True

    def __new__(cls, *obj):
        if len(obj)==1:
          obj=obj[0]
        if isinstance(obj, (int, long)):
            return Integer(obj)
        if isinstance(obj,tuple) and len(obj)==2:
            return Rational(*obj)
        if isinstance(obj, (str,float,mpmath.mpf,decimal.Decimal)):
            return Real(obj)
        if isinstance(obj, Number):
            return obj
        raise TypeError("expected str|int|long|float|Decimal|Number object but got %r" % (obj))

    def _as_mpf_val(self, prec):
        """Evaluate to mpf tuple accurate to at least prec bits"""
        raise NotImplementedError('%s needs ._as_mpf_val() method' % \
            (self.__class__.__name__))

    def _eval_evalf(self, prec):
        return Real._new(self._as_mpf_val(prec), prec)

    def _as_mpf_op(self, prec):
        prec = max(prec, self._prec)
        return self._as_mpf_val(prec), prec

    def __float__(self):
        return mlib.to_float(self._as_mpf_val(53))

    def _eval_derivative(self, s):
        return S.Zero

    def _eval_conjugate(self):
        return self

    def _eval_order(self, *symbols):
        # Order(5, x, y) -> Order(1,x,y)
        return C.Order(S.One, *symbols)

    def __eq__(self, other):
        raise NotImplementedError('%s needs .__eq__() method' % (self.__class__.__name__))
    def __ne__(self, other):
        raise NotImplementedError('%s needs .__ne__() method' % (self.__class__.__name__))
    def __lt__(self, other):
        raise NotImplementedError('%s needs .__lt__() method' % (self.__class__.__name__))
    def __le__(self, other):
        raise NotImplementedError('%s needs .__le__() method' % (self.__class__.__name__))

    def __gt__(self, other):
        return _sympify(other).__lt__(self)
    def __ge__(self, other):
        return _sympify(other).__le__(self)

    def as_coeff_terms(self, x=None):
        # a -> c * t
        return self, tuple()


class Real(Number):
    """
    Represents a floating point number. It is capable of representing
    arbitrary-precision floating-point numbers

    Usage:
    ======
        Real(3.5)   .... 3.5 (the 3.5 was converted from a python float)
        Real("3.0000000000000005")

        Real((1,3,0,2)) # mpmath tuple: (-1)**1 * 3 * 2**0; 3 has 2 bits
        -3.00000000000000

    Notes:
    ======
        - Real(x) with x being a Python int/long will return Integer(x)
    """
    is_real = True
    is_irrational = False
    is_integer = False

    __slots__ = ['_mpf_', '_prec']

    # mpz can't be pickled
    def __getnewargs__(self):
        return (mlib.to_pickable(self._mpf_),)

    def __getstate__(self):
        d = Basic.__getstate__(self).copy()
        del d["_mpf_"]
        return mlib.to_pickable(self._mpf_), d

    def __setstate__(self, state):
        _mpf_, d = state
        _mpf_ = mlib.from_pickable(_mpf_)
        self._mpf_ = _mpf_
        Basic.__setstate__(self, d)

    is_Real = True

    def floor(self):
        return C.Integer(int(mlib.to_int(mlib.mpf_floor(self._mpf_, self._prec))))

    def ceiling(self):
        return C.Integer(int(mlib.to_int(mlib.mpf_ceil(self._mpf_, self._prec))))

    @property
    def num(self):
        return mpmath.mpf(self._mpf_)

    def _as_mpf_val(self, prec):
        return self._mpf_

    def _as_mpf_op(self, prec):
        return self._mpf_, max(prec, self._prec)

    def __new__(cls, num, prec=15):
        prec = mpmath.settings.dps_to_prec(prec)
        if isinstance(num, (int, long)):
            return Integer(num)
        if isinstance(num, (str, decimal.Decimal)):
            _mpf_ = mlib.from_str(str(num), prec, rnd)
        elif isinstance(num, tuple) and len(num) == 4:
            if type(num[1]) is str:
                # it's a hexadecimal (coming from a pickled object)
                # assume that it is in standard form
                num = list(num)
                num[1] = long(num[1], 16)
                _mpf_ = tuple(num)
            else:
                _mpf_ = mpmath.mpf(
                    S.NegativeOne ** num[0] * num[1] * 2 ** num[2])._mpf_
        else:
            _mpf_ = mpmath.mpf(num)._mpf_
        if not num:
            return C.Zero()
        obj = Basic.__new__(cls)
        obj._mpf_ = _mpf_
        obj._prec = prec
        return obj

    @classmethod
    def _new(cls, _mpf_, _prec):
        if _mpf_ == mlib.fzero:
            return S.Zero
        obj = Basic.__new__(cls)
        obj._mpf_ = _mpf_
        obj._prec = _prec
        return obj

    def _hashable_content(self):
        return (self._mpf_, self._prec)

    def _eval_is_positive(self):
        return self.num > 0

    def _eval_is_negative(self):
        return self.num < 0

    def __neg__(self):
        return Real._new(mlib.mpf_neg(self._mpf_), self._prec)

    @_sympifyit('other', NotImplemented)
    def __mul__(self, other):
        if isinstance(other, Number):
            rhs, prec = other._as_mpf_op(self._prec)
            return Real._new(mlib.mpf_mul(self._mpf_, rhs, prec, rnd), prec)
        return Number.__mul__(self, other)

    @_sympifyit('other', NotImplemented)
    def __mod__(self, other):
        if isinstance(other, Number):
            rhs, prec = other._as_mpf_op(self._prec)
            return Real._new(mlib.mpf_mod(self._mpf_, rhs, prec, rnd), prec)
        return Number.__mod__(self, other)

    @_sympifyit('other', NotImplemented)
    def __add__(self, other):
        if (other is S.NaN) or (self is NaN):
            return S.NaN
        if isinstance(other, Number):
            rhs, prec = other._as_mpf_op(self._prec)
            return Real._new(mlib.mpf_add(self._mpf_, rhs, prec, rnd), prec)
        return Number.__add__(self, other)

    def _eval_power(b, e):
        """
        b is Real but not equal to rationals, integers, 0.5, oo, -oo, nan
        e is symbolic object but not equal to 0, 1

        (-p) ** r -> exp(r * log(-p)) -> exp(r * (log(p) + I*Pi)) ->
                  -> p ** r * (sin(Pi*r) + cos(Pi*r) * I)
        """
        if isinstance(e, Number):
            if isinstance(e, Integer):
                prec = b._prec
                return Real._new(mlib.mpf_pow_int(b._mpf_, e.p, prec, rnd), prec)
            e, prec = e._as_mpf_op(b._prec)
            b = b._mpf_
            try:
                y = mpf_pow(b, e, prec, rnd)
                return Real._new(y, prec)
            except mlib.ComplexResult:
                re, im = mlibc.mpc_pow((b, mlib.fzero), (e, mlib.fzero), prec, rnd)
                return Real._new(re, prec) + Real._new(im, prec) * S.ImaginaryUnit

    def __abs__(self):
        return Real._new(mlib.mpf_abs(self._mpf_), self._prec)

    def __int__(self):
        return int(mlib.to_int(self._mpf_))

    def __eq__(self, other):
        try:
            other = _sympify(other)
        except SympifyError:
            return False    # sympy != other  -->  not ==
        if isinstance(other, NumberSymbol):
            if other.is_irrational: return False
            return other.__eq__(self)
        if isinstance(other, FunctionClass): #cos as opposed to cos(x)
            return False
        if other.is_comparable: other = other.evalf()
        if isinstance(other, Number):
            return bool(mlib.mpf_eq(self._mpf_, other._as_mpf_val(self._prec)))
        return False    # Real != non-Number

    def __ne__(self, other):
        try:
            other = _sympify(other)
        except SympifyError:
            return True     # sympy != other
        if isinstance(other, NumberSymbol):
            if other.is_irrational: return True
            return other.__ne__(self)
        if isinstance(other, FunctionClass): #cos as opposed to cos(x)
            return True
        if other.is_comparable: other = other.evalf()
        if isinstance(other, Number):
            return bool(not mlib.mpf_eq(self._mpf_, other._as_mpf_val(self._prec)))
        return True     # Real != non-Number

    def __lt__(self, other):
        try:
            other = _sympify(other)
        except SympifyError:
            return False    # sympy > other
        if isinstance(other, NumberSymbol):
            return other.__ge__(self)
        if other.is_comparable: other = other.evalf()
        if isinstance(other, Number):
            return bool(mlib.mpf_lt(self._mpf_, other._as_mpf_val(self._prec)))
        return Basic.__lt__(self, other)

    def __le__(self, other):
        try:
            other = _sympify(other)
        except SympifyError:
            return False    # sympy > other  -->  ! <=
        if isinstance(other, NumberSymbol):
            return other.__gt__(self)
        if other.is_comparable: other = other.evalf()
        if isinstance(other, Number):
            return bool(mlib.mpf_le(self._mpf_, other._as_mpf_val(self._prec)))
        return Basic.__le__(self, other)

    def epsilon_eq(self, other, epsilon="10e-16"):
        return abs(self - other) < Real(epsilon)

    def _sage_(self):
        import sage.all as sage
        return sage.RealNumber(str(self))


# this is here to work nicely in Sage
RealNumber = Real


def _parse_rational(s):
    """Parse rational number from string representation"""
    # Simple fraction
    if "/" in s:
        p, q = s.split("/")
        return int(p), int(q)
    # Recurring decimal
    elif "[" in s:
        sign = 1
        if s[0] == "-":
            sign = -1
            s = s[1:]
        s, periodic = s.split("[")
        periodic = periodic.rstrip("]")
        offset = len(s) - s.index(".") - 1
        n1 = int(periodic)
        n2 = int("9" * len(periodic))
        r = Rational(*_parse_rational(s)) + Rational(n1, n2*10**offset)
        return sign*r.p, r.q
    # Ordinary decimal string. Use the Decimal class's built-in parser
    else:
        sign, digits, expt = decimal.Decimal(s).as_tuple()
        p = (1, -1)[sign] * int("".join(str(x) for x in digits))
        if expt >= 0:
            return p*(10**expt), 1
        else:
            return p, 10**-expt

class Rational(Number):
    """Represents integers and rational numbers (p/q) of any size.

    Examples
    ========
        >>> from sympy import Rational
        >>> Rational(3)
        3
        >>> Rational(1,2)
        1/2

    You can create a rational from a string:
        >>> Rational("3/5")
        3/5
        >>> Rational("1.23")
        123/100

    Use square brackets to indicate a recurring decimal:
        >>> Rational("0.[333]")
        1/3
        >>> Rational("1.2[05]")
        1193/990
        >>> float(Rational(1193,990))
        1.20505050505


    Low-level
    ---------

    Access nominator and denominator as .p and .q:
        >>> r = Rational(3,4)
        >>> r
        3/4
        >>> r.p
        3
        >>> r.q
        4

    """
    is_real = True
    is_integer = False
    is_rational = True

    __slots__ = ['p', 'q']

    is_Rational = True

    @cacheit
    def __new__(cls, p, q = None):
        if q is None:
            if isinstance(p, basestring):
                p, q = _parse_rational(p)
            elif isinstance(p, Rational):
                return p
            else:
                return Integer(p)
        q = int(q)
        p = int(p)
        if q==0:
            if p==0:
                if _errdict["divide"]:
                    raise ValueError("Indeterminate 0/0")
                else:
                    return S.NaN
            if p<0: return S.NegativeInfinity
            return S.Infinity
        if q<0:
            q = -q
            p = -p
        n = igcd(abs(p), q)
        if n>1:
            p //= n
            q //= n
        if q==1: return Integer(p)
        if p==1 and q==2: return S.Half
        obj = Basic.__new__(cls)
        obj.p = p
        obj.q = q
        #obj._args = (p, q)
        return obj

    def __getnewargs__(self):
        return (self.p, self.q)

    def _hashable_content(self):
        return (self.p, self.q)

    def _eval_is_positive(self):
        return self.p > 0

    def _eval_is_zero(self):
        return self.p == 0

    def __neg__(self): return Rational(-self.p, self.q)

    @_sympifyit('other', NotImplemented)
    def __mul__(self, other):
        if (other is S.NaN) or (self is S.NaN):
            return S.NaN
        if isinstance(other, Real):
            return other * self
        if isinstance(other, Rational):
            return Rational(self.p * other.p, self.q * other.q)
        return Number.__mul__(self, other)

    @_sympifyit('other', NotImplemented)
    def __mod__(self, other):
        if isinstance(other, Rational):
            n = (self.p*other.q) // (other.p*self.q)
            return Rational(self.p*other.q - n*other.p*self.q, self.q*other.q)
        if isinstance(other, Real):
            return self.evalf() % other
        return Number.__mod__(self, other)

    # TODO reorder
    @_sympifyit('other', NotImplemented)
    def __add__(self, other):
        if (other is S.NaN) or (self is S.NaN):
            return S.NaN
        if isinstance(other, Real):
            return other + self
        if isinstance(other, Rational):
            if self.is_unbounded:
                if other.is_bounded:
                    return self
                elif self==other:
                    return self
            else:
                if other.is_unbounded:
                    return other
            return Rational(self.p * other.q + self.q * other.p, self.q * other.q)
        return Number.__add__(self, other)

    def _eval_power(b, e):
        if (e is S.NaN): return S.NaN
        if isinstance(e, Number):
            if isinstance(e, Real):
                return b._eval_evalf(e._prec) ** e
            if e.is_negative:
                # (3/4)**-2 -> (4/3)**2
                ne = -e
                if (ne is S.One):
                    return Rational(b.q, b.p)
                if b < 0:
                    if e.q != 1:
                        return -(S.NegativeOne) ** ((e.p % e.q) / S(e.q)) * Rational(b.q, -b.p) ** ne
                    else:
                        return S.NegativeOne ** ne * Rational(b.q, -b.p) ** ne
                else:
                    return Rational(b.q, b.p) ** ne
            if (e is S.Infinity):
                if b.p > b.q:
                    # (3/2)**oo -> oo
                    return S.Infinity
                if b.p < -b.q:
                    # (-3/2)**oo -> oo + I*oo
                    return S.Infinity + S.Infinity * S.ImaginaryUnit
                return S.Zero
            if isinstance(e, Integer):
                # (4/3)**2 -> 4**2 / 3**2
                return Rational(b.p ** e.p, b.q ** e.p)
            if isinstance(e, Rational):
                if b.p != 1:
                    # (4/3)**(5/6) -> 4**(5/6) * 3**(-5/6)
                    return Integer(b.p) ** e * Integer(b.q) ** (-e)
                if b >= 0:
                    return Integer(b.q)**Rational(e.p * (e.q-1), e.q) / ( Integer(b.q) ** Integer(e.p))
                else:
                    return (-1)**e * (-b)**e

        c,t = b.as_coeff_terms()
        if e.is_even and isinstance(c, Number) and c < 0:
            return (-c * Mul(*t)) ** e

        return

    def _as_mpf_val(self, prec):
        return mlib.from_rational(self.p, self.q, prec, rnd)

    def _mpmath_(self, prec, rnd):
        return mpmath.make_mpf(mlib.from_rational(self.p, self.q, prec, rnd))

    def __abs__(self):
        return Rational(abs(self.p), self.q)

    def __int__(self):
        return int(self.p//self.q)

    def __eq__(self, other):
        try:
            other = _sympify(other)
        except SympifyError:
            return False    # sympy != other  -->  not ==
        if isinstance(other, NumberSymbol):
            if other.is_irrational: return False
            return other.__eq__(self)
        if isinstance(other, FunctionClass): #cos as opposed to cos(x)
            return False
        if other.is_comparable and not isinstance(other, Rational): other = other.evalf()
        if isinstance(other, Number):
            if isinstance(other, Real):
                return bool(mlib.mpf_eq(self._as_mpf_val(other._prec), other._mpf_))
            return bool(self.p==other.p and self.q==other.q)

        return False    # Rational != non-Number

    def __ne__(self, other):
        try:
            other = _sympify(other)
        except SympifyError:
            return True     # sympy != other
        if isinstance(other, NumberSymbol):
            if other.is_irrational: return True
            return other.__ne__(self)
        if isinstance(other, FunctionClass): #cos as opposed to cos(x)
            return True
        if other.is_comparable and not isinstance(other, Rational): other = other.evalf()
        if isinstance(other, Number):
            if isinstance(other, Real):
                return bool(not mlib.mpf_eq(self._as_mpf_val(other._prec), other._mpf_))
            return bool(self.p!=other.p or self.q!=other.q)

        return True     # Rational != non-Number

    def __lt__(self, other):
        try:
            other = _sympify(other)
        except SympifyError:
            return False    # sympy > other  --> not <
        if isinstance(other, NumberSymbol):
            return other.__ge__(self)
        if other.is_comparable and not isinstance(other, Rational): other = other.evalf()
        if isinstance(other, Number):
            if isinstance(other, Real):
                return bool(mlib.mpf_lt(self._as_mpf_val(other._prec), other._mpf_))
            return bool(self.p * other.q < self.q * other.p)
        return Basic.__lt__(self, other)

    def __le__(self, other):
        try:
            other = _sympify(other)
        except SympifyError:
            return False    # sympy > other  -->  not <=
        if isinstance(other, NumberSymbol):
            return other.__gt__(self)
        if other.is_comparable and not isinstance(other, Rational): other = other.evalf()
        if isinstance(other, Number):
            if isinstance(other, Real):
                return bool(mlib.mpf_le(self._as_mpf_val(other._prec), other._mpf_))
            return bool(self.p * other.q <= self.q * other.p)
        return Basic.__le__(self, other)

    def factors(self):
        f = factor_trial_division(self.p).copy()
        for p,e in factor_trial_division(self.q).items():
            try: f[p] += -e
            except KeyError: f[p] = -e

        if len(f)>1 and 1 in f: del f[1]
        return f

    def as_numer_denom(self):
        return Integer(self.p), Integer(self.q)

    def _sage_(self):
        import sage.all as sage
        #XXX: fixme, this should work:
        #return sage.Integer(self[0])/sage.Integer(self[1])
        return sage.Integer(self.p)/sage.Integer(self.q)

# int -> Integer
_intcache = {}


# TODO move this tracing facility to  sympy/core/trace.py  ?
def _intcache_printinfo():
    ints = sorted(_intcache.keys())
    nhit = _intcache_hits
    nmiss= _intcache_misses

    if nhit == 0 and nmiss == 0:
        print
        print 'Integer cache statistic was not collected'
        return

    miss_ratio = float(nmiss) / (nhit+nmiss)

    print
    print 'Integer cache statistic'
    print '-----------------------'
    print
    print '#items: %i' % len(ints)
    print
    print ' #hit   #miss               #total'
    print
    print '%5i   %5i (%7.5f %%)   %5i'    % (nhit, nmiss, miss_ratio*100, nhit+nmiss)
    print
    print ints

_intcache_hits   = 0
_intcache_misses = 0

def int_trace(f):
    import os
    if os.getenv('SYMPY_TRACE_INT', 'no').lower() != 'yes':
        return f

    def Integer_tracer(cls, i):
        global _intcache_hits, _intcache_misses

        try:
            _intcache_hits += 1
            return _intcache[i]
        except KeyError:
            _intcache_hits   -= 1
            _intcache_misses += 1

            return f(cls, i)


    # also we want to hook our _intcache_printinfo into sys.atexit
    import atexit
    atexit.register(_intcache_printinfo)

    return Integer_tracer




class Integer(Rational):

    q = 1
    is_integer = True

    is_Integer = True

    __slots__ = ['p']

    def _as_mpf_val(self, prec):
        return mlib.from_int(self.p)

    def _mpmath_(self, prec, rnd):
        return mpmath.make_mpf(self._as_mpf_val(prec))

    # TODO caching with decorator, but not to degrade performance
    @int_trace
    def __new__(cls, i):
        try:
            return _intcache[i]
        except KeyError:
            # We only work with well-behaved integer types. This converts, for
            # example, numpy.int32 instances.
            ival = int(i)
            if ival == 0: obj = S.Zero
            elif ival == 1: obj = S.One
            elif ival == -1: obj = S.NegativeOne
            else:
                obj = Basic.__new__(cls)
                obj.p = ival

            _intcache[i] = obj
            return obj

    def __getnewargs__(self):
        return (self.p,)

    # Arithmetic operations are here for efficiency
    def __int__(self):
        return self.p

    def __neg__(self):
        return Integer(-self.p)

    def __abs__(self):
        if self.p >= 0:
            return self
        else:
            return Integer(-self.p)

    def __mod__(self, other):
        return self.p % other

    def __rmod__(self, other):
        return other % self.p

    # TODO make it decorator + bytecodehacks?
    def __add__(a, b):
        if type(b) is int:
            return Integer(a.p + b)
        elif isinstance(b, Integer):
            return Integer(a.p + b.p)
        return Rational.__add__(a, b)   # a,b -not- b,a

    def __radd__(a, b):
        if type(b) is int:
            return Integer(b + a.p)
        elif isinstance(b, Integer):
            return Integer(b.p + a.p)
        return Rational.__add__(a, b)

    def __sub__(a, b):
        if type(b) is int:
            return Integer(a.p - b)
        elif isinstance(b, Integer):
            return Integer(a.p - b.p)
        return Rational.__sub__(a, b)

    def __rsub__(a, b):
        if type(b) is int:
            return Integer(b - a.p)
        elif isinstance(b, Integer):
            return Integer(b.p - a.p)
        return Rational.__rsub__(a, b)

    def __mul__(a, b):
        if type(b) is int:
            return Integer(a.p * b)
        elif isinstance(b, Integer):
            return Integer(a.p * b.p)
        return Rational.__mul__(a, b)

    def __rmul__(a, b):
        if type(b) is int:
            return Integer(b * a.p)
        elif isinstance(b, Integer):
            return Integer(b.p * a.p)
        return Rational.__mul__(a, b)

    # XXX __pow__ ?

    # XXX do we need to define __cmp__ ?
#   def __cmp__(a, b):

    def __eq__(a, b):
        if type(b) is int:
            return (a.p == b)
        elif isinstance(b, Integer):
            return (a.p == b.p)
        return Rational.__eq__(a, b)

    def __ne__(a, b):
        if type(b) is int:
            return (a.p != b)
        elif isinstance(b, Integer):
            return (a.p != b.p)
        return Rational.__ne__(a, b)

    def __gt__(a, b):
        if type(b) is int:
            return (a.p >  b)
        elif isinstance(b, Integer):
            return (a.p >  b.p)
        return Rational.__gt__(a, b)

    def __lt__(a, b):
        if type(b) is int:
            return (a.p <  b)
        elif isinstance(b, Integer):
            return (a.p <  b.p)
        return Rational.__lt__(a, b)

    def __ge__(a, b):
        if type(b) is int:
            return (a.p >= b)
        elif isinstance(b, Integer):
            return (a.p >= b.p)
        return Rational.__ge__(a, b)

    def __le__(a, b):
        if type(b) is int:
            return (a.p <= b)
        elif isinstance(b, Integer):
            return (a.p <= b.p)
        return Rational.__le__(a, b)

    ########################################

    def _eval_is_odd(self):
        return bool(self.p % 2)

    def _eval_power(b, e):
        """
        Tries to do some simplifications on b ** e, where b is
        an instance of Integer

        Returns None if no further simplifications can be done

        When exponent is a fraction (so we have for example a square root),
        we try to find the simplest possible representation, so that
          - 4**Rational(1,2) becomes 2
          - (-4)**Rational(1,2) becomes 2*I
        We will
        """
        if e is S.NaN: return S.NaN
        if b is S.One: return S.One
        if b is S.NegativeOne: return
        if e is S.Infinity:
            if b.p > S.One: return S.Infinity
            if b.p == -1: return S.NaN
            # cases 0, 1 are done in their respective classes
            return S.Infinity + S.ImaginaryUnit * S.Infinity
        if not isinstance(e, Number):
            # simplify when exp is even
            # (-2) ** k --> 2 ** k
            c,t = b.as_coeff_terms()
            if e.is_even and isinstance(c, Number) and c < 0:
                return (-c * Mul(*t)) ** e
        if not isinstance(e, Rational): return
        if e is S.Half and b < 0:
            # we extract I for this special case since everyone is doing so
            return S.ImaginaryUnit * Pow(-b, e)
        if e < 0:
            # invert base and change sign on exponent
            ne = -e
            if b < 0:
                if e.q != 1:
                    return -(S.NegativeOne) ** ((e.p % e.q) / S(e.q)) * Rational(1, -b) ** ne
                else:
                    return (S.NegativeOne) ** ne * Rational(1, -b) ** ne
            else:
                return Rational(1, b.p) ** ne
        # see if base is a perfect root, sqrt(4) --> 2
        x, xexact = integer_nthroot(abs(b.p), e.q)
        if xexact:
            # if it's a perfect root we've finished
            result = Integer(x ** abs(e.p))
            if b < 0: result *= (-1)**e
            return result
        # The following is an algorithm where we collect perfect roots
        # from the factors of base
        if b > 4294967296:
            # Prevent from factorizing too big integers
            return None
        dict = b.factors()
        out_int = 1
        sqr_int = 1
        sqr_gcd = 0
        sqr_dict = {}
        for prime,exponent in dict.iteritems():
            exponent *= e.p
            div_e = exponent // e.q
            div_m = exponent % e.q
            if div_e > 0:
                out_int *= prime**div_e
            if div_m > 0:
                sqr_dict[prime] = div_m
        for p,ex in sqr_dict.iteritems():
            if sqr_gcd == 0:
                sqr_gcd = ex
            else:
                sqr_gcd = igcd(sqr_gcd, ex)
        for k,v in sqr_dict.iteritems():
            sqr_int *= k**(v // sqr_gcd)
        if sqr_int == b.p and out_int == 1:
            result = None
        else:
            result = out_int * Pow(sqr_int , Rational(sqr_gcd, e.q))
        return result

    def _eval_is_prime(self):
        if self.p < 0:
            return False

    def as_numer_denom(self):
        return self, S.One

    def __floordiv__(self, other):
        return Integer(self.p // Integer(other).p)

    def __rfloordiv__(self, other):
        return Integer(Integer(other).p // self.p)

class Zero(Integer):
    __metaclass__ = SingletonMeta

    p = 0
    q = 1
    is_positive = False
    is_negative = False
    is_finite = False
    is_zero = True
    is_prime = False
    is_composite = False

    __slots__ = []

    @staticmethod
    def __abs__():
        return S.Zero

    @staticmethod
    def __neg__():
        return S.Zero

    def _eval_power(b, e):
        if e.is_negative:
            return S.Infinity
        if e.is_positive:
            return b
        d = e.evalf()
        if isinstance(d, Number):
            if d.is_negative:
                return S.Infinity
            return b
        coeff, terms = e.as_coeff_terms()
        if coeff.is_negative:
            return S.Infinity ** Mul(*terms)
        if coeff is not S.One:
            return b ** Mul(*terms)

    def _eval_order(self, *symbols):
        # Order(0,x) -> 0
        return self

class One(Integer):
    __metaclass__ = SingletonMeta

    p = 1
    q = 1

    is_prime = True

    __slots__ = []

    def _eval_evalf(self, prec):
        return self

    @staticmethod
    def __abs__():
        return S.One

    @staticmethod
    def __neg__():
        return S.NegativeOne

    def _eval_order(self, *symbols):
        return

    @staticmethod
    def factors():
        return {1: 1}

class NegativeOne(Integer):
    __metaclass__ = SingletonMeta

    p = -1
    q = 1

    __slots__ = []

    def _eval_evalf(self, prec):
        return self

    @staticmethod
    def __abs__():
        return S.One

    @staticmethod
    def __neg__():
        return S.One

    def _eval_power(b, e):
        if e.is_odd: return S.NegativeOne
        if e.is_even: return S.One
        if isinstance(e, Number):
            if isinstance(e, Real):
                return Real(-1.0) ** e
            if e is S.NaN:
                return S.NaN
            if e is S.Infinity  or  e is S.NegativeInfinity:
                return S.NaN
            if e is S.Half:
                return S.ImaginaryUnit
            if isinstance(e, Rational):
                if e.q == 2:
                    return S.ImaginaryUnit ** Integer(e.p)
                q = int(e)
                if q:
                    q = Integer(q)
                    return b ** q * b ** (e - q)
        return

class Half(Rational):
    __metaclass__ = SingletonMeta

    p = 1
    q = 2

    __slots__ = []

    @staticmethod
    def __abs__():
        return S.Half


class Infinity(Rational):
    __metaclass__ = SingletonMeta

    p = 1
    q = 0

    __slots__ = []

    is_commutative = True
    is_positive = True
    is_bounded = False
    is_finite   = False
    is_infinitesimal = False
    is_integer  = None
    is_rational = None
    is_odd = None

    @staticmethod
    def __abs__():
        return S.Infinity

    @staticmethod
    def __neg__():
        return S.NegativeInfinity

    def _eval_power(b, e):
        """
        e is symbolic object but not equal to 0, 1

        oo ** nan -> nan
        oo ** (-p) -> 0, p is number, oo
        """
        if e.is_positive:
            return S.Infinity
        if e.is_negative:
            return S.Zero
        if isinstance(e, Number):
            if e is S.NaN:
                return S.NaN
        d = e.evalf()
        if isinstance(d, Number):
            return b ** d
        return

    def _as_mpf_val(self, prec):
        return mlib.finf

    def _sage_(self):
        import sage.all as sage
        return sage.oo

    def __gt__(a, b):
        if b is S.Infinity:
            return False
        return True

    def __lt__(a, b):
        return False

    def __ge__(a, b):
        return True

    def __le__(a, b):
        if b is S.Infinity:
            return True
        return False

    def __mod__(self, other):
        return S.NaN

    __rmod__ = __mod__


class NegativeInfinity(Rational):
    __metaclass__ = SingletonMeta

    p = -1
    q = 0

    __slots__ = []

    is_commutative = True
    is_real = True
    is_positive = False
    is_bounded = False
    is_finite = False
    is_infinitesimal = False
    is_integer  = None
    is_rational = None

    @staticmethod
    def __abs__():
        return S.Infinity

    @staticmethod
    def __neg__():
        return S.Infinity

    def _eval_power(b, e):
        """
        e is symbolic object but not equal to 0, 1

        (-oo) ** nan -> nan
        (-oo) ** oo  -> nan
        (-oo) ** (-oo) -> nan
        (-oo) ** e -> oo, e is positive even integer
        (-oo) ** o -> -oo, o is positive odd integer

        """
        if isinstance(e, Number):
            if (e is S.NaN)  or  (e is S.Infinity)  or  (e is S.NegativeInfinity):
                return S.NaN
            if isinstance(e, Integer):
                if e.is_positive:
                    if e.is_odd:
                        return S.NegativeInfinity
                    return S.Infinity
            return S.NegativeOne**e * S.Infinity ** e
        return

    def _as_mpf_val(self, prec):
        return mlib.fninf

    def __gt__(a, b):
        return False

    def __lt__(a, b):
        if b is S.NegativeInfinity:
            return False
        return True

    def __ge__(a, b):
        if b is S.NegativeInfinity:
            return True
        return False

    def __le__(a, b):
        return True


class NaN(Rational):
    __metaclass__ = SingletonMeta

    p = 0
    q = 0

    is_commutative = True
    is_real = None
    is_rational = None
    is_integer  = None
    is_comparable = False
    is_finite   = None
    is_bounded = None
    #is_unbounded = False
    is_zero     = None
    is_prime    = None
    is_positive = None

    __slots__ = []

    def _as_mpf_val(self, prec):
        return mlib.fnan

    def _eval_power(b, e):
        if e is S.Zero:
            return S.One
        return b

    def _sage_(self):
        import sage.all as sage
        return sage.NaN

class ComplexInfinity(Atom):
    __metaclass__ = SingletonMeta

    is_commutative = True
    is_comparable = None
    is_bounded = False
    is_real = None

    __slots__ = []

    @staticmethod
    def __abs__():
        return S.Infinity

    @staticmethod
    def __neg__():
        return S.ComplexInfinity

    def _eval_power(b, e):
        if e is S.ComplexInfinity:
            return S.NaN

        if isinstance(e, Number):
            if e is S.Zero:
                return S.NaN
            else:
                if e.is_positive:
                    return S.ComplexInfinity
                else:
                    return S.Zero

class NumberSymbol(Atom):
    __metaclass__ = SingletonMeta

    is_commutative = True
    is_comparable = True
    is_bounded = True
    is_finite = True

    __slots__ = []

    is_NumberSymbol = True

    def approximation(self, number_cls):
        """ Return an interval with number_cls endpoints
        that contains the value of NumberSymbol.
        If not implemented, then return None.
        """

    def _eval_evalf(self, prec):
        return Real._new(self._as_mpf_val(prec), prec)

    def _eval_derivative(self, s):
        return S.Zero

    def __eq__(self, other):
        try:
            other = _sympify(other)
        except SympifyError:
            return False    # sympy != other  -->  not ==
        if self is other: return True
        if isinstance(other, Number) and self.is_irrational: return False

        return False    # NumberSymbol != non-(Number|self)

    def __ne__(self, other):
        try:
            other = _sympify(other)
        except SympifyError:
            return True     # sympy != other
        if self is other: return False
        if isinstance(other, Number) and self.is_irrational: return True

        return True     # NumberSymbol != non(Number|self)

    def __lt__(self, other):
        try:
            other = _sympify(other)
        except SympifyError:
            return False    # sympy > other  --> not <
        if self is other: return False
        if isinstance(other, Number):
            approx = self.approximation_interval(other.__class__)
            if approx is not None:
                l,u = approx
                if other < l: return False
                if other > u: return True
            return self.evalf()<other
        if other.is_comparable:
            other = other.evalf()
            return self.evalf()<other
        return Basic.__lt__(self, other)
    def __le__(self, other):
        try:
            other = _sympify(other)
        except SympifyError:
            return False    # sympy > other  --> not <=
        if self is other: return True
        if other.is_comparable: other = other.evalf()
        if isinstance(other, Number):
            return self.evalf()<=other
        return Basic.__le__(self, other)
    def __gt__(self, other):
        return (-self) < (-other)
    def __ge__(self, other):
        return (-self) <= (-other)


class Exp1(NumberSymbol):

    is_real = True
    is_positive = True
    is_negative = False # XXX Forces is_negative/is_nonnegative
    is_irrational = True

    __slots__ = []

    @staticmethod
    def __abs__():
        return S.Exp1

    def _as_mpf_val(self, prec):
        return mpf_e(prec)

    def approximation_interval(self, number_cls):
        if issubclass(number_cls,Integer):
            return (Integer(2),Integer(3))
        elif issubclass(number_cls,Rational):
            pass

    def _eval_power(self, exp):
        return C.exp(exp)

    def _sage_(self):
        import sage.all as sage
        return sage.e

class Pi(NumberSymbol):

    is_real = True
    is_positive = True
    is_negative = False
    is_irrational = True

    __slots__ = []

    @staticmethod
    def __abs__():
        return S.Pi

    def _as_mpf_val(self, prec):
        return mpf_pi(prec)

    def approximation_interval(self, number_cls):
        if issubclass(number_cls, Integer):
            return (Integer(3), Integer(4))
        elif issubclass(number_cls, Rational):
            return (Rational(223,71), Rational(22,7))

    def _sage_(self):
        import sage.all as sage
        return sage.pi

class GoldenRatio(NumberSymbol):

    is_real = True
    is_positive = True
    is_negative = False
    is_irrational = True

    __slots__ = []

    def _as_mpf_val(self, prec):
        return mlib.from_man_exp(phi_fixed(prec+10), -prec-10)

    def _eval_expand_func(self, deep=True, **hints):
        return S.Half + S.Half*S.Sqrt(5)

    def approximation_interval(self, number_cls):
        if issubclass(number_cls, Integer):
            return (S.One, Rational(2))
        elif issubclass(number_cls, Rational):
            pass

    def _sage_(self):
        import sage.all as sage
        return sage.golden_ratio

class EulerGamma(NumberSymbol):

    is_real = True
    is_positive = True
    is_negative = False
    is_irrational = None

    __slots__ = []

    def _as_mpf_val(self, prec):
        return mlib.from_man_exp(mpmath.gammazeta.euler_fixed(prec+10), -prec-10)

    def approximation_interval(self, number_cls):
        if issubclass(number_cls, Integer):
            return (S.Zero, S.One)
        elif issubclass(number_cls, Rational):
            return (S.Half, Rational(3, 5))

    def _sage_(self):
        import sage.all as sage
        return sage.euler_gamma

class Catalan(NumberSymbol):

    is_real = True
    is_positive = True
    is_negative = False
    is_irrational = None

    __slots__ = []

    def _as_mpf_val(self, prec):
        return mlib.from_man_exp(mpmath.gammazeta.catalan_fixed(prec+10), -prec-10)

    def approximation_interval(self, number_cls):
        if issubclass(number_cls, Integer):
            return (S.Zero, S.One)
        elif issubclass(number_cls, Rational):
            return (Rational(9, 10), S.One)

    def _sage_(self):
        import sage.all as sage
        return sage.catalan

class ImaginaryUnit(Atom):
    __metaclass__ = SingletonMeta

    is_commutative = True
    is_imaginary = True
    is_bounded = True
    is_finite = True

    __slots__ = []

    @staticmethod
    def __abs__():
        return S.One

    def _eval_evalf(self, prec):
        return self

    def _eval_conjugate(self):
        return -S.ImaginaryUnit

    def _eval_derivative(self, s):
        return S.Zero

    def _eval_power(b, e):
        """
        b is I = sqrt(-1)
        e is symbolic object but not equal to 0, 1

        I ** r -> (-1)**(r/2) -> exp(r/2 * Pi * I) -> sin(Pi*r/2) + cos(Pi*r/2) * I, r is decimal
        I ** 0 mod 4 -> 1
        I ** 1 mod 4 -> I
        I ** 2 mod 4 -> -1
        I ** 3 mod 4 -> -I
        """


        if isinstance(e, Number):
            #if isinstance(e, Decimal):
            #    a = decimal_math.pi() * exponent.num / 2
            #    return Decimal(decimal_math.sin(a) + decimal_math.cos(a) * ImaginaryUnit())
            if isinstance(e, Integer):
                e = e.p % 4
                if e==0: return S.One
                if e==1: return S.ImaginaryUnit
                if e==2: return -S.One
                return -S.ImaginaryUnit
            return (-S.One) ** (e * S.Half)
        return

    def as_base_exp(self):
        return -S.One, S.Half

    def _sage_(self):
        import sage.all as sage
        return sage.I

_intcache[0] = S.Zero
_intcache[1] = S.One
_intcache[-1]= S.NegativeOne

Basic.singleton['E'] = Exp1
Basic.singleton['pi'] = Pi
Basic.singleton['I'] = ImaginaryUnit
Basic.singleton['oo'] = Infinity
Basic.singleton['nan'] = NaN

Basic.singleton['zoo'] = ComplexInfinity

Basic.singleton['GoldenRatio'] = GoldenRatio
Basic.singleton['EulerGamma'] = EulerGamma
Basic.singleton['Catalan'] = Catalan

from basic import Basic, Atom, S, C, SingletonMeta
from cache import Memoizer
from sympify import _sympify, SympifyError
from function import FunctionClass
from power import Pow, integer_nthroot
from mul import Mul
