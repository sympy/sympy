import math
import decimal
import decimal_math
from basic import Basic, Atom, Singleton, S, Memoizer, MemoizerArg
from methods import RelMeths, ArithMeths
from power import integer_nthroot


@Memoizer((int, long), (int, long))
def gcd(a, b):
    '''Returns the Greatest Common Divisor,
    implementing Euclid\'s algorithm.'''
    while a:
        a, b = b%a, a
    return b

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


class Number(Atom, RelMeths, ArithMeths):
    """Represents any kind of number in sympy.


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

    def __new__(cls, *obj):
        if len(obj)==1: obj=obj[0]
        if isinstance(obj, (int, long)):
            return Integer(obj)
        if isinstance(obj,tuple) and len(obj)==2:
            return Rational(*obj)
        if isinstance(obj, (str,float,decimal.Decimal)):
            return Real(obj)
        if isinstance(obj, Number):
            return obj
        raise TypeError("expected str|int|long|float|Decimal|Number object but got %r" % (obj))

    def _eval_evalf(self):
        r = self._as_decimal()
        return Real(r)

    def __float__(self):
        return float(self._as_decimal())

    def _as_decimal(self):
        raise NotImplementedError('%s needs ._as_decimal() method' % (self.__class__.__name__))

    def _eval_derivative(self, s):
        return Zero()

    def _eval_conjugate(self):
        return self

    def _eval_apply(self, a):
        return self*a

    def _eval_order(self, *symbols):
        # Order(5, x, y) -> Order(1,x,y)
        return Basic.Order(Basic.One(),*symbols)

    def sqrt(self): return Real(decimal_math.sqrt(self._as_decimal()))
    def exp(self): return Real(decimal_math.exp(self._as_decimal()))
    def log(self): return Real(decimal_math.log(self._as_decimal()))
    def sin(self): return Real(decimal_math.sin(self._as_decimal()))
    def cos(self): return Real(decimal_math.cos(self._as_decimal()))
    def tan(self): return Real(decimal_math.tan(self._as_decimal()))
    def cot(self): return Real(decimal_math.cot(self._as_decimal()))
    def asin(self): return Real(decimal_math.asin(self._as_decimal()))
    def acos(self): return Real(decimal_math.acos(self._as_decimal()))
    def atan(self): return Real(decimal_math.atan(self._as_decimal()))
    def acot(self): return Real(decimal_math.acot(self._as_decimal()))
    def sinh(self): return Real(decimal_math.sinh(self._as_decimal()))
    def cosh(self): return Real(decimal_math.cosh(self._as_decimal()))
    def tanh(self): return Real(decimal_math.tanh(self._as_decimal()))
    def coth(self): return Real(decimal_math.coth(self._as_decimal()))
    def asinh(self): return Real(decimal_math.asinh(self._as_decimal()))
    def acosh(self): return Real(decimal_math.acosh(self._as_decimal()))
    def atanh(self): return Real(decimal_math.atanh(self._as_decimal()))
    def acoth(self): return Real(decimal_math.acoth(self._as_decimal()))
    def floor(self): return Real(decimal_math.floor(self._as_decimal()))
    def ceiling(self): return Real(decimal_math.ceiling(self._as_decimal()))

    def __eq__(self, other):
        raise NotImplementedError,'%s needs .__eq__() method' % (self.__class__.__name__)
    def __ne__(self, other):
        raise NotImplementedError,'%s needs .__ne__() method' % (self.__class__.__name__)
    def __lt__(self, other):
        raise NotImplementedError,'%s needs .__lt__() method' % (self.__class__.__name__)
    def __le__(self, other):
        raise NotImplementedError,'%s needs .__le__() method' % (self.__class__.__name__)

    def __gt__(self, other):
        return Basic.sympify(other).__lt__(self)
    def __ge__(self, other):
        return Basic.sympify(other).__le__(self)

    def as_coeff_terms(self, x=None):
        # a -> c * t
        return self, []

decimal_to_Number_cls = {
    decimal.Decimal('0').as_tuple():'Zero',
    decimal.Decimal('-0').as_tuple():'Zero',
    decimal.Decimal('1').as_tuple():'One',
    decimal.Decimal('-1').as_tuple():'NegativeOne',
    decimal.Decimal('Infinity').as_tuple():'Infinity',
    decimal.Decimal('-Infinity').as_tuple():'NegativeInfinity',
    decimal.Decimal('NaN').as_tuple():'NaN',
    }

def convert_to_Decimal(num):
    if isinstance(num, (str, int, long)):
        num = decimal.Decimal(num)
    elif isinstance(num, float):
        num = Real.float_to_decimal(num)
    return num

class Real(Number):
    """Represents a floating point number. It is capable of representing
    arbitrary-precision floating-point numbers

    Usage:

    Real(3.5)   .... 3.5 (the 3.5 was converted from a python float)
    Real("3.0000000000000005")

    """
    is_real = True
    is_irrational = False
    is_integer = False

    @Memoizer(type, MemoizerArg((str, int, long, float, decimal.Decimal), convert_to_Decimal))
    def __new__(cls, num):
        singleton_cls_name = decimal_to_Number_cls.get(num.as_tuple(), None)
        if singleton_cls_name is not None:
            return getattr(Basic, singleton_cls_name)()
        obj = Basic.__new__(cls)
        obj.num = num
        return obj

    @staticmethod
    def float_to_decimal(f):
        "Convert a floating point number to a Decimal with no loss of information"
        # Transform (exactly) a float to a mantissa (0.5 <= abs(m) < 1.0) and an
        # exponent.  Double the mantissa until it is an integer.  Use the integer
        # mantissa and exponent to compute an equivalent Decimal.  If this cannot
        # be done exactly, then retry with more precision.

        try:
            mantissa, exponent = math.frexp(f)
        except OverflowError:
            return decimal.Inf

        while mantissa != int(mantissa):
            mantissa *= 2.0
            exponent -= 1
        mantissa = int(mantissa)

        oldcontext = decimal.getcontext()
        decimal.setcontext(decimal.Context(traps=[decimal.Inexact]))
        try:
            while True:
                try:
                    return mantissa * decimal.Decimal(2) ** exponent
                except decimal.Inexact:
                    decimal.getcontext().prec += 1
        finally:
            decimal.setcontext(oldcontext)

    def _hashable_content(self):
        return (self.num,)

    def tostr(self, level=0):
        r = str(self.num.normalize())
        if self.precedence<=level:
            return '(%s)' % (r)
        return r

    def torepr(self):
        return '%s(%r)' % (self.__class__.__name__, str(self.num))

    def _eval_is_positive(self):
        return self.num.as_tuple()[0] == 0

    def _eval_is_negative(self):
        return self.num.as_tuple()[0] == 1

    def _eval_evalf(self):
        return self

    def _as_decimal(self):
        return self.num

    def __neg__(self):
        return Real(-self.num)

    def __mul__(self, other):
        other = Basic.sympify(other)
        if isinstance(other, Number):
            return Real(self.num * other._as_decimal())
        return Number.__mul__(self, other)

    def __add__(self, other):
        other = Basic.sympify(other)
        if isinstance(other, NaN) or isinstance(self, NaN):
            return NaN()
        if isinstance(other, Number):
            return Real(self.num + other._as_decimal())
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
                e = e.p
                return Real(decimal_math.pow(b.num, e))

            e2 = e._as_decimal()
            if b.is_negative and not e.is_integer:
                m = decimal_math.pow(-b.num, e2)
                a = decimal_math.pi() * e2
                s = m * decimal_math.sin(a)
                c = m * decimal_math.cos(a)
                return Real(s) + Real(c) * ImaginaryUnit()
            return Real(decimal_math.pow(b.num, e2))
        return

    def __abs__(self):
        return Real(abs(self.num))

    def __int__(self):
        return int(self.num)

    def __float__(self):
        return float(self.num)

    def __eq__(self, other):
        other = Basic.sympify(other)
        if isinstance(other, NumberSymbol):
            if other.is_irrational: return False
            return other.__eq__(self)
        if other.is_comparable: other = other.evalf()
        if isinstance(other, Number):
            return bool(self._as_decimal()==other._as_decimal())
        return RelMeths.__eq__(self, other)
    def __ne__(self, other):
        other = Basic.sympify(other)
        if isinstance(other, NumberSymbol):
            if other.is_irrational: return True
            return other.__ne__(self)
        if other.is_comparable: other = other.evalf()
        if isinstance(other, Number):
            return bool(self._as_decimal()!=other._as_decimal())
        return RelMeths.__ne__(self, other)
    def __lt__(self, other):
        other = Basic.sympify(other)
        if isinstance(other, NumberSymbol):
            return other.__ge__(self)
        if other.is_comparable: other = other.evalf()
        if isinstance(other, Number):
            return bool(self._as_decimal() < other._as_decimal())
        return RelMeths.__lt__(self, other)
    def __le__(self, other):
        other = Basic.sympify(other)
        if isinstance(other, NumberSymbol):
            return other.__gt__(self)
        if other.is_comparable: other = other.evalf()
        if isinstance(other, Number):
            return bool(self._as_decimal()<=other._as_decimal())
        return RelMeths.__le__(self, other)

    def epsilon_eq(self, other, epsilon="10e-16"):
        return abs(self - other) < Basic.Real(epsilon)

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
        1.2050505050505051
    """
    is_real = True
    is_integer = False
    is_rational = True

    @Memoizer(type, (int, long, str), MemoizerArg((int, long, type(None)), name="q"))
    def __new__(cls, p, q = None):
        if q is None:
            if isinstance(p, str):
                p, q = _parse_rational(p)
            else:
                return Integer(p)
        if q==0:
            if p==0: return NaN()
            if p<0: return NegativeInfinity()
            return Infinity()
        if q<0:
            q = -q
            p = -p
        n = gcd(abs(p), q)
        if n>1:
            p /= n
            q /= n
        if q==1: return Integer(p)
        if p==1 and q==2: return Half()
        obj = Basic.__new__(cls)
        obj.p = p
        obj.q = q
        return obj

    def _hashable_content(self):
        return (self.p, self.q)

    def tostr(self, level=0):
        if self.precedence<=level:
            return '(%s/%s)' % (self.p, self.q)
        return '%s/%s' % (self.p, self.q)

    def torepr(self):
        return '%s(%r, %r)' % (self.__class__.__name__, self.p, self.q)

    @property
    def precedence(self):
        if self.p < 0:
            return Basic.Add_precedence
        return Basic.Mul_precedence

    def _eval_is_positive(self):
        return self.p > 0

    def _eval_is_zero(self):
        return self.p == 0

    def __neg__(self): return Rational(-self.p, self.q)

    def __mul__(self, other):
        other = Basic.sympify(other)
        if isinstance(other, NaN) or isinstance(self, NaN):
            return NaN()
        if isinstance(other, Real):
            return Real(self._as_decimal() * other.num)
        if isinstance(other, Rational):
            return Rational(self.p * other.p, self.q * other.q)
        return Number.__mul__(self, other)

    def __add__(self, other):
        other = Basic.sympify(other)
        if isinstance(other, NaN) or isinstance(self, NaN):
            return NaN()
        if isinstance(other, Real):
            return Real(self._as_decimal() + other.num)
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
        if isinstance(e, Number):
            if isinstance(e, NaN): return NaN()
            if isinstance(e, Real):
                return Real(decimal_math.pow(b._as_decimal(), e.num))
            if e.is_negative:
                # (3/4)**-2 -> (4/3)**2
                ne = -e
                if isinstance(ne, One):
                    return Rational(b.q, b.p)
                return Rational(b.q, b.p) ** ne
            if isinstance(e, Infinity):
                if b.p > b.q:
                    # (3/2)**oo -> oo
                    return Infinity()
                if b.p < -b.q:
                    # (-3/2)**oo -> oo + I*oo
                    return Infinity() + Infinity() * ImaginaryUnit()
                return Zero()
            if isinstance(e, Integer):
                # (4/3)**2 -> 4**2 / 3**2
                return Rational(b.p ** e.p, b.q ** e.p)
            if isinstance(e, Rational):
                if b.p != 1:
                    # (4/3)**(5/6) -> 4**(5/6) * 3**(-5/6)
                    return Integer(b.p) ** e * Integer(b.q) ** (-e)
                if b >= 0:
                    x, xexact = integer_nthroot(b.p, e.q)
                    y, yexact = integer_nthroot(b.q, e.q)
                    if xexact and yexact:
                        res = Rational(x ** abs(e.p), y ** abs(e.p))
                        if e >= 0:
                            return res
                        else:
                            return 1/res
                    # Now check also devisors of the exponents denominator
                    # TODO: Check if this slows down to much.
                    for i in xrange(2, e.q/2 + 1):
                        if e.q % i == 0:
                            x, xexact = integer_nthroot(b.p, i)
                            y, yexact = integer_nthroot(b.q, i)
                            if xexact and yexact:
                                return Rational(x, y)**Rational(e.p, e.q/i)
                    else:
                        # Try to get some part of the base out, if exp > 1
                        if e.p > e.q:
                            i = e.p / e.q
                            r = e.p % e.q
                            return b**i * b**Rational(r, e.q)
                else:
                    return (-1)**e * (-b)**e

        c,t = b.as_coeff_terms()
        if e.is_even and isinstance(c, Basic.Number) and c < 0:
            return (-c * Basic.Mul(*t)) ** e

        return

    def _as_decimal(self):
        return decimal.Decimal(self.p) / decimal.Decimal(self.q)

    def __abs__(self):
        return Rational(abs(self.p), self.q)

    def __int__(self):
        return int(self.p//self.q)

    def __eq__(self, other):
        other = Basic.sympify(other)
        if isinstance(other, NumberSymbol):
            if other.is_irrational: return False
            return other.__eq__(self)
        if other.is_comparable and not isinstance(other, Rational): other = other.evalf()
        if isinstance(other, Number):
            if isinstance(other, Real):
                return bool(self._as_decimal()==other._as_decimal())
            return bool(self.p==other.p and self.q==other.q)
        return RelMeths.__eq__(self, other)
    def __ne__(self, other):
        other = Basic.sympify(other)
        if isinstance(other, NumberSymbol):
            if other.is_irrational: return True
            return other.__ne__(self)
        if other.is_comparable and not isinstance(other, Rational): other = other.evalf()
        if isinstance(other, Number):
            if isinstance(other, Real):
                return bool(self._as_decimal()!=other._as_decimal())
            return bool(self.p!=other.p or self.q!=other.q)
        return RelMeths.__ne__(self, other)
    def __lt__(self, other):
        other = Basic.sympify(other)
        if isinstance(other, NumberSymbol):
            return other.__ge__(self)
        if other.is_comparable and not isinstance(other, Rational): other = other.evalf()
        if isinstance(other, Number):
            if isinstance(other, Real):
                return bool(self._as_decimal() < other._as_decimal())
            return bool(self.p * other.q < self.q * other.p)
        return RelMeths.__lt__(self, other)
    def __le__(self, other):
        other = Basic.sympify(other)
        if isinstance(other, NumberSymbol):
            return other.__gt__(self)
        if other.is_comparable and not isinstance(other, Rational): other = other.evalf()
        if isinstance(other, Number):
            if isinstance(other, Real):
                return bool(self._as_decimal()<=other._as_decimal())
            return bool(self.p * other.q <= self.q * other.p)
        return RelMeths.__le__(self, other)

    def factors(self):
        f = factor_trial_division(self.p).copy()
        for p,e in factor_trial_division(self.q).items():
            try: f[p] += -e
            except KeyError: f[p] = -e
        fi = {}
        for p,e in f.items():
            if e==0:
                del f[p]
            else:
                try: fi[e] *= p
                except KeyError: fi[e] = p
        f = {}
        for e,p in fi.items():
            f[p] = e
        if len(f)>1 and f.has_key(1): del f[1]
        return f

    def as_numer_denom(self):
        return Integer(self.p), Integer(self.q)

class Integer(Rational):

    q = 1
    is_integer = True

    @Memoizer(type, (int, long))
    def __new__(cls, i):
        if isinstance(i, Integer):
            return i
        if i==0: return Zero()
        if i==1: return One()
        if i==-1: return NegativeOne()
        obj = Basic.__new__(cls)
        obj.p = i
        return obj

    def _eval_is_odd(self):
        return bool(self.p % 2)

    @property
    def precedence(self):
        if self.p < 0:
            return 40 # same as Add
        return Atom.precedence

    def tostr(self, level=0):
        if self.precedence<=level:
            return '(%s)' % (self.p)
        return str(self.p)

    def torepr(self):
        return '%s(%r)' % (self.__class__.__name__, self.p)

    def _eval_power(b, e):
        if isinstance(e, Number):
            if isinstance(e, NaN): return NaN()
            if isinstance(e, Real):
                return Real(decimal_math.pow(b._as_decimal(), e.num))
            if e.is_negative:
                # (3/4)**-2 -> (4/3)**2
                ne = -e
                if isinstance(ne, One):
                    return Rational(1, b.p)
                return Rational(1, b.p) ** ne
            if isinstance(e, Infinity):
                if b.p > 1:
                    # (3)**oo -> oo
                    return Infinity()
                if b.p < -1:
                    # (-3)**oo -> oo + I*oo
                    return Infinity() + Infinity() * ImaginaryUnit()
                return Zero()
            if isinstance(e, Integer):
                # (4/3)**2 -> 4**2 / 3**2
                return Integer(b.p ** e.p)
            if isinstance(e, Rational):
                if b == -1:  # any one has tested this ???
                    # calculate the roots of -1
                    r = cos(pi/e.q) + ImaginaryUnit()*sin(pi/e.q)
                    return r**e.p
                if b >= 0:
                    x, xexact = integer_nthroot(b.p, e.q)
                    if xexact:
                        res = Integer(x ** abs(e.p))
                        if e >= 0:
                            return res
                        else:
                            return 1/res
                    # Now check also devisors of the exponents denominator
                    # TODO: Check if this slows down to much.
                    for i in xrange(2, e.q/2 + 1):
                        if e.q % i == 0:
                            x, xexact = integer_nthroot(b.p, i)
                            if xexact:
                                return Integer(x)**Integer(b.p * i)
                    else:
                        # Try to get some part of the base out, if exp > 1
                        if e.p > e.q:
                            i = e.p / e.q
                            r = e.p % e.q
                            return b**i * b**Rational(r, e.q)
                else:
                    return (-1)**e * (-b)**e

        c,t = b.as_coeff_terms()
        if e.is_even and isinstance(c, Basic.Number) and c < 0:
            return (-c * Basic.Mul(*t)) ** e

        return

    def _eval_is_prime(self):
        if self.p < 0:
            return False

    def as_numer_denom(self):
        return self, One()

    def __floordiv__(self, other):
        return Integer(self.p // Integer(other).p)

    def __rfloordiv__(self, other):
        return Integer(Integer(other).p // self.p)

class Zero(Singleton, Integer):

    p = 0
    q = 1
    is_positive = False
    is_negative = False
    is_finite = False
    is_zero = True
    is_prime = False
    is_composite = True

    def _eval_power(b, e):
        if e.is_negative:
            return Infinity()
        if e.is_positive:
            return b
        d = e.evalf()
        if isinstance(d, Number):
            if d.is_negative:
                return Infinity()
            return b
        coeff, terms = e.as_coeff_terms()
        if coeff.is_negative:
            return Infinity() ** Basic.Mul(*terms)
        if not isinstance(coeff, Basic.One):
            return b ** Basic.Mul(*terms)

    def _eval_order(self, *symbols):
        # Order(0,x) -> 0
        return self

class One(Singleton, Integer):

    p = 1
    q = 1

    is_prime = True

    def _eval_power(b, e):
        return b

    def _eval_order(self, *symbols):
        return

class NegativeOne(Singleton, Integer):

    p = -1
    q = 1

    def _eval_power(b, e):
        if e.is_odd: return NegativeOne()
        if e.is_even: return One()
        if isinstance(e, Number):
            if isinstance(e, Real):
                a = e.num * decimal_math.pi()
                s = decimal_math.sin(a)
                c = decimal_math.cos(a)
                return Real(s) + Real(c) * ImaginaryUnit()
            if isinstance(e, NaN):
                return NaN()
            if isinstance(e, (Infinity, NegativeInfinity)):
                return NaN()
            if isinstance(e, Half):
                return ImaginaryUnit()
            if isinstance(e, Rational):
                if e.q == 2:
                    return ImaginaryUnit() ** Integer(e.p)
                q = int(e)
                if q:
                    q = Integer(q)
                    return b ** q * b ** (e - q)
        return

class Half(Singleton, Rational):

    p = 1
    q = 2

class Infinity(Singleton, Rational):

    p = 1
    q = 0

    is_commutative = True
    is_positive = True
    is_bounded = False
    is_finite = None
    is_odd = None

    def tostr(self, level=0):
        return 'oo'

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
            if isinstance(e, NaN):
                return NaN()
        d = e.evalf()
        if isinstance(d, Number):
            return b ** d
        return

    def _as_decimal(self):
        return decimal.Decimal('Infinity')

class NegativeInfinity(Singleton, Rational):

    p = -1
    q = 0

    is_commutative = True
    is_real = True
    is_positive = False
    is_bounded = False
    is_finite = False

    precedence = 40 # same as Add

    def tostr(self, level=0):
        return '-oo'

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
            if isinstance(e, (NaN, Infinity, NegativeInfinity)):
                return NaN()
            if isinstance(e, Integer):
                if e.is_positive:
                    if e.is_odd:
                        return NegativeInfinity()
                    return Infinity()
            return NegativeOne()**e * Infinity() ** e
        return

    def _as_decimal(self):
        return decimal.Decimal('-Infinity')

class NaN(Singleton, Rational):

    p = 0
    q = 0

    is_commutative = True
    is_real = None
    is_comparable = None
    is_bounded = None
    #is_unbounded = False

    def tostr(self, level=0):
        return 'nan'

    def _as_decimal(self):
        return decimal.Decimal('NaN')

    def _eval_power(b, e):
        if isinstance(e, Basic.Zero):
            return S.One
        return b

class NumberSymbol(Singleton, Atom, RelMeths, ArithMeths):

    is_commutative = True
    is_comparable = True
    is_bounded = True
    is_finite = True

    def approximation(self, number_cls):
        """ Return an interval with number_cls endpoints
        that contains the value of NumberSymbol.
        If not implemented, then return None.
        """

    def _eval_derivative(self, s):
        return Zero()
    def __eq__(self, other):
        other = Basic.sympify(other)
        if self is other: return True
        if isinstance(other, Number) and self.is_irrational: return False
        return RelMeths.__eq__(self, other)
    def __ne__(self, other):
        other = Basic.sympify(other)
        if self is other: return False
        if isinstance(other, Number) and self.is_irrational: return True
        return RelMeths.__ne__(self, other)
    def __lt__(self, other):
        other = Basic.sympify(other)
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
        return RelMeths.__lt__(self, other)
    def __le__(self, other):
        other = Basic.sympify(other)
        if self is other: return True
        if other.is_comparable: other = other.evalf()
        if isinstance(other, Number):
            return self.evalf()<=other
        return RelMeths.__le__(self, other)
    def __gt__(self, other):
        return (-self) < (-other)
    def __ge__(self, other):
        return (-self) <= (-other)


class Exp1(NumberSymbol):

    is_real = True
    is_positive = True
    is_negative = False # XXX Forces is_negative/is_nonnegative
    is_irrational = True

    def tostr(self, level=0):
        return 'E'

    def _eval_evalf(self):
        return Real(decimal_math.e())

    def approximation_interval(self, number_cls):
        if issubclass(number_cls,Integer):
            return (Integer(2),Integer(3))
        elif issubclass(number_cls,Rational):
            pass

    def _eval_power(self, exp):
        return Basic.Exp()(exp)

class Pi(NumberSymbol):

    is_real = True
    is_positive = True
    is_negative = False
    is_irrational = True

    def _eval_evalf(self):
        return Real(decimal_math.pi())

    def approximation_interval(self, number_cls):
        if issubclass(number_cls, Integer):
            return (Integer(3), Integer(4))
        elif issubclass(number_cls, Rational):
            return (Rational(223,71), Rational(22,7))

    def tostr(self, level=0):
        return 'Pi'

class GoldenRatio(NumberSymbol):

    is_real = True
    is_positive = True
    is_negative = False
    is_irrational = True

    def _eval_evalf(self):
        return Real(decimal_math.golden_ratio())

    def approximation_interval(self, number_cls):
        if issubclass(number_cls, Integer):
            return (S.One, Rational(2))
        elif issubclass(number_cls, Rational):
            pass

    def tostr(self, level=0):
        return 'GoldenRatio'

class EulerGamma(NumberSymbol):

    is_real = True
    is_positive = True
    is_negative = False
    is_irrational = None

    def _eval_evalf(self):
        return

    def approximation_interval(self, number_cls):
        if issubclass(number_cls, Integer):
            return (S.Zero, S.One)
        elif issubclass(number_cls, Rational):
            return (S.Half, Rational(3, 5))

    def tostr(self, level=0):
        return 'EulerGamma'

class Catalan(NumberSymbol):

    is_real = True
    is_positive = True
    is_negative = False
    is_irrational = None

    def _eval_evalf(self):
        return

    def approximation_interval(self, number_cls):
        if issubclass(number_cls, Integer):
            return (S.Zero, S.One)
        elif issubclass(number_cls, Rational):
            return (Rational(9, 10), S.One)

    def tostr(self, level=0):
        return 'Catalan'

class ImaginaryUnit(Singleton, Atom, RelMeths, ArithMeths):

    is_commutative = True
    is_imaginary = True
    is_bounded = True
    is_finite = True

    def tostr(self, level=0):
        return 'I'

    def _eval_conjugate(self):
        return -Basic.ImaginaryUnit()

    def _eval_derivative(self, s):
        return Zero()

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
                if e==0: return One()
                if e==1: return ImaginaryUnit()
                if e==2: return -One()
                return -ImaginaryUnit()
            return -One() ** (e * Half())
        return

    def as_base_exp(self):
        return -One(),Rational(1,2)

Basic.singleton['E'] = Exp1
Basic.singleton['pi'] = Pi
Basic.singleton['Pi'] = Pi
Basic.singleton['I'] = ImaginaryUnit
Basic.singleton['oo'] = Infinity
Basic.singleton['nan'] = NaN

Basic.singleton['GoldenRatio'] = GoldenRatio
Basic.singleton['EulerGamma'] = EulerGamma
Basic.singleton['Catalan'] = Catalan
