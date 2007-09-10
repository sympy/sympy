"""
This module implements a class Float for arbitrary-precision binary
floating-point arithmetic. It is typically 10-100 times faster
than Python's standard Decimals. For details on usage, refer to the
docstrings in the Float class. A ComplexFloat class is also provided.
"""

import math
_clog = math.log
_csqrt = math.sqrt

from sympy import Basic, Rational
from utils_ import bitcount, trailing_zeros

#----------------------------------------------------------------------
# Rounding modes
#

class _RoundingMode(int):
    def __new__(cls, level, name):
        a = int.__new__(cls, level)
        a.name = name
        return a
    def __repr__(self):
        return self.name

ROUND_DOWN    = _RoundingMode(1, 'ROUND_DOWN')
ROUND_UP      = _RoundingMode(2, 'ROUND_UP')
ROUND_FLOOR   = _RoundingMode(3, 'ROUND_FLOOR')
ROUND_CEILING = _RoundingMode(4, 'ROUND_CEILING')
ROUND_HALF_UP = _RoundingMode(5, 'ROUND_HALF_UP')
ROUND_HALF_DOWN = _RoundingMode(6, 'ROUND_HALF_DOWN')
ROUND_HALF_EVEN = _RoundingMode(7, 'ROUND_HALF_EVEN')


#----------------------------------------------------------------------
# Helper functions for bit manipulation
#

def rshift(x, n, mode):
    """
    Shift x n bits to the right (i.e., calculate x/(2**n)), and
    round to the nearest integer in accordance with the specified
    rounding mode. The exponent n may be negative, in which case x is
    shifted to the left (and no rounding is necessary).
    """
    if n == 0 or x == 0: return x
    if n < 0:  return x << -n

    # Bit-fiddling is relatively expensive in Python. To get away easily, we
    # can exploit the fact that Python rounds positive integers toward
    # zero and negative integers away from zero when dividing/shifting

    # These cases can be handled by simple shifting
    if mode < ROUND_HALF_UP:
        if mode == ROUND_DOWN:
            if x > 0: return x >> n
            else:     return -((-x) >> n)
        if mode == ROUND_UP:
            if x > 0: return -((-x) >> n)
            else:     return x >> n
        if mode == ROUND_FLOOR:
            return x >> n
        if mode == ROUND_CEILING:
            return -((-x) >> n)

    # Here we need to inspect the bits around the cutoff point
    if x > 0: t = x >> (n-1)
    else:     t = (-x) >> (n-1)
    if t & 1:
        if mode == ROUND_HALF_UP or \
           (mode == ROUND_HALF_DOWN and x & ((1<<(n-1))-1)) or \
           (mode == ROUND_HALF_EVEN and (t&2 or x & ((1<<(n-1))-1))):
            if x > 0:  return (t>>1)+1
            else:      return -((t>>1)+1)
    if x > 0: return t>>1
    else:     return -(t>>1)

def normalize(man, exp, prec, mode):
    """
    Normalize the binary floating-point number represented by
    man * 2**exp to the specified precision level, rounding
    according to the specified rounding mode if necessary.
    The mantissa is also stripped of trailing zero bits, and
    its bits are counted. The returned value is a tuple
    (man, exp, bc).
    """
    if not man:
        return 0, 0, 0
    bc = bitcount(man)
    if bc > prec:
        man = rshift(man, bc-prec, mode)
        exp += (bc - prec)
        bc = prec
    # Stripping trailing zeros permits faster equality testing
    if not man & 1:
        tr = trailing_zeros(man)
        if tr:
            man >>= tr
            exp += tr
            bc -= tr
    if not man:
        return 0, 0, 0
    return man, exp, bc

# Shortcut for constructing a Float from given mantissa and exponent, used by
# some of the arithmetic operators in the Float class
def makefloat(man, exp):
    return tuple.__new__(Float, normalize(man, exp, Float._prec, Float._mode))


#----------------------------------------------------------------------
# Other helper functions
#

def binary_to_decimal(man, exp, n):
    """Represent as a decimal string with at most n digits"""
    import decimal
    prec_ = decimal.getcontext().prec
    decimal.getcontext().prec = n
    if exp >= 0: d = decimal.Decimal(man) * (1<<exp)
    else:        d = decimal.Decimal(man) / (1<<-exp)
    a = str(d)
    decimal.getcontext().prec = prec_
    return a

_convratio = _clog(10,2) # 3.3219...


#---------------------------------------------------------------------------#
#                              Float class                                  #
#---------------------------------------------------------------------------#

class Float(tuple):
    """
    A Float is a rational number of the form

        man * 2**exp

    ("man" and "exp" are short for "mantissa" and "exponent"). Both man
    and exp are integers, possibly negative, and may be arbitrarily large.
    Essentially, a larger mantissa corresponds to a higher precision
    and a larger exponent corresponds to larger magnitude.

    A Float instance is represented by a tuple

        (man, exp, bc)

    where bc is the bitcount of the mantissa. The elements can be
    accessed as named properties:

        >>> x = Float(3)
        >>> x.man
        3
        >>> x.exp
        0
        >>> x.bc
        2

    When performing an arithmetic operation on two Floats, or creating a
    new Float from an existing numerical value, the result gets rounded
    to a fixed precision level, just like with ordinary Python floats.
    Unlike regular Python floats, however, the precision level can be
    set arbitrarily high. You can also change the rounding mode (all
    modes supported by Decimal are also supported by Float).

    The precision level and rounding mode are stored as properties of
    the Float class. (An individual Float instances does not have any
    precision or rounding data associated with it.) The precision level
    and rounding mode make up the current working context. You can
    change the working context through static methods of the Float
    class:

        Float.setprec(n)    -- set precision to n bits
        Float.extraprec(n)  -- increase precision by n bits
        Float.setdps(n)     -- set precision equivalent to n decimals
        Float.setmode(mode) -- set rounding mode

    Corresponding methods are available for inspection:

        Float.getprec()
        Float.getdps()
        Float.getmode()

    There are also two methods Float.store() and Float.revert(). If
    you call Float.store() before changing precision or mode, the
    old context can be restored with Float.revert(). (If Float.revert()
    is called one time too much, the default settings are restored.)
    You can nest multiple uses of store() and revert().

    (In the future, it will also be possible to use the 'with'
    statement to change contexts.)

    Note that precision is measured in bits. Since the ratio between
    binary and decimal precision levels is irrational, setprec and
    setdps work slightly differently. When you set the precision with
    setdps, the bit precision is set slightly higher than the exact
    corresponding precision to account for the fact that decimal
    numbers cannot generally be represented exactly in binary (the
    classical example is 0.1). The exact number given to setdps
    is however used by __str__ to determine number of digits to
    display. Likewise, when you set a bit precision, the decimal
    printing precision used for __str__ is set slightly lower.

    The following rounding modes are available:

        ROUND_DOWN       -- toward zero
        ROUND_UP         -- away from zero
        ROUND_FLOOR      -- towards -oo
        ROUND_CEILING    -- towards +oo
        ROUND_HALF_UP    -- to nearest; 0.5 to 1
        ROUND_HALF_DOWN  -- to nearest; 0.5 to 0
        ROUND_HALF_EVEN  -- to nearest; 0.5 to 0 and 1.5 to 2

    The rounding modes are available both as global constants defined
    in this module and as properties of the Float class, e.g.
    Float.ROUND_CEILING.

    The default precision level is 53 bits and the default rounding
    mode is ROUND_HALF_EVEN. In this mode, Floats should round exactly
    like regular Python floats (in the absence of bugs!).
    """

    #------------------------------------------------------------------
    # Static methods for context management
    #

    # Also make these constants available from the class
    ROUND_DOWN = ROUND_DOWN
    ROUND_UP = ROUND_UP
    ROUND_FLOOR = ROUND_FLOOR
    ROUND_CEILING = ROUND_CEILING
    ROUND_HALF_UP = ROUND_HALF_UP
    ROUND_HALF_DOWN = ROUND_HALF_DOWN
    ROUND_HALF_EVEN = ROUND_HALF_EVEN

    _prec = 53
    _dps = 15
    _mode = ROUND_HALF_EVEN
    _stack = []

    @staticmethod
    def store():
        """Store the current precision/rounding context. It can
        be restored by calling Float.revert()"""
        Float._stack.append((Float._prec, Float._dps, Float._mode))

    @staticmethod
    def revert():
        """Revert to last precision/rounding context stored with
        Float.store()"""
        if Float._stack:
            Float._prec, Float._dps, Float._mode = Float._stack.pop()
        else:
            Float._prec, Float._dps, Float._mode = 53, 15, ROUND_HALF_EVEN

    @staticmethod
    def setprec(n):
        """Set precision to n bits"""
        n = int(n)
        Float._prec = n
        Float._dps = int(round(n/_convratio)-1)

    @staticmethod
    def setdps(n):
        """Set the precision high enough to allow representing numbers
        with at least n decimal places without loss."""
        n = int(n)
        Float._prec = int(round((n+1)*_convratio))
        Float._dps = n

    @staticmethod
    def extraprec(n):
        Float.setprec(Float._prec + n)

    @staticmethod
    def setmode(mode):
        assert isinstance(mode, _RoundingMode)
        Float._mode = mode

    @staticmethod
    def getprec(): return Float._prec

    @staticmethod
    def getdps(): return Float._dps

    @staticmethod
    def getmode(): return Float._mode


    #------------------------------------------------------------------
    # Core object functionality
    #

    man = property(lambda self: self[0])
    exp = property(lambda self: self[1])
    bc = property(lambda self: self[2])

    def __new__(cls, x=0, prec=None, mode=None):
        """
        Float(x) creates a new Float instance with value x. The usual
        types are supported for x:

            >>> Float(3)
            Float('3')
            >>> Float(3.5)
            Float('3.5')
            >>> Float('3.5')
            Float('3.5')
            >>> Float(Rational(7,2))
            Float('3.5')

        You can also create a Float from a tuple specifying its
        mantissa and exponent:

            >>> Float((5, -3))
            Float('0.625')

        Use the prec and mode arguments to specify a custom precision
        level (in bits) and rounding mode. If these arguments are
        omitted, the current working precision is used instead.

            >>> Float('0.500001', prec=3, mode=ROUND_DOWN)
            Float('0.5')
            >>> Float('0.500001', prec=3, mode=ROUND_UP)
            Float('0.625')

        """
        prec = prec or cls._prec
        mode = mode or cls._mode
        if isinstance(x, tuple):
            return tuple.__new__(cls, normalize(x[0], x[1], prec, mode))
        elif isinstance(x, (int, long)):
            return tuple.__new__(cls, normalize(x, 0, prec, mode))
        elif isinstance(x, float):
            # We assume that a float mantissa has 53 bits
            m, e = math.frexp(x)
            return tuple.__new__(cls, normalize(int(m*(1<<53)), e-53, prec, mode))
        elif isinstance(x, (str, Rational)):
            if isinstance(x, str):
                x = Rational(x)
            n = prec + bitcount(x.q) + 2
            return tuple.__new__(cls, normalize((x.p<<n)//x.q, -n, prec, mode))
        else:
            raise TypeError

    def __hash__(s):
        try:
            # Try to be compatible with hash values for floats and ints
            return hash(float(s))
        except OverflowError:
            # We must unfortunately sacrifice compatibility with ints here. We
            # could do hash(man << exp) when the exponent is positive, but
            # this would cause unreasonable inefficiency for large numbers.
            return hash(self.man) + hash(self.exp)

    def __pos__(s):
        """s.__pos__() <==> +s

        Normalize s to the current working precision, rounding according
        to the current rounding mode."""
        return Float(s)

    def __float__(s):
        """Convert s to a Python float. OverflowError will be raised
        if the magnitude of s is too large."""
        try:
            return math.ldexp(s.man, s.exp)
        # Handle case when mantissa has too many bits (will still
        # overflow if exp is too large)
        except OverflowError:
            n = s.bc - 64
            m = s.man >> n
            return math.ldexp(m, s.exp + n)

    def __int__(s):
        """Convert to a Python int, using floor rounding"""
        return rshift(s.man, -s.exp, 0)

    def rational(s):
        """Convert to a SymPy Rational"""
        if s.exp > 0:
            return Rational(s.man * 2**s.exp, 1)
        else:
            return Rational(s.man, 2**(-s.exp))

    def __repr__(s):
        """Represent s as a decimal string, with sufficiently many
        digits included to ensure that Float(repr(s)) == s at the
        current working precision."""
        st = "Float('%s')"
        return st % binary_to_decimal(s.man, s.exp, Float._dps + 2)

    def __str__(s):
        """Print slightly more prettily than __repr__"""
        return binary_to_decimal(s.man, s.exp, Float._dps)


    #------------------------------------------------------------------
    # Comparisons
    #

    def __cmp__(s, t):
        """s.__cmp__(t) <==> cmp(s, t)

        If t is a Float, int or float, this returns -1 if s < t, 0 if
        s == t, and 1 if s > t. The comparison operators >, >=, <, <=
        are defined in terms of this function. If t is a SymPy Basic
        instance, s is converted into a Rational and Rational.__cmp__
        is called.

        Warning: in extreme cases, the truncation error resulting from
        calling Float(t) will result in an erroneous comparison: for
        example, Float(2**80) will compare as equal to 2**80+1. This
        problem can be circumvented by manually increasing the working
        precision or by converting numbers into Rationals for
        comparisons.
        """

        if not isinstance(t, Float):
            if isinstance(t, Basic):
                return s.rational().__cmp__(t)
            t = Float(t)

        # Note: the reason we call cmp(x,y) below instead of directly forming
        # the integer x-y is that the output from __cmp__ must be a machine
        # size integer in Python.

        # An inequality between two numbers s and t is determined by looking
        # at the value of s-t. A full floating-point subtraction is relatively
        # slow, so we first try to look at the exponents and signs of s and t.
        sm, se, sbc = s
        tm, te, tbc = t

        # Very easy cases: check for 0's and opposite signs
        if not tm: return cmp(sm, 0)
        if not sm: return cmp(0, tm)
        if sm > 0 and tm < 0: return 1
        if sm < 0 and tm > 0: return -1

        # In this case, the numbers likely have the same magnitude
        if se == te: return cmp(sm, tm)

        # The numbers have the same sign but different exponents. In this
        # case we try to determine if they are of different magnitude by
        # checking the position of the highest set bit in each number.
        a = sbc + se
        b = tbc + te
        if sm > 0:
            if a < b: return -1
            if a > b: return 1
        else:
            if a < b: return 1
            if a < b: return -1

        # The numbers have similar magnitude but different exponents.
        # So we subtract and check the sign of resulting mantissa.
        return cmp((s-t)[0], 0)


    # Since Float inherits from tuple, these functions must be defined
    # in addition to __cmp__ to override the tuple methods.
    __lt__ = lambda s, t: s.__cmp__(t) < 0
    __gt__ = lambda s, t: s.__cmp__(t) > 0
    __le__ = lambda s, t: s.__cmp__(t) <= 0
    __ge__ = lambda s, t: s.__cmp__(t) >= 0


    """
    We implement __eq__ separately from __cmp__. One reason for doing this
    performance: due to the way Floats are normalized, two Floats are
    mathematically equal iff all three parts (man, exp, bc) are equal. So we
    can do a direct tuple comparison and avoid the machinery in __cmp__.

    Another reason is to support == and != between Floats and ComplexFloat.
    ComplexFloats are not supported by __cmp__, since complex numbers are
    unordered.
    """
    def __eq__(s, t):
        """s.__eq__(t) <==> s == Float(t)

        Determine whether s and Float(t) are equal (see warning for
        __cmp__ about conversion between different types.)"""
        if isinstance(t, Basic):
            return s.rational() == t
        if isinstance(t, Float):
            return s[:] == t[:]
        if isinstance(t, ComplexFloat):
            return ComplexFloat(s) == t
        return not s.__cmp__(t)

    __ne__ = lambda s, t: not s.__eq__(t)

    def ae(s, t, rel_eps=None, abs_eps=None):
        """
        Determine whether the difference between s and t is smaller
        than a given epsilon ("ae" is short for "almost equal").

        Both a maximum relative difference and a maximum difference
        ('epsilons') may be specified. The absolute difference is
        defined as |s-t| and the relative difference is defined
        as |s-t|/max(|s|, |t|).

        If only one epsilon is given, both are set to the same value.
        If none is given, both epsilons are set to 2**(-prec+m) where
        prec is the current working precision and m is a small integer.
        """
        if not isinstance(t, Float):
            t = Float(t)
        if abs_eps is None and rel_eps is None:
            rel_eps = tuple.__new__(Float, (1, -s._prec+4, 1))
            rel_eps = Float((1, -s._prec+4))
        if abs_eps is None:
            abs_eps = rel_eps
        elif rel_eps is None:
            rel_eps = abs_eps
        diff = abs(s-t)
        if diff <= abs_eps:
            return True
        abss = abs(s)
        abst = abs(t)
        if abss < abst:
            err = diff/abst
        else:
            err = diff/abss
        return err <= rel_eps

    def almost_zero(s, prec):
        """Quick check if |s| < 2**-prec. May return a false negative
        if s is very close to the threshold."""
        return s.bc + s.exp < prec

    def __nonzero__(s):
        return bool(s[0])

    #------------------------------------------------------------------
    # Arithmetic
    #

    def __abs__(s):
        if s[0] < 0:
            return -s
        return s

    def __add__(s, t):
        if isinstance(t, Float):
            if t[1] > s[1]:
                s, t = t, s
            sman, sexp, sbc = s
            tman, texp, tbc = t
            if not tman: return +s
            if not sman: return +t
            if sexp - texp > 100:
                bitdelta = (sbc+sexp)-(tbc+texp)
                if bitdelta > s._prec+5:
                    # XXX: handle rounding
                    return +s
            return makefloat(tman+(sman<<(sexp-texp)), texp)
        if isinstance(t, (int, long)):
            # XXX: cancellation is possible here
            return s + Float(t)
        if isinstance(t, (ComplexFloat, complex)):
            return ComplexFloat(s) + t
        return s + Float(t)

    __radd__ = __add__

    def __neg__(s):
        return makefloat(-s[0], s[1])

    def __sub__(s, t):
        if isinstance(t, Float):
            return s + tuple.__new__(Float, (-t[0],) + t[1:])
        else:
            return s + (-t)

    def __rsub__(s, t):
        return (-s) + t

    def __mul__(s, t):
        if isinstance(t, Float):
            sman, sexp, sbc = s
            tman, texp, tbc = t
            return makefloat(sman*tman, sexp+texp)
        if isinstance(t, (int, long)):
            return makefloat(s[0]*t, s[1])
        if isinstance(t, (ComplexFloat, complex)):
            return ComplexFloat(s) * t
        return s * Float(t)

    __rmul__ = __mul__

    def __div__(s, t):
        if t == 0:
            raise ZeroDivisionError
        if isinstance(t, Float):
            sman, sexp, sbc = s
            tman, texp, tbc = t
            extra = max(0, s._prec - sbc + tbc + 4)
            return makefloat((sman<<extra)//tman, sexp-texp-extra)
        if isinstance(t, (int, long)):
            sman, sexp, sbc = s
            extra = s._prec - sbc + bitcount(t) + 4
            return makefloat((sman<<extra)//t, sexp-extra)
        if isinstance(t, (ComplexFloat, complex)):
            return ComplexFloat(s) / t
        return s / Float(t)

    def __rdiv__(s, t):
        if isinstance(t, (ComplexFloat, complex)):
            return t / ComplexFloat(s)
        return Float(t) / s

    def __pow__(s, n):
        """Calculate (man*2**exp)**n, currently for integral n only."""
        if isinstance(n, (int, long)):
            if n == 0: return Float((1, 0))
            if n == 1: return +s
            if n == 2: return s * s
            if n == -1: return Float(1) / s
            if n < 0:
                Float._prec += 2
                r = Float(1) / (s ** -n)
                Float._prec -= 2
                return +r
            else:
                prec2 = Float._prec + int(4*_clog(n, 2) + 4)
                man, exp, bc = normalize(s.man, s.exp, prec2, ROUND_FLOOR)
                pm, pe, bc = 1, 0, 1
                while n:
                    if n & 1:
                        pm, pe, bc = normalize(pm*man, pe+exp, prec2, ROUND_FLOOR)
                        n -= 1
                    man, exp, _ = normalize(man*man, exp+exp, prec2, ROUND_FLOOR)
                    n = n // 2
                #return Float((pm, pe))
                return makefloat(pm, pe)
        # TODO: support arbitrary powers through exp
        if n == 0.5:
            from functions import sqrt
            return sqrt(s)
        raise ValueError


class ComplexFloat(object):

    def __init__(s, real=0, imag=0):
        if isinstance(real, (complex, ComplexFloat)):
            real, imag = real.real, real.imag
        s.real = Float(real)
        s.imag = Float(imag)

    def __repr__(s):
        r = repr(s.real)[6:-1]
        i = repr(s.imag)[6:-1]
        return "ComplexFloat(real=%s, imag=%s)" % (r, i)

    def __str__(s):
        return "(%s + %s*I)" % (s.real, s.imag)

    def __complex__(s):
        return complex(float(s.real), float(s.imag))

    def __pos__(s):
        return ComplexFloat(s.real, s.imag)

    def __abs__(s):
        from functions import hypot
        return hypot(s.real, s.imag)

    def __eq__(s, t):
        if not isinstance(t, ComplexFloat):
            t = ComplexFloat(t)
        return s.real == t.real and s.imag == t.imag

    # TODO: refactor and merge with Float.ae
    def ae(s, t, rel_eps=None, abs_eps=None):
        if not isinstance(t, ComplexFloat):
            t = ComplexFloat(t)
        if abs_eps is None and rel_eps is None:
            rel_eps = Float((1, -s.real._prec+4))
        if abs_eps is None:
            abs_eps = rel_eps
        elif rel_eps is None:
            rel_eps = abs_eps
        diff = abs(s-t)
        if diff <= abs_eps:
            return True
        abss = abs(s)
        abst = abs(t)
        if abss < abst:
            err = diff/abst
        else:
            err = diff/abss
        return err <= rel_eps

    def __nonzero__(s):
        return s.real != 0 or s.imag != 0

    def conjugate(s):
        return ComplexFloat(s.real, -s.imag)

    def __add__(s, t):
        if not isinstance(t, ComplexFloat):
            t = ComplexFloat(t)
        return ComplexFloat(s.real+t.real, s.imag+t.imag)

    __radd__ = __add__

    def __neg__(s):
        return ComplexFloat(-s.real, -s.imag)

    def __sub__(s, t):
        if not isinstance(t, ComplexFloat):
            t = ComplexFloat(t)
        return ComplexFloat(s.real-t.real, s.imag-t.imag)

    def __rsub__(s, t):
        return (-s) + t

    def __mul__(s, t):
        if not isinstance(t, ComplexFloat):
            t = ComplexFloat(t)
        a = s.real; b = s.imag; c = t.real; d = t.imag
        if b == d == 0:
            return ComplexFloat(a*c, 0)
        else:
            return ComplexFloat(a*c-b*d, a*d+b*c)

    __rmul__ = __mul__

    def __div__(s, t):
        if not isinstance(t, ComplexFloat):
            t = ComplexFloat(t)
        a = s.real; b = s.imag; c = t.real; d = t.imag
        mag = c*c + d*d
        return ComplexFloat((a*c+b*d)/mag, (b*c-a*d)/mag)

    def __rdiv__(s, t):
        return ComplexFloat(t) / s

    def __pow__(s, n):
        if n == 0: return ComplexFloat(1)
        if n == 1: return +s
        if n == -1: return 1/s
        if n == 2: return s*s
        if isinstance(n, (int, long)) and n > 0:
            # TODO: should increase working precision here
            w = ComplexFloat(1)
            while n:
                if n & 1:
                    w = w*s
                    n -= 1
                s = s*s
                n //= 2
            return w
        raise NotImplementedError
