"""
This module implements a class Float for arbitrary-precision binary
floating-point arithmetic. It is typically 10-100 times faster
than Python's standard Decimals. For details on usage, refer to the
docstrings in the Float class.
"""

import math
_clog = math.log
_csqrt = math.sqrt

from sympy import Rational
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
    Return a tuple containing the new (man, exp).
    """
    if man == 0:
        return man, 0
    bc = bitcount(man)
    if bc > prec:
        man = rshift(man, bc-prec, mode)
        exp += (bc - prec)
    # It is not necessary to strip trailing zeros, but this
    # standardization permits faster equality testing of numbers
    # with the same exponent
    tr = trailing_zeros(man)
    if tr:
        man >>= tr
        exp += tr
    if man == 0:
        exp = 0
    return man, exp

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

class Float(object):
    """
    A Float is a rational number of the form

        man * 2**exp

    ("man" and "exp" are short for "mantissa" and "exponent"). Both man
    and exp are integers, possibly negative, and may be arbitrarily large.
    Essentially, a larger mantissa corresponds to a higher precision
    and a larger exponent corresponds to larger magnitude.

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

    __slots__ = ["man", "exp"]

    def __init__(s, x=0, prec=None, mode=None):
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
        prec = prec or s._prec
        mode = mode or s._mode
        if isinstance(x, tuple):
            s.man, s.exp = normalize(x[0], x[1], prec, mode)
        elif isinstance(x, Float):
            s.man, s.exp = normalize(x.man, x.exp, prec, mode)
        elif isinstance(x, (int, long)):
            s.man, s.exp = normalize(x, 0, prec, mode)
        elif isinstance(x, float):
            m, e = math.frexp(x)
            s.man, s.exp = normalize(int(m*2**53), e-53, prec, mode)
        elif isinstance(x, (str, Rational)):
            if isinstance(x, str):
                x = Rational(x)
            n = prec + bitcount(x.q) + 2
            s.man, s.exp = normalize((x.p<<n)//x.q, -n, prec, mode)
        else:
            raise TypeError

    def __pos__(s):
        """s.__pos__() <==> +s

        Normalize s to the current working precision, rounding according
        to the current rounding mode."""
        return Float(s)

    def __repr__(s):
        """Represent s as a decimal string, with sufficiently many
        digits included to ensure that Float(repr(s)) == s at the
        current working precision."""
        st = "Float('%s')"
        return st % binary_to_decimal(s.man, s.exp, Float._dps + 2)

    def __str__(s):
        """Print slightly more prettily than __repr__"""
        return binary_to_decimal(s.man, s.exp, Float._dps)

    def __float__(s):
        """Convert s to a Python float. OverflowError will be raised
        if the magnitude of s is too large."""
        try:
            return math.ldexp(s.man, s.exp)
        # Handle case when mantissa has too many bits (will still
        # overflow if exp is large)
        except OverflowError:
            n = bitcount(s.man) - 64
            m = s.man >> n
            return math.ldexp(m, s.exp + n)

    def __int__(s):
        return rshift(s.man, -s.exp, 0)

    #------------------------------------------------------------------
    # Comparison
    #

    def __cmp__(s, t):
        """__cmp__(s, t) <==> cmp(s, t)

        Returns -1 if s < t, 0 if s == t, and 1 if s > t"""
        if not isinstance(t, Float):
            t = Float(t)
        sm, se, tm, te = s.man, s.exp, t.man, t.exp
        if tm == 0: return cmp(sm, 0)
        if sm == 0: return cmp(0, tm)
        if sm > 0 and tm < 0: return 1
        if sm < 0 and tm > 0: return -1
        if se == te: return cmp(sm, tm)
        a = bitcount(sm) + se
        b = bitcount(tm) + te
        if sm > 0:
            if a < b: return -1
            if a > b: return 1
        else:
            if a < b: return 1
            if a < b: return -1
        return cmp((s-t).man, 0)

    def ae(s, t, rel_eps=None, abs_eps=None):
        """
        "ae" is short for "almost equal"

        Determine whether the difference between s and t is smaller
        than a given epsilon.

        Both a maximum relative difference and a maximum difference
        (or 'epsilons') may be specified. The absolute difference is
        defined as |s-t| and the relative difference is defined
        as |s-t|/max(|s|, |t|).

        If only one epsilon is given, both are set to the same value.
        If none is given, both epsilons are set to 2**(-prec+m) where
        prec is the current working precision and m is a small integer.
        """
        if not isinstance(t, Float):
            t = Float(t)
        if abs_eps is None and rel_eps is None:
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
            err = diff/t
        else:
            err = diff/s
        return err <= rel_eps

    def almost_zero(s, prec):
        """Quick check if |s| < 2**-prec"""
        return bitcount(s.man) + s.exp < prec

    def __nonzero__(s):
        return bool(s.man)

    #------------------------------------------------------------------
    # Arithmetic
    #

    def __abs__(s):
        if s.man < 0:
            return -s
        return s

    def __add__(s, t):
        if isinstance(t, Float):
            if t.exp > s.exp:
                s, t = t, s
            return Float((t.man+(s.man<<(s.exp-t.exp)), t.exp))
        if isinstance(t, (int, long)):
            # XXX: cancellation is possible here
            return s + Float(t)

    __radd__ = __add__

    def __neg__(s):
        return Float((-s.man, s.exp))

    def __sub__(s, t):
        return s + (-t)

    def __rsub__(s, t):
        return (-s) + t

    def __mul__(s, t):
        if isinstance(t, Float):
            return Float((s.man*t.man, s.exp+t.exp))
        if isinstance(t, (int, long)):
            return 

    __rmul__ = __mul__

    def __div__(s, t):
        if t == 0:
            raise ZeroDivisionError
        if isinstance(t, Float):
            extra = max(0, s._prec - bitcount(s.man) + bitcount(t.man) + 4)
            return Float(((s.man<<extra)//t.man, s.exp-t.exp-extra))
        if isinstance(t, (int, long)):
            extra = s._prec - bitcount(s.man) + bitcount(t) + 4
            return Float(((s.man<<extra)//t, s.exp-extra))

    def __pow__(s, n):
        """Calculate (man*2**exp)**n, currently for integral n only."""
        if isinstance(n, (int, long)):
            if n == 0: return Float((1, 0))
            if n == 1: return +s
            if n == 2: return s * s
            if n == -1: return 1 / s
            if n < 0:
                Float._prec += 2
                r = 1 / (s ** -n)
                Float._prec -= 2
                return +r
            else:
                prec2 = Float._prec + int(4*_clog(n, 2) + 4)
                man, exp = normalize(s.man, s.exp, prec2, ROUND_FLOOR)
                pm, pe = 1, 0
                while n:
                    if n & 1:
                        pm, pe = normalize(pm*man, pe+exp, prec2, ROUND_FLOOR)
                        n -= 1
                    man, exp = normalize(man*man, exp+exp, prec2, ROUND_FLOOR)
                    n = n // 2
                return Float((pm, pe))
        if n == 0.5:
            return s.sqrt()
        raise ValueError

