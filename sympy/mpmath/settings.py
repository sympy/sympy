import os
import sys

# All supported rounding modes
round_nearest = intern('n')
round_floor = intern('f')
round_ceiling = intern('c')
round_up = intern('u')
round_down = intern('d')

round_fast = round_down

def prec_to_dps(n):
    """Return number of accurate decimals that can be represented
    with a precision of n bits."""
    return max(1, int(round(int(n)/3.3219280948873626)-1))

def dps_to_prec(n):
    """Return the number of bits required to represent n decimals
    accurately."""
    return max(1, int(round((int(n)+1)*3.3219280948873626)))

def repr_dps(n):
    """Return the number of decimal digits required to represent
    a number with n-bit precision so that it can be uniquely
    reconstructed from the representation."""
    dps = prec_to_dps(n)
    if dps == 15:
        return 17
    return dps + 3

#----------------------------------------------------------------------------#
# Support GMPY for high-speed large integer arithmetic.                      #
#                                                                            #
# To allow an external module to handle arithmetic, we need to make sure     #
# that all high-precision variables are declared of the correct type. MP_BASE#
# is the constructor for the high-precision type. It defaults to Python's    #
# long type but can be assinged another type, typically gmpy.mpz.            #
#                                                                            #
# MP_BASE must be used for the mantissa component of an mpf and must be used #
# for internal fixed-point operations.                                       #
#                                                                            #
# Side-effects                                                               #
# 1) "is" cannot be used to test for special values. Must use "==".          #
# 2) There are bugs in GMPY prior to v1.02 so we must use v1.03 or later.    #
#----------------------------------------------------------------------------#

# So we can import it from this module
gmpy = None

MODE = 'python'
MP_BASE = long
if 'MPMATH_NOGMPY' not in os.environ:
    try:
        import gmpy
        if gmpy.version() >= '1.03':
            MODE = 'gmpy'
            MP_BASE = gmpy.mpz
    except:
        pass

if os.environ.has_key('MPMATH_STRICT'):
    STRICT = True
else:
    STRICT = False

MP_BASE_TYPE = type(MP_BASE(0))
MP_ZERO = MP_BASE(0)
MP_ONE = MP_BASE(1)
MP_TWO = MP_BASE(2)
MP_THREE = MP_BASE(3)
MP_FIVE = MP_BASE(5)

if MODE == 'gmpy':
    int_types = (int, long, MP_BASE_TYPE)
else:
    int_types = (int, long)


class Context(object):

    def __repr__(self):
        lines = ["Mpmath settings:",
            ("  mp.prec = %s" % self.prec).ljust(30) + "[default: 53]",
            ("  mp.dps = %s" % self.dps).ljust(30) + "[default: 15]",
            ("  mp.trap_complex = %s" % self.trap_complex).ljust(30) + "[default: False]",
        ]
        return "\n".join(lines)

    def default(self):
        self._prec = prec_rounding[0] = 53
        self._dps = 15
        self.trap_complex = False

    def set_prec(self, n):
        self._prec = prec_rounding[0] = max(1, int(n))
        self._dps = prec_to_dps(n)

    def set_dps(self, n):
        self._prec = prec_rounding[0] = dps_to_prec(n)
        self._dps = max(1, int(n))

    prec = property(lambda self: self._prec, set_prec)
    dps = property(lambda self: self._dps, set_dps)

# Hack for fast access
prec_rounding = [53, round_nearest]

mp = Context()
mp.default()



class PrecisionManager:

    def __init__(self, precfun, dpsfun, normalize_output=False):
        self.precfun = precfun
        self.dpsfun = dpsfun
        self.normalize_output = normalize_output

    def __call__(self, f):
        def g(*args, **kwargs):
            orig = mp.prec
            try:
                if self.precfun:
                    mp.prec = self.precfun(mp.prec)
                else:
                    mp.dps = self.dpsfun(mp.dps)
                if self.normalize_output:
                    v = f(*args, **kwargs)
                    if type(v) is tuple:
                        return tuple([+a for a in v])
                    return +v
                else:
                    return f(*args, **kwargs)
            finally:
                mp.prec = orig
        g.__name__ = f.__name__
        g.__doc__ = f.__doc__
        return g

    def __enter__(self):
        self.origp = mp.prec
        if self.precfun:
            mp.prec = self.precfun(mp.prec)
        else:
            mp.dps = self.dpsfun(mp.dps)

    def __exit__(self, exc_type, exc_val, exc_tb):
        mp.prec = self.origp
        return False

def extraprec(n, normalize_output=False):
    """
    The block

        with extraprec(n):
            <code>

    increases the precision n bits, executes <code>, and then
    restores the precision.

    extraprec(n)(f) returns a decorated version of the function f
    that increases the working precision by n bits before execution,
    and restores the parent precision afterwards. With
    normalize_output=True, it rounds the return value to the parent
    precision.
    """
    return PrecisionManager(lambda p: p + n, None, normalize_output)

def extradps(n, normalize_output=False):
    """
    This function is analogous to extraprec (see documentation)
    but changes the decimal precision instead of the number of bits.
    """
    return PrecisionManager(None, lambda d: d + n, normalize_output)

def workprec(n, normalize_output=False):
    """
    The block

        with workprec(n):
            <code>

    sets the precision to n bits, executes <code>, and then restores
    the precision.

    workprec(n)(f) returns a decorated version of the function f
    that sets the precision to n bits before execution,
    and restores the precision afterwards. With normalize_output=True,
    it rounds the return value to the parent precision.
    """
    return PrecisionManager(lambda p: n, None, normalize_output)

def workdps(n, normalize_output=False):
    """
    This function is analogous to workprec (see documentation)
    but changes the decimal precision instead of the number of bits.
    """
    return PrecisionManager(None, lambda d: n, normalize_output)

