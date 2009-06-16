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
sage = None

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

if 'MPMATH_NOSAGE' not in os.environ:
    try:
        import sage.all
        if hasattr(sage.all.Integer, "trailing_zero_bits"):
            sage = sage.all
            MODE = 'sage'
            MP_BASE = sage.Integer
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

