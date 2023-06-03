import os
import sys
from typing import Tuple as tTuple, Type

import mpmath.libmp as mlib

from sympy.external import import_module

__all__ = [
    # GROUND_TYPES is either 'gmpy' or 'python' depending on which is used. If
    # gmpy is installed then it will be used unless the environment variable
    # SYMPY_GROUND_TYPES is set to something other than 'auto', 'gmpy', or
    # 'gmpy2'.
    'GROUND_TYPES',

    # If HAS_GMPY is 0, no supported version of gmpy is available. Otherwise,
    # HAS_GMPY will be 2 for gmpy2 if GROUND_TYPES is 'gmpy'. It used to be
    # possible for HAS_GMPY to be 1 for gmpy but gmpy is no longer supported.
    'HAS_GMPY',

    # SYMPY_INTS is a tuple containing the base types for valid integer types.
    # This is either (int,) or (int, type(mpz(0))) depending on GROUND_TYPES.
    'SYMPY_INTS',

    # MPQ is either gmpy.mpq or the Python equivalent from
    # sympy.external.pythonmpq
    'MPQ',

    # MPZ is either gmpy.mpz or int.
    'MPZ',

    # Either the gmpy or the mpmath function
    'factorial',

    # isqrt from gmpy or mpmath
    'sqrt',

    # gcd from gmpy or math
    'gcd',

    # lcm from gmpy or math
    'lcm',

    # invert from gmpy or pow
    'invert',

    # legendre from gmpy or sympy
    'legendre',

    # jacobi from gmpy or sympy
    'jacobi',
]


#
# SYMPY_GROUND_TYPES can be gmpy, gmpy2, python or auto
#
GROUND_TYPES = os.environ.get('SYMPY_GROUND_TYPES', 'auto').lower()


#
# Try to import gmpy2 by default. If gmpy or gmpy2 is specified in
# SYMPY_GROUND_TYPES then warn if gmpy2 is not found. In all cases there is a
# fallback based on pure Python int and PythonMPQ that should still work fine.
#
if GROUND_TYPES in ('auto', 'gmpy', 'gmpy2'):

    # Actually import gmpy2
    gmpy = import_module('gmpy2', min_module_version='2.0.0',
                module_version_attr='version', module_version_attr_call_args=())

    # Warn if user explicitly asked for gmpy but it isn't available.
    if gmpy is None and GROUND_TYPES in ('gmpy', 'gmpy2'):
        from warnings import warn
        warn("gmpy library is not installed, switching to 'python' ground types")

elif GROUND_TYPES == 'python':

    # The user asked for Python so ignore gmpy2 module.
    gmpy = None

else:

    # Invalid value for SYMPY_GROUND_TYPES. Ignore the gmpy2 module.
    from warnings import warn
    warn("SYMPY_GROUND_TYPES environment variable unrecognised. "
         "Should be 'python', 'auto', 'gmpy', or 'gmpy2'")
    gmpy = None


#
# At this point gmpy will be None if gmpy2 was not successfully imported or if
# the environment variable SYMPY_GROUND_TYPES was set to 'python' (or some
# unrecognised value). The two blocks below define the values exported by this
# module in each case.
#
SYMPY_INTS: tTuple[Type, ...]

if gmpy is not None:

    HAS_GMPY = 2
    GROUND_TYPES = 'gmpy'
    SYMPY_INTS = (int, type(gmpy.mpz(0)))
    MPZ = gmpy.mpz
    MPQ = gmpy.mpq

    factorial = gmpy.fac
    sqrt = gmpy.isqrt
    is_square = gmpy.is_square
    sqrtrem = gmpy.isqrt_rem
    gcd = gmpy.gcd
    lcm = gmpy.lcm
    invert = gmpy.invert
    legendre = gmpy.legendre
    jacobi = gmpy.jacobi

else:
    from .pythonmpq import PythonMPQ
    import math

    HAS_GMPY = 0
    GROUND_TYPES = 'python'
    SYMPY_INTS = (int,)
    MPZ = int
    MPQ = PythonMPQ

    factorial = lambda x: int(mlib.ifac(x))
    sqrt = lambda x: int(mlib.isqrt(x))
    is_square = lambda x: x >= 0 and mlib.sqrtrem(x)[1] == 0
    sqrtrem = lambda x: tuple(int(r) for r in mlib.sqrtrem(x))
    if sys.version_info[:2] >= (3, 9):
        gcd = math.gcd
        lcm = math.lcm
    else:
        # Until python 3.8 is no longer supported
        from functools import reduce
        gcd = lambda *args: reduce(math.gcd, args, 0)

        def lcm(*args):
            if 0 in args:
                return 0
            return reduce(lambda x, y: x*y//math.gcd(x, y), args, 1)

    def invert(x, m):
        """ Return y such that x*y == 1 modulo m.

        Uses ``math.pow`` but reproduces the behaviour of ``gmpy2.invert``
        which raises ZeroDivisionError if no inverse exists.
        """
        try:
            return pow(x, -1, m)
        except ValueError:
            raise ZeroDivisionError("invert() no inverse exists")

    def legendre(x, y):
        """ Return Legendre symbol (x / y).

        Following the implementation of gmpy2,
        the error is raised only when y is an even number.
        """
        if y <= 0 or not y % 2:
            raise ValueError("y should be an odd prime")
        x %= y
        if not x:
            return 0
        if pow(x, (y - 1) // 2, y) == 1:
            return 1
        return -1

    def jacobi(x, y):
        """ Return Jacobi symbol (x / y)."""
        if y <= 0 or not y % 2:
            raise ValueError("y should be an odd positive integer")
        x %= y
        if not x:
            return int(y == 1)
        if y == 1 or x == 1:
            return 1
        if gcd(x, y) != 1:
            return 0
        j = 1
        while x != 0:
            while x % 2 == 0 and x > 0:
                x >>= 1
                if y % 8 in [3, 5]:
                    j = -j
            x, y = y, x
            if x % 4 == y % 4 == 3:
                j = -j
            x %= y
        return j
