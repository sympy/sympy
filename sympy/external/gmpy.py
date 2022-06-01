import os
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

else:
    from .pythonmpq import PythonMPQ

    HAS_GMPY = 0
    GROUND_TYPES = 'python'
    SYMPY_INTS = (int,)
    MPZ = int
    MPQ = PythonMPQ

    factorial = lambda x: int(mlib.ifac(x))
    sqrt = lambda x: int(mlib.isqrt(x))
