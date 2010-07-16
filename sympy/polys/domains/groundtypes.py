"""Ground types for various mathematical domains in SymPy. """

HAS_FRACTION = True

try:
    import fractions
except ImportError:
    HAS_FRACTION = False

HAS_GMPY = True

try:
    import gmpy
except ImportError:
    HAS_GMPY = False

from __builtin__ import (
    int     as PythonIntegerType,
    float   as PythonRealType,
    complex as PythonComplexType,
)

if HAS_FRACTION:
    from fractions import (
        Fraction as PythonRationalType,
    )
else:
    class PythonRationalType(object):
        def __init__(self, obj):
            pass

from sympy.core.numbers import (
    ifactorial as python_factorial,
    igcdex     as python_gcdex,
    igcd       as python_gcd,
    ilcm       as python_lcm,
)

from sympy import (
    Real     as SymPyRealType,
    Integer  as SymPyIntegerType,
    Rational as SymPyRationalType,
)

from sympy.core.numbers import (
    ifactorial as sympy_factorial,
    igcdex     as sympy_gcdex,
    igcd       as sympy_gcd,
    ilcm       as sympy_lcm,
)

if HAS_GMPY:
    from gmpy import (
        mpz    as GMPYIntegerType,
        mpq    as GMPYRationalType,
        numer  as gmpy_numer,
        denom  as gmpy_denom,
        gcdext as gmpy_gcdex,
        gcd    as gmpy_gcd,
        lcm    as gmpy_lcm,
        sqrt   as gmpy_sqrt,
        fac    as gmpy_factorial,
    )
else:
    class GMPYIntegerType(object):
        def __init__(self, obj):
            pass

    class GMPYRationalType(object):
        def __init__(self, obj):
            pass

    gmpy_factorial   = None
    gmpy_numer       = None
    gmpy_denom       = None
    gmpy_gcdex       = None
    gmpy_gcd         = None
    gmpy_lcm         = None
    gmpy_sqrt        = None
    gmpy_factorial   = None

from sympy.mpmath import (
    mpf as MPmathRealType,
    mpc as MPmathComplexType,
    mpi as MPmathIntervalType,
)

from sympy.mpmath.libmp.libmpf import isqrt

def python_sqrt(a):
    return int(isqrt(a))
