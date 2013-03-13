"""Ground types for various mathematical domains in SymPy. """

from sympy.external import import_module

# The logic for detecting if a compatible version of gmpy/gmpy2 is present is
# done in sympy.core.compatibility.

from sympy.core.compatibility import HAS_GMPY

from __builtin__ import (
    int as PythonIntegerType,
    float as PythonRealType,
    complex as PythonComplexType,
)

from pythonrational import PythonRational

from sympy.core.numbers import (
    igcdex as python_gcdex,
    igcd as python_gcd,
    ilcm as python_lcm,
)

from sympy import (
    Float as SymPyRealType,
    Integer as SymPyIntegerType,
    Rational as SymPyRationalType,
)

if HAS_GMPY == 1:
    from gmpy import (
        mpz as GMPYIntegerType,
        mpq as GMPYRationalType,
        fac as gmpy_factorial,
        numer as gmpy_numer,
        denom as gmpy_denom,
        gcdext as gmpy_gcdex,
        gcd as gmpy_gcd,
        lcm as gmpy_lcm,
        sqrt as gmpy_sqrt,
        qdiv as gmpy_qdiv,
    )
elif HAS_GMPY == 2:
    from gmpy2 import (
        mpz as GMPYIntegerType,
        mpq as GMPYRationalType,
        fac as gmpy_factorial,
        numer as gmpy_numer,
        denom as gmpy_denom,
        gcdext as gmpy_gcdex,
        gcd as gmpy_gcd,
        lcm as gmpy_lcm,
        isqrt as gmpy_sqrt,
        qdiv as gmpy_qdiv,
    )
else:
    class GMPYIntegerType(object):
        def __init__(self, obj):
            pass

    class GMPYRationalType(object):
        def __init__(self, obj):
            pass

    gmpy_factorial = None
    gmpy_numer = None
    gmpy_denom = None
    gmpy_gcdex = None
    gmpy_gcd = None
    gmpy_lcm = None
    gmpy_sqrt = None
    gmpy_qdiv = None

from sympy.mpmath import (
    mpf as MPmathRealType,
    mpc as MPmathComplexType,
    mpi as MPmathIntervalType,
)

import sympy.mpmath.libmp as mlib


def python_sqrt(n):
    return int(mlib.isqrt(n))


def python_factorial(n):
    return int(mlib.ifac(n))
