"""Implementation of mathematical domains. """

__all__ = []

from sympy.core.compatibility import GROUND_TYPES

from . import algebraicfield, complexfield, domain, expressiondomain, \
    finitefield, fractionfield, gmpyfinitefield, gmpyintegerring, \
    gmpyrationalfield, integerring, polynomialring, pythonfinitefield, \
    pythonintegerring, pythonrationalfield, rationalfield, realfield
from .algebraicfield import *
from .complexfield import *
from .domain import *
from .expressiondomain import *
from .finitefield import *
from .fractionfield import *
from .gmpyfinitefield import *
from .gmpyintegerring import *
from .gmpyrationalfield import *
from .integerring import *
from .polynomialring import *
from .pythonfinitefield import *
from .pythonintegerring import *
from .pythonrational import PythonRational
from .pythonrationalfield import *
from .rationalfield import *
from .realfield import *

__all__.extend(domain.__all__)

__all__.extend(finitefield.__all__)

__all__.extend(integerring.__all__)

__all__.extend(rationalfield.__all__)

__all__.extend(realfield.__all__)

__all__.extend(complexfield.__all__)

__all__.extend(pythonfinitefield.__all__)

__all__.extend(gmpyfinitefield.__all__)

__all__.extend(pythonintegerring.__all__)

__all__.extend(gmpyintegerring.__all__)

__all__.extend(pythonrationalfield.__all__)

__all__.extend(gmpyrationalfield.__all__)

__all__.extend(algebraicfield.__all__)

__all__.extend(polynomialring.__all__)

__all__.extend(fractionfield.__all__)

__all__.extend(expressiondomain.__all__)

FF_python = PythonFiniteField
FF_gmpy = GMPYFiniteField

ZZ_python = PythonIntegerRing
ZZ_gmpy = GMPYIntegerRing

QQ_python = PythonRationalField
QQ_gmpy = GMPYRationalField

RR = RealField()
CC = ComplexField()



_GROUND_TYPES_MAP = {
    'gmpy': (FF_gmpy, ZZ_gmpy(), QQ_gmpy()),
    'python': (FF_python, ZZ_python(), QQ_python()),
}

try:
    FF, ZZ, QQ = _GROUND_TYPES_MAP[GROUND_TYPES]
except KeyError:
    raise ValueError("invalid ground types: %s" % GROUND_TYPES)

GF = FF

EX = ExpressionDomain()

__all__.extend([
    "FF_python", "FF_gmpy",
    "ZZ_python", "ZZ_gmpy",
    "QQ_python", "QQ_gmpy",
    "GF", "FF", "ZZ", "QQ", "RR", "CC", "EX",
])
