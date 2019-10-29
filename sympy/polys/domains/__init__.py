"""Implementation of mathematical domains. """

__all__ = ['Domain', 'FiniteField', 'IntegerRing', 'RationalField',
        'RealField', 'ComplexField', 'PythonFiniteField', 'GMPYFiniteField',
        'PythonIntegerRing', 'GMPYIntegerRing', 'PythonRational',
        'GMPYRationalField', 'AlgebraicField', 'PolynomialRing',
        'FractionField', 'ExpressionDomain', 'PythonRational']

from .domain import Domain
from .finitefield import FiniteField
from .integerring import IntegerRing
from .rationalfield import RationalField
from .realfield import RealField
from .complexfield import ComplexField
from .pythonfinitefield import PythonFiniteField
from .gmpyfinitefield import GMPYFiniteField
from .pythonintegerring import PythonIntegerRing
from .gmpyintegerring import GMPYIntegerRing
from .pythonrationalfield import PythonRationalField
from .gmpyrationalfield import GMPYRationalField
from .algebraicfield import AlgebraicField
from .polynomialring import PolynomialRing
from .fractionfield import FractionField
from .expressiondomain import ExpressionDomain
from .pythonrational import PythonRational

FF_python = PythonFiniteField
FF_gmpy = GMPYFiniteField

ZZ_python = PythonIntegerRing
ZZ_gmpy = GMPYIntegerRing

QQ_python = PythonRationalField
QQ_gmpy = GMPYRationalField

RR = RealField()
CC = ComplexField()

from sympy.core.compatibility import GROUND_TYPES

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
