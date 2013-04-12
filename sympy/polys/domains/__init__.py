"""Implementation of mathematical domains. """

__all__ = []

import domain
__all__.extend(domain.__all__)
from domain import Domain

import finitefield
__all__.extend(finitefield.__all__)
from finitefield import FiniteField

import integerring
__all__.extend(integerring.__all__)
from integerring import IntegerRing

import rationalfield
__all__.extend(rationalfield.__all__)
from rationalfield import RationalField

import realfield
__all__.extend(realfield.__all__)
from realfield import RealField

import complexfield
__all__.extend(complexfield.__all__)
from complexfield import ComplexField

import pythonfinitefield
__all__.extend(pythonfinitefield.__all__)
from pythonfinitefield import PythonFiniteField

import gmpyfinitefield
__all__.extend(gmpyfinitefield.__all__)
from gmpyfinitefield import GMPYFiniteField

import pythonintegerring
__all__.extend(pythonintegerring.__all__)
from pythonintegerring import PythonIntegerRing

import gmpyintegerring
__all__.extend(gmpyintegerring.__all__)
from gmpyintegerring import GMPYIntegerRing

import pythonrationalfield
__all__.extend(pythonrationalfield.__all__)
from pythonrationalfield import PythonRationalField

import gmpyrationalfield
__all__.extend(gmpyrationalfield.__all__)
from gmpyrationalfield import GMPYRationalField

import algebraicfield
__all__.extend(algebraicfield.__all__)
from algebraicfield import AlgebraicField

import polynomialring
__all__.extend(polynomialring.__all__)
from polynomialring import PolynomialRing

import fractionfield
__all__.extend(fractionfield.__all__)
from fractionfield import FractionField

import expressiondomain
__all__.extend(expressiondomain.__all__)
from expressiondomain import ExpressionDomain

FF_python = PythonFiniteField
FF_gmpy = GMPYFiniteField

ZZ_python = PythonIntegerRing
ZZ_gmpy = GMPYIntegerRing

QQ_python = PythonRationalField
QQ_gmpy = GMPYRationalField

RR = RealField()
CC = ComplexField()

from pythonrational import PythonRational

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
