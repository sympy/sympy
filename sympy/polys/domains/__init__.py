"""Implementation of mathematical domains. """

from domain import Domain
from ring import Ring
from field import Field

from simpledomain import SimpleDomain
from compositedomain import CompositeDomain
from characteristiczero import CharacteristicZero

from finitefield import FiniteField
from integerring import IntegerRing
from rationalfield import RationalField
from realdomain import RealDomain

from pythonfinitefield import PythonFiniteField
from gmpyfinitefield import GMPYFiniteField

from pythonintegerring import PythonIntegerRing
from gmpyintegerring import GMPYIntegerRing

from pythonrationalfield import PythonRationalField
from gmpyrationalfield import GMPYRationalField

from mpmathrealdomain import MPmathRealDomain

from algebraicfield import AlgebraicField

from old_polynomialring import PolynomialRing
from old_fractionfield import FractionField

from expressiondomain import ExpressionDomain

from quotientring import QuotientRing

FF_python = PythonFiniteField
FF_gmpy = GMPYFiniteField

ZZ_python = PythonIntegerRing
ZZ_gmpy = GMPYIntegerRing

QQ_python = PythonRationalField
QQ_gmpy = GMPYRationalField

RR_mpmath = MPmathRealDomain

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

RR = RR_mpmath()

EX = ExpressionDomain()
