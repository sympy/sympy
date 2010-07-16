"""Implementation of mathematical domains. """

from domain import Domain
from ring import Ring
from field import Field

from simpledomain import SimpleDomain
from compositedomain import CompositeDomain

from integerring import IntegerRing
from rationalfield import RationalField
from realdomain import RealDomain

from pythonintegerring import PythonIntegerRing
from sympyintegerring import SymPyIntegerRing
from gmpyintegerring import GMPYIntegerRing

from pythonrationalfield import PythonRationalField
from sympyrationalfield import SymPyRationalField
from gmpyrationalfield import GMPYRationalField

from sympyrealdomain import SymPyRealDomain
from pythonrealdomain import PythonRealDomain
from mpmathrealdomain import MPmathRealDomain

from pythoncomplexdomain import PythonComplexDomain
from mpmathcomplexdomain import MPmathComplexDomain

from algebraicfield import AlgebraicField

from polynomialring import PolynomialRing
from fractionfield import FractionField

from expressiondomain import ExpressionDomain

ZZ_python = PythonIntegerRing
ZZ_sympy = SymPyIntegerRing
ZZ_gmpy = GMPYIntegerRing

QQ_python = PythonRationalField
QQ_sympy = SymPyRationalField
QQ_gmpy = GMPYRationalField

RR_sympy = SymPyRealDomain
RR_python = PythonRealDomain
RR_mpmath = MPmathRealDomain

CC_python = PythonComplexDomain
CC_mpmath = MPmathComplexDomain

from groundtypes import HAS_FRACTION, HAS_GMPY

def _getenv(key, default=None):
    from os import getenv
    return getenv(key, default)

GROUND_TYPES = _getenv('SYMPY_GROUND_TYPES', 'gmpy').lower()

if GROUND_TYPES == 'python':  # XXX: needs 2.6 or better (at least for now)
    ZZ = ZZ_python()

    if HAS_FRACTION:
        QQ = QQ_python()
    elif HAS_GMPY:
        QQ = QQ_gmpy()
        GROUND_TYPES = 'python/gmpy'
    else:
        QQ = QQ_sympy()
        GROUND_TYPES = 'python/sympy'
elif GROUND_TYPES == 'sympy': # XXX: this is *very* slow, guess why ;)
    ZZ = ZZ_sympy()
    QQ = QQ_sympy()
elif GROUND_TYPES == 'gmpy':  # XXX: should be fine? sorry, but no, try -Qnew, damn
    if HAS_GMPY:
        ZZ = ZZ_gmpy()
        QQ = QQ_gmpy()
    else:
        ZZ = ZZ_python()

        if HAS_FRACTION:
            QQ = QQ_python()
            GROUND_TYPES = 'python'
        else:
            QQ = QQ_sympy()
            GROUND_TYPES = 'python/sympy'
else:
    raise ValueError("invalid ground types: %s" % GROUND_TYPES)

RR = RR_mpmath()
CC = CC_mpmath()

EX = ExpressionDomain()
