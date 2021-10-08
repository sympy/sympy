"""

Module to define exceptions to be used in sympy.polys.matrices modules and
classes.

Ideally all exceptions raised in these modules would be defined and documented
here and not e.g. imported from matrices. Also ideally generic exceptions like
ValueError/TypeError would not be raised anywhere.

"""

from sympy.matrices.common import (NonInvertibleMatrixError,
    NonSquareMatrixError, ShapeError)


class DMError(Exception):
    """Base class for errors raised by DDM"""
    pass


class DMBadInputError(DMError):
    """list of lists is inconsistent with shape"""
    pass


class DMDomainError(DMError):
    """domains do not match"""
    pass


class DMFormatError(DMError):
    """mixed dense/sparse not supported"""
    pass


class DMRankError(DMError):
    """matrix does not have expected rank"""
    pass


class DMShapeError(DMError):
    """shapes are inconsistent"""
    pass


__all__ = [
    'DMError', 'DMBadInputError', 'DMDomainError', 'DMFormatError',
    'DMRankError', 'DMShapeError',

    'NonSquareMatrixError', 'NonInvertibleMatrixError', 'ShapeError',
]
