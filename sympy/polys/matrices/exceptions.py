"""

Module to define exceptions to be used in sympy.polys.matrices modules and
classes.

Ideally all exceptions raised in these modules would be defined and documented
here and not e.g. imported from matrices. Also ideally generic exceptions like
ValueError/TypeError would not be raised anywhere.

"""

from sympy.matrices.common import (NonInvertibleMatrixError,
    NonSquareMatrixError, ShapeError)


class DDMError(Exception):
    """Base class for errors raised by DDM"""
    pass


class DDMBadInputError(DDMError):
    """list of lists is inconsistent with shape"""
    pass


class DDMDomainError(DDMError):
    """domains do not match"""
    pass


class DDMShapeError(DDMError):
    """shapes are inconsistent"""
    pass


class DDMFormatError(DDMError):
    """mixed dense/sparse not supported"""
    pass


__all__ = [
    'DDMError', 'DDMShapeError', 'DDMDomainError', 'DDMFormatError',

    'NonSquareMatrixError', 'NonInvertibleMatrixError', 'ShapeError',
]
