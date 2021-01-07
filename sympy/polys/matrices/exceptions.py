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


__all__ = [
    'DDMError', 'DDMShapeError', 'DDMDomainError',

    'NonSquareMatrixError', 'NonInvertibleMatrixError', 'ShapeError',
]
