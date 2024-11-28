#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from sympy.external import import_module
django = import_module('django')

from django.core.exceptions   import ValidationError
from django.utils.deconstruct import deconstructible
from django.utils.translation import gettext_lazy as _
from sympy                    import Basic

class SympyValidationError(ValidationError):
    """
    Sympy Validation Error.
    """
    message = _('Enter a Sympy Object.')
    
    def __init__(self, message=message, code=None, params=None):
        super().__init__(message, code, params)

@deconstructible
class SympyValidator:
    """
    Validate if the input is an instance of Sympy, otherwise raise SympyValidationError.
    """
    
    def __init__(self):
        pass
    
    def __call__(self, value):
        if not isinstance(value, Basic):
            raise SympyValidationError()
    
    def __eq__(self, other):
        return (isinstance(other, self.__class__))

# 
# Unfortunately, it is impossible to implement the validator below because
# Sympy did not implement 
# raising the TypeError:
# "Type <class 'sympy.physics.units.quantities.Quantity'> not implemented for get_dimensional_dependencies"
# 
# If you need this dimensional verification use
# Unum instead https://unum.readthedocs.io/en/latest/
# 
# from sympy.physics.units.systems.si import dimsys_SI
# 
# @deconstructible
# class SympyDimensionValidator(SympyValidator):
#     """
#     Validate if the input has the same Sympy dimension, otherwise raise SympyValidationError.
#     """
#     messages = {
#         'invalid': _('Enter a Sympy Object.'),
#         'dimension': _(
#             'Ensure dimension is equal than %(dim).',
#             'Ensure dimension is equal than %(dim).',
#             'dim'
#         ),
#     }
#     
#     def __init__(self, dimension=None, **kwargs):
#         self.dimension = dimension
#         super().__init__(**kwargs)
#     
#     def __call__(self, value):
#         super().__call__(value)
#         if dimsys_SI.equivalent_dims(self.dimension, value):
#             raise ValidationError(
#                 self.messages['dimension'],
#                 code='dimension',
#                 params={'dim': (self.dimension)},
#             )
# 
#     def __eq__(self, other):
#         return (
#             super().__eq__(other) and
#             self.dimension == other.dimension
#         )

