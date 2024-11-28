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
