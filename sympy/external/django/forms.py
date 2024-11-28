#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from sympy.external import import_module
django = import_module('django')

from django.core.validators           import MaxValueValidator, MinValueValidator
from django.forms                     import TextInput, CharField, Field
from django.utils.translation         import gettext as _
from sympy.external.django            import str2sympy
from sympy.external.django.validators import SympyValidator

class SympyCharField(CharField):
    description = _("Sympy Char Field")
    widget      = TextInput

    def __init__(self, *, max_value=None, min_value=None, max_length=None,
                 min_length=None, strip=True, empty_value=0, **kwargs):
        # It is not possible to use MaxLengthValidator with Sympy
        # because validators are run together and MaxLengthValidator
        # needs a string and SympyValidator a Sympy object
        # Sympy object doesn't len() function
        self.max_length  = max_length
        self.min_length  = min_length
        self.strip       = strip
        self.empty_value = empty_value
        Field.__init__(self, **kwargs)

        self.validators.append(SympyValidator())
        if min_value is not None:
            self.validators.append(MinValueValidator(str2sympy(min_value)))
        if max_value is not None:
            self.validators.append(MaxValueValidator(str2sympy(max_value)))

    def to_python(self, value):
        """Return a Sympy object."""
        return str2sympy(super().to_python(value))
