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

# from django.forms.fields              import IntegerField
# from sympycharfield.validators        import SympyDimensionValidator
# 
# class SympyDimensionCharField(CharField, IntegerField):
#     description = "Sympy Char Field with dimensional verification"
#     
#     def __init__(self, *, dimension=None, **kwargs):
#         if kwargs.get('min_value'):
#             kwargs['min_value'] = str2sympy(kwargs.get('min_value'))
#         if kwargs.get('max_value'):
#             kwargs['max_value'] = str2sympy(kwargs.get('max_value'))
#         super().__init__(**kwargs)
#         
#         self.validators.append(SympyValidator())
#         
#         self.dimension = dimension
#         if self.dimension:
#             self.validators.append(SympyDimensionValidator(dimension))
#     
#     def to_python(self, value):
#         """Return a Sympy object."""
#         
#         if value in self.empty_values:
#             return None
#         value = str(value)
#         if self.strip:
#             value = value.strip()
#         if self.localize:
#             value = formats.sanitize_separators(value)
#         try:
#             value = str2sympy(value)
#         except (ValueError, TypeError):
#             raise ValidationError(self.error_messages['invalid'], code='invalid')
#         return value
#     
#     def widget_attrs(self, widget):
#         attrs = super().widget_attrs(widget)
#         if isinstance(widget, TextInput):
#             if self.min_value is not None:
#                 attrs['min'] = self.min_value
#             if self.max_value is not None:
#                 attrs['max'] = self.max_value
#         return attrs
#     def clean(self, value):
#         return Field.clean(self, value)
#     def validate(self, value):
#         super().validate(value)
