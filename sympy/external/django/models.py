#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from sympy.external import import_module
django = import_module('django')

from django.db.models                 import CharField, Field
from django.utils.translation         import gettext as _
from sympy.external.django            import str2sympy
from sympy.external.django.validators import SympyValidator

class SympyCharField(CharField):
    description = _("Sympy Char Field")

    def __init__(self, *args, db_collation=None, **kwargs):
        Field.__init__(self, *args, **kwargs)
        self.db_collation = db_collation
        # It is not possible to use MaxLengthValidator with Sympy
        # because validators are run together and MaxLengthValidator
        # needs a string and SympyValidator a Sympy object
        # Sympy object doesn't len() function
        self.validators.append(SympyValidator())

    def from_db_value(self, value, expression, connection, context=None):
        return str2sympy(value)

    def to_python(self, value):
        return str2sympy(super().to_python(value))

    def get_prep_value(self, value):
        value = super().get_prep_value(value)
        return str(value)
