# -*- coding: utf-8 -*-
from sympy.core.decorators import deprecated
from .sympycode import SymPyPrinter, sympy_code, print_sympy_code

@deprecated(
    last_supported_version='1.2',
    useinstead="SymPyPrinter",
    issue=12807,
    deprecated_since_version='1.1')
class PythonPrinter(SymPyPrinter):
    pass


@deprecated(
    last_supported_version='1.2',
    useinstead="SymPyPrinter",
    issue=12807,
    deprecated_since_version='1.1')
def python(expr, **settings):
    return sympy_code(expr, **settings)


@deprecated(
    last_supported_version='1.2',
    useinstead="SymPyPrinter",
    issue=12807,
    deprecated_since_version='1.1')
def print_python(expr, **settings):
    return print_sympy_code(expr, **settings)
