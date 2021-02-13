"""
This is a shim file to provide backwards compatibility (ccode.py was renamed
to c.py in SymPy 1.7).
"""

from sympy.utilities.exceptions import SymPyDeprecationWarning

SymPyDeprecationWarning(
    feature="importing from sympy.printing.ccode",
    useinstead="Import from sympy.printing.c",
    issue=20256,
    deprecated_since_version="1.7").warn()

from .c import (ccode, print_ccode, known_functions_C89, known_functions_C99, # noqa:F401
                reserved_words, reserved_words_c99, get_math_macros,
                C89CodePrinter, C99CodePrinter, C11CodePrinter,
                c_code_printers)
