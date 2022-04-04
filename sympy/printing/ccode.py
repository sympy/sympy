"""
.. deprecated:: 1.7

   ccode.py was deprecated and renamed to c.py. This is a shim file to provide
   backwards compatibility.

"""

from sympy.utilities.exceptions import sympy_deprecation_warning

sympy_deprecation_warning(
    """
    The sympy.printing.ccode submodule is deprecated. It has been renamed to
    sympy.printing.c.
    """,
    deprecated_since_version="1.7",
    active_deprecations_target="deprecated-printing-code-submodules",
)

from .c import (ccode, print_ccode, known_functions_C89, known_functions_C99, # noqa:F401
                reserved_words, reserved_words_c99, get_math_macros,
                C89CodePrinter, C99CodePrinter, C11CodePrinter,
                c_code_printers)
