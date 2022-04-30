"""
.. deprecated:: 1.7

   fcode.py was deprecated and renamed to fortran.py. This is a shim file to
   provide backwards compatibility.

"""

from sympy.utilities.exceptions import sympy_deprecation_warning

sympy_deprecation_warning(
    """
    The sympy.printing.fcode submodule is deprecated. It has been renamed to
    sympy.printing.fortran.
    """,
    deprecated_since_version="1.7",
    active_deprecations_target="deprecated-printing-code-submodules",
)

from .fortran import fcode, print_fcode, known_functions, FCodePrinter # noqa:F401
