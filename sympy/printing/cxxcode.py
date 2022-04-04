"""
.. deprecated:: 1.7

   cxxcode.py was deprecated and renamed to cxx.py. This is a shim file to
   provide backwards compatibility.

"""

from sympy.utilities.exceptions import sympy_deprecation_warning

sympy_deprecation_warning(
    """
    The sympy.printing.cxxcode submodule is deprecated. It has been renamed to
    sympy.printing.cxx.
    """,
    deprecated_since_version="1.7",
    active_deprecations_target="deprecated-printing-code-submodules",
)

from .cxx import (cxxcode, reserved, CXX98CodePrinter, # noqa:F401
                  CXX11CodePrinter, CXX17CodePrinter, cxx_code_printers)
