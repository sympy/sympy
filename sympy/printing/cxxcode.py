"""
This is a shim file to provide backwards compatibility (cxxcode.py was renamed
to cxx.py in SymPy 1.7).
"""

from sympy.utilities.exceptions import SymPyDeprecationWarning

SymPyDeprecationWarning(
    feature="importing from sympy.printing.cxxcode",
    useinstead="Import from sympy.printing.cxx",
    issue=20256,
    deprecated_since_version="1.7").warn()

from .cxx import (cxxcode, reserved, CXX98CodePrinter, # noqa:F401
                  CXX11CodePrinter, CXX17CodePrinter, cxx_code_printers)
