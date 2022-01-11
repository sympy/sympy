"""
This is a shim file to provide backwards compatibility (fcode.py was renamed
to fortran.py in SymPy 1.7).
"""

from sympy.utilities.exceptions import SymPyDeprecationWarning

SymPyDeprecationWarning(
    feature="importing from sympy.printing.fcode",
    useinstead="Import from sympy.printing.fortran",
    issue=20256,
    deprecated_since_version="1.7").warn()

from .fortran import fcode, print_fcode, known_functions, FCodePrinter # noqa:F401
