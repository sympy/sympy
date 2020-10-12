"""
This is a shim file to provide backwards compatibility (fcode.py was renamed
to fortran.py in SymPy 1.7).
"""

from .fortran import fcode, print_fcode, known_functions, FCodePrinter # noqa:F401
