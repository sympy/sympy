from sympy.utilities.exceptions import SymPyDeprecationWarning

SymPyDeprecationWarning(
    feature="Import sympy.utilities.tmpfiles",
    useinstead="Import from sympy.testing.tmpfiles",
    issue=18095,
    deprecated_since_version="1.6").warn()

from sympy.testing.tmpfiles import *  # noqa:F401
