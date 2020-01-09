from sympy.utilities.exceptions import SymPyDeprecationWarning

SymPyDeprecationWarning(
    feature="Import sympy.utilities.quality_unicode",
    useinstead="Import from sympy.testing.quality_unicode",
    issue=18095,
    deprecated_since_version="1.6").warn()

from sympy.testing.quality_unicode import *  # noqa:F401
