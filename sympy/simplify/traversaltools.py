from sympy.core.decorators import deprecated
from sympy.core.traversal import use as _use

use = deprecated(
    useinstead="sympy.core.traversal.use",
    deprecated_since_version="1.10", issue=22288)(_use)
