from sympy.core.traversal import use as _use
from sympy.utilities.decorator import deprecated

use = deprecated(
    useinstead="sympy.core.traversal.use",
    deprecated_since_version="1.10", issue=22288)(_use)
