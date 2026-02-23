"""
Utility functions for safe Float comparison in SymPy.

SymPy's Float type uses mpmath binary floating-point internally, so
``Float('0.1') + Float('0.2') != Float('0.3')`` due to representation
error. These utilities provide tolerance-based comparison.
"""

from sympy.core.numbers import Float
from sympy.core.singleton import S


def float_eq(a, b, tol=1e-15):
    """Test approximate equality of two SymPy Float values.

    Parameters
    ==========
    a, b : Float or numeric
        Values to compare.
    tol : float, optional
        Absolute tolerance for comparison. Default is 1e-15.

    Returns
    =======
    bool
        True if ``|a - b| <= tol``.

    Examples
    ========

    >>> from sympy import Float
    >>> from sympy.core.float_utils import float_eq
    >>> float_eq(Float('0.1') + Float('0.2'), Float('0.3'))
    True
    >>> float_eq(Float('1.0'), Float('2.0'))
    False
    """
    from sympy.core.sympify import _sympify
    a = _sympify(a)
    b = _sympify(b)
    diff = a - b
    # Use evalf to get a numeric value for comparison
    try:
        diff_val = float(diff.evalf())
    except (TypeError, ValueError, AttributeError):
        return False
    return abs(diff_val) <= tol


def float_close(a, b, rel_tol=1e-9, abs_tol=1e-15):
    """Test approximate equality using both relative and absolute tolerance.

    Similar to ``math.isclose`` but works with SymPy Float objects.

    Parameters
    ==========
    a, b : Float or numeric
        Values to compare.
    rel_tol : float, optional
        Relative tolerance. Default is 1e-9.
    abs_tol : float, optional
        Absolute tolerance. Default is 1e-15.

    Returns
    =======
    bool
        True if values are close within specified tolerances.

    Examples
    ========

    >>> from sympy import Float
    >>> from sympy.core.float_utils import float_close
    >>> float_close(Float('0.1') + Float('0.2'), Float('0.3'))
    True
    """
    from sympy.core.sympify import _sympify
    a = _sympify(a)
    b = _sympify(b)
    try:
        a_val = float(a.evalf())
        b_val = float(b.evalf())
    except (TypeError, ValueError, AttributeError):
        return False
    diff = abs(a_val - b_val)
    return diff <= max(rel_tol * max(abs(a_val), abs(b_val)), abs_tol)
