import os

from sympy.core.function import expand_mul
from sympy.simplify.simplify import dotprodsimp as _dotprodsimp


# The following is an internal variable for controlling the recently introduced
# dotprodsimp intermediate simplification step in matrix operations in one
# place. It is intended as an emergency switch in cases where user code does not
# like the different structure of results that comes from this simplification
# and can not be adapted for some reason. When the intermediate simplification
# step is considered fully compatible with user code and this mechanism is no
# longer needed in can be removed.

# The default value of `None` specifies that dotprodsimp be used in a few
# selected low-level functions but not in others. Setting this global variable
# to `False` will turn off the dotprodsimp intermediate simplifications
# everywhere and setting to `True` will turn it on everywhere in matrices where
# it can be applied.

# There are a few other places in the matrix code where dotprodsimp can probably
# help, these are places where a call is made to:
#
# dps = _get_intermediate_simp()
#
# To determine whether dotprodsimp helps in these places testing needs to be
# done, to turn dotprodsimp on in these places by default replace this call with:
#
# from sympy.simplify.simplify import dotprodsimp as _dotprodsimp
# dps = _get_intermediate_simp(_dotprodsimp)
#
# Or possibly:
#
# from sympy import expand_mul
# dps = _get_intermediate_simp(expand_mul, expand_mul)
#
# This second form uses lighter simplification by default but may still do
# better than nothing.
#
# The other place where dotprodsimp may be added is any place where matrices are
# multiplied via:
#
# A.multiply(B) -> A.multiply(B, dotprodsimp=True)

# True, False or None
_DOTPRODSIMP_MODE = False if os.environ.get('SYMPY_DOTPRODSIMP', '').lower() in \
        ('false', 'off', '0') else None

def _get_intermediate_simp(deffunc=lambda x: x, offfunc=lambda x: x,
        onfunc=_dotprodsimp, dotprodsimp=None):
    """Support function for controlling intermediate simplification. Returns a
    simplification function according to the global setting of dotprodsimp
    operation.

    ``deffunc``     - Function to be used by default.
    ``offfunc``     - Function to be used if dotprodsimp has been turned off.
    ``onfunc``      - Function to be used if dotprodsimp has been turned on.
    ``dotprodsimp`` - True, False or None. Will be overriden by global
                      _DOTPRODSIMP_MODE if that is not None.
    """

    if dotprodsimp is False or _DOTPRODSIMP_MODE is False:
        return offfunc
    if dotprodsimp is True or _DOTPRODSIMP_MODE is True:
        return onfunc

    return deffunc # None, None

def _get_intermediate_simp_bool(default=False, dotprodsimp=None):
    """Same as ``_get_intermediate_simp`` but returns bools instead of functions
    by default."""

    return _get_intermediate_simp(default, False, True, dotprodsimp)


def _iszero(x):
    """Returns True if x is zero."""
    return getattr(x, 'is_zero', None)


def _is_zero_after_expand_mul(x):
    """Tests by expand_mul only, suitable for polynomials and rational
    functions."""
    return expand_mul(x) == 0
