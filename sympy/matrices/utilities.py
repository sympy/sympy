from __future__ import division, print_function

from sympy.core.function import expand_mul


def _toselfclass(self, M):
    return M if isinstance(M, self.__class__) else self.__class__(M)


def _iszero(x):
    """Returns True if x is zero."""
    return getattr(x, 'is_zero', None)


def _is_zero_after_expand_mul(x):
    """Tests by expand_mul only, suitable for polynomials and rational
    functions."""
    return expand_mul(x) == 0


