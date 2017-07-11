"""Calculus-related methods."""

from .euler import euler_equations
from .finite_diff import apply_finite_diff, as_finite_diff, \
    differentiate_finite, finite_diff_weights
from .singularities import is_decreasing, is_increasing, is_monotonic, \
    is_strictly_decreasing, is_strictly_increasing, singularities
from .util import AccumBounds, not_empty_in, periodicity

__all__ = [
'euler_equations',

'singularities', 'is_increasing',
'is_strictly_increasing', 'is_decreasing',
'is_strictly_decreasing', 'is_monotonic',

'finite_diff_weights', 'apply_finite_diff', 'as_finite_diff', 'differentiate_finite',

'periodicity', 'not_empty_in', 'AccumBounds',
]
