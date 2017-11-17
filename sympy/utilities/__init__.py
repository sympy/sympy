"""This module contains some general purpose utilities that are used across
SymPy.
"""

__all__ = []

from .iterables import (
    flatten, group, take, subsets,
    variations, numbered_symbols, cartes, capture, dict_merge,
    postorder_traversal, interactive_traversal,
    prefixes, postfixes, sift, topological_sort, unflatten,
    has_dups, has_variety, reshape, default_sort_key, ordered
)
__all__ += [
    "flatten", "group", "take", "subsets",
    "variations", "numbered_symbols", "cartes", "capture", "dict_merge",
    "postorder_traversal", "interactive_traversal",
    "prefixes", "postfixes", "sift", "topological_sort", "unflatten",
    "has_dups", "has_variety", "reshape", "default_sort_key", "ordered"
]

from .misc import filldedent
__all__ += ["filldedent"]

from .lambdify import lambdify
__all__ += ["lambdify"]

from .source import source
__all__ += ["source"]

from .decorator import threaded, xthreaded, public, memoize_property
__all__ += ["threaded", "xthreaded", "public", "memoize_property"]

from .runtests import test, doctest
__all__ += ["test", "doctest"]

from .timeutils import timed
__all__ += ["timed"]
