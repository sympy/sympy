"""This module contains some general purpose utilities that are used across
SymPy.
"""
from __future__ import annotations
from .iterables import (flatten, group, take, subsets,
    variations, numbered_symbols, cartes, capture, dict_merge,
    prefixes, postfixes, sift, topological_sort, unflatten,
    has_dups, has_variety, reshape, rotations)

from .misc import filldedent

from .lambdify import lambdify

from .decorator import threaded, xthreaded, public, memoize_property

from .timeutils import timed

from .optimize_derivatives import ( generate_optimized_derivatives, count_operations )

__all__ = [
    'flatten', 'group', 'take', 'subsets', 'variations', 'numbered_symbols',
    'cartes', 'capture', 'dict_merge', 'prefixes', 'postfixes', 'sift',
    'topological_sort', 'unflatten', 'has_dups', 'has_variety', 'reshape',
    'rotations',

    'filldedent',

    'lambdify',

    'threaded', 'xthreaded', 'public', 'memoize_property',

    'timed',

    'generate_optimized_derivatives', 'count_operations'
]
