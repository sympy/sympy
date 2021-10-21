"""This module contains some general purpose utilities that are used across
SymPy.
"""
from .iterables import (flatten, group, take, subsets,
    variations, numbered_symbols, cartes, capture, dict_merge,
    prefixes, postfixes, sift, topological_sort, unflatten,
    has_dups, has_variety, reshape, default_sort_key, ordered,
    rotations)

from .misc import filldedent

from .source import source

from .decorator import threaded, xthreaded, public, memoize_property

from .timeutils import timed

__all__ = [
    'flatten', 'group', 'take', 'subsets', 'variations', 'numbered_symbols',
    'cartes', 'capture', 'dict_merge', 'prefixes', 'postfixes', 'sift',
    'topological_sort', 'unflatten', 'has_dups', 'has_variety', 'reshape',
    'default_sort_key', 'ordered', 'rotations',

    'filldedent',

    'source',

    'threaded', 'xthreaded', 'public', 'memoize_property',

    'timed',
]
