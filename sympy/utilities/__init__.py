"""This module contains some general purpose utilities that are used across
SymPy.
"""
from .decorator import memoize_property, public, threaded, xthreaded
from .iterables import capture, cartes, default_sort_key, dict_merge, \
    flatten, group, has_dups, has_variety, interactive_traversal, \
    numbered_symbols, ordered, postfixes, postorder_traversal, prefixes, \
    reshape, sift, subsets, take, topological_sort, unflatten, variations
from .lambdify import lambdify
from .misc import filldedent
from .runtests import doctest, test
from .source import source
from .timeutils import timed
