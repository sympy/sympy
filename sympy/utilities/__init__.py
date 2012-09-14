"""This module contains some general purpose utilities that are used across
SymPy.
"""
from iterables import (flatten, group, take, subsets,
    variations, numbered_symbols, cartes, capture, dict_merge,
    postorder_traversal, interactive_traversal,
    prefixes, postfixes, sift, topological_sort, unflatten,
    quick_sort, has_dups, has_variety, reshape)

from lambdify import lambdify
from source import source

from decorator import threaded, xthreaded

from runtests import test, doctest

from cythonutils import cythonized
from timeutils import timed

from misc import default_sort_key
