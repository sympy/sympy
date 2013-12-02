"""This module contains some general purpose utilities that are used across
SymPy.
"""
from .iterables import (flatten, group, take, subsets,
    variations, numbered_symbols, cartes, capture, dict_merge,
    postorder_traversal, interactive_traversal,
    prefixes, postfixes, sift, topological_sort, unflatten,
    has_dups, has_variety, reshape, default_sort_key, ordered)

from .lambdify import lambdify
from .source import source

from .decorator import threaded, xthreaded, public

from .runtests import test, doctest

from .timeutils import timed

from .solution import add_step, add_comment, add_eq, add_exp, reset_solution, last_solution, start_subroutine, cancel_subroutine, commit_subroutine
