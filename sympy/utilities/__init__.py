"""Some utilities that may help.
"""
import sys

from iterables import (iff, make_list, flatten, group, split,
    subsets, variations, numbered_symbols, take, dict_merge)

if sys.version_info[0] <= 2 and sys.version_info[1] < 5:
    from iterables import any, all
else:
    any = any
    all = all

from lambdify import lambdify
from source import source

from decorator import threaded, xthreaded, deprecated, wraps

from runtests import test, doctest

from pytest import raises

from cythonutils import cythonized
from timeutils import timed

