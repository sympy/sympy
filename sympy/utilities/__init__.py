"""Some utilities that may help.
"""
import sys

from iterables import make_list, flatten, subsets, numbered_symbols

if sys.version_info[0] <= 2 and sys.version_info[1] < 5:
    from iterables import any, all
else:
    any = any
    all = all

from lambdify import lambdify
from source import source

from decorator import threaded

from runtests import test, doctest

from pytest import raises

