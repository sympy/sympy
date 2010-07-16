"""Some utilities that may help.
"""
from iterables import (iff, flatten, group, split, take, subsets,
    variations, numbered_symbols, cartes, capture, any, all, dict_merge)

from lambdify import lambdify
from source import source

from decorator import threaded, deprecated, wraps

from cythonutils import cythonized
from timeutils import timed
