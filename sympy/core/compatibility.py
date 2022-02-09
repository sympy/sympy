"""
Reimplementations of constructs introduced in later versions of Python than
we support. Also some functions that are needed SymPy-wide and are located
here for easy import.
"""


from .sorting import ordered as _ordered, _nodes as __nodes, default_sort_key as _default_sort_key
from sympy.utilities.decorator import deprecated
from sympy.utilities.misc import as_int as _as_int
from sympy.utilities.iterables import (iterable as _iterable,
    is_sequence as _is_sequence, NotIterable as _notiterable)


default_sort_key = deprecated(useinstead="sympy.core.sorting.default_sort_key",
    deprecated_since_version="1.10", issue=22352)(_default_sort_key)


ordered = deprecated(useinstead="sympy.core.sorting.ordered",
    deprecated_since_version="1.10", issue=22352)(_ordered)


_nodes = deprecated(useinstead="sympy.core.sorting._nodes",
    deprecated_since_version="1.10", issue=22352)(__nodes)


as_int = deprecated(useinstead="sympy.utilities.misc.as_int",
    deprecated_since_version="1.10", issue=22352)(_as_int)


is_sequence = deprecated(useinstead="sympy.utilities.iterables.is_sequence",
    deprecated_since_version="1.10", issue=22352)(_is_sequence)


iterable = deprecated(useinstead="sympy.utilities.iterables.iterable",
    deprecated_since_version="1.10", issue=22352)(_iterable)

NotIterable = deprecated(useinstead="sympy.utilities.iterables.NotIterable",
    deprecated_since_version="1.10", issue=22352)(_notiterable)
