"""
Replacement rules.
"""

from __future__ import print_function, division


class Transform(object):
    """
    Immutable mapping that can be used as a generic transformation rule.

    Parameters
    ----------
    transform : callable
        Computes the value corresponding to any key.
    filter : callable, optional
        If supplied, specifies which objects are in the mapping.

    Examples
    ========

    >>> from sympy.core.rules import Transform
    >>> from sympy.abc import x

    This Transform will return, as a value, one more than the key:

    >>> add1 = Transform(lambda x: x + 1)
    >>> add1[1]
    2
    >>> add1[x]
    x + 1

    By default, all values are considered to be in the dictionary. If a filter
    is supplied, only the objects for which it returns True are considered as
    being in the dictionary:

    >>> add1_odd = Transform(lambda x: x + 1, lambda x: x%2 == 1)
    >>> 2 in add1_odd
    False
    >>> add1_odd.get(2, 0)
    0
    >>> 3 in add1_odd
    True
    >>> add1_odd[3]
    4
    >>> add1_odd.get(3, 0)
    4
    """

    def __init__(self, transform, filter=lambda x: True):
        self._transform = transform
        self._filter = filter

    def __contains__(self, item):
        return self._filter(item)

    def __getitem__(self, key):
        if self._filter(key):
            return self._transform(key)
        else:
            raise KeyError(key)

    def get(self, item, default=None):
        if item in self:
            return self[item]
        else:
            return default


class MapMatcher(object):
    """
    Dictionary-like data structure to handle multiple pattern matching statements.
    This object is assigned pattern matching expressions and corresponding functions.
    The functions take one parameter, i.e. the dictionary of substitutions.

    Data can be passed on costruction, or later with item assignment.
    It is important to notice that the matched pattern for a given expression
    is the first pattern that matches, in the order given upon construction,
    any following pattern is disregarded.

    Examples
    ========

    >>> from sympy.core.rules import MapMatcher
    >>> from sympy import symbols, Wild, sin, S
    >>> from sympy import Tuple
    >>> x, y = symbols('x y')
    >>> a = Wild('a', exclude=[x, y])
    >>> b = Wild('b', exclude=[x])
    >>> map_matcher = MapMatcher()
    >>> map_matcher[a/b + x] = Tuple(a, b)
    >>> map_matcher[a + y] = Tuple(a, b)
    >>> map_matcher[3/y + x]
    (3, y)
    >>> map_matcher[y + 2]
    (2, b_)

    If no match is successful, a ``KeyError`` exception is raised.

    >>> # map_matcher[S.One]

    """
    def __init__(self, initial_data=()):
        self._pattern_func_list = []
        for i, j in initial_data:
            self[i] = j

    def __setitem__(self, key, value):
        self._pattern_func_list.append( (key, value) )

    def __getitem__(self, key):
        # the first match is returned.
        for pattern, value in self._pattern_func_list:
            d = key.match(pattern)
            if d:
                return value.xreplace(d)
        raise KeyError()
