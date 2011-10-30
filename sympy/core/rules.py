"""
Replacement rules.
"""

class Transform(object):
    """
    Generic transformation rule.

    The Transform object is hybrid function/dictionary that is created
    with a ``transform`` function and a ``filter``. When accessed with
    an index Transform()[key] it behaves like a default dictionary (or
    function) and when accessed with get() it behaves like a dictionary.

    >>> from sympy.core.rules import Transform
    >>> from sympy.abc import x

    This Transform will return, as a value, one more than the key:

    >>> add1 = Transform(lambda x: x + 1)
    >>> add1[1]
    2
    >>> add1[x]
    x + 1

    By default, all values are considered to be in the dictionary. If a filter
    is supplied, it will limit the scope of the application of the transform
    when accessing the Transform with the get() method; the filter is ignored
    when accessing it via index notation (and in this sense the Transform is
    like a function or defaultdict):

    >>> add1_odd = Transform(lambda x: x + 1, lambda x: x%2 == 1)
    >>> add1_odd[2]
    3
    >>> add1_odd.get(2) is None
    True
    >>> add1_odd.get(2, 0)
    0

    The filter is also used to answer a query as to whether the Transform
    contains a given value or not:

    >>> 2 in add1_odd
    False
    >>> 3 in add1_odd
    True
    """

    def __init__(self, transform, filter=lambda x: True):
        self.transform = transform
        self.filter = filter

    def __contains__(self, item):
        return self.filter(item)

    def __getitem__(self, key):
        return self.transform(key)

    def get(self, item, default=None):
        if item in self:
            print 'yo'
            return self[item]
        else:
            return default
