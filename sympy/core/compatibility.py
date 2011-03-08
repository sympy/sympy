"""
Reimplementations of constructs introduced in Python 2.5 for compatibility with
Python 2.4
"""
#XXX: When we drop Python 2.4 support, replace minkey, iff, all, and any
# with their builtin equivalents.
def minkey(sequence, key):
    """
    Implements the min() function with the key argument.

    This is because the key argument isn't supported in Python 2.4.  The
    argument works the same as the key argument to sorted. `key` should be a
    function on the elements of `sequence` that returns objects that are
    comparable with each other, such as integers.  If more than one element
    of the sequence is the smallest with respect to key, then the first one
    will be returned.  This only supports the sequence version of min().

    == Example ==

    >>> from sympy.utilities.iterables import minkey
    >>> minkey(['ab', 'a', 'abc', 'x'], key=len)
    'a'

    """
    smallest = sequence[0]
    smallestkey = key(sequence[0])
    for i in xrange(1, len(sequence)):
        keyi = key(sequence[i])
        if keyi < smallestkey:
            smallest = sequence[i]
            smallestkey = keyi
    return smallest

def iff(condition, result1, result2):
    """
    Return result1 if condition else result2

    This is a replacement for the conditional if statement that is part of
    python 2.5+. If the condition must should not be called unless the
    condition is met, then wrap the result in a lambda; it will be called
    to return the result:

    >>> from sympy import iff
    >>> x = 0.5
    >>> iff(x == 0, x, lambda: 1/x)
    2.0
    >>> x = 0
    >>> iff(x == 0, x, lambda: 1/x)
    0
    """

    if condition:
        rv = result1
    else:
        rv = result2
    # XXX this is fragile; is there a better way to tell if it's a lambda?
    if '<lambda>' in str(rv):
        return rv()
    else:
        return rv

try:
    from __builtins__ import all
except ImportError:
    def all(iterable):
        """
        Return True if all elements are set to True. This
        function does not support predicates explicitly,
        but this behavior can be simulated easily using
        list comprehension.

        >>> from sympy import all
        >>> all( [True, True, True] )
        True
        >>> all( [True, False, True] )
        False
        >>> all( [ x % 2 == 0 for x in [2, 6, 8] ] )
        True
        >>> all( [ x % 2 == 0 for x in [2, 6, 7] ] )
        False

        NOTE: Starting from Python 2.5 this a built-in.
        """
        for item in iterable:
            if not item:
                return False
        return True

try:
    from __builtins__ import all
except ImportError:
    def any(iterable):
        """
        Return True if at least one element is set to True.
        This function does not support predicates explicitly,
        but this behavior can be simulated easily using
        list comprehension.

        >>> from sympy import any
        >>> any( [False, False, False] )
        False
        >>> any( [False, True, False] )
        True
        >>> any( [ x % 2 == 1 for x in [2, 6, 8] ] )
        False
        >>> any( [ x % 2 == 1 for x in [2, 6, 7] ] )
        True

        NOTE: Starting from Python 2.5 this a built-in.
        """
        for item in iterable:
            if item:
                return True
        return False

