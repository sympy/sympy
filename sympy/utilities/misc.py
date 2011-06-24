"""Miscellaneous stuff that doesn't really fit anywhere else."""

from sympy.core import sympify

def default_sort_key(item):
    """
    A default sort key for lists of SymPy objects to pass to functions like sorted().

    This uses the default ordering. If you want a nonstandard ordering, you will
    have to create your own sort key using the sort_key() method of the object.

    **Examples**

    >>> from sympy import Basic, S, I, default_sort_key
    >>> from sympy.abc import x

    >>> sorted([S(1)/2, I, -I], key=default_sort_key)
    [1/2, -I, I]
    >>> a = [S(1)/2, I, -I]
    >>> a.sort(key=default_sort_key)
    >>> a
    [1/2, -I, I]

    >>> b = S("[x, 1/x, 1/x**2, x**2, x**(1/2), x**(1/4), x**(3/2)]")
    >>> b.sort(key=default_sort_key)

    The built-in functions min() and max() also take a key function (in Python
    2.5 or higher), that this can be used for.
    """

    #XXX: The following should also be in the docstring, but orders do not
    # actually work at the moment.

    # To use a nonstandard order, you must create your own sort key.  The default
    # order is lex.

    # >>> from sympy import sympify
    # >>> mykey = lambda item: sympify(item).sort_key(order='rev-lex')
    # >>> sorted([x, x**2, 1], key=default_sort_key)
    # [x**2, x, 1]
    # >>> sorted([x, x**2, 1], key=mykey)
    # [1, x, x**2]

    return sympify(item).sort_key()

import sys
size = getattr(sys, "maxint", None)
if size is None: #Python 3 doesn't have maxint
    size = sys.maxsize
if size > 2**32:
    ARCH = "64-bit"
else:
    ARCH = "32-bit"
