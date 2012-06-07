"""Miscellaneous stuff that doesn't really fit anywhere else."""

from textwrap import fill, dedent

# if you use
# filldedent('''
#             the text''')
# a space will be put before the first line because dedent will
# put a \n as the first line and fill replaces \n with spaces
# so we strip off any leading and trailing \n since printed wrapped
# text should not have leading or trailing spaces.
filldedent = lambda s, w=70: '\n' + fill(dedent(str(s)).strip('\n'), width=w)

def default_sort_key(item, order=None):
    """A default sort key for lists of SymPy objects to pass to functions
    like sorted().

    This uses the default ordering. If you want a nonstandard ordering,
    you will have to create your own sort key using the sort_key() method
    of the object.

    Examples
    ========

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

    from sympy.core import S, Basic
    from sympy.core.sympify import sympify, SympifyError
    from sympy.core.compatibility import iterable

    if isinstance(item, Basic):
        return item.sort_key(order=order)

    if iterable(item, exclude=basestring):
        if isinstance(item, dict):
            args = item.items()
        else:
            args = list(item)

        args = [default_sort_key(arg, order=order) for arg in args]

        if isinstance(item, dict):
            args = sorted(args)

        cls_index, args = 10, (len(args), tuple(args))
    else:
        if not isinstance(item, basestring):
            try:
                item = sympify(item)
            except SympifyError:
                pass

        if isinstance(item, Basic):
            return item.sort_key(order=order)

        cls_index, args = 0, (1, (str(item),))

    return (cls_index, 0, item.__class__.__name__), args, S.One.sort_key(), S.One

def rawlines(s):
    """Return a cut-and-pastable string that, when printed, is equivalent
    to the input. The string returned is formatted so it can be indented
    nicely within tests; in some cases it is wrapped in the dedent
    function which has to be imported from textwrap.

    Examples
    ========

    Note: because there are characters in the examples below that need
    to be escaped because they are themselves within a triple quoted
    docstring, expressions below look more complicated than they would
    be if they were printed in an interpreter window.

    >>> from sympy.utilities.misc import rawlines
    >>> from sympy import TableForm
    >>> s = str(TableForm([[1, 10]], headings=(None, ['a', 'bee'])))
    >>> print rawlines(s) # the \\ appears as \ when printed
    (
        'a bee\\n'
        '-----\\n'
        '1 10 '
    )
    >>> print rawlines('''this
    ... that''')
    dedent('''\\
        this
        that''')

    >>> print rawlines('''this
    ... that
    ... ''')
    dedent('''\\
        this
        that
        ''')

    >>> s = \"\"\"this
    ... is a triple '''
    ... \"\"\"
    >>> print rawlines(s)
    dedent(\"\"\"\\
        this
        is a triple '''
        \"\"\")

    >>> print rawlines('''this
    ... that
    ...     ''')
    (
        'this\\n'
        'that\\n'
        '    '
    )
    """
    lines = s.split('\n')
    if len(lines) == 1:
        return repr(lines[0])
    triple = ["'''" in s, '"""' in s]
    if any(li.endswith(' ') for li in lines) or '\\' in s or all(triple):
        rv = ["("]
        # add on the newlines
        trailing = s.endswith('\n')
        n = last = len(lines) - 1
        for i, li in enumerate(lines):
            if i != last or trailing:
                rv.append(repr(li)[:-1] + '\\n\'')
            else:
                rv.append(repr(li))
        return '\n    '.join(rv) + '\n)'
    else:
        rv = '\n    '.join(lines)
        if triple[0]:
            return 'dedent("""\\\n    %s""")' % rv
        else:
            return "dedent('''\\\n    %s''')" % rv

import sys
size = getattr(sys, "maxint", None)
if size is None: #Python 3 doesn't have maxint
    size = sys.maxsize
if size > 2**32:
    ARCH = "64-bit"
else:
    ARCH = "32-bit"

def debug(*args):
    """
    Print ``*args`` if SYMPY_DEBUG is True, else do nothing.
    """
    from sympy import SYMPY_DEBUG
    if SYMPY_DEBUG:
        for a in args:
            print a,
        print
