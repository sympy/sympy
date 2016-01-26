"""Miscellaneous stuff that doesn't really fit anywhere else."""

from __future__ import print_function, division

import sys
import os
from textwrap import fill, dedent
from sympy.core.compatibility import get_function_name, range


def filldedent(s, w=70):
    """
    Strips leading and trailing empty lines from a copy of `s`, then dedents,
    fills and returns it.

    Empty line stripping serves to deal with docstrings like this one that
    start with a newline after the initial triple quote, inserting an empty
    line at the beginning of the string."""
    return '\n' + fill(dedent(str(s)).strip('\n'), width=w)


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
    >>> print(rawlines(s)) # the \\ appears as \ when printed
    (
        'a bee\\n'
        '-----\\n'
        '1 10 '
    )
    >>> print(rawlines('''this
    ... that'''))
    dedent('''\\
        this
        that''')

    >>> print(rawlines('''this
    ... that
    ... '''))
    dedent('''\\
        this
        that
        ''')

    >>> s = \"\"\"this
    ... is a triple '''
    ... \"\"\"
    >>> print(rawlines(s))
    dedent(\"\"\"\\
        this
        is a triple '''
        \"\"\")

    >>> print(rawlines('''this
    ... that
    ...     '''))
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
        last = len(lines) - 1
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

size = getattr(sys, "maxint", None)
if size is None:  # Python 3 doesn't have maxint
    size = sys.maxsize
if size > 2**32:
    ARCH = "64-bit"
else:
    ARCH = "32-bit"


# XXX: PyPy doesn't support hash randomization
HASH_RANDOMIZATION = getattr(sys.flags, 'hash_randomization', False)

_debug_tmp = []
_debug_iter = 0

def debug_decorator(func):
    """If SYMPY_DEBUG is True, it will print a nice execution tree with
    arguments and results of all decorated functions, else do nothing.
    """
    from sympy import SYMPY_DEBUG

    if not SYMPY_DEBUG:
        return func

    def maketree(f, *args, **kw):
        global _debug_tmp
        global _debug_iter
        oldtmp = _debug_tmp
        _debug_tmp = []
        _debug_iter += 1

        def tree(subtrees):
            def indent(s, type=1):
                x = s.split("\n")
                r = "+-%s\n" % x[0]
                for a in x[1:]:
                    if a == "":
                        continue
                    if type == 1:
                        r += "| %s\n" % a
                    else:
                        r += "  %s\n" % a
                return r
            if len(subtrees) == 0:
                return ""
            f = []
            for a in subtrees[:-1]:
                f.append(indent(a))
            f.append(indent(subtrees[-1], 2))
            return ''.join(f)

        # If there is a bug and the algorithm enters an infinite loop, enable the
        # following lines. It will print the names and parameters of all major functions
        # that are called, *before* they are called
        #from sympy.core.compatibility import reduce
        #print("%s%s %s%s" % (_debug_iter, reduce(lambda x, y: x + y, \
        #    map(lambda x: '-', range(1, 2 + _debug_iter))), get_function_name(f), args))

        r = f(*args, **kw)

        _debug_iter -= 1
        s = "%s%s = %s\n" % (get_function_name(f), args, r)
        if _debug_tmp != []:
            s += tree(_debug_tmp)
        _debug_tmp = oldtmp
        _debug_tmp.append(s)
        if _debug_iter == 0:
            print((_debug_tmp[0]))
            _debug_tmp = []
        return r

    def decorated(*args, **kwargs):
        return maketree(func, *args, **kwargs)

    return decorated


def debug(*args):
    """
    Print ``*args`` if SYMPY_DEBUG is True, else do nothing.
    """
    from sympy import SYMPY_DEBUG
    if SYMPY_DEBUG:
        print(*args, file=sys.stderr)


def find_executable(executable, path=None):
    """Try to find 'executable' in the directories listed in 'path' (a
    string listing directories separated by 'os.pathsep'; defaults to
    os.environ['PATH']).  Returns the complete filename or None if not
    found
    """
    if path is None:
        path = os.environ['PATH']
    paths = path.split(os.pathsep)
    extlist = ['']
    if os.name == 'os2':
        (base, ext) = os.path.splitext(executable)
        # executable files on OS/2 can have an arbitrary extension, but
        # .exe is automatically appended if no dot is present in the name
        if not ext:
            executable = executable + ".exe"
    elif sys.platform == 'win32':
        pathext = os.environ['PATHEXT'].lower().split(os.pathsep)
        (base, ext) = os.path.splitext(executable)
        if ext.lower() not in pathext:
            extlist = pathext
    for ext in extlist:
        execname = executable + ext
        if os.path.isfile(execname):
            return execname
        else:
            for p in paths:
                f = os.path.join(p, execname)
                if os.path.isfile(f):
                    return f
    else:
        return None


def func_name(x):
    '''return function name of `x` (if defined) else the `type(x)`.'''
    return getattr(getattr(x, 'func', x), '__name__', type(x))
