"""Miscellaneous stuff that doesn't really fit anywhere else."""

from textwrap import fill, dedent

# if you use
# filldedent('''
#             the text''')
# a space will be put before the first line because dedent will
# put a \n as the first line and fill replaces \n with spaces
# so we strip off any leading and trailing \n since printed wrapped
# text should not have leading or trailing spaces.
filldedent = lambda s: '\n' + fill(dedent(s).strip('\n'))

def default_sort_key(item, order=None):
    """
    A default sort key for lists of SymPy objects to pass to functions like sorted().

    This uses the default ordering. If you want a nonstandard ordering, you will
    have to create your own sort key using the sort_key() method of the object.

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

from functools import wraps
from inspect import getargspec

def __unchanged(x):
    return x

def normalize_args(*normalizers):
    """A decorator that representationally normalizes function arguments.

    "Representational normalization" means converting the value without
    making it unequal to the original.
    This is not useful for operations like rounding an input value to
    the nearest integer (since that changes the value), but it is useful
    for accepting 2.0 (a float) where the function requires an integer.

    Parameters to normalize_args must be "normalizer functions" which
    take a single parameter and return a normalized result.
    normalize_args will then check the results for equality with the
    original value.

    normalize_args will take great care to construct an exception message
    that is useful to end users: it will contain the name of the parameter,
    the name of the normalizer function, the input value, and the unequal
    result value.

    If no normalizer or None is provided for a parameter, the argument is
    passed through unchanged.
    In particular, you will want to use None for the first (self) argument
    of a class function.

    **Purpose**

    normalize_args is useful to narrow the set of potential values down to
    what the decorated function actually wants to handle.

    For example, with

    >>> from sympy.utilities.misc import normalize_args
    >>> @normalize_args(int)
    ... def some_function(n):
    ...     print type(n), n

    some_function() will, for N, never see anything other than an ordinary
    Python integer (at least for those objects that follow the Python
    convention of having an int() function that either returns a Python
    integer or fails):

    >>> some_function(42)
    <type 'int'> 42
    >>> some_function(42.0)
    <type 'int'> 42

    some_function() also does not need to worry about whether the caller
    passed in an int; 2 will do just fine, as will 2.0, or a symbolic
    SymPy expression that evalutes to 2 or 2.0.

    In the code above, the conditions on the argument value are:
    - Since the decorator specified int, it must be possible to call int()
      for it.
    - The result must compare equal to the original value.

    **Usage tips**

    A lambda will work as a normalizer, but it will be reported as
    "<lambda>" in the exception text, which is usually not helpful.
    If you need an ad-hoc normalizer, just wrap it in a function;
    if you make that function public, this will also give your callers
    the ability to inspect the normalizer's result directly.

    Most of the time, you will want to convert stuff to a standard
    data type. In a SymPy context, this usually means int().

    **Known traps**

    "Comparing equal" is not the same as "being equal".
    Dictionary key comparison is one major area where this can bite.

    Since Python ints have unlimited precision, normalizing using
    ``Integer`` is not necessary to do arithmetic, ``int`` is enough.
    Besides, Python ints "are equal" as dictionary keys; Integers are
    not always (in particular, they "are not equal" to Python ints).

    Avoid ``float``. Python ints have unlimited precision, Python
    floats do not, so normalizing via ``float`` will fail with a
    "does not compare equal" diagnostic if the int is too large to
    be represented as a float.
    The language manual explicitly states this on Python floats:
    "You are at the mercy of the underlying machine architecture [...]
    for the accepted range [...]."

    Equality for irrational numbers is undecidable, so normalizing
    real numbers is not a sensible thing to try.

    Converting rationals to float is possible but will work as
    expected only if the denominator is a power of 2. In all other
    cases, rounding errors will make the normalization fail due to
    inequality (which is by design).
    """
    def wrap(f):

        # Get all relevant metadata about f.
        param_names = getargspec(f)[0]

        # Build the paramNormalizers dictionary.
        # It maps parameter names to pertinent normalizer functions.
        # If no normalizer is given or the normalizer is None,
        # the parameter name is mapped to the __unchanged function.
        # We start with mapping all parameters to __unchanged:
        param_normalizers = dict.fromkeys(param_names, __unchanged);
        # Then fill normalizers into their proper slots:
        for param_name, normalizer in zip(param_names, normalizers):
            if normalizer != None:
                param_normalizers[param_name] = normalizer

        # Now we can define the wrapped function:
        @wraps(f)
        def wrapped_f(*args, **kwargs):

            # Merge args into kwargs
            for name, value in zip(param_names, args):
                if name in kwargs:
                    # Generate the same error as
                    # getcallargs in Python 2.7
                    raise TypeError(
                        "%(funcname)s() got multiple values "
                        "for keyword argument '%(paramname)s'"
                        % {
                            "funcname": f.func_name,
                            "paramname": name
                        })
                kwargs[name] = value

            # Generate a dictionary of parameter names with
            # normalized argument values
            normalized_args = {}
            for name, value in kwargs.iteritems ():
                normalized_value = param_normalizers[name](value)
                normalized_args[name] = normalized_value
            # Check that normalization returned values that compare equal.
            if kwargs == normalized_args:
                # call f with normalized args
                return f(**normalized_args)

            # We reach here only if at least one argument turned out
            # unequal with its normalized equivalent.
            # We now pick any argument where this happened and report that
            # one; in theory, we could try and report all of them, but that
            # would generate overly long error messages.
            for name, value in kwargs.iteritems():
                normalizer = param_normalizers[name]
                normalized_value = normalizer(value)
                if value != normalized_value:
                    raise ValueError(
                        "Invalid value for parameter %(paramname)s: "
                        "%(funcname)s(%(value)s) yielded %(result)s, "
                        "which is not equal to %(value)s"
                        % {
                            "funcname": f.func_name,
                            "paramname": name,
                            "value": value,
                            "result": normalized_value
                        }
                    )
            raise ValueError(
                "INTERNAL ERROR: "
                "Dicts are different but all key-value pairs match."
            )

        return wrapped_f
    return wrap
