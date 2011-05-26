"""
SymPy core decorators.

The purpose of this module is to expose decorators without any other
dependencies, so that they can be easily imported anywhere in sympy/core.
"""
import re
import operator
import warnings

from sympify import SympifyError, sympify

try:
    from functools import wraps
except ImportError:
    def wraps(old_func):
        """Copy private data from ``old_func`` to ``new_func``. """
        def decorate(new_func):
            new_func.__dict__.update(old_func.__dict__)
            new_func.__module__ = old_func.__module__
            new_func.__name__   = old_func.__name__
            new_func.__doc__    = old_func.__doc__
            return new_func
        return decorate

def deprecated(func):
    """This is a decorator which can be used to mark functions
    as deprecated. It will result in a warning being emitted
    when the function is used."""
    @wraps(func)
    def new_func(*args, **kwargs):
        warnings.warn("Call to deprecated function %s." % func.__name__,
                      category=DeprecationWarning)
        return func(*args, **kwargs)
    return new_func

def _sympifyit(arg, retval=None):
    """decorator to smartly _sympify function arguments

       @_sympifyit('other', NotImplemented)
       def add(self, other):
           ...

       In add, other can be thought of as already being a SymPy object.

       If it is not, the code is likely to catch an exception, then other will
       be explicitly _sympified, and the whole code restarted.

       if _sympify(arg) fails, NotImplemented will be returned

       see: __sympifyit
    """
    def deco(func):
        return __sympifyit(func, arg, retval)

    return deco

def __sympifyit(func, arg, retval=None):
    """decorator to _sympify `arg` argument for function `func`

       don't use directly -- use _sympifyit instead
    """

    # we support f(a,b) only
    assert func.func_code.co_argcount
    # only b is _sympified
    assert func.func_code.co_varnames[1] == arg

    if retval is None:
        @wraps(func)
        def __sympifyit_wrapper(a, b):
            return func(a, sympify(b, strict=True))

    else:
        @wraps(func)
        def __sympifyit_wrapper(a, b):
            try:
                return func(a, sympify(b, strict=True))
            except SympifyError:
                return retval

    return __sympifyit_wrapper

_re_direct = re.compile(r"__([^r]\w+)__")
_re_reverse = re.compile(r"__r(\w+)__")
_op_special = {'__rshift__': (operator.rshift, False),
        '__rrshift__': (operator.rshift, True)}


def _get_operator(name):
    """Find the operator corresponding to a magic method name.

    Returns (operator, reversed) where reversed is a flag specifying whether
    the input is a reverse method, like '__radd__'.
    Returns (None, False) if `name` is not a magic method name.
    """
    try:
        return _op_special[name]
    except KeyError:
        pass
    op_name = None
    m = _re_direct.match(name)
    if m:
        reverse = False
        op_name = m.groups()[0]
    else:
        m = _re_reverse.match(name)
        if m:
            reverse = True
            op_name = m.groups()[0]
    if op_name:
        try:
            return getattr(operator, op_name), reverse
        except AttributeError:
            pass
    return None, False


def sympify_other(func):
    from sympy.core.basic import Basic
    op, reverse = _get_operator(func.__name__)
    if op is None:
        op = func
    if not reverse:
        @wraps(func)
        def decorated(self, other):
            if not isinstance(other, Basic):
                try:
                    return op(self, sympify(other, strict=True))
                except SympifyError:
                    return NotImplemented
            return func(self, other)
    else:
        @wraps(func)
        def decorated(self, other):
            if not isinstance(other, Basic):
                try:
                    return op(sympify(other, strict=True), self)
                except SympifyError:
                    return NotImplemented
            return func(self, other)
    return decorated
