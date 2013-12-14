from __future__ import print_function, division

__all__ = ["dispatch"]


class DispatchedFunction(object):
    """
    Class to handle dispatched functions and methods.
    """
    _fnames_ = dict()
    __slots__ = ("_ftypes_", "_ftypes_with_varargs_", "_fcache_", "_last_added_function_",)

    def __new__(cls, types, func, prefix="", varargs=None):
        if prefix:
            func_name = prefix + "." + func.__name__
        else:
            func_name = func.__name__

        if func_name in cls._fnames_:
            obj = DispatchedFunction._fnames_[func_name]
            obj._ftypes_[types] = func
            # destroy any previous cache when add a new dispatching rule,
            # otherwise errors with subtypes may occur.
            obj._fcache_ = dict()
            obj._last_added_function_ = func
            return obj

        obj = super(DispatchedFunction, cls).__new__(cls)
        obj._ftypes_ = dict()
        obj._ftypes_[types] = func
        obj._fcache_ = dict()
        obj._last_added_function_ = func
        DispatchedFunction._fnames_[func_name] = obj
        return obj

    def __call__(self, *args, **kw_args):
        arg_types = tuple(type(_) for _ in args)
        # check if such argument types are already stored in the cache,
        # in that case, retrieve the function call there, in order to avoid
        # to perform the whole search-algorithm:
        if arg_types in self._fcache_:
            return self._fcache_[arg_types](*args, **kw_args)

        mfunc = None
        mtypes = None
        for types, func in self._ftypes_.items():
            if len(types) != len(args):
                continue
            if not all(map(lambda x: isinstance(*x), zip(args, types))):
                # incompatible types, skip to next.
                continue

            if mtypes is not None:
                # this means that a match has already been found,
                # select the one with the most specific subclasses.
                for arg_type, _type, m_type in zip(arg_types, types, mtypes):
                    if _type not in arg_type.__mro__:
                        _type = object
                    if m_type not in arg_type.__mro__:
                        m_type = object
                    ord1 = arg_type.__mro__.index(_type)
                    ord2 = arg_type.__mro__.index(m_type)
                    if ord1 == ord2:
                        continue
                    elif ord1 < ord2:
                        mtypes = types
                        mfunc = func
                        break
                    else:
                        break
            else:
                mtypes = types
                mfunc = func

        if mfunc is None:
            raise TypeError("no match found for signature {}".format(types))

        # store the result into _fcache_, so next time it will be
        # faster to retrieve it:
        self._fcache_[arg_types] = mfunc

        return mfunc(*args, **kw_args)

    @property
    def methods(self):
        return self._ftypes_.items()

    def __len__(self):
        return len(self._ftypes_.keys())


def dispatch(*types, **kwargs):
    """
    Create dispatched functions and methods.

    This decorator replaces the function
    with a ``DispatchedFunction`` object, which handles links according to
    the types of arguments.

    Example
    =======

    >>> from sympy import Add, Mul, symbols, Basic
    >>> from sympy.core.dispatch import dispatch
    >>> x, y = symbols('x, y')
    >>> @dispatch(Add)
    ... def f(a):
    ...     return "addition"
    >>> @dispatch(Mul)
    ... def f(a):
    ...     return "multiplication"
    >>> f(x+y)
    'addition'
    >>> f(x*y)
    'multiplication'

    The most specific subclass is selected, in this example the type of ``x``
    is an instance of ``Basic``, so it will match ``Basic`` rather than
    ``object``. The type of number ``1`` instead does not subclass ``Basic``,
    so it will match ``object`` (the most generic match):

    >>> @dispatch(Basic)
    ... def g(a):
    ...     return str(a) + " is an instance of Basic"
    >>> @dispatch(object)
    ... def g(a):
    ...     return str(a) + " is an instance of object"
    >>> g(x)
    'x is an instance of Basic'
    >>> g(1)
    '1 is an instance of object'
    """
    prefix = kwargs.get("prefix", "")
    varargs = kwargs.get("varargs", None)

    def decorator_f(f):
        if isinstance(f, DispatchedFunction):
            return DispatchedFunction(types, f._last_added_function_, prefix, varargs)
        else:
            return DispatchedFunction(types, f, prefix, varargs)
    return decorator_f
