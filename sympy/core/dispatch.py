from __future__ import print_function, division

__all__ = ["dispatch"]


class DispatchedFunction(object):
    """
    Class to handle dispatched functions and methods.
    """
    _fnames_ = dict()
    __slots__ = ("_ftypes_", "_ftypes_with_varargs_", "_fcache_", "_last_added_function_",)

    def __new__(cls, types, func, prefix="", varargs=None):
        if varargs is None:
            ftypes = "_ftypes_"
        else:
            ftypes = "_ftypes_with_varargs_"
            types = types + (varargs,)
        if prefix:
            func_name = prefix + "." + func.__name__
        else:
            func_name = func.__name__

        if func_name in cls._fnames_:
            obj = DispatchedFunction._fnames_[func_name]
            getattr(obj, ftypes)[types] = func
            # destroy any previous cache when add a new dispatching rule,
            # otherwise errors with subtypes may occur.
            obj._fcache_ = dict()
            obj._last_added_function_ = func
            return obj

        obj = super(DispatchedFunction, cls).__new__(cls)
        obj._ftypes_ = dict()
        obj._ftypes_with_varargs_ = dict()
        setattr(obj, ftypes, dict({types: func}))
        obj._fcache_ = dict()
        obj._last_added_function_ = func
        DispatchedFunction._fnames_[func_name] = obj
        return obj

    def __call__(self, *args, **kw_args):
        arg_types = tuple(type(arg) for arg in args)
        # check if such argument types are already stored in the cache,
        # in that case, retrieve the function call there, in order to avoid
        # to perform the whole search-algorithm:
        if arg_types in self._fcache_:
            return self._fcache_[arg_types](*args, **kw_args)

        def check_match(types, func):
            # this function first checks that types are indeed instances of the
            # arguments. After that it performs a for-loop over all
            # previously-found matches, and checks if it is possible to
            # determine which match is more specific (i.e. has lowest subclasses
            # inside its `__mro__`, method resolution order, Python's default
            # method searching order in classes). If a match is more specific
            # than another match, the less specific match is dropped.

            # Matches are stored as a static variable of this function,
            # namely, `check_match.mdata`. It is a python `list` instance.

            # In the end, there will be three cases. Either no matches have
            # been found, or just one, or multiple ones. No matches and more
            # than ones raise an error, while a single match passes execution
            # over to the matched function.

            if not all(map(lambda x: isinstance(*x), zip(args, types))):
                # incompatible types, return so that it will skip to the
                # next one.
                return

            # determine the method resolution order index of the current
            # types inside the `__mro__` field of the function arguments:
            mroord = []
            for arg_type, _type in zip(arg_types, types):
                if _type not in arg_type.__mro__:
                    _type = object
                mroord.append(arg_type.__mro__.index(_type))
            mroord = tuple(mroord)

            # Now compare this with all previous matches, more specific
            # matches are dropped:
            pop_list = []

            for indx, comp_data in enumerate(check_match.mdata):
                comp_types, comp_func, comp_mroord = comp_data
                if len(comp_types) != len(types):
                    1/0  # should never be here
                    continue
                less_eq_comp = [i <= j for i, j in zip(mroord, comp_mroord)]
                if all(less_eq_comp):
                    # the current match is more specific than the previous one,
                    # drop the previous one:
                    pop_list.append(indx)
                    continue
                if all([i >= j for i, j in zip(mroord, comp_mroord)]):
                    # a previous match is more specific than the current one,
                    # drop this match:
                    return

            for i in reversed(pop_list):
                check_match.mdata.pop(i)
            # append data about current match to a static variable:
            check_match.mdata.append((types, func, mroord))

        check_match.mdata = []

        for types, func in self._ftypes_.items():
            # for-loop to match types without varargs (variable arguments)
            if len(types) != len(args):
                continue
            check_match(types, func)

        for types, func in self._ftypes_with_varargs_.items():
            # for-loop to match types with varargs (variable arguments),
            # last type in `types` is repeated as many times as the number
            # of missing arguments.
            vararg_type = types[-1]
            fixed_types = types[:-1]

            number_of_varargs = len(args) - len(types[:-1])
            if number_of_varargs < 0:
                # too few arguments, impossible to match
                continue
            # extend the types to match by the vararg type, multiplied by
            # the number of missing arguments to match the function args:
            types_to_match = fixed_types + (vararg_type,)*number_of_varargs

            check_match(types_to_match, func)

        if not check_match.mdata:
            # no matches have been found, raise an error
            raise TypeError("dispatch: no match found for signature {0}".format(types))
        elif len(check_match.mdata) > 1:
            # ambiguous matches, multiple matches have been found, but it was
            # not possible to determine which one was the more specific, raise
            # an error:
            raise TypeError("dispatch: ambiguous type resolution for signature {0}".format(types))

        # store the result into _fcache_, so that next time the dispatched
        # function is called with the same argument types, it will be
        # faster to retrieve it (the current function, which has just
        # determined the match, will not be called anymore):
        matched_func = check_match.mdata[0][1]
        self._fcache_[arg_types] = matched_func

        return matched_func(*args, **kw_args)

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
