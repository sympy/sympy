"""
SymPy core decorators.

The purpose of this module is to expose decorators without any other
dependencies, so that they can be easily imported anywhere in sympy/core.
"""

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

    def __sympifyit_wrapper(a, b):
        # our task is to call `func` with `b` _sympified.
        #
        # if we _sympify from the beginning, we'll get unneccesary overhead,
        # because _sympify has far non-zero cost even for Basic instances.
        #
        # the idea is to first run `func` with `b` as is, catch any error, and
        # try to rerun with b _sympified.
        #
        # so for Basic instances we'll get almost no overhead, and for other
        # objects we agree to take additional overhead because `func` has to be
        # run first, and only when it raises we can proceed with the second
        # phase.
        #
        # however there is one important exception -- python ints.
        # ints are used heavily, e.g. in sum([x**i for i in range(n)]) and
        # other places, so it is important to sympify ints as fast as possible
        # too.


        # python ints are used frequently -- it is important to convert them as
        # fast as possible
        #
        # %timeit type(1) is int            ->  1.43 us
        # %timeit type('abc') is int        ->  1.48 us
        # %timeit isinstance(1, int)        ->  1.29 us
        # %timeit isinstance('abc', int)    ->  2.23 us
        # %timeit isinstance(x, int)        ->  4.28 us
        # z = S.Half
        # %timeit isinstance(z, int)        ->  5.2 us
        #
        # so we use:
        if type(b) is int:
            from numbers import Integer
            b = Integer(b)

        try:
            # fast-path: let's hope b is already SymPy object
            return func(a, b)

        except Exception, e:
            from sympify import SympifyError
            # we've got an exception.
            # maybe it's from nested __sympifyit? then we have to quit.
            if isinstance(e, SympifyError):
                #print 'double deep sympify'
                if retval is not None:
                    return retval
                else:
                    raise

            # slow-path: b seems to be not SymPy object -- let's _sympify it
            from sympify import _sympify
            try:
                b = _sympify(b)
                #print 'deep sympify'
            except SympifyError:
                # sympify failed, let's return requested value
                if retval is not None:
                    return retval
                else:
                    # or pass exception through
                    raise

            # b successfully _sympified, lets call func again.
            # if it raises here -- exception goes to caller
            return func(a, b)

    return __sympifyit_wrapper
