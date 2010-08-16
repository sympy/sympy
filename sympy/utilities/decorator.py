import warnings

def threaded_factory(func, use_add):
    """A factory for ``threaded`` decorators. """
    from sympy.core import sympify, Add
    from sympy.matrices import Matrix

    def threaded_decorator(expr, *args, **kwargs):
        if isinstance(expr, Matrix):
            return expr.applyfunc(lambda f: func(f, *args, **kwargs))
        elif hasattr(expr, '__iter__'):
            return expr.__class__([ func(f, *args, **kwargs) for f in expr ])
        else:
            expr = sympify(expr)

            if use_add and expr.is_Add:
                return expr.__class__(*[ func(f, *args, **kwargs) for f in expr.args ])
            elif expr.is_Relational:
                return expr.__class__(func(expr.lhs, *args, **kwargs),
                                      func(expr.rhs, *args, **kwargs))
            else:
                return func(expr, *args, **kwargs)

    return threaded_decorator

def threaded(func):
    """Apply ``func`` to sub--elements of an object, including :class:`Add`.

       This decorator is intended to make it uniformly possible to apply a
       function to all elements of composite objects, e.g. matrices, lists,
       tuples and other iterable containers, or just expressions.

       This version of :func:`threaded` decorator allows threading over
       elements of :class:`Add` class. If this behavior is not desirable
       use :func:`xthreaded` decorator.

       Functions using this decorator must have the following signature::

          @threaded
          def function(expr, *args, **kwargs):

    """
    return wraps(func, threaded_factory(func, True))

def xthreaded(func):
    """Apply ``func`` to sub--elements of an object, excluding :class:`Add`.

       This decorator is intended to make it uniformly possible to apply a
       function to all elements of composite objects, e.g. matrices, lists,
       tuples and other iterable containers, or just expressions.

       This version of :func:`threaded` decorator disallows threading over
       elements of :class:`Add` class. If this behavior is not desirable
       use :func:`threaded` decorator.

       Functions using this decorator must have the following signature::

          @xthreaded
          def function(expr, *args, **kwargs):

    """
    return wraps(func, threaded_factory(func, False))

def deprecated(func):
    """This is a decorator which can be used to mark functions
    as deprecated. It will result in a warning being emitted
    when the function is used."""
    def new_func(*args, **kwargs):
        warnings.warn("Call to deprecated function %s." % func.__name__,
                      category=DeprecationWarning)
        return func(*args, **kwargs)
    new_func.__name__ = func.__name__
    new_func.__doc__ = func.__doc__
    new_func.__dict__.update(func.__dict__)
    return new_func

def wraps(old_func, new_func):
    """Copy private data from ``old_func`` to ``new_func``. """
    new_func.__dict__.update(old_func.__dict__)

    new_func.__module__ = old_func.__module__
    new_func.__name__   = old_func.__name__
    new_func.__doc__    = old_func.__doc__

    return new_func

