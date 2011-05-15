from sympy.core.decorators import wraps


def threaded_factory(func, use_add):
    """A factory for ``threaded`` decorators. """
    from sympy.core import sympify, Add
    from sympy.matrices import Matrix

    @wraps(func)
    def threaded_func(expr, *args, **kwargs):
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

    return threaded_func

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
    return threaded_factory(func, True)

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
    return threaded_factory(func, False)
