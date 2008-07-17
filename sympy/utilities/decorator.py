
from sympy.core.add import Add
from sympy.core.sympify import sympify
from sympy.core.relational import Relational

def threaded(**flags):
    """Call a function on all elements of composite objects.

       This decorator is intended to make it uniformly possible  to apply
       functions to all elements of composite or iterable objects or just
       to expressions. Currently supported, so called, threadable classes
       are Add, Matrix, Relational and all other which implement __iter__
       attribute.  If the input expression is not of this type,  then the
       function will be applied to this object, as to a singleton.

       Functions which will use this decorator must have the following
       signature:

          @threaded()
          def function(expr, *args, **kwargs):

       where 'expr' is obligatory threadable object (listed above). Other
       arguments, if available, are transferred without any change. As an
       example see functions in sympy.simplify module.

       By default threading is done on elements of Add instance. To avoid
       this behaviour set 'use_add' flag  with False in keyword arguments
       (see integrate() for details), e.g:

          @threaded(use_add=False)
          def function(expr, *args, **kwargs):

    """
    from sympy.matrices import Matrix

    use_add = flags.get('use_add', True)

    def threaded_proxy(func):
        def threaded_decorator(expr, *args, **kwargs):
            if isinstance(expr, Matrix):
                return expr.applyfunc(lambda f: func(f, *args, **kwargs))
            elif hasattr(expr, '__iter__'):
                return expr.__class__([ func(f, *args, **kwargs) for f in expr ])
            else:
                expr = sympify(expr)

                if isinstance(expr, Relational):
                    lhs = func(expr.lhs, *args, **kwargs)
                    rhs = func(expr.rhs, *args, **kwargs)

                    return Relational(lhs, rhs, expr.rel_op)
                elif expr.is_Add and use_add:
                    return Add(*[ func(f, *args, **kwargs) for f in expr.args ])
                else:
                    return func(expr, *args, **kwargs)

        threaded_decorator.__doc__  = func.__doc__
        threaded_decorator.__name__ = func.__name__

        return threaded_decorator

    return threaded_proxy
