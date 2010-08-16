"""Tools for applying functions to specified parts of expressions. """

from sympy.core import sympify

def use(expr, func, level=0, args=(), kwargs={}):
    """
    Use ``func`` to transform ``expr`` at the given level.

    Example
    =======

    >>> from sympy import use, expand
    >>> from sympy.abc import x, y

    >>> f = (x + y)**2*x + 1

    >>> use(f, expand, level=2)
    1 + x*(2*x*y + x**2 + y**2)

    >>> expand(f)
    1 + x*y**2 + 2*y*x**2 + x**3

    """
    def _use(expr, level):
        if not level:
            return func(expr, *args, **kwargs)
        else:
            if expr.is_Atom:
                return expr
            else:
                level -= 1
                _args = []

                for arg in expr.args:
                    _args.append(_use(arg, level))

                return expr.__class__(*_args)

    return _use(sympify(expr), level)
