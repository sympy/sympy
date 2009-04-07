from sympy.core import sympify, Basic

def series(expr, x, point=0, n=6, dir="+"):
    """
    Series expansion of expr(x) around "point".

    This is just a wrapper around Basic.series(), see its docstring of for a
    thourough docstring for this function. In isympy you can just type
    Basic.series? and enter.
    """
    return expr.series(x, point, n, dir)
