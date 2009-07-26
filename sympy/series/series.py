from sympy.core import sympify, Basic

def series(expr, x, point=0, n=6, dir="+"):
    """Series expansion of expr around point `x=point`.

    See the doctring of Basic.series() for complete details of this wrapper.
    """
    return expr.series(x, point, n, dir)
