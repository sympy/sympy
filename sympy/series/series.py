def series(expr, x, x0=0, n=6, dir="+"):
    """Series expansion of expr around point `x = x0`.

    See the doctring of Expr.series() for complete details of this wrapper.
    """
    return expr.series(x, x0, n, dir)
