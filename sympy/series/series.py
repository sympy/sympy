from sympy.core.sympify import sympify
from sympy.polys.ltaylor import taylor

def series(expr, x=None, x0=0, n=6, dir="+",pol_pars=[]):
    """Series expansion of expr around point `x = x0`.

    See the doctring of taylor and
    the doctring of Expr.series() for complete details of this wrapper.
    """
    return taylor(expr,x,x0,n,dir,pol_pars=pol_pars)
