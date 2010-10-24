"""
This module implements the Residue function and related tools for working with
residues.
"""

from sympy import Wild, sympify, Integer, Add

def residue(expr, x, x0):
    """
    Finds the residue of ``expr`` at the point x=x0.

    The residue is defined as the coefficient of 1/(x-x0) in the power series
    expansion about x=x0. The motivation behind this is that the limit(expr, x,
    x0) tells you the value of expr when x goes to x0. If the limit is
    infinite, it is often useful to know the residue.

    Examples:

    >>> from sympy import Symbol, residue, sin
    >>> x = Symbol("x")
    >>> residue(1/x, x, 0)
    1
    >>> residue(1/x**2, x, 0)
    0
    >>> residue(2/sin(x), x, 0)
    2

    """
    expr = sympify(expr)
    if x0 != 0:
        expr = expr.subs(x, x+x0)
    s = expr.series(x, 0, 0).removeO()
    # TODO: this sometimes helps, but the series expansion should rather be
    # fixed, see #1627:
    if s == 0:
        s = expr.series(x, 0, 6).removeO()
    if x0 != 0:
        s = s.subs(x, x-x0)
    a = Wild("r", exclude=[x])
    c = Wild("c", exclude=[x])
    r = s.match(a/(x-x0)+c)
    if r:
        return r[a]
    elif isinstance(s, Add):
        # TODO: this is to overcome a bug in match (#1626)
        for t in s.args:
            r = t.match(a/(x-x0)+c)
            if r:
                return r[a]
    return Integer(0)
