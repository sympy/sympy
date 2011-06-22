"""
This module implements the Residue function and related tools for working with
residues.
"""

from sympy import Wild, sympify, Integer, Add

def residue(expr, x, x0):
    """
    Finds the residue of ``expr`` at the point x=x0.

    The residue is defined as the coefficient of 1/(x-x0) in the power series
    expansion about x=x0.

    Examples:

    >>> from sympy import Symbol, residue, sin
    >>> x = Symbol("x")
    >>> residue(1/x, x, 0)
    1
    >>> residue(1/x**2, x, 0)
    0
    >>> residue(2/sin(x), x, 0)
    2

    This function is essential for the Residue Theorem [1].

    The current implementation uses series expansion to calculate it. A more
    general implementation is explained in the section 5.6 of the Bronstein's
    book [2]. For purely rational functions, the algorithm is much easier. See
    sections 2.4, 2.5, and 2.7 (this section actually gives an algorithm for
    computing any Laurent series coefficient for a rational function). The
    theory in section 2.4 will help to understand why the resultant works in
    the general algorithm. For the definition of a resultant, see section 1.4
    (and any previous sections for more review).

    [1] http://en.wikipedia.org/wiki/Residue_theorem
    [2] M. Bronstein: Symbolic Integration I, Springer Verlag (2005)

    """
    from sympy import collect
    expr = sympify(expr)
    if x0 != 0:
        expr = expr.subs(x, x+x0)
    s = expr.series(x, 0, 0).removeO()
    # TODO: this sometimes helps, but the series expansion should rather be
    # fixed, see #1627:
    if s == 0:
        s = expr.series(x, 0, 6).removeO()
    s = collect(s, x)
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
