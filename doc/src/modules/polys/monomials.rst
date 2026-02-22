================
Monomials
================

This section documents utilities for working with monomials in SymPy.

itermonomials
=============

.. autofunction:: sympy.polys.monomials.itermonomials

Examples
--------

Generate all monomials in variables ``x`` and ``y`` of total degree up to 2:

    >>> from sympy import symbols
    >>> from sympy.polys.monomials import itermonomials
    >>> x, y = symbols('x y')
    >>> list(itermonomials([x, y], 2))
    [1, x, y, x**2, x*y, y**2]

Generate monomials with degree between 1 and 3:

    >>> list(itermonomials([x, y], 3, 1))
    [x, y, x**2, x*y, y**2, x**3, x**2*y, x*y**2, y**3]
