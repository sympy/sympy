=========
 Solvers
=========

    >>> from sympy import *
    >>> x, y, z = symbols('x y z')
    >>> init_printing(use_unicode=True)

A Note about Equations
======================

Recall from the :ref:`gotchas section <tutorial_gotchas_equals>` of this
tutorial that symbolic equations in SymPy are not represented by ``=`` or
``==``, but by ``Eq``.


    >>> Eq(x, y)
    x = y


However, there is an even easier way.  In SymPy, any expression is not in an
``Eq`` is automatically assumed to equal 0 by the solving functions.  Since `a
= b` if and only if `a - b = 0`, this means that instead of using ``x == y``,
you can just use ``x - y``.  For example

    >>> solve(Eq(x**2, 1), x)
    [-1, 1]
    >>> solve(Eq(x**2 - 1, 0), x)
    [-1, 1]
    >>> solve(x**2 - 1, x)
    [-1, 1]

This is particularly useful if the equation you wish to solve is already equal
to 0.  Instead of typing ``solve(Eq(expr, 0), x)``, you can just use
``solve(expr, x)``.

Solving Algebraic Equations
===========================

The main function for solving algebraic equations, as we saw above, is
``solve``.  The syntax is ``solve(equations, variables)``, where, as we saw
above, ``equations`` may be in the form of ``Eq`` instances or expressions
that are assumed to be equal to zero.

When solving a single equation, the output of solve is a list of the
solutions.

    >>>

If no solutions are found, an empty list is returned, or
``NotImplementedError`` is raised.

    >>> solve


.. note::

   If solve returns ``[]`` or raises ``NotImplementedError``, it doesn't mean
   that the equation has no solutions.  It just means that it couldn't find
   any.  Often this means that the solutions cannot be represented
   symbolically

   >>>
