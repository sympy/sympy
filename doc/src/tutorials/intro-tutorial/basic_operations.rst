.. _tutorial-basic:

==================
 Basic Operations
==================

Here we discuss some of the most basic operations needed for expression
manipulation in SymPy.  Some more advanced operations will be discussed later
in the :ref:`advanced expression manipulation <tutorial-manipulation>` section.

    >>> from sympy import *
    >>> x, y, z = symbols("x y z")

Substitution
============

One of the most common things you might want to do with a mathematical
expression is substitution.  Substitution replaces all instances of something
in an expression with something else.  It is done using the ``subs`` method.
For example

    >>> expr = cos(x) + 1
    >>> expr.subs(x, y)
    cos(y) + 1

Substitution is usually done for one of two reasons:

1. Evaluating an expression at a point. For example, if our expression is
   ``cos(x) + 1`` and we want to evaluate it at the point ``x = 0``, so that
   we get ``cos(0) + 1``, which is 2.

   >>> expr.subs(x, 0)
   2

2. Replacing a subexpression with another subexpression.  There are two
   reasons we might want to do this.  The first is if we are trying to build
   an expression that has some symmetry, such as `x^{x^{x^x}}`.  To build
   this, we might start with ``x**y``, and replace ``y`` with ``x**y``.  We
   would then get ``x**(x**y)``.  If we replaced ``y`` in this new expression
   with ``x**x``, we would get ``x**(x**(x**x))``, the desired expression.

   >>> expr = x**y
   >>> expr
   x**y
   >>> expr = expr.subs(y, x**y)
   >>> expr
   x**(x**y)
   >>> expr = expr.subs(y, x**x)
   >>> expr
   x**(x**(x**x))

   The second is if we want to perform a very controlled simplification, or
   perhaps a simplification that SymPy is otherwise unable to do.  For
   example, say we have `\sin(2x) + \cos(2x)`, and we want to replace
   `\sin(2x)` with `2\sin(x)\cos(x)`.  As we will learn later, the function
   ``expand_trig`` does this.  However, this function will also expand
   `\cos(2x)`, which we may not want.  While there are ways to perform such
   precise simplification, and we will learn some of them in the
   :ref:`advanced expression manipulation <tutorial-manipulation>` section, an
   easy way is to just replace `\sin(2x)` with `2\sin(x)\cos(x)`.

   >>> expr = sin(2*x) + cos(2*x)
   >>> expand_trig(expr)
   2*sin(x)*cos(x) + 2*cos(x)**2 - 1
   >>> expr.subs(sin(2*x), 2*sin(x)*cos(x))
   2*sin(x)*cos(x) + cos(2*x)

There are two important things to note about ``subs``.  First, it returns a
new expression.  SymPy objects are immutable.  That means that ``subs`` does
not modify it in-place.  For example

   >>> expr = cos(x)
   >>> expr.subs(x, 0)
   1
   >>> expr
   cos(x)
   >>> x
   x

.. sidebar:: Quick Tip

   SymPy expressions are immutable.  No function will change them in-place.

Here, we see that performing ``expr.subs(x, 0)`` leaves ``expr`` unchanged.
In fact, since SymPy expressions are immutable, no function will change them
in-place.  All functions will return new expressions.

To perform multiple substitutions at once, pass a list of ``(old, new)`` pairs
to ``subs``.

    >>> expr = x**3 + 4*x*y - z
    >>> expr.subs([(x, 2), (y, 4), (z, 0)])
    40

It is often useful to combine this with a list comprehension to do a large set
of similar replacements all at once.  For example, say we had `x^4 - 4x^3 + 4x^2 -
2x + 3` and we wanted to replace all instances of `x` that have an even power
with `y`, to get `y^4 - 4x^3 + 4y^2 - 2x + 3`.

    >>> expr = x**4 - 4*x**3 + 4*x**2 - 2*x + 3
    >>> replacements = [(x**i, y**i) for i in range(5) if i % 2 == 0]
    >>> expr.subs(replacements)
    -4*x**3 - 2*x + y**4 + 4*y**2 + 3

Converting Strings to SymPy Expressions
=======================================

The ``sympify`` function (that's ``sympify``, not to be confused with
``simplify``) can be used to convert strings into SymPy expressions.

For example

    >>> str_expr = "x**2 + 3*x - 1/2"
    >>> expr = sympify(str_expr)
    >>> expr
    x**2 + 3*x - 1/2
    >>> expr.subs(x, 2)
    19/2

.. warning:: ``sympify`` uses ``eval``.  Don't use it on unsanitized input.

``evalf``
=========

To evaluate a numerical expression into a floating point number, use
``evalf``.

    >>> expr = sqrt(8)
    >>> expr.evalf()
    2.82842712474619

SymPy can evaluate floating point expressions to arbitrary precision.  By
default, 15 digits of precision are used, but you can pass any number as the
argument to ``evalf``.  Let's compute the first 100 digits of `\pi`.

    >>> pi.evalf(100)
    3.141592653589793238462643383279502884197169399375105820974944592307816406286208998628034825342117068

To numerically evaluate an expression with a Symbol at a point, we might use
``subs`` followed by ``evalf``, but it is more efficient and numerically
stable to pass the substitution to ``evalf`` using the ``subs`` flag, which
takes a dictionary of ``Symbol: point`` pairs.

    >>> expr = cos(2*x)
    >>> expr.evalf(subs={x: 2.4})
    0.0874989834394464

Sometimes there are roundoff errors smaller than the desired precision that
remain after an expression is evaluated. Such numbers can be removed at the
user's discretion by setting the ``chop`` flag to True.

    >>> one = cos(1)**2 + sin(1)**2
    >>> (one - 1).evalf()
    -0.e-124
    >>> (one - 1).evalf(chop=True)
    0

``lambdify``
============

``subs`` and ``evalf`` are good if you want to do simple evaluation, but if
you intend to evaluate an expression at many points, there are more efficient
ways.  For example, if you wanted to evaluate an expression at a thousand
points, using SymPy would be far slower than it needs to be, especially if you
only care about machine precision.  Instead, you should use libraries like
`NumPy <http://www.numpy.org/>`_ and `SciPy <http://www.scipy.org/>`_.

The easiest way to convert a SymPy expression to an expression that can be
numerically evaluated is to use the ``lambdify`` function.  ``lambdify`` acts
like a ``lambda`` function, except it converts the SymPy names to the names of
the given numerical library, usually NumPy.  For example

    >>> import numpy # doctest:+SKIP
    >>> a = numpy.arange(10) # doctest:+SKIP
    >>> expr = sin(x)
    >>> f = lambdify(x, expr, "numpy") # doctest:+SKIP
    >>> f(a) # doctest:+SKIP
    [ 0.          0.84147098  0.90929743  0.14112001 -0.7568025  -0.95892427
     -0.2794155   0.6569866   0.98935825  0.41211849]

.. warning:: ``lambdify`` uses ``eval``.  Don't use it on unsanitized input.

You can use other libraries than NumPy. For example, to use the standard
library math module, use ``"math"``.

    >>> f = lambdify(x, expr, "math")
    >>> f(0.1)
    0.0998334166468

To use lambdify with numerical libraries that it does not know about, pass a
dictionary of ``sympy_name:numerical_function`` pairs.  For example

    >>> def mysin(x):
    ...     """
    ...     My sine. Note that this is only accurate for small x.
    ...     """
    ...     return x
    >>> f = lambdify(x, expr, {"sin":mysin})
    >>> f(0.1)
    0.1

.. TODO: Write an advanced numerics section

