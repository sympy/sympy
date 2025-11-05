.. _tips:


Tips
====

This page provides practical tips to help you use SymPy more effectively.
These tips cover common use cases, performance optimizations, and shortcuts
that can make your work with SymPy easier and more efficient.

.. contents:: Contents
   :local:

Quick Start Tips
================

Use isympy for Interactive Sessions
-----------------------------------

The :command:`isympy` command launches an enhanced IPython session with
SymPy already imported and common symbols predefined:

.. code-block:: bash

    $ isympy

This automatically runs:

    >>> from sympy import *
    >>> x, y, z, t = symbols('x y z t')
    >>> k, m, n = symbols('k m n', integer=True)
    >>> f, g, h = symbols('f g h', cls=Function)

This saves you from having to type these imports and definitions every time.

Use Tab Completion
------------------

In IPython or Jupyter notebooks, use the Tab key to:

- See available methods on an object: ``expr.<TAB>``
- Complete function names: ``simp<TAB>`` shows ``simplify``, ``simpify``, etc.
- See function signatures: ``expand(<TAB>``

Use ? and ?? for Help
---------------------

In IPython or Jupyter, you can get quick help on any function:

    >>> expand?        # Shows docstring
    >>> expand??       # Shows source code

Alternatively, use Python's built-in help:

    >>> help(expand)

Working with Symbols
====================

Define Multiple Symbols at Once
-------------------------------

Use the :func:`~.symbols` function to define multiple symbols efficiently:

    >>> a, b, c = symbols('a b c')
    >>> x1, x2, x3, x4, x5 = symbols('x1:6')  # Creates x1 through x5

Import Common Symbols from sympy.abc
------------------------------------

For quick scripts, import symbols from :mod:`sympy.abc`:

    >>> from sympy.abc import x, y, z, alpha, beta, gamma

This is especially useful for Greek letters and common variable names.

Add Assumptions When Defining Symbols
-------------------------------------

Adding assumptions helps SymPy simplify expressions more effectively:

    >>> x = symbols('x', positive=True)
    >>> sqrt(x**2)
    x
    >>> n = symbols('n', integer=True)
    >>> cos(n*pi)
    (-1)**n

Common assumptions include: ``real``, ``positive``, ``negative``,
``integer``, ``rational``, ``nonzero``, ``nonnegative``, ``finite``.

Working with Numbers
====================

Use Rational() for Exact Fractions
----------------------------------

Avoid using Python floats when you want exact arithmetic:

**Don't:**

    >>> x + 1/3
    x + 0.333333333333333

**Do:**

    >>> x + Rational(1, 3)
    x + 1/3

Or use the ``S`` shorthand:

    >>> x + S(1)/3
    x + 1/3

Use S() to Create SymPy Objects
-------------------------------

The ``S`` function converts Python objects to SymPy objects:

    >>> S(1)/2        # Creates Rational(1, 2)
    1/2
    >>> S('1/2')      # Parses string to Rational
    1/2
    >>> S(0.5)        # Converts float to exact representation
    0.5

Use SymPy Constants
-------------------

Always use SymPy's mathematical constants, not Python's:

**Don't:**

    >>> import math
    >>> sin(math.pi)
    1.22464679914735e-16

**Do:**

    >>> sin(pi)
    0

Other useful constants: :obj:`~.E` (Euler's number), :obj:`~.I`
(imaginary unit), :obj:`~.oo` (infinity), :obj:`~.zoo` (complex infinity).

Simplification and Manipulation
===============================

Choose the Right Simplification Function
----------------------------------------

Don't always reach for :func:`~.simplify`. Use targeted functions:

- :func:`~.expand` - Expand products and powers
- :func:`~.factor` - Factor polynomials
- :func:`~.cancel` - Cancel common factors in rational expressions
- :func:`~.collect` - Collect terms with same power
- :func:`~.trigsimp` - Simplify trigonometric expressions
- :func:`~.powsimp` - Simplify powers
- :func:`~.radsimp` - Rationalize denominators
- :func:`~.simplify` - General purpose (tries everything)

Use expand() Strategically
--------------------------

The :func:`~.expand` function has many options for selective expansion:

    >>> expand((x + 1)*(x + 2), mul=True)
    x**2 + 3*x + 2
    >>> expand(sin(x + y), trig=True)
    sin(x)*cos(y) + sin(y)*cos(x)
    >>> expand((x*y)**n, power_base=True)
    x**n*y**n

Use .rewrite() to Change Function Forms
---------------------------------------

The :meth:`~sympy.core.basic.Basic.rewrite` method converts between equivalent forms:

    >>> tan(x).rewrite(sin)
    2*sin(x)**2/sin(2*x)
    >>> factorial(n).rewrite(gamma)
    gamma(n + 1)
    >>> exp(I*x).rewrite(cos)
    cos(x) + I*sin(x)

Use together() for Combining Fractions
--------------------------------------

The :func:`~.together` function combines multiple fractions:

    >>> together(1/x + 1/y)
    (x + y)/(x*y)
    >>> together(1/(x + 1) + 1/(x + 2))
    (2*x + 3)/((x + 1)*(x + 2))

Substitution and Evaluation
===========================

Use .subs() for Substitution
----------------------------

Substitute values or expressions with :meth:`~sympy.core.basic.Basic.subs`:

    >>> expr = x**2 + 2*x + 1
    >>> expr.subs(x, 2)
    9
    >>> expr.subs(x, y + 1)
    (y + 1)**2 + 2*(y + 1) + 1

You can substitute multiple values at once:

    >>> expr = x + y + z
    >>> expr.subs({x: 1, y: 2})
    z + 3
    >>> expr.subs([(x, 1), (y, 2)])  # Order matters for sequential subs
    z + 3

Use .evalf() for Numerical Evaluation
-------------------------------------

Convert symbolic expressions to numerical values:

    >>> pi.evalf()
    3.14159265358979
    >>> pi.evalf(50)  # 50 digits of precision
    3.1415926535897932384626433832795028841971693993751

You can evaluate with substitution:

    >>> expr = sin(x)
    >>> expr.evalf(subs={x: 1})
    0.841470984807897

Use N() as Shorthand for evalf()
--------------------------------

:func:`~.N` is equivalent to :meth:`~sympy.core.evalf.EvalfMixin.evalf`:

    >>> N(pi)
    3.14159265358979
    >>> N(pi, 30)
    3.14159265358979323846264338328

Solving Equations
=================

Use solve() for Equations
-------------------------

The :func:`~.solve` function solves equations symbolically:

    >>> solve(x**2 - 4, x)
    [-2, 2]
    >>> solve(Eq(x + 2, 5), x)
    [3]

For systems of equations, pass multiple equations:

    >>> solve([x + y - 2, x - y], [x, y])
    {x: 1, y: 1}

Use solveset() for Better Control
---------------------------------

:func:`~.solveset` is more powerful and returns a :class:`~.FiniteSet`:

    >>> solveset(x**2 - 4, x)
    {-2, 2}
    >>> solveset(sin(x), x, domain=S.Reals)
    {2*n*pi | n in Integers} ∪ {2*n*pi + pi | n in Integers}

Check Your Domain
-----------------

Solutions depend on the domain. Specify when needed:

    >>> solveset(exp(x) - 1, x)
    {2*n*I*pi | n in Integers}
    >>> solveset(exp(x) - 1, x, domain=S.Reals)
    {0}

Calculus
========

Use diff() for Derivatives
--------------------------

Compute derivatives with :func:`~.diff`:

    >>> diff(x**3, x)
    3*x**2
    >>> diff(x**3, x, 2)  # Second derivative
    6*x
    >>> diff(x**2*y**3, x, y)  # Mixed partial derivative
    6*x*y**2

Use integrate() for Integration
-------------------------------

Compute integrals with :func:`~.integrate`:

    >>> integrate(x**2, x)
    x**3/3
    >>> integrate(x**2, (x, 0, 1))  # Definite integral
    1/3
    >>> integrate(exp(-x**2), (x, -oo, oo))  # Infinite bounds
    sqrt(pi)

Use limit() for Limits
----------------------

Compute limits with :func:`~.limit`:

    >>> limit(sin(x)/x, x, 0)
    1
    >>> limit(1/x, x, oo)
    0
    >>> limit(1/x, x, 0, '+')  # One-sided limit from right
    oo

Use series() for Taylor Series
------------------------------

Expand expressions as series with :func:`~.series`:

    >>> series(exp(x), x, 0, 6)
    1 + x + x**2/2 + x**3/6 + x**4/24 + x**5/120 + O(x**6)
    >>> series(1/(1 - x), x, 0, 5)
    1 + x + x**2 + x**3 + x**4 + O(x**5)

Remove O() Terms with .removeO()
--------------------------------

Remove order terms from series:

    >>> s = series(sin(x), x, 0, 6)
    >>> s
    x - x**3/6 + x**5/120 + O(x**6)
    >>> s.removeO()
    x**5/120 - x**3/6 + x

Matrices
========

Create Matrices Easily
----------------------

Use :class:`~.Matrix` for symbolic matrices:

    >>> M = Matrix([[1, 2], [3, 4]])
    >>> M
    Matrix([
    [1, 2],
    [3, 4]])

Access elements and properties:

    >>> M[0, 1]  # Element at row 0, column 1
    2
    >>> M.det()  # Determinant
    -2
    >>> M.inv()  # Inverse
    Matrix([
    [-2, 1],
    [3/2, -1/2]])

Use Eye for Identity Matrices
-----------------------------

Create identity matrices with :func:`~.eye`:

    >>> eye(3)
    Matrix([
    [1, 0, 0],
    [0, 1, 0],
    [0, 0, 1]])

Use zeros() and ones()
----------------------

Create matrices filled with zeros or ones:

    >>> zeros(2, 3)
    Matrix([
    [0, 0, 0],
    [0, 0, 0]])
    >>> ones(2, 2)
    Matrix([
    [1, 1],
    [1, 1]])

Matrix Operations
-----------------

Matrices support standard operations:

    >>> M = Matrix([[1, 2], [3, 4]])
    >>> M.T  # Transpose
    Matrix([
    [1, 3],
    [2, 4]])
    >>> M.eigenvals()  # Eigenvalues
    {5/2 - sqrt(33)/2: 1, sqrt(33)/2 + 5/2: 1}
    >>> M.eigenvects()  # Eigenvectors
    [(...), (...)]

Performance Tips
================

Use evaluate=False to Skip Simplification
-----------------------------------------

When building large expressions, disable automatic evaluation:

    >>> expr = Add(x, y, evaluate=False)
    >>> big_sum = Add(*range(1000), evaluate=False)

Then simplify once at the end:

    >>> expr.doit()

Use lambdify() for Fast Numerical Evaluation
--------------------------------------------

Convert SymPy expressions to fast numerical functions:

    >>> from sympy import lambdify
    >>> import numpy as np
    >>> f = lambdify(x, x**2 + 2*x + 1)
    >>> f(5)
    36
    >>> f(np.array([1, 2, 3]))  # Works with NumPy arrays
    array([ 4,  9, 16])

Use cse() for Common Subexpression Elimination
----------------------------------------------

Optimize expressions by eliminating repeated subexpressions:

    >>> from sympy import cse
    >>> expr1 = sin(x) + cos(x)
    >>> expr2 = sin(x) - cos(x)
    >>> replacements, reduced = cse([expr1, expr2])
    >>> replacements
    [(x0, sin(x)), (x1, cos(x))]
    >>> reduced
    [x0 + x1, x0 - x1]

Display and Printing
====================

Use pprint() for Pretty Printing
--------------------------------

The :func:`~.pprint` function displays expressions nicely:

    >>> pprint(Integral(exp(-x**2), (x, -oo, oo)))
      ∞
      ⌠
      ⎮   -x²
      ⎮  ℯ    dx
      ⌡
     -∞

Use init_printing() for Automatic Pretty Printing
-------------------------------------------------

Enable automatic pretty printing in IPython/Jupyter:

    >>> init_printing()
    >>> Integral(exp(-x**2), (x, -oo, oo))
    # Displays in nice mathematical notation

Use latex() to Generate LaTeX Code
----------------------------------

Convert expressions to LaTeX for papers and presentations:

    >>> from sympy import latex
    >>> latex(Integral(exp(-x**2), (x, -oo, oo)))
    '\\int\\limits_{-\\infty}^{\\infty} e^{- x^{2}}\\, dx'

Use mathematica_code() and other code generation
------------------------------------------------

Convert to other systems:

    >>> from sympy import mathematica_code, octave_code
    >>> mathematica_code(x**2 + y**2)
    'x^2 + y^2'

Testing and Debugging
=====================

Use .equals() to Test Mathematical Equality
-------------------------------------------

The ``==`` operator tests structural equality, not mathematical equality.
Use :meth:`~sympy.core.expr.Expr.equals` for mathematical equality:

    >>> (x + 1)**2 == x**2 + 2*x + 1
    False
    >>> ((x + 1)**2).equals(x**2 + 2*x + 1)
    True

Check if Expression is Zero
---------------------------

To check if an expression is mathematically zero:

    >>> expr = (x + 1)**2 - (x**2 + 2*x + 1)
    >>> simplify(expr)
    0
    >>> simplify(expr) == 0
    True

Use .is_* Properties
--------------------

Check properties of expressions:

    >>> x = symbols('x', real=True)
    >>> (x**2).is_positive  # Returns True, False, or None (unknown)
    >>> (x**2).is_real
    True
    >>> I.is_imaginary
    True

Common Pitfalls and Solutions
=============================

Remember: ** not ^ for Exponentiation
-------------------------------------

Python uses ``**`` for powers, not ``^``:

    >>> x**2  # Correct
    x**2
    >>> # x^2 gives x XOR 2, not x squared

Always Use * for Multiplication
-------------------------------

Implicit multiplication doesn't work:

    >>> 2*x  # Correct
    2*x
    >>> # 2x gives SyntaxError

Variables Must Be Defined Before Use
------------------------------------

Unlike Mathematica or Maple, symbols must be created first:

    >>> x = symbols('x')  # Must define first
    >>> y**2  # Error if y not defined

Changing a Variable Doesn't Change Expressions
----------------------------------------------

Once an expression is created, changing variables doesn't affect it:

    >>> x = symbols('x')
    >>> expr = x + 1
    >>> x = 2  # This doesn't change expr
    >>> expr
    x + 1  # Still symbolic

Use ``.subs()`` to change values in expressions.

Further Resources
=================

- :ref:`Introductory Tutorial <intro-tutorial>` - Start here if you're new
- :ref:`Gotchas and Pitfalls <gotchas>` - Common mistakes to avoid
- :ref:`Best Practices <best-practices>` - Guidelines for writing good SymPy code
- :ref:`API Reference <reference>` - Complete function documentation
- `SymPy Tutorial <https://docs.sympy.org/latest/tutorial/index.html>`_ - Official tutorial
- `SymPy Examples <https://github.com/sympy/sympy/wiki/Quick-examples>`_ - Code examples on the wiki
