.. module:: sympy.holonomic.holonomic

Uses and Current limitations
============================

Integration
-----------

One can perform integrations using holonomic functions by following these steps:

1. Convert the integrand to a holonomic function.
2. Now integrate the holonomic representation of the function.
3. Convert the integral back to expressions.

Examples
^^^^^^^^

>>> from sympy.abc import x, a
>>> from sympy import sin
>>> from sympy.holonomic import expr_to_holonomic
>>> expr_to_holonomic(1/(x**2+a), x).integrate(x).to_expr()
atan(x/sqrt(a))/sqrt(a)
>>> expr_to_holonomic(sin(x)/x).integrate(x).to_expr()
Si(x)

As you can see in the first example we converted the function to holonomic, integrated the result
and then converted back to symbolic expression.

Limitations
-----------

1. Converting to expressions is not always possible. The holonomic function
should have a hypergeometric series at ``x0``.

2. Implementation of converting to holonomic sequence currently doesn't support
``Frobenius method`` when the solutions need to have `\log` terms. This happens
when at least one pair of the roots of the indicial equation differ by an integer and
frobenius method yields linearly dependent series solutions. Since we use this while converting
to expressions, sometimes :func:`~HolonomicFunction.to_expr()` fails.

3. There doesn't seem to be a way for computing indefinite integrals, so :func:`~HolonomicFunction.integrate()`
basically computes `\int_{x_0}^{x} f(x)dx` if no limits are given, where `x_0` is the point at
which initial conditions for the integrand are stored. Sometimes this gives an additional constant in the result.
For instance:

>>> expr_to_holonomic(sin(x)).integrate(x).to_expr()
-cos(x) + 1
>>> sin(x).integrate(x)
-cos(x)

The indefinite integral of `\sin(x)` is `-\cos(x)`. But the output is `-\cos(x) + 1`
which is `\int_{0}^{x} sin(x)dx`. Although both are considered correct but `-\cos(x)`
is simpler.