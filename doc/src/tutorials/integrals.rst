.. _integrate-partial-fractions:

=========================================================
Example: Integration of Rational Functions using SymPy
=========================================================

SymPy can automatically perform integration on rational functions using
partial fraction decomposition. This is a common operation when dealing
with symbolic calculus, and SymPy makes it very simple to execute.

This example demonstrates how to integrate a rational function and verify
the result.

.. code-block:: python

   from sympy import symbols, integrate

   # Define the variable
   x = symbols('x')

   # Define the rational function
   expr = (2*x + 3) / (x**2 + 3*x + 2)

   # Perform the integration
   result = integrate(expr, x)
   print(result)

The output will be:

.. code-block:: text

   log(x + 1) + log(x + 2)

SymPy automatically decomposes the rational expression into partial fractions
and integrates them.

You can verify that the result is correct by differentiating the integral:

.. code-block:: python

   from sympy import diff, simplify

   simplify(diff(result, x) - expr)

This returns ``0``, confirming that the integration result is correct.

---

**Explanation**

- ``integrate(expr, x)`` computes the indefinite integral of ``expr`` with respect to ``x``.
- For rational functions, SymPy internally uses *partial fraction decomposition* when appropriate.
- You can use ``diff(result, x)`` to differentiate the result and confirm it matches the original expression.

This simple example demonstrates how SymPy’s integration engine handles rational
functions efficiently while allowing users to symbolically verify results.
