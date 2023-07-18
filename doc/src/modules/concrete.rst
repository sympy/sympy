========
Concrete
========

Hypergeometric terms
--------------------

Hypergeometric terms play a prominent role in the theory of recurrences and
summation. A sequence `a(n)` is a hypergeometric term provided that the
consecutive term ratio `a(n + 1) / a(n)` is a rational function in `n`. Such
sequences are also called just "hypergeometric." Hypergeometric terms include
polynomials, rational functions, factorials, exponential functions, and
binomial coefficients, along with products and quotients of any of these.

The ``is_hypergeometric`` method checks if an expression is hypergeometric.
Here are some examples:

    >>> from sympy import *
    >>> n, k = symbols('n,k')
    >>> (n**2 + 1).is_hypergeometric(n)
    True
    >>> factorial(n).is_hypergeometric(n)
    True
    >>> binomial(n, k).is_hypergeometric(n)
    True
    >>> rf(n, k).is_hypergeometric(n)
    True
    >>> ff(n, k).is_hypergeometric(n)
    True
    >>> gamma(n).is_hypergeometric(n)
    True
    >>> (2**n).is_hypergeometric(n)
    True
    >>> binomial(n, k).is_hypergeometric(k)
    True
    >>> binomial(n, k).is_hypergeometric(n)
    True
    >>> rf(n, k).is_hypergeometric(k)
    True
    >>> ff(n, k).is_hypergeometric(k)
    True
    >>> (2**k).is_hypergeometric(k)
    True

Compositions of hypergeometric terms are usually not hypergeometric.

    >>> factorial(n**2).is_hypergeometric(n)
    False
    >>> (2**(n**3 + 1)).is_hypergeometric(n)
    False

The function ``hypersimp()`` will try to compute the simplified term ratio of a
hypergeometric term:

    >>> hypersimp(factorial(2*n), n)
    2*(n + 1)*(2*n + 1)
    >>> hypersimp(factorial(n**2), n)

Sums involving hypergeometric terms are extremely common in many areas of
discrete mathematics, and there are powerful algorithms that help evaluate
them. The ones SymPy implements are documented below.

Concrete Class Reference
------------------------
.. autoclass:: sympy.concrete.summations.Sum
   :members:

.. autoclass:: sympy.concrete.products.Product
   :members:

.. autoclass:: sympy.concrete.expr_with_intlimits.ExprWithIntLimits
   :members:

Concrete Functions Reference
----------------------------

.. autofunction:: sympy.concrete.summations.summation

.. autofunction:: sympy.concrete.products.product

.. autofunction:: sympy.concrete.gosper.gosper_normal

.. autofunction:: sympy.concrete.gosper.gosper_term

.. autofunction:: sympy.concrete.gosper.gosper_sum
