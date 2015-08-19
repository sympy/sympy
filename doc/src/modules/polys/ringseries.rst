.. _polys-ringseries:

=====================================
Series Manipulation using Polynomials
=====================================

Any finite Taylor series, for all practical purposes is, in fact a polynomial.
This module makes use of the efficient representation and operations of sparse
polynomials for very fast multivariate series manipulations. Typical speedups
compared to SymPy's ``series`` method are in the range 20-100, with the gap
widening as the series being handled gets larger. Though the definition of a
polynomial limits the use of Polynomial module to Taylor series, we extend it
to allow Laurent and even Puiseux series (with fractional exponents).

All the functions expand any given series on some ring specified by the user.
Thus, the coefficients of the calculated series depend on the ring being used.
For example::

    >>> R, x, y = ring('x, y', QQ)
    >>> rs_sin(x*y, x, 5)
    -1/6*x**3*y**3 + x*y

``QQ`` stands for the Rational domain. Here all coefficients are rationals.
Similarly, if a Real domain is used::

    >>> R, x, y = ring('x, y', RR)
    >>> rs_sin(x*y, x, 5)
    -0.166666666666667*x**3*y**3 + x*y

All series returned by the functions of this module are objects of the
``PolyElement`` class. To use them with other SymPy types, convert them  to ``Expr``::

    >>> a = symbols('a')
    >>> series = rs_exp(x, x, 5)
    >>> a + series.as_expr()
    a + x**4/24 + x**3/6 + x**2/2 + x + 1

Manipulation of power series
****************************************************************************
.. currentmodule:: sympy.polys.ring_series

Functions in this module carry the prefix ``rs_``, standing for "ring series".
They manipulate finite power series in the sparse representation provided
by ``polys.ring.ring``.

**Elementary functions**

.. autofunction:: rs_log
.. autofunction:: rs_LambertW
.. autofunction:: rs_exp
.. autofunction:: rs_atan
.. autofunction:: rs_asin
.. autofunction:: rs_tan
.. autofunction:: rs_cot
.. autofunction:: rs_sin
.. autofunction:: rs_cos
.. autofunction:: rs_cos_sin
.. autofunction:: rs_atanh
.. autofunction:: rs_sinh
.. autofunction:: rs_cosh
.. autofunction:: rs_tanh
.. autofunction:: rs_hadamard_exp

**Operations**

.. autofunction:: rs_mul
.. autofunction:: rs_square
.. autofunction:: rs_pow
.. autofunction:: rs_series_inversion
.. autofunction:: rs_series_reversion
.. autofunction:: rs_nth_root
.. autofunction:: rs_trunc
.. autofunction:: rs_subs
.. autofunction:: rs_diff
.. autofunction:: rs_integrate
.. autofunction:: rs_newton
.. autofunction:: rs_compose_add

**Utility functions**

.. autofunction:: rs_is_puiseux
.. autofunction:: rs_puiseux
.. autofunction:: rs_puiseux2
.. autofunction:: rs_series_from_list
.. autofunction:: rs_fun
.. autofunction:: mul_xin
.. autofunction:: pow_xin
