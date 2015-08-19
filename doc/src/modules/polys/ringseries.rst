.. _polys-ringseries:

=====================================
Series Manipulation using Polynomials
=====================================

Any finite Taylor series, for all practical purposes is, in fact a polynomial.
This module makes use of the efficient representation and operations of sparse
polynomials for very fast multivariate series manipulations. Typical speedups
compared to SymPy's ``series`` method are in the range 20-100, with the gap
widening as the series being handled gets larger.

All the functions expand any given series on some ring specified by the user.
Thus, the coefficients of the calculated series depend on the ring being used.
For example::

    >>> R, x, y = ring('x, y', QQ)
    >>> rs_sin(x*y, x, 5)
    -1/6*x**3*y**3 + x*y

``QQ`` stands for the Rational domain. Here all coefficients are rationals. It
is recommended to use ``QQ`` with ring series as it automatically chooses the
fastest Rational type.

Similarly, if a Real domain is used::

    >>> R, x, y = ring('x, y', RR)
    >>> rs_sin(x*y, x, 5)
    -0.166666666666667*x**3*y**3 + x*y

Though the definition of a polynomial limits the use of Polynomial module to
Taylor series, we extend it to allow Laurent and even Puiseux series (with
fractional exponents)::

    >>> R, x, y = ring('x, y', QQ)
    >>> rs_cos(x**QQ(1, 3) + x*y, x, 2)
    -x**(4/3)*y + 1/24*x**(4/3) - 1/2*x**(2/3) + 1
    >>> rs_tan(x**QQ(2, 5)*y**QQ(1, 2), x, 2)
    1/3*x**(6/5)*y**(3/2) + x**(2/5)*y**(1/2)

All series returned by the functions of this module are instances of the
``PolyElement`` class. To use them with other SymPy types, convert them  to
``Expr``::

    >>> a = symbols('a')
    >>> series = rs_exp(x, x, 5)
    >>> a + series.as_expr()
    a + x**4/24 + x**3/6 + x**2/2 + x + 1

rs_series
=========

Direct use of elementary ring series functions does give more control, but is
limiting at the same time. Creating an appropriate ring for the desired series
expansion and knowing which ring series function to call, are things not
everyone might be familiar with.

`rs\_series` is a function that takes an arbitrary ``Expr`` and returns its
expansion by calling the appropriate ring series functions. The returned series
is over the simplest (almost) possible ring that does the job. It recursively
builds the ring as it parses the given expression, adding generators to the
ring when it needs them. Some examples::

    >>> from sympy.abc import a, b, c
    >>> rs_series(sin(a + b), a, 5)
    1/24*sin(b)*a**4 - 1/2*sin(b)*a**2 + sin(b) - 1/6*cos(b)*a**3 + cos(b)*a
    >>> rs_series(sin(exp(a*b) + cos(a + c)), a, 2)
    -sin(c)*cos(cos(c) + 1)*a + cos(cos(c) + 1)*a*b + sin(cos(c) + 1)
    >>> rs_series(sin(a + b)*cos(a + c)*tan(a**2 + b), a, 2)
    cos(b)*cos(c)*tan(b)*a - sin(b)*sin(c)*tan(b)*a + sin(b)*cos(c)*tan(b)

`rs\_series` can expand complicated multivariate expressions involving multiple
functions and most importantly, it does so blazingly fast::

    >>> %timeit ((sin(a) + cos(a))**10).series(a, 0, 5)
    1 loops, best of 3: 1.33 s per loop
    >>> %timeit rs_series((sin(a) + cos(a))**10, a, 5)
    100 loops, best of 3: 4.13 ms per loop

`rs\_series` is over 300 times faster. Given an expression to expand, there is
some fixed overhead to parse it. Thus, for larger orders, the speed
improvement becomes more prominent::

    >>> %timeit rs_series((sin(a) + cos(a))**10, a, 100)
    10 loops, best of 3: 32.8 ms per loop

Contribute
==========

`rs\_series` is not fully implemented yet. As of now, it supports only
multivariate Taylor expansions of expressions involving ``sin``, ``cos``,
``exp`` and ``tan``. Adding the remaining functions is not at all difficult and
they will be gradually added. If you are interested in helping, read the
comments in ``ring_series.py``. Currently, it does not support Puiseux series
(though the elementary functions do). This is expected to be fixed soon.

You can also add more functions to ``ring_series.py``. Only elementary
functions are supported currently.

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
