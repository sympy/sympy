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

    >>> from sympy.polys import ring, QQ, RR
    >>> from sympy.polys.ring_series import rs_sin
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

    >>> from sympy.polys.ring_series import rs_cos, rs_tan
    >>> R, x, y = ring('x, y', QQ)

    >>> rs_cos(x + x*y, x, 3)/x**3
    -1/2*x**(-1)*y**2 - x**(-1)*y - 1/2*x**(-1) + x**(-3)

    >>> rs_tan(x**QQ(2, 5)*y**QQ(1, 2), x, 2)
    1/3*x**(6/5)*y**(3/2) + x**(2/5)*y**(1/2)

By default, ``PolyElement`` did not allow non-natural numbers as exponents. It
converted a fraction to an integer and raised an error on getting negative
exponents. The goal of the ``ring series`` module is fast series expansion, and
not to use the ``polys`` module. The reason we use it as our backend is simply
because it implements a sparse representation and most of the basic functions
that we need. However, this default behaviour of ``polys`` was limiting for
``ring series``.

Note that there is no such constraint (in having rational exponents) in the
data-structure used by ``polys``- ``dict``. Sparse polynomials
(``PolyElement``) use the Python dict to store a polynomial term by term, where
a tuple of exponents is the key and the coefficient of that term is the value.
There is no reason why we can't have rational values in the ``dict`` so as to
support rational exponents.

So the approach we took was to modify sparse ``polys`` to allow non-natural
exponents. And it turned out to be quite simple. We only had to delete the
conversion to ``int`` of exponents in the ``__pow__`` method of
``PolyElement``. So::

    >>> x**QQ(3, 4)
    x**(3/4)

and not ``1`` as was the case earlier.

Though this change violates the definition of a polynomial, it doesn't break
anything yet.  Ideally, we shouldn't modify ``polys`` in any way. But to have
all the ``series`` capabilities we want, no other simple way was found. If need
be, we can separate the modified part of ``polys`` from core ``polys``. It
would be great if any other elegant solution is found.

All series returned by the functions of this module are instances of the
``PolyElement`` class. To use them with other SymPy types, convert them  to
``Expr``::

    >>> from sympy.polys.ring_series import rs_exp
    >>> from sympy.abc import a, b, c
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
is a polynomial over the simplest (almost) possible ring that does the job. It
recursively builds the ring as it parses the given expression, adding
generators to the ring when it needs them. Some examples::

    >>> rs_series(sin(a + b), a, 5) # doctest: +SKIP
    1/24*sin(b)*a**4 - 1/2*sin(b)*a**2 + sin(b) - 1/6*cos(b)*a**3 + cos(b)*a

    >>> rs_series(sin(exp(a*b) + cos(a + c)), a, 2) # doctest: +SKIP
    -sin(c)*cos(cos(c) + 1)*a + cos(cos(c) + 1)*a*b + sin(cos(c) + 1)

    >>> rs_series(sin(a + b)*cos(a + c)*tan(a**2 + b), a, 2) # doctest: +SKIP
    cos(b)*cos(c)*tan(b)*a - sin(b)*sin(c)*tan(b)*a + sin(b)*cos(c)*tan(b)

It can expand complicated multivariate expressions involving multiple functions
and most importantly, it does so blazingly fast::

    >>> %timeit ((sin(a) + cos(a))**10).series(a, 0, 5) # doctest: +SKIP
    1 loops, best of 3: 1.33 s per loop

    >>> %timeit rs_series((sin(a) + cos(a))**10, a, 5) # doctest: +SKIP
    100 loops, best of 3: 4.13 ms per loop

`rs\_series` is over 300 times faster. Given an expression to expand, there is
some fixed overhead to parse it. Thus, for larger orders, the speed
improvement becomes more prominent::

    >>> %timeit rs_series((sin(a) + cos(a))**10, a, 100) # doctest: +SKIP
    10 loops, best of 3: 32.8 ms per loop

To figure out the right ring for a given expression, `rs\_series` uses the
``sring`` function, which in turn uses other functions of ``polys``. As
explained above, non-natural exponents are not allowed. But the restriction is
on exponents and not generators. So, ``polys`` allows all sorts of symbolic
terms as generators to make sure that the exponent is a natural number::

    >>> from sympy.polys.rings import sring
    >>> R, expr = sring(1/a**3 + a**QQ(3, 7)); R
    Polynomial ring in 1/a, a**(1/7) over ZZ with lex order

In the above example, `1/a` and `a**(1/7)` will be treated as completely
different atoms. For all practical purposes, we could let `b = 1/a` and `c =
a**(1/7)` and do the manipulations. Effectively, expressions involving `1/a`
and `a**(1/7)` (and their powers) will never simplify::

    >>> expr*R(1/a) # doctest: +SKIP
    (1/a)**2 + (1/a)*(a**(1/7))**3

This leads to similar issues with manipulating Laurent and Puiseux series as
faced earlier. Fortunately, this time we have an elegant solution and are able
to isolate the ``series`` and ``polys`` behaviour from one another. We
introduce a boolean flag ``series`` in the list of allowed ``Options`` for
polynomials (see :class:`sympy.polys.polyoptions.Options`). Thus, when we want
``sring`` to allow rational exponents we supply a ``series=True`` flag to
``sring``::

    >>> rs_series(sin(a**QQ(1, 2)), a, 3) # doctest: +SKIP
    -1/5040*a**(7/3) + 1/120*a**(5/3) - 1/6*a + a**(1/3)

Contribute
==========

`rs\_series` is not fully implemented yet. As of now, it supports only
multivariate Taylor expansions of expressions involving ``sin``, ``cos``,
``exp`` and ``tan``. Adding the remaining functions is not at all difficult and
they will be gradually added. If you are interested in helping, read the
comments in ``ring_series.py``. Currently, it does not support Puiseux series
(though the elementary functions do). This is expected to be fixed soon.

You can also add more functions to ``ring_series.py``. Only elementary
functions are supported currently. The long term goal is to replace SymPy's
current ``series`` method with ``rs_series``.

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
