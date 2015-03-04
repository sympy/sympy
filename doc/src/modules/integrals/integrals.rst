Symbolic Integrals
==================

.. module:: sympy.integrals

The ``integrals`` module in SymPy implements methods to calculate definite and indefinite integrals of expressions.

Principal method in this module is :func:`integrate`

  - ``integrate(f, x)`` returns the indefinite integral :math:`\int f\,dx`
  - ``integrate(f, (x, a, b))`` returns the definite integral :math:`\int_{a}^{b} f\,dx`

Examples
--------
SymPy can integrate a vast array of functions. It can integrate polynomial functions::

    >>> from sympy import *
    >>> init_printing(use_unicode=False, wrap_line=False, no_global=True)
    >>> x = Symbol('x')
    >>> integrate(x**2 + x + 1, x)
     3    2
    x    x
    -- + -- + x
    3    2

Rational functions::

	>>> integrate(x/(x**2+2*x+1), x)
	               1
	log(x + 1) + -----
	             x + 1


Exponential-polynomial functions. These multiplicative combinations of polynomials and the functions ``exp``, ``cos`` and ``sin`` can be integrated by hand using repeated integration by parts, which is an extremely tedious process. Happily, SymPy will deal with these integrals.

::

    >>> integrate(x**2 * exp(x) * cos(x), x)
     2  x           2  x                         x           x
    x *e *sin(x)   x *e *cos(x)      x          e *sin(x)   e *cos(x)
    ------------ + ------------ - x*e *sin(x) + --------- - ---------
         2              2                           2           2



even a few nonelementary integrals (in particular, some integrals involving the error function) can be evaluated::

	>>> integrate(exp(-x**2)*erf(x), x)
	  ____    2
	\/ pi *erf (x)
	--------------
	      4


Integral Transforms
-------------------

.. module:: sympy.integrals.transforms

SymPy has special support for definite integrals, and integral transforms.

.. autofunction:: mellin_transform
.. autofunction:: inverse_mellin_transform
.. autofunction:: laplace_transform
.. autofunction:: inverse_laplace_transform
.. autofunction:: fourier_transform
.. autofunction:: inverse_fourier_transform
.. autofunction:: sine_transform
.. autofunction:: inverse_sine_transform
.. autofunction:: cosine_transform
.. autofunction:: inverse_cosine_transform
.. autofunction:: hankel_transform
.. autofunction:: inverse_hankel_transform


Internals
---------

There is a general method for calculating antiderivatives of elementary functions, called the *Risch algorithm*. The Risch algorithm is a decision procedure that can determine whether an elementary solution exists, and in that case calculate it. It can be extended to handle many nonelementary functions in addition to the elementary ones.

SymPy currently uses a simplified version of the Risch algorithm, called the *Risch-Norman algorithm*. This algorithm is much faster, but may fail to find an antiderivative, although it is still very powerful. SymPy also uses pattern matching and heuristics to speed up evaluation of some types of integrals, e.g. polynomials.

For non-elementary definite integrals, SymPy uses so-called Meijer G-functions.
Details are described here:

.. toctree::
   :maxdepth: 1

   g-functions.rst

API reference
-------------

.. autofunction:: sympy.integrals.integrate
.. autofunction:: sympy.integrals.line_integrate
.. autofunction:: sympy.integrals.deltaintegrate
.. autofunction:: sympy.integrals.ratint
.. autofunction:: sympy.integrals.heurisch
.. autofunction:: sympy.integrals.heurisch_wrapper
.. autofunction:: sympy.integrals.trigintegrate
.. autofunction:: sympy.integrals.manualintegrate
.. autofunction:: sympy.integrals.manualintegrate.integral_steps

The class `Integral` represents an unevaluated integral and has some methods that help in the integration of an expression.

.. autoclass:: sympy.integrals.Integral
   :members:

   .. data:: is_commutative

      Returns whether all the free symbols in the integral are commutative.

TODO and Bugs
-------------
There are still lots of functions that SymPy does not know how to integrate. For bugs related to this module, see https://github.com/sympy/sympy/issues?labels=Integration

Numeric Integrals
=================

SymPy has functions to calculate points and weights for Gaussian quadrature of
any order and any precision:

.. autofunction:: sympy.integrals.quadrature.gauss_legendre

.. autofunction:: sympy.integrals.quadrature.gauss_laguerre

.. autofunction:: sympy.integrals.quadrature.gauss_hermite

.. autofunction:: sympy.integrals.quadrature.gauss_gen_laguerre

.. autofunction:: sympy.integrals.quadrature.gauss_chebyshev_t

.. autofunction:: sympy.integrals.quadrature.gauss_chebyshev_u

.. autofunction:: sympy.integrals.quadrature.gauss_jacobi
