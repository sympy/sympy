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

The class `Integral` represents an unevaluated integral and has some methods that help in the integration of an expression.

.. autoclass:: sympy.integrals.Integral
   :members:

   .. data:: is_commutative

      Returns whether all the free symbols in the integral are commutative.

TODO and Bugs
-------------
There are still lots of functions that SymPy does not know how to integrate. For bugs related to this module, see https://github.com/sympy/sympy/issues?q=label%3AIntegration

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

Integration over Polytopes
==========================

.. module:: sympy.integrals.intpoly

The ``intpoly`` module in SymPy implements methods to calculate the integral of a polynomial over 2/3-Polytopes.
Uses evaluation techniques as described in Chin et al. (2015) [1].

The input for 2-Polytope or Polygon uses the already existing ``Polygon`` data structure in SymPy. See
:mod:`sympy.geometry.polygon` for how to create a polygon.

For the 3-Polytope or Polyhedron, the most economical representation
is to specify a list of vertices and then to provide each constituting face(Polygon) as a list of vertex indices.

For example, consider the unit cube. Here is how it would be represented.

``unit_cube = [[(0, 0, 0), (0, 0, 1), (0, 1, 0), (0, 1, 1), (1, 0, 0),(1, 0, 1), (1, 1, 0), (1, 1, 1)],``
            ``[3, 7, 6, 2], [1, 5, 7, 3], [5, 4, 6, 7], [0, 4, 5, 1], [2, 0, 1, 3], [2, 6, 4, 0]]``

Here, the first sublist is the list of vertices. The other smaller lists such as ``[3, 7, 6, 2]`` represent a 2D face
of the polyhedra with vertices having index ``3, 7, 6 and 2`` in the first sublist(in that order).

Principal method in this module is :func:`polytope_integrate`

  - ``polytope_integrate(Polygon((0, 0), (0, 1), (1, 0)), x)`` returns the integral of :math:`x` over the triangle with vertices (0, 0), (0, 1) and (1, 0)

  - ``polytope_integrate(unit_cube, x + y + z)`` returns the integral of :math:`x + y + z` over the unit cube.

References
----------
[1] : Chin, Eric B., Jean B. Lasserre, and N. Sukumar. "Numerical integration of homogeneous functions on convex and nonconvex polygons and polyhedra." Computational Mechanics 56.6 (2015): 967-981

PDF link : http://dilbert.engr.ucdavis.edu/~suku/quadrature/cls-integration.pdf

Examples
--------

For 2D Polygons
---------------

Single Polynomial::

    >>> from sympy.integrals.intpoly import *
    >>> init_printing(use_unicode=False, wrap_line=False, no_global=True)
    >>> polytope_integrate(Polygon((0, 0), (0, 1), (1, 0)), x)
    1/6
    >>> polytope_integrate(Polygon((0, 0), (0, 1), (1, 0)), x + x*y + y**2)
    7/24

List of specified polynomials::

    >>> polytope_integrate(Polygon((0, 0), (0, 1), (1, 0)), [3, x*y + y**2, x**4], max_degree=4)
              4               2
    {3: 3/2, x : 1/30, x*y + y : 1/8}
    >>> polytope_integrate(Polygon((0, 0), (0, 1), (1, 0)), [1.125, x, x**2, 6.89*x**3, x*y + y**2, x**4], max_degree=4)
                           2              3  689    4               2
    {1.125: 9/16, x: 1/6, x : 1/12, 6.89*x : ----, x : 1/30, x*y + y : 1/8}
                                             2000

Computing all monomials up to a maximum degree::

    >>> polytope_integrate(Polygon((0, 0), (0, 1), (1, 0)),max_degree=3)
                            2         3                 2         3                      2         2
    {0: 0, 1: 1/2, x: 1/6, x : 1/12, x : 1/20, y: 1/6, y : 1/12, y : 1/20, x*y: 1/24, x*y : 1/60, x *y: 1/60}


For 3-Polytopes/Polyhedra
-------------------------

Single Polynomial::

    >>> from sympy.integrals.intpoly import *
    >>> cube = [[(0, 0, 0), (0, 0, 5), (0, 5, 0), (0, 5, 5), (5, 0, 0), (5, 0, 5), (5, 5, 0), (5, 5, 5)], [2, 6, 7, 3], [3, 7, 5, 1], [7, 6, 4, 5], [1, 5, 4, 0], [3, 1, 0, 2], [0, 4, 6, 2]]
    >>> polytope_integrate(cube, x**2 + y**2 + z**2 + x*y + y*z + x*z)
    -21875/4
    >>> octahedron = [[(S(-1) / sqrt(2), 0, 0), (0, S(1) / sqrt(2), 0), (0, 0, S(-1) / sqrt(2)), (0, 0, S(1) / sqrt(2)), (0, S(-1) / sqrt(2), 0), (S(1) / sqrt(2), 0, 0)], [3, 4, 5], [3, 5, 1], [3, 1, 0], [3, 0, 4], [4, 0, 2], [4, 2, 5], [2, 0, 1], [5, 2, 1]]
    >>> polytope_integrate(octahedron, x**2 + y**2 + z**2 + x*y + y*z + x*z)
      ___
    \/ 2
    -----
      20

List of specified polynomials::

    >>> polytope_integrate(Polygon((0, 0), (0, 1), (1, 0)), [3, x*y + y**2, x**4], max_degree=4)
              4               2
    {3: 3/2, x : 1/30, x*y + y : 1/8}
    >>> polytope_integrate(Polygon((0, 0), (0, 1), (1, 0)), [1.125, x, x**2, 6.89*x**3, x*y + y**2, x**4], max_degree=4)
                           2              3  689    4               2
    {1.125: 9/16, x: 1/6, x : 1/12, 6.89*x : ----, x : 1/30, x*y + y : 1/8}
                                             2000

Computing all monomials up to a maximum degree::

    >>> polytope_integrate(Polygon((0, 0), (0, 1), (1, 0)),max_degree=3)
                            2         3                 2         3                      2         2
    {0: 0, 1: 1/2, x: 1/6, x : 1/12, x : 1/20, y: 1/6, y : 1/12, y : 1/20, x*y: 1/24, x*y : 1/60, x *y: 1/60}

API reference
-------------

.. autofunction:: sympy.integrals.intpoly.polytope_integrate
