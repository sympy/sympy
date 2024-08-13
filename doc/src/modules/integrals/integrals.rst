=========
Integrals
=========

.. module:: sympy.integrals

The ``integrals`` module in SymPy implements methods to calculate definite and indefinite integrals of expressions.

Principal method in this module is :func:`~.integrate`

  - ``integrate(f, x)`` returns the indefinite integral :math:`\int f\,dx`
  - ``integrate(f, (x, a, b))`` returns the definite integral :math:`\int_{a}^{b} f\,dx`

Examples
--------
SymPy can integrate a vast array of functions. It can integrate polynomial functions::

    >>> from sympy import *
    >>> init_printing(use_unicode=False)
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
.. autoclass:: MellinTransform
   :members:
.. autofunction:: inverse_mellin_transform
.. autoclass:: InverseMellinTransform
   :members:
.. autofunction:: laplace_transform
.. autofunction:: laplace_correspondence
.. autofunction:: laplace_initial_conds
.. autoclass:: LaplaceTransform
   :members:
.. autofunction:: inverse_laplace_transform
.. autoclass:: InverseLaplaceTransform
   :members:
.. autofunction:: fourier_transform
.. autofunction:: _fourier_transform
.. autoclass:: FourierTransform
   :members:
.. autofunction:: inverse_fourier_transform
.. autoclass:: InverseFourierTransform
   :members:
.. autofunction:: sine_transform
.. autoclass:: SineTransform
   :members:
.. autofunction:: inverse_sine_transform
.. autoclass:: InverseSineTransform
   :members:
.. autofunction:: cosine_transform
.. autoclass:: CosineTransform
   :members:
.. autofunction:: inverse_cosine_transform
.. autoclass:: InverseCosineTransform
   :members:
.. autofunction:: hankel_transform
.. autoclass:: HankelTransform
   :members:
.. autofunction:: inverse_hankel_transform
.. autoclass:: InverseHankelTransform
   :members:
.. autoclass:: IntegralTransform
   :members:
.. autoexception:: IntegralTransformError

Internals
---------

SymPy uses a number of algorithms to compute integrals. Algorithms are tried
in order until one produces an answer. Most of these algorithms can be enabled
or disabled manually using various flags to :func:`~.integrate` or :meth:`~.Integral.doit`.

SymPy first applies several heuristic algorithms, as these are the fastest:

1. If the function is a rational function, there is a complete algorithm for
   integrating rational functions called the Lazard-Rioboo-Trager and the
   Horowitz-Ostrogradsky algorithms. They are implemented in :func:`.ratint`.

   .. autofunction:: sympy.integrals.rationaltools::ratint
   .. autofunction:: sympy.integrals.rationaltools::ratint_ratpart
   .. autofunction:: sympy.integrals.rationaltools::ratint_logpart

2. :func:`.trigintegrate` solves integrals of trigonometric functions using
   pattern matching

   .. autofunction:: sympy.integrals.trigonometry::trigintegrate

3. :func:`.deltaintegrate` solves integrals with :class:`~.DiracDelta` objects.

   .. autofunction:: sympy.integrals.deltafunctions::deltaintegrate

4. :func:`.singularityintegrate` is applied if the function contains a :class:`~.SingularityFunction`

   .. autofunction:: sympy.integrals.singularityfunctions::singularityintegrate

5. If the heuristic algorithms cannot be applied, :func:`.risch_integrate` is
   tried next. The *Risch algorithm* is a general method for calculating
   antiderivatives of elementary functions. The Risch algorithm is a decision
   procedure that can determine whether an elementary solution exists, and in
   that case calculate it. It can be extended to handle many nonelementary
   functions in addition to the elementary ones. However, the version implemented
   in SymPy only supports a small subset of the full algorithm, particularly, on
   part of the transcendental algorithm for exponentials and logarithms is
   implemented. An advantage of :func:`.risch_integrate` over other methods is
   that if it returns an instance of :class:`.NonElementaryIntegral`, the
   integral is proven to be nonelementary by the algorithm, meaning the integral
   cannot be represented using a combination of exponentials, logarithms, trig
   functions, powers, rational functions, algebraic functions, and function
   composition.

   .. autofunction:: sympy.integrals.risch::risch_integrate
   .. autoclass:: sympy.integrals.risch::NonElementaryIntegral
      :members:

6. For non-elementary definite integrals, SymPy uses so-called Meijer G-functions.
   Details are described in :ref:`g-functions`.

7. All the algorithms mentioned thus far are either pattern-matching based
   heuristic, or solve integrals using algorithms that are much different from
   the way most people are taught in their calculus courses. SymPy also
   implements a method that can solve integrals in much the same way you would in
   calculus. The advantage of this method is that it is possible to extract the
   integration steps from, so that one can see how to compute the integral "by
   hand". This is used by `SymPy Gamma <https://sympygamma.com>`_. This is
   implemented in the :func:`.manualintegrate` function. The steps for an integral
   can be seen with the :func:`.integral_steps` function.

   .. autofunction:: sympy.integrals.manualintegrate::manualintegrate
   .. autofunction:: sympy.integrals.manualintegrate::integral_steps

8. Finally, if all the above fail, SymPy also uses a simplified version of the
   Risch algorithm, called the *Risch-Norman algorithm*. This algorithm is tried
   last because it is often the slowest to compute. This is implemented in
   :func:`.heurisch`:

   .. autofunction:: sympy.integrals.heurisch::heurisch
   .. autofunction:: sympy.integrals.heurisch::components

API reference
-------------

.. autofunction:: sympy.integrals.integrals::integrate
.. autofunction:: sympy.integrals.integrals::line_integrate

The class :class:`~.Integral` represents an unevaluated integral and has some methods that help in the integration of an expression.

.. autoclass:: sympy.integrals.integrals::Integral
   :members:

   .. data:: is_commutative

      Returns whether all the free symbols in the integral are commutative.

:class:`~.Integral` subclasses from :class:`~.ExprWithLimits`, which is a
common superclass of :class:`~.Integral` and :class:`~.Sum`.

.. autoclass:: sympy.concrete.expr_with_limits::ExprWithLimits
   :members:

TODO and Bugs
-------------
There are still lots of functions that SymPy does not know how to integrate. For bugs related to this module, see https://github.com/sympy/sympy/issues?q=is%3Aissue+is%3Aopen+label%3Aintegrals

Numeric Integrals
-----------------

SymPy has functions to calculate points and weights for Gaussian quadrature of
any order and any precision:

.. autofunction:: sympy.integrals.quadrature::gauss_legendre

.. autofunction:: sympy.integrals.quadrature::gauss_laguerre

.. autofunction:: sympy.integrals.quadrature::gauss_hermite

.. autofunction:: sympy.integrals.quadrature::gauss_gen_laguerre

.. autofunction:: sympy.integrals.quadrature::gauss_chebyshev_t

.. autofunction:: sympy.integrals.quadrature::gauss_chebyshev_u

.. autofunction:: sympy.integrals.quadrature::gauss_jacobi

.. autofunction:: sympy.integrals.quadrature::gauss_lobatto

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
    >>> init_printing(use_unicode=False)
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

.. autofunction:: sympy.integrals.intpoly::polytope_integrate
