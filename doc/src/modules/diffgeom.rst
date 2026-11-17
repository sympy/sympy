=====================
Differential Geometry
=====================

.. module:: sympy.diffgeom

Introduction
------------

The module models differential geometry symbolically. A :class:`Manifold` carries
:class:`Patch` objects, a patch carries the :class:`CoordSystem` charts defined on it, and a
chart knows how its coordinates relate to those of another chart. Everything below is exact:
the transformations are SymPy expressions, not floating point routines.

Two charts on a two-dimensional patch, cartesian and polar, with the relation between them
given in both directions:

    >>> from sympy import symbols, sqrt, atan2, cos, sin
    >>> from sympy.diffgeom import Manifold, Patch, CoordSystem
    >>> m = Manifold('M', 2)
    >>> p = Patch('P', m)
    >>> x, y = symbols('x y', real=True)
    >>> r, theta = symbols('r theta', nonnegative=True)
    >>> relation_dict = {
    ... ('Car2D', 'Pol'): [(x, y), (sqrt(x**2 + y**2), atan2(y, x))],
    ... ('Pol', 'Car2D'): [(r, theta), (r*cos(theta), r*sin(theta))]
    ... }
    >>> Car2D = CoordSystem('Car2D', p, (x, y), relation_dict)
    >>> Pol = CoordSystem('Pol', p, (r, theta), relation_dict)

``transform`` gives the coordinates of the other chart, symbolically or at a point:

    >>> Car2D.transform(Pol)
    Matrix([
    [sqrt(x**2 + y**2)],
    [      atan2(y, x)]])
    >>> Car2D.transform(Pol, [1, 2])
    Matrix([
    [sqrt(5)],
    [atan(2)]])

Numerical evaluation
~~~~~~~~~~~~~~~~~~~~

To convert many points, turn the transformation into a numerical function with
:func:`~.lambdify` once and hand it whole arrays. Passing the entries as a tuple rather than
the matrix keeps the results one-dimensional:

    >>> import numpy
    >>> from sympy import lambdify
    >>> to_polar = lambdify(Car2D.symbols, tuple(Car2D.transform(Pol)), modules='numpy')
    >>> radius, angle = to_polar(numpy.array([1.0, 0.0, -1.0]),
    ...                         numpy.array([0.0, 2.0, 1.0]))
    >>> [round(float(value), 3) for value in radius]
    [1.0, 2.0, 1.414]
    >>> [round(float(value), 3) for value in angle]
    [0.0, 1.571, 2.356]

Call the function on the arrays themselves, as above, rather than once per point: the
vectorised call spends its time inside NumPy, while a Python loop over points pays the call
overhead for every one of them. See :ref:`numeric_computation` for the general picture.

Base Class Reference
--------------------
.. autoclass:: Manifold
   :members:

.. autoclass:: Patch
   :members:

.. autoclass:: CoordSystem
   :members:

.. autoclass:: CoordinateSymbol
   :members:

.. autoclass:: Point
   :members:

.. autoclass:: BaseScalarField
   :members:

.. autoclass:: BaseVectorField
   :members:

.. autoclass:: Commutator
   :members:

.. autoclass:: Differential
   :members:

.. autoclass:: TensorProduct
   :members:

.. autoclass:: WedgeProduct
   :members:

.. autoclass:: LieDerivative
   :members:

.. autoclass:: BaseCovarDerivativeOp
   :members:

.. autoclass:: CovarDerivativeOp
   :members:

.. autofunction:: intcurve_series

.. autofunction:: intcurve_diffequ

.. autofunction:: vectors_in_basis

.. autofunction:: twoform_to_matrix

.. autofunction:: metric_to_Christoffel_1st

.. autofunction:: metric_to_Christoffel_2nd

.. autofunction:: metric_to_Riemann_components

.. autofunction:: metric_to_Ricci_components
