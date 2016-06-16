===============
Geometry Module
===============


Introduction
------------

The geometry module for SymPy allows one to create two-dimensional geometrical
entities, such as lines and circles, and query for information about these
entities. This could include asking the area of an ellipse, checking for
collinearity of a set of points, or finding the intersection between two lines.
The primary use case of the module involves entities with numerical values, but
it is possible to also use symbolic representations.

Available Entities
------------------

The following entities are currently available in the geometry module:

* ``Point``
* ``Line``, ``Ray``, ``Segment``
* ``Ellipse``, ``Circle``
* ``Polygon``, ``RegularPolygon``, ``Triangle``

Most of the work one will do will be through the properties and methods of
these entities, but several global methods exist:

* ``intersection(entity1, entity2)``
* ``are_similar(entity1, entity2)``
* ``convex_hull(points)``

For a full API listing and an explanation of the methods and their return
values please see the list of classes at the end of this document.

Example Usage
-------------

The following Python session gives one an idea of how to work with some of the
geometry module.

    >>> from sympy import *
    >>> from sympy.geometry import *
    >>> x = Point(0, 0)
    >>> y = Point(1, 1)
    >>> z = Point(2, 2)
    >>> zp = Point(1, 0)
    >>> Point.is_collinear(x, y, z)
    True
    >>> Point.is_collinear(x, y, zp)
    False
    >>> t = Triangle(zp, y, x)
    >>> t.area
    1/2
    >>> t.medians[x]
    Segment2D(Point2D(0, 0), Point2D(1, 1/2))
    >>> Segment(Point(1, S(1)/2), Point(0, 0))
    Segment2D(Point2D(0, 0), Point2D(1, 1/2))
    >>> m = t.medians
    >>> intersection(m[x], m[y], m[zp])
    [Point2D(2/3, 1/3)]
    >>> c = Circle(x, 5)
    >>> l = Line(Point(5, -5), Point(5, 5))
    >>> c.is_tangent(l) # is l tangent to c?
    True
    >>> l = Line(x, y)
    >>> c.is_tangent(l) # is l tangent to c?
    False
    >>> intersection(c, l)
    [Point2D(-5*sqrt(2)/2, -5*sqrt(2)/2), Point2D(5*sqrt(2)/2, 5*sqrt(2)/2)]

Intersection of medians
-----------------------
::

    >>> from sympy import symbols
    >>> from sympy.geometry import Point, Triangle, intersection

    >>> a, b = symbols("a,b", positive=True)

    >>> x = Point(0, 0)
    >>> y = Point(a, 0)
    >>> z = Point(2*a, b)
    >>> t = Triangle(x, y, z)

    >>> t.area
    a*b/2

    >>> t.medians[x]
    Segment2D(Point2D(0, 0), Point2D(3*a/2, b/2))

    >>> intersection(t.medians[x], t.medians[y], t.medians[z])
    [Point2D(a, b/3)]

An in-depth example: Pappus' Hexagon Theorem
--------------------------------------------

From Wikipedia ([WikiPappus]_):

  Given one set of collinear points `A`, `B`, `C`, and another set of collinear
  points `a`, `b`, `c`, then the intersection points `X`, `Y`, `Z` of line pairs `Ab` and
  `aB`, `Ac` and `aC`, `Bc` and `bC` are collinear.

::

    >>> from sympy import *
    >>> from sympy.geometry import *
    >>>
    >>> l1 = Line(Point(0, 0), Point(5, 6))
    >>> l2 = Line(Point(0, 0), Point(2, -2))
    >>>
    >>> def subs_point(l, val):
    ...    """Take an arbitrary point and make it a fixed point."""
    ...    t = Symbol('t', real=True)
    ...    ap = l.arbitrary_point()
    ...    return Point(ap.x.subs(t, val), ap.y.subs(t, val))
    ...
    >>> p11 = subs_point(l1, 5)
    >>> p12 = subs_point(l1, 6)
    >>> p13 = subs_point(l1, 11)
    >>>
    >>> p21 = subs_point(l2, -1)
    >>> p22 = subs_point(l2, 2)
    >>> p23 = subs_point(l2, 13)
    >>>
    >>> ll1 = Line(p11, p22)
    >>> ll2 = Line(p11, p23)
    >>> ll3 = Line(p12, p21)
    >>> ll4 = Line(p12, p23)
    >>> ll5 = Line(p13, p21)
    >>> ll6 = Line(p13, p22)
    >>>
    >>> pp1 = intersection(ll1, ll3)[0]
    >>> pp2 = intersection(ll2, ll5)[0]
    >>> pp3 = intersection(ll4, ll6)[0]
    >>>
    >>> Point.is_collinear(pp1, pp2, pp3)
    True

References
~~~~~~~~~~

.. [WikiPappus] "Pappus's Hexagon Theorem" Wikipedia, the Free Encyclopedia.
        Web. 26 Apr. 2013.
        <http://en.wikipedia.org/wiki/Pappus's_hexagon_theorem>

Miscellaneous Notes
-------------------

* The area property of ``Polygon`` and ``Triangle`` may return a positive or
  negative value, depending on whether or not the points are oriented
  counter-clockwise or clockwise, respectively. If you always want a
  positive value be sure to use the ``abs`` function.
* Although ``Polygon`` can refer to any type of polygon, the code has been
  written for simple polygons. Hence, expect potential problems if dealing
  with complex polygons (overlapping sides).
* Since SymPy is still in its infancy some things may not simplify
  properly and hence some things that should return ``True`` (e.g.,
  ``Point.is_collinear``) may not actually do so. Similarly, attempting to find
  the intersection of entities that do intersect may result in an empty
  result.

Future Work
-----------

Truth Setting Expressions
~~~~~~~~~~~~~~~~~~~~~~~~~

When one deals with symbolic entities, it often happens that an assertion
cannot be guaranteed. For example, consider the following code:

    >>> from sympy import *
    >>> from sympy.geometry import *
    >>> x,y,z = map(Symbol, 'xyz')
    >>> p1,p2,p3 = Point(x, y), Point(y, z), Point(2*x*y, y)
    >>> Point.is_collinear(p1, p2, p3)
    False

Even though the result is currently ``False``, this is not *always* true. If the
quantity `z - y - 2*y*z + 2*y**2 == 0` then the points will be collinear. It
would be really nice to inform the user of this because such a quantity may be
useful to a user for further calculation and, at the very least, being nice to
know. This could be potentially done by returning an object (e.g.,
GeometryResult) that the user could use. This actually would not involve an
extensive amount of work.

Three Dimensions and Beyond
~~~~~~~~~~~~~~~~~~~~~~~~~~~

Currently a limited subset of the geometry module has been extended to
three dimensions, but it certainly would be a good addition to extend
more. This would probably involve a fair amount of work since many of
the algorithms used are specific to two dimensions.

Geometry Visualization
~~~~~~~~~~~~~~~~~~~~~~

The plotting module is capable of plotting geometric entities. See
:ref:`Plotting Geometric Entities <plot_geom>` in
the plotting module entry.

Submodules
~~~~~~~~~~

.. toctree::
    :maxdepth: 3

    entities.rst
    utils.rst
    points.rst
    lines.rst
    curves.rst
    ellipses.rst
    polygons.rst
    plane.rst
