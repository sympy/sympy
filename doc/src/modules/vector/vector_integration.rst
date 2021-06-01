================================
Applications of Vector Integrals
================================

To integrate a scalar or vector field over a region, we have to first define a region. SymPy provides three methods for defining a region:

1. Using Parametric Equations with :class:`~sympy.vector.parametricregion.ParametricRegion`.

2. Using Implicit Equation with :class:`~sympy.vector.implicitregion.ImplicitRegion`.

3. Using objects of geometry module.

The :func:`~sympy.vector.integrals.vector_integrate` function is used to integrate scalar or vector field over any type of region. It automatically determines the type of integration (line, surface, or volume) depending on the nature of the object.

We define a coordinate system and make necesssary imports for examples.

>>> from sympy import sin, cos, exp, pi, symbols
>>> from sympy.vector import CoordSys3D, ParametricRegion, ImplicitRegion, vector_integrate
>>> from sympy.abc import r, x, y, z, theta, phi
>>> C = CoordSys3D('C')

Calculation of Perimeter, Surface Area, and Volume
==================================================

To calculate the perimeter of a circle, we need to define it. Let's define it using its parametric equation.

>>> param_circle = ParametricRegion((4*cos(theta), 4*sin(theta)), (theta, 0, 2*pi))

We can also define a circle using its implicit equation.

>>> implicit_circle = ImplicitRegion((x, y), x**2 + y**2 - 4)

The perimeter of a figure is equal to the absolute value of its integral over a unit scalar field.

>>> vector_integrate(1, param_circle)
8*pi
>>> vector_integrate(1, implicit_circle)
4*pi

Suppose a user wants to calculate the perimeter of a triangle. Determining the parametric representation of a triangle can be difficult. Instead, the user can use an object of :class:`~sympy.geometry.polygon.Polygon` class in the geometry module.

>>> from sympy.geometry import Point, Polygon
>>> triangle = Polygon(Point(1, 2), (3, 5), (1,6))
>>> vector_integrate(1, triangle)
sqrt(5) + sqrt(13) + 4

To define a solid sphere, we need to use three parameters (r, theta and phi). For :class:`~sympy.vector.parametricregion.ParametricRegion` obextj, the order of limits determine the sign of the integral.

>>> solidsphere = ParametricRegion((r*sin(phi)*cos(theta),r*sin(phi)*sin(theta), r*cos(phi)),
...                             (phi, 0, pi), (theta, 0, 2*pi), (r, 0, 3))
>>> vector_integrate(1, solidsphere)
36*pi

Calculation of mass of a body
=============================

Consider a triangular lamina 𝑅  with vertices (0,0), (0, 5), (5,0) and with density :math:`\rho(x, y) = xy\:kg/m^2`. Find the total mass.

>>> triangle = ParametricRegion((x, y), (x, 0, 5), (y, 0, 5 - x))
>>> vector_integrate(C.x*C.y, triangle)
625/24

Find the mass of a cylinder centered on the z-axis which has height h, radius a, and density :math:`\rho = x^2 + y^2\:kg/m^2`.

>>> a, h = symbols('a h', positive=True)
>>> cylinder = ParametricRegion((r*cos(theta), r*sin(theta), z),
...                     (theta, 0, 2*pi), (z, 0, h), (r, 0, a))
>>> vector_integrate(C.x**2 + C.y**2, cylinder)
pi*a**4*h/2

Calculation of Flux
===================

1. Consider a region of space in which there is a constant vectorfield
:math:`E(x, y, z) = a\mathbf{\hat{k}}`.
A  hemisphere of radius r  lies on the x-y plane. What is the flux of the field through the sphere?

>>> semisphere = ParametricRegion((r*sin(phi)*cos(theta), r*sin(phi)*sin(theta), r*cos(phi)),\
...                             (phi, 0, pi/2), (theta, 0, 2*pi))
>>> flux = vector_integrate(a*C.k, semisphere)
>>> flux
pi*a*r**2

2. Consider  a  region  of  space  in  which  there  is  a  vector  field
:math:`E(x, y, z) = x^2 \mathbf{\hat{k}}` above the x-y plane, and a field
:math:`E(x, y, z) = y^2 \mathbf{\hat{k}}` below the x-y plane. What is the flux of that vector field through a cube of side length L with its center at the origin?”

The field is parallel to the z-axis so only the top and bottom face of the box will contribute to flux.

>>> L = symbols('L', positive=True)
>>> top_face = ParametricRegion((x, y, L/2), (x, -L/2, L/2), (y, -L/2, L/2))
>>> bottom_face = ParametricRegion((x, y, -L/2), (x, -L/2, L/2), (y, -L/2, L/2))
>>> flux = vector_integrate(C.x**2*C.k, top_face) + vector_integrate(C.y**2*C.k, bottom_face)
>>> flux
L**4/6

Verifying Stoke's Theorem
=========================

See https://en.wikipedia.org/wiki/Stokes%27_theorem

Example 1
    >>> from sympy.vector import curl
    >>> curve = ParametricRegion((cos(theta), sin(theta)), (theta, 0, pi/2))
    >>> surface = ParametricRegion((r*cos(theta), r*sin(theta)), (r, 0, 1), (theta, 0, pi/2))
    >>> F = C.y*C.i + C.z*C.k + C.x*C.k
    >>>
    >>> vector_integrate(F, curve)
    -pi/4
    >>> vector_integrate(curl(F), surface)
    -pi/4

Example 2
    >>> circle = ParametricRegion((cos(theta), sin(theta), 1), (theta, 0, 2*pi))
    >>> cone = ParametricRegion((r*cos(theta), r*sin(theta), r), (r, 0, 1), (theta, 0, 2*pi))
    >>> cone = ParametricRegion((r*cos(theta), r*sin(theta), r), (r, 0, 1), (theta, 0, 2*pi))
    >>> f = (-C.y**3/3 + sin(C.x))*C.i + (C.x**3/3 + cos(C.y))*C.j + C.x*C.y*C.z*C.k
    >>> vector_integrate(f,  circle)
    pi/2
    >>> vector_integrate(curl(f),  cone)
    pi/2


Verifying Divergence Theorem
============================

See https://en.wikipedia.org/wiki/Divergence_theorem

Example 1
    >>> from sympy.vector import divergence
    >>> sphere = ParametricRegion((4*sin(phi)*cos(theta),4*sin(phi)*sin(theta), 4*cos(phi)),
    ...                         (phi, 0, pi), (theta, 0, 2*pi))
    >>> solidsphere = ParametricRegion((r*sin(phi)*cos(theta),r*sin(phi)*sin(theta), r*cos(phi)),
    ...     (r, 0, 4),(phi, 0, pi), (theta, 0, 2*pi))
    >>> field = C.x**3*C.i + C.y**3*C.j + C.z**3*C.k
    >>> vector_integrate(field, sphere)
    12288*pi/5
    >>> vector_integrate(divergence(field), solidsphere)
    12288*pi/5

Example 2
    >>> cube = ParametricRegion((x, y, z), (x, 0, 1), (y, 0, 1), (z, 0, 1))
    >>> field = 2*C.x*C.y*C.i + 3*C.x*C.y*C.j + C.z*exp(C.x + C.y)*C.k
    >>> vector_integrate(divergence(field), cube)
    -E + 7/2 + E*(-1 + E)
