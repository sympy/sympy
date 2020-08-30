================================
Applications of Vector Integrals
================================

To integrate a scalar or vector field over a region, we have to first define a region. SymPy provides three methods for defining a region. 

1. Using Parametric Equations.

2. Using Implicit Equation.

3. Using objects of geometry module.

The vector_integrate method is used to integrate scalar or vector field over any type of region. It automatically determines the type of integration (line, surface or volume) depending on the nature of the object.

Importing necessary classes and functions:

>>> from sympy.vector import pi
>>> from sympy.vector import ParametricRegion, vector_integrate
>>> from sympy.abc import r, x, y, z, theta, phi
>>> C = CoordSys3D('C')
    
Calculation of Perimeter, Surface Area and Volume
=================================================
To calculate the perimeter of circle, we need to define a circle, Lets define it using its parametric equation.

>>> circle = ParametricRegion((4*cos(theta), 4*sin(theta)), (theta, 0, 2*pi))

To calculate the perimeter of the circle, we can integrate it over a scalar field of unit magnitude.

>>> vector_integrate(1, circle)
8*pi

Suppose a user wants to calculate the perimeter of triangle. Determining the parametric representation of triangle can be difficult. Instead, the user can use an object of Polygon class in the geometry module.

>>> from sympy.geometry import Point, Polygon. 
>>> triangle = Polygon(Point(1, 2), (3, 5), (1,6))
>>> vector_integrate(1, triangle)
sqrt(5) + sqrt(13) + 4

To define a solid sphere, we need to use three parameters (r, theta and phi).

>>> sphere = ParametricRegion((r*sin(phi)*cos(theta),r*sin(phi)*sin(theta), r*cos(phi)),
... (theta, 0, pi), (phi, 0, pi), (r, 0, 3))
>>> vector_integrate(1, sphere)
-18*pi

Calculation of charge on a body
===============================


Calculation of Flux
===================


Verifying Stoke's Theorem
=========================


Example 1
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
    >>> f
    (-C.y**3/3 + sin(C.x))*C.i + (C.x**3/3 + cos(C.y))*C.j + C.x*C.y*C.z*C.k
    >>> vector_integrate(f,  circle)
    pi/2
    >>> vector_integrate(curl(f),  cone)
    pi/2


Verifying Divergence Theorem
============================
Example 1
    >>> sphere = ParametricRegion((4*sin(phi)*cos(theta),4*sin(phi)*sin(theta), 4*cos(phi)),
    ...     (theta, 0, 2*pi), (phi, 0, pi))
    >>> solidsphere = ParametricRegion((r*sin(phi)*cos(theta),r*sin(phi)*sin(theta), r*cos(phi)),
    ...     (r, 0, 4), (theta, 0, 2*pi), (phi, 0, pi))
    >>> field = C.x**3*C.i + C.y**3*C.j + C.z**3*C.k
    >>> vector_integrate(field, sphere)
    12288*pi/5
    >>> vector_integrate(divergence(field), solidsphere)
    -12288*pi/5

Example 2
    >>> cube = ParametricRegion((x, y, z), (x, 0, 1), (y, 0, 1), (z, 0, 1))
    >>> field = 2*C.x*C.y*C.i + 3*C.x*C.y*C.j + C.z*exp(C.x + C.y)*C.k
    >>> vector_integrate(divergence(field), cube)
    -E + 7/2 + E*(-1 + E)
