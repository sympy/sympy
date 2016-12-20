============================
Basic Implementation details
============================

Coordinate Systems and Vectors
==============================

As of now, :mod:`sympy.vector` only deals with the Cartesian (also called
rectangular) coordinate systems. A 3D Cartesian coordinate system can
be initialized in :mod:`sympy.vector` as

  >>> from sympy.vector import CoordSysCartesian
  >>> N = CoordSysCartesian('N')

The string parameter to the constructor denotes the name assigned to the
system, and will primarily be used for printing purposes.

Once a coordinate system (in essence, a ``CoordSysCartesian`` instance)
has been defined, we can access the orthonormal unit vectors (i.e. the
:math:`\mathbf{\hat{i}}`, :math:`\mathbf{\hat{j}}` and
:math:`\mathbf{\hat{k}}` vectors) and coordinate variables/base
scalars (i.e. the :math:`\mathbf{x}`, :math:`\mathbf{y}` and
:math:`\mathbf{z}` variables) corresponding to it. We will talk
about coordinate variables in the later sections.

The basis vectors for the :math:`X`, :math:`Y` and :math:`Z`
axes can be accessed using the ``i``, ``j`` and ``k``
properties respectively.

  >>> N.i
  N.i
  >>> type(N.i)
  <class 'sympy.vector.vector.BaseVector'>

As seen above, the basis vectors are all instances of a class called
``BaseVector``.

When a ``BaseVector`` is multiplied by a scalar (essentially any
SymPy ``Expr``), we get a ``VectorMul`` - the product of
a base vector and a scalar.

  >>> 3*N.i
  3*N.i
  >>> type(3*N.i)
  <class 'sympy.vector.vector.VectorMul'>

Addition of ``VectorMul`` and ``BaseVectors`` gives rise to
formation of ``VectorAdd`` - except for special cases, ofcourse.

  >>> v = 2*N.i + N.j
  >>> type(v)
  <class 'sympy.vector.vector.VectorAdd'>
  >>> v - N.j
  2*N.i
  >>> type(v - N.j)
  <class 'sympy.vector.vector.VectorMul'>

What about a zero vector? It can be accessed using the ``zero``
attribute assigned to class ``Vector``. Since the notion of a zero
vector remains the same regardless of the coordinate system in
consideration, we use ``Vector.zero`` wherever such a quantity is
required.

  >>> from sympy.vector import Vector
  >>> Vector.zero
  0
  >>> type(Vector.zero)
  <class 'sympy.vector.vector.VectorZero'>
  >>> N.i + Vector.zero
  N.i
  >>> Vector.zero == 2*Vector.zero
  True

All the classes shown above - ``BaseVector``, ``VectorMul``,
``VectorAdd`` and ``VectorZero`` are subclasses of ``Vector``.

You should never have to instantiate objects of any of the
subclasses of ``Vector``. Using the ``BaseVector`` instances assigned to a
``CoordSysCartesian`` instance and (if needed) ``Vector.zero``
as building blocks, any sort of vectorial expression can be constructed
with the basic mathematical operators ``+``, ``-``, ``*``.
and ``/``.

  >>> v = N.i - 2*N.j
  >>> v/3
  1/3*N.i + (-2/3)*N.j
  >>> v + N.k
  N.i + (-2)*N.j + N.k
  >>> Vector.zero/2
  0
  >>> (v/3)*4
  4/3*N.i + (-8/3)*N.j


In addition to the elementary mathematical operations, the vector
operations of ``dot`` and ``cross`` can also be performed on
``Vector``.

  >>> v1 = 2*N.i + 3*N.j - N.k
  >>> v2 = N.i - 4*N.j + N.k
  >>> v1.dot(v2)
  -11
  >>> v1.cross(v2)
  (-1)*N.i + (-3)*N.j + (-11)*N.k
  >>> v2.cross(v1)
  N.i + 3*N.j + 11*N.k

The ``&`` and ``^`` operators have been overloaded for the
``dot`` and ``cross`` methods respectively.

  >>> v1 & v2
  -11
  >>> v1 ^ v2
  (-1)*N.i + (-3)*N.j + (-11)*N.k

However, this is not the recommended way of performing these operations.
Using the original methods makes the code clearer and easier to follow.

In addition to these operations, it is also possible to compute the
outer products of ``Vector`` instances in :mod:`sympy.vector`. More
on that in a little bit.


SymPy operations on Vectors
===========================

The SymPy operations of ``simplify``, ``trigsimp``, ``diff``,
and ``factor`` work on ``Vector`` objects, with the standard SymPy API.

In essence, the methods work on the measure numbers(The coefficients
of the basis vectors) present in the provided vectorial expression.

  >>> from sympy.abc import a, b, c
  >>> from sympy import sin, cos, trigsimp, diff
  >>> v = (a*b + a*c + b**2 + b*c)*N.i + N.j
  >>> v.factor()
  ((a + b)*(b + c))*N.i + N.j
  >>> v = (sin(a)**2 + cos(a)**2)*N.i - (2*cos(b)**2 - 1)*N.k
  >>> trigsimp(v)
  N.i + (-cos(2*b))*N.k
  >>> v.simplify()
  N.i + (-cos(2*b))*N.k
  >>> diff(v, b)
  (4*sin(b)*cos(b))*N.k
  >>> from sympy import Derivative
  >>> Derivative(v, b).doit()
  (4*sin(b)*cos(b))*N.k

``Integral`` also works with ``Vector`` instances, similar to
``Derivative``.

  >>> from sympy import Integral
  >>> v1 = a*N.i + sin(a)*N.j - N.k
  >>> Integral(v1, a)
  (Integral(a, a))*N.i + (Integral(sin(a), a))*N.j + (Integral(-1, a))*N.k
  >>> Integral(v1, a).doit()
  a**2/2*N.i + (-cos(a))*N.j + (-a)*N.k

Points
======

As mentioned before, every coordinate system corresponds to a unique origin
point. Points, in general, have been implemented in :mod:`sympy.vector` in the
form of the ``Point`` class.

To access the origin of system, use the ``origin`` property of the
``CoordSysCartesian`` class.

  >>> from sympy.vector import CoordSysCartesian
  >>> N = CoordSysCartesian('N')
  >>> N.origin
  N.origin
  >>> type(N.origin)
  <class 'sympy.vector.point.Point'>

You can instantiate new points in space using the ``locate_new``
method of ``Point``. The arguments include the name(string) of the
new ``Point``, and its position vector with respect to the
'parent' ``Point``.

  >>> from sympy.abc import a, b, c
  >>> P = N.origin.locate_new('P', a*N.i + b*N.j + c*N.k)
  >>> Q = P.locate_new('Q', -b*N.j)

Like ``Vector``, a user never has to expressly instantiate an object of
``Point``. This is because any location in space (albeit relative) can be
pointed at by using the ``origin`` of a ``CoordSysCartesian`` as the
reference, and then using ``locate_new`` on it and subsequent
``Point`` instances.

The position vector of a ``Point`` with respect to another ``Point`` can
be computed using the ``position_wrt`` method.

  >>> P.position_wrt(Q)
  b*N.j
  >>> Q.position_wrt(N.origin)
  a*N.i + c*N.k

Additionally, it is possible to obtain the :math:`X`, :math:`Y` and :math:`Z`
coordinates of a ``Point`` with respect to a ``CoordSysCartesian``
in the form of a tuple. This is done using the ``express_coordinates``
method.

  >>> Q.express_coordinates(N)
  (a, 0, c)


Dyadics
=======

A dyadic, or dyadic tensor, is a second-order tensor formed by the
juxtaposition of pairs of vectors. Therefore, the outer products of vectors
give rise to the formation of dyadics. Dyadic tensors have been implemented
in :mod:`sympy.vector` in the ``Dyadic`` class.

Once again, you never have to instantiate objects of ``Dyadic``.
The outer products of vectors can be computed using the ``outer``
method of ``Vector``. The ``|`` operator has been overloaded for
``outer``.

  >>> from sympy.vector import CoordSysCartesian
  >>> N = CoordSysCartesian('N')
  >>> N.i.outer(N.j)
  (N.i|N.j)
  >>> N.i|N.j
  (N.i|N.j)

Similar to ``Vector``, ``Dyadic`` also has subsequent subclasses like
``BaseDyadic``, ``DyadicMul``, ``DyadicAdd``. As with ``Vector``,
a zero dyadic can be accessed from ``Dyadic.zero``.

All basic mathematical operations work with ``Dyadic`` too.

  >>> dyad = N.i.outer(N.k)
  >>> dyad*3
  3*(N.i|N.k)
  >>> dyad - dyad
  0
  >>> dyad + 2*(N.j|N.i)
  (N.i|N.k) + 2*(N.j|N.i)

``dot`` and ``cross`` also work among ``Dyadic`` instances as well as
between a ``Dyadic`` and ``Vector`` (and also vice versa) - as per the
respective mathematical definitions. As with ``Vector``, ``&`` and
``^`` have been overloaded for ``dot`` and ``cross``.

  >>> d = N.i.outer(N.j)
  >>> d.dot(N.j|N.j)
  (N.i|N.j)
  >>> d.dot(N.i)
  0
  >>> d.dot(N.j)
  N.i
  >>> N.i.dot(d)
  N.j
  >>> N.k ^ d
  (N.j|N.j)
