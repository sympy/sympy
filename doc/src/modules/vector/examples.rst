=========================
General examples of usage
=========================

This section details the solution of two basic problems in vector
math/calculus using the :mod:`sympy.vector` package.

Quadrilateral problem
=====================

The Problem
-----------

*OABC is any quadrilateral in 3D space. P is the
midpoint of OA, Q is the midpoint of AB, R is the midpoint of BC
and S is the midpoint of OC. Prove that PQ is parallel to SR*

Solution
--------

The solution to this problem demonstrates the usage of ``Point``,
and basic operations on ``Vector``.

Define a coordinate system

  >>> from sympy.vector import CoordSysCartesian
  >>> Sys = CoordSysCartesian('Sys')

Define point O to be Sys' origin. We can do this without
loss of generality

  >>> O = Sys.origin

Define point A with respect to O

  >>> from sympy import symbols
  >>> a1, a2, a3 = symbols('a1 a2 a3')
  >>> A = O.locate_new('A', a1*Sys.i + a2*Sys.j + a3*Sys.k)

Similarly define points B and C

  >>> b1, b2, b3 = symbols('b1 b2 b3')
  >>> B = O.locate_new('B', b1*Sys.i + b2*Sys.j + b3*Sys.k)
  >>> c1, c2, c3 = symbols('c1 c2 c3')
  >>> C = O.locate_new('C', c1*Sys.i + c2*Sys.j + c3*Sys.k)

P is the midpoint of OA. Lets locate it with respect to O
(you could also define it with respect to A).

  >>> P = O.locate_new('P', A.position_wrt(O) + (O.position_wrt(A) / 2))

Similarly define points Q, R and S as per the problem definitions.

  >>> Q = A.locate_new('Q', B.position_wrt(A) / 2)
  >>> R = B.locate_new('R', C.position_wrt(B) / 2)
  >>> S = O.locate_new('R', C.position_wrt(O) / 2)

Now compute the vectors in the directions specified by PQ and SR.

  >>> PQ = Q.position_wrt(P)
  >>> SR = R.position_wrt(S)

Compute cross product

  >>> PQ.cross(SR)
  0

Since the cross product is a zero vector, the two vectors have to be
parallel, thus proving that PQ || SR.


Third product rule for Del operator
===================================

See
---

.. [WikiDel] http://en.wikipedia.org/wiki/Del

The Problem
-----------

Prove the third rule -
:math:`\nabla \cdot (f \vec v) = f (\nabla \cdot \vec v) + \vec v \cdot (\nabla f)`

Solution
--------

Start with a coordinate system

  >>> from sympy.vector import CoordSysCartesian
  >>> C = CoordSysCartesian('C')

The scalar field :math:`f` and the measure numbers of the vector field
:math:`\vec v` are all functions of the coordinate variables of the
coordinate system in general.
Hence, define SymPy functions that way.

  >>> from sympy import symbols
  >>> v1, v2, v3, f = symbols('v1 v2 v3 f', type="Function")

``v1``, ``v2`` and ``v3`` are the :math:`X`, :math:`Y` and :math:`Z`
components of the vector field respectively.

Define the vector field as ``vfield`` and the scalar field as ``sfield``.

  >>> vfield = v1(C.x, C.y, C.z)*C.i + v2(C.x, C.y, C.z)*C.j + v3(C.x, C.y, C.z)*C.k
  >>> ffield = f(C.x, C.y, C.z)

Construct the expression for the LHS of the equation using ``C.delop``.

  >>> lhs = (C.delop.dot(ffield * vfield)).doit()

Similarly, the RHS would be defined.

  >>> rhs = ((vfield.dot(C.delop(ffield))) + (ffield * (C.delop.dot(vfield)))).doit()

Now, to prove the product rule, we would just need to equate the expanded and
simplified versions of the lhs and the rhs, so that the SymPy expressions match.

  >>> lhs.expand().simplify() == rhs.expand().simplify()
  True

Thus, the general form of the third product rule mentioned above can be proven
using :mod:`sympy.vector`.
