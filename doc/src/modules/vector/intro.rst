=======================================================
Introduction to Coordinate Systems, Vectors and Dyadics
=======================================================

A coordinate system is an abstract mathematical entity used to define
the notion of directions and locations in n-dimensional spaces. This
module deals with 3-dimensional spaces, with the conventional `\mathbf{X}`, 
`\mathbf{Y}` and `\mathbf{Z}` directions (or axes) defined with respect 
to each coordinate system.

Each coordinate system also has a special reference point called the 
'origin' defined for it. This point is used either while pointing to 
locations in 3D space, or while calculating the coordinates of 
pre-defined points with respect to the system in question.

As of now, :mod:`sympy.vector` only deals with the Cartesian (also called 
rectangular) coordinate systems. A 3D Cartesian coordinate system can
be initialized in :mod:`sympy.vector` as

  >>> from sympy.vector import CoordSysCartesian
  >>> N = CoordSysCartesian('N')

The string parameter to the constructor denotes the name assigned to the
system, and will primarily be used for printing purposes.

Once a coordinate system (in essence, a :mod:`CoordSysCartesian` instance)
has been defined, we can access the basis/base vectors (i.e. the 
`\mathbf{\hat{i}}`, `\mathbf{\hat{j}}` and `\mathbf{\hat{k}}` vectors) 
and coordinate variables/base scalars (i.e. the `\mathbf{x}`, 
`\mathbf{y}` and `\mathbf{z}` variables) corresponding to it. We will talk
about coordinate variables in the later section on fields.

The basis vectors for the `\mathbf{X}`, `\mathbf{Y}` and `\mathbf{Z}` 
directions can be accessed using the :mod:`i`, :mod:`j` and :mod:`k` 
:mod:`property`'s respectively.

  >>> N.i
  N.i
  >>> type(N.i)
  <class 'sympy.vector.vector.BaseVector'>

As seen above, the basis vectors are all instances of a class called 
:mod:`BaseVector`.

When a :mod:`BaseVector` is multiplied by a scalar (essentially any
:mod:`SymPy` :mod:`Expr`), we get a :mod:`VectorMul` - the product of
a base vector and a scalar.

  >>> 3*N.i
  3*N.i
  >>> type(3*N.i)
  <class 'sympy.vector.vector.VectorMul'>

Addition of :mod:`VectorMul`s and :mod:`BaseVectors`s gives rise to
formation of :mod:`VectorAdd`s - except for special cases, ofcourse.

  >>> v = 2*N.i + N.j
  >>> type(v)
  <class 'sympy.vector.vector.VectorAdd'>
  >>> v - N.j
  2*N.i
  >>> type(v - N.j)
  <class 'sympy.vector.vector.VectorMul'>

What about a zero vector? It can be accessed using the :mod:`zero`
attribute assigned to class :mod:`Vector`. Since the notion of a zero
vector remains the same regardless of the coordinate system in 
consideration, we use :mod:`Vector.zero` whereever such a quantity is
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

Few points worth noting about the :mod:`Vector` architecture in :mod:`sympy.vector`
-----------------------------------------------------------------------------------

1. All the classes shown above - :mod:`BaseVector`, :mod:`VectorMul`, 
:mod:`VectorAdd` and :mod:`VectorZero` are subclasses of :mod:`Vector`.

2. The user should never have to instantiate objects of any of the
subclasses of :mod:`Vector` - using the base vectors assigned to a
:mod:`CoordSysCartesian` instance and (if needed) :mod:`Vector.zero`
as building blocks, any sort of vectorial expression can be constructed
with the basic mathematical operators :mod:`+`, :mod:`-`, :mod:`*`
and :mod:`/`.

  >>> v = N.i - 2*N.j
  >>> v/3
  1/3*N.i + (-2/3)*N.j
  >>> v + N.k
  N.i + (-2)*N.j + N.k
  >>> Vector.zero/2
  0
  >>> (v/3)*4
  4/3*N.i + (-8/3)*N.j


In addition to the elementary mathematical operations, the vectorial 
operations of :mod:`dot` and :mod:`cross` can also be performed on 
:mod:`Vector`s.

  >>> 
