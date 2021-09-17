.. _matrices-random:

Generate Specific Random Matrices
=================================

In order to setup distinct mathematical problems, like

 - systems of linear equations,
 - invertible matrices,
 - triangular matrices
 - differential equations with given eigenvalues
 - linear isometries
 - non-standard geometries by non-standard scalar products
 - orthogonal basis of vector space,

specific types of matrices are required.

A purely random generation of matrix entries does not grant any
specific shape or property of the matrix.
So, those matrices have to be created carefully using
classical classification theorems of linear algebra.

Generators on different types of matrices are presented below.

.. note::

   Even for carefully sampled matrices which are elements in topological matrix groups,
   the sample distribution does not claim to meet any uniform distribution
   (with respect to the groups Haar measure).
   So the sampling procedures may be used to produce educational or academic problems
   but do not generate proper pre-defined distribution.
   These matrices can be useful for testing algorithms

There are three different types of matrix generators:

   - the *base matrices* which form the atomic building block of the later
   - the *compound matrices* which are simply products of base matrices of given types
   - the *conjugate matrices* which are usually given as base or compound matrix $\mathbf{A}$
     conjugate by another compound matrix $\mathbf{S}$,  i.e.

      .. math::

         \mathbf{B} = \mathbf{S} \cdot \mathbf{A} \cdot \mathbf{S}^{-1}
         \text{ or }
         \mathbf{C} = \mathbf{S} \cdot \mathbf{A} \cdot \mathbf{S}^{T}


 1. *base matrices* are simple matrix types with only a few non trivial entries.
    Those entries can set explicit by arguments or, if not given, are chosen randomly.

    - :py:func:`jordan <sympy.matrices.random.jordan>`
    - :py:func:`elementary <sympy.matrices.random.elementary>`

    - :py:func:`transposition <sympy.matrices.random.transposition>`
    - :py:func:`permutation <sympy.matrices.random.permutation>`

    - :py:func:`projection <sympy.matrices.random.projection>`

    - :py:func:`rotation <sympy.matrices.random.rotation>`
    - :py:func:`reflection <sympy.matrices.random.reflection>`

    as well as specific *normal form* matrices defined by their spectrum of eigenvalues
    given via the **spec** argument.

    - :py:func:`diagonal_normal <sympy.matrices.random.diagonal_normal>`
    - :py:func:`jordan_normal <sympy.matrices.random.jordan_normal>`
    - :py:func:`isometry_normal <sympy.matrices.random.isometry_normal>`

 2. *compound matrices* are build as a multiplication of several base matrices.
    Since by this the entries of the base matrices are chosen randomly
    but can be controlled by **domain** arguments.

    - :py:func:`permutation <sympy.matrices.random.permutation>` as product of
      :py:func:`transposition <sympy.matrices.random.transposition>` matrices

    - :py:func:`orthogonal <sympy.matrices.random.orthogonal>` as product of
      :py:func:`rotation <sympy.matrices.random.rotation>` and
      :py:func:`reflection <sympy.matrices.random.reflection>` matrices

    - :py:func:`invertible <sympy.matrices.random.invertible>` as product of
      :py:func:`elementary <sympy.matrices.random.elementary>` matrices

    - :py:func:`triangular <sympy.matrices.random.triangular>` as product of *triangular*
      :py:func:`elementary <sympy.matrices.random.elementary>` matrices

 3. *conjugate matrices* that arise as a product $\mathbf{SAS}^{-1}$
    of specific matrices $\mathbf{A}$
    - usually given as a *normal form* - conjugate by an
    :py:func:`invertible <sympy.matrices.random.invertible>` $\mathbf{S}$
    such as

    - :py:func:`diagonalizable <sympy.matrices.random.diagonalizable>` with
      :py:func:`diagonal_normal <sympy.matrices.random.diagonal_normal>`

    - :py:func:`idempotent <sympy.matrices.random.idempotent>` with
      :py:func:`projection <sympy.matrices.random.projection>` (with only zero and one entries)

    - :py:func:`trigonalizable <sympy.matrices.random.trigonalizable>` with
      :py:func:`jordan_normal <sympy.matrices.random.jordan_normal>`

    - :py:func:`nilpotent <sympy.matrices.random.nilpotent>` with
      :py:func:`jordan_normal <sympy.matrices.random.jordan_normal>` (with zero eigenvalues)

    *isometry conjugate matrices* that arise as a product $\mathbf{OAO}^{-1} = \mathbf{OAO}^{T}$
    of specific matrices $A$
    where $\mathbf{A}$ sets the (complex) eigenvalue spectrum.

    - :py:func:`orthogonal <sympy.matrices.random.orthogonal>`
      can be a conjugate of
      :py:func:`isometry_normal <sympy.matrices.random.isometry_normal>`
      by
      :py:func:`orthogonal <sympy.matrices.random.orthogonal>`

    - same holds for :py:func:`unitary <sympy.matrices.random.unitary>`
      which can be a conjugate of
      :py:func:`isometry_normal <sympy.matrices.random.isometry_normal>`
      by
      :py:func:`unitary <sympy.matrices.random.unitary>`

    - :py:func:`normal <sympy.matrices.random.normal>`
      is a conjugate of
      :py:func:`diagonal_normal <sympy.matrices.random.diagonal_normal>`
      by either
      :py:func:`orthogonal <sympy.matrices.random.orthogonal>`
      or
      :py:func:`unitary <sympy.matrices.random.unitary>`

    finally *symmetric matrices* given as a as product $\mathbf{MM}^T$
    resp. $\mathbf{MM}^H$ of $\mathbf{M}$, such as

    - :py:func:`symmetric <sympy.matrices.random.symmetric>` $\mathbf{MM}^T$ with
      :py:func:`invertible <sympy.matrices.random.invertible>` $\mathbf{M}$

    - :py:func:`hermite <sympy.matrices.random.hermite>` $\mathbf{MM}^H$ with
      :py:func:`invertible <sympy.matrices.random.orthogonal>` $\mathbf{M}$


In addition to the type of matrix, also the type of entries (as a commutative ring with one)
to be able to define values (**domain**) can be specified,
from which the value entries (**scalar**) of the basic matrices are randomly generated.

Because the complexity and amount of entries in the generated compound matrices
in addition to the **domain**, also by the number of base matrices multiplied for generation
This can be set using the argument **length**.

Normal forms as well as conjugate types have the arg **spec** to provide a spectrum of eigenvalues.

The product of *elementary* matrices

    >>> from sympy.matrices.random import elementary

    >>> A = elementary(3, (0, 2), -2)
    >>> A
    Matrix([
    [1, 0, -2],
    [0, 1,  0],
    [0, 0,  1]])

    >>> B = elementary(3, (1, 2), 1)
    >>> B
    Matrix([
    [1, 0, 0],
    [0, 1, 1],
    [0, 0, 1]])

    >>> C = elementary(3, (1, 0), -1)
    >>> C
    Matrix([
    [ 1, 0, 0],
    [-1, 1, 0],
    [ 0, 0, 1]])

gives an *invertible* matrix $\mathbf{M}$

    >>> M = A*B*C
    >>> M
    Matrix([
    [ 1, 0, -2],
    [-1, 1,  1],
    [ 0, 0,  1]])

    >>> M.inv()
    Matrix([
    [1, 0, 2],
    [1, 1, 1],
    [0, 0, 1]])


Seeding Random Number Generators
--------------------------------

There are Three different ways to generate randomly matrix entries.

1. By default the module uses its own random generator ``rand`` which is a instance of ``random.Random()``.
Hence, the state can be set by invoking ``rand.seed()``.

    >>> from sympy.matrices.random import rand, orthogonal

    >>> rand.seed(1)
    >>> orthogonal(3)
    Matrix([
    [-1,          0,          0],
    [ 0,  sqrt(2)/2, -sqrt(2)/2],
    [ 0, -sqrt(2)/2, -sqrt(2)/2]])

    >>> rand.seed(2)
    >>> orthogonal(3)
    Matrix([
    [0, -1,  0],
    [0,  0, -1],
    [1,  0,  0]])

    >>> rand.seed(1)
    >>> orthogonal(3)
    Matrix([
    [-1,          0,          0],
    [ 0,  sqrt(2)/2, -sqrt(2)/2],
    [ 0, -sqrt(2)/2, -sqrt(2)/2]])

2. Provide a custom seed
as in ``random.Random().seed()``

    >>> from sympy.matrices.random import orthogonal
    >>> orthogonal(3, seed=1)
    Matrix([
    [-1,          0,          0],
    [ 0,  sqrt(2)/2, -sqrt(2)/2],
    [ 0, -sqrt(2)/2, -sqrt(2)/2]])

    >>> orthogonal(3, seed=2)
    Matrix([
    [0, -1,  0],
    [0,  0, -1],
    [1,  0,  0]])

    >>> orthogonal(3, seed=1)
    Matrix([
    [-1,          0,          0],
    [ 0,  sqrt(2)/2, -sqrt(2)/2],
    [ 0, -sqrt(2)/2, -sqrt(2)/2]])

3. Provide a custom *random number generator* which must implement the same instance method ``sample``
as ``random.Random()`` does.

    >>> from random import Random
    >>> from sympy.matrices.random import orthogonal
    >>> rng = Random()

    >>> rng.seed(1)
    >>> orthogonal(3, seed=rng)
    Matrix([
    [-1,          0,          0],
    [ 0,  sqrt(2)/2, -sqrt(2)/2],
    [ 0, -sqrt(2)/2, -sqrt(2)/2]])

    >>> rng.seed(2)
    >>> orthogonal(3, seed=rng)
    Matrix([
    [0, -1,  0],
    [0,  0, -1],
    [1,  0,  0]])

    >>> rng.seed(1)
    >>> orthogonal(3, seed=rng)
    Matrix([
    [-1,          0,          0],
    [ 0,  sqrt(2)/2, -sqrt(2)/2],
    [ 0, -sqrt(2)/2, -sqrt(2)/2]])

4. For sampling matrix entries from ``spec``, ``domain`` or ``units`` arguments
all these can be *random number generator*, too.
At least they should provide some method ``sample``
drawing a list of **k** random items (*without putting back*) via ``sample(k)``.

    >>> class Domain(list):
    ...     rng = Random()
    ...     def _sample(self, k):
    ...         return self.rng.sample(self, k)

    >>> from sympy import sqrt
    >>> from sympy.matrices.random import rand, orthogonal
    >>> domain = Domain([(sqrt(2) / 2, sqrt(2) / 2), (0, -1)])

    >>> # to obtain the same result, seed must use the same random number generator
    >>> # since seed is used for drawing indices, too.
    >>> domain.rng.seed(1)
    >>> orthogonal(3, domain=domain, seed=domain.rng)
    Matrix([
    [-1,          0,          0],
    [ 0,  sqrt(2)/2, -sqrt(2)/2],
    [ 0, -sqrt(2)/2, -sqrt(2)/2]])

    >>> domain.rng.seed(2)
    >>> orthogonal(3, domain=domain, seed=domain.rng)
    Matrix([
    [0, -1,  0],
    [0,  0, -1],
    [1,  0,  0]])

    >>> domain.rng.seed(1)
    >>> orthogonal(3, domain=domain, seed=domain.rng)
    Matrix([
    [-1,          0,          0],
    [ 0,  sqrt(2)/2, -sqrt(2)/2],
    [ 0, -sqrt(2)/2, -sqrt(2)/2]])


Random Matrix Generation Functions Reference
--------------------------------------------

.. automodule:: sympy.matrices.random
   :members:
   :member-order: bysource
