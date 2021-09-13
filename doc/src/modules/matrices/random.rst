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

require specific types of matrices.
A purely random generation of matrix entries does not grant any
specific shape or property of the matrix.
So, those matrices have to be created carefully using
classical classification theorems of linear algebra.

Generators on different types of matrices are presented below.

.. note::

   Even for carefully sampled matrices which are elements in topological matrix groups,
   the sample distribution does not claim to meet any uniform distribution
   (with respect to the groups haar measure). The sample procedures focuses
   mainly on educational and academic problems than a proper pre-defined
   distribution.

There are three different types of matrix generators:

   - the *base matrices* which from the atomic building block of the later
   - the *compound matrices* which are simply products of base matrices of given types
   - the *conjugate matrices* which are usually given as base or compound matrix $A$
     conjugate by another compound matrix $S$, i.e.

      .. math::

         B = S \cdot A \cdot S^{-1} \text{ or } C = S \cdot A \cdot S^{t}


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
    but can be controlled by **scalar_set** arguments.

    - :py:func:`permutation <sympy.matrices.random.permutation>` as product of
      :py:func:`transposition <sympy.matrices.random.transposition>` matrices

    - :py:func:`orthogonal <sympy.matrices.random.orthogonal>` as product of
      :py:func:`rotation <sympy.matrices.random.rotation>` and
      :py:func:`reflection <sympy.matrices.random.reflection>` matrices

    - :py:func:`invertible <sympy.matrices.random.invertible>` as product of
      :py:func:`elementary <sympy.matrices.random.elementary>` matrices

    - :py:func:`triangular <sympy.matrices.random.triangular>` as product of *triangular*
      :py:func:`elementary <sympy.matrices.random.elementary>` matrices

 3. *conjugate matrices* that arise as a product $SAS^{-1}$ of specific matrices $A$
    - usually given as a *normal form* - conjugate by an
    :py:func:`invertible <sympy.matrices.random.invertible>` $S$
    such as

    - :py:func:`diagonalizable <sympy.matrices.random.diagonalizable>` with
      :py:func:`diagonal_normal <sympy.matrices.random.diagonal_normal>`

    - :py:func:`idempotent <sympy.matrices.random.idempotent>` with
      :py:func:`projection <sympy.matrices.random.projection>` (with only zero and one entries)

    - :py:func:`trigonalizable <sympy.matrices.random.trigonalizable>` with
      :py:func:`jordan_normal <sympy.matrices.random.jordan_normal>`

    - :py:func:`nilpotent <sympy.matrices.random.nilpotent>` with
      :py:func:`jordan_normal <sympy.matrices.random.jordan_normal>` (with zero eigenvalues)

    *isometry conjugate matrices* that arise as a product $OAO^{-1} = OAO^{t}$ of specific matrices $A$
    where $A$ sets the (complex) eigenvalue spectrum.

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

    finally *symmetric matrices* given as a as product $SS^t$ resp. $S\bar{S}^t$ of $S$, such as

    - :py:func:`symmetric <sympy.matrices.random.symmetric>` $SS^t$ with
      :py:func:`invertible <sympy.matrices.random.invertible>` $S$

    - :py:func:`hermite <sympy.matrices.random.hermite>` $S\bar{S}^t$ with
      :py:func:`invertible <sympy.matrices.random.orthogonal>` $S$


In addition to the type of matrix, also the type of entries (as a commutative ring with one)
to be able to define values (**scalar_set**) can be specified,
from which the value entries (**scalar**) of the basic matrices are randomly generated.

Because the complexity and amount of entries in the generated compound matrices
in addition to the **scalar_set**, also by the number of base matrices multiplied for generation
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

gives an *invertible* matrix ``M``

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


Matrix Functions Reference
--------------------------

.. automodule:: sympy.matrices.random
   :members:
   :member-order: bysource
