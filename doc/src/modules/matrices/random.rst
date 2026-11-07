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

A purely random generation of matrix entries, as
:py:func:`randMatrix <sympy.matrices.dense.randMatrix>`
does not guarantee any specific shape or property of the matrix.
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

Generating Random Matrices
--------------------------

There are different ways of matrix generation.

   - the *base matrices* which form the atomic building block of the later
   - the *compound matrices* which are simply products of base matrices of given types

 1. *base matrices* are simple matrix types with only a few non trivial entries.
    Those entries can set explicit by arguments or, if not given, are chosen randomly.

    - :py:func:`elementary <sympy.matrices.random.elementary>`

 2. *compound matrices* are build as a multiplication of several base matrices.
    Since by this the entries of the base matrices are chosen randomly
    but can be controlled by **scalars** arguments.

    - :py:func:`invertible <sympy.matrices.random.square>` as product of
      :py:func:`elementary <sympy.matrices.random.elementary>` matrices

    - :py:func:`triangular <sympy.matrices.random.triangular>` as product of *triangular*
      :py:func:`elementary <sympy.matrices.random.elementary>` matrices

In addition to the type of matrix, also the type of entries (as a commutative ring with one)
to be able to define values (**scalars**) can be specified,
from which the value entries (**scalar**) of the basic matrices are randomly generated.

Because the complexity and amount of entries in the generated compound matrices
in addition to the **scalars**, also by the number of base matrices multiplied for generation
This can be set using the argument **length**.

The product of *elementary* matrices

    >>> from sympy.matrices.random import elementary

    >>> A = elementary(3, index=(0, 2), scalar=-2)
    >>> A
    Matrix([
    [1, 0, -2],
    [0, 1,  0],
    [0, 0,  1]])

    >>> B = elementary(3, index=(1, 2), scalar=1)
    >>> B
    Matrix([
    [1, 0, 0],
    [0, 1, 1],
    [0, 0, 1]])

    >>> C = elementary(3, index=(1, 0), scalar=-1)
    >>> C
    Matrix([
    [ 1, 0, 0],
    [-1, 1, 0],
    [ 0, 0, 1]])

gives an invertible *square* matrix $\mathbf{M}$

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

By default the module uses the global sympy random number generator
function ``sample()`` which is a methode of an instance of ``random.Random()``.
The random state can be set by invoking ``seed()``.

    .. ..testsetup::

       >>> from sympy.core.random import rng
       >>> _rng_state = rng.getstate()

    >>> from sympy.core.random import seed
    >>> from sympy.matrices.random import square

    >>> seed(1)
    >>> square(3)
    Matrix([
    [-1, -1, -1],
    [-2, -1, -1],
    [ 0,  0,  1]])

    >>> seed(2)
    >>> square(3)
    Matrix([
    [ 1,  1, 0],
    [ 1,  2, 0],
    [-1, -1, 1]])

    >>> seed(1)
    >>> square(3)
    Matrix([
    [-1, -1, -1],
    [-2, -1, -1],
    [ 0,  0,  1]])

    .. ..testcleanup::

       >>> assert not rng.getstate() == _rng_state
       >>> rng.setstate(_rng_state)
       >>> assert rng.getstate() == _rng_state

Random Matrix Generation Functions Reference
============================================

Base Matrices
-------------

.. autofunction:: sympy.matrices.random.elementary

Compound Matrices
-----------------

.. autofunction:: sympy.matrices.random.triangular
.. autofunction:: sympy.matrices.random.square


Useful Functions
----------------

.. autofunction:: sympy.matrices.random.regular_to_singular
