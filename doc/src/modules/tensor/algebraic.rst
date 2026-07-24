.. _tensor-algebraic:

Algebraic Tensor Products
=========================

.. module:: sympy.tensor.algebraic

The algebraic tensor module provides support for tensor products of
matrix-like objects and their linear combinations. Unlike
:mod:`~sympy.tensor.array` which operates on concrete indexed components,
algebraic tensors work with symbolic matrix expressions and preserve
the tensor product structure.

This module was used to calculate all results of the mathematical
physics paper `arXiv:2511.08159 <https://arxiv.org/abs/2511.08159>`_
on finite noncommutative geometry.

A pure tensor is a tensor product of matrix (or vector) factors:

    >>> from sympy.matrices.expressions import MatrixSymbol
    >>> from sympy.tensor.algebraic import AlgebraicPureTensor
    >>> A = MatrixSymbol("A", 3, 4)
    >>> B = MatrixSymbol("B", 4, 5)
    >>> T = AlgebraicPureTensor(A, B)
    >>> T
    A ⊗ B

Pure tensors support addition, forming algebraic tensors (linear
combinations of pure tensors):

    >>> C = MatrixSymbol("C", 3, 4)
    >>> D = MatrixSymbol("D", 4, 5)
    >>> S = AlgebraicPureTensor(C, D)
    >>> S + T
    A ⊗ B + C ⊗ D

Composition of pure tensors performs factor-wise matrix multiplication:

    >>> from sympy.tensor.algebraic import compose_algebraic_pure_tensors
    >>> X = MatrixSymbol("X", 4, 3)
    >>> Y = MatrixSymbol("Y", 5, 4)
    >>> compose_algebraic_pure_tensors(T, AlgebraicPureTensor(X, Y))
    A*X ⊗ B*Y

The :func:`~sympy.tensor.algebraic.simplify.tensorsimplify` function
can simplify algebraic tensor expressions by combining like terms and
reducing zero factors.

Worked Example: Dirac Commutator in Finite Geometry
----------------------------------------------------

The following example demonstrates the simplification power of this
module. We define a Dirac operator :math:`D` and an algebra element
:math:`a` as algebraic tensors in a finite noncommutative geometry
model (from `arXiv:2511.08159 <https://arxiv.org/abs/2511.08159>`_).
Their commutator :math:`[D, a] = D \cdot a - a \cdot D` expands
to many pure tensor terms but simplifies to only 8.

First, define the fermion mass symbols (noncommutative) and the
complex scalar symbols for the algebra element:

    >>> from sympy import symbols, Matrix, eye
    >>> from sympy.tensor.algebraic import (
    ...     AlgebraicPureTensor, AlgebraicTensor, tensorsimplify)

    >>> (upsilon_R, upsilonc_nu, cupsilon_nu, upsilon_nu,
    ...  upsilont_nu, upsilonc_R) = symbols(
    ...     r"\Upsilon_R, \Upsilon^*_\nu, \overline{\Upsilon}_\nu, "
    ...     r"\Upsilon_\nu, \Upsilon^t_\nu, \Upsilon^*_R",
    ...     commutative=False)
    >>> (upsilonc_u, cupsilon_u, upsilon_u, upsilont_u) = symbols(
    ...     r"\Upsilon^*_u, \overline{\Upsilon}_u, \Upsilon_u, "
    ...     r"\Upsilon^t_u", commutative=False)
    >>> (upsilonc_e, cupsilon_e, upsilon_e, upsilont_e) = symbols(
    ...     r"\Upsilon^*_e, \overline{\Upsilon}_e, \Upsilon_e, "
    ...     r"\Upsilon^t_e", commutative=False)
    >>> (upsilonc_d, cupsilon_d, upsilon_d, upsilont_d) = symbols(
    ...     r"\Upsilon^*_d, \overline{\Upsilon}_d, \Upsilon_d, "
    ...     r"\Upsilon^t_d", commutative=False)

    >>> (z, w, alpha, beta, gamma, delta) = symbols(
    ...     r"z, w, \alpha, \beta, \gamma, \delta",
    ...     complex=True)
    >>> (m11, m12, m13, m21, m22, m23,
    ...  m31, m32, m33) = symbols(
    ...     r"m_{11}, m_{12}, m_{13}, m_{21}, m_{22}, "
    ...     r"m_{23}, m_{31}, m_{32}, m_{33}", complex=True)

Define the Dirac operator as a sum of four pure tensors, each with
three matrix factors (2x2, 4x4, 4x4):

    >>> proj_L = Matrix([[1, 0], [0, 0]])
    >>> proj_R = Matrix([[0, 0], [0, 1]])
    >>> id4 = eye(4)
    >>> pL_f = Matrix([[1,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0]])
    >>> pR_f = Matrix([[0,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,1]])
    >>> m_nu = Matrix([
    ...     [0, 0, upsilon_R, upsilonc_nu],
    ...     [0, 0, cupsilon_nu, 0],
    ...     [upsilonc_R, upsilont_nu, 0, 0],
    ...     [upsilon_nu, 0, 0, 0]])
    >>> m_u = Matrix([
    ...     [0, 0, 0, upsilonc_u],
    ...     [0, 0, cupsilon_u, 0],
    ...     [0, upsilont_u, 0, 0],
    ...     [upsilon_u, 0, 0, 0]])
    >>> m_e = Matrix([
    ...     [0, 0, 0, upsilonc_e],
    ...     [0, 0, cupsilon_e, 0],
    ...     [0, upsilont_e, 0, 0],
    ...     [upsilon_e, 0, 0, 0]])
    >>> m_d = Matrix([
    ...     [0, 0, 0, upsilonc_d],
    ...     [0, 0, cupsilon_d, 0],
    ...     [0, upsilont_d, 0, 0],
    ...     [upsilon_d, 0, 0, 0]])

    >>> D1 = AlgebraicPureTensor(proj_L, pL_f, m_nu)
    >>> D2 = AlgebraicPureTensor(proj_L, pR_f, m_u)
    >>> D3 = AlgebraicPureTensor(proj_R, pL_f, m_e)
    >>> D4 = AlgebraicPureTensor(proj_R, pR_f, m_d)
    >>> Dirac = D1 + D2 + D3 + D4

Define the algebra element :math:`a` as a sum of three pure tensors
(the representation :math:`\pi(a)`):

    >>> a_f1 = Matrix([[z, 0], [0, w]])
    >>> a_f2 = Matrix([[alpha, beta], [gamma, delta]])
    >>> a_m  = Matrix([
    ...     [z, 0, 0, 0],
    ...     [0, m11, m12, m13],
    ...     [0, m21, m22, m23],
    ...     [0, m31, m32, m33]])
    >>> a_f3 = Matrix([[0,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,0]])

    >>> a_pi1 = AlgebraicPureTensor(a_f1, id4, pL_f)
    >>> a_pi2 = AlgebraicPureTensor(a_f2, id4, pR_f)
    >>> a_pi3 = AlgebraicPureTensor(eye(2), a_m, a_f3)
    >>> a = a_pi1 + a_pi2 + a_pi3

Compute the commutator, expand, and simplify:

    >>> da = (Dirac * a - a * Dirac).expand()
    >>> len(da.args)
    18
    >>> da_simp = tensorsimplify(da)
    >>> len(da_simp.args)
    8

The simplification reduces the expression from 18 pure tensor terms to
just 8, by combining like terms and eliminating zero factors.

Classes
-------

.. autoclass:: AlgebraicPureTensor
   :members:

.. autoclass:: AlgebraicTensor
   :members:

.. autoclass:: AlgebraicZeroTensor
   :members:

.. autoclass:: ShapeMismatchError
   :members:

Functions
---------

.. autofunction:: algebraic_tensor_product

.. autofunction:: compose_algebraic_pure_tensors

.. autofunction:: compose_algebraic_tensors

.. autofunction:: tensorsimplify
