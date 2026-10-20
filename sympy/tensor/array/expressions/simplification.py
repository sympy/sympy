"""
Opt-in simplification functions for array expressions.

This module provides :func:`collect_tensor_products`, which collects sums
of tensor products whose factors coincide in all slots except at most one
by using the multilinearity of the tensor product:

    A x B + A x C  ->  A x (B + C)

Unlike the automatic collection of proportional terms performed by
``ArrayAdd`` canonicalization, this transformation can grow intermediate
factors, so it is not applied automatically.
"""
from __future__ import annotations

from sympy.core.basic import Basic
from sympy.core.singleton import S
from sympy.matrices.expressions.matexpr import MatrixExpr
from sympy.tensor.array.expressions.array_expressions import (
    ArrayAdd, _array_add, _array_term_as_coeff_arrays,
    _array_term_from_coeff_arrays, get_shape)


def collect_tensor_products(expr):
    """
    Collect sums of tensor products differing in at most one factor.

    Addends of an ``ArrayAdd`` whose tensor-product factors coincide in
    all slots except at most one are merged by linearity in the differing
    slot: ``A x B + A x C`` becomes ``A x (B + C)``.  Scalar coefficients
    are absorbed into the differing slot.  The transformation is applied
    recursively to the whole expression tree.

    Examples
    ========

    >>> from sympy import MatrixSymbol
    >>> from sympy.tensor.array.expressions import ArrayAdd, ArrayTensorProduct
    >>> from sympy.tensor.array.expressions import collect_tensor_products
    >>> A = MatrixSymbol("A", 2, 2)
    >>> B = MatrixSymbol("B", 2, 2)
    >>> C = MatrixSymbol("C", 2, 2)

    >>> collect_tensor_products(ArrayAdd(ArrayTensorProduct(A, B), ArrayTensorProduct(A, C)))
    ArrayTensorProduct(A, B + C)

    Scalar coefficients are absorbed into the merged slot:

    >>> expr = ArrayAdd(ArrayTensorProduct(2, A, B), ArrayTensorProduct(3, C, B))
    >>> collect_tensor_products(expr)
    ArrayTensorProduct(2*A + 3*C, B)

    Terms differing in more than one slot are left alone:

    >>> D = MatrixSymbol("D", 2, 2)
    >>> expr = ArrayAdd(ArrayTensorProduct(A, B), ArrayTensorProduct(C, D))
    >>> collect_tensor_products(expr) == expr
    True
    """
    if not isinstance(expr, Basic):
        return expr
    if expr.args:
        newargs = [collect_tensor_products(arg) for arg in expr.args]
        if list(newargs) != list(expr.args):
            expr = expr.func(*newargs)
    if isinstance(expr, ArrayAdd):
        return _collect_array_add(expr)
    return expr


def _scale_slot(coeff, arg):
    """Multiply the array factor *arg* by the scalar *coeff*."""
    if coeff is S.One:
        return arg
    return _array_term_from_coeff_arrays(coeff, (arg,))


def _merge_single_slot(term1, term2):
    """Merge two ``(coefficient, factors)`` pairs into a single pair
    representing their sum, or return None if they cannot be merged.

    The sum of two tensor products is itself a tensor product only when
    the factors coincide in all slots except at most one:
    ``c1*(A x B) + c2*(C x B) = (c1*A + c2*C) x B``.
    """
    (c1, arrays1), (c2, arrays2) = term1, term2
    if len(arrays1) != len(arrays2) or len(arrays1) == 0:
        return None
    if any(get_shape(a) != get_shape(b) for a, b in zip(arrays1, arrays2)):
        return None
    diff = [i for i, (a, b) in enumerate(zip(arrays1, arrays2)) if a != b]
    if len(diff) == 0:
        return (c1 + c2, arrays1)
    if len(diff) != 1:
        return None
    i = diff[0]
    a, b = arrays1[i], arrays2[i]
    sa, sb = _scale_slot(c1, a), _scale_slot(c2, b)
    if isinstance(sa, MatrixExpr) and isinstance(sb, MatrixExpr):
        slot = sa + sb
    else:
        slot = _array_add(sa, sb)
    new_arrays = list(arrays1)
    new_arrays[i] = slot
    return (S.One, tuple(new_arrays))


def _collect_array_add(add: ArrayAdd):
    terms = [_array_term_as_coeff_arrays(arg) for arg in add.args]

    # Merge greedily until no pair of terms can be merged any more (a
    # merged slot can enable further merges with other terms):
    while True:
        merged: list[tuple] = []
        for term in terms:
            for i, other in enumerate(merged):
                combined = _merge_single_slot(other, term)
                if combined is not None:
                    merged[i] = combined
                    break
            else:
                merged.append(term)
        if len(merged) == len(terms):
            break
        terms = merged

    if len(terms) == len(add.args):
        return add
    return _array_add(*[_array_term_from_coeff_arrays(c, arrays)
                        for c, arrays in terms])
