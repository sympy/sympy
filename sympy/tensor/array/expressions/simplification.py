"""
Opt-in simplification functions for array expressions.

This module provides two simplification functions for sums of tensor
products:

- :func:`collect_tensor_products` collects sums of tensor products whose
  factors coincide in all slots except at most one by using the
  multilinearity of the tensor product:

      A x B + A x C  ->  A x (B + C)

- :func:`commutativity_simplify` additionally decomposes explicit
  commutative factors into the elementary-array basis, so that terms
  which distribute their scalar entries differently across the factors
  can be compared, combined and cancelled.

Unlike the automatic collection of proportional terms performed by
``ArrayAdd`` canonicalization, these transformations can grow
intermediate expressions, so they are not applied automatically.
"""
from __future__ import annotations

import itertools
from typing import TYPE_CHECKING

from sympy.core.basic import Basic
from sympy.core.function import expand
from sympy.core.singleton import S
from sympy.matrices.matrixbase import MatrixBase
from sympy.matrices.expressions.matexpr import MatrixExpr
from sympy.tensor.array.dense_ndim_array import ImmutableDenseNDimArray
from sympy.tensor.array.ndim_array import NDimArray
from sympy.tensor.array.expressions.array_expressions import (
    ArrayAdd, ZeroArray, _array_add, _array_term_as_coeff_arrays,
    _array_term_from_coeff_arrays, get_shape)

if TYPE_CHECKING:
    from sympy.core.expr import Expr


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
    if isinstance(a, (NDimArray, MatrixBase)) and \
            isinstance(b, (NDimArray, MatrixBase)) and \
            c1.is_commutative is not False and \
            c2.is_commutative is not False:
        # Two explicit factors: perform the sum explicitly.
        aa = ImmutableDenseNDimArray(a) if isinstance(a, MatrixBase) else a
        bb = ImmutableDenseNDimArray(b) if isinstance(b, MatrixBase) else b
        slot = aa*c1 + bb*c2
    else:
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


def _is_explicit_commutative(factor):
    """Return True if *factor* is an explicit array whose entries contain
    no noncommutative symbols."""
    if isinstance(factor, (NDimArray, MatrixBase)):
        return all(getattr(s, "is_commutative", True) is not False
                   for s in factor.free_symbols)
    return False


def _iter_nonzero_entries(factor):
    """Iterate over the (index, entry) pairs of the nonzero entries of an
    explicit array or matrix."""
    shape = get_shape(factor)
    for index in itertools.product(*[range(int(d)) for d in shape]):
        entry = factor[index if len(index) > 1 else index[0]]
        if entry == S.Zero or getattr(entry, "is_zero", None) is True:
            continue
        yield index, entry


def _basis_array(shape, index):
    """Return the elementary array of the given *shape* with a single 1
    at position *index*."""
    total = 1
    for d in shape:
        total *= int(d)
    flat_position = 0
    for d, i in zip(shape, index):
        flat_position = flat_position*int(d) + i
    data = [S.Zero]*total
    data[flat_position] = S.One
    return ImmutableDenseNDimArray(data).reshape(*[int(d) for d in shape])


def commutativity_simplify(expr):
    """
    Simplify a sum of tensor products by decomposing explicit commutative
    factors into the elementary-array basis.

    Each addend of an ``ArrayAdd`` (whose addends must be tensor products
    with the same number of factors and the same factor shapes) is
    processed as follows: every factor which is an explicit array with no
    noncommutative entries is decomposed into elementary arrays, with its
    entries collected into the scalar prefactor.  Terms with the same
    elementary basis and the same remaining (noncommutative or symbolic)
    factors are then combined by summing their prefactors, cancelling
    terms whose sum of prefactors is zero.  Finally the surviving terms
    are recombined with :func:`collect_tensor_products`.

    This makes it possible to combine and cancel terms which distribute
    scalar entries differently across their factors, a situation
    :func:`collect_tensor_products` alone cannot handle since the factors
    differ in more than one slot.  The transformation is applied
    recursively to the whole expression tree.

    Examples
    ========

    The two addends below represent the same tensor, with the scalar
    ``z`` placed in a different factor:

    >>> from sympy import Array, Symbol
    >>> from sympy.tensor.array.expressions import ArrayAdd, ArrayTensorProduct
    >>> from sympy.tensor.array.expressions import commutativity_simplify
    >>> z = Symbol("z")
    >>> u = Symbol("u", commutative=False)
    >>> X = Array([[u, 0], [0, u]])
    >>> E00 = Array([[1, 0], [0, 0]])
    >>> zE00 = Array([[z, 0], [0, 0]])
    >>> expr = ArrayAdd(ArrayTensorProduct(zE00, E00, X), ArrayTensorProduct(E00, zE00, X))
    >>> commutativity_simplify(expr)
    ArrayTensorProduct(2*z, [[1, 0], [0, 0]], [[1, 0], [0, 0]], [[u, 0], [0, u]])

    Cancelling terms are detected:

    >>> from sympy.tensor.array.expressions import ArrayTensorProduct
    >>> expr = ArrayAdd(ArrayTensorProduct(zE00, E00, X), ArrayTensorProduct(-1, E00, zE00, X))
    >>> commutativity_simplify(expr)
    ZeroArray(2, 2, 2, 2, 2, 2)
    """
    if not isinstance(expr, Basic):
        return expr
    if expr.args and not isinstance(expr, (NDimArray, MatrixBase)):
        newargs = [commutativity_simplify(arg) for arg in expr.args]
        if list(newargs) != list(expr.args):
            expr = expr.func(*newargs)
    if isinstance(expr, ArrayAdd):
        return _commutativity_simplify_add(expr)
    return expr


def _commutativity_simplify_add(add: ArrayAdd):
    terms = [_array_term_as_coeff_arrays(arg) for arg in add.args]

    # The addends must be tensor products with aligned factor shapes:
    if any(not arrays for _, arrays in terms):
        return _collect_array_add(add)
    n = len(terms[0][1])
    if any(len(arrays) != n for _, arrays in terms):
        return _collect_array_add(add)
    shapes = [get_shape(a) for a in terms[0][1]]
    if any([get_shape(a) for a in arrays] != shapes for _, arrays in terms):
        return _collect_array_add(add)
    if any(None in (s or (None,)) or not all(getattr(d, "is_Integer", isinstance(d, int)) for d in s)
           for s in shapes):
        # Symbolic dimensions cannot be decomposed entry by entry:
        return _collect_array_add(add)

    comm_slots = [i for i in range(n)
                  if all(_is_explicit_commutative(arrays[i]) for _, arrays in terms)]
    if not comm_slots:
        return _collect_array_add(add)
    nc_slots = [i for i in range(n) if i not in comm_slots]

    # Decompose every term into elementary-basis components and group the
    # prefactors by (elementary basis, remaining factors):
    groups: dict[tuple, dict[tuple, Expr]] = {}
    for coeff, arrays in terms:
        decomposed: list[tuple[Expr, dict[int, tuple]]] = [(coeff, {})]
        for i in comm_slots:
            new_decomposed: list[tuple[Expr, dict[int, tuple]]] = []
            for prefactor, basis_map in decomposed:
                for index, entry in _iter_nonzero_entries(arrays[i]):
                    new_basis_map = dict(basis_map)
                    new_basis_map[i] = index
                    new_decomposed.append((prefactor*entry, new_basis_map))
            decomposed = new_decomposed
        nc_tuple = tuple(arrays[i] for i in nc_slots)
        for prefactor, basis_map in decomposed:
            key = tuple(basis_map[i] for i in comm_slots)
            ncmap = groups.setdefault(key, {})
            ncmap[nc_tuple] = ncmap.get(nc_tuple, S.Zero) + prefactor

    # Rebuild the surviving terms:
    new_terms = []
    for key, ncmap in groups.items():
        basis_arrays = {i: _basis_array(shapes[i], index)
                        for i, index in zip(comm_slots, key)}
        for nc_tuple, coeff in ncmap.items():
            coeff = expand(coeff)
            if coeff.is_zero is True:
                continue
            nc_iter = iter(nc_tuple)
            arrays = tuple(basis_arrays[i] if i in basis_arrays else next(nc_iter)
                           for i in range(n))
            new_terms.append(_array_term_from_coeff_arrays(coeff, arrays))

    if not new_terms:
        return ZeroArray(*[d for s in shapes for d in s])
    if len(new_terms) == 1:
        result = new_terms[0]
    else:
        result = _array_add(*new_terms)
    return collect_tensor_products(result)
