"""Transforms that are always applied to quantum expressions.

This module uses the kind and _constructor_postprocessor_mapping APIs
to transform different combinations of Operators, Bras, and Kets into
Inner/Outer/TensorProducts. These transformations are registered
with the postprocessing API of core classes like `Mul` and `Pow` and
are always applied to any expression involving Bras, Kets, and
Operators. This API replaces the custom `__mul__` and `__pow__`
methods of the quantum classes, which were found to be inconsistent.

THIS IS EXPERIMENTAL.
"""
from sympy.core.basic import Basic
from sympy.core.expr import Expr
from sympy.core.mul import Mul
from sympy.core.singleton import S
from sympy.multipledispatch.dispatcher import (
    Dispatcher, ambiguity_register_error_ignore_dup
)

from sympy.physics.quantum.innerproduct import InnerProduct
from sympy.physics.quantum.kind import KetKind, BraKind
from sympy.physics.quantum.operator import (
    OuterProduct, IdentityOperator, Operator
)
from sympy.physics.quantum.state import BraBase, KetBase, StateBase
from sympy.physics.quantum.tensorproduct import TensorProduct


#-----------------------------------------------------------------------------
# Multipledispatch based transformed for Mul and Pow
#-----------------------------------------------------------------------------


_transform_state_pair = Dispatcher('_transform_state_pair')

@_transform_state_pair.register(Expr, Expr)
def _transform_expr(a, b):
    return None

@_transform_state_pair.register(IdentityOperator, Expr)
def _transform_id_any(a, b):
    return (b,)

_transform_state_pair.add(
    (IdentityOperator, Expr),
    lambda x, y: (y,),
    on_ambiguity=ambiguity_register_error_ignore_dup
)
_transform_state_pair.add(
    (Expr, IdentityOperator),
    lambda x, y: (x,),
    on_ambiguity=ambiguity_register_error_ignore_dup
)
_transform_state_pair.add(
    (IdentityOperator, IdentityOperator),
    lambda x, y: S.One,
    on_ambiguity=ambiguity_register_error_ignore_dup
)

@_transform_state_pair.register(BraBase, KetBase)
def _transform_bra_ket(a, b):
    return (InnerProduct(a, b),)

@_transform_state_pair.register(KetBase, BraBase)
def _transform_ket_bra(a, b):
    return (OuterProduct(a, b),)

@_transform_state_pair.register(KetBase, KetBase)
def _transform_ket_ket(a, b):
    return (TensorProduct(a, b),)

@_transform_state_pair.register(BraBase, BraBase)
def _transform_bra_bra(a, b):
    return (TensorProduct(a, b),)

@_transform_state_pair.register(OuterProduct, KetBase)
def _transform_op_ket(a, b):
    return (a.ket, InnerProduct(a.bra, b))

@_transform_state_pair.register(BraBase, OuterProduct)
def _transform_bra_op(a, b):
    return (InnerProduct(a, b.ket), b.bra)

@_transform_state_pair.register(TensorProduct, KetBase)
def _transform_tp_ket(a, b):
    if a.kind == KetKind:
        return (TensorProduct(*(a.args + (b,))),)

@_transform_state_pair.register(KetBase, TensorProduct)
def _transform_ket_tp(a, b):
    if b.kind == KetKind:
        return (TensorProduct(*((a,) + b.args)),)

@_transform_state_pair.register(TensorProduct, BraBase)
def _transform_tp_bra(a, b):
    if a.kind == BraKind:
        return (TensorProduct(*(a.args + (b,))),)

@_transform_state_pair.register(BraBase, TensorProduct)
def _transform_bra_tp(a, b):
    if b.kind == BraKind:
        return (TensorProduct(*((a,) + b.args)),)

@_transform_state_pair.register(TensorProduct, TensorProduct)
def _transform_tp_tp(a, b):
    if a.kind == BraKind and b.kind == KetKind:
        if len(a.args) == len(b.args):
            return tuple([InnerProduct(i, j) for (i, j) in zip(a.args, b.args)])

@_transform_state_pair.register(OuterProduct, OuterProduct)
def _transform_op_op(a, b):
    return (a.ket, InnerProduct(a.bra, b.ket), b.bra)


#-----------------------------------------------------------------------------
# Postprocessing transforms for Mul and Pow
#-----------------------------------------------------------------------------


def _postprocess_state_mul(expr):
    args = list(expr.args)
    result = []

    # Continue as long as we have at least 2 elements
    while len(args) > 1:
        # Get first two elements
        first = args.pop(0)
        second = args[0]  # Look at second element without popping yet

        transformed = _transform_state_pair(first, second)

        if transformed is None:
            # If transform returns None, append first element
            result.append(first)
        else:
            args.pop(0) # This item was transformed, pop and discard
            args.insert(0, transformed[-1])
            result.extend(transformed[:-1])

    # Append any remaining element
    if args:
        result.append(args[0])

    return Mul._from_args(result, is_commutative=False)


def _postprocess_state_pow(expr):
    base, exp = expr.as_base_exp()
    if exp.is_integer and exp.is_positive:
        return TensorProduct(*(base for i in range(int(exp))))


#-----------------------------------------------------------------------------
# Register the transformers with Basic._constructor_postprocessor_mapping
#-----------------------------------------------------------------------------


Basic._constructor_postprocessor_mapping[StateBase] = {
    "Mul": [_postprocess_state_mul],
    "Pow": [_postprocess_state_pow]
}

Basic._constructor_postprocessor_mapping[TensorProduct] = {
    "Mul": [_postprocess_state_mul]
}

Basic._constructor_postprocessor_mapping[Operator] = {
    "Mul": [_postprocess_state_mul]
}
