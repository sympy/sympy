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
from sympy.utilities.misc import debug

from sympy.physics.quantum.innerproduct import InnerProduct
from sympy.physics.quantum.kind import KetKind, BraKind, OperatorKind
from sympy.physics.quantum.operator import (
    OuterProduct, IdentityOperator, Operator
)
from sympy.physics.quantum.slidingtransform import DipatchingSlidingTransform
from sympy.physics.quantum.state import BraBase, KetBase, StateBase
from sympy.physics.quantum.tensorproduct import TensorProduct


#-----------------------------------------------------------------------------
# Multipledispatch based transformed for Mul and Pow
#-----------------------------------------------------------------------------

_postprocess_state_mul = DipatchingSlidingTransform()
"""Transform a pair of expression in a Mul to their canonical form.

All functions that are registered with this dispatcher need to take
two inputs and return either tuple of transformed outputs, or None if no
transform is applied. The output tuple is inserted into the right place
of the ``Mul`` that is being put into canonical form. It works something like
the following:

``Mul(a, b, c, d, e, f) -> Mul(*(_transform_state_pair(a, b) + (c, d, e, f))))``

The transforms here are always applied when quantum objects are multiplied.

THIS IS EXPERIMENTAL.

However, users of ``sympy.physics.quantum`` can import this dispatcher and
register their own transforms to control the canonical form of products
of quantum expressions.
"""


@_postprocess_state_mul.binary.register(Expr, Expr)
def _transform_expr(a, b):
    """Default transformer that does nothing for base types."""
    return None


# The identity times anything is the anything.
_postprocess_state_mul.binary.add(
    (IdentityOperator, Expr),
    lambda x, y: (y,),
    on_ambiguity=ambiguity_register_error_ignore_dup
)
_postprocess_state_mul.binary.add(
    (Expr, IdentityOperator),
    lambda x, y: (x,),
    on_ambiguity=ambiguity_register_error_ignore_dup
)
_postprocess_state_mul.binary.add(
    (IdentityOperator, IdentityOperator),
    lambda x, y: S.One,
    on_ambiguity=ambiguity_register_error_ignore_dup
)

@_postprocess_state_mul.binary.register(BraBase, KetBase)
def _transform_bra_ket(a, b):
    """Transform a bra*ket -> InnerProduct(bra, ket)."""
    return (InnerProduct(a, b),)

@_postprocess_state_mul.binary.register(KetBase, BraBase)
def _transform_ket_bra(a, b):
    """Transform a keT*bra -> OuterProduct(ket, bra)."""
    return (OuterProduct(a, b),)

@_postprocess_state_mul.binary.register(KetBase, KetBase)
def _transform_ket_ket(a, b):
    """Raise a TypeError if a user tries to multiply two kets.

    Multiplication based on `*` is not a shorthand for tensor products.
    """
    raise TypeError(
        'Multiplication of two kets is not allowed. Use TensorProduct instead.'
    )

@_postprocess_state_mul.binary.register(BraBase, BraBase)
def _transform_bra_bra(a, b):
    """Raise a TypeError if a user tries to multiply two bras.

    Multiplication based on `*` is not a shorthand for tensor products.
    """
    raise TypeError(
        'Multiplication of two bras is not allowed. Use TensorProduct instead.'
    )

@_postprocess_state_mul.binary.register(OuterProduct, KetBase)
def _transform_op_ket(a, b):
    return (InnerProduct(a.bra, b), a.ket)

@_postprocess_state_mul.binary.register(BraBase, OuterProduct)
def _transform_bra_op(a, b):
    return (InnerProduct(a, b.ket), b.bra)

@_postprocess_state_mul.binary.register(TensorProduct, KetBase)
def _transform_tp_ket(a, b):
    """Raise a TypeError if a user tries to multiply TensorProduct(*kets)*ket.

    Multiplication based on `*` is not a shorthand for tensor products.
    """
    if a.kind == KetKind:
        raise TypeError(
            'Multiplication of TensorProduct(*kets)*ket is invalid.'
        )

@_postprocess_state_mul.binary.register(KetBase, TensorProduct)
def _transform_ket_tp(a, b):
    """Raise a TypeError if a user tries to multiply ket*TensorProduct(*kets).

    Multiplication based on `*` is not a shorthand for tensor products.
    """
    if b.kind == KetKind:
        raise TypeError(
            'Multiplication of ket*TensorProduct(*kets) is invalid.'
        )

@_postprocess_state_mul.binary.register(TensorProduct, BraBase)
def _transform_tp_bra(a, b):
    """Raise a TypeError if a user tries to multiply TensorProduct(*bras)*bra.

    Multiplication based on `*` is not a shorthand for tensor products.
    """
    if a.kind == BraKind:
        raise TypeError(
            'Multiplication of TensorProduct(*bras)*bra is invalid.'
        )

@_postprocess_state_mul.binary.register(BraBase, TensorProduct)
def _transform_bra_tp(a, b):
    """Raise a TypeError if a user tries to multiply bra*TensorProduct(*bras).

    Multiplication based on `*` is not a shorthand for tensor products.
    """
    if b.kind == BraKind:
        raise TypeError(
            'Multiplication of bra*TensorProduct(*bras) is invalid.'
        )

@_postprocess_state_mul.binary.register(TensorProduct, TensorProduct)
def _transform_tp_tp(a, b):
    """Combine a product of tensor products if their number of args matches."""
    debug('_transform_tp_tp', a, b)
    if len(a.args) == len(b.args):
        if a.kind == BraKind and b.kind == KetKind:
            return tuple([InnerProduct(i, j) for (i, j) in zip(a.args, b.args)])
        else:
            return (TensorProduct(*(i*j for (i, j) in zip(a.args, b.args))), )

@_postprocess_state_mul.binary.register(OuterProduct, OuterProduct)
def _transform_op_op(a, b):
    """Extract an inner produt from a product of outer products."""
    return (InnerProduct(a.bra, b.ket), OuterProduct(a.ket, b.bra))


#-----------------------------------------------------------------------------
# Postprocessing transforms for Mul and Pow
#-----------------------------------------------------------------------------


def _postprocess_state_pow(expr):
    """Handle bras and kets raised to powers.

    Under ``*`` multiplication this is invalid. Users should use a
    TensorProduct instead.
    """
    base, exp = expr.as_base_exp()
    if base.kind == KetKind or base.kind == BraKind:
        raise TypeError(
            'A bra or ket to a power is invalid, use TensorProduct instead.'
        )


def _postprocess_tp_pow(expr):
    """Handle TensorProduct(*operators)**(positive integer).

    This handles a tensor product of operators, to an integer power.
    The power here is interpreted as regular multiplication, not
    tensor product exponentiation. The form of exponentiation performed
    here leaves the space and dimension of the object the same.

    This operation does not make sense for tensor product's of states.
    """
    base, exp = expr.as_base_exp()
    debug('_postprocess_tp_pow: ', base, exp, expr.args)
    if isinstance(base, TensorProduct) and exp.is_integer and exp.is_positive and base.kind == OperatorKind:
        new_args = [a**exp for a in base.args]
        return TensorProduct(*new_args)


#-----------------------------------------------------------------------------
# Register the transformers with Basic._constructor_postprocessor_mapping
#-----------------------------------------------------------------------------


Basic._constructor_postprocessor_mapping[StateBase] = {
    "Mul": [_postprocess_state_mul],
    "Pow": [_postprocess_state_pow]
}

Basic._constructor_postprocessor_mapping[TensorProduct] = {
    "Mul": [_postprocess_state_mul],
    "Pow": [_postprocess_tp_pow]
}

Basic._constructor_postprocessor_mapping[Operator] = {
    "Mul": [_postprocess_state_mul]
}
