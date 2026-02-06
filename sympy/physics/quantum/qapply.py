"""Quantum operator application.

This module provides the qapply() function which applies quantum operators
to quantum states, handling distribution over linear combinations and
complex quantum mechanical expressions symbolically.
"""

from __future__ import annotations
from typing import TYPE_CHECKING, Any,  TypedDict,   cast

from sympy.concrete import Sum
from sympy.core.add import Add
from sympy.core.basic import Basic
from sympy.core.expr import Expr
from sympy.core.kind import NumberKind
from sympy.core.mul import Mul
from sympy.core.power import Pow
from sympy.core.singleton import S
from sympy.core.sympify import sympify, _sympify
from sympy.functions.elementary.complexes import Abs
from sympy.integrals import Integral
from sympy.multipledispatch import Dispatcher
from sympy.multipledispatch.dispatcher import ambiguity_register_error_ignore_dup


from sympy.physics.quantum.anticommutator import AntiCommutator
from sympy.physics.quantum.commutator import Commutator
from sympy.physics.quantum.dagger import Dagger
from sympy.physics.quantum.density import Density
from sympy.physics.quantum.innerproduct import InnerProduct
from sympy.physics.quantum.operator import (
    OuterProduct, Operator
)
from sympy.physics.quantum.state import KetBase, BraBase, Wavefunction
from sympy.physics.quantum.tensorproduct import TensorProduct
from sympy.physics.quantum.slidingtransform import SlidingTransform


if TYPE_CHECKING:
    from typing import Unpack

__all__ = [
    'qapply',
    'QApplyOptions'
]


class QApplyOptions(TypedDict, total=False):
    """Type definition for qapply options dictionary."""
    ip_doit: bool
    sum_doit: bool
    dagger: bool


#-----------------------------------------------------------------------------
# Utilities
#-----------------------------------------------------------------------------


def _ip_doit_func(e):
    """Transform the inner products in an expression by calling ``.doit()``."""
    return e.replace(InnerProduct, lambda *args: InnerProduct(*args).doit())


def _sum_doit_func(e):
    """Transform the sums in an expression by calling ``.doit()``."""
    result = e.replace(Sum, lambda *args: Sum(*args).doit())
    result = result.replace(Integral, lambda *args: Integral(*args).doit())
    return result


def _handle_doit_unary(f):
    """Decorator that applies doit() to inner products and sums after transformation."""

    def apply_doit(expr, **options):
        sum_doit = options.get('sum_doit', False)
        ip_doit = options.get('ip_doit', True)
        result = f(expr, **options)
        result = _ip_doit_func(result) if ip_doit else result
        result = _sum_doit_func(result) if sum_doit else result
        return result

    return apply_doit


def _flatten_mul(f):
    """Decorator that flattens Mul expressions in transformation results."""

    def g(*args, **options):
        seq = f(*args, **options)
        if seq is None:
            return None
        result = []
        for item in seq:
            if isinstance(item, Mul):
                result.extend(item.args)
            else:
                result.append(item)
        return tuple(result)

    return g


#-----------------------------------------------------------------------------
# _qapply_unary
#-----------------------------------------------------------------------------

_qapply_unary = Dispatcher('_qapply_unary')


@_qapply_unary.register(Dagger)
def _qapply_unary_dagger(expr: Dagger, **options: Unpack[QApplyOptions]) -> Dagger:
    """Apply qapply to the argument of Dagger."""
    return Dagger(qapply(expr.args[0], **options))


@_qapply_unary.register(Density)
def _qapply_unary_density(expr: Density, **options: Unpack[QApplyOptions]) -> Density:
    """Apply qapply to states in a Density matrix."""
    # For a Density operator call qapply on its state
    new_args = [(qapply(state, **options), prob) for (state,  # type: ignore[misc,has-type]
                     prob) in expr.args]
    return Density(*new_args)


@_qapply_unary.register(TensorProduct)
def _qapply_unary_tp(expr: TensorProduct, **options: Unpack[QApplyOptions]) -> TensorProduct:
    """Apply qapply to each factor in a TensorProduct."""
    new_args = [qapply(t, **options) for t in expr.args]
    return TensorProduct(*new_args)


@_qapply_unary.register(Add)
def _qapply_unary_state(expr: Add, **options: Unpack[QApplyOptions]) -> Expr:
    """Apply qapply to each term in a sum."""
    result: Expr = S.Zero
    for arg in expr.args:
        result += qapply(arg, **options)
    return result


@_qapply_unary.register(Pow)
def _qapply_unary_pow(expr: Pow, **options: Unpack[QApplyOptions]) -> Expr:
    """Apply qapply to the base of a power expression."""
    # For a Pow, call qapply on its base.
    base, exp = expr.as_base_exp()
    return qapply(base, **options)**exp


@_qapply_unary.register(Mul)
@_handle_doit_unary  # type: ignore[arg-type,misc]
def _qapply_unary_mul(expr: Mul, **options: Unpack[QApplyOptions]) -> Expr:
    """We have a Mul where there might be actual operators to apply."""
    dagger = options.get('dagger', False)
    result = qapply_Mul(expr, **options)
    if dagger:
        result = Dagger(qapply_Mul(Dagger(result), **options))   # type: ignore[arg-type]
    return result


@_qapply_unary.register(Sum)
@_handle_doit_unary  # type: ignore[arg-type,misc]
def _qapply_unary_sum(expr: Sum, **options: Unpack[QApplyOptions]) -> Sum:
    """For a Sum, call qapply on its function."""
    result = Sum(qapply(expr.function, **options), *expr.limits)
    return result


@_qapply_unary.register(Integral)
@_handle_doit_unary  # type: ignore[arg-type,misc]
def _qapply_unary_integral(expr: Integral, **options: Unpack[QApplyOptions]) -> Integral:
    """For a Sum, call qapply on its function."""
    result = Integral(qapply(expr.function, **options), *expr.limits)
    return result


@_qapply_unary.register(Abs)
def _qapply_unary_abs(expr: Abs, **options: Unpack[QApplyOptions]) -> Abs:
    """Apply qapply to the argument of Abs."""
    return Abs(qapply(expr.args[0], **options))



@_qapply_unary.register(Expr)
def _qapply_unary_expr(expr: Expr, **options: Unpack[QApplyOptions]) -> None:
    """Default handler for Expr - no transformation."""
    return None


@_qapply_unary.register((int, float, complex))
def _qapply_unary_number(expr: int | float | complex, **options: Unpack[QApplyOptions]) -> Expr:
    """Convert numbers to SymPy expressions."""
    return sympify(expr)


#-----------------------------------------------------------------------------
# qapply
#-----------------------------------------------------------------------------

def qapply(e: Basic | int | float | complex, **options: Unpack[QApplyOptions]) -> Expr:
    """Apply quantum operators to quantum states in expressions.

    Explanation
    ===========

    This function applies quantum operators to quantum states, handling the
    symbolic computation of operator-state interactions. It automatically
    distributes operators over linear combinations (sums) of states, evaluates
    operator products, and processes complex quantum mechanical expressions
    while preserving quantum mechanical properties like linearity.

    Parameters
    ==========

    e : Expr
        A SymPy expression containing quantum operators and states.

    **options : QApplyOptions
        ip_doit : bool (default True)
            Automatically call ``.doit()`` on inner products to evaluate them.
        sum_doit : bool (default False)
            Automatically call ``.doit()`` on sums to expand them.
        dagger : bool (default False)
            Apply operators to bras from the left using dagger transformation.

    Returns
    =======

    result : Expr
        The expression with operators applied to states.

    Examples
    ========

    Basic operator application to quantum states:

        >>> from sympy.physics.quantum import qapply
        >>> from sympy.physics.quantum.spin import Jz, JzKet
        >>> qapply(Jz * JzKet(1, 1))
        hbar*|1,1>

    Bosonic creation and annihilation operators:

        >>> from sympy.physics.quantum.boson import BosonOp, BosonFockKet
        >>> a = BosonOp("a")
        >>> qapply(a * BosonFockKet(3))
        sqrt(3)*|2>

    Automatic distribution over sums of states:

        >>> from sympy.physics.quantum.state import Ket
        >>> from sympy.physics.quantum.operator import Operator
        >>> A = Operator('A')
        >>> ket1 = Ket('0')
        >>> ket2 = Ket('1')
        >>> qapply(A * (ket1 + ket2))
        A*|0> + A*|1>

    Operator powers are expanded and applied:

        >>> from sympy.physics.quantum.dagger import Dagger
        >>> a = BosonOp("a")
        >>> qapply(a**2 * BosonFockKet(3))
        sqrt(6)*|1>
        >>> qapply(Dagger(a)**2 * BosonFockKet(1))
        sqrt(6)*|3>

    Commutators are evaluated before application:

        >>> from sympy.physics.quantum.commutator import Commutator
        >>> from sympy.physics.quantum.spin import Jx, Jy, JzKet
        >>> qapply(Commutator(Jx, Jy) * JzKet(1, 1))
        hbar**2*I*|1,1>

    Tensor products of operators and states:

        >>> from sympy.physics.quantum.tensorproduct import TensorProduct
        >>> a = BosonOp("a")
        >>> b = BosonOp("b")
        >>> qapply(TensorProduct(a, b) * TensorProduct(BosonFockKet(1), BosonFockKet(2)))
        sqrt(2)*|0>|1>

    Inner products with the ip_doit option:

        >>> qapply(JzKet(1, 1).dual * Jz * JzKet(1, 1))
        hbar
        >>> qapply(JzKet(1, 1).dual * Jz * JzKet(1, 1), ip_doit=False)
        hbar*<1,1|1,1>

    The dagger option applies operators to bras:

        >>> from sympy.physics.quantum.qubit import QubitBra
        >>> from sympy.physics.quantum.gate import H
        >>> qapply(QubitBra(0) * Dagger(H(0)), dagger=True)
        sqrt(2)*<0|/2 + sqrt(2)*<1|/2

    Complex nested expressions with double sums:

        >>> from sympy import Sum, symbols
        >>> p, q = symbols('p q', integer=True, nonnegative=True)
        >>> a = BosonOp("a")
        >>> expr = Sum(p * a, (p, 0, 1)) * Sum(BosonFockKet(q), (q, 1, 2))
        >>> qapply(expr)
        Sum(p*sqrt(q)*|q - 1>, (p, 0, 1), (q, 1, 2))

    See Also
    ========

    sympy.physics.quantum.operator.Operator
    sympy.physics.quantum.state.Ket, sympy.physics.quantum.state.Bra
    sympy.physics.quantum.tensorproduct.TensorProduct

    References
    ==========

    .. [1] https://en.wikipedia.org/wiki/Quantum_state
    .. [2] https://en.wikipedia.org/wiki/Operator_(physics)

    """
    ip_doit = options.get('ip_doit', True)

    # After this, e_basic is guaranteed to be a SymPy Basic
    e_basic = _sympify(e)

    # Using the kind API here helps us to narrow what types of expressions
    # we call ``_ip_doit_func`` on. Here, we are covering the case where there
    # is an inner product in a scalar input.
    if e_basic.kind == NumberKind:
        return _ip_doit_func(e_basic) if ip_doit else e_basic

    result = _qapply_unary(e_basic, **options)

    return e_basic if (result is None) else result


#-----------------------------------------------------------------------------
# qapply_Mul
#-----------------------------------------------------------------------------

_qapply_mul_unary = Dispatcher('_qapply_mul_unary')
_qapply_mul_binary = Dispatcher('_qapply_mul_binary')

class _QApplyMulUnary:
    def __call__(self, expr: Any, **kwargs: Unpack[QApplyOptions]):
        pass
class _QApplyMulBinary:
    def __call__(self, lhs: Any, rhs: Any, **kwargs: Unpack[QApplyOptions]):
        pass

qapply_Mul = SlidingTransform(
    unary=cast(_QApplyMulUnary, _qapply_mul_unary),
    binary=cast(_QApplyMulBinary, _qapply_mul_binary),
    reverse=True,
    from_args=False
)


#-----------------------------------------------------------------------------
# qapply_Mul: SlidingTransform unary
#-----------------------------------------------------------------------------


@_qapply_mul_unary.register(OuterProduct)
def _qapply_mul_unary_op(expr: OuterProduct, **options: Unpack[QApplyOptions]) -> tuple[Expr, ...]:
    """Split OuterProduct into ket and bra factors."""
    return (expr.ket, expr.bra)


@_qapply_mul_unary.register((Commutator, AntiCommutator))
@_flatten_mul
def _qapply_mul_unary_comm(expr: Commutator | AntiCommutator, **options: Unpack[QApplyOptions]) -> tuple[Expr, ...]:
    """Evaluate commutators and anticommutators."""
    return (expr.doit(),)


@_qapply_mul_unary.register(TensorProduct)
@_flatten_mul
def _qapply_mul_unary_tp(expr: TensorProduct, **options: Unpack[QApplyOptions]) -> tuple[Expr, ...]:
    """Apply qapply to each factor in TensorProduct."""
    new_args = [qapply(t, **options) for t in expr.args]
    return (TensorProduct(*new_args),)


@_qapply_mul_unary.register(Density)
def _qapply_mul_unary_density(expr: Density, **options: Unpack[QApplyOptions]) -> tuple[Expr, ...]:
    """Apply qapply to states in Density matrix."""
    # For a Density operator call qapply on its state
    new_args = [(qapply(state, **options), prob) for (state,  # type: ignore[misc,has-type]
                     prob) in expr.args]
    return (Density(*new_args),)


@_qapply_mul_unary.register(Pow)
@_flatten_mul
def _qapply_mul_unary_pow(expr: Pow, **options: Unpack[QApplyOptions]) -> tuple[Expr, ...]:
    """Expand integer powers or apply qapply to power expressions."""
    base, exp = expr.as_base_exp()
    if exp.is_Integer and exp > 0:
        return tuple(expr.base for i in range(int(expr.exp)))
    else:
        return (qapply(expr, **options),)


@_qapply_mul_unary.register((int, float, complex))
def _qapply_mul_unary_number(expr: int | float | complex, **options: Unpack[QApplyOptions]) -> tuple[Expr, ...]:
    """Convert numbers to SymPy expressions in Mul context."""
    return (sympify(expr),)


@_qapply_mul_unary.register(Expr)
def qapply_mul_unary_expr(expr: Expr, **options: Unpack[QApplyOptions]) -> None:
    """Default handler for Expr in Mul context - no transformation."""
    return None


#-----------------------------------------------------------------------------
# qapply_Mul: SlidingTransform binary
#-----------------------------------------------------------------------------


@_qapply_mul_binary.register(Operator, (KetBase, TensorProduct, Wavefunction))
@_flatten_mul
def _qapply_mul_binary_op_wavefunction(lhs: Operator, rhs: KetBase | TensorProduct | Wavefunction, **options: Unpack[QApplyOptions]) -> tuple[Expr, ...] | None:
    """Apply operator to ket, tensor product, or wavefunction."""
    _apply = getattr(lhs, '_apply_operator', None)
    if _apply is not None:
        try:
            result = _apply(rhs, **options)
        except NotImplementedError:
            result = None
    else:
        result = None

    if result is None:
        _apply_right = getattr(rhs, '_apply_from_right_to', None)
        if _apply_right is not None:
            try:
                result = _apply_right(lhs, **options)
            except NotImplementedError:
                result = None

    return None if result is None else (result,)


@_qapply_mul_binary.register(BraBase, Operator)
def _qapply_mul_binary_bra_op(lhs: BraBase, rhs: Operator, **options: Unpack[QApplyOptions]) -> None:
    """Register this to remind us that we only do this when dagger=True."""
    return None


@_qapply_mul_binary.register(Sum, Sum)
def _qapply_mul_binary_sum_sum(lhs: Sum, rhs: Sum, **options: Unpack[QApplyOptions]) -> tuple[Expr, ...] | None:
    """Combine two sums into a single sum."""
    # Require that both are OperatorKind, BraKind, KetKind, or UndefinedKind
    if lhs.kind != NumberKind and rhs.kind != NumberKind:
        if set(lhs.variables).intersection(set(rhs.variables)):
            raise ValueError(
                'Duplicated dummy indices in separate sums in qapply.'
            )
        limits = lhs.limits + rhs.limits
        result = Sum(
            qapply(lhs.function*rhs.function, **options),
            *limits
        )
        return (result,)
    return None


@_qapply_mul_binary.register(Integral, Integral)
def _qapply_mul_binary_integral_integral(lhs: Integral, rhs: Integral, **options: Unpack[QApplyOptions]) -> tuple[Expr, ...] | None:
    """Combine two integrals into a single integral."""
    # Require that both are OperatorKind, BraKind, KetKind, or UndefinedKind
    if lhs.kind != NumberKind and rhs.kind != NumberKind:
        if set(lhs.variables).intersection(set(rhs.variables)):
            raise ValueError(
                'Duplicated dummy indices in separate sums in qapply.'
            )
        limits = lhs.limits + rhs.limits
        result = Integral(
            qapply(lhs.function*rhs.function, **options),
            *limits
        )
        return (result,)
    return None


@_qapply_mul_binary.register(Add, Add)
def _qapply_mul_binary_add_add(lhs: Add, rhs: Add, **options: Unpack[QApplyOptions]) -> tuple[Expr, ...] | None:
    """Distribute product of two sums using distributive property."""
    # Require that both are OperatorKind, BraKind, KetKind, or UndefinedKind
    if lhs.kind != NumberKind and rhs.kind != NumberKind:
        terms = []
        for lhs_term in lhs.args:
            for rhs_term in rhs.args:
                terms.append(qapply(lhs_term * rhs_term, **options))
        result = Add(*terms)
        return (result,)
    return None


@_qapply_mul_binary.register(Sum, Add)
def _qapply_mul_binary_sum_add(lhs: Sum, rhs: Add, **options: Unpack[QApplyOptions]) -> tuple[Expr, ...] | None:
    """Apply sum to add combination."""
    # Require that both are OperatorKind, BraKind, KetKind, or UndefinedKind
    if lhs.kind != NumberKind and rhs.kind != NumberKind:
        terms = (qapply(lhs.function*term, **options) for term in rhs.args)
        result = Sum(Add(*terms), *lhs.limits)
        return (result,)
    return None


@_qapply_mul_binary.register(Add, Sum)
def _qapply_mul_binary_add_sum(lhs: Add, rhs: Sum, **options: Unpack[QApplyOptions]) -> tuple[Expr, ...] | None:
    """Apply add to sum combination."""
    # Require that both are OperatorKind, BraKind, KetKind, or UndefinedKind
    if lhs.kind != NumberKind and rhs.kind != NumberKind:
        terms = (qapply(term*rhs.function, **options) for term in lhs.args)
        result = Sum(Add(*terms), *rhs.limits)
        return (result,)
    return None


@_qapply_mul_binary.register(Integral, Add)
def _qapply_mul_binary_integral_add(lhs: Integral, rhs: Add, **options: Unpack[QApplyOptions]) -> tuple[Expr, ...] | None:
    """Apply sum to add combination."""
    # Require that both are OperatorKind, BraKind, KetKind, or UndefinedKind
    if lhs.kind != NumberKind and rhs.kind != NumberKind:
        terms = (qapply(lhs.function*term, **options) for term in rhs.args)
        result = Integral(Add(*terms), *lhs.limits)
        return (result,)
    return None


@_qapply_mul_binary.register(Add, Integral)
def _qapply_mul_binary_add_integral(lhs: Add, rhs: Integral, **options: Unpack[QApplyOptions]) -> tuple[Expr, ...] | None:
    """Apply add to sum combination."""
    # Require that both are OperatorKind, BraKind, KetKind, or UndefinedKind
    if lhs.kind != NumberKind and rhs.kind != NumberKind:
        terms = (qapply(term*rhs.function, **options) for term in lhs.args)
        result = Integral(Add(*terms), *rhs.limits)
        return (result,)
    return None


@_qapply_mul_binary.register(Sum, Integral)
def _qapply_mul_binary_sum_integral(lhs: Sum, rhs: Integral, **options: Unpack[QApplyOptions]) -> tuple[Expr, ...] | None:
    """Apply sum to integral combination."""
    # Require that both are OperatorKind, BraKind, KetKind, or UndefinedKind
    if lhs.kind != NumberKind and rhs.kind != NumberKind:
        result = Sum(
            Integral(
                qapply(lhs.function*rhs.function, **options),
                *rhs.limits
            ),
            *lhs.limits
        )
        return (result,)
    return None


@_qapply_mul_binary.register(Integral, Sum)
def _qapply_mul_binary_integral_sum(lhs: Integral, rhs: Sum, **options: Unpack[QApplyOptions]) -> tuple[Expr, ...] | None:
    """Apply integral to sum combination."""
    # Require that both are OperatorKind, BraKind, KetKind, or UndefinedKind
    if lhs.kind != NumberKind and rhs.kind != NumberKind:
        result = Integral(
            Sum(
                qapply(lhs.function*rhs.function, **options),
                *rhs.limits
            ),
            *lhs.limits
        )
        return (result,)
    return None


@_qapply_mul_binary.register(Add, Expr, on_ambiguity=ambiguity_register_error_ignore_dup)
def _qapply_mul_binary_add_expr(lhs: Add, rhs: Expr, **options: Unpack[QApplyOptions]) -> tuple[Expr, ...] | None:
    """Distribute expression over sum terms."""
    # Require that both are OperatorKind, BraKind, KetKind, or UndefinedKind
    if lhs.kind != NumberKind and rhs.kind != NumberKind:
        result = Add(*(qapply(item*rhs, **options) for item in lhs.args))
        return (result,)
    return None


@_qapply_mul_binary.register(Expr, Add, on_ambiguity=ambiguity_register_error_ignore_dup)
def _qapply_mul_binary_expr_add(lhs: Expr, rhs: Add, **options: Unpack[QApplyOptions]) -> tuple[Expr, ...] | None:
    """Distribute expression over sum terms."""
    # Require that both are OperatorKind, BraKind, KetKind, or UndefinedKind
    if lhs.kind != NumberKind and rhs.kind != NumberKind:
        result = Add(*(qapply(lhs*item, **options) for item in rhs.args))
        return (result,)
    return None


@_qapply_mul_binary.register((Sum, Integral), Expr, on_ambiguity=ambiguity_register_error_ignore_dup)
def _qapply_mul_binary_sum_expr(lhs: Sum | Integral, rhs: Expr, **options: Unpack[QApplyOptions]) -> tuple[Expr, ...] | None:
    """Apply expression to sum/integral function."""
    # Require that both are OperatorKind, BraKind, KetKind, or UndefinedKind
    if lhs.kind != NumberKind and rhs.kind != NumberKind:
        result = lhs.func(qapply(lhs.function*rhs, **options), *lhs.limits)
        return (result,)
    return None


@_qapply_mul_binary.register(Expr, (Sum, Integral), on_ambiguity=ambiguity_register_error_ignore_dup)
def _qapply_mul_binary_expr_sum(lhs: Expr, rhs: Sum | Integral, **options: Unpack[QApplyOptions]) -> tuple[Expr, ...] | None:
    """Apply expression to sum/integral function."""
    # Require that both are OperatorKind, BraKind, KetKind, or UndefinedKind
    if lhs.kind != NumberKind and rhs.kind != NumberKind:
        result = rhs.func(qapply(lhs*rhs.function, **options), *rhs.limits)
        return (result,)
    return None


@_qapply_mul_binary.register(Expr, Expr, on_ambiguity=ambiguity_register_error_ignore_dup)
def _qapply_mul_binary_expr_expr(lhs: Expr, rhs: Expr, **options: Unpack[QApplyOptions]) -> None:
    """Default binary handler for Expr pairs - no transformation."""
    return None
