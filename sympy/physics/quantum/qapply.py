"""
Quantum operator application system for symbolic quantum mechanics.

The qapply function provides a comprehensive system for applying quantum operators
to quantum states, handling the complex interactions between operators, states,
and mathematical expressions in symbolic quantum mechanics.

The system operates through a hierarchical dispatch mechanism:

1. **Unary Transformations**: Individual expressions (Add, Mul, Pow, etc.) are
   processed through type-specific handlers that recursively apply qapply to
   sub-expressions while preserving mathematical structure.

2. **Multiplication Processing**: For Mul expressions containing operators and
   states, a specialized SlidingTransform processes the multiplication from
   right-to-left, applying binary transformation rules to adjacent pairs of
   factors.

3. **Operator-State Application**: When an operator is adjacent to a compatible
   state (ket, bra, wavefunction), the system attempts to apply the operator
   by calling the operator's `_apply_operator` method or the state's
   `_apply_from_right_to` method.

4. **Expression Distribution**: The system handles distribution of operators
   over sums and integrals, automatically propagating applications through
   linear combinations of states.

The qapply function preserves quantum mechanical properties like linearity
and handles complex expressions involving tensor products, density matrices,
commutators, and other quantum constructs while maintaining symbolic exactness.

This system is also extensible so developers can add handlers for new operators
and states through the multiple dispatch system.
"""

from sympy.concrete import Sum
from sympy.core.add import Add
from sympy.core.expr import Expr
from sympy.core.kind import NumberKind
from sympy.core.mul import Mul
from sympy.core.power import Pow
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


__all__ = [
    'qapply'
]


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
            return
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


@_qapply_unary.register(Expr)
def _qapply_unary_expr(expr, **options):
    """Default handler for Expr - no transformation."""
    return None


@_qapply_unary.register((int, float, complex))
def _qapply_unary_number(expr, **options):
    """Convert numbers to SymPy expressions."""
    return sympify(expr)


@_qapply_unary.register(Abs)
def _qapply_unary_abs(expr, **options):
    """Apply qapply to the argument of Abs."""
    return Abs(qapply(expr.args[0], **options))


@_qapply_unary.register(Add)
def _qapply_unary_state(expr, **options):
    """Apply qapply to each term in a sum."""
    result = 0
    for arg in expr.args:
        result += qapply(arg, **options)
    return result


@_qapply_unary.register(Pow)
def _qapply_unary_pow(expr, **options):
    """Apply qapply to the base of a power expression."""
    # For a Pow, call qapply on its base.
    base, exp = expr.as_base_exp()
    return qapply(base, **options)**exp


@_qapply_unary.register(Mul)
@_handle_doit_unary
def _qapply_unary_mul(expr, **options):
    """We have a Mul where there might be actual operators to apply."""
    dagger = options.get('dagger', False)
    result = qapply_Mul(expr, **options)
    if dagger:
        result = Dagger(qapply_Mul(Dagger(result), **options))
    return result


@_qapply_unary.register(Sum)
@_handle_doit_unary
def _qapply_unary_sum(expr, **options):
    """For a Sum, call qapply on its function."""
    result = Sum(qapply(expr.function, **options), *expr.limits)
    return result


@_qapply_unary.register(Integral)
@_handle_doit_unary
def _qapply_unary_integral(expr, **options):
    """For a Sum, call qapply on its function."""
    result = Integral(qapply(expr.function, **options), *expr.limits)
    return result


@_qapply_unary.register(Dagger)
def _qapply_unary_dagger(expr, **options):
    """Apply qapply to the argument of Dagger."""
    return Dagger(qapply(expr, **options))


@_qapply_unary.register(Density)
def _qapply_unary_density(expr, **options):
    """Apply qapply to states in a Density matrix."""
    # For a Density operator call qapply on its state
    new_args = [(qapply(state, **options), prob) for (state,
                     prob) in expr.args]
    return Density(*new_args)


@_qapply_unary.register(TensorProduct)
def _qapply_unary_tp(expr, **options):
    """Apply qapply to each factor in a TensorProduct."""
    new_args = [qapply(t, **options) for t in expr.args]
    return TensorProduct(*new_args)


#-----------------------------------------------------------------------------
# qapply
#-----------------------------------------------------------------------------

def qapply(e, **options):
    """Apply quantum operators to states and expressions."""
    ip_doit = options.get('ip_doit', True)

    e = _sympify(e)

    # Using the kind API here helps us to narrow what types of expressions
    # we call ``_ip_doit_func`` on. Here, we are covering the case where there
    # is an inner product in a scalar input.
    if e.kind == NumberKind:
        return _ip_doit_func(e) if ip_doit else e

    result = _qapply_unary(e, **options)

    return e if (result is None) else result


#-----------------------------------------------------------------------------
# qapply_Mul
#-----------------------------------------------------------------------------


qapply_Mul = SlidingTransform(
    unary=Dispatcher('_qapply_mul_unary'),
    binary=Dispatcher('_qapply_mul_binary'),
    reverse=True,
    from_args=False
)

#-----------------------------------------------------------------------------
# qapply_Mul: Unary
#-----------------------------------------------------------------------------


@qapply_Mul.unary.register((int, float, complex))
def _qapply_mul_unary_number(expr, **options):
    """Convert numbers to SymPy expressions in Mul context."""
    return (sympify(expr),)


@qapply_Mul.unary.register(Expr)
def qapply_mul_unary_expr(expr, **options):
    """Default handler for Expr in Mul context - no transformation."""
    return None


@qapply_Mul.unary.register(Pow)
@_flatten_mul
def _qapply_mul_unary_pow(expr, **options):
    """Expand integer powers or apply qapply to power expressions."""
    base, exp = expr.as_base_exp()
    if exp.is_Integer and exp > 0:
        return tuple(expr.base for i in range(int(expr.exp)))
    else:
        return (qapply(expr, **options),)


@qapply_Mul.unary.register(OuterProduct)
def _qapply_mul_unary_op(expr, **options):
    """Split OuterProduct into ket and bra factors."""
    return (expr.ket, expr.bra)


@qapply_Mul.unary.register((Commutator, AntiCommutator))
@_flatten_mul
def _qapply_mul_unary_comm(expr, **options):
    """Evaluate commutators and anticommutators."""
    return (expr.doit(),)


@qapply_Mul.unary.register(TensorProduct)
@_flatten_mul
def _qapply_mul_unary_tp(expr, **options):
    """Apply qapply to each factor in TensorProduct."""
    new_args = [qapply(t, **options) for t in expr.args]
    return (TensorProduct(*new_args),)


@qapply_Mul.unary.register(Density)
def _qapply_mul_unary_density(expr, **options):
    """Apply qapply to states in Density matrix."""
    # For a Density operator call qapply on its state
    new_args = [(qapply(state, **options), prob) for (state,
                     prob) in expr.args]
    return (Density(*new_args),)


#-----------------------------------------------------------------------------
# qapply_Mul: Binary
#-----------------------------------------------------------------------------


@qapply_Mul.binary.register(Sum, Sum)
def _qapply_mul_binary_sum_sum(lhs, rhs, **options):
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


@qapply_Mul.binary.register(Integral, Integral)
def _qapply_mul_binary_integral_integral(lhs, rhs, **options):
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


@qapply_Mul.binary.register(Add, Add)
def _qapply_mul_binary_add_add(lhs, rhs, **options):
    """Distribute product of two sums using distributive property."""
    # Require that both are OperatorKind, BraKind, KetKind, or UndefinedKind
    if lhs.kind != NumberKind and rhs.kind != NumberKind:
        terms = []
        for lhs_term in lhs.args:
            for rhs_term in rhs.args:
                terms.append(qapply(lhs_term * rhs_term, **options))
        result = Add(*terms)
        return (result,)


@qapply_Mul.binary.register(Operator, (KetBase, TensorProduct, Wavefunction))
@_flatten_mul
def _qapply_mul_binary_op_wavefunction(lhs, rhs, **options):
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


@qapply_Mul.binary.register(BraBase, Operator)
def _qapply_mul_binary_bra_op(lhs, rhs, **options):
    """Register this to remind us that we only do this when dagger=True."""
    return None


@qapply_Mul.binary.register(Sum, Add)
def _qapply_mul_binary_sum_add(lhs, rhs, **options):
    """Apply sum to add combination."""
    # Require that both are OperatorKind, BraKind, KetKind, or UndefinedKind
    if lhs.kind != NumberKind and rhs.kind != NumberKind:
        terms = (qapply(lhs.function*term, **options) for term in rhs.args)
        return Sum(Add(*terms), *lhs.limits)


@qapply_Mul.binary.register(Add, Sum)
def _qapply_mul_binary_add_sum(lhs, rhs, **options):
    """Apply add to sum combination."""
    # Require that both are OperatorKind, BraKind, KetKind, or UndefinedKind
    if lhs.kind != NumberKind and rhs.kind != NumberKind:
        terms = (qapply(term*rhs.function, **options) for term in lhs.args)
        return Sum(Add(*terms), *rhs.limits)


@qapply_Mul.binary.register(Integral, Add)
def _qapply_mul_binary_integral_add(lhs, rhs, **options):
    """Apply sum to add combination."""
    # Require that both are OperatorKind, BraKind, KetKind, or UndefinedKind
    if lhs.kind != NumberKind and rhs.kind != NumberKind:
        terms = (qapply(lhs.function*term, **options) for term in rhs.args)
        return Integral(Add(*terms), *lhs.limits)


@qapply_Mul.binary.register(Add, Integral)
def _qapply_mul_binary_add_integral(lhs, rhs, **options):
    """Apply add to sum combination."""
    # Require that both are OperatorKind, BraKind, KetKind, or UndefinedKind
    if lhs.kind != NumberKind and rhs.kind != NumberKind:
        terms = (qapply(term*rhs.function, **options) for term in lhs.args)
        return Integral(Add(*terms), *rhs.limits)


@qapply_Mul.binary.register(Sum, Integral)
def _qapply_mul_binary_sum_integral(lhs, rhs, **options):
    """Apply sum to integral combination."""
    # Require that both are OperatorKind, BraKind, KetKind, or UndefinedKind
    if lhs.kind != NumberKind and rhs.kind != NumberKind:
        return Sum(
            Integral(
                qapply(lhs.function*rhs.function, **options),
                *rhs.limits
            ),
            *lhs.limits
        )


@qapply_Mul.binary.register(Integral, Sum)
def _qapply_mul_binary_integral_sum(lhs, rhs, **options):
    """Apply integral to sum combination."""
    # Require that both are OperatorKind, BraKind, KetKind, or UndefinedKind
    if lhs.kind != NumberKind and rhs.kind != NumberKind:
        return Integral(
            Sum(
                qapply(lhs.function*rhs.function, **options),
                *rhs.limits
            ),
            *lhs.limits
        )


@qapply_Mul.binary.register(Add, Expr, on_ambiguity=ambiguity_register_error_ignore_dup)
def _qapply_mul_binary_add_expr(lhs, rhs, **options):
    """Distribute expression over sum terms."""
    # Require that both are OperatorKind, BraKind, KetKind, or UndefinedKind
    if lhs.kind != NumberKind and rhs.kind != NumberKind:
        result = Add(*(qapply(item*rhs, **options) for item in lhs.args))
        return (result,)


@qapply_Mul.binary.register(Expr, Add, on_ambiguity=ambiguity_register_error_ignore_dup)
def _qapply_mul_binary_expr_add(lhs, rhs, **options):
    """Distribute expression over sum terms."""
    # Require that both are OperatorKind, BraKind, KetKind, or UndefinedKind
    if lhs.kind != NumberKind and rhs.kind != NumberKind:
        result = Add(*(qapply(lhs*item, **options) for item in rhs.args))
        return (result,)


@qapply_Mul.binary.register((Sum, Integral), Expr, on_ambiguity=ambiguity_register_error_ignore_dup)
def _qapply_mul_binary_sum_expr(lhs, rhs, **options):
    """Apply expression to sum/integral function."""
    # Require that both are OperatorKind, BraKind, KetKind, or UndefinedKind
    if lhs.kind != NumberKind and rhs.kind != NumberKind:
        result = lhs.func(qapply(lhs.function*rhs, **options), *lhs.limits)
        return (result,)


@qapply_Mul.binary.register(Expr, (Sum, Integral), on_ambiguity=ambiguity_register_error_ignore_dup)
def _qapply_mul_binary_expr_sum(lhs, rhs, **options):
    """Apply expression to sum/integral function."""
    # Require that both are OperatorKind, BraKind, KetKind, or UndefinedKind
    if lhs.kind != NumberKind and rhs.kind != NumberKind:
        result = rhs.func(qapply(lhs*rhs.function, **options), *rhs.limits)
        return (result,)


@qapply_Mul.binary.register(Expr, Expr, on_ambiguity=ambiguity_register_error_ignore_dup)
def _qapply_mul_binary_expr_expr(lhs, rhs, **options):
    """Default binary handler for Expr pairs - no transformation."""
    return None
