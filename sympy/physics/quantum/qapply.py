"""Logic for applying operators to states.

TODO:

- Finish slidingtransform
- Finish slidingtransform tests
- Add extra type handlers to get rid of warnings in qapply
- Optimize performance
- Finish tests for qapply
- Fix all other tests

"""

from sympy.concrete import Sum
from sympy.core.add import Add
from sympy.core.expr import Expr
from sympy.core.kind import NumberKind
from sympy.core.mul import Mul
from sympy.core.power import Pow
from sympy.core.singleton import S
from sympy.core.sympify import sympify, _sympify
from sympy.functions.elementary.complexes import Abs
from sympy.integrals import Integral
from sympy.multipledispatch import Dispatcher

from sympy.utilities.misc import debug

from sympy.physics.quantum.anticommutator import AntiCommutator
from sympy.physics.quantum.commutator import Commutator
from sympy.physics.quantum.dagger import Dagger
from sympy.physics.quantum.density import Density
from sympy.physics.quantum.innerproduct import InnerProduct
from sympy.physics.quantum.operator import OuterProduct, Operator
from sympy.physics.quantum.state import State, KetBase, BraBase, Wavefunction
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

    def apply_doit(expr, **options):
        sum_doit = options.get('sum_doit', False)
        ip_doit = options.get('ip_doit', True)
        result = f(expr, **options)
        result = _ip_doit_func(result) if ip_doit else result
        result = _sum_doit_func(result) if sum_doit else result
        return result

    return apply_doit


def _flatten_mul(f):
    
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
# Main code
#-----------------------------------------------------------------------------

# Level 0: int, float, complex
# Level 1: Expr
# Level 2: Mul, Add, Pow, Abs
# Level 3: Sum, Integral
# Level 3: Bra, Ket, Operator, Dagger
# Level 4: InnerProduct, OuterProduct
# Level 5: TensorProduct, Density

#-----------------------------------------------------------------------------
# _qapply_unary
#-----------------------------------------------------------------------------

_qapply_unary = Dispatcher('_qapply_unary')


@_qapply_unary.register(Expr)
def _qapply_unary_expr(expr, **options):
    return None


@_qapply_unary.register((int, float, complex))
def _qapply_unary_number(expr, **options):
    return sympify(expr)


@_qapply_unary.register(Abs)
def _qapply_unary_abs(expr, **options):
    debug('Abs: ', expr, **options)
    return Abs(qapply(expr.args[0], **options))


@_qapply_unary.register(Add)
def _qapply_unary_state(expr, **options):
    result = 0
    for arg in expr.args:
        result += qapply(arg, **options)
    return result


@_qapply_unary.register(Pow)
def _qapply_unary_pow(expr, **options):
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
    return Dagger(qapply(expr, **options))


@_qapply_unary.register(Density)
def _qapply_unary_density(expr, **options):
    # For a Density operator call qapply on its state
    new_args = [(qapply(state, **options), prob) for (state,
                     prob) in expr.args]
    return Density(*new_args)


@_qapply_unary.register(TensorProduct)
def _qapply_unary_tp(expr, **options):
    new_args = [qapply(t, **options) for t in expr.args]
    return TensorProduct(*new_args)


#-----------------------------------------------------------------------------
# qapply
#-----------------------------------------------------------------------------

def qapply(e, **options):
    debug('qapply: ', e)
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
    reverse=True
)

# Unary transformations for qapply_Mul


@qapply_Mul.unary.register((int, float, complex))
def _qapply_mul_unary_number(expr, **options):
    return (sympify(expr),)


@qapply_Mul.unary.register(Expr)
def qapply_mul_unary_expr(expr, **options):
    return None


@qapply_Mul.unary.register(Pow)
@_flatten_mul
def _qapply_mul_unary_pow(expr, **options):
    base, exp = expr.as_base_exp()
    if exp.is_Integer and exp > 0:
        return tuple(expr.base for i in range(int(expr.exp)))
    else:
        return (qapply(expr, **options),)


@qapply_Mul.unary.register(OuterProduct)
def _qapply_mul_unary_op(expr, **options):
    return (expr.ket, expr.bra)


@qapply_Mul.unary.register((Commutator, AntiCommutator))
@_flatten_mul
def _qapply_mul_unary_comm(expr, **options):
    return (expr.doit(),)


@qapply_Mul.unary.register(TensorProduct)
@_flatten_mul
def _qapply_mul_unary_tp(expr, **options):
    new_args = [qapply(t, **options) for t in expr.args]
    return (TensorProduct(*new_args),)


@qapply_Mul.unary.register(Density)
def _qapply_mul_unary_density(expr, **options):
    # For a Density operator call qapply on its state
    new_args = [(qapply(state, **options), prob) for (state,
                     prob) in expr.args]
    return (Density(*new_args),)


# Binary transformations for qapply_Mul


@qapply_Mul.binary.register(Expr, Expr)
def _qapply_mul_binary_expr_expr(lhs, rhs, **options):
    return None


@qapply_Mul.binary.register(Operator, KetBase)
@_flatten_mul
def _qapply_mul_binary_op_ket(lhs, rhs, **options):
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


@qapply_Mul.binary.register(Add, Expr)
def _qapply_mul_binary_add_expr(lhs, rhs, **options):
    # Require that both are OperatorKind, BraKind, KetKind, or UndefinedKind
    if lhs.kind != NumberKind and rhs.kind != NumberKind:
        result = Add(*(qapply(item*rhs, **options) for item in lhs.args))
        return (result,)


@qapply_Mul.binary.register(Expr, Add)
def _qapply_mul_binary_expr_add(lhs, rhs, **options):
    # Require that both are OperatorKind, BraKind, KetKind, or UndefinedKind
    if lhs.kind != NumberKind and rhs.kind != NumberKind:
        result = Add(*(qapply(lhs*item, **options) for item in rhs.args))
        return (result,)


@qapply_Mul.binary.register((Sum, Integral), Expr)
def _qapply_mul_binary_sum_expr(lhs, rhs, **options):
    # Require that both are OperatorKind, BraKind, KetKind, or UndefinedKind
    if lhs.kind != NumberKind and rhs.kind != NumberKind:
        result = lhs.func(qapply(lhs.function*rhs, **options), *lhs.limits)
        return (result,)


@qapply_Mul.binary.register(Expr, (Sum, Integral))
def _qapply_mul_binary_expr_sum(lhs, rhs, **options):
    # Require that both are OperatorKind, BraKind, KetKind, or UndefinedKind
    if lhs.kind != NumberKind and rhs.kind != NumberKind:
        result = rhs.func(qapply(lhs*rhs.function, **options), *rhs.limits)
        return (result,)


@qapply_Mul.binary.register(Sum, Sum)
def _qapply_mul_binary_sum_sum(lhs, rhs, **options):
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
