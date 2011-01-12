"""Logic for applying operators to states."""

import copy

from sympy import S, Add, Mul, Pow

from sympy.physics.quantum.qexpr import split_commutative_parts
from sympy.physics.quantum.dagger import Dagger
from sympy.physics.quantum.operator import Operator
from sympy.physics.quantum.state import KetBase, BraBase

__all__ = [
    'apply_operators'
]


#-----------------------------------------------------------------------------
# Main code
#-----------------------------------------------------------------------------

def apply_operators(e, **options):
    """Apply operators to states in a quantum expression.

    By default, operators acting on states (O|psi>) are left in symbolic,
    unevaluated form. This function uses various special methods to attempt
    and apply operators to states. When this happens there are two possible
    outcomes: i) it is not known how the operator acts on the state, which
    will result in the expression being unchanged and ii) it is known how the
    operator acts on the sate, which will result in the action being carried
    out.

    Parameters
    ==========
    e : Expr
        The expression containing operators and states. This expression tree
        will be walked to find operators acting on states symbolically.
    options : dict
        A dict of key/value pairs that determine how the operator actions
        are carried out.

    Returns
    =======
    e : Expr
        The original expression, but with the operators applied to states.
    """

    # TODO: Fix early 0 in apply_operators.

    # This may be a bit aggressive but ensures that everything gets expanded
    # to its simplest form before trying to apply operators. This includes
    # things like (A+B+C)*|a> and A*(|a>+|b>) and all Commutators and
    # TensorProducts. The only problem with this is that if we can't apply
    # all the Operators, we have just expanded everything.
    # TODO: don't expand the scalars in front of each Mul.
    e = e.expand(commutator=True, tensorproduct=True).doit()

    # If we just have a raw ket, return it.
    # TODO: make this acts_like_KetBase
    if isinstance(e, KetBase):
        return e

    # We have an Add(a, b, c, ...) and compute 
    # Add(apply_operators(a), apply_operators(b), ...)
    elif isinstance(e, Add):
        result = 0
        for arg in e.args:
            result += apply_operators(arg, **options)
        return result

    # We have a Mul where there might be actual operators to apply to kets.
    elif isinstance(e, Mul):
        return apply_operators_Mul(e, **options)

    # In all other cases (State, Operator, Pow, Commutator, InnerProduct,
    # OuterProduct) we won't ever have operators to apply to kets.
    else:
        return e


def apply_operators_Mul(e, **options):

    # The general form of e at this point is:
    #
    # ZeroOrMore(scalars)*ZeroOrMore(Bra)*ZeroOrMore(Operators)*Ket*
    # ZeroOrMore(Bra*ZeroMore(Operators)*Ket)
    # 
    # We want to pick out pieces that look like:
    #
    # Bra*ZeroOrMore(Operators)*Ket

    if not isinstance(e, Mul):
        return e

    # Split the Mul into a Mul of scalars and a list of non-commutative parts.
    c_part, nc_part = split_commutative_parts(e)
    c_part = Mul(*c_part)
    n_nc = len(nc_part)
    if n_nc == 0 or n_nc == 1:
        return c_part

    terms = []
    last_ket = 0
    for i in range(n_nc):
        # Make this acts_like_KetBase
        if isinstance(nc_part[i], KetBase):
            terms.append(nc_part[last_ket:i+1])
            last_ket = i+1
    last_part = nc_part[last_ket:]
    if last_part:
        terms.append(last_part)

    result = 1
    for term in terms:
        # print 'term', term

        bra = 1
        # TODO: make this acts_like_BraBase
        if isinstance(term[0], BraBase):
            bra = term.pop(0)

        ket = 1
        # TODO: make this acts_like(term[-1], KetBase)
        if isinstance(term[-1], KetBase):
            ket = term.pop(-1)

        # print 'bra', bra
        # print 'ket', ket

        if bra == 1 and ket == 1:
            tresult = Mul(*term)

        elif len(term) == 0:
            # This cover the cases where either bra or ket is 1 and will also
            # automatically create an InnerProduct if neither are 1.
            tresult = bra*ket

        # The <bra|*A*B case
        elif ket == 1:
            try:
                new_term = [Dagger(i) for i in reversed(term)]
            except NotImplementedError:
                tresult = bra*Mul(*term)
            else:
                # TODO: what is the best way of handling copy.copy in all this
                tresult = Dagger(
                    apply_list_of_ops(1, copy.copy(new_term), Dagger(bra))
                )

        # The full <bra|*A*B*|ket> case
        # TODO: Get scalar<bra|*|ket> into an InnerProduct.
        else:
            try:
                tresult = bra*apply_list_of_ops(1, copy.copy(term), ket)
            except NotImplementedError:
                tresult = bra*Mul(*term)*ket

        # print 'tresult', tresult
        # print
        result *= tresult

    return c_part*result


def apply_list_of_ops(scalar, ops, ket):
    # scalar is a number
    # ops is a list of Operators or Powers of Operators.
    # ket can be a*|alpha> + b*|beta> + ... but cannot include any Operators.

    if scalar == 0 or ket == 0:
        return S.Zero

    # print scalar, ops, ket

    if not ops:
        return scalar*ket

    if isinstance(ket, Add):
        result = 0
        for arg in ket.args:
            if isinstance(arg, KetBase):
                result += apply_list_of_ops(scalar, ops, arg)
            elif isinstance(arg, Mul):
                k = arg.args[-1]
                s = Mul(scalar, *arg.args[:-1])
                result += apply_list_of_ops(s, copy.copy(ops), k)
            else:
                raise TypeError('Ket or Mul expected, got: %r' % arg)
        return result

    if isinstance(ket, Mul):
        scalar = Mul(scalar, *ket.args[:-1])
        ket = ket.args[-1]

    next_op = ops.pop(-1)
    if isinstance(next_op, Pow):
        ops.append(next_op.base**(next_op.exp-1))
        next_op = next_op.base

    # If unsuccessful, this will raise NotImplementedError which we let
    # propagate.
    new_ket = apply_single_op(next_op, ket)

    if new_ket == 0:
        return 0

    # TODO: Try is without this expand, I don't think we will need it.
    # We definitely need it!!!
    new_ket = new_ket.expand()

    return apply_list_of_ops(scalar, ops, new_ket)


def apply_single_op(op, ket):
    # Apply a single Operator to a Ket, at this point, we must have only an
    # Operator and a Ket and nothing else!
    if not isinstance(ket, KetBase):
        raise TypeError('Ket expected, got: %r' % ket)
    if not isinstance(op, Operator):
        raise TypeError('Operator expected, got: %r' % op)
    try:
        result = op._apply_operator(ket)
    except NotImplementedError:
        try:
            result = ket._apply_operator(op)
        except NotImplementedError:
            raise
    return result
