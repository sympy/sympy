"""
Evaluates quantum expressions by applying operators and powers to operators
and states and general expressions by multiplying factors and distributing
over summands. ``qapply`` is based on the generic ``lapply`` ("linear apply"). This
module just provides the ``qapply()`` function interface with additional options
and registers the ``multipledispatch`` handlers required for the quantum objects.
"""

from typing import cast, Any, List, Tuple, Union # For Python 3.8

from sympy.core.expr import Expr
from sympy.core.mul import Mul

from sympy.physics.quantum.lapply import (lapply,                 #-- correct path once location of lapply is final
    # lapply is object agnostic. All object specific information has been
    # isolated in five Dispatcher handler families:
    lapply1_type, lapply1_pow_type, c_nc_ncef,
    lapply_mul2types_A, lapply_mul2types_B,
    # import lapply support functions that may be used in the handlers here
    lapply_mul2, lapply_mul2elems, lapply_no_act)


__all__ = [
    'qapply'
]

#--------+---------+---------+---------+---------+---------+---------+---------
# qapply function header and doc string
#
# Options of lapply are listed for completeness. Options taken from expand()
# have same names (mul, power_xxx, tensorproduct)
#------------------------------------------------------------------------------

def qapply(e:Expr, ip_doit=True, dagger=False, op_join=True, tensorproduct=None,
                   mul=True, power_base=None, power_exp=False, nested_exp=None,
                   apply_exp=False, cache_for_pow=None, **options ) -> Expr:
    """
    Evaluates expressions of quantum operators and applies quantum operators,
    powers and factors to each other. Knows about quantum objects and invokes
    their methods ``._apply_operator`` and ``._apply_from_right_to``. Is based
    on ``lapply`` and may be extended by ``multipledispatch`` handlers.

    Parameters
    ==========

    e : Expr
        The expression containing operators, states resp. factors, summands and
        powers. The expression tree will be walked and operators resp.
        factors and powers applied and distributed over summands.

    dagger : Boolean, optional
        If True then if ``a * b`` doesn't compute, also try
        ``Dagger(Dagger(b)*Dagger(a))``. Defaults to False.

    ip_doit : Boolean, optional
        If True then call ``.doit()`` on inner products when they are
        encountered. Defaults to True.

    op_join : Boolean, optional
        If True then rewrite ``ket * bra`` as outer product ``OuterProduct(ket, bra)``
        and avoid breaking them up to ``ket * bra``. Defaults to True.

    tensorproduct : Boolean, optional
        If True expand sums in factors of ``TensorProduct`` into sums of
        ``TensorProduct``. Same option as in :obj:`sympy.core.function.expand()`.
        Defaults to value of option ``mul``.

    mul : Boolean, optional
        If True (default) distribute commutative factors over sums. Corresponds to option
        ``mul`` of :obj:`sympy.core.function.expand()`. ``mul=False`` will do the minimal
        application only, so use only when expansion is slow.

    power_base : Boolean
        Refers to commutative factors. Same option as in :obj:`sympy.core.function.expand()`:
        If True split powers of multiplied bases. Defaults to value of option ``mul``.

    power_exp : Boolean, optional
        Refers to commutative factors. Same option as in :obj:`sympy.core.function.expand()`:
        It True expand addition in exponents into multiplied bases. Defaults to False.

    nested_exp : Boolean, optional
        If True, then if powers are nested and if their exponents may be multiplied,
        explicitly expand the product of the exponents. Defaults to value of option ``mul``.

    apply_exp : Boolean, optional
        If True apply ``lapply`` to exponents of a power before applying the exponent
        to the base. Default is False and to apply exponents as provided.

    cache_for_pow : dict, optional
        a dictionary that provides cache values to speed up the evaluation
        of powers. If provided it is updated in place, so it can be passed on to the
        next invocation of ``lapply`` as long as options and attributes of symbols
        remain unaltered. Defaults to None, in which case an internal cache is used.

    options : dict, optional
        A dict of keyword/value pairs that are passed to the handlers and the operator
        methods like ``._apply_operator_<rhs>()`` or ``._apply_from_right_to_<lhs>()``
        as specific options.



    Returns
    =======

    e : Expr
        The original expression, but with the operators applied to operators and states
        and with moderate expansion done by distributing commutative factors over summands.

    Examples
    ========

        >>> from sympy.physics.quantum import qapply, Ket, Bra
        >>> from sympy.physics.quantum.qubit import QubitBra
        >>> from sympy.physics.quantum.gate import XGate
        >>> from sympy.core.symbol import symbols
        >>> b = Bra('b')
        >>> k = Ket('k')
        >>> A = k * b
        >>> A
        |k><b|
        >>> qapply(A * b.dual / (b * b.dual))
        |k>
        >>> qapply(k.dual * A / (k.dual * k))
        <b|
        >>> qapply(QubitBra(1) * XGate(0), dagger = True)
        <0|
        >>> qapply(A * k * b * A, op_join = False)
        <b|k>**2*|k>*<b|
        >>> n = symbols("n", integer=True, nonnegative=True)
        >>> qapply(A ** (n + 2))
        <b|k>**(n + 1)*|k><b|

    See Also
    ========

    sympy.physics.quantum.lapply.lapply: Generic evaluation engine to apply linear operators, powers and factors

"""                                  # ----- adapt path to lapply above once its position in code tree is determined !! ---

    # set quantum specific options for qapply
    # if a*b doesn't compute, try Dagger(Dagger(b)*Dagger(a))
    options['dagger'       ] = dagger
    # evaluate ip.doit() on all InnerProduct objects
    options['ip_doit'      ] = ip_doit
    # try to join Ket*Bra to OuterProduct, avoid breaking OuterProduct
    options['op_join'      ] = op_join
    # expand sums in tensor product factors (= hint tensorproduct of expand())
    if tensorproduct is None: tensorproduct = mul
    options['tensorproduct'] = tensorproduct

    # invoke lapply passing lapply options on
    res = lapply(e, mul=mul, power_base=power_base, power_exp=power_exp,
        nested_exp=nested_exp, apply_exp=apply_exp,
        cache_for_pow=cache_for_pow, **options)
    return res


# quantum package imports at this point avoid circular import with Density
from sympy.physics.quantum.qexpr import QExpr
from sympy.physics.quantum.density import Density
from sympy.physics.quantum.operator import Operator, OuterProduct
from sympy.physics.quantum.innerproduct import InnerProduct
from sympy.physics.quantum.tensorproduct import TensorProduct
from sympy.physics.quantum.commutator import Commutator
from sympy.physics.quantum.anticommutator import AntiCommutator

from sympy.functions.elementary.complexes import adjoint
from sympy.physics.quantum.dagger import Dagger

##    ##  ---------------------------------------------------------------------
##    ##
##    ##  @dispatch Handlers for lapply_mul2types_*(lhs, rhs, **options)
##    ##
########  See lapply.py: ------------------------------------------------------
########  Handlers may return result expr, or None, or, if an interaction of
##    ##  lhs and rhs is confirmed, a list of factors suitable as input to
##    ##  first argument of lapply_mul2(). If the handler is able to determine
##    ##  that lhs and rhs will not interact for mathematical reasons is should
##    ##  return Mul(lhs, rhs, evaluate=False)

@lapply_mul2types_A.register(TensorProduct, TensorProduct)
def hdlrTT_for(lhs:TensorProduct, rhs:TensorProduct, **options) \
                                                -> Union[Expr, List, None]:
    # class TensorProduct doesn't define .__mul__  or *
    if len(lhs.args) == len(rhs.args): # we can't check factor dim's !!
        # Do the *-operator for TensorProduct*TensorProduct. If there are
        # c-factors in the tensor factors TensorProduct will pull them
        # up front automatically and return a Mul:
        res = TensorProduct(*[ cast(Expr, lf) * cast(Expr, rf)
                                for (lf, rf) in zip(lhs.args, rhs.args)])
        return res
    else: # Currently we cannot handle product with different factor numbers
        return None  # as dim's of factor spaces are not known # type: ignore[no-redef]

@lapply_mul2types_A.register((Operator, OuterProduct), Density)
def hdlrOOD_for(lhs:Operator, rhs:Density, **options) \
                                            -> Union[Expr, List, None]:
    # class Density only defines Op*Density:=Density.apply_op(Op)
    # We allow only (unitary) operators as lhs as .apply_op is unable
    # to cope with anything else. Rewrite this function if needed.
    return rhs.apply_op(lhs)

# This dispatch types are put in layer B to resolve type ambiguities
# with the type pairs in layer A, prioritising layer A:
# e.g. (Pow, Expr) and (Expr, OuterProduct) would be ambiguous for
# (Pow, OuterProduct) and require a very special rule for
# (Pow, OuterProduct) which has no real meaning.
# For input types (QExpr, OuterProduct) dispatch chooses
# (QExpr, (QExpr, Expr)) over (Expr, OuterProduct) without warning!
# So use dispatch ((QExpr, Expr), OuterProduct) to fix that.

@lapply_mul2types_B.register((QExpr, Expr), OuterProduct)
def hdlrEO_for(lhs:Expr, rhs:OuterProduct, **options) \
                                            -> Union[Expr, List, None]:
    res = lapply_mul2elems(lhs, rhs.ket, **options)
    if type(res) is list:   # assume interaction then
        return res + [rhs.bra]
    elif lapply_no_act(lhs, rhs.ket, res): # type: ignore[arg-type]
        return Mul(lhs, rhs, evaluate = False) # return unevaluated
    else:
        return [res, rhs.bra] # chose return factors as list # type: ignore[no-redef]

@lapply_mul2types_B.register(OuterProduct, (OuterProduct, QExpr, Expr))
def hdlrOE_for(lhs:OuterProduct, rhs:Expr, **options) \
                                            -> Union[Expr, List, None]:
    res = lapply_mul2elems(lhs.bra, rhs, **options)
    if type(res) is list:   # assume interaction then
        return [lhs.ket] + res
    elif lapply_no_act(lhs.bra, rhs, res): # type: ignore[arg-type]
        return Mul(lhs, rhs, evaluate = False) # return unevaluated
    else:
        return [lhs.ket, res] # chose return factors as list # type: ignore[no-redef]


# The catch-all for quantum object * Expr or Expr * quantum object:
# Invokes the _apply_operator_xx and _apply_from_right_to methods
# that are defined for many quantum objects. We cannot restrict to
# QExpr x QExpr here, as e.g. Wavefunction is a quantum object but
# no QExpr. Also interaction between QExpr and Expr should be
# supported in general.

@lapply_mul2types_B.register(QExpr, (QExpr, Expr))
def hdlrQQ_for(lhs:QExpr, rhs:Expr, **options) -> Union[Expr, List, None]:
    # The quantum package uses methods _apply_operator and _apply_from_right
    # methods for many of its objects to describe application. So for quantum
    # objects try to apply _apply_operator and _apply_from_right_to as a
    # catch-all. multipledispatch will give priority to more specific rules.
    res = cast(Any, None)
    if hasattr(lhs, '_apply_operator'):
        try: # for historical reasons, _apply_operator is tried first
            res = cast(Any, lhs._apply_operator(rhs, **options)) # type: ignore[attr-defined]
        except (NotImplementedError): # action lhs not implemented for rhs
            pass
    if res is None and isinstance(rhs, QExpr):
        # Invoke the catch-all function for Expr*QExpr
        res = (lapply_mul2types_B.dispatch(Expr, QExpr))(lhs, rhs, **options)
    if res is None and options["dagger"]:
        res = lapply_mul2daggered(lhs, rhs, **options)
    return res # either some result or None # type: ignore[no-redef]

@lapply_mul2types_B.register(Expr, QExpr) # sides (QExpr, Expr) switched
def hdlrEQ_for(lhs:Expr, rhs:QExpr, **options) -> Union[Expr, List, None]:
    res = cast(Any, None)
    if hasattr(rhs, '_apply_from_right_to'):
        try:
            res = cast(Any, rhs._apply_from_right_to(lhs, **options)) # type: ignore[attr-defined]
        except (NotImplementedError): # action rhs not implemented for lhs
            pass
    if res is None and options["dagger"]:
        res = lapply_mul2daggered(lhs, rhs, **options)
    return res # either some result or None # type: ignore[no-redef]

# If Layer A had no result for QExpr, QExpr, and options["dagger"]==True,
# try identity lhsR*rhsL = Dagger(Dagger(lhsR*rhsL))
#                        = Dagger(Dagger(rhsL)*Dagger(lhsR))
# as the last resort. Returns Mul(lhs, rhs, evaluate=False) if no success.
# Dagger is built on class functions.elementary.complexes.adjoint, which
# is also used in matrices.Adjoint, so it is actually a general class and
# independent of package quantum where it is located! However the option
# "dagger" makes sense only in packages like quantum which has one-sided
# definition of 'apply' resp. *, i.e. where A*B may be undefined, but
# Dagger(B)*Dagger(A) is defined.

def lapply_mul2daggered(lhs:Expr, rhs:Expr, **options) \
                                                 -> Union[Expr, None]:
    #assert options["dagger"] # precondition for this function
    res0 = lhs * rhs # try the standard last resort first
    if lapply_no_act(lhs, rhs, res0):
        res0 = Mul(lhs, rhs, evaluate=False)
    else: return res0
    # Return if the daggered values are just formal adjoints
    rhs_d = Dagger(rhs)
    if isinstance(rhs_d, adjoint): return res0
    lhs_d = Dagger(lhs)
    if isinstance(lhs_d, adjoint): return res0
    # If daggered values both exist try to mul them
    options["dagger"] = False # don't recurse infinitely.
    res_ud = lapply_mul2elems(rhs_d, lhs_d, **options)
    options["dagger"] = True
    # check for interaction before doing un-dagger as this
    # may disturb recognition of no interaction
    if type(res_ud) is list: # so assume interaction ok
        return Dagger(lapply_mul2(res_ud, [], **options))
    elif lapply_no_act(rhs_d, lhs_d, res_ud): # type: ignore[arg-type]
        return Mul(lhs, rhs, evaluate=False)
    else: # successful
        return Dagger(res_ud)


##    ##  ---------------------------------------------------------------------
##    ##
##    ##  @dispatch Handlers for lapply1_pow_type(base, exp, **options)
##    ##
########  see lapply.py: ------------------------------------------------------
########  This handlers are invoked by the lapply1_type(Pow) to handle a power
########  Pow(base, exp) with a particular type of base and/or exponent.
##    ##  Handlers should take care of non-evaulated objects and recurse into
##    ##  inner structure of e, if any, and process it so that the result is
##    ##  flattened and c_nc_ncf() may extract factors from the result in a
##    ##  next step.

@lapply1_pow_type.register()
def hdlrTP_for(base:TensorProduct, exp:Expr, **options): # -> Expr:
    # quantum TensorProduct doesn't provide .__pow__ or ._pow, and
    # we don't use tensor_product_simp_Pow(e).
    # Note: Using ** will pull out commutative factors and return a
    # Mul(Pow(c-factors, exp), Pow(nc-part, exp)).
    # And TensorProduct will also pull out all commutative factors from
    # its elements, returning Mul in that case
    return TensorProduct(*[lapply_mul2([b**exp], [], **options)
                                    for b in base.args]        )


##    ##  ---------------------------------------------------------------------
##    ##
##    ##  Handlers for lapply1_type(Expr, to_L=..., **options)
##    ##
##    ##  See lapply.py: ------------------------------------------------------
########  Handlers should take care of non-evaulated objects and recurse into
########  inner structure of e, if any, and process it so that the result is
##    ##  flattened and c_nc_ncf() may extract factors from the result in a
##    ##  next step.
##    ##  Note: These handlers are also called for commutative factors! In
##    ##  most cases no action is wanted as lapply leaves c factors to Mul.
##    ##  So put "if e.is_commutative: return e" on top.


@lapply1_type.register(InnerProduct) #call as lapply1_type(e, **options)
def hdlrI_for(e:InnerProduct, **options) -> Expr:
    return cast(Expr, e.doit() if options["ip_doit"] else e)

@lapply1_type.register(Density) #call as lapply1_type(e, **options)
def hdlrD_for(e:Density, **options) -> Expr:
    new_args = [(lapply_mul2([cast(Expr, state)],[], **options),
                 cast(Expr, prob))  for (state, prob) in cast(Tuple, e.args)]
    return cast(Expr, Density(*new_args)) # cast for MyPy # type: ignore[no-redef]

@lapply1_type.register(TensorProduct) #call as lapply1_type(e, **options)
def hdlrT_for(e:TensorProduct, **options) -> Expr:
    # call lapply on args of TensorProduct, may return TP or Mul
    res = cast(Expr, TensorProduct(*[lapply_mul2([t], [], **options)
                                              for t in e.args ]))
    # if res is TensorProduct and option tensorproduct set, expand is
    if options["tensorproduct"] and isinstance(res, TensorProduct):
        res = res._eval_expand_tensorproduct() # may become Mul
    # TensorProduct is treated as elementary even if it is of the form
    # (A*B)x(C*D) that would allow factorisation (AxC)*(BxD)
    return cast(Expr, res) # cast for MyPy  # type: ignore[no-redef]

@lapply1_type.register((Commutator, AntiCommutator))
def hdlrCAC_for(e:Union[Commutator, AntiCommutator], **options) -> Expr:
    # .doit() in any case and loose structure. May return any type
    return cast(Expr, e.doit()) # cast for MyPy # type: ignore[no-redef]


##    ##  ---------------------------------------------------------
##    ##
##    ##  Handlers for c_nc_ncef(Expr, to_L=..., **options)
##    ##
##    ##  See lapply.py: ------------------------------------------
########  Split up expression e in commutative factors and non-
########  commutative factors. From the nc factors pick the leftmost
##    ##  (if to_L==True) or rightmost (to_L==False) elementary factor
##    ##  ncef and return (c_list, nc_list_without_ncef, ncef). If
##    ##  e is commutative, return (c_list, [], S.One). If e cannot
##    ##  be split up, return ([], [], e).

@c_nc_ncef.register(OuterProduct) # Call as c_nc_ncef(e, to_L=..., **options)
def hdlrOP_for(e:OuterProduct, to_L:bool, **options) -> Tuple[List, List, Expr]:
    if options["op_join"]:  # option to not break up OuterProducts
        return ([], [], e)  # so consider OuterProduct elementary
    else:    # else break it up into ket * bra
        return (     ([], [e.bra], e.ket) if to_L
                else ([], [e.ket], e.bra)        )
