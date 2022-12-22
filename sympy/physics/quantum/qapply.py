# qapply - evaluate expressions of quantum expressions
#
# See e.g. https://github.com/sympy/sympy/issues/17665 for a discussion of qapply
# See e.g. https://github.com/sympy/sympy/issues/19540 for a discussion of InnerProduct and TensorProduct

# Design Principles of this implementation:

#  1. The objects and notation of quantum mechanics (Ket, Bra, ...) are just special cases of
#     general objects (vector, linear operator,..) and the notation of standard Linear Algebra.
#     So I built the algorithm work with generic/symbolic/abstract vectors, linear operators etc.
#     from Linear Algebra in mind. It should be possible extended it to standard Linear Algebra
#     with little effort just by providing object type specific rules in 3. and 4.
#     As can be seen from the imports, the only knowledge required from quantum is about
#     its composite types. No knowledge about states, operators, wavefunctions, qubits etc.

#  2. Centralize the algorithm that walks the expression tree and determines to pair factors a,b,c,d
#     as (a*b)*(c*d) or a*((b*c)*d) etc. in one place (qapply_Mul2) with no object type specifics
#     in it. I decided to implement (a*(b*(c*d))) in general.
#
# This requires to
#  3. have all rules how to break up composite objects in factors made explicit and grouped
#     in one place (split_c_nc_atom)
#  4. have all inbuilt rules how to multiply two object types made explicit and grouped in one place
#     (qapply_Mul2Atoms)
#
# Moreover:
#  5. do not use expand() indiscriminately: Only do an expansion of Add if required to apply operators, so
#     keeping expansion to the bare minimum. As default, distribute commutative factors over sums as this allows
#     terms to cancel out thus simplifying the expression. In order to simplify commutative factors qapply()
#     employs the three mayor hint functions of expand() Mul._eval_expand_mul(), Pow._eval_expand_power_base()
#     and Pow._eval_expand_power_exp(). See options mul, power_base and power_exp.
#  6. apply recursion along the structure of the expression but don't use recursion for iterations:
#     this avoids hitting the Python recursion limit in very complex expressions or powers with
#     large exponents
#  7. No need to try to maintain subclass types of objects and expressions, as we simply cannot know whether
#     a product is still in a subclass of Mul etc. So use the top types.

# Most rules might at a later time be implemented as methods of the object classes, however I decided
# to be minimally invasive into the existing code base and to implement each set of rules locally
# in this file, so this file is a 1:1 replacement of the previous version of qapply.py

# Note: The types OuterProduct resp. InnerProduct are sympy.physics.quantum.operator.OuterProduct/InnerProduct
#       and by definition are restricted to hold subclasses of KetBase and BraBase. Not to be confused
#       with other implementations of same or similar mathematical concepts elsewhere in SymPy (e.g.
#       .matrices or .physics.vector deal with vectors/matrices described by coefficients.
#       Classical outer/inner product is implemented by sympy.physics.vector.functions.outer resp.
#       sympy.matrices.dot/sympy.physics.vector.functions.dot).


"""Evaluates quantum expressions by applying operators to operators and states.
The expressions are expanded as far as required to apply factors to summands.
"""

from typing import cast, Any, List, Tuple # keep compatibility with Python 3.8
from sympy.core.expr import Expr
from sympy.core.add import Add
from sympy.core.mul import Mul
from sympy.core.power import Pow
from sympy.core.singleton import S
from sympy.core.numbers import Number
from sympy.functions.elementary.integers import floor

from sympy.physics.quantum.dagger import Dagger

# Composite types from .quantum that qapply needs to break up or compose:
from sympy.physics.quantum.anticommutator import AntiCommutator
from sympy.physics.quantum.commutator import Commutator
from sympy.physics.quantum.innerproduct import InnerProduct
from sympy.physics.quantum.operator import OuterProduct
from sympy.physics.quantum.tensorproduct import TensorProduct
#from sympy.physics.quantum.density import Density # moved to qapply body
#otherwise generates circular import reference qapply->Density->represent->qapply



__all__ = [
    'qapply'
]


#-----------------------------------------------------------------------------
# Main code
#-----------------------------------------------------------------------------

def qapply(e, **options):
    """Apply operators to operators and states in a quantum expression.

    Parameters
    ==========

    e : Expr
        The expression containing operators and states. This expression tree
        will be walked to find operators acting on operators and states symbolically.
    options : dict
        A dict of key/value pairs that determine how the operator actions
        are carried out. The options are passed on to invocations of operator methods
        <lhs>._apply_operator_<rhs>() and <rhs>._apply_from_right_to_<lhs>().

        The following options are defined for qapply itself:

        * ``dagger``: if a * b doesn't compute, also try ``Dagger(Dagger(b)*Dagger(a))``
          (default: False).
        * ``ip_doit``: call ``.doit()`` in inner products when they are
          encountered (default: True).
        * ``op_join``: rewrite ket * bra as outer product ``OuterProduct(ket, bra)``
          and avoid breaking them up to ket * bra (default: True).
        * ``mul``: distribute commutative factors over sums. Corresponds to option ``mul``
          of ``expand``. ``mul=False`` will return application only (default: True).
        * ``power_base``: for commutative factors, split powers of multiplied bases
          (same option as in ``expand``) (default: value of option ``mul``).
        * ``power_exp``: for cummutative factors, expand addition in exponents into
          multiplied bases (same option as in ``expand``) (default: value of option ``mul``).
        * ``tensorproduct``: Expand sums in factors of ``TensorProducts`` into sums of
          ``TensorProducts`` (same option as in ``expand``) (default: True).



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
        >>> qapply(A ** (n + 2), power_exp=False)
        <b|k>**(n + 1)*|k><b|

    """
    # moved import here to avoid circular reference on inport qapply->Density->represent->qapply
    from sympy.physics.quantum.density import Density

    # Note: The variables EmptyMul, options, dagger, ip_doit, op_join
    # are assigned and in read-only scope: see bottom of file for body of qapply()

    """
    Remarks on quantum type constructors and Mul
    --------------------------------------------

    1.) The code relies on the auto-simplifications that Mul and all the composite objects constructors like
    InnerProduct, OuterProduct, TensorProduct, Pow etc do: Extraction of commutative factors, contracting
    of same factors to Powers etc, so these constructors may return their type, Mul or even other types.

    2.) Mul(a,b) != a * b:
    Mul(a, b, c) processes the arguments a, b, c (e.g. factoring out scalars, grouping same factors to powers)
    but does not try to actually multiply the arguments, i.e. is not trying (a*(b*c)), (a*b) etc.
    On the other hand a*b is executed using *-operator __mul__ for type(a), type(b):
    E.g. for Kets and Bras, __mul__ is overloaded to form InnerProducts resp. OuterProducts types, and
    e.g. IdentityOperator()*Ket(1)-> Ket(1)
    (Examples: Mul(Bra(1), Ket(1)) returns unchanged. Bra(1)*Mul(Ket(1), Bra(1)) returns
    Mul(Bra(1), Ket(1), Bra(1)), while Bra(1)*Ket(1)*Bra(1) returns Mul(InnerProduct(Bra(1),Ket(1)), Bra(1)).
    So while the return type of Mul(a,b) will be Mul (or Zero) the return type of a*b completely depends
    on the input types.

    3.) Performance of Mul: refer to the remarks in code of Mul.flatten:
    - avoid a*b*c if you don't need the multiplication and use Mul(a,b,c) instead
    - instead of Mul(a,b)*Mul(c) use Mul(a,b,c)
    - results of Mul(a,b,c) are cached


    Remarks on qapply vs. Mul and "atomic" factors:
    ----------------------------------------------

    I will use the verbs 'to apply' and 'to multiply' equivalently as in this context they
    mean the same and are both denoted by *.

    So the purpose of qapply is to multiply factors that Mul doesn't know how to multiply, and
    leave those that Mul knows to Mul. This necessarily includes the distribution of factors
    over Add summands as far as required, which forces qapply to be a kind of workalike of expand().

    For this purpose we have to classify all types of factors in two ways:

    1. Does Mul know how to multiply this type with all other types, or does it not?

       For efficiency, qapply would try to only handle those types that Mul doesn't
       know how to multiply, and will leave all others to Mul. This also helps to avoid
       contradicting results of Mul and qapply.
       Of course passing factors to Mul makes only sense if there are more than one
       to pass. So commutative factors are best. So it fits just nicely that in the
       quantum context the only types Mul knows how to multiply are the scalars, i.e.
       the complex numbers that may be pulled out.

       So we would best use 'is_scalar' as predicate to send factors to Mul. Unfor-
       tunately that doesn't exist, and the closest approximation 'is_complex' can
       not be used: on most symbols, is_complex is undefined. The only predicate that
       SymPy guarantees is 'is_commutative'!

       So we use 'is_commutative', and we can even profit from args_cnc(). Currently
       there is no need to filter out commutative multiples of the identity operators,
       as Identity.is_commutative==False and IdentityOperator.is_commutative==False.
       But if this quirks should vanish, a better filter is required.

    2. Is the factor a compound type (e.g. OuterProduct, Commutator) that is actually
       an expression that should be evaluated and broken up into its components to be
       multiplied one by one, or is it 'elementary'?
       So qapply implements rules which compound types to break up ((Anti-)Commutator,
       Pow, Mul, Add) and which to break as an option (OuterProduct). Those that will
       not be broken up are called "atomic" (note that this is similar in meaning but
       not identical to the sympy definition of atom/atomic).

       Note that "atomic" refers to the general algorithm (qapply_Mul2) that handles
       walking the expression tree and the association of factors. qapply still may in-
       spect the inner contents of factors to do the multiplication (in qapply_Mul2Atoms).
    """

    def no_interaction(lhs:Expr, rhs:Expr, res:Expr) -> bool:
        """Returns True if lhs, rhs didn't interact so lhs*rhs just gave res=Mul(lhs,rhs)
           or Pow(lhs, 2) in case lhs == rhs""" # assumes lhs, rhs non-commutative!
        return ((res.is_Mul and res.args == (lhs, rhs)) or
                (lhs == rhs and res.is_Pow and res.base == lhs))

    # -----------
    # qapply_Mul2
    # -----------
    # Multiply all factors from fL to fR from the left
    # - fL is typically derived by list(Mul.args), but may contain any expression incl. Mul, Add and c arguments
    # - fR may not contain Add, Mul or factors Mul should multiply, may only contain "atomic" factors (see above)
    # - Important: fL, fR are modified during execution!

    def qapply_Mul2(fL:List, fR:List) -> Expr:
        """Multiplies the left hand list of factors into the right hand list of factors"""

        # define shorthand to distribute eadd: cl*fL*(eadd.args[0]+eadd.args[1]+..)*fR
        def qapply_DoAddL(eadd:Add): # Important: arguments to qapply_Mul2 must be fresh copies!
                                     # TBD: implement a cache to avoid frequent evaluation of fL?
            cle = qapply_expand_powers(cl)  # do expansion of power_base / power_exp
            if expand_mul: # distribute factors over summands, but multiply each summand with cl
                c = qapply_expand_mul(cle)  # deeply distribute factors over sums in cle
                return Add(*[qapply_Mul2(fL + [c * arg], fR.copy()) for arg in eadd.args])
            else: # distribute factors over summands, but keep commutative factors up front
                return Mul(*(cle + [Add(*[qapply_Mul2(fL + [arg], fR.copy()) for arg in eadd.args])] ))

        cl : List[Expr] = [] #gathers factors for Mul. Avoids calling Mul for any new factor.

        lhsR = EmptyMul # init leftmost  left  hand side factor
        rhsL = EmptyMul # init rightmost right hand side factor

        while(True):         # replacing recursion on same level by iteration

            if lhsR is EmptyMul:     # Need to pick factor from left hand side
                if fL == []:         # but there aren't any more
                    if rhsL is not EmptyMul:
                        fR.insert(0, rhsL) # push rhsL factor not yet in fR
                    cl = qapply_expand_powers(cl) # do expansion of power_base / power_exp
                    return Mul(*( ([qapply_expand_mul(cl)] if expand_mul else cl) + fR)) # return final result
                # else
                lhs0 = fL.pop()      # get right-most factor from fL
                # split lhs0 into c factors list lhs_c, left factor list lhsL, right atom factor lhsR
                # so that Mul(*lhsL)*lhsR == lhs  and  Mul(*lhs_c)*lhs == lhs0, also applying .doit()s!
                lhs_c, lhsL, lhsR = split_c_nc_atom(lhs0, to_the_right)
                cl.extend(lhs_c)     # save c factors, if any
                if lhsR is EmptyMul: # lhs was all for Mul, so went into lhs_c and cl completely
                    continue         # so pick next factor
                if lhsR is S.Zero:
                    return S.Zero
                fL.extend(lhsL)      # push back lhsL factors, if any. A Mul is flattened into fL

                if lhsR.is_Add:      # do expansion of Add
                    if rhsL is not EmptyMul:
                        fR.insert(0, rhsL)  # push back rhsL factor if any
                    return qapply_DoAddL(lhsR) # distribute factors over summands of lhsR

            # Now lhsR is a non-trivial nc factor and no Add. Get rhsL:
            if rhsL is EmptyMul:
                if fR == []:        # no factors in fR
                    rhsL = lhsR     # so shift lhsR over to become the rhsL, is atomic
                    lhsR = EmptyMul # request new lhsR instead
                    continue
                # else: get next left most factor from fR
                # fR is built by atomic nc factors, no Add, Mul or c factors
                rhsL = fR.pop(0)    # atomic, as only atoms are pushed to fR

            # the heart: we multiply two atomic factors with no obvious scalar factors within
            res0 = qapply_Mul2Atoms(lhsR, rhsL)  # this may return ANY type!
            # try dagger if the two factors have not interacted at all:
            nia = no_interaction(lhsR, rhsL, res0)
            if nia and dagger: #Identity lhsR*rhsL = Dagger(Dagger(lhsR*rhsL)) = Dagger(Dagger(rhsL)*Dagger(lhsR))
                res0 = Dagger(rud := qapply_Mul2Atoms(rd := Dagger(rhsL), ld := Dagger(lhsR)))
                nia = no_interaction(rd, ld, rud) # check before un-dagger that may do transformations

            if nia:  # the two factors lhsR*rhsL have not interacted
                # so continue with a fresh factor from the left: fL*lhsR
                fR.insert(0, rhsL) # push rhsL (is atomic) to fR
                rhsL = lhsR        # continue with previous lhsR as new rhsL (is atomic)
                lhsR = EmptyMul    # request a new lhsR from fL
            else: # Multiplication had interaction lhsR*rhsL->res0, so continue with res0*fR
                res_c, resL, resR = split_c_nc_atom(res0, to_the_right) # incl. ip_doit, Commutator->Add, -> ANY type in resL
                if resR is S.Zero: # the case resR==1 handled below
                    return S.Zero
                cl.extend(res_c)   # save away any factors for Mul

                # if resR is Add (implies resR != rhsL), distribute it
                if resR.is_Add:
                    fL.extend(resL)  # push resL factors
                    return qapply_DoAddL(resR) # distribute factors over summands of resR

                if resR == rhsL:
                    if resL == [lhsR]:  # implies we have broken up res0 into its original factors
                        # so push rhsL to fR and continue with fL*lhsR
                        fR.insert(0, rhsL) # = resR, push rhsL (is atomic) to fR
                        rhsL = lhsR        # = resL, continue with previous lhsR as new rhsL (is atomic)
                        lhsR = EmptyMul    # request a new lhsR
                    else: # resL must be equal to lhsR, but in different representation, so compute resL*rhsL
                        fL.extend(resL)    # Note: if lhsR was an Identity, resL == 1 will be eliminated in push
                        lhsR = EmptyMul    # get new factor from fL
                else: # lhsR*rhsL -> res_c*resL*resR, so continue with resR*fR
                    fL.extend(resL)
                    lhsR = resR            # atomic; if resR==1==EmptyMul, request new factor from fL
                    rhsL = EmptyMul        # get new factor from fR


    # qapply_expand_powers() and qapply_expand_mul() implement subsets of expand(), applied
    # here to commutative factors only. In order to simplify some expressions qapply needs
    # to deeply expand and multiply out all commutative Add-factors to make
    # summands cancel out. Also expansion of powers via expand(power_base=True) or
    # expand(power_exp=True) is enabled by default, as it is in expand().
    # On the other hand we do not want to call expand(deep=True), as this also
    # aggressively expands function arguments like sqrt(j*(j+1)+m*(m+1)).
    # So we use the parts of expand that have greatest effect: the hint functions
    # Mul._eval_expand_mul and Pow._eval_expand_power_xxx:

    def qapply_expand_powers(cl:List) -> List:
        """Expands commutative factors in cl according to options power_xxx"""
        if cl != []: # shortcut for frequent case
            # else: Expand powers in cl according to options power_base and power_exp:
            if expand_power_base:
                cle = []      # apply _eval.expand_power_base() to all c-factors that are Pow's and pass other factors
                for c in cl:  # If result is a Mul flatten it into the result list cle
                    cle.extend(c.args) if (c.is_Pow and (c := c._eval_expand_power_base()).is_Mul) else cle.append(c)
                cl = cle      # resulting list of c-factors, Muls flattened in
            if expand_power_exp:
                cl = [(c._eval_expand_power_exp()  if c.is_Pow else c) for c in cl] # no hints applicable
                # no need to flatten Muls in cl as cl may contain Muls
        return cl

    def qapply_expand_mul(cl:List) -> Expr:
        """Expands list of commutative factors by key option 'mul' of expand(), i.e.
        does a deep distribution of factors over summands, and returns result"""
        return mc._eval_expand_mul() if (mc := Mul(*cl)).is_Mul else mc

    #-----------------
    # qapply_Mul2Atoms
    #-----------------

    # qapply_Mul2Atoms applies a set of application rules that Mul doesn't:
    # The main issue here is to deal with the plentitude of options that SymPy and the
    # quantum package provide how und where such rules for
    # 'lhs times rhs' may be defined and to predict the result that would have:
    # 1. lhs*rhs invokes lhs.__mul__ and then rhs.__rmul__  if defined in class or parents
    # 2. Mul(lhs, rhs) may behave differently from lhs*rhs (e.g. Bra(1)*Ket(1) vs. Mul(Bra(1),Ket(1)))
    # 3. lhs._apply_operator_<type(rhs)>(rhs), rhs._apply_from_right_to_<type(lhs)>(lhs) if defined
    # #  (but note: only for perfect fit of types; no searching up in parent classes)
    # 4. lhs.__pow__(rhs), lhs._pow(rhs) if defined in class or parents
    # Moreover:
    # 5. some classes provide xxx_simp or .flatten methods or other
    #    special handling routines for mul and power, e.g. TensorProduct

    def qapply_Mul2Atoms(lhs:Expr, rhs:Expr) -> Expr:
        """Applies atomic lhs and rhs to each other"""

        # Special case: lhs or rhs is PowHold type:
        # PowHold objects come here only if we have declared them atomic, because we don't know how
        # to safely extract factors from them, or base*base doesn't interact.
        # But joining powers on same bases is always safe and shortens the overall expression tree.
        # As PowHold is our artificial type we cannot use the logic in Pow() here, so mimic it here.
        # Note that this type generic rules cannot be expressed using _apply_operator / _apply_from_right_to:
        if isinstance(lhs, PowHold):
            cl, ncl, al = lhs.split_c_nc_atom(to_the_right) # extract factors of base to the right
            if isinstance(rhs, PowHold):
                if lhs.base == rhs.base:
                    return lhs.add_exp(rhs.exp)    # Note PowHold does A**0->1, A**1->A
                elif (lhs.exp - 1).is_positive and (rhs.exp - 1).is_positive:  # both powers can extract one base factor
                    cr, ncr, ar = rhs.split_c_nc_atom(to_the_left)  # extract factors of base to the left
                    alar = qapply_Mul2([al, ar], [])
                    if no_interaction(al, ar, alar):
                        return Mul(lhs, rhs, evaluate = False) # return unevaluated
                    else:
                        return Mul(*(cl + cr + ncl + [alar] + ncr))
                else: # lhs or rhs have exponents <= 1, so we cannot extract factors
                    return Mul(lhs, rhs, evaluate = False) # return unevaluated
            elif lhs.base == rhs: # catch special case rhs == lhs.base == PowHold(lhs.base,1)
                return lhs.add_exp(1)
            elif (lhs.exp - 1).is_positive: # (PowHold with 1 < exp) * anything
                alrhs = qapply_Mul2([al, rhs], [])
                if no_interaction(al, rhs, alrhs):
                    pass # maybe further down rhs can act on lhs, e.g. if rhs==Id
                else:
                    return Mul(*(cl + ncl + [alrhs]))
        elif isinstance(rhs, PowHold):
            if lhs == rhs.base: # catch special case lhs == rhs.base = PowHold(rhs.base,1)
                return rhs.add_exp(1)
            elif (rhs.exp - 1).is_positive:
                cr, ncr, ar = rhs.split_c_nc_atom(to_the_left) # extract factors of base to the left
                lhsar = qapply_Mul2([lhs, ar], [])
                if no_interaction(lhs, ar, lhsar):
                    pass # maybe further down lhs can act on rhs, e.g. if lhs==Id
                else:
                    return Mul(*(cr + [lhsar] + ncr))

        # case TensorProduct * TensorProduct: no *-Operator, ie. no __mul__ defined
        if isinstance(lhs, TensorProduct) and isinstance(rhs, TensorProduct) and len(lhs.args) == len(rhs.args):
            # Do the *-operator for TensorProduct*TensorProduct. If there are
            # c-factors in the tensor factors these will be pulled up front automatically, making it Mul
            result = TensorProduct(*[ cast(Expr, lf) * cast(Expr, rf) for (lf, rf) in zip(lhs.args, rhs.args)])
            #return qapply_Mul2([result], []) #qapply'ies to the tensor product factors
            return result # the tensor factors will be qapply'ed in slit_c_nc_atom_factor

        # add more cases type1 * type2 go here if they require recursive calls to qapply()
        # ...

        # Now try to apply the _apply_operator and _apply_from_right_to methods, so
        # these take precedence over the built-in methods below.
        result = cast(Any, None)
        try:
            result = cast(Any, lhs)._apply_operator(rhs, **options)
        except (NotImplementedError, AttributeError):
            try:
                result = cast(Any, rhs)._apply_from_right_to(lhs, **options)
            except (NotImplementedError, AttributeError):
                pass
        if result is not None:
            return result

        # If no succes so far: Try the following built-in rules.
        # Unfortunately the dispatch function behind _apply_operator/_apply_from_right_to
        # accepts perfect matches of types only and doesn't search up in class parents, so
        # we put those general rules here:
        if op_join:
            if isinstance(lhs, OuterProduct): # provide generic OuterProduct._apply_operator_*(rhs)
                res = qapply_Mul2Atoms(lhs.bra, rhs)
                if no_interaction(lhs.bra, rhs, res):
                    return Mul(lhs, rhs, evaluate = False) # return unevaluated
                else:
                    return lhs.ket * res
            if isinstance(rhs, OuterProduct): # provide generic OuterProduct._apply_from_right_to_*(lhs)
                res = qapply_Mul2Atoms(lhs, rhs.ket)
                if no_interaction(lhs, rhs.ket, res):
                    return Mul(lhs, rhs, evaluate = False) # return unevaluated
                else:
                    return res * rhs.bra

        # Still no result, so return lhs*rhs using *-operator, which
        # may be overloaded, so e.g. bra*ket->InnerProduct, ket*bra->OuterProduct etc.
        return lhs*rhs

    """
    Remark: qapply and powers Pow(base, power)
    1) e as a Pow object is already result of the constructor Pow.__new__(base, power) that
       contains at lot of logic to compute this power, including calling base._eval_power(power).
       Pow also has pulled (recognizable) commutative factors from its base into a separate Pow object.
       Some classes like OuterProduct and TensorProduct lack adequate methods .__pow__, ._pow
       or ._eval_power to apply powers that are required for qapply to work. I decided to add these
       here instead of in their classes. So if someone adds them to the classes, remove them here.
    2) qapply will try some 'tricks' to simplify certain powers that Pow doesn't use, nor 
       Pow._eval_expand_power_xxx, e.g. on idempotent operators like dyads.
    3) But the general approach is to unroll positive powers of an operator, i.e. apply it one by one.
       In order to make this more efficient, the base is evaluated by qapply only once and then cached
       (as a PowHold object). Nevertheless that may be costly.
    4) But: If higher powers of an operator op occur, the **best approach** is to make sure
       that Pow(op, exp) evaluates this, so the user should provide appropriate methods
       in the class definition of the operator op.
    """
    #--------------------
    # split_pow_c_nc_atom
    #--------------------

    def split_pow_c_nc_atom(e:Pow, split_left:bool) -> Tuple[List, List, Expr]:
        """Prepare expression Pow(base, exp) for further application; return (c-factor list if known, nc-list, atomR)"""

        # TensorProduct doesn't provide .__pow__ or ._pow. We don't use tensor_product_simp_Pow(e).
        if isinstance(e.base, TensorProduct):
            # Note ** will pull out commutative factors and return Mul(Pow(c-factors, exp), Pow(nc-part, exp))
            # And TensorProduct will also pull out all commutative factors from its elements, returning Mul in that case
            return split_c_nc_atom(
                TensorProduct(*[qapply_Mul2([b**e.exp], []) for b in e.base.args]), split_left)

        # try to turn a negative exponent positive, if base has an explicit inverse, so we may unroll the power
        # Pow() doesn't try this on creation of a Pow object. Requires exp to be integer to be always correct mathematically.
        if e.exp.is_negative and e.exp.is_Integer and hasattr(e.base, 'inv'): # < operator may be undecidable
            baseinv = e.base.inv()
            if not (baseinv.is_Pow and baseinv.exp.is_negative): # explicit inverse, not e.base**-1 etc.
                e = baseinv ** (-e.exp)    # rewrite e with a positive exponent (but may return anything, even neg. exp!)
                if not e.is_Pow: # if the power e was transformed into sth else
                    return split_c_nc_atom(e, split_left) # return it
                # else continue with the new e with whatever exponent


        # dispatch method depending on type of e.base, return or fall through with e as Power
        # ...

        # if no methods has given a result so far, at least recurse into e.base
        base = qapply_Mul2([e.base], [])
        # let Pow do all it can do using class methods ._pow, ._eval_power etc.
        e = cast(Pow, base ** e.exp)  # may return any Expr; cast() is for MyPy
        if e.is_Pow: # if it's still a Pow, try some
            # check for commutative elements in base
            base = e.base
            cl, ncl, ar = split_c_nc_atom(base, to_the_right) # cl*ncl*ar == base
            if ar is EmptyMul: # implies ncl==[] and base is all commutative
                return [e], [], EmptyMul
            if split_left:
                raise NotImplementedError('Splitting Powers to the left should not occur and has not been implemented.')
            if cl != []:       # strip off commutative factors from base
                base = qapply_Mul2(ncl + [ar], [])  # ncl*ar == base. Note that ar may by type Add.

            # Split exponent in part expi and "frational" part: e.exp = expi + expf. expi is not necessarily integer!
            expi = (floor(e.exp) if (e.exp.is_positive) else (floor(e.exp) - 1 if (e.exp.is_negative) else e.exp))
            expf = e.exp - expi # fractional part, 0 if we couldn't identify a fractional part

            # if expi can be be determined to be integer exponent >= 2, we try harder:
            if expi.is_integer and 2 <= expi:  # Note: use .is_integer; .is_Integer checks for type Integer
                if ncl != []: # base has two or more factors
                    # if base=a*b*c and c*a is commutative, then base**expi=(c*a)**(expi-1) * a * b**expi * c.
                    # We don't check all possible factors, just the simplest case so that
                    # cases like (|k>*<b|)**exp or (A*B*A.inv)**exp are simplified
                    c2, ncr2, al = split_c_nc_atom(Mul(*ncl), to_the_left) # al*ncr2*ar == base
                    if c2 != []: # should not happen, since we split ncl. But done here for safety
                        base = qapply_Mul2(ncr2 + [al], [])
                    p = PowHold(base, e.exp, True, ([], ncr2 + [ar], al), ([], [al] + ncr2, ar))
                    pfl = [] if expf == 0 else [ p.copy_exp_atomic(expf, True) ] # fractional power
                    ia = qapply_Mul2([ar, al], [])  # note ar, al may by type Add
                    if no_interaction(ar, al, ia):  # factors didn't interact, so it makes
                        # no sense for the power pi to unroll, as factors will just line up, so keep atomic
                        # Case 1: expi>=2, component base, base*base no interaction
                        return [c ** e.exp for c in cl+c2], [], p
                    else: # factors did interact: check whether ia is commutative
                        if (getattr(ia, 'is_commutative', False) or isinstance(ia, Number)):
                            # (No need to check instance(ia, IdentityOperator), as ia would be S.One instead)
                            # base = a*b*c and c*a is commutative, so: base**expi = (c*a)**expi * a * b**expi * c
                            # Case 2: expi>=2, component base, base*base interacts, and commutative
                            return [c ** e.exp for c in cl+c2] + [ia ** (expi - 1)],\
                                    pfl + [al, Pow(Mul(*ncr2), expi)], ar
                        else: # it makes at least sense to permit the power to unroll itself
                            #Case 3: expi>=2, component base, base*base interacts, but non-c
                            return [c ** e.exp for c in cl+c2], [p.copy_exp_atomic(e.exp - 1, False)] + ncl, ar
                else: # ncl == [], base is atomic, expi is Integer >= 2
                    p   = PowHold(base, e.exp, atomic = True , b_sp_l = ([], [], ar), b_sp_r = ([], [], ar))
                    pfl = [] if expf == 0 else [ p.copy_exp_atomic(expf, True) ] # fractional power
                    ar2 = qapply_Mul2([ar, ar], [])  # does the base interact with itself?
                    if no_interaction(ar, ar, ar2):           # factors didn't interact, so it makes
                        # no sense for the power to unroll, as factors will just line up, so make atomic
                        # Case 4: expi>=2, base atomic and base*base no interaction
                        return ([c ** e.exp for c in cl], [], p)
                    else: # Check for important simplifications if base is a dyad, "idempotent" or
                        # "involutoric" operator and permit the power to unroll:
                        ar2c, ar2ncl, ar2ar = split_c_nc_atom(ar2, to_the_right)
                        if ar2ar == EmptyMul: # ar2 is commutative, and ncl==[], so base = cl*ar.
                            # Case 5a: expi >= 2, base atomic + "involutoric"
                            # simplify: base**e.exp = base**expf * cl**e.expi * ar**e.expi-1 * ar
                            c_list = [c ** e.expi for c in cl] + [c ** (expi // 2) for c in ar2c]
                            return (c_list, pfl, ar) if (expi % 2 == 1) else (c_list, [], Mul(*pfl))
                        if ar2ar == ar and ar2ncl == []: # Case 5b: expi >= 2, base atomic + "idempotent"
                            # kind of idempotency: a*a = c*a -> a**expi = c**(expi-1) * a
                            return [c ** e.exp for c in cl] + [c ** (expi - 1) for c in ar2c], pfl, ar
                        # else: Case 6: base atomic but not "idempotent"
                        # return base**(exp-2) * ar * ar = base**(exp-2)*ar2c*ar2ncl*ar2ar
                        return [c ** e.exp for c in cl] + ar2c, [p.copy_exp_atomic(e.exp - 2, False)] + ar2ncl, ar2ar
            else: # Case 7: exponent is non-integer or < 2
                return [c ** e.exp for c in cl], [], \
                    PowHold(base, e.exp, atomic = True, b_sp_l=tuple(), b_sp_r = ([], ncl, ar))
        else: # Case 8: e is no longer a power
            return split_c_nc_atom(e, split_left)


    #-------------------
    # Factor Helper functions for qaaply_Mul2 - variables ip_doit, dagger, options etc are in scope
    #-------------------
    #
    # Code relies on the sympy convention that all commutative factors have been pulled out
    # from objects like OuterProduct, TensorProduct, InnerProduct, Pow etc by their constructors,
    # so if an an expression e is an instance of those objects e does not contain such factors

    def split_c_nc_atom(e:Expr, split_left:bool) -> Tuple[List, List, Expr]:
        """returns (c_factors list of e, e/(c_factors*nc_atomic_factor) list, nc_atomic_factor)"""
        if  e == EmptyMul or e == S.Zero:  # speed up for special cases
            return [], [], e

        elif e.is_Mul: # by definition of Mul has at least 2 factors in it
            # Currently .is_commutative ensures 'Mul can handle it'. Should non-scalar factors
            # with .is_commutative=True arrive, a filter like hasattr('_apply_operator') is required
            c1, nc = e.args_cnc() # depends on .is_commutative
            if ip_doit: # evaluate InnerProducts if ip_doit==True. As .is_commutative they are in c.
                c1 = [(f.doit() if isinstance(f, InnerProduct) else f) for f in c1]
            #extract atomic factor from nc:
            if  nc == []:   # if e was all commutative
                return c1, nc, EmptyMul
            else: # nc is list with either 1 or more elements
                # pick the left- resp. right most factor f from nc
                nc1, f = (nc[1:], nc[0]) if split_left else (nc[:-1], nc[-1])
                # break up f further. c2 has ip_doit done.
                c2, nc2, a = split_c_nc_atom(f, split_left)
                # Note that x.extend(y) returns None, so can't use c1.extend(c2) etc.
                return (c1+c2, (nc2+nc1 if split_left else nc1+nc2), a)

        elif e.is_Pow: # by definition represents at least 2 factors, except A**-1:
            # e.exp == 0 or e.exp == 1 will not occur since Pow(A,0)->1, Pow(A,1)->A
            # Do simplification of Powers, incl. ip_doit, and wraps Powers in PowHold class
            return split_pow_c_nc_atom(e, split_left)
        elif isinstance(e, PowHold):
            if getattr(e, "atomic"): # so MyPy won't complain about e.atomic
                return [], [], e  # PowHold won't roll off factors in split
            else: # PowHold rolls itself off providing factors
                return e.split_c_nc_atom(split_left)

        elif isinstance(e, InnerProduct):
            return [e.doit() if ip_doit else e], [], EmptyMul # is commutative in any case

        elif isinstance(e, OuterProduct):
            if op_join:
                return ([], [], e)  # consider OuterProduct atomic
            else:                   # break it up into ket * bra
                return ([], [e.bra], e.ket) if split_left else ([], [e.ket], e.bra)

        elif isinstance(e, Density):    # For a Density need to call qapply on its state vectors
            new_args = [(qapply_Mul2([cast(Expr, state)],[]), cast(Expr, prob)) for (state, prob) in cast(Tuple, e.args)]
            res = cast(Expr, Density(*new_args)) # might return Mul or Density (for MyPy)
            if isinstance(res, Density):
                return ([], [], res)  # Density is considered atomic
            else: # e.g. Mul, OuterProduct, Add of OuterProduct etc..
                return split_c_nc_atom(res, split_left)

        elif isinstance(e, TensorProduct): # For a raw TensorProduct, call qapply on its args first
            res = cast(Expr, TensorProduct(*[qapply_Mul2([t], []) for t in e.args])) # may return Mul or TP (for MyPy)
            if expand_tensorproduct and isinstance(res, TensorProduct):
                res = res._eval_expand_tensorproduct() # no hints for expansion
            if isinstance(res, TensorProduct):    
                return ([], [], res)  # Tensorproduct itself is atomic
            else: # e.g. has become Mul, scalar, ..
                return split_c_nc_atom(res, split_left)

        elif isinstance(e, (Commutator, AntiCommutator)): #requires .doit() first
            res = e.doit()            # may return any type including Add, S.Zero, unmodified,..
            if type(res) is type(e):
                return ([], [], res)  # recursion would be infinite, so nothing we can do with it
            else:
                return split_c_nc_atom(res, split_left)

        else: # generic rules that are not type dependendent:
            if (getattr(e, 'is_commutative', False) or isinstance(e, Number)):
                # remaining e is commutative and needs no .doit() etc.
                return ([e], [], EmptyMul)
            else: # e is nc, doesn't contain commuting factors, needs no .doit() etc., e.g. Ket, Bra, Add
                return ([], [], e)


    # Minimalistic wrapper for Pow to hold algorithm from iterated inspection of base in Pow(base,exp)
    # Derived from Pow.__new__:
    class PowHold(Expr): # formally make it an Expr, so it looks more SymPy
        # see definition of class Pow: this permits modification of obj.is_commutative
        __slots__ = ('is_commutative', '_b_sp_l', '_b_sp_r', 'atomic')

        def __new__(cls, base:Expr, exp:Expr, atomic:bool, b_sp_l:Tuple, b_sp_r:Tuple):
            if exp is S.Zero:  return S.One    # by sympy convention, see Pow.
            elif exp is S.One: return base
            obj = Expr.__new__(cls, base, exp) # sets up self.func, self.args
            obj.is_commutative = (base.is_commutative and exp.is_commutative) # see Pow.__new__
            if base.is_commutative:
                obj._b_sp_l = obj._b_sp_r = ([base], [], S.One)
            else:
                obj._b_sp_l = b_sp_l # base split to left:  (c, ncr, al) with c*al*ncr = base
                obj._b_sp_r = b_sp_r # base split to right: (c, ncl, ar) with c*ncl*ar = base
            obj.atomic = atomic  # true if PowHold behaves like atomic
            return obj

        def adjoint(self): # required to make Dagger(PowHold) work
            def dagger_sp(b_sp:tuple) -> tuple:
                return (b_sp if (b_sp == tuple()) else
                        ([Dagger(f) for f in b_sp[0] ], [Dagger(f) for f in reversed(b_sp[1])], Dagger(b_sp[2])))
            (b, e) = (b.base, b.exp*self.exp) if ((b := Dagger(self.base)).is_Pow) else (b, self.exp) # e.g. if Dagger(Unitary) returns Unitary**-1            
            return PowHold(b, e, self.atomic, dagger_sp(self._b_sp_l), dagger_sp(self._b_sp_r))

        def add_exp(self, exp2): # return PowHold with exponent incremented by exp2
            return PowHold(self.base, self.exp + exp2, self.atomic, self._b_sp_l, self._b_sp_r)
        def copy_exp_atomic(self, exp, atomic):
            return PowHold(self.base,      exp,             atomic, self._b_sp_l, self._b_sp_r)
        def split_c_nc_atom(self, split_to_left):
            # split up by split_left; ignores .atomic, hands factor out as long as possible
            if (self.exp - 1).is_nonnegative:
                if split_to_left:    # extract base to the left:
                    if self._b_sp_l == tuple():
                        self._b_sp_l = split_c_nc_atom(self.base, to_the_left)
                    res = self._b_sp_l[0], self._b_sp_l[1] + [self.add_exp(-1)], self._b_sp_l[2]
                    return res
                else:          # extract base to the right:
                    if self._b_sp_r == tuple():
                        self._b_sp_r = split_c_nc_atom(self.base, to_the_right)
                    res = self._b_sp_r[0], [self.add_exp(-1)] + self._b_sp_r[1], self._b_sp_r[2]
                    return res
            else:
                return [], [], self.copy_exp_atomic(self.exp, True)

        @property
        def base(self): return self._args[0]
        @property
        def exp(self) : return self._args[1]

    #-----------------------------------
    # qapply(e:Expr, **options) -> Expr:
    #-----------------------------------

    # if e is no sympy expression, return it unmodified (e.g. Numbers)
    if not isinstance(e, Expr):
        return e

    # set options for qapply_* and helper functions:
    dagger               = options.get('dagger', False)       # if a*b doesn't compute, try Dagger(Dagger(b)*Dagger(a))
    ip_doit              = options.get('ip_doit', True)       # evaluate ip.doit() on all InnerProduct objects
    op_join              = options.get('op_join', True)       # try to join Ket*Bra to OuterProduct, avoid breaking OuterProduct
    expand_mul           = options.get('mul', True)           # distribute commutative factors and apply hint 'mul' of expand() 
    expand_power_base    = options.get('power_base', expand_mul) # for commutative factors: apply hint 'power_base' of expand() 
    expand_power_exp     = options.get('power_exp', expand_mul)  # for commutative factors: apply hint 'power_exp' of expand()
    expand_tensorproduct = options.get('tensorproduct', True) # expand sums in tensor product factors (= hint tensorproduct of expand())

    # some constants for convenience
    EmptyMul = Mul() # is Integer(1), is S.One, but more intuitive in this context
    to_the_left  = True  # option for split_c_nc_atom: split atomic factor to the left
    to_the_right = False # option for split_c_nc_atom: split atomic factor to the right

    # call the workhorse, with all variables from this namespace (dagger, ip_doit, ...) in scope
    res = qapply_Mul2([e], [])

    # clean up the result, if necessary
    if hasattr(res, 'replace'):
        res = res.replace(PowHold, Pow)

    """ debugging aid: show various representations and expansion of result
    from sympy.printing import srepr
    try:    # in case e contains objects (e.g. functions) that crash print()
        ep = str(e)
    except:
        ep = f"{str(type(e))}({str(e.args)})" # strange enough, printing e.args works
    print(f"    1.qapply({ep},{options}) -> {srepr(res)}")
    print(f"    2.qapply({ep},{options}) -> {res}")
    if hasattr(res, 'expand'):
        print(f"    3.qapply({ep},{options}).expand(TP=True) = {res.expand(tensorproduct=True)}")
    #"""
    return res
