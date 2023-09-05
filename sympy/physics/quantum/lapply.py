"""
``lapply`` (from "linear apply") is a generic evaluation engine to apply linear
operators, powers and factors and distribute them over summands, with focus on
non-commutative operators.

Its primary purpose is evaluate larger expression trees of SymPy objects whose
multiplication has not been defined with :obj:`~.Mul()` and without having to
overload the ``*`` operator. To this end ``lapply`` may be extended by
``multipledispatch`` handlers for multiplication or to break up compound
objects into elementary objects.

``lapply`` is also useful to compute expressions that grow too large without
intermediate simplification and would otherwise need to be manually split up
and computed piecewise. The amount of expansion of terms may be tuned by
options. ``lapply`` also contains simplifications for powers with symbolic
exponents and involutoric or idempotent factors.

Handlers for SymPy matrix objects are included so ``lapply`` is available to
evaluate matrix expressions that contain ``MatrixSymbols``.
"""
#
# lapply evolved as the generic kernel required to apply the objects of the
# quantum package (see sympy.physics.quantum.qapply).
# Handlers for matrices are included in the generic lapply module as a
# convenience, also as matrices are loaded with "from sympy import *"
#
# Design Principles:
#
# 1. The objects and notation of the quantum package (Ket, Bra, ...) are just
#   special cases of general objects (vector, linear operator, ...) and
#   the notation of standard Linear Algebra. The same is true for the
#   matrices package, quaternion package etc. Their addition is derived from
#   Add, and their multiplication from Mul, their fields of commutative
#   scalars being the commutative scalars of SymPy. However the rules for
#   multiplication/application of this objects are not all known to Mul (
#   and perhaps should not for performance reasons), and the coding differs
#   per package.
#   So lapply is a generic evaluation engine for associative and distributive
#   law trying the rules for multiplication provided by handlers, defaulting
#   to * resp. Mul. Special rules have been provided for handling of powers
#   with idempotent or involutoric factors in the base.
#
# Object type specifics have completely been encapsuled in @dispatch handlers:
# a. The algorithm that walks the expression tree and determines to pair factors
#   a,b,c,d as (a*b)*(c*d) or a*((b*c)*d) etc. is completely object type agnostic
#   (lapply_Mul2). Walking of the tree is dependant on the success of pairing:
#   I decided to implement (a*(b*(c*d))) association from the right in general
#   but depending on whether factors interacted to new factors and re-tries to
#   the right side in that case.
# b. Recursing into inner structure of an object (lapply1_type handlers)
# c. Splitting up composite objects in commutative (c) and non-commutative
#   (nc) factors (c_nc_ncef handlers)
# d. All rules how to multiply two object types (lapply_mul2elems with
#   lapply_mul2types_x handlers)
#
# Moreover:
# 5. expand() is not used indiscriminately: Only do an expansion of Add if
#   required to apply operators, so keeping expansion to the bare minimum. As
#   default, distribute commutative factors over sums as this allows terms to
#   cancel out thus simplifying the expression. In order to simplify commuta-
#   tive factors lapply() supports the three major hints of expand() "mul",
#   "power_base" and "power_exp" (see options mul, power_base and power_exp).
# 6. apply recursion along the structure of the expression but don't use
#   recursion for iterations: this avoids hitting the Python recursion limit
#   in very complex expressions or powers with exponents > 1200 or so.
#7. No need to try to maintain subclass types of objects and expressions, as
#   we simply cannot know whether a product is still in a subclass of Mul etc.
#   So use the top types.

from typing import cast, Any, List, Tuple, Callable, Union  # For Python 3.8
from enum import Enum

from sympy.core.expr import Expr
from sympy.core.add import Add
from sympy.core.mul import Mul
from sympy.core.power import Pow
from sympy.core.singleton import S
from sympy.core.numbers import Number
from sympy.functions.elementary.integers import floor
from sympy.multipledispatch.dispatcher import Dispatcher

# lapply is object agnostic. All object specific information has been
# isolated in handlers that are dispatched by object type.
# Imports required to define object types for handlers are coded
# just in front of the handlers.

__all__ = [
    'lapply'
]

#--------+---------+---------+---------+---------+---------+---------+---------
# lapply function header and doc string
#
# Options taken from expand() have same names (mul, power_xxx, tensorproduct)
#------------------------------------------------------------------------------

def lapply(e:Expr, mul=True, Add=Add,
                   power_base=None, power_exp=False, nested_exp=None,
                   apply_exp=False, cache_for_pow=None, **options) -> Expr:
    """
    Generic evaluation engine to apply linear operators, powers and factors
    and distribute them over summands with focus on non-commutative operators.
    Multiplication of arbitrary SymPy objects can be defined by
    ``multipledispatch`` handlers (see :obj:`sympy.physics.quantum.qapply()`
    for an example).

    ``lapply`` incorporates simplification rules for powers with symbolic
    exponents and involutoric or idempotent factors in their base, so it can
    help to simplify terms in a way like combinations of
    :obj:`sympy.core.function.expand()`, :obj:`sympy.simplify.simplify()`,
    :obj:`sympy.simplify.powsimp()` and ``.doit()``.

    The level of expansion may be tuned by options especially for commutative
    factors. The handling of matrices is included for matrices that are SymPy
    expressions, i.e. derived from :obj:`~.MatrixExpr` and immutable, and is
    best used with terms that contain :obj:`~.MatrixSymbol` (see examples).

    Parameters
    ==========

    e : Expr
        The expression containing operators, states resp. factors, summands
        and powers. The expression tree will be walked and operators resp.
        factors and powers applied and distributed over summands.

    mul : Boolean, optional
        If True (default) distribute commutative factors over sums. Corresponds
        to option ``mul`` of :obj:`sympy.core.function.expand()`. ``mul=False``
        will do the minimal application only, so use only when expansion is slow.

    power_base : Boolean, optional
        Refers to commutative factors. Same option as in
        :obj:`sympy.core.function.expand()`:. If True split powers of multiplied
        bases. Defaults to value of option ``mul``.

    power_exp : Boolean, optional
        Refers to commutative factors. Same option as in
        :obj:`sympy.core.function.expand()`:.If True expand addition in exponents
        into multiplied bases. Defaults to False.

    nested_exp : Boolean, optional
        If True, then if powers are nested and if their exponents may be
        multiplied, explicitly expand the product of the exponents. Defaults
        to value of option ``mul``.

    apply_exp : Boolean, optional
        If True apply ``lapply`` to exponents of a power before applying the
        exponent to the base. Default is False and to apply exponents as provided.

    cache_for_pow : dict, optional
        a dictionary that provides cache values to speed up the evaluation
        of powers. If provided it is updated in place, so it can be passed on
        to the next invocation of ``lapply`` as long as options and attributes
        of symbols remain unaltered.
        Defaults to None, in which case an internal cache is used.

    Add : Callable, optional
        This function is used to add summands. It must accept the option
        ``evaluate=True`` and return type :obj:`~.Expr`. If set to None, the
        individual .func argument of each addition expression will be used.
        Defaults to standard :obj:`~.Add()`.

    options : dict, optional
        A dict of keyword/value pairs that are passed to the handlers and the
        operator methods as operator and handler specific options.


    Returns
    =======

    e : Expr
        The original expression, but with the operators applied to operators
        and states and with moderate expansion done by distributing commutative
        factors over summands.


    Examples
    ========

        >>> from sympy.physics.quantum.lapply import lapply, c_nc_ncef #------------Adapt path to final location of lapply
        >>> from sympy import symbols, Mul, sqrt, S
        >>> from sympy.matrices import (ImmutableMatrix, MatrixSymbol, MatMul, Identity)
        >>> from sympy.algebras.quaternion import Quaternion

        >>> a, b, c, d, x = symbols("a b c d x", commutative=True)
        >>> M1 = ImmutableMatrix([[a, 2], [3, 4]])
        >>> A = MatrixSymbol("A", 2, 2)
        >>> B = MatrixSymbol("B", 2, 2)
        >>> o = symbols("o", integer=True, nonnegative=True, even=True)

        ``lapply`` passes explicit matrices on to the matrix package, but
        may achieve strong simplification of symbolic matrix terms:

        >>> Z = M1.T * (A * M1.I + B * M1.I).T
        >>> lapply(Z)
        A.T + B.T
        >>> p = Mul(*[(k*Identity(2) + Z - A.T - B.T) for k in range(1,3)])
        >>> lapply(p * A)
        2*A
        >>> (p * A).expand().simplify().doit() == 2*A
        False
        >>> lapply((p * A).expand().simplify().doit())
        2*A

        Also with powers of matrices with ``MatrixSymbols``:

        >>> M3 = ImmutableMatrix([[0, a], [b, 0]])
        >>> Z = ((M3 * A * M3)**(o + 10000)).simplify()
        >>> (Z.base, Z.exp) == (M3 * A * M3, o + 10000)
        True
        >>> lapply(Z) == MatMul(a**(o+9999)*b**(o+9999), M3, A**(o+10000),M3)
        True

        >>> Idb = ImmutableMatrix([[b, 0], [0, b]])
        >>> Z = ((M1 * A * M1.I * Idb)**10000).simplify()
        >>> (Z.base, Z.exp) == (M1 * A * (b * M1.I.simplify()), 10000)
        True
        >>> lapply(Z) == MatMul(b**9999, M1, A**10000, (b*M1.I.simplify()))
        True

        >>> M2 = ImmutableMatrix([[0, sqrt(a*(a+1))], [sqrt(a*(a-1)), 0]])
        >>> F = M2 * Idb * M1
        >>> Z = ((M1 * A * M2 * Idb)**4).simplify()
        >>> (Z.base, Z.exp) == (M1 * A * (b * M2.simplify()), 4)
        True
        >>> lapply(Z) == M1 * A*F* A*F* A*F* A * M2*Idb
        True

        Another example: a matrix annuls it characteristic polynomial:

        >>> Mq = ImmutableMatrix([[a, -b, -c, -d], [b, a, -d, c],
        ...                       [c,  d,  a, -b], [d, -c, b, a]])
        >>> cp = Mq.charpoly().eval(x)
        >>> lapply(cp.subs(x, Mq))
        0

        As ``Add`` doesn't handle Quaternions, this example shows how
        to provide an proprietary addition function to ``lapply`` (note: don't
        use complex numbers with Quaternions as SymPy and so ``lapply`` consider
        complex numbers as commutative which is not correct with Quaternions):

        >>> q = Quaternion(a, b, c, d)
        >>> quatAdd = (lambda *summands, evaluate=True: sum(summands).expand())
        >>> lapply(cp.subs(x, q), Add=quatAdd)
        0 + 0*i + 0*j + 0*k

        If ``lapply`` shall handle real Quaternions as commutative scalars,
        add a handler for Quaternion objects to the ``c_nc_ncef`` dispatcher:

        >>> c_nc_ncef.add((Quaternion,), (lambda q, to_L, **options:
        ...   ([q.scalar_part()], [], S.One) if q.vector_part().is_zero_quaternion()
        ...              else ([],[], q) ))
        >>> lapply(q * q.conjugate())
        a**2 + b**2 + c**2 + d**2

        See :obj:`sympy.physics.quantum.qapply()` as example how to define handlers
        for a collection of objects.
    """
    # if e is no SymPy Expr, return it unmodified (e.g. Numbers).
    # It is useless to call lapply with Basic or even less object types, as
    # the purpose of lapply requires that at least the arithmetic operation
    # "apply" is defined on all objects in e in some form. Since in SymPy
    # "Everything that requires arithmetic operations to be defined should
    # subclass Expr" (see definition of class Expr), and since this holds
    # for numbers and symbols too, no need to consider type Basic anywhere.
    if not isinstance(e, Expr) or e.is_Atom:
        return e

    # set options for lapply_* and helper functions:
    # distribute commutative factors and apply hint 'mul' of expand()
    options['mul'          ] = opt_mul = mul
    # for commutative factors: apply hint 'power_base' of expand()
    if power_base    is None: power_base    = opt_mul
    options['power_base'   ] = power_base
    # for commutative factors: apply hint 'power_exp' of expand()
    options['power_exp'    ] = power_exp
    # if powers are nested, expand_mul() the product of exponents
    if nested_exp     is None: nested_exp   = opt_mul
    options['nested_exp'   ] = nested_exp
    # apply lapply to exponents, otherwise exponents are used "as is"
    options['apply_exp'    ] = apply_exp
    # Add function, or .func per addition expression
    options['Add'          ] = Add  # no check done

    # The cache for Pow_with_Cache() is just a dictionary that is
    # either passed as option "cache_for_pow" or provided here, so
    # it is local to this particular invocation of lapply. This makes
    # it threadsafe and also **bound to the values of the options**.
    # If a dict is provided, entries will be added and the caller may
    # thus re-use it in the next invocation.
    options["cache_for_pow"] = cache_for_pow \
                 if isinstance(cache_for_pow, dict) else {}

    # call the workhorse
    res = lapply_mul2([e], [], **options)
    # clean up if any, and return
    return res                              # End of lapply


# some modul local constants for convenience
EmptyMul = Mul()     # is Integer(1), is S.One, but more intuitive here
to_the_left  = True  # option for qrec_x: split elem factor to the left
to_the_right = False # option for qrec_x: split elem factor to the right


# Internal Remarks on object type constructors and Mul
# -----------------------------------------------------
#
# 1.) The code relies on the auto-simplifications that Mul and the composite
# objects constructors like Pow, Add, InnerProduct, OuterProduct, TensorProduct
# etc do: Extraction of commutative factors, contracting of same factors to
# Powers etc, so these constructors may return their type, Mul or other types.
#
# 2.) Mul(a,b) != a * b:
# Mul(a, b, c) processes the arguments a, b, c (e.g. factoring out scalars,
# grouping same factors to powers) but actually multiplies the
# arguments only in simple cases or if a postprocessor is defined.
# On the other hand a*b is executed using *-operator __mul__ for type(a),
# type(b): E.g. for Kets and Bras, __mul__ is overloaded to form InnerProducts
# resp. OuterProducts types, and e.g. IdentityOperator()*Ket(1)-> Ket(1)
# (Examples: Mul(Bra(1), Ket(1)) returns unchanged. Bra(1)*Mul(Ket(1), Bra(1))
# returns Mul(Bra(1), Ket(1), Bra(1)), while Bra(1)*Ket(1)*Bra(1) returns Mul(
# InnerProduct(Bra(1),Ket(1)), Bra(1)). So while the return type of Mul(a,b)
# will be Mul (or Zero) the return type of a*b completely depends on the input.
#
# 3.) Performance of Mul: refer to the remarks in code of Mul.flatten:
# - avoid a*b*c if you don't need the multiplication and use Mul(a,b,c) instead
# - instead of Mul(a,b)*Mul(c) use Mul(a,b,c)
# - results of Mul(a,b,c) are cached

# Remarks on lapply vs. Mul and expand, and "elementary" factors:
# ---------------------------------------------------------------
#
# I will use the verbs 'to apply' and 'to multiply' equivalently as in this
# context of linear operators they mean the same and are both denoted by *.
#
# So the purpose of lapply is to multiply factors that Mul doesn't know
# how to multiply, and leave those that Mul knows to Mul. This necessarily
# includes the distribution of factors over Add summands as far as
# required which forces lapply to partly appear similar to expand().
# Since lapply also has to handle powers (which come in naturally as
# repetitions of operators), nested powers (which come in as repetitions of
# repetitions of operators) and fractional powers (which come in through
# fractional powers of operators), this is where I see the difference to
# expand(): lapply will do only expansions that are required: E.g. it won't
# expand arguments of functions, won't lapply exponents of powers (except
# if explicitly optioned), and focusses on getting factors multiplied.
# On the other hand expand's primary task is to form terms and products, its
# hints are on how to expand a given type, but not on how to multiply
# two factors. lapply is trying hard to apply factors to each other, handle
# powers and get a result, and expansion is a byproduct.
#
# Now, to draw the border between lapply and Mul, we have to classify all
# factors in two ways:
#
# 1. Does Mul know how to multiply this type with other types, or doesn't it?
#
#   For efficiency, lapply would try to only handle those types that Mul
#   doesn't know how to multiply, and will leave all others to Mul. This
#   also helps to avoid contradicting results of Mul and lapply.
#   Of course passing factors to Mul makes only sense if there are more
#   than one to pass. So commutative factors are best as they can be
#   grouped together up front. So it fits just nicely that many relevant
#   mathematical objects are non-commutative and require multiplication
#   methodes that Mul probably doesn't know.
#   So I choose to leave commutative objects to Mul assuming that Mul
#   knows how to handle them. So lapply pulls out commutative factors. All
#   it does is to apply lapply1_type handlers once (to process rare cases
#   like InnerProduct) before passing them on.
#
#   Scalars would have been a even more simpler group to leave to Mul.
#   Unfortunately .is_scalar doesn't exist, and the closest approximation
#   'is_complex' can not be used: on most symbols, is_complex is undefined.
#   The only predicate that SymPy guarantees is 'is_commutative'!
#
#   So we use 'is_commutative' to pass factors on to Mul. This might raise
#   issues with commutative operators like Identity or IdentityOperator
#   that are unknown to Mul. So lapply adapts the usance that commutative
#   operators should be reduced to the scalar they are. And luckily
#   Identity.is_commutative==False and IdentityOperator.is_commutative==False.
#   But when this quirks should vanish action is required.
#
# 2. Is the factor a compound type (e.g. OuterProduct, Commutator) that is
#   actually an expression that should be evaluated and broken up into
#   its components to be multiplied one by one, or is it 'elementary'?
#   So lapply implements rules (Pow, Mul) and the handler family c_nc_ncef
#   to break up compound types (e.g. (Anti-)Commutator) or which to break as
#   an option (e.g. OuterProduct). Those that will not be broken up are called
#   "elementary" (note that this is similar in meaning but not identical
#   to the sympy definition of Atom or is_Atom!).
#
#   Note that "elementary" refers to the general algorithm (lapply_mul2)
#   that handles walking the expression tree and the association of factors.
#   lapply_mul2elems may still inspect the inner contents of factors to do
#   the multiplication (see e.g. handler for OuterProduct).


def lapply_no_act(lhs:Expr, rhs:Expr, res:Expr) -> bool:
    """ Assumes lhs, rhs elementary and nc and res=lapply_mul2([lhs, rhs], []).
    Returns True if lhs, rhs didn't interact so lhs*rhs just gave
    res=Mul(lhs,rhs) or Pow(lhs, 2)."""
    return ((res.is_Mul and res.args == (lhs, rhs)) or
            (lhs == rhs and res.is_Pow and
             res.base == (lhs.base if lhs.is_Pow else lhs)) # type: ignore[attr-defined]
            )
    # Remark: This heuristics accepts A**3 * A**4 = A**7 as interaction,
    # but rejects A**3 * A**3 = A**6. This is since Mul(X,X)=Pow(X,2) is
    # rejected as just a shorthand in any case. Combining A**3 * A**7
    # takes more mathematics, so this is accepted, and consequently
    # additions of power exponents is implemented in lapply_mul2elems().
    # Along the same lines lhs=(AxB), rhs=(CxD) and res=(A*C)x(B*D)
    # returns False ("interaction"), otherwise lapply wouldn't multiply
    # tensor products per factor, quite unexpectedly (same action as
    # tensor_product_simp_Mul())


# Focus of lapply_mul2 is on non-commutative operators. Handling of
# commutative factors (i.e. scalars!) is therefore delegated to Mul.
# However in order to simplify some expressions lapply needs to deeply
# expand and multiply out all commutative Add-factors to find whether
# they cancel out or at least become much simpler. Also expansion of
# powers via expand(power_base=True) or expand(power_exp=True) is enabled
# by default, as it is in expand(). On the other hand we do not want to
# call expand(deep=True), as this also aggressively expands function
# arguments like sqrt(j*(j+1)+m*(m+1)). So we use the parts of expand that
# have greatest effect: the hint function Mul._eval_expand_mul and
# expand() with options power_xxx, which are invoked in the helper functions
# lapply_expand_powers() and lapply_expand_mul():

# Helper functions to expand commutative factors (wrapping expand())
def lapply_expand_powers(cl:List, **options) -> List:
    """Expands commutative factors in cl according to options power_xxx"""
    if cl == []: return cl # shortcut for frequent case
    # else: Expand powers in cl according to options power_base/_exp:
    popts = { key: options[key] for key in ["power_base", "power_exp"] }
    if popts == {}: return cl # nothing to be done
    # set the additional options of expand, especially mul=False!
    popts.update({"deep": True, "mul":False, "log":False,
                                "multinomial":False, "basic":False })
    cle = [] # Though cl may contain Muls, we flatten the expansion..
    for c in cl:  # ..results, just for consistency. It's not required.
        c = c.expand(**popts)
        cle.extend(c.args) if c.is_Mul else cle.append(c)
    return cle

# This is a "less intrusive" version of expand_mul():
def lapply_expand_mul(e:Expr) -> Expr:
    """Like expand_mul(), does a deep distribution of factors over summands.
    But other than expand_mul() doesn't expand into arguments of functions."""
    # Reduced version of sympy.core.expr.expand() and Expr._expand_hint()
    # relies on Mul._eval_expand_mul() that expands into Add's only.
    chg = True
    def lapply_expand_chg(e:Expr) -> Expr:
        nonlocal chg
        if e.is_Mul:  # e._eval_expand_mul returns e if no change
            e1 = e._eval_expand_mul(deep=True) # type: ignore[attr-defined]
            chg |= (e1 is not e)
        elif e.is_Add:
            sargs = [lapply_expand_chg(ae) for ae in e.args] # type: ignore[arg-type]
            e1 = e.func(*sargs) if chg else e
        else: return e
        return e1
    while chg:
        chg, e = False, lapply_expand_chg(e)
    return e


#--------+---------+---------+---------+---------+---------+---------+---------
# lapply_mul2 / lapply_mul2c_nc
# -----------------------------------------------------------------------------
# Multiply all factors from fL to fR from the left
# - fL is typically derived by list(Mul.args), but may contain any
#      expression incl. Mul, Add and commutative arguments.
# - fR may not contain Add, Mul or factors Mul should multiply, i.e.
#      may only contain "elementary" factors (see above)
# - options MUST contain settings for all lapply option keywords!
# - Important: fL, fR are modified during execution!
#
# lapply_mul2_c_nc returns (comm. factor list, nc factors list), and
# lapply_mul2 the full Mul of that.

def lapply_mul2(fL:List, fR:List, **options) -> Expr:
    """Multiplies the left hand list of factors into the right hand list
    of factors. Returns Mul'ed result."""
    (cl, ncl) = lapply_mul2_c_nc(fL, fR, **options)
    if options["mul"]:      # expand(cl) and use * instead of Mul
        return lapply_expand_mul(Mul(*cl)) * Mul(*ncl)
    return Mul(*(cl + ncl))

def lapply_mul2_c_nc(fL:List, fR:List, **options) -> Tuple[List, List]:
    """Multiplies the left hand list of factors into the right hand list
    of factors. Returns (comm. factors list, nc factors list as tuples)."""

    # Shorthand to distribute cl*fL*(eadd.args[0]+eadd.args[1]+..)*fR
    def lapply_DoAddL(eadd:Add, fadd:Callable):
        # Important: arguments to lapply_mul2 must be fresh copies
        # TBD: implement a cache to avoid frequent evaluation of fL?
        # if options["Add"] is set, use this as the add function, else
        # use the original function of the eadd argument. Set option
        # evaluate=True explicitly in case default is False (e.g. MatAdd)
        fadd = options.get("Add", fadd)
        if options["mul"]: # multiply **each summand** with cl
            cle = [lapply1_type(cf, **options) for cf in cl]
            cle = lapply_expand_powers(cle, **options) # by power_base/_exp
            c = lapply_expand_mul(Mul(*cle))
            cle, r = [], fadd(*[lapply_mul2(fL + [c * arg], fR.copy(),
                            **options) for arg in eadd.args], evaluate=True)
        else: # distribute factors, but keep commutative factors up front
            cle, r = cl, fadd(*[lapply_mul2(fL + [    arg], fR.copy(),
                            **options) for arg in eadd.args], evaluate=True)
        # r is a sum of products that are output of lapply_mul2 and consists
        # of elementary factors that don't interact. If r has become a Mul,
        # the nc factors are thus considered non-interacting elems
        clr, rl = r.args_cnc() if r.is_Mul else ([], [r])
        return [lapply1_type(cf, **options) for cf in cle+clr], rl


    cl : List[Expr] = [] # gathers commutative factors

    lhsR = EmptyMul # init leftmost  left  hand side factor
    rhsL = EmptyMul # init rightmost right hand side factor

    while(True):         # replacing recursion on same level by iteration
        if lhsR is EmptyMul:     # Need to pick factor from left hand side
            if fL == []:         # but there aren't any more
                if rhsL is not EmptyMul:
                    fR.insert(0, rhsL) # push rhsL factor not yet in fR
                cl = [lapply1_type(cf, **options) for cf in cl]
                cl = lapply_expand_powers(cl, **options)
                return (cl if cl != [S.One] else []), fR
            # else
            lhs0 = fL.pop()      # get right-most factor from fL
            # split lhs0 up into (also applying .doit()):
            # c factors lhs_c, left factor list lhsL, right elem factor lhsR
            # so that Mul(*lhsL)*lhsR == lhs  and  Mul(*lhs_c)*lhs == lhs0
            lhs_c, lhsL, lhsR = lapply_c_nc_ncef(lhs0, to_the_right, **options)
            cl.extend(lhs_c)     # save c factors, if any
            if lhsR is EmptyMul: # lhs0 went into lhs_c and cl completely
                continue         # so pick next factor
            if lhsR is S.Zero:
                return [S.Zero], []
            fL.extend(lhsL)      # push back lhsL factors, if any

            if lhsR.is_Add:      # do expansion of Add over cl,fL,R
                if rhsL is not EmptyMul:
                    fR.insert(0, rhsL)  # push back rhsL factor if any
                return lapply_DoAddL(cast(Add, lhsR), lhsR.func)

        # Now lhsR is a non-trivial nc factor and no Add. Get rhsL:
        if rhsL is EmptyMul:
            if fR == []:    # if no factors in fR:
                rhsL = lhsR     # shift lhsR (elementary) over to rhsL
                lhsR = EmptyMul # request new lhsR instead
                continue
            # else: get next left most factor from fR
            # fR is built by elementary nc factors, no Add, Mul or c factors
            rhsL = fR.pop(0)    # elem, as only elems are pushed to fR

        # the heart of lapply: we multiply two elementary factors with no
        # obvious scalar factors within them:
        res0 = lapply_mul2elems(lhsR, rhsL, **options)
        nia = isinstance(res0, Expr) and lapply_no_act(lhsR, rhsL, res0)
        if nia: # the two factors lhsR*rhsL have not interacted
            # push rhsL to fR & continue with a fresh fL factor: fL*lhsR
            fR.insert(0, rhsL) # push rhsL (is elem) to fR
            rhsL = lhsR        # continue with elementary lhsR as new rhsL
            lhsR = EmptyMul    # request a new lhsR from fL
        else: # since interaction lhsR*rhsL->res0, continue res0*fR
            if type(res0) is list:
                fL.extend(res0[:-1]) # push back left factors
                res0 = res0[-1]      # focus on right most -> Expr
            res_c, resL, resR = \
                    lapply_c_nc_ncef(res0, to_the_right, **options) # type: ignore[arg-type]
            if resR is S.Zero: # the case resR==1 handled below
                return [S.Zero], []
            cl.extend(res_c)   # save away any factors for Mul

            # if resR is Add (implies resR != rhsL), distribute it
            if resR.is_Add:      # distribute cl,fL,fR over resR
                fL.extend(resL)  # push resL factors
                return lapply_DoAddL(cast(Add, resR), resR.func)

            if resR == rhsL: # (right factor of res0) == rhsL
                if resL == [lhsR]: # we have regained lhsR*rhsL from res0!
                    # so push rhsL to fR and continue with fL*lhsR
                    fR.insert(0, rhsL) # = resR, push rhsL (is elem) to fR
                    rhsL = lhsR        # = resL = lhsR (elem)
                    lhsR = EmptyMul    # request a new lhsR
                else: # resL == lhsR but in different representation,
                    # so compute resL*rhsL
                    fL.extend(resL)    # resL == 1 will be eliminated in Mul
                    lhsR = EmptyMul    # get new factor from fL
            else: # lhsR*rhsL -> res_c*resL*resR, so continue with resR*fR
                fL.extend(resL)
                lhsR = resR            # if resR==EmptyMul gets fL factor
                rhsL = EmptyMul        # get new factor from fR

#--------+---------+---------+---------+---------+---------+---------+---------
# lapply_c_nc_ncef traverses into the structure of object, if required, and
# returns the c_nc_ncef factorisation for elementary factor left/right.
# Code relies on the sympy convention that all commutative factors have been
# pulled out from objects like Pow, OuterProduct, TensorProduct, InnerProduct
# etc by their constructors. So if an an expression e is an instance of those
# objects e does not contain such factors.

def lapply_c_nc_ncef(e:Expr, to_L:bool, **options) -> Tuple[List, List, Expr]:
    """returns (list c factors of e, list nc factors of e without elementary
    factor, nc elementary factor according to left/right indication to_L)"""

    # First prepare expression e by factor extraction. May involve doit(),
    # recursion into inner structure by lapply, flattening etc.
    # e should then be ready for c_nc_ncef(e, to_L) to generate factor lists
    # and to pick the left-most resp. right-most elementary nc factor:
    while(True):
        if e.is_Atom: # just a shortcut for frequent cases, Atoms are elems
            return ([e], [], EmptyMul) if e.is_commutative else ([], [], e)
        # Check cache for base of a Pow before checking cache in general
        if e.is_Pow and (cached := get_cached_for_base(e.base, **options)): # type: ignore[attr-defined]
            if cached[0] != RO.UNKN: # base analyzed; return split if feasible
                return lapply_pow_roll1_c_nc_ncef(e.base, e.exp, to_L,      # type: ignore[attr-defined]
                                           False, cached, **options)
            e1 = e
            break
        # is e itself in cache - which would prove it has been processed?
        cached = options["cache_for_pow"].get(e, None)
        if cached == "elem":
            return [], [], e  # shortcut for elements cached as elementary
        elif cached:          # rare case:
            e1 = e            # e has hit the cache for exponents and bases
            break

        # Not in cache: regular processing of e: clean up and recurse into e
        e1 = lapply1_type(e, **options) # dispatch on type e
        if isinstance(e1, type(e)):     # e1 is same or derived from type(e)
            break
        else:
            e = e1  # we must dig deeper into e
    # Now get the factorisation into c, nc, and left/right-most elementary
    res = c_nc_ncef(e1, to_L=to_L, **options)  # dispatch on type of e1
    return res


#--------+---------+---------+---------+---------+---------+---------+---------
# lapply_mul2elems
#------------------------------------------------------------------------------

# lapply_mul2elems applies a set of application rules that Mul doesn't:
# The main issue here is to deal with the plentitude of options that SymPy
# and the packages provide how und where such rules for
# 'lhs times rhs' may be defined and to predict the result that would have:
# 1. lhs*rhs invokes lhs.__mul__ and then rhs.__rmul__  if defined in class
#    or parents
# 2. Mul(lhs, rhs) may behave differently from lhs*rhs (e.g. Bra(1)*Ket(1)
#    vs. Mul(Bra(1),Ket(1))), or Quaternion*Quaternion vs. Mul(Q,Q)
# 3. In quantum: lhs._apply_operator_<type(rhs)>(rhs) or
#    rhs._apply_from_right_to_<type(lhs)>(lhs) if defined
#    (but note: only for perfect fit of types; no searching up in parents)
# 4. lhs.__pow__(rhs), lhs._pow(rhs) if defined in class or parents
# Moreover:
# 5. some classes provide xxx_simp or .flatten methods or other
#    special handling routines for mul and power, e.g. TensorProduct

def lapply_mul2elems(lhs:Expr, rhs:Expr, **options) -> Union[Expr, List]:
    """Applies elementary lhs and rhs to each other. May return result expr
    or, in case an interaction of lhs and rhs is confirmed, a list of
    factors suitable as first argument to lapply_mul2."""
    res = cast(Any, None)
    # Try layer A of rules lapply_mul2types dispatched by types:
    if (handler := lapply_mul2types_A.dispatch(type(lhs), type(rhs))):
        if (res := handler(lhs, rhs, **options)) is not None:
            return res

    # Try layer B of rules lapply_mul2types dispatched by types.
    # These are rules of lower priority, also rules that would be
    # flagged as ambiguous to layer A rules
    if (handler := lapply_mul2types_B.dispatch(type(lhs), type(rhs))):
        if (res := handler(lhs, rhs, **options)) is not None:
            return res

    # If still no result, default to lhs*rhs using *-operator, which may be
    # overloaded, so e.g. bra*ket->InnerProduct, ket*bra->OuterProduct,
    # Quaterion*Quaternion -> Quaternion, etc. etc. At least Mul(lhs,rhs).
    return (lhs * rhs) if res is None else res

##    ##  ---------------------------------------------------------------------
##    ##
##    ##  Handlers for lapply_mul2types_*(lhs, rhs, **options)
##    ##
########  ---------------------------------------------------------------------
########  Handlers may return result expr, or None, or, if an interaction of
##    ##  lhs and rhs is confirmed, a list of factors suitable as input to
##    ##  first argument of lapply_mul2(). If the handler is able to determine
##    ##  that lhs and rhs will not interact for mathematical reasons is should
##    ##  return Mul(lhs, rhs, evaluate=False)

lapply_mul2types_A = Dispatcher("lapply_mul2types_A")

@lapply_mul2types_A.register(Pow, (Pow, Expr))
def hdlrPE_for(lhs:Pow, rhs:Expr, **options) -> Union[Expr, List, None]:
    # Note: Handler may assume that lhs has been cached
    # Note: Matrix exponents or nc exponents in general not implemented.
    if not lhs.exp.is_commutative:
        return None

    if lhs.base == rhs: # special case rhs = lhs.base
        return Pow(lhs.base, lhs.exp + 1) # just increment exponent
    elif rhs.is_Pow and lhs.base == rhs.base and rhs.exp.is_commutative: # type: ignore[attr-defined]
        # Both exponents are commutative, so adding is possible:
        return Pow(lhs.base, lhs.exp + rhs.exp) # type: ignore[attr-defined]

    # Additional feature over powsimp(lhs*rhs): check interaction
    if (lhs.exp - 1).is_positive: # and also finite
        cached = get_cached_for_base(lhs.base, **options)
        if not cached: raise NotImplementedError("lhs not cached")
        cl, ncl, al = lapply_pow_roll1_c_nc_ncef(lhs.base, lhs.exp,
                to_the_right, force=True, cached=cached, **options)

        alrhs = lapply_mul2([al, rhs], [], **options)
        if not lapply_no_act(al, rhs, alrhs):
            # interaction; return resulting factors as list..
            return cl + ncl + [alrhs]    # ..instead of Mul()
    return None

@lapply_mul2types_A.register(Expr, Pow)
def hdlrEP_for(lhs:Expr, rhs:Pow, **options) -> Union[Expr, List, None]:
    # Note: Handler may assume that rhs has been cached
    # Note: Matrix exponents or nc exponents in general not implemented.
    if not rhs.exp.is_commutative:
        return None
    if rhs.base == lhs: # special case lhs = rhs.base
        return Pow(rhs.base, rhs.exp + 1) # just increment exponent

    # Additional feature over powsimp(lhs*rhs): check interaction
    # of lhs*rhs.base
    if (rhs.exp - 1).is_positive: # and also finite
        cached = get_cached_for_base(rhs.base, **options)
        if not cached: raise NotImplementedError("rhs not cached")
        cr, ncr, ar = lapply_pow_roll1_c_nc_ncef(rhs.base, rhs.exp,
                to_the_left, force=True, cached=cached, **options)
        lhsar = lapply_mul2([lhs, ar], [], **options)
        if not lapply_no_act(lhs, ar, lhsar):
            # interaction; return resulting factors as list..
            return cr + [lhsar] + ncr   # ..instead of Mul()
    return None

# Create (empty) Dispatcher for layer B of rules:
# Layer B contains handlers that either have lower priority
# than layer A, or would create a type ambiguity in Dispatcher.
lapply_mul2types_B = Dispatcher("lapply_mul2types_B")
# See sympy.physics.quantum.qapply for examples for layer B.


# -----------------------------------------------------------------------------
# Cache for Pow
# -----------------------------------------------------------------------------
# The purpose of caching is to avoid frequent recomputation if a Pow is
# rolled off to Pow(base,exp-1)*base -> Pow(base,exp-2)*base*base etc.
#
# Principles of Caching for base and exponents of Powers:
# - Only cache exprs that will probably be used multiple times. Do not litter the
#   the cache with one-time results. So bases of powers are cached. Exponents are
#   cached only in case of option["apply_exp"] and non-simple exponents.
# - Exprs being cached must have been returned in ncl from lapply_mul2_c_nc, so
#   all its factors have been processed fully.
# - The cache returns three types of information on an expr (besides None):
#     i.   "elem": expr is elementary (in the sense defined above)
#     ii.  "expr": expr has been processed fully, all factors elementary
#     iii. (rolloff, cl, encl) expr is base of a power
#   which roughly correspond to output of c_nc_ncef, lapply_mul2, lapply_mul2_c_nc.
# - Base and exponents of Pows are cached separately. If a power is rolled off
#   the exponents counted down are removed from cache to avoid litter. This may
#   lead to re-computation of exponents in situations as Pow(A,exp)*B+Pow(A,exp)*C
#   which is unavoidable as the cache only has one instance.

class RO(Enum): # the roll-off type of a base:
    UNKN = 0    # Roll-off type is unknown
    NORO = 1    # base*base doesn't interact, so don't roll off base**exp
    ROLL = 2    # base*base interacts, so roll off base**exp off

# For a base, cache the rolloff decision assigned to it (rolloff) and the
# splitting of base into cl and elementary ncl factors as derived e.g.
# from lapply_mul2_c_nc
def cache_for_base(base:Expr,
                   rolloff:RO, cl:List, encl:List, **options) -> Tuple:
    """Generates the info to be cached for base, updates the cache with
    it and returns the tuple that has been cached."""
    if base.is_commutative: # in this case overwrite any other rolloff
        # note that factors in encl may be nc but give a c product!
        rolloff, cl, encl = RO.ROLL, [base], [] # let Mul handle c base.

    to_cache = [rolloff, cl.copy(), encl.copy()]
    cache = options["cache_for_pow"]
    # The cache is indexed on base, as only attributes of base are stored
    cached = cache.setdefault(base, tuple(to_cache))
    if cached == to_cache:          # cache is now what is should be
        return cached
    if type(cached) is tuple:       # cache had another full base info
        if to_cache[0] == RO.UNKN:  # Existing Rolloff in cache may not
            to_cache[0] = cached[0] # be downgraded to RO.UNKN
    # else: cache[base] is "elem" or "expr". This will be updated
    # to tuple(to_cache) as it gives more information.

    cached = cache[base] = tuple(to_cache) # (rolloff, to_cache[1], to_cache[2])
    return cached

# Pow_with_Cache packs a call cache_for_base() in front of Pow(base, exp)
def Pow_with_cache(base:Expr, exp:Expr,
                   rolloff:RO, cl:List, encl:List, **options) -> Expr:
    """Returns Pow(base, exp), but also calls cache_for_base() to
    cache factorization data (rolloff, cl, encl) of base."""

    cache_for_base(base, rolloff, cl, encl, **options)

    if   exp is S.Zero: return S.One
    elif exp is S.One : return base
    else              : return Pow(base, exp, evaluate=False)


def get_cached_for_base(base:Expr, **options) \
                            -> Union[Tuple[RO, List, List], None]:
    """Get cached contents (RO, cl, encl) for base, None if none"""
    cached = options["cache_for_pow"].get(base, None)
    if type(cached) is tuple:
        return (cached[0], cached[1].copy(), cached[2].copy())
    elif cached == "elem": # happens to have been cached as elementary
        return (RO.UNKN, [], [base])
    return None

# Check if expr is in lapply_cache; if yes, return it, if no
# apply lapply_mul2, cache the result and return it.
def lapply_mul2_cache(expr:Expr, **options) -> Expr:
    """Compute lapply_mul2([expr],[],**options) and cache the result"""
    if options["cache_for_pow"].get(expr, None) is None:
        expr =  lapply_mul2([expr], [], **options)
        options["cache_for_pow"].setdefault(expr, "expr")
    return expr

def lapply_pow_roll1_c_nc_ncef(base:Expr, exp:Expr, to_L:bool, force:bool,
                         cached:Tuple[RO, List, List], **options) \
                             -> Tuple[List, List, Expr]:
    """Roll-off one base factor from Pow(base, exp) and return the c_nc_ncef
    factor representation of Pow(p.base, p.exp-1)*p.base resp. of p.base*\
    Pow(p.base, p.exp-1) from the cached data"""

    if options["apply_exp"] and not exp.is_Atom:
        exp = lapply_mul2_cache(exp, **options)
    if base.is_commutative and exp.is_commutative:
        return [Pow(base, exp)], [], EmptyMul

    # Return if exp is not verified >=1 or (type=NORO and not force).
    # Note that we first must assure that exp - S.One is defined at all
    # before computing exp1 = exp - 1:
    if ((not exp.is_positive) or (not (exp1 := exp - S.One).is_nonnegative) # type: ignore[attr-defined]
         or (cached[0]==RO.NORO and not force)):
        return [], [], Pow(base, exp, evaluate=False)

    ncl = cached[2] # ncl has at least 1 factor, as base is nc
    pl = []         # result power in list
    if   exp1 == S.One : ncl += ncl # exp == 2
    elif exp1 == S.Zero: pass       # exp == 1
    else:
        pl = [Pow(base, exp1, evaluate=False)]
        # if exp1 is complex to compute, cache it and remove
        # the previous exp so cache is not littered
        if options["apply_exp"] and not exp1.is_Atom:
            cache = options["cache_for_pow"]
            if cache.get(exp, None) == "expr":
                del cache[exp] #don't delete in other cases!
            cache.setdefault(exp1, "expr") # register exp1 as lapplied

    # hand out c-factors, nc-factors and nc elementary factor
    ncl2, a = (ncl[1:], ncl[0]) if to_L else (ncl[:-1], ncl[-1])
    res = cached[1], (ncl2 + pl) if to_L else (pl + ncl2), a
    return res

##    ##  ---------------------------------------------------------------------
##    ##
##    ##  @dispatch Handlers for lapply1_pow_type(base, exp, **options)
##    ##
########  see lapply.py: ------------------------------------------------------
########  This handlers are invoked by the lapply1_type(Pow) to handle a power
########  Pow(base, exp) with a particular type of base and/or exponent.
##    ##  Handlers should take care of non-evaluated objects and recurse into
##    ##  inner structure of e, if any, and process it so that the result is
##    ##  flattened and c_nc_ncf() may extract factors from the result in a
##    ##  next step.

# No handler present here. See sympy.physics.quantum.qapply for an example:
#
#@lapply1_pow_type.register(TensorProduct, Expr)
#def hdlr_TE(base:TensorProduct, exp:Expr, **options) -> Expr:
#    return TensorProduct(*[lapply_mul2([b**exp], [], **options)
#                                    for b in base.args]        )

lapply1_pow_type = Dispatcher('lapply1_pow_type') # create empty Dispatcher


##    ##  ---------------------------------------------------------------------
##    ##
##    ##  Handlers for lapply1_type(Expr, **options)
##    ##
##    ##  ---------------------------------------------------------------------
########  Handlers should take care of non-evaulated objects and recurse into
########  inner structure of e, if any, and process it so that the result is
##    ##  flattened and c_nc_ncf() may extract factors from the result in a
##    ##  next step.
##    ##  Note: These handlers are also called for commutative factors! In
##    ##  most cases no action is wanted as lapply leaves c factors to Mul.
##    ##  So put "if e.is_commutative: return e" on top.

lapply1_type = Dispatcher("lapply1_type")

@lapply1_type.register(Mul)
def hdlrM_for(e:Mul, **options) -> Expr:
    # Unevaluated Muls exist by explicit user option only, so it is for
    # a reason. Therefore I decided not to call .doit() indiscriminately.
    # This means that Mul may be in non-canonical form and its args
    # may need sorting for c, nc. Calling args_cnc() might fail on non-
    # evaluated Muls (issue #24480), so we mimick args_cnc() here:
    cl  : List[Expr] = []
    ncl : List[Expr] = []
    for a in e.args:
        if a.is_commutative: cl  += a,  # quicker than cl.append(a)
        else:                ncl += a,
    return Mul(*(cl+ncl), evaluate=False)


# Note on Sum() and Product():
# Do not process them, even if the ranges are numeric. This explicit symbolic
# representations as Sum() or Product() exists by user decision with a purpose
# (see e.g. quantum JxKet etc.). If unintended, user may apply .doit() before
# passing them to lapply.


# Internal Remark: lapply and powers Pow(base, power)
# ---------------------------------------------------
#
# 1) e as a Pow object is already result of the constructor Pow.__new__(base,
#    power) that contains at lot of logic to compute this power, including
#    calling base._eval_power(power). Pow also has pulled (recognizable)
#    commutative factors from its base into a separate Pow object.
#    Some classes (e.g. in quantum OuterProduct and TensorProduct) lack
#    adequate methods .__pow__, ._pow or ._eval_power to apply powers.
#    So where lapply needs to recursively enter into powers of these types,
#    add a dispatched method lapply1_pow_type for the type of base.
# 2) lapply will try some 'tricks' to simplify certain powers that Pow doesn't
#    use, nor Pow._eval_expand_power_xxx, e.g. on idempotent operators/dyads.
# 3) But the last resort of lapply is to unroll positive powers of an operator,
#    i.e. apply it one by one to a vector v: For dense A, to compute Pow(A,k)*v
#    via M:=Pow(A,k) and then M*v is much(!) more expensive than (A*..*(A*(A*v))),
#    as long as log2(k)*dim(A)>k - which should be the case considering that A
#    annuls its characteristic polynomial. So, if A*v computes to some b which
#    is not the expression Mul(A,v), we compute (A*..*(A*(A*v))). Only if A*v
#    doesn't compute, we try to compute Pow(A,k) via A*A*..*A*A. If we find that
#    A*A doesn't compute, Pow(A,k) is left alone.
#    In order to make this more efficient, the base is evaluated by lapply only
#    once and then cached (calling Pow_with_Cache(), see caching above).
# 4) But as this is an heuristic only: If higher powers of an operator op occur,
#    the **best approach** is that the user provides an appropriate method in
#    the class definition of op and makes sure that Pow(op, exp) evaluates this.

@lapply1_type.register(Pow)
def hdlrP_for(e:Pow, **options) -> Expr:
    """Prepare expression Pow(base, exp) for further application by
    going into base and exp."""
    if e.is_commutative: return e

    cached = get_cached_for_base(e.base, **options)
    # dispatch on the type of base and exp of the power, if handler:
    if not cached:
        handler = lapply1_pow_type.dispatch(type(e.base), type(e.exp))
        if handler:
            if (res := handler(e.base, e.exp, **options)) is not None:
                return res

    # Note: we don't call lapply or expand on e.exp by default, as there seems
    # little use for it in the context of lapply: Powers indicate repetitions
    # of operators, and what sense would operators in exponents make? Of course
    # terms like e**A do have meaning (e.g. Hamiltonial evolution in quantum
    # mechanics), but in many fields this terms would not occur. So I make
    # that an optional feature (option "apply_exp") with default = False.
    exp =      lapply_mul2_cache(e.exp, **options) if options["apply_exp"] \
          else e.exp

    # try to turn a neg exponent positive, if base has an explicit inverse,
    # so we may unroll the power. Pow() doesn't try this on creation of a Pow
    # object. Requires exp to be integer to be always correct mathematically.
    if e.exp.is_negative and e.exp.is_integer and hasattr(e.base, 'inv'): # type: ignore[attr-defined]
        binv = e.base.inv().doit() # type: ignore[attr-defined] #for inv()
        # We require an explicit inverse of base, not just e.base**-1. Also
        # detect fake inverses that are not Pow instances (e.g. MatPow):
        e1 = binv ** (-e.exp) # so, give binv a try
        if e1 != e and not \
            (hasattr(e1, "args") and e1.args and e1.args[0] == e.base):
            # Detect if no longer a power. Should be "not e1.is_Pow", but use
            # a hack to detect MatPow that is not is_Pow since PR #1636
            if not hasattr(e1, "base") or not hasattr(e1, "exp"):
                return e1     # return it with new type.
            if e1.exp.is_nonnegative:  # so if exp turned nonnegative
                # continue with the new e=e1 and new positive exponent
                e = e1; exp = e.exp
                cached = get_cached_for_base(e.base, **options)
    base = e.base
    if not cached:
        # if no methods has given a result so far, at least recurse into e.base.
        cl, encl = lapply_mul2_c_nc([e.base], [], **options)
        for a in encl: options["cache_for_pow"].setdefault(a, "elem")
        cml = [lapply_expand_mul(Mul(*cl))] if options["mul"] else cl
        base = Mul(*(cml + encl))
        cache_for_base(base, RO.UNKN, cl, encl, **options)

    if base != e.base or exp != e.exp:
        # Let ** do all it can do using class methods ._pow, ._eval_power etc.
        e = cast(Pow, base ** exp)  # may return any Expr; cast() is for MyPy
        if e.is_Pow:  # repeat until base and exp no longer change.
            return lapply1_type(e, **options)
            # also ensures that e.base is cached and factorized
    return e


# catch-all of lapply1_type for an arbitrary class e where we know
# nothing about its structure and factors. So return e unmodified.
# dispatch will throw "NotImplemented" if a non-Expr appears - which
# is actually impossible for lapply to handle as a non-Expr don't have
# arithmetic operations (in SymPy, everything that has arithmetic
# operations must be derived from Expr).
@lapply1_type.register(Expr)
def hdlrE_for(e:Expr, **options) -> Expr:
    return e


##    ##  ---------------------------------------------------------------------
##    ##
##    ##  Handlers for c_nc_ncef(Expr, to_L=..., **options)
##    ##
##    ##  ---------------------------------------------------------------------
########  Split up expression e in commutative factors and non-commutative
########  factors. From the nc factors pick the leftmost (if to_L==True) or
##    ##  rightmost (to_L==False) elementary factor ncef and return (c_list,
##    ##  nc_list_without_ncef, ncef). If  e is commutative, return (c_list, [],
##    ##  S.One). If e cannot be split up, return ([], [], e).
##    ##  Important: Call with argument to_L= explictly named otherwise the
##    ##  Dispatcher will try to dispatch on first and 2nd argument!

c_nc_ncef = Dispatcher("c_nc_ncef")

@c_nc_ncef.register(Mul) # Call as c_nc_ncef(e, to_L=..., **options)
def hadlrM_for(e:Mul, to_L:bool, **options) -> Tuple[List, List, Expr]:
    # !!Code relies on e.args being sorted commutative, non-commutative,
    # i.e. canonical form of Mul. May fail if Mul(.., evaluate=False)!
    if e.args[0].is_commutative:    # not all non-commutative?
        # search backwards for last commutative factor
        for (i, a) in enumerate(reversed(e.args)):
            if a.is_commutative:
                i = len(e.args) - i  # i is index of first nc factor
                break
    else: # no commutative element in e.args at all
        i = 0
    cl, ncl, a = list(e.args[0:i]), list(e.args[i:]), EmptyMul
    while ncl != [] and a is EmptyMul:
        # extract elementary factor from left or right side of ncl:
        # pick the left- resp. right most factor f from ncl
        ncl, f = (ncl[1:], ncl[0]) if to_L else (ncl[:-1], ncl[-1])
        # get leftmost/rightmost factor of f
        fcl, fncl, a = lapply_c_nc_ncef(f, to_L=to_L, **options)
        cl.extend(fcl)
        if to_L: ncl[:0] = fncl         # insert fncl to the left
        else   : ncl.extend(fncl)       # append fncl to the right
    return cl, ncl, a


# ### Simplification Rules for Powers ####
# same as used by expand() or by Pow() with numeric exponents:
# (Rules differ from https://docs.sympy.org/latest/tutorials/intro-tutorial/simplification.html#powers):
# 1.  for arbitrary A, and a, b commutative:
#     A**a * A**b == A**(a+b) == A**(b+a)
# 2a. for c commutative, A (non-)commutative, n integer (positive or
#     negative assuming A**-1, c**-1 exists):
#     (c*A)**n == c**n * A**n
# 2b. for c nonnegative, A (non-)commutative, a commutative:
#     (c*A)**a == c**a * A**a
# 2c. for 0<=r<1 (from 1 and 2a):
#     (c*A)**(n+r) == (c*A)**r * (c*A)**n == c**n * (c*A)**r * A**n
# 2d. for c and d nonnegative, a (non-)commutative:
#     (c*d)**a == c**a * d**a  (is commutative if a is)
# 3a. (A**b)**a = A**(a*b) if a integer, A,b arbitrary
# 3b. if c is nonnegative, b is nonnegative and a (non-)commutative:
#     (c**b)**a = c**(a*b)  (is commutative if a is) (see 2d)
# 3c. for 0<=r<1 and n integer (from 1 and 3a):
#     (A**a)**(n+r) == (A**a)**r * (A**(a*n))


@c_nc_ncef.register(Pow) # Call as c_nc_ncef(e, to_L=..., **options)
def hadlrP_for(e:Pow, to_L:bool, **options) -> Tuple[List, List, Expr]:
    """Try to simplify power e it by probing factors of base and of base*base.
    Return it as factors c-factors (if known), nc factors,
    right hand elementary nc factor"""
    # if e.is_commutative: return [e], [], EmptyMul
    if e.is_commutative: # don't touch it. But c Pows that come here are the only
        # factor in the expr, probably for a reason. So make an exception:
        return [Pow(lapply_expand_mul(e.base), e.exp, evaluate=False) \
                if options["mul"] else e], [], EmptyMul

    # Requires that e.base has been cached and factorized in lapply1_type(Pow)
    if (cached := get_cached_for_base(e.base, **options)) is None:
        raise NotImplementedError(f"c_nc_ncef(Pow({e.base}, {e.exp})) not cached.")

    # Two Helper functions pow_x, pow_i to implement option "nested_exp":
    #   Raise elements c of a list cl to power expx: If c happens to be a
    #   Pow(c.base, c.exp) SymPy will automatically write c**exp as
    #   Pow(c.base, c.exp*exp) if that is feasible. However it won't expand
    #   c.exp*exp. So if option "nested_exp" is on, do it here (delaying it
    #   to lapply_expand_power() would require another check and another
    #   call of the Pow constructor).
    if options["nested_exp"]:
        # expand the product of exponents if condition cond is met;
        # cond describes feasibility for multiplication with expx:
        pow_x = (lambda cl, expx, cond:
            [( (c.base ** lapply_expand_mul(c.exp * expx))
               if c.is_Pow and cond(c)  else  c ** expx)     for c in cl])
        pow_i = (lambda cl, expi: pow_x(cl, expi, (lambda c: True)))
    else: # else return c**exp unsimplified
        pow_x = (lambda cl, expx, cond: [c ** expx for c in cl])
        pow_i = (lambda cl, expi:       [c ** expi for c in cl])

    ##### End Helper functions ###### Start function body ###### ---------------

    if not e.exp.is_commutative: # so rules 2a, 2b don't apply
        # and we cannot extract any factor from e.base**exp
        return [], [], e

    # e.exp is commutative, so we may apply 2b or even 2a.
    # get factors of e.base = cl*encl
    (cl, encl) = cached[1:]
    if encl == []: # implies e.base is commutative -> rolloff=RO.ROLL
        cache_for_base(e.base, RO.ROLL, cl, encl, **options)
        return [e], [], EmptyMul # Power is commutative
    if to_L:     # this function can handle split_to_right only
        raise NotImplementedError('Splitting Powers to the left should' +
                                ' not occur and is not implemented.')

    if cl == [] and e.exp.is_integer:        # type: ignore[attr-defined]
        return c_nc_pow_nc_i(e.base, encl, e.exp, **options)

    if not e.exp.is_real:
        # case 8: no knowledge about e.exp except that it is commutative
        # so we may apply 2b and pull out c**e.exp for nonneg factors
        clnonneg = [c for c in cl if c.is_nonnegative]
        clelse   = [c for c in cl if not c in clnonneg]
        # ep = (clelse*base)**e.exp
        ep = cast(Expr, Pow_with_cache(Mul(*(clelse + encl)), e.exp,
                                   RO.UNKN, clelse, encl, **options))
        # we can only apply Rule 3b here and not 3a as exp is no int:
        return pow_x(clnonneg, e.exp, (lambda d:
                d.base.is_nonnegative and d.exp.is_nonnegative)), [], ep

    # Split real exponent e.exp in part expi and "fractional" part:
    # e.exp = expi + expf. expi = largest integer <= e.exp that is
    # known to be integer, expf is real, but probably fractional.
    expi = floor(e.exp) # makes expi.is_integer == True
    # (comment this if you want symbolic calculation with floor(e.exp):)
    if expi.has(floor): # there is a floor(x) with non-integer x
        if    expi.func == floor: expi = S.Zero # no tangible integer
        elif  expi.is_Add: # extract explicit numeric integers only and
            # leave symbolic integers to expf for more natural results
            expi = Add(*[a for a in expi.args if a.is_Integer])
        else: raise NotImplementedError(f"Cannot extract int part of {expi}")
    # part with potential fraction and unknown attr. 0 if e.exp.is_integer.
    expf = e.exp - expi

    # Describe the fractional power part e.base**expf:
    # create efpcl, efpl so that e.base**expf == Mul(*efpcl) * Mul(*efpl)
    # Even if expf == e.exp we need to split cl up
    efpl = efpcl = []
    if expf != S.Zero: # and cl != []
        # pull non-neg factors from cl up front (rule 2b)
        clnonneg = [c for c in cl if c.is_nonnegative]
        if clnonneg != []:
            efpcl  = pow_x(clnonneg, expf, # [c**expf for c in clnonneg]
                (lambda c: c.base.is_nonnegative and c.exp.is_positive))
            clelse = [c for c in cl if not c in clnonneg]
            # efp = (clelse*base)**expf holds remaining factors from cl
            efp    = Pow_with_cache(Mul(*(clelse + encl)),
                                expf, RO.UNKN, clelse, encl, **options)
            efpl   = [] if efp is S.One else [efp]
        else:
            efpl   = [Pow(e.base, expf, evaluate = False)]
    # From here on, the fractional power e.base**expf is accounted for
    # by adding efpcl to the commutative factors, and efp to the nc-factors
    if expf == e.exp: # implies expi == 0, so e == clefpl
        return efpcl, [], Mul(*efpl)

    # handle Pow(encl, expi) with encl factors and expi int:
    base = Mul(*encl) # e.base without cl
    c2, nc2, a2 = c_nc_pow_nc_i(base, encl, expi, **options)
    # Join c2, nc2, a2 with cl, efpcl, efpl factors:
    if a2 == EmptyMul: # implies nc2 == []
        return pow_i(cl, expi) + efpcl + c2, [], Mul(*efpl)
    else:
        return pow_i(cl, expi) + efpcl + c2, efpl + nc2, a2


# Do c_nc_ncef (to_the_right) for Power where base = encl has elementary
# nc-factors only and exponent expi is verified expi.is_integer==True:
def c_nc_pow_nc_i(base:Expr, encl:List, expi:Expr, **options):
    """c_nc_ncef (to the right) for Pow(ncl, expi) with ncl elementary
    nc factors and integer expi"""

    ncl, ar = encl[0:-1], encl[-1] # split encl to the right

    # Handle special cases expi == 0 or expi == 1:
    if expi is S.Zero:
        cache_for_base(base, RO.UNKN, cl = [], encl=encl, **options)
        return [], [], S.One
    elif expi is S.One :
        cache_for_base(base, RO.UNKN, cl = [], encl=encl, **options)
        return [], ncl, ar

    #### Now: expi is integer and != 0 and != 1 (secured by if's above)
    # Code below relies on this condition to be mathematically correct!
    assert expi.is_integer   # type: ignore[attr-defined]

    def pb(expi, rolloff):   # shorthand
        return Pow_with_cache(base, expi, rolloff=rolloff,
                                        cl = [], encl=encl, **options)

    # Handle nested powers (see rule 3a, 3c): base=ar, ncl = []
    # (powdenest() can't handle symbolic exponents and doesn't expand exp.)
    if ar.is_Pow and ncl == []:
        # ar**expi = ar.base**(ar.exp*expi):
        # De-nest nested power and split it up into factors:
        return lapply_c_nc_ncef(
            Pow(ar.base, lapply_expand_mul(ar.exp * expi)),
                                            to_the_right, **options)

    if not ((expi - 1).is_nonnegative): # expi not confirmed >= 1
        # Case 7: int expi, not 0, 1, no sign confirmed, base=ncl*ar
        return [], [], pb(expi, RO.UNKN)

    #### Case: expi is integer >= 1 (secured by If's above)
    ####
    # Code below relies on this condition to be mathematically correct!
    assert (expi - 1).is_nonnegative # Superfluous. Optical reminder only

    # pow_i is shorthand to raise comm. factors in cl to expi:
    if options["nested_exp"]:
        pow_i = (lambda cl, expi:
            [( (c.base ** lapply_expand_mul(c.exp * expi))
               if c.is_Pow  else c ** expi)                  for c in cl])
    else:
        pow_i = (lambda cl, expi: [c ** expi for c in cl])

    # if rolloff state of base has already been determined, use that:
    if (cached := get_cached_for_base(base, **options)) is not None:
        if cached[0] != RO.UNKN:
            return lapply_pow_roll1_c_nc_ncef(base, expi,
                            to_the_right, False, cached, **options)

    # Else determine rolloff state of base and potential simplifications:
    if ncl == []: # base=ar is elementary, expi Int >= 1
        ar2cl, ar2ncl = lapply_mul2_c_nc([ar, ar], [], **options)
        if len(ar2ncl) >= 2: # Elementary factors remain.
            # Unrolling makes no sense, as factor count will just grow.
            # Case 4: expi>=1, base=ar elementary and no base*base
            return ([], [], pb(expi, RO.NORO))

        # is ar "involutoric" operator resp. ar*ar is commutative?
        if ar2ncl == []: # ar*ar=ar2cl is even commutative!
            # Case 5a: expi >= 1, base=ar elementary + "involutoric"
            # ar**expi=(ar2c**(expi//2)) * (ar ? (expi%2==1) : 1)
            # Note: If expi is no explicit integer then expi % 2 is
            # symbolic and the formal Pow(ar, expi%2) cannot unroll.
            return pow_i(ar2cl, expi // 2), [], \
                        S.One if expi.is_even else Pow(ar, expi % 2) # type: ignore[attr-defined]

        # Check if ar is dyad / "idempotent" and let power unroll:
        if ar2ncl[0] == ar :
            # Case 5b: expi >= 1, base=ar elementary + "idempotent"
            # "idempotency": ar*ar=c*ar -> ar**expi=c**(expi-1)*ar
            return pow_i(ar2cl, expi - 1), [], ar

        if ar2cl == [] and lapply_no_act(ar, ar, ar2ncl[0]):
            return ([], [], pb(expi, RO.NORO))

        # else: Case 6: base=ar elementary but not "idempotent"
        if (expi - 2).is_nonnegative: # Idea: split off base**2 at once
            # base**(exp-2) * ar*ar = base**(exp-2)*ar2c*ar2ncl*ar2ar
            # so return ar2cl, [pb(expi - 2, RO.ROLL)]+ar2ncl[:-1],ar2ncl[-1].
            # But: makes result inconsistent from result if handed out by
            # lapply_pow_roll1_c_nc_ncef() for same base, exp. So:
            return [], [pb(expi - 1, RO.ROLL)], ar
        else: # we may split off at least one base factor ar:
            return [], [pb(expi - 1, RO.ROLL)], ar

    #### Case: expi integer >= 1 and ncl!=(), so base has >=2 nc-factors
    ####
    # if base=a*b*c and c*a is commutative, then
    # base**expi=(c*a)**(expi-1) * a * b**expi * c.
    # We don't check all possible factors, just the simplest cases
    # like (A*B*A.inv)**exp or (|k>*<b|)**exp are simplified
    al, ncr2 = ncl[0], ncl[1:] # split ncl to the left

    ia = lapply_mul2([ar, al], [], **options) # ar, al may be Add
    if lapply_no_act(ar, al, ia): # factors didn't interact, so
        # it makes no sense for the power p to unroll, as factors
        # will just line up, so choose rolloff=NORO.
        # Tbd: Also cases ar=CxD, al=AxB, ia=(C*A)x(D*B) and
        # ar=(C1*C2)x(D1*D2),al=(A1*A2)x(B1*B2) should go here!
        # That will need special nia rule here for tensor factors.
        # Case 1: expi>=1, component base, base*base not interact
        return [], [], pb(expi, RO.NORO)

    # factors did interact: check whether ia is commutative
    if (getattr(ia, 'is_commutative', False) or isinstance(ia, Number)):
        # (No need to check for ia == Identity or IdentityOperator
        # or IdentityGate etc, as ia would be S.One instead.
        # base = a*b*c and c*a is commutative, so:
        # base**expi = (c*a)**expi * a * b**expi * c
        # Case 2: expi>=1, comp. base, base*base interact/comm.
        base2 = Mul(*ncr2)
        cache_for_base(base2,
                       rolloff=RO.UNKN, cl=[], encl=ncr2, **options)
        return [ia ** (expi - 1)], [al, Pow(base2, expi)], ar

    # it makes sense to permit the power to unroll.
    # Case 3: expi>=1, component base, base*base interacts
    # but is non-c. base**expi = base**(expi-1) * (ncl * ar)
    return [], [pb(expi-1, RO.ROLL)] + ncl, ar


# catch-all of c_nc_ncef for an arbitrary Expr
# where we know nothing about its structure and factors
@c_nc_ncef.register(Expr) # Call as c_nc_ncef(e, to_L=..., **options)
def hadrlE_for(e:Expr, to_L, **options) -> Tuple[List, List, object]:
    # generic rules that are not type dependendent:
    if (getattr(e, 'is_commutative', False) or isinstance(e, Number)):
        # remaining e is commutative: scalars, InnerProduct, ..
        return ([e], [], S.One)
    else: # e is non-commutative, e.g. in quantum: Ket, Bra, Add, ..
        return ([], [], e)


# ---------+---------+---------+---------+---------+---------+---------+-------
# MM          MM
# MMM        MMM   lapply Handlers for Matrices
# MMMM      MMMM
# MM MM    MM MM --------------------------------------------------------------
# MM  MM  MM  MM  Reminder: SymPy Exprs must be immutable for the SymPy internal
# MM   MMMM   MM  cache to handle them. So the matrices package has type Matrix
# MM          MM  that is mutable and co-exists with SymPy, type MatrixExpr for
# MM          MM  SymPy expressions and ImmutableMatrix that should be used to
# MM          MM  create explicit matrices in SymPy.
#
# Only the handler c_nc_ncef(MatrixExpr) would required to make lapply work with
# matrices, to turn multiples of the identity to scalars which helps a lot in
# multiplication and simplification of powers. It might generate issues with
# terms like (1 - Identity(2)) which however vanish as soon as this term is
# multiplied with something - so this is accepted. Note that MatAdd cannot
# handle scalar summands, but lapply by default uses Add which can.
#
# Two further handlers are required as workaround for MatPow not being derived
# from Pow but MatMul being derived from Mul (see comments below).

from sympy.matrices import (ZeroMatrix, Identity, MatrixExpr,
                            ImmutableMatrix, MatMul, MatAdd)

# c_nc_ncef(MatrixExpr) tries to identify multiples of the Identity.
# It also does .simplify() on the matrix - costly, but necessary for
# meaningful detection of diagonal matrices and multiples of the
# Identity. lapply adapts the convention of quantum package that
# multiples of Identity are reduced to scalars. Also see comment above.

# Handle MatMul as a Mul and not as MatrixExpr
c_nc_ncef.add((MatMul,), c_nc_ncef.dispatch(Mul))

@c_nc_ncef.register(MatrixExpr) # Call as c_nc_ncef(e, to_L=..., **options)
def hdlrME_for(e:MatrixExpr, to_L:bool, **options) -> Tuple[List,List,Expr]:
    if isinstance(e, (MatAdd, MatrixSymbol)): return [], [], e
    # Evaluate e, turn into explicit type, and simplify:
    e = e.doit().simplify()
    # ZeroMatrix and Identity don't have .is_square and .is_diagonal!
    if isinstance(e, ZeroMatrix) or getattr(e, "is_zero_matrix", False):
        return [S.Zero],[], S.Zero
    if isinstance(e, Identity): return [],[], EmptyMul
    # we don't seek a common factor of all matrix elements via
    # gcd(e.args[2]), as that seems to overdo it. Just find k*Identity:
    # .is_square and .is_diagonal available only if e is explicit matrix:
    if isinstance(e, ImmutableMatrix) and e.is_square and \
        e.is_diagonal() and all(e[0,0] == d for d in e.diagonal(0)[1:]):
        if e[0,0].is_commutative: # return (commutative) scalar
            return [e[0,0]], [], EmptyMul
        else: # matrices don't support nc scalars, but maybe in the future
            return [], [], e[0,0] # return non-commutative scalar
    return [], [], e


#######  Workaround for MatPow not being derived from Pow  #######
#
# 1) PR #1636 in 2012 removed inheritance of MatPow from Pow!
# This broke SymPy consistency that Pow is the shorthand of Mul for identical
# factors, and breaks generic code like lapply that must rely on .is_Pow to
# identify powers as result of Mul, and to treat Pow just as a particular
# occurrence of a Mul. See conversation at PR #1636.
# Example: A=MatrixSymbol("A", 2,2): Mul(A,A) -> MatPow(A,2),
# even Pow(A, 2)->MatPow(A,2), but MatPow(A,2).is_Pow==False!

# So a workaround is required make lapply work with matrices:
# a) MatPow is changed to Pow in lapply1_type(MatPow) and piped through
# lapply1_type(Pow), except if its base is symbolic (e.g. MatrixSymbol,
# Transpose, Inverse, DiagMatrix, Adjoint etc.) and cannot be evaluated anyway
# b) Same trick in c_nc_ncef(MatPow); it is processed as Pow and switched back
# to MatPow.
# See comment at the end for further details.

from sympy.matrices import MatPow, MatrixSymbol

# Workaround for MatPow.is_Pow == False for lapply1_type:
@lapply1_type.register(MatPow)
def hdlrMP_for(e:MatPow, **options) -> Expr:
    e = e.doit() # process and get final type, rely on native pow!
    if not isinstance(e, MatPow):
        return e
    if isinstance(e.base, (ImmutableMatrix, MatrixSymbol, MatMul, MatAdd)):
        e = lapply1_type(Pow(e.base, e.exp, evaluate=False), **options)
        if e.is_Pow:     # convert to MatPow if matrices in base
            e = e.doit().doit()  # evaluate; see comment below.
    return e


# Workaround for MatPow.is_Pow == False for c_nc_ncef:
@c_nc_ncef.register(MatPow)
def hdlrMaP_for(e:MatPow, to_L:bool, **options) -> Tuple[List,List,Expr]:
    if isinstance(e.base, (ImmutableMatrix, MatMul, MatAdd)):
        cl, ncl, a = c_nc_ncef(Pow(e.base, e.exp, evaluate=False),
                                                to_L=to_L, **options)
        # Eval pows that have matrices in their base back to MatPows
        if a.is_Pow: a = a.doit().doit()  # see comment below
        ncl = [(f.doit().doit() if f.is_Pow else f) for f in ncl]
        return cl, ncl, a
    return [], [], e

# Since Pow(m, n) converts to MatPow(m, n) if there is Matrix in base m, I
# use Pow.doit() for that conversion. Since MatPow(m, n) has evaluate=False,
# I use another .doit() to trigger evaluation, computation of the power if
# possible and to make types consistent (e.g. converts MatPow(m, -1) to
# Inverse(m) type). Type consistency is crucial for lapply_no_act() to work.
# As every invocation of Pow(Matrix, n) turns Pow into MatPow statements
# "if x.is_Pow ..." won't fire. But since these are mostly needed for
# caching during roll off of powers which is a back up feature if the native
# Pow (i.e. MatPow) doesn't handle powers, so it is rarely needed here.
