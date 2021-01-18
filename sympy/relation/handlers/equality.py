from sympy.assumptions import ask, Q
from sympy.core import Add, Basic, Expr, S, Symbol, Tuple
from sympy.core.kind import NumberKind
from sympy.core.function import AppliedUndef
from sympy.core.logic import fuzzy_bool, fuzzy_and, fuzzy_xor, fuzzy_not
from sympy.functions import arg
from sympy.logic.boolalg import BooleanAtom
from sympy.multipledispatch import Dispatcher
from sympy.simplify.simplify import clear_coefficients
from sympy.utilities.iterables import sift


EqualHandler = Dispatcher(
    "EqualHandler",
    doc="""
    Handler for Q.eq.
    Test that two expressions are equal.
    """
)

@EqualHandler.register(Basic, Basic)
def _(lhs, rhs, assumptions):
    return None

@EqualHandler.register(Tuple, Expr)
def _(lhs, rhs, assumptions):
    return False

@EqualHandler.register(Tuple, AppliedUndef)
def _(lhs, rhs, assumptions):
    return None

@EqualHandler.register(Tuple, Symbol)
def _(lhs, rhs, assumptions):
    return None

@EqualHandler.register(Tuple, Tuple)
def _(lhs, rhs, assumptions):
    if len(lhs) != len(rhs):
        return False
    return fuzzy_and(fuzzy_bool(ask(Q.eq(s, o), assumptions)) for s, o in zip(lhs, rhs))

@EqualHandler.register(BooleanAtom, BooleanAtom)
def _(lhs, rhs, assumptions):
    return lhs is rhs

@EqualHandler.register(Expr, Expr)
def _(lhs, rhs, assumptions):
    # 1. final check that lhs - rhs is zero
    try:
        diff = lhs - rhs
        if diff.kind is NumberKind:
            ret = ask(Q.zero(diff), assumptions)
            if ret is not None:
                return ret
    except TypeError:
        pass

    # 2. check infinite
    lhs_inf, rhs_inf = ask(Q.infinite(lhs), assumptions), ask(Q.infinite(rhs), assumptions)
    if lhs_inf or rhs_inf:
        # 2-1. check if both is infinite
        if fuzzy_xor([lhs_inf, rhs_inf]):
            return False
        lhs_exreal = ask(Q.extended_real(lhs), assumptions)
        rhs_exreal = ask(Q.extended_real(rhs), assumptions)
        if fuzzy_xor([lhs_exreal, rhs_exreal]):
            return False
        if fuzzy_and([lhs_exreal, rhs_exreal]):
            lhs_expos = ask(Q.positive(lhs, assumptions))
            rhs_expos = ask(Q.positive(rhs, assumptions))
            return fuzzy_xor([lhs_expos, fuzzy_not(rhs_expos)])

        # 2-2. Try to split real/imaginary parts and equate them
        I = S.ImaginaryUnit

        def split_real_imag(expr):
            real_imag = lambda t: (
                'real' if ask(Q.extended_real(t), assumptions) else
                'imag' if ask(Q.extended_real(I*t), assumptions) else None)
            return sift(Add.make_args(expr), real_imag)

        lhs_ri = split_real_imag(lhs)
        if not lhs_ri[None]:
            rhs_ri = split_real_imag(rhs)
            if not rhs_ri[None]:
                eq_real = ask(Q.eq(Add(*lhs_ri['real']), Add(*rhs_ri['real'])), assumptions)
                eq_imag = ask(Q.eq((I * Add(*lhs_ri['imag']), I * Add(*rhs_ri['imag']))), assumptions)
                return fuzzy_and(map(fuzzy_bool, [eq_real, eq_imag]))

        # Compare e.g. zoo with 1+I*oo by comparing args
        arglhs = arg(lhs)
        argrhs = arg(rhs)
        # Guard against Eq(nan, nan) -> Falsesymp
        if not (arglhs == S.NaN and argrhs == S.NaN):
            return fuzzy_bool(ask(Q.eq(arglhs, argrhs), assumptions))

    # 3. see if the ratio evaluates
    if all(isinstance(i, Expr) for i in (lhs, rhs)):
        dif = lhs - rhs
        n, d = dif.as_numer_denom()
        rv = None
        if ask(Q.zero(n), assumptions):
            rv = ask(Q.nonzero(d), assumptions)
        elif ask(Q.finite(n), assumptions):
            if  ask(Q.infinite(d), assumptions):
                rv = True
            elif ask(Q.zero(n), assumptions) is False:
                rv = ask(Q.infinite(d), assumptions)
                if rv is None:
                    # if the condition that makes the denominator
                    # infinite does not make the original expression
                    # True then False can be returned
                    l, r = clear_coefficients(d, S.Infinity)
                    args = [_.subs(l, r) for _ in (lhs, rhs)]
                    if args != [lhs, rhs]:
                        rv = fuzzy_bool(ask(Q.eq(*args), assumptions))
                        if rv is True:
                            rv = None
        elif any(ask(Q.infinite(a), assumptions) for a in Add.make_args(n)):
            # (inf or nan)/x != 0
            rv = False
        if rv is not None:
            return rv

    return None
