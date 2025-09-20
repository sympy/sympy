from sympy.assumptions.assume import recursive_ask
'\nHandlers for predicates related to set membership: integer, rational, etc.\n'
from sympy.assumptions import Q
from sympy.core import Add, Basic, Expr, Mul, Pow, S
from sympy.core.numbers import AlgebraicNumber, ComplexInfinity, Exp1, Float, GoldenRatio, ImaginaryUnit, Infinity, Integer, NaN, NegativeInfinity, Number, NumberSymbol, Pi, pi, Rational, TribonacciConstant, E
from sympy.core.logic import fuzzy_bool, fuzzy_not
from sympy.functions import Abs, acos, acot, asin, atan, cos, cot, exp, im, log, re, sin, tan
from sympy.core.numbers import I
from sympy.core.relational import Eq
from sympy.functions.elementary.complexes import conjugate
from sympy.matrices import Determinant, MatrixBase, Trace
from sympy.matrices.expressions.matexpr import MatrixElement
from sympy.multipledispatch import MDNotImplementedError
from .common import test_closed_group, ask_all, ask_any
from ..predicates.sets import IntegerPredicate, RationalPredicate, IrrationalPredicate, RealPredicate, ExtendedRealPredicate, HermitianPredicate, ComplexPredicate, ImaginaryPredicate, AntihermitianPredicate, AlgebraicPredicate, TranscendentalPredicate

def _IntegerPredicate_number(expr, assumptions):
    try:
        i = int(expr.round())
        if not (expr - i).equals(0):
            raise TypeError
        return True
    except TypeError:
        return False

@IntegerPredicate.register_many(int, Integer)
def _(expr, assumptions, rec):
    return True

@IntegerPredicate.register_many(Exp1, GoldenRatio, ImaginaryUnit, Infinity, NegativeInfinity, Pi, Rational, TribonacciConstant)
def _(expr, assumptions, rec):
    return False

@IntegerPredicate.register(Expr)
def _(expr, assumptions, rec):
    ret = expr.is_integer
    if ret is None:
        raise MDNotImplementedError
    return ret

@IntegerPredicate.register(Add)
def _(expr, assumptions, rec):
    """
    * Integer + Integer       -> Integer
    * Integer + !Integer      -> !Integer
    * !Integer + !Integer -> ?
    """
    if expr.is_number:
        return _IntegerPredicate_number(expr, assumptions)
    return test_closed_group(expr, assumptions, Q.integer, rec=rec)

@IntegerPredicate.register(Pow)
def _(expr, assumptions, rec):
    if expr.is_number:
        return _IntegerPredicate_number(expr, assumptions)
    if ask_all(~Q.zero(expr.base), Q.finite(expr.base), Q.zero(expr.exp), assumptions=assumptions, rec=rec):
        return True
    if ask_all(Q.integer(expr.base), Q.integer(expr.exp), assumptions=assumptions, rec=rec):
        if ask_any(Q.positive(expr.exp), Q.nonnegative(expr.exp) & ~Q.zero(expr.base), Q.zero(expr.base - 1), Q.zero(expr.base + 1), assumptions=assumptions, rec=rec):
            return True

@IntegerPredicate.register(Mul)
def _(expr, assumptions, rec):
    """
    * Integer*Integer      -> Integer
    * Integer*Irrational   -> !Integer
    * Odd/Even             -> !Integer
    * Integer*Rational     -> ?
    """
    if expr.is_number:
        return _IntegerPredicate_number(expr, assumptions)
    _output = True
    for arg in expr.args:
        if not recursive_ask(Q.integer(arg), assumptions=assumptions, rec=rec):
            if arg.is_Rational:
                if arg.q == 2:
                    return recursive_ask(Q.even(2 * expr), assumptions=assumptions, rec=rec)
                if ~(arg.q & 1):
                    return None
            elif recursive_ask(Q.irrational(arg), assumptions=assumptions, rec=rec):
                if _output:
                    _output = False
                else:
                    return
            else:
                return
    return _output

@IntegerPredicate.register(Abs)
def _(expr, assumptions, rec):
    if recursive_ask(Q.integer(expr.args[0]), assumptions=assumptions, rec=rec):
        return True

@IntegerPredicate.register_many(Determinant, MatrixElement, Trace)
def _(expr, assumptions, rec):
    return recursive_ask(Q.integer_elements(expr.args[0]), assumptions=assumptions, rec=rec)

@RationalPredicate.register(Rational)
def _(expr, assumptions, rec):
    return True

@RationalPredicate.register(Float)
def _(expr, assumptions, rec):
    return None

@RationalPredicate.register_many(Exp1, GoldenRatio, ImaginaryUnit, Infinity, NegativeInfinity, Pi, TribonacciConstant)
def _(expr, assumptions, rec):
    return False

@RationalPredicate.register(Expr)
def _(expr, assumptions, rec):
    ret = expr.is_rational
    if ret is None:
        raise MDNotImplementedError
    return ret

@RationalPredicate.register_many(Add, Mul)
def _(expr, assumptions, rec):
    """
    * Rational + Rational     -> Rational
    * Rational + !Rational    -> !Rational
    * !Rational + !Rational   -> ?
    """
    if expr.is_number:
        if expr.as_real_imag()[1]:
            return False
    return test_closed_group(expr, assumptions, Q.rational, rec=rec)

@RationalPredicate.register(Pow)
def _(expr, assumptions, rec):
    """
    * Rational ** Integer      -> Rational
    * Irrational ** Rational   -> Irrational
    * Rational ** Irrational   -> ?
    """
    if expr.base == E:
        x = expr.exp
        if recursive_ask(Q.rational(x), assumptions=assumptions, rec=rec):
            return recursive_ask(Q.zero(x), assumptions=assumptions, rec=rec)
        return
    is_exp_integer = recursive_ask(Q.integer(expr.exp), assumptions=assumptions, rec=rec)
    if is_exp_integer:
        is_base_rational = recursive_ask(Q.rational(expr.base), assumptions=assumptions, rec=rec)
        if is_base_rational:
            is_base_zero = recursive_ask(Q.zero(expr.base), assumptions=assumptions, rec=rec)
            if is_base_zero is False:
                return True
            if is_base_zero and recursive_ask(Q.positive(expr.exp), assumptions=assumptions, rec=rec):
                return True
        if recursive_ask(Q.algebraic(expr.base), assumptions=assumptions, rec=rec) is False:
            return recursive_ask(Q.zero(expr.exp), assumptions=assumptions, rec=rec)
        if recursive_ask(Q.irrational(expr.base), assumptions=assumptions, rec=rec) and recursive_ask(Q.zero(expr.exp - -1), assumptions=assumptions, rec=rec):
            return False
        return
    elif recursive_ask(Q.rational(expr.exp), assumptions=assumptions, rec=rec):
        if recursive_ask(Q.prime(expr.base), assumptions=assumptions, rec=rec) and is_exp_integer is False:
            return False
        if recursive_ask(Q.zero(expr.base), assumptions=assumptions, rec=rec) and recursive_ask(Q.positive(expr.exp), assumptions=assumptions, rec=rec):
            return True
        if recursive_ask(Q.zero(expr.base - 1), assumptions=assumptions, rec=rec):
            return True

@RationalPredicate.register_many(asin, atan, cos, sin, tan)
def _(expr, assumptions, rec):
    x = expr.args[0]
    if recursive_ask(Q.rational(x), assumptions=assumptions, rec=rec):
        return recursive_ask(~Q.nonzero(x), assumptions=assumptions, rec=rec)

@RationalPredicate.register(exp)
def _(expr, assumptions, rec):
    x = expr.exp
    if recursive_ask(Q.rational(x), assumptions=assumptions, rec=rec):
        return recursive_ask(~Q.nonzero(x), assumptions=assumptions, rec=rec)

@RationalPredicate.register_many(acot, cot)
def _(expr, assumptions, rec):
    x = expr.args[0]
    if recursive_ask(Q.rational(x), assumptions=assumptions, rec=rec):
        return False

@RationalPredicate.register_many(acos, log)
def _(expr, assumptions, rec):
    x = expr.args[0]
    if recursive_ask(Q.rational(x), assumptions=assumptions, rec=rec):
        return recursive_ask(~Q.nonzero(x - 1), assumptions=assumptions, rec=rec)

@IrrationalPredicate.register(Expr)
def _(expr, assumptions, rec):
    ret = expr.is_irrational
    if ret is None:
        raise MDNotImplementedError
    return ret

@IrrationalPredicate.register(Basic)
def _(expr, assumptions, rec):
    _real = recursive_ask(Q.real(expr), assumptions=assumptions, rec=rec)
    if _real:
        _rational = recursive_ask(Q.rational(expr), assumptions=assumptions, rec=rec)
        if _rational is None:
            return None
        return not _rational
    else:
        return _real

def _RealPredicate_number(expr, assumptions):
    i = expr.as_real_imag()[1].evalf(2)
    if i._prec != 1:
        return not i

@RealPredicate.register_many(Abs, Exp1, Float, GoldenRatio, im, Pi, Rational, re, TribonacciConstant)
def _(expr, assumptions, rec):
    return True

@RealPredicate.register_many(ImaginaryUnit, Infinity, NegativeInfinity)
def _(expr, assumptions, rec):
    return False

@RealPredicate.register(Expr)
def _(expr, assumptions, rec):
    ret = expr.is_real
    if ret is None:
        raise MDNotImplementedError
    return ret

@RealPredicate.register(Add)
def _(expr, assumptions, rec):
    """
    * Real + Real              -> Real
    * Real + (Complex & !Real) -> !Real
    """
    if expr.is_number:
        return _RealPredicate_number(expr, assumptions)
    return test_closed_group(expr, assumptions, Q.real, rec=rec)

@RealPredicate.register(Mul)
def _(expr, assumptions, rec):
    """
    * Real*Real               -> Real
    * Real*Imaginary          -> !Real
    * Imaginary*Imaginary     -> Real
    """
    if expr.is_number:
        return _RealPredicate_number(expr, assumptions)
    result = True
    for arg in expr.args:
        if recursive_ask(Q.real(arg), assumptions=assumptions, rec=rec):
            pass
        elif recursive_ask(Q.imaginary(arg), assumptions=assumptions, rec=rec):
            result = result ^ True
        else:
            break
    else:
        return result

@RealPredicate.register(Pow)
def _(expr, assumptions, rec):
    """
    * Real**Integer              -> Real
    * Positive**Real             -> Real
    * Negative**Real             -> ?
    * Real**(Integer/Even)       -> Real if base is nonnegative
    * Real**(Integer/Odd)        -> Real
    * Imaginary**(Integer/Even)  -> Real
    * Imaginary**(Integer/Odd)   -> not Real
    * Imaginary**Real            -> ? since Real could be 0 (giving real)
                                    or 1 (giving imaginary)
    * b**Imaginary               -> Real if log(b) is imaginary and b != 0
                                    and exponent != integer multiple of
                                    I*pi/log(b)
    * Real**Real                 -> ? e.g. sqrt(-1) is imaginary and
                                    sqrt(2) is not
    """
    if expr.is_number:
        return _RealPredicate_number(expr, assumptions)
    if expr.base == E:
        return recursive_ask(Q.integer(expr.exp / I / pi) | Q.real(expr.exp), assumptions=assumptions, rec=rec)
    if expr.base.func == exp or (expr.base.is_Pow and expr.base.base == E):
        if recursive_ask(Q.imaginary(expr.base.exp), assumptions=assumptions, rec=rec):
            if recursive_ask(Q.imaginary(expr.exp), assumptions=assumptions, rec=rec):
                return True
        i = expr.base.exp / I / pi
        if recursive_ask(Q.integer(2 * i), assumptions=assumptions, rec=rec):
            return recursive_ask(Q.real((S.NegativeOne ** i) ** expr.exp), assumptions=assumptions, rec=rec)
        return
    if recursive_ask(Q.imaginary(expr.base), assumptions=assumptions, rec=rec):
        if recursive_ask(Q.integer(expr.exp), assumptions=assumptions, rec=rec):
            odd = recursive_ask(Q.odd(expr.exp), assumptions=assumptions, rec=rec)
            if odd is not None:
                return not odd
            return
    if recursive_ask(Q.imaginary(expr.exp), assumptions=assumptions, rec=rec):
        imlog = recursive_ask(Q.imaginary(log(expr.base)), assumptions=assumptions, rec=rec)
        if imlog is not None:
            return imlog
    if recursive_ask(Q.real(expr.base), assumptions=assumptions, rec=rec):
        if recursive_ask(Q.real(expr.exp), assumptions=assumptions, rec=rec):
            if recursive_ask(Q.zero(expr.base), assumptions=assumptions, rec=rec) is not False:
                if recursive_ask(Q.positive(expr.exp), assumptions=assumptions, rec=rec):
                    return True
                return
            if expr.exp.is_Rational and recursive_ask(Q.even(expr.exp.q), assumptions=assumptions, rec=rec):
                return recursive_ask(Q.positive(expr.base), assumptions=assumptions, rec=rec)
            elif recursive_ask(Q.integer(expr.exp), assumptions=assumptions, rec=rec):
                return True
            elif recursive_ask(Q.positive(expr.base), assumptions=assumptions, rec=rec):
                return True

@RealPredicate.register_many(cos, sin)
def _(expr, assumptions, rec):
    if recursive_ask(Q.real(expr.args[0]), assumptions=assumptions, rec=rec):
        return True

@RealPredicate.register(exp)
def _(expr, assumptions, rec):
    return recursive_ask(Q.integer(expr.exp / I / pi) | Q.real(expr.exp), assumptions=assumptions, rec=rec)

@RealPredicate.register(log)
def _(expr, assumptions, rec):
    return recursive_ask(Q.positive(expr.args[0]), assumptions=assumptions, rec=rec)

@RealPredicate.register_many(Determinant, MatrixElement, Trace)
def _(expr, assumptions, rec):
    return recursive_ask(Q.real_elements(expr.args[0]), assumptions=assumptions, rec=rec)

@ExtendedRealPredicate.register(object)
def _(expr, assumptions, rec):
    return recursive_ask(Q.negative_infinite(expr) | Q.negative(expr) | Q.zero(expr) | Q.positive(expr) | Q.positive_infinite(expr), assumptions=assumptions, rec=rec)

@ExtendedRealPredicate.register_many(Infinity, NegativeInfinity)
def _(expr, assumptions, rec):
    return True

@ExtendedRealPredicate.register_many(Add, Mul, Pow)
def _(expr, assumptions, rec):
    return test_closed_group(expr, assumptions, Q.extended_real, rec=rec)

@HermitianPredicate.register(object)
def _(expr, assumptions, rec):
    if isinstance(expr, MatrixBase):
        return None
    return recursive_ask(Q.real(expr), assumptions=assumptions, rec=rec)

@HermitianPredicate.register(Add)
def _(expr, assumptions, rec):
    """
    * Hermitian + Hermitian  -> Hermitian
    * Hermitian + !Hermitian -> !Hermitian
    """
    if expr.is_number:
        raise MDNotImplementedError
    return test_closed_group(expr, assumptions, Q.hermitian, rec=rec)

@HermitianPredicate.register(Mul)
def _(expr, assumptions, rec):
    """
    As long as there is at most only one noncommutative term:

    * Hermitian*Hermitian         -> Hermitian
    * Hermitian*Antihermitian     -> !Hermitian
    * Antihermitian*Antihermitian -> Hermitian
    """
    if expr.is_number:
        raise MDNotImplementedError
    nccount = 0
    result = True
    for arg in expr.args:
        if recursive_ask(Q.antihermitian(arg), assumptions=assumptions, rec=rec):
            result = result ^ True
        elif not recursive_ask(Q.hermitian(arg), assumptions=assumptions, rec=rec):
            break
        if recursive_ask(~Q.commutative(arg), assumptions=assumptions, rec=rec):
            nccount += 1
            if nccount > 1:
                break
    else:
        return result

@HermitianPredicate.register(Pow)
def _(expr, assumptions, rec):
    """
    * Hermitian**Integer -> Hermitian
    """
    if expr.is_number:
        raise MDNotImplementedError
    if expr.base == E:
        if recursive_ask(Q.hermitian(expr.exp), assumptions=assumptions, rec=rec):
            return True
        raise MDNotImplementedError
    if recursive_ask(Q.hermitian(expr.base), assumptions=assumptions, rec=rec):
        if recursive_ask(Q.integer(expr.exp), assumptions=assumptions, rec=rec):
            return True
    raise MDNotImplementedError

@HermitianPredicate.register_many(cos, sin)
def _(expr, assumptions, rec):
    if recursive_ask(Q.hermitian(expr.args[0]), assumptions=assumptions, rec=rec):
        return True
    raise MDNotImplementedError

@HermitianPredicate.register(exp)
def _(expr, assumptions, rec):
    if recursive_ask(Q.hermitian(expr.exp), assumptions=assumptions, rec=rec):
        return True
    raise MDNotImplementedError

@HermitianPredicate.register(MatrixBase)
def _(mat, assumptions, rec):
    rows, cols = mat.shape
    ret_val = True
    for i in range(rows):
        for j in range(i, cols):
            cond = recursive_ask(Q.zero(mat[i, j] - conjugate(mat[j, i])), assumptions=assumptions, rec=rec)
            if cond is None:
                ret_val = None
            if cond == False:
                return False
    if ret_val is None:
        raise MDNotImplementedError
    return ret_val

@ComplexPredicate.register_many(Abs, cos, exp, im, ImaginaryUnit, log, Number, NumberSymbol, re, sin)
def _(expr, assumptions, rec):
    return True

@ComplexPredicate.register_many(Infinity, NegativeInfinity)
def _(expr, assumptions, rec):
    return False

@ComplexPredicate.register(Expr)
def _(expr, assumptions, rec):
    ret = expr.is_complex
    if ret is None:
        raise MDNotImplementedError
    return ret

@ComplexPredicate.register_many(Add, Mul)
def _(expr, assumptions, rec):
    return test_closed_group(expr, assumptions, Q.complex, rec=rec)

@ComplexPredicate.register(Pow)
def _(expr, assumptions, rec):
    if expr.base == E:
        return True
    return test_closed_group(expr, assumptions, Q.complex, rec=rec)

@ComplexPredicate.register_many(Determinant, MatrixElement, Trace)
def _(expr, assumptions, rec):
    return recursive_ask(Q.complex_elements(expr.args[0]), assumptions=assumptions, rec=rec)

@ComplexPredicate.register(NaN)
def _(expr, assumptions, rec):
    return None

def _Imaginary_number(expr, assumptions):
    r = expr.as_real_imag()[0].evalf(2)
    if r._prec != 1:
        return not r

@ImaginaryPredicate.register(ImaginaryUnit)
def _(expr, assumptions, rec):
    return True

@ImaginaryPredicate.register(Expr)
def _(expr, assumptions, rec):
    ret = expr.is_imaginary
    if ret is None:
        raise MDNotImplementedError
    return ret

@ImaginaryPredicate.register(Add)
def _(expr, assumptions, rec):
    """
    * Imaginary + Imaginary -> Imaginary
    * Imaginary + Complex   -> ?
    * Imaginary + Real      -> !Imaginary
    """
    if expr.is_number:
        return _Imaginary_number(expr, assumptions)
    reals = 0
    for arg in expr.args:
        if recursive_ask(Q.imaginary(arg), assumptions=assumptions, rec=rec):
            pass
        elif recursive_ask(Q.real(arg), assumptions=assumptions, rec=rec):
            reals += 1
        else:
            break
    else:
        if reals == 0:
            return True
        if reals in (1, len(expr.args)):
            return False

@ImaginaryPredicate.register(Mul)
def _(expr, assumptions, rec):
    """
    * Real*Imaginary      -> Imaginary
    * Imaginary*Imaginary -> Real
    """
    if expr.is_number:
        return _Imaginary_number(expr, assumptions)
    result = False
    reals = 0
    for arg in expr.args:
        if recursive_ask(Q.imaginary(arg), assumptions=assumptions, rec=rec):
            result = result ^ True
        elif not recursive_ask(Q.real(arg), assumptions=assumptions, rec=rec):
            break
    else:
        if reals == len(expr.args):
            return False
        return result

@ImaginaryPredicate.register(Pow)
def _(expr, assumptions, rec):
    """
    * Imaginary**Odd        -> Imaginary
    * Imaginary**Even       -> Real
    * b**Imaginary          -> !Imaginary if exponent is an integer
                               multiple of I*pi/log(b)
    * Imaginary**Real       -> ?
    * Positive**Real        -> Real
    * Negative**Integer     -> Real
    * Negative**(Integer/2) -> Imaginary
    * Negative**Real        -> not Imaginary if exponent is not Rational
    """
    if expr.is_number:
        return _Imaginary_number(expr, assumptions)
    if expr.base == E:
        a = expr.exp / I / pi
        return recursive_ask(Q.integer(2 * a) & ~Q.integer(a), assumptions=assumptions, rec=rec)
    if expr.base.func == exp or (expr.base.is_Pow and expr.base.base == E):
        if recursive_ask(Q.imaginary(expr.base.exp), assumptions=assumptions, rec=rec):
            if recursive_ask(Q.imaginary(expr.exp), assumptions=assumptions, rec=rec):
                return False
            i = expr.base.exp / I / pi
            if recursive_ask(Q.integer(2 * i), assumptions=assumptions, rec=rec):
                return recursive_ask(Q.imaginary((S.NegativeOne ** i) ** expr.exp), assumptions=assumptions, rec=rec)
    if recursive_ask(Q.imaginary(expr.base), assumptions=assumptions, rec=rec):
        if recursive_ask(Q.integer(expr.exp), assumptions=assumptions, rec=rec):
            odd = recursive_ask(Q.odd(expr.exp), assumptions=assumptions, rec=rec)
            if odd is not None:
                return odd
            return
    if recursive_ask(Q.imaginary(expr.exp), assumptions=assumptions, rec=rec):
        imlog = recursive_ask(Q.imaginary(log(expr.base)), assumptions=assumptions, rec=rec)
        if imlog is not None:
            return False
    if recursive_ask(Q.real(expr.base) & Q.real(expr.exp), assumptions=assumptions, rec=rec):
        if recursive_ask(Q.positive(expr.base), assumptions=assumptions, rec=rec):
            return False
        else:
            rat = recursive_ask(Q.rational(expr.exp), assumptions=assumptions, rec=rec)
            if not rat:
                return rat
            if recursive_ask(Q.integer(expr.exp), assumptions=assumptions, rec=rec):
                return False
            else:
                half = recursive_ask(Q.integer(2 * expr.exp), assumptions=assumptions, rec=rec)
                if half:
                    return recursive_ask(Q.negative(expr.base), assumptions=assumptions, rec=rec)
                return half

@ImaginaryPredicate.register(log)
def _(expr, assumptions, rec):
    if recursive_ask(Q.real(expr.args[0]), assumptions=assumptions, rec=rec):
        if recursive_ask(Q.positive(expr.args[0]), assumptions=assumptions, rec=rec):
            return False
        return
    if expr.args[0].func == exp or (expr.args[0].is_Pow and expr.args[0].base == E):
        if expr.args[0].exp in [I, -I]:
            return True
    im = recursive_ask(Q.imaginary(expr.args[0]), assumptions=assumptions, rec=rec)
    if im is False:
        return False

@ImaginaryPredicate.register(exp)
def _(expr, assumptions, rec):
    a = expr.exp / I / pi
    return recursive_ask(Q.integer(2 * a) & ~Q.integer(a), assumptions=assumptions, rec=rec)

@ImaginaryPredicate.register_many(Number, NumberSymbol)
def _(expr, assumptions, rec):
    return not expr.as_real_imag()[1] == 0

@ImaginaryPredicate.register(NaN)
def _(expr, assumptions, rec):
    return None

@AntihermitianPredicate.register(object)
def _(expr, assumptions, rec):
    if isinstance(expr, MatrixBase):
        return None
    if recursive_ask(Q.zero(expr), assumptions=assumptions, rec=rec):
        return True
    return recursive_ask(Q.imaginary(expr), assumptions=assumptions, rec=rec)

@AntihermitianPredicate.register(Add)
def _(expr, assumptions, rec):
    """
    * Antihermitian + Antihermitian  -> Antihermitian
    * Antihermitian + !Antihermitian -> !Antihermitian
    """
    if expr.is_number:
        raise MDNotImplementedError
    return test_closed_group(expr, assumptions, Q.antihermitian, rec=rec)

@AntihermitianPredicate.register(Mul)
def _(expr, assumptions, rec):
    """
    As long as there is at most only one noncommutative term:

    * Hermitian*Hermitian         -> !Antihermitian
    * Hermitian*Antihermitian     -> Antihermitian
    * Antihermitian*Antihermitian -> !Antihermitian
    """
    if expr.is_number:
        raise MDNotImplementedError
    nccount = 0
    result = False
    for arg in expr.args:
        if recursive_ask(Q.antihermitian(arg), assumptions=assumptions, rec=rec):
            result = result ^ True
        elif not recursive_ask(Q.hermitian(arg), assumptions=assumptions, rec=rec):
            break
        if recursive_ask(~Q.commutative(arg), assumptions=assumptions, rec=rec):
            nccount += 1
            if nccount > 1:
                break
    else:
        return result

@AntihermitianPredicate.register(Pow)
def _(expr, assumptions, rec):
    """
    * Hermitian**Integer  -> !Antihermitian
    * Antihermitian**Even -> !Antihermitian
    * Antihermitian**Odd  -> Antihermitian
    """
    if expr.is_number:
        raise MDNotImplementedError
    if recursive_ask(Q.hermitian(expr.base), assumptions=assumptions, rec=rec):
        if recursive_ask(Q.integer(expr.exp), assumptions=assumptions, rec=rec):
            return False
    elif recursive_ask(Q.antihermitian(expr.base), assumptions=assumptions, rec=rec):
        if recursive_ask(Q.even(expr.exp), assumptions=assumptions, rec=rec):
            return False
        elif recursive_ask(Q.odd(expr.exp), assumptions=assumptions, rec=rec):
            return True
    raise MDNotImplementedError

@AntihermitianPredicate.register(MatrixBase)
def _(mat, assumptions, rec):
    rows, cols = mat.shape
    ret_val = True
    for i in range(rows):
        for j in range(i, cols):
            cond = recursive_ask(Q.zero(mat[i, j] - -conjugate(mat[j, i])), assumptions=assumptions, rec=rec)
            if cond is None:
                ret_val = None
            if cond == False:
                return False
    if ret_val is None:
        raise MDNotImplementedError
    return ret_val

@AlgebraicPredicate.register_many(AlgebraicNumber, GoldenRatio, ImaginaryUnit, TribonacciConstant)
def _(expr, assumptions, rec):
    return True

@AlgebraicPredicate.register(Float)
def _(expr, assumptions, rec):
    return None

@AlgebraicPredicate.register_many(ComplexInfinity, Exp1, Infinity, NegativeInfinity, Pi)
def _(expr, assumptions, rec):
    return False

@AlgebraicPredicate.register_many(Add, Mul)
def _(expr, assumptions, rec):
    return test_closed_group(expr, assumptions, Q.algebraic, rec=rec)

@AlgebraicPredicate.register(Pow)
def _(expr, assumptions, rec):
    if expr.base == E:
        if recursive_ask(Q.algebraic(expr.exp), assumptions=assumptions, rec=rec):
            return recursive_ask(~Q.nonzero(expr.exp), assumptions=assumptions, rec=rec)
        return
    if expr.base == pi:
        if recursive_ask(Q.integer(expr.exp), assumptions=assumptions, rec=rec) and recursive_ask(Q.positive(expr.exp), assumptions=assumptions, rec=rec):
            return False
        return
    exp_rational = recursive_ask(Q.rational(expr.exp), assumptions=assumptions, rec=rec)
    base_algebraic = recursive_ask(Q.algebraic(expr.base), assumptions=assumptions, rec=rec)
    exp_algebraic = recursive_ask(Q.algebraic(expr.exp), assumptions=assumptions, rec=rec)
    if base_algebraic and exp_algebraic:
        if exp_rational:
            return True
        if recursive_ask(~Q.zero(expr.base - 0) & ~Q.zero(expr.base - 1), assumptions=assumptions, rec=rec) and exp_rational is False:
            return False
    exp_integer = recursive_ask(Q.integer(expr.exp), assumptions=assumptions, rec=rec)
    if base_algebraic is False and exp_integer:
        if expr.exp > 0:
            return False

@AlgebraicPredicate.register(Rational)
def _(expr, assumptions, rec):
    return expr.q != 0

@AlgebraicPredicate.register_many(asin, atan, cos, sin, tan)
def _(expr, assumptions, rec):
    x = expr.args[0]
    if recursive_ask(Q.algebraic(x), assumptions=assumptions, rec=rec):
        return recursive_ask(~Q.nonzero(x), assumptions=assumptions, rec=rec)

@AlgebraicPredicate.register(exp)
def _(expr, assumptions, rec):
    x = expr.exp
    if recursive_ask(Q.algebraic(x), assumptions=assumptions, rec=rec):
        return recursive_ask(~Q.nonzero(x), assumptions=assumptions, rec=rec)

@AlgebraicPredicate.register_many(acot, cot)
def _(expr, assumptions, rec):
    x = expr.args[0]
    if recursive_ask(Q.algebraic(x), assumptions=assumptions, rec=rec):
        return False

@AlgebraicPredicate.register_many(acos, log)
def _(expr, assumptions, rec):
    x = expr.args[0]
    if recursive_ask(Q.algebraic(x), assumptions=assumptions, rec=rec):
        return recursive_ask(~Q.nonzero(x - 1), assumptions=assumptions, rec=rec)

@TranscendentalPredicate.register(Expr)
def _(expr, assumptions, rec):
    ret = expr.is_transcendental
    if ret is not None:
        return ret
    is_complex = recursive_ask(Q.complex(expr), assumptions=assumptions, rec=rec)
    if is_complex:
        is_algebraic = recursive_ask(Q.algebraic(expr), assumptions=assumptions, rec=rec)
        return fuzzy_not(is_algebraic)
    return is_complex