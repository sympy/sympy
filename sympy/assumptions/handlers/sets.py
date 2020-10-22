"""
Handlers for predicates related to set membership: integer, rational, etc.
"""

from sympy.assumptions import Q, ask
from sympy.assumptions.handlers import CommonHandler, test_closed_group
from sympy.core import Basic, Expr, Add, Mul, Pow
from sympy.core.numbers import (
    pi, Integer, Pi, Exp1, GoldenRatio, TribonacciConstant, Infinity,
    NegativeInfinity, ImaginaryUnit, Rational, Float, Number, NumberSymbol,
    ComplexInfinity, AlgebraicNumber
)
from sympy.core.logic import fuzzy_bool
from sympy.functions import (
    exp, log, Abs, cot, sin, cos, tan, asin, atan, acos, acot,
    re, im
)
from sympy.matrices import (
    MatrixBase, Determinant, Trace
)
from sympy.matrices.expressions.matexpr import MatrixElement
from sympy import I, Eq, conjugate


### AskIntegerHandler ###

AskIntegerHandler = CommonHandler.copy(
    'AskIntegerHandler',
    doc="""
    Handler for Q.integer
    Test that an expression belongs to the field of integer numbers
    """
)

for sig in (int, Integer):
    AskIntegerHandler.register(sig)(AskIntegerHandler.AlwaysTrue)

for sig in (
Pi, Exp1, GoldenRatio, TribonacciConstant, Infinity, NegativeInfinity, ImaginaryUnit
):
    AskIntegerHandler.register(sig)(AskIntegerHandler.AlwaysFalse)

def _number_integer(expr, assumptions):
    # helper method
    try:
        i = int(expr.round())
        if not (expr - i).equals(0):
            raise TypeError
        return True
    except TypeError:
        return False

@AskIntegerHandler.register(Expr)
def _(expr, assumptions):
    return expr.is_integer

@AskIntegerHandler.register(Pow)
@AskIntegerHandler.register(Add)
def _(expr, assumptions):
    """
    Integer + Integer       -> Integer
    Integer + !Integer      -> !Integer
    !Integer + !Integer -> ?
    """
    if expr.is_number:
        return _number_integer(expr, assumptions)
    return test_closed_group(expr, assumptions, Q.integer)

@AskIntegerHandler.register(Mul)
def _(expr, assumptions):
    """
    Integer*Integer      -> Integer
    Integer*Irrational   -> !Integer
    Odd/Even             -> !Integer
    Integer*Rational     -> ?
    """
    if expr.is_number:
        return _number_integer(expr, assumptions)
    _output = True
    for arg in expr.args:
        if not ask(Q.integer(arg), assumptions):
            if arg.is_Rational:
                if arg.q == 2:
                    return ask(Q.even(2*expr), assumptions)
                if ~(arg.q & 1):
                    return None
            elif ask(Q.irrational(arg), assumptions):
                if _output:
                    _output = False
                else:
                    return
            else:
                return

    return _output

@AskIntegerHandler.register(Rational)
def _(expr, assumptions):
    # rationals with denominator one get
    # evaluated to Integers
    return False

@AskIntegerHandler.register(Abs)
def _(expr, assumptions):
    return ask(Q.integer(expr.args[0]), assumptions)

@AskIntegerHandler.register(MatrixBase)
@AskIntegerHandler.register(Determinant)
@AskIntegerHandler.register(Trace)
@AskIntegerHandler.register(MatrixElement)
def _(expr, assumptions):
    return ask(Q.integer_elements(expr.args[0]), assumptions)


### AskRationalHandler ###

AskRationalHandler = CommonHandler.copy(
    'AskRationalHandler',
    doc="""
    Handler for Q.rational
    Test that an expression belongs to the field of rational numbers
    """
)

AskRationalHandler.register(Rational)(AskRationalHandler.AlwaysTrue)

AskRationalHandler.register(Float)(AskRationalHandler.AlwaysNone)

for sig in (
ImaginaryUnit, Infinity, NegativeInfinity, Pi, Exp1, GoldenRatio, TribonacciConstant
):
    AskRationalHandler.register(sig)(AskRationalHandler.AlwaysFalse)

@AskRationalHandler.register(Expr)
def _(expr, assumptions):
    return expr.is_rational

@AskRationalHandler.register(Add)
@AskRationalHandler.register(Mul)
def _(expr, assumptions):
    """
    Rational + Rational     -> Rational
    Rational + !Rational    -> !Rational
    !Rational + !Rational   -> ?
    """
    if expr.is_number:
        if expr.as_real_imag()[1]:
            return False
    return test_closed_group(expr, assumptions, Q.rational)

@AskRationalHandler.register(Pow)
def _(expr, assumptions):
    """
    Rational ** Integer      -> Rational
    Irrational ** Rational   -> Irrational
    Rational ** Irrational   -> ?
    """
    if ask(Q.integer(expr.exp), assumptions):
        return ask(Q.rational(expr.base), assumptions)
    elif ask(Q.rational(expr.exp), assumptions):
        if ask(Q.prime(expr.base), assumptions):
            return False

@AskRationalHandler.register(exp)
@AskRationalHandler.register(sin)
@AskRationalHandler.register(cos)
@AskRationalHandler.register(tan)
@AskRationalHandler.register(asin)
@AskRationalHandler.register(atan)
def _(expr, assumptions):
    x = expr.args[0]
    if ask(Q.rational(x), assumptions):
        return ask(~Q.nonzero(x), assumptions)

@AskRationalHandler.register(cot)
@AskRationalHandler.register(acot)
def _(expr, assumptions):
    x = expr.args[0]
    if ask(Q.rational(x), assumptions):
        return False

@AskRationalHandler.register(log)
@AskRationalHandler.register(acos)
def _(expr, assumptions):
    x = expr.args[0]
    if ask(Q.rational(x), assumptions):
        return ask(~Q.nonzero(x - 1), assumptions)


### AskIrrationalHandler ###

AskIrrationalHandler = CommonHandler.copy(
    'AskIrrationalHandler',
    doc="""
    Handler for Q.irrational
    """
)

@AskIrrationalHandler.register(Expr)
def _(expr, assumptions):
    return expr.is_irrational

@AskIrrationalHandler.register(Basic)
def _(expr, assumptions):
    _real = ask(Q.real(expr), assumptions)
    if _real:
        _rational = ask(Q.rational(expr), assumptions)
        if _rational is None:
            return None
        return not _rational
    else:
        return _real


### AskRealHandler ###

AskRealHandler = CommonHandler.copy(
    'AskRealHandler',
    doc="""
    Handler for Q.real
    Test that an expression belongs to the field of real numbers
    """
)

for sig in (
Rational, Float, Pi, Exp1, GoldenRatio, TribonacciConstant, Abs, re, im
):
    AskRealHandler.register(sig)(AskRealHandler.AlwaysTrue)

for sig in (ImaginaryUnit, Infinity, NegativeInfinity):
    AskRealHandler.register(sig)(AskRealHandler.AlwaysFalse)

def _number_real(expr, assumptions):
    # let as_real_imag() work first since the expression may
    # be simpler to evaluate
    i = expr.as_real_imag()[1].evalf(2)
    if i._prec != 1:
        return not i
    # allow None to be returned if we couldn't show for sure
    # that i was 0

@AskRealHandler.register(Expr)
def _(expr, assumptions):
    return expr.is_real

@AskRealHandler.register(Add)
def _(expr, assumptions):
    """
    Real + Real              -> Real
    Real + (Complex & !Real) -> !Real
    """
    if expr.is_number:
        return _number_real(expr, assumptions)
    return test_closed_group(expr, assumptions, Q.real)

@AskRealHandler.register(Mul)
def _(expr, assumptions):
    """
    Real*Real               -> Real
    Real*Imaginary          -> !Real
    Imaginary*Imaginary     -> Real
    """
    if expr.is_number:
        return _number_real(expr, assumptions)
    result = True
    for arg in expr.args:
        if ask(Q.real(arg), assumptions):
            pass
        elif ask(Q.imaginary(arg), assumptions):
            result = result ^ True
        else:
            break
    else:
        return result

@AskRealHandler.register(Pow)
def _(expr, assumptions):
    """
    Real**Integer              -> Real
    Positive**Real             -> Real
    Real**(Integer/Even)       -> Real if base is nonnegative
    Real**(Integer/Odd)        -> Real
    Imaginary**(Integer/Even)  -> Real
    Imaginary**(Integer/Odd)   -> not Real
    Imaginary**Real            -> ? since Real could be 0 (giving real) or 1 (giving imaginary)
    b**Imaginary               -> Real if log(b) is imaginary and b != 0 and exponent != integer multiple of I*pi/log(b)
    Real**Real                 -> ? e.g. sqrt(-1) is imaginary and sqrt(2) is not
    """
    if expr.is_number:
        return _number_real(expr, assumptions)

    if expr.base.func == exp:
        if ask(Q.imaginary(expr.base.args[0]), assumptions):
            if ask(Q.imaginary(expr.exp), assumptions):
                return True
        # If the i = (exp's arg)/(I*pi) is an integer or half-integer
        # multiple of I*pi then 2*i will be an integer. In addition,
        # exp(i*I*pi) = (-1)**i so the overall realness of the expr
        # can be determined by replacing exp(i*I*pi) with (-1)**i.
        i = expr.base.args[0]/I/pi
        if ask(Q.integer(2*i), assumptions):
            return ask(Q.real(((-1)**i)**expr.exp), assumptions)
        return

    if ask(Q.imaginary(expr.base), assumptions):
        if ask(Q.integer(expr.exp), assumptions):
            odd = ask(Q.odd(expr.exp), assumptions)
            if odd is not None:
                return not odd
            return

    if ask(Q.imaginary(expr.exp), assumptions):
        imlog = ask(Q.imaginary(log(expr.base)), assumptions)
        if imlog is not None:
            # I**i -> real, log(I) is imag;
            # (2*I)**i -> complex, log(2*I) is not imag
            return imlog

    if ask(Q.real(expr.base), assumptions):
        if ask(Q.real(expr.exp), assumptions):
            if expr.exp.is_Rational and \
                    ask(Q.even(expr.exp.q), assumptions):
                return ask(Q.positive(expr.base), assumptions)
            elif ask(Q.integer(expr.exp), assumptions):
                return True
            elif ask(Q.positive(expr.base), assumptions):
                return True
            elif ask(Q.negative(expr.base), assumptions):
                return False

@AskRealHandler.register(sin)
@AskRealHandler.register(cos)
def _(expr, assumptions):
    if ask(Q.real(expr.args[0]), assumptions):
        return True

@AskRealHandler.register(exp)
def _(expr, assumptions):
    return ask(Q.integer(expr.args[0]/I/pi) | Q.real(expr.args[0]), assumptions)

@AskRealHandler.register(log)
def _(expr, assumptions):
    return ask(Q.positive(expr.args[0]), assumptions)

@AskRealHandler.register(MatrixElement)
@AskRealHandler.register(Determinant)
@AskRealHandler.register(Trace)
def _(expr, assumptions):
    return ask(Q.real_elements(expr.args[0]), assumptions)


### AskExtendedRealHandler ###

AskExtendedRealHandler = AskRealHandler.copy(
    'AskExtendedRealHandler',
    doc="""
    Handler for Q.extended_real
    Test that an expression belongs to the field of extended real numbers,
    that is real numbers union {Infinity, -Infinity}
    """
)

for sig in (Infinity, NegativeInfinity):
    AskExtendedRealHandler.register(sig)(AskExtendedRealHandler.AlwaysTrue)

@AskExtendedRealHandler.register(Add)
@AskExtendedRealHandler.register(Mul)
@AskExtendedRealHandler.register(Pow)
def _(expr, assumptions):
    return test_closed_group(expr, assumptions, Q.extended_real)


### AskHermitianHandler ###

AskHermitianHandler = AskRealHandler.copy(
    'AskHermitianHandler',
    doc="""
    Handler for Q.hermitian
    Test that an expression belongs to the field of Hermitian operators
    """
)

@AskHermitianHandler.register(Expr)
def _(expr, assumptions):
    if isinstance(expr, MatrixBase):
        return None

@AskHermitianHandler.register(Add)
def _(expr, assumptions):
    """
    Hermitian + Hermitian  -> Hermitian
    Hermitian + !Hermitian -> !Hermitian
    """
    if expr.is_number:
        return _number_real(expr, assumptions)
    return test_closed_group(expr, assumptions, Q.hermitian)

@AskHermitianHandler.register(Mul)
def _(expr, assumptions):
    """
    As long as there is at most only one noncommutative term:
    Hermitian*Hermitian         -> Hermitian
    Hermitian*Antihermitian     -> !Hermitian
    Antihermitian*Antihermitian -> Hermitian
    """
    if expr.is_number:
        return _number_real(expr, assumptions)
    nccount = 0
    result = True
    for arg in expr.args:
        if ask(Q.antihermitian(arg), assumptions):
            result = result ^ True
        elif not ask(Q.hermitian(arg), assumptions):
            break
        if ask(~Q.commutative(arg), assumptions):
            nccount += 1
            if nccount > 1:
                break
    else:
        return result

@AskHermitianHandler.register(Pow)
def _(expr, assumptions):
    """
    Hermitian**Integer -> Hermitian
    """
    if expr.is_number:
        return _number_real(expr, assumptions)
    if ask(Q.hermitian(expr.base), assumptions):
        if ask(Q.integer(expr.exp), assumptions):
            return True

@AskHermitianHandler.register(sin)
@AskHermitianHandler.register(cos)
@AskHermitianHandler.register(exp)
def _(expr, assumptions):
    if ask(Q.hermitian(expr.args[0]), assumptions):
        return True

@AskHermitianHandler.register(MatrixBase)
def _(mat, assumptions):
    rows, cols = mat.shape
    ret_val = True
    for i in range(rows):
        for j in range(i, cols):
            cond = fuzzy_bool(Eq(mat[i, j], conjugate(mat[j, i])))
            if cond == None:
                ret_val = None
            if cond == False:
                return False
    return ret_val


### AskComplexHandler ###

AskComplexHandler = CommonHandler.copy(
    'AskComplexHandler',
    doc="""
    Handler for Q.complex
    Test that an expression belongs to the field of complex numbers
    """
)

for sig in (
Number, sin, cos, log, exp, re, im, NumberSymbol, Abs, ImaginaryUnit
):
    # they are all complex functions or expressions
    AskComplexHandler.register(sig)(AskComplexHandler.AlwaysTrue)

for sig in (Infinity, NegativeInfinity):
    AskComplexHandler.register(sig)(AskComplexHandler.AlwaysFalse)

@AskComplexHandler.register(Expr)
def _(expr, assumptions):
    return expr.is_complex

@AskComplexHandler.register(Add)
@AskComplexHandler.register(Mul)
@AskComplexHandler.register(Pow)
def _(expr, assumptions):
    return test_closed_group(expr, assumptions, Q.complex)

@AskComplexHandler.register(MatrixElement)
@AskComplexHandler.register(Determinant)
@AskComplexHandler.register(Trace)
def _(expr, assumptions):
    return ask(Q.complex_elements(expr.args[0]), assumptions)


### AskImaginaryHandler ###

AskImaginaryHandler = CommonHandler.copy(
    'AskImaginaryHandler',
    doc="""
    Handler for Q.imaginary
    Test that an expression belongs to the field of imaginary numbers,
    that is, numbers in the form x*I, where x is real
    """
)

AskImaginaryHandler.register(ImaginaryUnit)(AskImaginaryHandler.AlwaysTrue)

def _number_imaginary(expr, assumptions):
    # let as_real_imag() work first since the expression may
    # be simpler to evaluate
    r = expr.as_real_imag()[0].evalf(2)
    if r._prec != 1:
        return not r
    # allow None to be returned if we couldn't show for sure
    # that r was 0

@AskImaginaryHandler.register(Expr)
def _(expr, assumptions):
    return expr.is_imaginary

@AskImaginaryHandler.register(Add)
def _(expr, assumptions):
    """
    Imaginary + Imaginary -> Imaginary
    Imaginary + Complex   -> ?
    Imaginary + Real      -> !Imaginary
    """
    if expr.is_number:
        return _number_imaginary(expr, assumptions)

    reals = 0
    for arg in expr.args:
        if ask(Q.imaginary(arg), assumptions):
            pass
        elif ask(Q.real(arg), assumptions):
            reals += 1
        else:
            break
    else:
        if reals == 0:
            return True
        if reals == 1 or (len(expr.args) == reals):
            # two reals could sum 0 thus giving an imaginary
            return False

@AskImaginaryHandler.register(Mul)
def _(expr, assumptions):
    """
    Real*Imaginary      -> Imaginary
    Imaginary*Imaginary -> Real
    """
    if expr.is_number:
        return _number_imaginary(expr, assumptions)
    result = False
    reals = 0
    for arg in expr.args:
        if ask(Q.imaginary(arg), assumptions):
            result = result ^ True
        elif not ask(Q.real(arg), assumptions):
            break
    else:
        if reals == len(expr.args):
            return False
        return result

@AskImaginaryHandler.register(Pow)
def _(expr, assumptions):
    """
    Imaginary**Odd        -> Imaginary
    Imaginary**Even       -> Real
    b**Imaginary          -> !Imaginary if exponent is an integer multiple of I*pi/log(b)
    Imaginary**Real       -> ?
    Positive**Real        -> Real
    Negative**Integer     -> Real
    Negative**(Integer/2) -> Imaginary
    Negative**Real        -> not Imaginary if exponent is not Rational
    """
    if expr.is_number:
        return _number_imaginary(expr, assumptions)

    if expr.base.func == exp:
        if ask(Q.imaginary(expr.base.args[0]), assumptions):
            if ask(Q.imaginary(expr.exp), assumptions):
                return False
            i = expr.base.args[0]/I/pi
            if ask(Q.integer(2*i), assumptions):
                return ask(Q.imaginary(((-1)**i)**expr.exp), assumptions)

    if ask(Q.imaginary(expr.base), assumptions):
        if ask(Q.integer(expr.exp), assumptions):
            odd = ask(Q.odd(expr.exp), assumptions)
            if odd is not None:
                return odd
            return

    if ask(Q.imaginary(expr.exp), assumptions):
        imlog = ask(Q.imaginary(log(expr.base)), assumptions)
        if imlog is not None:
            return False  # I**i -> real; (2*I)**i -> complex ==> not imaginary

    if ask(Q.real(expr.base) & Q.real(expr.exp), assumptions):
        if ask(Q.positive(expr.base), assumptions):
            return False
        else:
            rat = ask(Q.rational(expr.exp), assumptions)
            if not rat:
                return rat
            if ask(Q.integer(expr.exp), assumptions):
                return False
            else:
                half = ask(Q.integer(2*expr.exp), assumptions)
                if half:
                    return ask(Q.negative(expr.base), assumptions)
                return half

@AskImaginaryHandler.register(log)
def _(expr, assumptions):
    if ask(Q.real(expr.args[0]), assumptions):
        if ask(Q.positive(expr.args[0]), assumptions):
            return False
        return
    # XXX it should be enough to do
    # return ask(Q.nonpositive(expr.args[0]), assumptions)
    # but ask(Q.nonpositive(exp(x)), Q.imaginary(x)) -> None;
    # it should return True since exp(x) will be either 0 or complex
    if expr.args[0].func == exp:
        if expr.args[0].args[0] in [I, -I]:
            return True
    im = ask(Q.imaginary(expr.args[0]), assumptions)
    if im is False:
        return False

@AskImaginaryHandler.register(exp)
def _(expr, assumptions):
    a = expr.args[0]/I/pi
    return ask(Q.integer(2*a) & ~Q.integer(a), assumptions)

@AskImaginaryHandler.register(Number)
@AskImaginaryHandler.register(NumberSymbol)
def _(expr, assumptions):
    return not (expr.as_real_imag()[1] == 0)


### AskAntiHermitianHandler ###

AskAntiHermitianHandler = AskImaginaryHandler.copy(
    'AskAntiHermitianHandler',
    doc="""
    Handler for Q.antihermitian
    Test that an expression belongs to the field of anti-Hermitian operators,
    that is, operators in the form x*I, where x is Hermitian
    """
)

@AskAntiHermitianHandler.register(Expr)
def _(expr, assumptions):
    if isinstance(expr, MatrixBase):
        return None

@AskAntiHermitianHandler.register(Add)
def _(expr, assumptions):
    """
    Antihermitian + Antihermitian  -> Antihermitian
    Antihermitian + !Antihermitian -> !Antihermitian
    """
    if expr.is_number:
        return _number_imaginary(expr, assumptions)
    return test_closed_group(expr, assumptions, Q.antihermitian)

@AskAntiHermitianHandler.register(Mul)
def _(expr, assumptions):
    """
    As long as there is at most only one noncommutative term:
    Hermitian*Hermitian         -> !Antihermitian
    Hermitian*Antihermitian     -> Antihermitian
    Antihermitian*Antihermitian -> !Antihermitian
    """
    if expr.is_number:
        return _number_imaginary(expr, assumptions)
    nccount = 0
    result = False
    for arg in expr.args:
        if ask(Q.antihermitian(arg), assumptions):
            result = result ^ True
        elif not ask(Q.hermitian(arg), assumptions):
            break
        if ask(~Q.commutative(arg), assumptions):
            nccount += 1
            if nccount > 1:
                break
    else:
        return result

@AskAntiHermitianHandler.register(Pow)
def _(expr, assumptions):
    """
    Hermitian**Integer  -> !Antihermitian
    Antihermitian**Even -> !Antihermitian
    Antihermitian**Odd  -> Antihermitian
    """
    if expr.is_number:
        return _number_imaginary(expr, assumptions)
    if ask(Q.hermitian(expr.base), assumptions):
        if ask(Q.integer(expr.exp), assumptions):
            return False
    elif ask(Q.antihermitian(expr.base), assumptions):
        if ask(Q.even(expr.exp), assumptions):
            return False
        elif ask(Q.odd(expr.exp), assumptions):
            return True

@AskAntiHermitianHandler.register(MatrixBase)
def _(mat, assumptions):
    rows, cols = mat.shape
    ret_val = True
    for i in range(rows):
        for j in range(i, cols):
            cond = fuzzy_bool(Eq(mat[i, j], -conjugate(mat[j, i])))
            if cond == None:
                ret_val = None
            if cond == False:
                return False
    return ret_val


### AskAlgebraicHandler ###

AskAlgebraicHandler = CommonHandler.copy(
    'AskAlgebraicHandler',
    doc="""Handler for Q.algebraic key. """
)

for sig in (
Float, GoldenRatio, TribonacciConstant, ImaginaryUnit, AlgebraicNumber
):
    AskAlgebraicHandler.register(sig)(AskAlgebraicHandler.AlwaysTrue)

for sig in (
Infinity, NegativeInfinity, ComplexInfinity, Pi, Exp1
):
    AskAlgebraicHandler.register(sig)(AskAlgebraicHandler.AlwaysFalse)

@AskAlgebraicHandler.register(Add)
def _(expr, assumptions):
    return test_closed_group(expr, assumptions, Q.algebraic)

@AskAlgebraicHandler.register(Mul)
def _(expr, assumptions):
    return test_closed_group(expr, assumptions, Q.algebraic)

@AskAlgebraicHandler.register(Pow)
def _(expr, assumptions):
    return expr.exp.is_Rational and ask(
        Q.algebraic(expr.base), assumptions)

@AskAlgebraicHandler.register(Rational)
def _(expr, assumptions):
    return expr.q != 0

@AskAlgebraicHandler.register(exp)
@AskAlgebraicHandler.register(sin)
@AskAlgebraicHandler.register(cos)
@AskAlgebraicHandler.register(tan)
@AskAlgebraicHandler.register(asin)
@AskAlgebraicHandler.register(atan)
def _(expr, assumptions):
    x = expr.args[0]
    if ask(Q.algebraic(x), assumptions):
        return ask(~Q.nonzero(x), assumptions)

@AskAlgebraicHandler.register(cot)
@AskAlgebraicHandler.register(acot)
def _(expr, assumptions):
    x = expr.args[0]
    if ask(Q.algebraic(x), assumptions):
        return False

@AskAlgebraicHandler.register(log)
@AskAlgebraicHandler.register(acos)
def _(expr, assumptions):
    x = expr.args[0]
    if ask(Q.algebraic(x), assumptions):
        return ask(~Q.nonzero(x - 1), assumptions)
