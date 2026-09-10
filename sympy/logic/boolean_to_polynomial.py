from sympy.core import S
from sympy.logic.boolalg import And, Or, Not, Implies, Equivalent, Xor, Nand, BooleanAtom
from sympy.logic.boolalg import Nor, Xnor, ITE
from sympy.core.cache import cacheit


def _simplify_polynomial(poly):
    """
    Internal helper to simplify polynomials:
    - Reduces powers (e.g., a^2 -> a)
    - Flattens nested Add/Mul
    - Collects terms by symbol
    """
    from sympy.core import Add, Mul, Pow
    from sympy.simplify import collect

    if isinstance(poly, Mul):
        args = list(poly.args)
        seen = set()
        simplified_args = []
        for arg in args:
            if isinstance(arg, Pow):
                base, _ = arg.args
                simplified_args.append(base)
            elif arg in seen:
                continue
            else:
                simplified_args.append(arg)
                seen.add(arg)
        return Mul(*simplified_args)

    elif isinstance(poly, Add):
        symbols = poly.free_symbols
        if symbols:
            poly = collect(poly, symbols)
        return poly

    elif isinstance(poly, Pow):
        base, exp = poly.args
        return base if exp >= 1 else S.One

    return poly


@cacheit
def boolean_to_polynomial(expr, simplify=False):
    """
    Convert a Boolean expression to its equivalent Boolean polynomial
    (Zhegalkin polynomial) over GF(2).

    Supported operators:
    - And
    - Or
    - Not
    - Implies
    - Equivalent (multi-argument supported)
    - Xor (multi-argument supported)
    - Nand
    - Nor
    - Xnor
    - ITE

    The output is an arithmetic expression using {0,1} with +, *, -.

    Parameters
    ----------
    expr : Boolean expression (SymPy expression)
    simplify : bool or str, optional
        Whether to apply simplifications ('deep' for recursive simplification).

    Returns
    -------
    SymPy expression (arithmetic polynomial)

    Examples
    --------
    >>> from sympy import symbols
    >>> from sympy.logic.boolalg import And, Or, Not
    >>> from sympy.logic.boolean_to_polynomial import boolean_to_polynomial

    >>> a, b = symbols('a b')

    >>> boolean_to_polynomial(And(a, b)) == a * b
    True

    >>> boolean_to_polynomial(Or(a, b)) == a + b - a*b
    True

    >>> boolean_to_polynomial(Not(a)) == 1 - a
    True
    """
    from sympy.core.symbol import Symbol
    from sympy.core import Mul

    if isinstance(expr, BooleanAtom):
        result = S.One if expr == S.true else S.Zero

    elif isinstance(expr, Symbol):
        result = expr

    elif isinstance(expr, Not):
        arg = expr.args[0]
        result = S.One - boolean_to_polynomial(arg, simplify)

    elif isinstance(expr, And):
        args = [boolean_to_polynomial(arg, simplify) for arg in expr.args]
        result = Mul(*args)
        if simplify:
            result = _simplify_polynomial(result)

    elif isinstance(expr, Or):
        args = [boolean_to_polynomial(arg, simplify) for arg in expr.args]
        result = S.Zero
        for arg in args:
            result += arg - result * arg
        if simplify:
            result = _simplify_polynomial(result)

    elif isinstance(expr, Implies):
        a = boolean_to_polynomial(expr.args[0], simplify)
        b = boolean_to_polynomial(expr.args[1], simplify)
        result = S.One - a + a * b

    elif isinstance(expr, Equivalent):
        args = [boolean_to_polynomial(arg, simplify) for arg in expr.args]
        result = S.One
        for i in range(len(args) - 1):
            ai = args[i]
            aj = args[i + 1]
            eq_pair = S.One - ai - aj + 2 * ai * aj
            result *= eq_pair
        if simplify:
            result = _simplify_polynomial(result)

    elif isinstance(expr, Xor):
        args = [boolean_to_polynomial(arg, simplify) for arg in expr.args]
        result = S.Zero
        for arg in args:
            result += arg - 2 * result * arg
        if simplify:
            result = _simplify_polynomial(result)

    elif isinstance(expr, Nand):
        a = boolean_to_polynomial(expr.args[0], simplify)
        b = boolean_to_polynomial(expr.args[1], simplify)
        result = S.One - a * b
        if simplify:
            result = _simplify_polynomial(result)

    elif isinstance(expr, Nor):
        a = boolean_to_polynomial(expr.args[0], simplify)
        b = boolean_to_polynomial(expr.args[1], simplify)
        result = S.One - (a + b - a * b)
        if simplify:
            result = _simplify_polynomial(result)

    elif isinstance(expr, Xnor):
        a = boolean_to_polynomial(expr.args[0], simplify)
        b = boolean_to_polynomial(expr.args[1], simplify)
        result = S.One - a - b + 2 * a * b
        if simplify:
            result = _simplify_polynomial(result)

    elif isinstance(expr, ITE):
        a = boolean_to_polynomial(expr.args[0], simplify)
        b = boolean_to_polynomial(expr.args[1], simplify)
        c = boolean_to_polynomial(expr.args[2], simplify)
        result = a * b + (S.One - a) * c
        if simplify == "deep":
            result = _simplify_polynomial(result)

    else:
        raise TypeError(f"Unsupported expression type: {type(expr)}")

    return result


def is_boolean_polynomial(expr):
    from sympy.core.symbol import Symbol
    from sympy.core import Add, Mul, Pow

    def _check(e):
        if e.is_Number:
            return e in (S.Zero, S.One)
        elif isinstance(e, Symbol):
            return True
        elif isinstance(e, (Add, Mul)):
            return all(_check(arg) for arg in e.args)
        elif isinstance(e, Pow):
            base, exp = e.args
            return _check(base) and exp == 1
        else:
            return False

    return _check(expr)


def degree(expr):
    from sympy.core import Add, Mul, Pow

    def _deg(e):
        if e.is_Number or isinstance(e, Pow):
            return 0
        elif isinstance(e, Mul):
            return sum(_deg(arg) for arg in e.args)
        elif isinstance(e, Add):
            return max(_deg(arg) for arg in e.args)
        elif e.is_Symbol:
            return 1
        return 0

    return _deg(expr)
