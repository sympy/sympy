from sympy.core import S
from sympy.logic.boolalg import And, Or, Not, Implies, Equivalent, Xor, Nand, BooleanAtom
from sympy.logic.boolalg import Nor, Xnor, ITE

def boolean_to_polynomial(expr, simplify=False, normalize=False, track_steps=False, _cache=None, _verbose=False, _steps=None):
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
    normalize : bool, optional
        Whether to normalize terms (canonical ordering).
    track_steps : bool, optional
        If True, also return a list of transformation steps.
    _cache : dict, internal
        Cache for memoization.
    _verbose : bool, optional
        If True, print debug info.
    _steps : list, internal
        Internal list of steps (for track_steps).

    Returns
    -------
    SymPy expression (arithmetic polynomial), or (polynomial, steps) if track_steps=True
    """
    from sympy.core.symbol import Symbol
    from sympy.core import Add, Mul, Pow
    from sympy.simplify import collect

    if _cache is None:
        _cache = {}

    if track_steps and _steps is None:
        _steps = []

    # Check cache first
    if expr in _cache:
        result = _cache[expr]
        if _verbose:
            print(f"Cache hit for {expr}")
        if track_steps:
            _steps.append((expr, result))
            return (result, _steps)
        return result

    # Helper function to simplify polynomials
    def simplify_polynomial(poly):
        """
        Simplifies polynomials:
        - Quadratic terms: a^2 = a
        - Nested Adds/Muls flattened and collected

        Parameters
        ----------
        poly : SymPy expression

        Returns
        -------
        Simplified SymPy expression
        """
        if isinstance(poly, Mul):
            args = list(poly.args)
            seen = set()
            simplified_args = []
            for arg in args:
                if isinstance(arg, Pow):
                    base, exp = arg.args
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

    # Base case: Boolean atoms
    if isinstance(expr, BooleanAtom):
        result = S.One if expr == S.true else S.Zero

    # Symbolic variables
    elif isinstance(expr, Symbol):
        result = expr

    # NOT operation
    elif isinstance(expr, Not):
        arg = expr.args[0]
        result = S.One - boolean_to_polynomial(arg, simplify, normalize, track_steps, _cache, _verbose, _steps)

    # AND operation
    elif isinstance(expr, And):
        args = [boolean_to_polynomial(arg, simplify, normalize, track_steps, _cache, _verbose, _steps) for arg in expr.args]
        result = Mul(*args)
        if simplify:
            result = simplify_polynomial(result)

    # OR operation
    elif isinstance(expr, Or):
        args = [boolean_to_polynomial(arg, simplify, normalize, track_steps, _cache, _verbose, _steps) for arg in expr.args]
        result = S.Zero
        for arg in args:
            result += arg - result * arg
        if simplify:
            result = simplify_polynomial(result)

    # IMPLIES operation
    elif isinstance(expr, Implies):
        a = boolean_to_polynomial(expr.args[0], simplify, normalize, track_steps, _cache, _verbose, _steps)
        b = boolean_to_polynomial(expr.args[1], simplify, normalize, track_steps, _cache, _verbose, _steps)
        result = S.One - a + a * b

    # EQUIVALENT operation
    elif isinstance(expr, Equivalent):
        args = [boolean_to_polynomial(arg, simplify, normalize, track_steps, _cache, _verbose, _steps) for arg in expr.args]
        result = S.One
        for i in range(len(args) - 1):
            ai = args[i]
            aj = args[i+1]
            eq_pair = S.One - ai - aj + 2 * ai * aj
            result *= eq_pair
        if simplify:
            result = simplify_polynomial(result)

    # XOR operation
    elif isinstance(expr, Xor):
        args = [boolean_to_polynomial(arg, simplify, normalize, track_steps, _cache, _verbose, _steps) for arg in expr.args]
        result = S.Zero
        for arg in args:
            result += arg - 2 * result * arg
        if simplify:
            result = simplify_polynomial(result)

    # NAND operation
    elif isinstance(expr, Nand):
        a = boolean_to_polynomial(expr.args[0], simplify, normalize, track_steps, _cache, _verbose, _steps)
        b = boolean_to_polynomial(expr.args[1], simplify, normalize, track_steps, _cache, _verbose, _steps)
        result = S.One - a * b
        if simplify:
            result = simplify_polynomial(result)

    # NOR operation
    elif isinstance(expr, Nor):
        a = boolean_to_polynomial(expr.args[0], simplify, normalize, track_steps, _cache, _verbose, _steps)
        b = boolean_to_polynomial(expr.args[1], simplify, normalize, track_steps, _cache, _verbose, _steps)
        result = S.One - (a + b - a * b)
        if simplify:
            result = simplify_polynomial(result)

    # XNOR operation
    elif isinstance(expr, Xnor):
        a = boolean_to_polynomial(expr.args[0], simplify, normalize, track_steps, _cache, _verbose, _steps)
        b = boolean_to_polynomial(expr.args[1], simplify, normalize, track_steps, _cache, _verbose, _steps)
        result = S.One - a - b + 2 * a * b
        if simplify:
            result = simplify_polynomial(result)

    # ITE operation
    elif isinstance(expr, ITE):
        a = boolean_to_polynomial(expr.args[0], simplify, normalize, track_steps, _cache, _verbose, _steps)
        b = boolean_to_polynomial(expr.args[1], simplify, normalize, track_steps, _cache, _verbose, _steps)
        c = boolean_to_polynomial(expr.args[2], simplify, normalize, track_steps, _cache, _verbose, _steps)
        result = a * b + (S.One - a) * c
        if simplify == "deep":
            result = simplify_polynomial(result)

    # Unsupported types
    else:
        raise TypeError(f"Unsupported expression type: {type(expr)}")

    if normalize:
        result = simplify_polynomial(result)

    _cache[expr] = result

    if _verbose:
        print(f"Computed polynomial for {expr}: {result}")

    if track_steps:
        _steps.append((expr, result))
        return (result, _steps)

    return result

def is_boolean_polynomial(expr):
    """
    Checks if an expression is a valid Boolean polynomial:
    - Only involves +, *, -, integers 0/1 and symbols
    - No powers >1

    Parameters
    ----------
    expr : SymPy expression

    Returns
    -------
    bool
    """
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
    """
    Computes the maximum degree of a Boolean polynomial.

    Parameters
    ----------
    expr : SymPy expression

    Returns
    -------
    int
    """
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
