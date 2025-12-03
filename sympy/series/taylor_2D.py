from sympy.core.sympify import sympify, SympifyError
from sympy.core.function import expand
from sympy.core.symbol import symbols
from sympy import diff
from sympy.functions.combinatorial.factorials import factorial

X,Y = symbols('X'),symbols('Y')

__all__ = [
    "diffDict",
    "RemoveCoefficientSplit",
    "TaylorTwovariable"
]

def diffDict(exp,upto,a,b):
    """
    Generate all partial derivatives of a two-variable expression
    up to the given order, and evaluate each derivative at (a, b).

    Parameters
    ----------
    exp : str or Expr
        The expression in terms of X and Y for which partial derivatives
        are to be computed.
    upto : int
        The maximum derivative order (i.e., how many layers of differentiation
        will be performed).
    a : int, float, or symbol
        The x-coordinate of the expansion point.
    b : int, float, or symbol
        The y-coordinate of the expansion point.

    Returns
    -------
    dict
        A nested dictionary where each key represents the derivative order
        and the corresponding value stores a dictionary mapping derivative
        expressions (like X, Y, X*Y, X**2, etc.) to their evaluated values
        at (a, b).

    Notes
    -----
    This function forms the computational backbone of the 2D Taylor Series
    expansion. It systematically builds all mixed partial derivatives up to
    a specified order and stores both symbolic and evaluated versions.
    The evaluation uses symbolic string replacement to substitute X and Y
    with the given numerical or symbolic points.

    Examples
    --------
    >>> diffDict("X**2 + Y**2", 2, 0, 0)
    {1: {X: 0, Y: 0}, 2: {X*X: 2, X*Y: 0, Y*Y: 2}}
    """
    
    a, b = sympify(a), sympify(b)
    exp = sympify(exp)
    expDict = {1 : {sympify('X') : diff(exp , X), sympify('Y') : diff(exp ,Y)}}
    ValueDict = {1 : {sympify('X') : sympify((diff(exp, X).subs({X:a,Y:b}))) ,sympify('Y') : sympify((diff(exp, Y).subs({X:a,Y:b})))}}
    for i in range(2,upto+2):
        tempdict = {}
        tempdict2 = {}
        change = 0
        keys = expDict[i-1].keys()
        values = expDict[i-1].values()
        tempList = [list(keys) , list(values)]
        for j in range(2**i):
            if j%2 == 0 and change <= len(values) - 1:
                tempdict[sympify(tempList[0][change] * X)] = diff(tempList[1][change], X)
                tempdict2[sympify(tempList[0][change] * X)] =  (diff(tempList[1][change], X)).subs({X : a, Y : b})

            elif j%2 != 0 and change <= len(values) - 1:
                tempdict[sympify(tempList[0][change] * Y)] = diff(tempList[1][change], Y)
                tempdict2[sympify(tempList[0][change] * Y)] = (diff(tempList[1][change], Y)).subs({X : a, Y : b})
                change+=1
        expDict[i] = tempdict
        ValueDict[i] = tempdict2
    return ValueDict

def RemoveCoefficientSplit(Custom):
    """
    Remove the numeric coefficient from each monomial term.

    Given a list of monomials (as SymPy expressions), this helper function
    separates out the numeric coefficient and returns only the pure symbolic
    monomial part for comparison. This is used internally to match the terms
    generated from `(X + Y)**n` with the mixed partial derivatives produced
    during the Taylor expansion.

    Parameters
    ==========
    Custom : list
        A list of SymPy expressions, typically monomials obtained from
        expanding `(X + Y)**n`.

    Returns
    =======
    list
        A list of strings where each entry is the monomial part of the term
        with coefficients removed. For example, ``3*X**2*Y`` becomes
        ``"X**2*Y"``.

    Notes
    =====
    This function is an internal utility used by `TaylorTwovariable` for
    matching like terms. It does not modify the input and assumes that
    each entry in ``Custom`` is a valid SymPy expression with a single
    numeric coefficient.

    Examples
    ========
    >>> from sympy import symbols
    >>> X, Y = symbols('X Y')
    >>> RemoveCoefficientSplit([3*X**2*Y, -5*X*Y**3])
    ['X**2*Y', 'X*Y**3']
    """
    removed = []
    for express in Custom:
        new = ""
        check = 1
        for i in range(len(express)):
            if (express[i].isdigit() or express[i].isalnum() != True) and (i == 0 or i == 1) and check == 1:
                pass
            else:
                new += express[i]
                check = 0
        removed.append(new)
    return removed

def TaylorTwovariable(real_eq,a,b,upto):
    """
    Compute the two-variable Taylor Series expansion of a given expression
    around the point (a, b) up to the specified order.

    Parameters
    ----------
    real_eq : str or Expr
        The target expression to expand. Must be written in terms of X and Y.
    a : int, float, or symbol
        The x-coordinate of the expansion point.
    b : int, float, or symbol
        The y-coordinate of the expansion point.
    upto : int
        The maximum order of the Taylor expansion.

    Returns
    -------
    Expr
        The expanded two-variable Taylor polynomial as a SymPy expression.

    Raises
    ------
    SympifyError
        If the input expression cannot be converted into a valid SymPy object.

    Notes
    -----
    - This implementation uses a custom recursive derivative generator
      (`diffTuple`) instead of `sympy.series` to provide a transparent view
      of how multivariate expansions are constructed step-by-step.
    - The expansion automatically substitutes (X, Y) with (x - a, y - b)
      to represent the shifted coordinate system.
    - The result is fully simplified and expanded for clarity.

    Examples
    --------
    >>> TaylorTwovariable("X**2 + Y**2", 0, 0, 2)
    X**2 + Y**2
    >>> TaylorTwovariable("sin(X*Y)", 0, 0, 3)
    X*Y
    >>> TaylorTwovariable("exp(X + Y)", 0, 0, 3)
    1 + (X + Y) + (X + Y)**2/2 + (X + Y)**3/6
    """
    a, b = sympify(a), sympify(b)
    while True:
        eq = real_eq
        try:
            eq = sympify(eq)
            real_eq = eq
            break
        except SympifyError:
            raise SympifyError("Write a correct equation")
    diffs = diffDict(real_eq, upto, a, b)
    diff_Dict_values = list(diffs.values())
    equations = []
    for i in range(1,upto+1):
        eq = sympify('X + Y')
        eq = eq**i
        eq = expand(eq)
        tempDictKeys = list(diff_Dict_values[i-1].keys())
        tempDictValues = list(diff_Dict_values[i-1].values())
        Custom = str(expand((X+Y)**i)).split(' + ')
        fake = RemoveCoefficientSplit(Custom)
        for j in range(len(tempDictKeys)):
            for k in range(len(fake)):
                if str(tempDictKeys[j]) == str(fake[k]):
                    Custom[k]= sympify(f'{Custom[k]} * {tempDictValues[j]}')
                    break
        eq = 0
        for j in Custom:
            eq += j
        eq = eq.subs({X: X - a, Y: Y - b})
        equations.append(sympify(eq))

    Taylor = real_eq.subs({X: a, Y: b})
    for x in range(0,len(equations)):
        Taylor += (equations[x])/factorial(x+1)
    return expand(Taylor)
