from sympy import Wild, solve, simplify, log, exp


def pod(expr, sym):
    # pod: Points of discontinuty, (will replace with a better name)
    """
    Find the points of Discontinuity of a real univariate function

    Example
    =========

    >>> from sympy.calculus.discontinuity import pod
    >>> from sympy import exp, log, Symbol
    >>> x = Symbol('x', real=True)
    >>> pod(log((x-2)**2) + x**2, x)
    [2]
    >>> pod(exp(1/(x-2)**2), x)
    [2]

    """
    #  For now this a hacky implementation.

    # not trying for trig function because they have infinitely
    # many solutions, and this can turn out to be problematic
    # for example solve(tan(x), x) returns only 0 for now
    if expr.is_polynomial():
        return []
    pods = []
    pods = pods + solve(simplify(1/expr), sym)
    p = Wild("p")
    q = Wild("q")
    r = Wild("r")

    # check the condition for log
    expr_dict = expr.match(r*log(p) + q)
    if not expr_dict[r].is_zero:
        pods += solve(expr_dict[p], sym)
        pods += pod(expr_dict[p], sym)
        pods += pod(expr_dict[r], sym)

    # check the condition for exp
    expr = expr.rewrite(exp)
    expr_dict = expr.match(r*exp(p) + q)
    if not expr_dict[r].is_zero:
        pods += solve(simplify(1/expr_dict[p]), sym)
        pods += pod(expr_dict[p], sym)
        pods += pod(expr_dict[r], sym)

    return list(set(pods)) # remove dublications
