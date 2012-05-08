from sympy import Symbol
from sympy.functions.elementary.exponential import exp, log
from sympy.functions.elementary.trigonometric import cos, sin, tan, acos, asin, atan
from sympy.functions.elementary.hyperbolic import cosh, sinh, tanh, acosh, asinh, atanh

class TYPE(int):
    def __repr__(self):
            return self.name
    __str__ = __repr__

LINEAR = TYPE(1); LINEAR.name = 'LINEAR'
POLYNOMIAL = TYPE(2); POLYNOMIAL.name = 'POLYNOMIAL'
RATIONAL = TYPE(3); RATIONAL.name = 'RATIONAL'
ROOT = TYPE(4); ROOT.name = 'ROOT'
EXP = TYPE(5); EXP.name = 'EXP'
LOG = TYPE(6); LOG.name = 'LOG'
HYPERBOLIC = TYPE(7); HYPERBOLIC.name = 'HYPERBOLIC'
TRIGONOMETRIC = TYPE(8); TRIGONOMETRIC.name = 'TRIGONOMETRIC'
INVHYPERBOLIC = TYPE(9); INVHYPERBOLIC.name = 'INVHYPERBOLIC'
INVTRIGONOMETRIC = TYPE(10); INVTRIGONOMETRIC.name = 'INVTRIGONOMETRIC'
NONELEMENTARY = TYPE(11); NONELEMENTARY.name = 'NONELEMENTARY'

f = function_complexity = {}
f[exp] = EXP
f[log] = LOG
f[cos] = f[sin] = f[tan] = TRIGONOMETRIC
f[cosh] = f[sinh] = f[tanh] = HYPERBOLIC
f[acos] = f[asin] = f[atan] = INVTRIGONOMETRIC
f[acosh] = f[asinh] = f[atanh] = INVHYPERBOLIC

def complexity(d):
    return sorted(d, reverse=True)

def max_dep(deps):
    deps = [d for d in deps if d]
    if not deps:
        return deps
    return max(deps, key=complexity)

def classify(expr, x):
    """
    Determining "class" of formula

    The function classify(expr, x) walks the expression top-down and at each
    level classifies it as a function of x.

    Examples
    ========
    >>> from sympy import classify
    >>> from sympy import Symbol,sin
    >> x = Symbol('x')
    >>> classify(x + 1, x)
    [LINEAR]
    >>> classify(sin(x +1), x)
    [TRIGONOMETRIC, LINEAR]
    >>> classify((sin(x) + 1)**x, x)
    [EXP, LOG, LINEAR, TRIGONOMETRIC, LINEAR]
    """
    if not expr.has(x):
        return []
    if expr == x:
        return [LINEAR]
    if expr.is_Add:
        d = max_dep(classify(arg, x) for arg in expr.args)
        if d:
            # [linear, rational, ...] -> [rational, ...]
            if d[0] in (LINEAR, POLYNOMIAL, RATIONAL):
                return d
        return [LINEAR] + d
    if expr.is_Mul:
        deps = filter(None, (classify(arg, x) for arg in expr.args))
        # Only multiplying by a constant
        # [linear, rational, ...] -> [rational, ...]
        mul_dep = LINEAR
        if len(deps) > 1:
            mul_dep = POLYNOMIAL
            # Mul absorbs polynomials and rational functions
            for d in deps:
                if d[0] in (LINEAR, POLYNOMIAL, RATIONAL):
                    mul_dep = max(mul_dep, d[0])
        d = max_dep(deps)
        if d and d[0] <= mul_dep:
            return [mul_dep] + d[1:]
        return [mul_dep] + d
    if expr.is_Pow:
        base, expt = expr.args
        base_dep = classify(base, x)
        expt_dep = classify(expt, x)
        # f(x)^g(x)
        if base_dep and expt_dep:
            return [EXP] + max_dep([[LOG] + base_dep, expt_dep])
        # f(x)^const
        if not expt_dep:
            d = base_dep
            if expt.is_integer:
                if expt.is_positive:
                    if d and (d[0] in (LINEAR, POLYNOMIAL)):
                        d = d[1:]
                    return [POLYNOMIAL] + d
                if d and (d[0] in (LINEAR, POLYNOMIAL, RATIONAL)):
                    d = d[1:]
                return [RATIONAL] + d
            if expt.is_rational:
                return [ROOT] + d
            # f(x)^const = exp(log(f(x)) * const)
            return [EXP, LINEAR, LOG] + base_dep
        # const^f(x)
        if expt_dep and (expt_dep[0] in (LINEAR, POLYNOMIAL, RATIONAL)):
            return [EXP] + expt_dep
        return [EXP, LINEAR] + expt_dep
    fdep = function_complexity.get(expr.func, NONELEMENTARY)
    return [fdep] + max_dep(classify(arg, x) for arg in expr.args)

def is_elementary(expr, x):
    return max(classify(expr, x)) < NONELEMENTARY
