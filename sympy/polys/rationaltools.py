"""Tools for manipulation of rational expressions. """


from sympy.core import Basic, Add, sympify, symbols, nan, zoo
from sympy.core.exprtools import gcd_terms
from sympy.utilities import public
from sympy.utilities.iterables import iterable


@public
def together(expr, deep=False, fraction=True):
    """
    Denest and combine rational expressions using symbolic methods.

    This function takes an expression or a container of expressions
    and puts it (them) together by denesting and combining rational
    subexpressions. No heroic measures are taken to minimize degree
    of the resulting numerator and denominator. To obtain completely
    reduced expression use :func:`~.cancel`. However, :func:`~.together`
    can preserve as much as possible of the structure of the input
    expression in the output (no expansion is performed).

    A wide variety of objects can be put together including lists,
    tuples, sets, relational objects, integrals and others. It is
    also possible to transform interior of function applications,
    by setting ``deep`` flag to ``True``.

    By definition, :func:`~.together` is a complement to :func:`~.apart`,
    so ``apart(together(expr))`` should return expr unchanged. Note
    however, that :func:`~.together` uses only symbolic methods, so
    it might be necessary to use :func:`~.cancel` to perform algebraic
    simplification and minimize degree of the numerator and denominator.

    Examples
    ========

    >>> from sympy import together, exp
    >>> from sympy.abc import x, y, z

    >>> together(1/x + 1/y)
    (x + y)/(x*y)
    >>> together(1/x + 1/y + 1/z)
    (x*y + x*z + y*z)/(x*y*z)

    >>> together(1/(x*y) + 1/y**2)
    (x + y)/(x*y**2)

    >>> together(1/(1 + 1/x) + 1/(1 + 1/y))
    (x*(y + 1) + y*(x + 1))/((x + 1)*(y + 1))

    >>> together(exp(1/x + 1/y))
    exp(1/y + 1/x)
    >>> together(exp(1/x + 1/y), deep=True)
    exp((x + y)/(x*y))

    >>> together(1/exp(x) + 1/(x*exp(x)))
    (x + 1)*exp(-x)/x

    >>> together(1/exp(2*x) + 1/(x*exp(3*x)))
    (x*exp(x) + 1)*exp(-3*x)/x

    """
    def _together(expr):
        if isinstance(expr, Basic):
            if expr.is_Atom or (expr.is_Function and not deep):
                return expr
            elif expr.is_Add:
                return gcd_terms(list(map(_together, Add.make_args(expr))), fraction=fraction)
            elif expr.is_Pow:
                base = _together(expr.base)

                if deep:
                    exp = _together(expr.exp)
                else:
                    exp = expr.exp

                return expr.func(base, exp)
            else:
                return expr.func(*[ _together(arg) for arg in expr.args ])
        elif iterable(expr):
            return expr.__class__([ _together(ex) for ex in expr ])

        return expr

    return _together(sympify(expr))


def thiele_interpolate(u, v, var=symbols('x'), simplify=True):
    """
    Build a rational function from a finite set of inputs in vector u
    and their function values in vector v (both vectors must be of equal
    lengths).

    This algorithm might encounter division by zero
    when values are equally spaced. For such cases, a better algorithm is
    sympy.polys.polyfuncs.rational_interpolate.

    Furthermore, the simplest solution is not always returned (see the last
    example below).

    An arbitrary symbol can optionally be provided as var
    (default being 'x').

    At each step of the algorithm, some simplification is done, which can
    optionally be prevented by setting simplify as False.

    Examples
    ========

    >>> from sympy.polys.rationaltools import thiele_interpolate as thiele

    >>> thiele([1, 2, 5, 6], [10, 12, 11, 13])
    (9*x**2 + 29*x - 238)/(8*x - 28)

    >>> from sympy import S
    >>> thiele([1, 2, 3, 4], [S.One, S.One/2, S.One/3, S.One/4])
    1/x

    >>> thiele([1, 2, 3, 4], [S.One, S.One/4, S.One/9, S.One/16])
    (x**2 - 10*x + 35)/(50*x - 24)

    See Also
    ========

    sympy.polys.polyfuncs.rational_interpolate : another algorithm

    """
    n = len(u)
    assert len(v) == n
    u = sympify(u)
    v = sympify(v)
    rho = [v, [(u[i] - u[i + 1])/(v[i] - v[i + 1]) for i in range(n - 1)]]
    for i in range(n - 2):
        r = []
        for j in range(n - i - 2):
            r.append((u[j] - u[j + i + 2])/(rho[-1][j] - rho[-1][j + 1]) + rho[-2][j + 1])
        rho.append(r)
    a = 0
    for i in range(n - 1, 1, -1):
        a = (var - u[i - 1]) / (rho[i][0] - rho[i - 2][0] + a)
        if simplify:
            a = a.cancel()
        if a.has(zoo, nan):
            raise ZeroDivisionError("division by zero")
    a = v[0] + (var - u[0]) / (rho[1][0] + a)
    if simplify:
        a = a.cancel()
    if a.has(zoo, nan):
        raise ZeroDivisionError("division by zero")
    return a
