from sympy.core.function import Function
from sympy.core.sympify import _sympify
from sympy.core.add import Add
from sympy.core.mul import Mul
from sympy.core.power import Pow
from sympy.functions.elementary.hyperbolic import (
    sinh, cosh, tanh, csch, sech, coth)
from sympy.functions.elementary.trigonometric import (
    cos, sec, sin, csc, tan, cot)
from sympy.functions.special.error_functions import erf, fresnelc, fresnels
from sympy.functions import Abs, sign


_even_univariate_functions = (
    Abs, cos, sec, cosh, sech
)

_odd_univariate_functions = (
    sign, sin, csc, tan, cot, sinh, csch, tanh, coth,
    erf, fresnels, fresnelc
)

_any_functions = (
    Function
)


def is_even_function(expr, *symbols):
    r"""Tests whether a function is even with respect to the given
    symbols.

    Explanation
    ===========

    A univariate function $f$ is an even function with respect to
    $x$ if the following property holds.

    .. math::
        f(-x) = f(x)

    A multivariate function $f$ is an even function with respect to
    $x_1, x_2, ..., x_n$ if the following property holds.

    .. math::
        f(-x_1, -x_2, ..., -x_n) = f(x_1, x_2, ..., x_n)

    Parameters
    ==========

    expr : Expr
        The expression to test the function parity.

    symbols : Symbol
        The symbols to test the function parity with respect to.

    Notes
    =====

    A function $f$ is even with respect to $x$:

    1. If $f(x) = \sum_{i=1}^{n} f_i(x)$ and
       $f_1, f_2, ..., f_n$ are even functions.

    2. If $f(x) = \prod_{i=1}^{m} f_i(x) \prod_{j=1}^{n} g_j(x)$ and
       $f_1, f_2, ..., f_m$ are even functions,
       $g_1, g_2, ..., g_n$ are odd functions,
       $n$ is an even integer.

    3. If $f(x) = g(x)^{h(x)}$ where $g$ is an even function, $h$ is an
       even function.

    4. If $f(x) = g(x)^{n}$ where $g$ is an odd function, $n$ is an even
       integer.

    5. If $f(x) = g(h_1(x), h_2(x), ..., h_n(x))$ where
       $h_1, h_2, ..., h_n$ are even functions.

    6. If $f(x) = g(h_1(x), h_2(x), ..., h_n(x))$ where
       $g$ is an even function and
       $h_1, h_2, ..., h_n$ are odd functions.

    The properties for univariate function carries over the multivariate
    function naturally as you replace $x$ with $x_1, x_2, ..., x_n$ and
    $-x$ with $-x_1, -x_2, ..., -x_n$

    References
    ==========

    .. [#] https://en.wikipedia.org/wiki/Even_and_odd_functions
    """
    expr = _sympify(expr)
    symbols = [_sympify(x) for x in symbols]

    if all([not expr.has(x) for x in symbols]):
        return True

    if isinstance(expr, Add):
        if all([is_even_function(arg, *symbols) for arg in expr.args]):
            return True
        return None

    if isinstance(expr, Mul):
        count_even = 0
        count_odd = 0
        for arg in expr.args:
            if is_even_function(arg, *symbols):
                count_even += 1
                continue
            if is_odd_function(arg, *symbols):
                count_odd += 1
                continue
            return None

        if count_odd % 2 == 0:
            return True
        return None

    if isinstance(expr, Pow):
        base, exp = expr.args
        if is_even_function(base, *symbols) and is_even_function(exp, *symbols):
            return True

        if not exp.is_integer:
            return None

        if is_odd_function(base, *symbols) and exp.is_even:
            return True
        return None

    if isinstance(expr, _even_univariate_functions):
        arg = expr.args[0]
        if is_even_function(arg, *symbols):
            return True
        if is_odd_function(arg, *symbols):
            return True
        return None

    if isinstance(expr, _odd_univariate_functions):
        arg = expr.args[0]
        if is_even_function(arg, *symbols):
            return True
        return None

    if isinstance(expr, _any_functions):
        if all([is_even_function(arg, *symbols) for arg in expr.args]):
            return True
        return None

    return None


def is_odd_function(expr, *symbols):
    r"""Tests whether a function is odd with respect to the given
    symbols.

    Explanation
    ===========

    A univariate function $f$ is an odd function with respect to
    $x$ if the following property holds.

    .. math::
        f(-x) = -f(x)

    A multivariate function $f$ is an odd function with respect to
    $x_1, x_2, ..., x_n$ if the following property holds.

    .. math::
        f(-x_1, -x_2, ..., -x_n) = -f(x_1, x_2, ..., x_n)

    Parameters
    ==========

    expr : Expr
        The expression to test the function parity.

    symbols : Symbol
        The symbols to test the function parity with respect to.

    Notes
    =====

    A function $f$ is odd with respect to $x$:

    1. If $f(x) = \sum_{i=1}^{n} f_i(x)$ and
       $f_1, f_2, ..., f_n$ are odd functions.

    2. If $f(x) = \prod_{i=1}^{m} f_i(x) \prod_{j=1}^{n} g_j(x)$ and
       $f_1, f_2, ..., f_m$ are even functions,
       $g_1, g_2, ..., g_n$ are odd functions,
       $n$ is an odd integer.

    3. If $f(x) = g(x)^{n}$ where $g$ is an odd function, $n$ is an odd
       integer.

    4. If $f(x) = g(h_1(x), h_2(x), ..., h_n(x))$ where
       $g$ is an odd function and
       $h_1, h_2, ..., h_n$ are odd functions.

    The properties for univariate function carries over the multivariate
    function naturally as you replace $x$ with $x_1, x_2, ..., x_n$ and
    $-x$ with $-x_1, -x_2, ..., -x_n$

    References
    ==========

    .. [#] https://en.wikipedia.org/wiki/Even_and_odd_functions
    """
    expr = _sympify(expr)
    symbols = [_sympify(x) for x in symbols]

    if expr in symbols:
        return True

    if all([not expr.has(x) for x in symbols]):
        if expr.is_zero:
            return True
        return None

    if isinstance(expr, Add):
        if all([is_odd_function(arg, *symbols) for arg in expr.args]):
            return True
        return None

    if isinstance(expr, Mul):
        count_even = 0
        count_odd = 0
        for arg in expr.args:
            if is_even_function(arg, *symbols):
                count_even += 1
                continue
            if is_odd_function(arg, *symbols):
                count_odd += 1
                continue
            return None

        if count_odd % 2 == 1:
            return True
        return None

    if isinstance(expr, Pow):
        base, exp = expr.args
        if is_odd_function(base, *symbols) and exp.is_odd:
            return True

    if isinstance(expr, _odd_univariate_functions):
        arg = expr.args[0]
        if is_odd_function(arg, *symbols):
            return True

    return None
