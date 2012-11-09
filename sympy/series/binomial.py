from sympy.core import Mul, Add, Pow, sympify
from sympy.core.symbol import symbols, Symbol


def binomial_expand(function):
    """
    Returns a valid binomial expansion of an expression, if there is any.
    Otherwise, returns the function itself. The output is according to
    that given by WolframAlpha, for functions with non-integer values.
    Output depends on the value of exponent.

    Examples
    ========

    >>> from sympy import binomial_expand, symbols
    >>> x,y,z = symbols('x,y,z')
    >>> binomial_expand((x + y + z) ** 2)
        x**2 + 2*x*y + 2*x*z + y**2 + 2*y*z + z**2
    >>> binomial_expand((1 + x) ** 0.3)
        1.0 + 0.3*x - 0.105*x**2 + 0.0595*x**3 - 0.0401625*x**4 +
        0.02972025*x**5 + O(x**6)

    """
    function = sympify(function).evalf()
    if function.is_Pow:
        xin, yin = symbols('xin,yin')
        if not(function.base.is_Add):
            return function
        exp = function.exp
        args_list = list(function.base.args)
        req = args_list.pop(1)
        other = Add(*args_list)
        args_list = list((((1 + xin) ** exp).series()).args)
        i = 0
        for i in range(len(args_list)):
            if args_list[i].is_Mul or args_list[i].is_Pow:
                args_list[i] = args_list[i].subs(xin, (xin / yin))
            i += 1
        i = 0
        for i in range(len(args_list)):
            if args_list[i].is_Mul:
                temp = list(args_list[i].args)
                j = 0
                for j in range(len(temp)):
                    if temp[j].is_Pow:
                        if temp[j].base == yin:
                            new = Pow(yin, (temp[j].exp + exp))
                            temp[j] = new
                    j += 1
                args_list[i] = Mul(*temp)
                i += 1
            else:
                args_list[i] = (args_list[i] * (other ** exp)).evalf()
                i += 1
        step1 = ((Add(*args_list)).subs(xin, req))
        step2 = (step1.subs(yin, other))
        return step2.expand()
    else:
        return function
