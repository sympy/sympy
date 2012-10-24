from sympy.core import Mul, Add, Pow
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
    if function.__class__ == Pow:
        xin, yin = symbols('xin,yin')
        exp = function.exp
        args_list = list(function.base.args)
        if len(args_list) == 0 or len(args_list) == 1:
            return function
        req = args_list.pop(1)
        other = Add(*args_list)
        args_list = list((((1 + xin) ** exp).series()).args)
        counter = 0
        while counter <= (len(args_list) - 1):
            if args_list[counter].__class__ == Mul or args_list[counter].__class__ == Pow:
                args_list[counter] = args_list[counter].subs(xin, (xin / yin))
            counter += 1
        counter = 0
        while counter <= (len(args_list) - 1):
            if args_list[counter].__class__ == Mul:
                temp = list(args_list[counter].args)
                counter2 = 0
                while counter2 <= (len(temp) - 1):
                    if temp[counter2].__class__ == Pow:
                        if temp[counter2].base == yin:
                            new = Pow(yin, (temp[counter2].exp + exp))
                            temp[counter2] = new
                    counter2 += 1
                args_list[counter] = Mul(*temp)
                counter += 1
            else:
                args_list[counter] = (args_list[counter] * (other ** exp)).evalf()
                counter += 1
        step1 = ((Add(*args_list)).subs(xin, req))
        step2 = (step1.subs(yin, other))
        return step2.expand()
    else:
        return function
