from sympy import S, C, sin, cos, tan, cot, exp, Add
from sympy.abc import x, y, z, a, b, c
from sympy.simplify.simplify_utils import replace_add_fgfg
from sympy.simplify.simplify import TR10_inv
from sympy.core.function import expand_multinomial, expand_mul
from sympy.core.compatibility import ordered
from sympy.core.numbers import pi

def _mexpand(expr):
    return expand_mul(expand_multinomial(expr))


def TR10(rv):
    """
    Examples
    ========

    >>> from sympy.simplify.fu import TR10
    >>> from sympy.abc import a, b, c
    >>> from sympy import cos, sin
    >>> TR10(cos(a+b))
    -sin(a)*sin(b) + cos(a)*cos(b)
    >>> TR10(sin(a+b))
    sin(a)*cos(b) + sin(b)*cos(a)
    >>> TR10(sin(a+b+c))
    (-sin(a)*sin(b) + cos(a)*cos(b))*sin(c) + (sin(a)*cos(b) + sin(b)*cos(a))*cos(c)
    """
    if rv.func not in (cos, sin):
        return rv.replace(lambda x: x.func in (cos, sin), lambda x: TR10(x))
    f = rv.func
    arg = rv.args[0]  # should expand_mul be used?
    if arg.is_Add:
        args = list(ordered(arg.args))
        a = args.pop()
        b = Add(*args)
        if b.is_Add:
            if f == sin:
                return sin(a)*TR10(cos(b)) + cos(a)*TR10(sin(b))
            else:
                return cos(a)*TR10(cos(b)) - sin(a)*TR10(sin(b))
        else:
            if f == sin:
                return sin(a)*cos(b) + cos(a)*sin(b)
            else:
                return cos(a)*cos(b) - sin(a)*sin(b)
    return rv


def test_TR10_inv():
    expr1 = (15*sin(b + c) + 11*cos(1)) * \
       (15*sin(x + 2) + 13*cos(y - 1) + 4*cos(z - 1) + 2*sin(3))
    expr1e = _mexpand(expr1)
    expr = _mexpand(TR10(expr1))
    res = TR10_inv(expr)
    assert expr1e == res

    expr = _mexpand(expr/cos(c))
    expr1e = _mexpand(expr1e/cos(c))
    res = TR10_inv(expr)
    assert expr1e == res

    expr1 = (11*sin(a + b) + 11*sin(b + 4) + sin(b + c) + 15*cos(a - 3) + \
            20*cos(b - 4) + 4*sin(7) + 31*cos(1)) * \
           (15*sin(x + 2) + 13*sin(y + 1) + 4*sin(z + 1) + 7*cos(x - z) + \
            16*cos(y - z) + 2*sin(4) + 11)
    expr1e = _mexpand(expr1)
    expr = _mexpand(TR10(expr1))
    assert len(expr.args) == 144 and len(expr1e.args) == 49
    res = TR10_inv(expr)
    assert expr1e == res

    expr = _mexpand((sin(1) + cos(1) + sin(2) + cos(2))*(sin(3) + cos(3)))
    res = TR10_inv(expr)
    # res can take a few different forms, so we do the following test
    res1 = _mexpand(res.rewrite((sin,cos),exp).expand())
    expr1 = _mexpand(expr.rewrite((sin,cos),exp).expand())
    assert res1 == expr1
    assert expr.count(C.TrigonometricFunction) == 16
    assert res.count(C.TrigonometricFunction) == 4
