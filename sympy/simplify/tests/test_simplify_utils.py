from sympy import S, C, sin, cos, tan, cot, exp, Add, sinh, cosh
from sympy.abc import x, y, z, a, b, c
from sympy.simplify.simplify_utils import replace_add_fgfg, replace_mul_f2g
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
    assert res == 2*sin(pi/4 + 1)*sin(pi/4 + 3) + 2*sin(pi/4 + 2)*sin(pi/4 + 3)

def test_replace_mul_f2g():
    rep = lambda x: replace_mul_f2g(x, sin, cos)
    assert rep(sin(x)*cos(x)) == sin(2*x)/2
    assert rep(sin(2*x)/sin(x)) == 2*cos(x)
    assert rep(sin(2*x)/cos(x)) == 2*sin(x)
    assert rep(sin(2*x)**2/sin(x)) == 2*sin(2*x)*cos(x)
    assert rep(sin(x)**2*cos(x)**2/sin(2*x)) == sin(2*x)/4
    assert rep(cos(x)*cos(2*x)/sin(4*x)) == 1/(4*sin(x))
    assert rep(cos(x)*cos(2*x)/(sin(x)*sin(4*x))) == 1/(4*sin(x)**2)
    assert rep(sinh(1)*sin(4*x)/(cos(x)*cos(2*x))) == 4*sinh(1)*sin(x)
    assert rep(cos(2*x)/(sin(x)*cos(x)*sin(4*x))) == 1/sin(2*x)**2

    rep = lambda x: replace_mul_f2g(x, sinh, cosh)
    assert rep(sinh(x)*cosh(x)) == sinh(2*x)/2
    assert rep(sinh(2*x)/sinh(x)) == 2*cosh(x)
    assert rep(sinh(2*x)/cosh(x)) == 2*sinh(x)
    assert rep(sinh(2*x)**2/sinh(x)) == 2*sinh(2*x)*cosh(x)
    assert rep(sinh(x)**2*cosh(x)**2/sinh(2*x)) == sinh(2*x)/4
    assert rep(cosh(x)*cosh(2*x)/sinh(4*x)) == 1/(4*sinh(x))
    assert rep(cosh(x)*cosh(2*x)/(sinh(x)*sinh(4*x))) == 1/(4*sinh(x)**2)
    assert rep(sinh(1)*sinh(4*x)/(cosh(x)*cosh(2*x))) == 4*sinh(1)*sinh(x)
    assert rep(cosh(2*x)/(sinh(x)*cosh(x)*sinh(4*x))) == 1/sinh(2*x)**2
