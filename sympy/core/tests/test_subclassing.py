"""
Test the precedence of implementations of arithmetic operators.

Subclasses of Expr need to be allowed to override Expr's implementation
of the builtin arithmetic operators.
"""
from sympy import Expr, Symbol, Integer, Rational, Float, S, sin
from sympy.utilities.pytest import XFAIL

x = Symbol('x')
y = Symbol('y', commutative=False)
_Expr_objects = [Integer(3), Rational(4, 3), Float(2.4), x, x+1, 3*x, S.Pi,
        S.ImaginaryUnit, y, x*y, y**3, sin(x+1)]

class Higher(Expr):

    def _direct(op_name):
        def __op__(self, other):
            if isinstance (other, Expr):
                return 'high'
            return getattr(Expr, op_name)(self, other)
        return __op__

    def _reverse(op_name):
        def __rop__(self, other):
            if isinstance(other, Expr):
                return 'high'
            return getattr(Expr, op_name)(self, other)
        return __rop__

    for name in ["add", "mul", "sub", "div", "pow"]:
        dname = "__%s__" % name
        locals()[dname] =  _direct(dname)
        rname = "__r%s__" % name
        locals()[rname] = _reverse(rname)


class Lower(Expr):

    def _direct(op_name):
        def __op__(self, other):
            if type(other) is type(self):
                return 'low'
            return getattr(Expr, op_name)(self, other)
        return __op__

    def _reverse(op_name):
        def __rop__(self, other):
            if isinstance(other, Lower):
                return 'low'
            return getattr(Expr, op_name)(self, other)
        return __rop__

    for name in ["add", "mul", "sub", "div", "pow"]:
        dname = "__%s__" % name
        locals()[dname] =  _direct(dname)
        rname = "__r%s__" % name
        locals()[rname] = _reverse(rname)

@XFAIL
def test_mul():
    h = Higher()
    l = Lower()
    assert h*h == l*h == h*l == 'high'
    assert l*l == 'low'
    for x in _Expr_objects:
        assert x*h == h*x == 'high'
        assert x*l != 'low'
        assert l*x != 'low'

@XFAIL
def test_add():
    h = Higher()
    l = Lower()
    assert h+h == l+h == h+l == 'high'
    assert l+l == 'low'
    for x in _Expr_objects:
        assert x+h == h+x == 'high'
        assert l+x == x+l != 'low'

@XFAIL
def test_sub():
    h = Higher()
    l = Lower()
    assert h-h == l-h == h-l == 'high'
    assert l -l == 'low'
    for x in _Expr_objects:
        assert x-h == h-x == 'high'
        assert l-x != 'low'
        assert x-l != 'low'

@XFAIL
def test_pow():
    h = Higher()
    l = Lower()
    assert h**h == l**h == h**l == 'high'
    assert l ** l == 'low'
    for x in _Expr_objects:
        assert x**h == h**x == 'high'
        assert l**x != 'low'
        assert x**l != 'low'

@XFAIL
def test_div():
    h = Higher()
    l = Lower()
    assert h/h == l/h == h/l == 'high'
    assert l/l == 'low'
    for x in _Expr_objects:
        assert x/h == h/x == 'high'
        assert l/x != 'low'
        assert x/l != 'low'
