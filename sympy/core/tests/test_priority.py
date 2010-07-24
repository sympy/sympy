from sympy import Expr, Symbol
from sympy.core.expr import call_other

class Higher(Expr):

    _op_priority = 20.0
    result = 'high'

    def __mul__(self, other):
        if call_other(self, other):
            try:
                return other.__rmul__(self)
            except AttributeError:
                pass
        return self.result

    def __rmul__(self, other):
        if call_other(self, other):
            try:
                return other.__mul__(self)
            except AttributeError:
                pass
        return self.result

    def __add__(self, other):
        if call_other(self, other):
            try:
                return other.__radd__(self)
            except AttributeError:
                pass
        return self.result

    def __radd__(self, other):
        if call_other(self, other):
            try:
                return other.__add__(self)
            except AttributeError:
                pass
        return self.result

    def __sub__(self, other):
        if call_other(self, other):
            try:
                return other.__rsub__(self)
            except AttributeError:
                pass
        return self.result

    def __rsub__(self, other):
        if call_other(self, other):
            try:
                return other.__sub__(self)
            except AttributeError:
                pass
        return self.result

    def __pow__(self, other):
        if call_other(self, other):
            try:
                return other.__rpow__(self)
            except AttributeError:
                pass
        return self.result

    def __rpow__(self, other):
        if call_other(self, other):
            try:
                return other.__pow__(self)
            except AttributeError:
                pass
        return self.result

    def __div__(self, other):
        if call_other(self, other):
            try:
                return other.__rdiv__(self)
            except AttributeError:
                pass
        return self.result

    def __rdiv__(self, other):
        if call_other(self, other):
            try:
                return other.__rdiv__(self)
            except AttributeError:
                pass
        return self.result

class Lower(Higher):

    _op_priority = 5.0
    result = 'low'

def test_mul():
    x = Symbol('x')
    h = Higher()
    l = Lower()
    assert l*h == h*l == 'high'
    assert x*h == h*x == 'high'
    assert l*x == x*l != 'low'

def test_add():
    x = Symbol('x')
    h = Higher()
    l = Lower()
    assert l+h == h+l == 'high'
    assert x+h == h+x == 'high'
    assert l+x == x+l != 'low'

def test_sub():
    x = Symbol('x')
    h = Higher()
    l = Lower()
    assert l-h == h-l == 'high'
    assert x-h == h-x == 'high'
    assert l-x == -(x-l) != 'low'

def test_pow():
    x = Symbol('x')
    h = Higher()
    l = Lower()
    assert l**h == h**l == 'high'
    assert x**h == h**x == 'high'
    assert l**x != 'low'
    assert x**l != 'low'

def test_div():
    x = Symbol('x')
    h = Higher()
    l = Lower()
    assert l/h == h/l == 'high'
    assert x/h == h/x == 'high'
    assert l/x != 'low'
    assert x/l != 'low'
