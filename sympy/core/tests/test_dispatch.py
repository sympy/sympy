from sympy.core.add import Add
from sympy.core.basic import Basic
from sympy.core.expr import Expr
from sympy.core.mul import Mul
from sympy.core.numbers import Integer
from sympy.core.singleton import S
from sympy.core.symbol import symbols
from sympy.utilities.pytest import raises

from sympy.core.dispatch import dispatch


def test_dispatch_on_sympy_objects():
    @dispatch(str)
    def f(a):
        return a*3

    @dispatch(int)
    def f(b):
        return "integer value: {}".format(b)

    @dispatch(int, int)
    def f(a, b):
        return "two integers: {}, {}".format(a, b)

    @dispatch(Basic)
    def g(a):
        return "Basic"

    @dispatch(Expr)
    def g(a):
        return "Expr"

    @dispatch(Integer)
    def g(a):
        return "Integer"

    assert g(S.One) == "Integer"
    assert g(Integer(2)/2) == "Integer"
    x = symbols('x')
    assert g(x) == "Expr"
    # assert g(x**2+1) == "Expr"

    @dispatch(Add)
    def g(a):
        return "Add"

    assert g(x**2+1) == "Add"
    assert f("triple this ") == "triple this triple this triple this "
    assert f(34) == "integer value: 34"
    assert f(3, 444) == "two integers: 3, 444"
    raises(TypeError, lambda: f("her", 32))


def test_intertwined_inheritances():
    class DispTestA(object): pass
    class DispTestB(object): pass
    class DispTestA1(DispTestA): pass
    class DispTestA2(DispTestA): pass
    class DispTestB1(DispTestB): pass
    class DispTestAB(DispTestA, DispTestB): pass
    class DispTestAB1(DispTestAB): pass
    class DispTestAB2(DispTestAB): pass
    class DispTestAB12(DispTestAB1, DispTestAB2): pass

    @dispatch(DispTestAB1, object)
    def dispfu(a, b):
        return (DispTestAB1, object)

    @dispatch(object, DispTestAB1)
    def dispfu(a, b):
        return (object, DispTestAB1)

    @dispatch(DispTestA, DispTestB)
    def dispfu(a, b):
        return (DispTestA, DispTestB)

    @dispatch(DispTestA)
    def dispfu(a):
        return (DispTestA)

    @dispatch(DispTestB)
    def dispfu(a):
        return (DispTestB)

    def assert_list1():
        assert dispfu(DispTestAB1(), DispTestAB1()) == (DispTestAB1, object)
        assert dispfu(DispTestAB(), DispTestAB()) == (DispTestA, DispTestB)
        assert dispfu(DispTestAB1(), DispTestAB1()) == (DispTestAB1, object)
        assert dispfu(DispTestAB()) == (DispTestA)

    assert_list1()

    @dispatch(object)
    def dispfu(a):
        return (object)

    assert_list1()


def test_dispatch_class_methods():
    class A(object):
        @dispatch(Basic, prefix="A")
        def method(self):
            return (Basic,)

        @dispatch(Mul, prefix="A")
        def method(self):
            return (Mul,)

        @dispatch(Add, prefix="A")
        def method(self):
            return (Add,)

    class B(object):
        @dispatch(Add, prefix="B")
        def method(self):
            return (B, Add)
    a = A()
    x, y = symbols('x, y')
    assert a.method(x*y) == (Mul,)
    assert a.method(x+y) == (Add,)
    assert a.method(x**y) == (Basic,)
    b = B()
    assert b.method(x+y) == (B, Add)

test_intertwined_inheritances()
test_dispatch_class_methods()