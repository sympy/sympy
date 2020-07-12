from sympy import symbols, S, FiniteSet
from sympy.map import Map

a, b, c, x = symbols('a b c x')

def test_nargs():
    class f1(Map):
        def eval(self, *args):
            return
    class f2(Map):
        def eval(self, x, y=3):
            return
    class f3(Map):
        def eval(self, x, y):
            return

    assert f1().nargs == S.Naturals0
    assert f2().nargs == FiniteSet(1, 2)
    assert f3().nargs == FiniteSet(2)

def test_parameters():

    class QuadraticFunction(Map):
        def eval(self, x):
            a, b, c = self.parameters
            return a*x**2 + b*x + c
    f = QuadraticFunction((a, b, c))

    assert f(x, evaluate=True) == f(x).doit() == a*x**2 + b*x + c
    assert f(x).free_symbols == {a, b, c, x}

    assert f.subs(a, 2)(x).doit() == f(x).subs(a, 2).doit() == 2*x**2 + b*x + c
    assert f(x).subs(a, 2).free_symbols == {b, c, x}
