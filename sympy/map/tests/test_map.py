from sympy import symbols, S, FiniteSet
from sympy.map import Map

a, b, c, x = symbols('a b c x')

def test_nargs():

    class f1(Map):
        def eval(self, *args):
            return
    assert f1().nargs == S.Naturals0

    class f2(Map):
        def eval(self, x, y=3):
            return
    assert f2().nargs == FiniteSet(1, 2)

    class f3(Map):
        def eval(self, x, y):
            return
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

def test_inversemap():

    class f1(Map):
        def eval(self, x):
            return x+1
    assert f1().inv(evaluate=False) == f1().inv(evaluate=True)

    class f2(Map):
        def eval(self, x):
            return 2*x
        def _eval_inverse(self):
            return f2_inv(
                parameters=self.parameters,
                domain=self.codomain, codomain=self.domain
            )
    class f2_inv(Map):
        def eval(self, x):
            return x/2
        def _eval_inverse(self):
            return f2(
                parameters=self.parameters,
                domain=self.codomain, codomain=self.domain
            )
    assert f2().inv() != f2().inv(evaluate=True) == f2_inv()
    assert f2_inv().inv() != f2_inv().inv(evaluate=True) == f2()
    assert f2().inv(evaluate=True)(x, evaluate=True) == f2().inv()(x).doit() == x/2
