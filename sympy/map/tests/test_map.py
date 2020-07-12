from sympy import symbols, S, FiniteSet
from sympy.map import Map, InverseMap, AppliedMap

a, b, c, x = symbols('a b c x')

def test_nargs():

    assert Map(domain=S.Reals).nargs == 1
    assert Map(domain=S.Reals**1).nargs == 1
    assert Map(domain=S.Reals*S.Integers).nargs == 2
    assert Map(domain=S.Reals**4).nargs == 4

def test_parameters():

    f1 = Map(parameters=(a, b, c))
    assert f1.parameters == (a, b, c)
    assert f1(x).free_symbols == {a, b, c, x}

    f2 = f1.subs(b, 2)
    assert f2.parameters == (a, 2, c)
    assert f2(x).free_symbols == {a, c, x}

def test_inversemap():

    class f1(Map):
        def eval(self, x):
            return x+1
    assert f1().inv() == f1().inv().doit() == f1().inv(evaluate=True) == InverseMap(f1())

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
    assert f2().inv() != f2().inv().doit() == f2().inv(evaluate=True) == f2_inv()
    assert f2_inv().inv() != f2_inv().inv().doit() == f2_inv().inv(evaluate=True) == f2()
