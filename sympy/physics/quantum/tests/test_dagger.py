from sympy import I, Matrix, symbols, conjugate, Expr

from sympy.physics.quantum.dagger import Dagger


def test_scalars():
    x,y,z = symbols('xyz')
    i,j,k = symbols('ijk',real=True)
    assert Dagger(x) == conjugate(x)
    assert Dagger(i) == i
    assert Dagger(I*x) == -I*conjugate(x)


def test_matrix():
    x = symbols('x')
    m = Matrix([[I,x*I],[2,4]])
    assert Dagger(m) == m.H


class Foo(Expr):

    def _eval_dagger(self):
        return I


def test_eval_dagger():
    f = Foo()
    d = Dagger(f)
    assert d == I

