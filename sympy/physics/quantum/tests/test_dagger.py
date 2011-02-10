from sympy import I, Matrix, symbols, conjugate, Expr, Integer

from sympy.physics.quantum.dagger import Dagger


def test_scalars():
    x = symbols('x',complex=True)
    assert Dagger(x) == conjugate(x)
    assert Dagger(I*x) == -I*conjugate(x)

    i = symbols('i',real=True)
    assert Dagger(i) == i

    p = symbols('p')
    assert isinstance(Dagger(p), Dagger)

    i = Integer(3)
    assert Dagger(i) == i


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

try:
    import numpy as np
except ImportError:
    pass
else:
    def test_numpy_dagger():
        a = np.matrix([[1.0,2.0j],[-1.0j,2.0]])
        adag = a.copy().transpose().conjugate()
        assert (Dagger(a) == adag).all()


try:
    from scipy import sparse
    import numpy as np
except ImportError:
    pass
else:
    def test_scipy_sparse_dagger():
        a = sparse.csr_matrix([[1.0+0.0j,2.0j],[-1.0j,2.0+0.0j]])
        adag = a.copy().transpose().conjugate()
        assert np.linalg.norm((Dagger(a) - adag).todense()) == 0.0
