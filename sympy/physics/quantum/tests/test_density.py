from sympy import pprint, latex, symbols, S, log
from sympy.matrices.matrices import Matrix
from sympy.core.trace import Tr
from sympy.external import import_module
from sympy.physics.quantum.density import Density, entropy
from sympy.physics.quantum.state import Ket, Bra
from sympy.physics.quantum.qubit import Qubit
from sympy.physics.quantum.qapply import qapply
from sympy.physics.quantum.gate import HadamardGate
from sympy.physics.quantum.represent import represent
from sympy.physics.quantum.dagger import Dagger
from sympy.physics.quantum.cartesian import XKet, PxKet, PxOp, XOp
from sympy.physics.quantum.spin import JzKet, Jz
from sympy.functions import sqrt
from sympy.utilities.pytest import raises
from sympy.physics.quantum.matrixutils import scipy_sparse_matrix


def test_eval_args():
    # check instance created
    assert isinstance(Density([Ket(0), 0.5], [Ket(1), 0.5]), Density)
    assert isinstance(Density([Qubit('00'), 1/sqrt(2)],
                              [Qubit('11'), 1/sqrt(2)]), Density)

    #test if Qubit object type preserved
    d = Density([Qubit('00'), 1/sqrt(2)], [Qubit('11'), 1/sqrt(2)])
    for (state, prob) in d.args:
        assert isinstance(state, Qubit)

    # check for value error, when prob is not provided
    raises(ValueError, lambda: Density([Ket(0)], [Ket(1)]))

def test_doit():
    x,y = symbols('x y')
    d = Density([XKet(),0.5], [PxKet(),0.5])
    assert (0.5*(PxKet()*Dagger(PxKet())) +
            0.5*(XKet()*Dagger(XKet()))) == d.doit()

    # check for kets with expr in them
    d_with_sym = Density([XKet(x*y),0.5], [PxKet(x*y),0.5])
    assert (0.5*(PxKet(x*y)*Dagger(PxKet(x*y))) +
            0.5*(XKet(x*y)*Dagger(XKet(x*y)))) == d_with_sym.doit()

def test_apply_op():
    d = Density([Ket(0), 0.5], [Ket(1), 0.5])
    assert d.apply_op(XOp()) == Density([XOp()*Ket(0), 0.5],
                                          [XOp()*Ket(1), 0.5])

def test_represent():
    x,y = symbols('x y')
    d = Density([XKet(),0.5], [PxKet(),0.5])
    assert (represent(0.5*(PxKet()*Dagger(PxKet()))) +
            represent(0.5*(XKet()*Dagger(XKet())))) == represent(d)

    # check for kets with expr in them
    d_with_sym = Density([XKet(x*y),0.5], [PxKet(x*y),0.5])
    assert (represent(0.5*(PxKet(x*y)*Dagger(PxKet(x*y)))) +
            represent(0.5*(XKet(x*y)*Dagger(XKet(x*y))))) == \
        represent(d_with_sym)

    # check when given explicit basis
    assert (represent(0.5*(XKet()*Dagger(XKet())), basis=PxOp()) +
            represent(0.5*(PxKet()*Dagger(PxKet())), basis=PxOp())) == \
        represent(d, basis=PxOp())

def test_states():
    d = Density([Ket(0), 0.5], [Ket(1), 0.5])
    states = d.states()
    assert states[0] == Ket(0) and states[1] == Ket(1)

def test_probs():
    d = Density([Ket(0), .75], [Ket(1), 0.25])
    probs = d.probs()
    assert probs[0] == 0.75 and probs[1] == 0.25

    #probs can be symbols
    x,y = symbols('x y')
    d = Density([Ket(0), x], [Ket(1), y])
    probs = d.probs()
    assert probs[0] == x and probs[1] == y

def test_get_state():
    x,y = symbols('x y')
    d = Density([Ket(0), x], [Ket(1), y])
    states = (d.get_state(0), d.get_state(1))
    assert states[0] == Ket(0) and states[1] == Ket(1)

def test_get_prob():
    x,y = symbols('x y')
    d = Density([Ket(0), x], [Ket(1), y])
    probs = (d.get_prob(0), d.get_prob(1))
    assert probs[0] == x and probs[1] == y

def test_entropy():
    up = JzKet(S(1)/2,S(1)/2)
    down = JzKet(S(1)/2,-S(1)/2)
    d = Density((up,0.5),(down,0.5))

    # test for density object
    ent = entropy(d)
    assert  ent.real == 0.69314718055994529 and ent.imag == 0

    np = import_module('numpy', min_python_version=(2, 6))
    if np:
        #do this test only if 'numpy' is available on test machine
        mat = represent(d)
        assert isinstance(mat, Matrix) and entropy(mat) == 0.5*log(2)

    scipy = import_module('scipy', __import__kwargs={'fromlist':['sparse']})
    if scipy and np:
        #do this test only if numpy and scipy are available
        mat = represent(d, format="scipy.sparse")
        assert isinstance(mat, scipy_sparse_matrix) and ent.real == \
                      0.69314718055994529 and ent.imag == 0

def test_eval_trace():
    up = JzKet(S(1)/2,S(1)/2)
    down = JzKet(S(1)/2,-S(1)/2)
    d = Density((up,0.5),(down,0.5))

    t = Tr(d)
    assert t.doit() == 1

    #TODO: partial trace
