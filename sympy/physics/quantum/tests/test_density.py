from sympy import pprint, latex
from sympy.physics.quantum.density import Density
from sympy.physics.quantum.state import Ket, Bra
from sympy.physics.quantum.qubit import Qubit
from sympy.physics.quantum.qapply import qapply
from sympy.physics.quantum.gate import HadamardGate
from sympy.physics.quantum.represent import represent
from sympy.physics.quantum.dagger import Dagger
from sympy.physics.quantum.operator import *
from sympy.functions import sqrt
from sympy.utilities.pytest import raises

def test_eval_args():
    # check instance created
    assert(isinstance(Density([Ket(0), 0.5], [Ket(1), 0.5]), Density))

    # check for value error, when prob is not provided
    raises(ValueError, 'Density([Ket(0)], [Ket(1)])')

    #check for valid state
    raises(ValueError, 'Density(1,1)')
    raises(ValueError, 'Density([Ket(0), 0.5], (Ket(1), 0.25), (1, 0.25))')

    #TODO: Need to implement Qubit based Density before
    #this test is done.
    #assert(isinstance(Density([Qubit('01'),0.5], [Qubit('01'),0.5]), Density))


def test_doit():
    pass

def test_represent():
    pass

def test_entropy():
    pass

def entropy():
    pass

def test_reduced_density():
    pass

def test_entropy_of_entanglement():
    pass

def test_latex():
    d = Density([Ket(0),0.5],[Ket(1),0.5])
    result =  (r'\rho\left(\begin{pmatrix}{\left|0\right\rangle }, & '
           r'0.5\end{pmatrix},\begin{pmatrix}{\left|1\right\rangle }, & '
           r'0.5\end{pmatrix}\right)')
    assert latex(d) == result


if __name__ == '__main__':
    #test_latex()
    #d = Density([Ket(),0.5],[Ket(),0.5])
    test_eval_args()
    #print d
    #pprint(d)
    #print latex(d)
    #pprint(d)
    #k = Ket()
    #pprint(k)
    #print k


    #o = OuterProduct(Ket(),Bra())
    #print latex(o)

    #k = Ket('k')
    #b = Bra('b')
    #op = OuterProduct(k, b)
    #print latex(op)

#def test_representDensity():
#    assert represent(Density([Qubit('01'),1]), nqubits=2) == mat*Dagger(mat)

#def test_apply_operators_density():
#    assert apply_operators(Density([Qubit('00'),1]).operate_on(HadamardGate(0)\
#           *HadamardGate(1)))\
#        == Density([apply_operators(HadamardGate(0)*HadamardGate(1)*Qubit('00')),1])

