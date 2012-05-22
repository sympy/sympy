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

def test_eval_args():
    pass

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
