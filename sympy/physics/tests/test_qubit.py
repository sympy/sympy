from sympy.physics.qubit import *
from sympy.physics.quantum import *
from sympy import symbols, Rational
from sympy.core.numbers import *
from sympy.functions.elementary import *
from sympy.physics.shor import *
from sympy.core.containers import Tuple
import random
x, y = symbols('xy')

epsilon = .000001

def test_Qubit():
    array = [0,0,1,1,0]
    qb = Qubit('00110')
    assert qb.flip(0) == Qubit('00111')
    assert qb.flip(1) == Qubit('00100')
    assert qb.flip(4) == Qubit('10110')
    assert qb.dimension == 5
    for i in range(5):
        assert qb[i] == array[4-i]
    assert len(qb) == 5
    qb = Qubit('110')
    """assert qb._represent_ZGate(0)**(ZGate(0)**(3)) == \
    Matrix([0,0,0,0,0,0,1,0])
    """

def test_Gate():
    c = CNOTGate(0,3)
    t = ToffoliGate(0,1,6)
    h = HadamardGate(2)
    assert c.minimum_dimension == 3
    assert t.minimum_dimension == 6
    assert h.minimum_dimension == 2
    assert c.input_number == 2
    assert t.input_number == 3
    assert h.input_number == 1

def test_Fourier():
    assert QFT(0,3).decompose() == SwapGate(0,2)*HadamardGate(0)\
    *RkGate(1,0,2)*HadamardGate(1)*RkGate(2,0,3)*RkGate(2,1,2)*HadamardGate(2)
    assert QFT(0,3).input_number == 2
    assert IQFT(0,3).decompose() == HadamardGate(2)*IRkGate(2,1,2)\
    *IRkGate(2,0,3)*HadamardGate(1)*IRkGate(1,0,2)*HadamardGate(0)*SwapGate(0,2)

def test_represent_HilbertSpace():
    import numpy as np
    a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p = symbols('abcdefghijklmnop')
    gateMat = Matrix([[a,b,c,d],[e,f,g,h],[i,j,k,l],[m,n,o,p]])
    assert represent_hilbert_space(gateMat, 3, (0,1)) == \
    Matrix([[a,c,b,d,0,0,0,0],[i,k,j,l,0,0,0,0],[e,g,f,h,0,0,0,0],\
    [m,o,n,p,0,0,0,0],[0,0,0,0,a,c,b,d],[0,0,0,0,i,k,j,l],\
    [0,0,0,0,e,g,f,h],[0,0,0,0,m,o,n,p]])
    assert type(represent_hilbert_space(gateMat, 2, \
    (0,1), format = 'numpy')) == type(np.matrix(1))

def test_represent_Hadamard_():
    circuit = HadamardGate(0)*Qubit('00')
    answer = represent(circuit, ZGate(0)**2)
    # check that the answers are same to within an epsilon
    assert answer == Matrix([1/sqrt(2),1/sqrt(2), 0, 0])

def test_represent_XGate_():
    circuit = XGate(0)*Qubit('00')
    answer = represent(circuit, ZGate(0)**2)
    assert Matrix([0, 1, 0, 0]) == answer

def test_represent_YGate_():
    circuit = YGate(0)*Qubit('00')
    answer = represent(circuit, ZGate(0)**2)
    assert answer[0] == 0 and answer[1] == ImaginaryUnit() and \
    answer[2] == 0 and answer[3] == 0

def test_represent_ZGate_():
    circuit = ZGate(0)*Qubit('00')
    answer = represent(circuit, ZGate(0)**2)
    assert Matrix([1, 0, 0, 0]) == answer

def test_represent_PhaseGate_():
    circuit = PhaseGate(0)*Qubit('01')
    answer = represent(circuit, ZGate(0)**2)
    assert Matrix([0, ImaginaryUnit(),0,0]) == answer

def test_represent_TGate_():
    circuit = TGate(0)*Qubit('01')
    assert Matrix([0, exp(I*Pi()/4), 0, 0]) == represent(circuit, ZGate(0)**2)

def test_CompoundGates_():
    circuit = YGate(0)*ZGate(0)*XGate(0)*HadamardGate(0)*Qubit('00')
    answer = represent(circuit, ZGate(0)**2)
    assert Matrix([.5*ImaginaryUnit()*sqrt(2),ImaginaryUnit()/sqrt(2), 0, 0])\
    == answer

def test_CNOTGate():
    circuit = CNOTGate(1,0)
    assert represent(circuit, ZGate(0)**2) == \
    Matrix([[1,0,0,0],[0,1,0,0],[0,0,0,1],[0,0,1,0]])
    circuit = circuit*Qubit('111')
    assert matrix_to_Qubits(represent(circuit, ZGate(0)**3)) == \
    apply_operators(circuit)

def test_ToffoliGate():
    circuit = ToffoliGate(2,1,0)
    assert represent(circuit, ZGate(0)**3) == \
    Matrix([[1,0,0,0,0,0,0,0],[0,1,0,0,0,0,0,0],[0,0,1,0,0,0,0,0],\
    [0,0,0,1,0,0,0,0],[0,0,0,0,1,0,0,0],[0,0,0,0,0,1,0,0],[0,0,0,0,0,0,0,1],\
    [0,0,0,0,0,0,1,0]])

    circuit = ToffoliGate(3,0,1)
    assert apply_operators(circuit*Qubit('1001')) == \
    matrix_to_Qubits(represent(circuit*Qubit('1001'), ZGate(0)**4))
    assert apply_operators(circuit*Qubit('0000')) == \
    matrix_to_Qubits(represent(circuit*Qubit('0000'), ZGate(0)**4))

def test_SwapGate():
    assert apply_operators(SwapGate(0,1)*Qubit('10')) == Qubit('01')
    assert Qubit('010') == apply_operators(SwapGate(1,0)*SwapGate(0,1)*Qubit('010'))
    assert matrix_to_Qubits(represent(SwapGate(0,1)*Qubit('10'), ZGate(0)**2))\
     == Qubit('01')
    assert Qubit('010') == matrix_to_Qubits(represent(SwapGate(1,0)\
    *SwapGate(0,1)*Qubit('010'), ZGate(0)**3))

def test_ControlledZ_Gate():
    assert apply_operators(CZGate(0,1)*Qubit('11')) == -Qubit('11')
    assert matrix_to_Qubits(represent(CZGate(0,1)*Qubit('11'),\
     ZGate(0)**2)) == -Qubit('11')

def test_CPhase_Gate():
    assert apply_operators(CPhaseGate(0,1)*Qubit('11')) == ImaginaryUnit()*Qubit('11')
    assert matrix_to_Qubits(represent(CPhaseGate(0,1)*Qubit('11'),\
     ZGate(0)**2)) == ImaginaryUnit()*Qubit('11')

def test_gateSort():
    assert gate_sort(XGate(1)*HadamardGate(0)**2*CNOTGate(0,1)*XGate(1)*XGate(0))\
     == HadamardGate(0)**2*XGate(1)*CNOTGate(0,1)*XGate(0)*XGate(1)

def test_gate_simp():
     assert gate_simp(HadamardGate(0)*XGate(1)*HadamardGate(0)**2*CNOTGate(0,1)\
     *XGate(1)**3*XGate(0)*ZGate(3)**2*PhaseGate(4)**3) == HadamardGate(0)*\
     XGate(1)*CNOTGate(0,1)*XGate(0)*XGate(1)*ZGate(4)*PhaseGate(4)

def test_gate_Qubit_strings():
    assert sstr(Qubit('01')) == "|01>"
    assert sstr(HadamardGate(3)) == "H(3)"
    assert sstr(XGate(2)) == "NOT(2)"
    assert sstr(ZGate(6)) == "Z(6)"
    assert sstr(YGate(6)) == "Y(6)"
    assert sstr(CNOTGate(1,0)) == "CNOT(1,0)"

def test_ArbMat4_apply():
    a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p = symbols('abcdefghijklmnop')
    class Arb(Gate):
        @property
        def matrix(self):
            return Matrix([[a,b,c,d],[e,f,g,h],[i,j,k,l],[m,n,o,p]])

    assert apply_operators(Arb(1,0)*Qubit('00101')) == b*Qubit('00100')\
     + f*Qubit('00101') + j*Qubit('00110') + n*Qubit('00111')
    assert apply_operators(Arb(2,4)*Qubit('00101')) == c*Qubit('00001')\
     + g*Qubit('10001') + k*Qubit('00101') + o*Qubit('10101')
    assert apply_operators(Arb(3,0)*Qubit('11111')) == d*Qubit('10110')\
     + h*Qubit('10111') + l*Qubit('11110') + p*Qubit('11111')
    assert apply_operators(Arb(6,9)*Qubit('0110110101')) ==\
     a*Qubit('0110110101') + e*Qubit('1110110101') +\
      i*Qubit('0111110101') + m*Qubit('1111110101')

def test_ArbMat8_apply():
    a,b,c,d,e,f,g,h = symbols('abcdefgh')
    class Arb(Gate):
        @property
        def matrix(self):
            symlist = [a,b,c,d,e,f,g,h]
            lout = []
            for i in range(8):
                lin = []
                for j in range(8):
                    lin.append(symlist[i]**j)
                lout.append(lin)
            return Matrix(lout)

    assert apply_operators(Arb(2,1,0)*Qubit('01101')) == \
    a**5*Qubit('01000') + b**5*Qubit('01001') + c**5*Qubit('01010') +\
     d**5*Qubit('01011') + e**5*Qubit('01100') + f**5*Qubit('01101') +\
      g**5*Qubit('01110') + h**5*Qubit('01111')

    assert apply_operators(Arb(0,4,3)*Qubit('11010')) == \
    a**3*Qubit('00010') + b**3*Qubit('01010') + c**3*Qubit('10010')\
     + d**3*Qubit('11010') + e**3*Qubit('00011') + f**3*Qubit('01011')\
      + g**3*Qubit('10011') + h**3*Qubit('11011')

    assert apply_operators(Arb(4,1,3)*Qubit('010010')) ==\
     a**6*Qubit('000000') + b**6*Qubit('001000') + c**6*Qubit('000010') \
     + d**6*Qubit('001010') + e**6*Qubit('010000') + f**6*Qubit('011000')\
     + g**6*Qubit('010010') + h**6*Qubit('011010')

    assert apply_operators(Arb(3,1,4)*Qubit('010100101')) == \
    Qubit('010100101') + Qubit('010110101') + Qubit('010100111')\
     + Qubit('010110111') + Qubit('010101101') +\
     Qubit('010111101') + Qubit('010101111') + Qubit('010111111')

    assert apply_operators(Arb(8,10,9)*Qubit('11101010101')) ==\
    a**7*Qubit('00001010101') + b**7*Qubit('01001010101')\
    + c**7*Qubit('10001010101') + d**7*Qubit('11001010101')\
    + e**7*Qubit('00101010101') + f**7*Qubit('01101010101') +\
     g**7*Qubit('10101010101') + h**7*Qubit('11101010101')

    assert apply_operators(Arb(9,2,3)*Qubit('0111111011'))\
     == a*Qubit('0111110011') + b*Qubit('0111111011')\
      + c*Qubit('0111110111') + d*Qubit('0111111111')\
       + e*Qubit('1111110011') + f*Qubit('1111111011')\
        + g*Qubit('1111110111') + h*Qubit('1111111111')

    assert apply_operators(Arb(2,1,0)*Qubit('010')) ==\
    a**2*Qubit('000') + b**2*Qubit('001') + c**2*Qubit('010') + d**2*Qubit('011')\
    + e**2*Qubit('100') + f**2*Qubit('101') + g**2*Qubit('110') + h**2*Qubit('111')

def test_ArbMat4_Equality():

    class Arb(Gate):
        @property
        def matrix(self):
            a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p = symbols('abcdefghijklmnop')
            return Matrix([[a,b,c,d],[e,f,g,h],[i,j,k,l],[m,n,o,p]])

    for i in range(4):
        for j in range(4):
            if j != i:
                assert apply_operators(Arb(i,j)*(Qubit('10110'))) ==\
                matrix_to_Qubits(represent(Arb(i,j)*Qubit('10110'),\
                ZGate(0)**5))

def test_Arb8_Matrix_Equality():
    class Arb(Gate):
        @property
        def matrix(self):
            a,b,c,d,e,f,g,h = symbols('abcdefgh')
            symlist = [a,b,c,d,e,f,g,h]
            lout = []
            for i in range(8):
                lin = []
                for j in range(8):
                    lin.append(symlist[i]**j)
                lout.append(lin)
            return Matrix(lout)

    for i in range(1):
        for j in range(4):
            for k in range(4):
                if j != i and k != i and k != j:
                    assert apply_operators(Arb(i,j,k)*(Qubit((0,1,1,1,0)))) ==\
                     matrix_to_Qubits(represent(Arb(i,j,k)*Qubit((0,1,1,1,0)),\
                     ZGate(0)**5))

def test_superposition_of_states():
    assert apply_operators(CNOTGate(0,1)*HadamardGate(0)*(1/sqrt(2)*Qubit((0,1))\
     + 1/sqrt(2)*Qubit((1,0)))).expand() == (Qubit((0,1))/2 + Qubit((0,0))/2 - Qubit((1,1))/2 +\
     Qubit((1,0))/2)

    assert matrix_to_Qubits(represent(CNOTGate(0,1)*HadamardGate(0)\
    *(1/sqrt(2)*Qubit((0,1)) + 1/sqrt(2)*Qubit((1,0))), ZGate(0)**(2)))\
     == (Qubit((0,1))/2 + Qubit((0,0))/2 - Qubit((1,1))/2 + Qubit((1,0))/2)

def test_tensor_product():
    try:
        import numpy as np
    except ImportError:
        print 'import error numpy'
        return
    l1 = zeros(4)
    for i in range(16):
        l1[i] = 2**i
    l2 = zeros(4)
    for i in range(16):
        l2[i] = i
    l3 = zeros(2)
    for i in range(4):
        l3[i] = i
    vec = Matrix([1,2,3])

    #test for Matrix known 4x4 matricies
    numpyl1 = np.matrix(l1.tolist())
    numpyl2 = np.matrix(l2.tolist())
    numpy_product = np.kron(numpyl1,numpyl2)
    args = [l1, l2]
    sympy_product = tensor_product(*args)
    assert numpy_product.tolist() == sympy_product.tolist()
    numpy_product = np.kron(numpyl2,numpyl1)
    args = [l2, l1]
    sympy_product = tensor_product(*args)
    assert numpy_product.tolist() == sympy_product.tolist()

    #test for other known matrix of different dimensions
    numpyl2 = np.matrix(l3.tolist())
    numpy_product = np.kron(numpyl1,numpyl2)
    args = [l1, l3]
    sympy_product = tensor_product(*args)
    assert numpy_product.tolist() == sympy_product.tolist()
    numpy_product = np.kron(numpyl2,numpyl1)
    args = [l3, l1]
    sympy_product = tensor_product(*args)
    assert numpy_product.tolist() == sympy_product.tolist()

    #test for non square matrix
    numpyl2 = np.matrix(vec.tolist())
    numpy_product = np.kron(numpyl1,numpyl2)
    args = [l1, vec]
    sympy_product = tensor_product(*args)
    assert numpy_product.tolist() == sympy_product.tolist()
    numpy_product = np.kron(numpyl2,numpyl1)
    args = [vec, l1]
    sympy_product = tensor_product(*args)
    assert numpy_product.tolist() == sympy_product.tolist()

    #test for random matrix with random values that are floats
    random_matrix1 = np.random.rand(np.random.rand()*5+1,np.random.rand()*5+1)
    random_matrix2 = np.random.rand(np.random.rand()*5+1,np.random.rand()*5+1)
    numpy_product = np.kron(random_matrix1,random_matrix2)
    args = [Matrix(random_matrix1.tolist()),Matrix(random_matrix2.tolist())]
    sympy_product = tensor_product(*args)
    assert not (sympy_product - Matrix(numpy_product.tolist())).tolist() > \
    (ones((sympy_product.rows,sympy_product.cols))*epsilon).tolist()

    #test for three matrix kronecker
    sympy_product = tensor_product(l1,vec,l2)
    npl1 = np.matrix(l1.tolist())
    npl2 = np.matrix(l2.tolist())
    npvec = np.matrix(vec.tolist())

    numpy_product = np.kron(l1,np.kron(vec,l2))
    assert numpy_product.tolist() == sympy_product.tolist()

#test apply methods
def test_apply_represent_equality():
    gates = [HadamardGate(int(3*random.random())),\
     XGate(int(3*random.random())), ZGate(int(3*random.random())),\
      YGate(int(3*random.random())), ZGate(int(3*random.random())),\
       PhaseGate(int(3*random.random()))]

    circuit = Qubit((int(random.random()*2),int(random.random()*2),\
    int(random.random()*2),int(random.random()*2),int(random.random()*2),\
    int(random.random()*2)))
    for i in range(int(random.random()*6)):
        circuit = gates[int(random.random()*6)]*circuit


    mat = represent(circuit, ZGate(0)**(6))
    states = apply_operators(circuit)
    state_rep = matrix_to_Qubits(mat)
    states = states.expand()
    state_rep = state_rep.expand()
    assert state_rep == states

def test_matrix_to_Qubits():
    assert matrix_to_Qubits(Matrix([1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]))\
    == Qubit((0,0,0,0))
    assert Qubits_to_matrix(Qubit((0,0,0,0))) ==\
    Matrix([1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0])
    assert matrix_to_Qubits(sqrt(2)*2*Matrix([1,1,1,1,1,1,1,1])) ==\
    (2*sqrt(2)*(Qubit((0,0,0)) + Qubit((0,0,1)) + Qubit((0,1,0)) + Qubit((0,1,1))\
    + Qubit((1,0,0)) + Qubit((1,0,1)) + Qubit((1,1,0)) + Qubit((1,1,1)))).expand()
    assert Qubits_to_matrix(2*sqrt(2)*(Qubit((0,0,0)) + Qubit((0,0,1)) + Qubit((0,1,0))\
    + Qubit((0,1,1)) + Qubit((1,0,0)) + Qubit((1,0,1)) + Qubit((1,1,0)) + Qubit((1,1,1))))\
    == sqrt(2)*2*Matrix([1,1,1,1,1,1,1,1])

def test_RkGate_and_inverse():
    assert RkGate(1,2,x).k == x
    assert RkGate(1,2,x).args[0] == Tuple(1,2)
    assert IRkGate(1,2,x).k == x
    assert IRkGate(1,2,x).args[0] == Tuple(1,2)

    assert represent(RkGate(0,1,2), ZGate(0)**2) ==\
    Matrix([[1,0,0,0], [0,1,0,0], [0,0,1,0], [0,0,0,\
    exp(2*ImaginaryUnit()*Pi()/2**2)]])

    assert represent(IRkGate(0,1,3), ZGate(0)**(2)) ==\
    Matrix([[1,0,0,0], [0,1,0,0], [0,0,1,0], [0,0,0,\
    exp(-2*ImaginaryUnit()*Pi()/2**3)]])

def test_quantum_fourier():
    assert QFT(0,3).decompose() == SwapGate(0,2)*HadamardGate(0)*RkGate(1,0,2)\
    *HadamardGate(1)*RkGate(2,0,3)*RkGate(2,1,2)*HadamardGate(2)
    assert IQFT(0,3).decompose() == HadamardGate(2)*IRkGate(2,1,2)*IRkGate(2,0,3)\
    * HadamardGate(1)*IRkGate(1,0,2)*HadamardGate(0)*SwapGate(0,2)
    assert represent(QFT(0,3).decompose()*IQFT(0,3).decompose(), ZGate(0)**(3))\
     == eye(8)
    assert apply_operators(QFT(0,3).decompose()*Qubit((0,0,0))) ==\
     apply_operators(HadamardGate(0)*HadamardGate(1)*HadamardGate(2)*Qubit((0,0,0))).expand()
