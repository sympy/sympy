from sympy.core.mul import Mul
from sympy.core.numbers import (I, Integer, Rational, pi)
from sympy.core.symbol import (Wild, symbols)
from sympy.functions.elementary.exponential import exp
from sympy.functions.elementary.miscellaneous import sqrt
from sympy.matrices import Matrix, ImmutableMatrix

from sympy.physics.quantum.gate import (XGate, YGate, ZGate, random_circuit,
        CNOT, IdentityGate, H, X, Y, S, T, Z, SwapGate, gate_simp, gate_sort,
        CNotGate, TGate, HadamardGate, PhaseGate, UGate, CGate,
        SdgGate, TdgGate, SXGate, PhaseShiftGate, RXGate, RYGate, RZGate,
        Sdg, Tdg, SX, P, RX, RY, RZ)
from sympy.physics.quantum.commutator import Commutator
from sympy.physics.quantum.anticommutator import AntiCommutator
from sympy.physics.quantum.represent import represent
from sympy.physics.quantum.qapply import qapply
from sympy.physics.quantum.qubit import Qubit, IntQubit, qubit_to_matrix, \
    matrix_to_qubit
from sympy.physics.quantum.matrixutils import matrix_to_zero
from sympy.physics.quantum.matrixcache import sqrt2_inv
from sympy.physics.quantum import Dagger


def test_gate():
    """Test a basic gate."""
    h = HadamardGate(1)
    assert h.min_qubits == 2
    assert h.nqubits == 1

    i0 = Wild('i0')
    i1 = Wild('i1')
    h0_w1 = HadamardGate(i0)
    h0_w2 = HadamardGate(i0)
    h1_w1 = HadamardGate(i1)

    assert h0_w1 == h0_w2
    assert h0_w1 != h1_w1
    assert h1_w1 != h0_w2

    cnot_10_w1 = CNOT(i1, i0)
    cnot_10_w2 = CNOT(i1, i0)
    cnot_01_w1 = CNOT(i0, i1)

    assert cnot_10_w1 == cnot_10_w2
    assert cnot_10_w1 != cnot_01_w1
    assert cnot_10_w2 != cnot_01_w1


def test_UGate():
    a, b, c, d = symbols('a,b,c,d')
    uMat = Matrix([[a, b], [c, d]])

    # Test basic case where gate exists in 1-qubit space
    u1 = UGate((0,), uMat)
    assert represent(u1, nqubits=1) == uMat
    assert qapply(u1*Qubit('0')) == a*Qubit('0') + c*Qubit('1')
    assert qapply(u1*Qubit('1')) == b*Qubit('0') + d*Qubit('1')

    # Test case where gate exists in a larger space
    u2 = UGate((1,), uMat)
    u2Rep = represent(u2, nqubits=2)
    for i in range(4):
        assert u2Rep*qubit_to_matrix(IntQubit(i, 2)) == \
            qubit_to_matrix(qapply(u2*IntQubit(i, 2)))


def test_cgate():
    """Test the general CGate."""
    # Test single control functionality
    CNOTMatrix = Matrix(
        [[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 0, 1], [0, 0, 1, 0]])
    assert represent(CGate(1, XGate(0)), nqubits=2) == CNOTMatrix

    # Test multiple control bit functionality
    ToffoliGate = CGate((1, 2), XGate(0))
    assert represent(ToffoliGate, nqubits=3) == \
        Matrix(
            [[1, 0, 0, 0, 0, 0, 0, 0], [0, 1, 0, 0, 0, 0, 0, 0], [0, 0, 1, 0, 0, 0, 0, 0],
    [0, 0, 0, 1, 0, 0, 0, 0], [0, 0, 0, 0, 1, 0, 0, 0], [0, 0, 0, 0, 0,
        1, 0, 0], [0, 0, 0, 0, 0, 0, 0, 1],
    [0, 0, 0, 0, 0, 0, 1, 0]])

    ToffoliGate = CGate((3, 0), XGate(1))
    assert qapply(ToffoliGate*Qubit('1001')) == \
        matrix_to_qubit(represent(ToffoliGate*Qubit('1001'), nqubits=4))
    assert qapply(ToffoliGate*Qubit('0000')) == \
        matrix_to_qubit(represent(ToffoliGate*Qubit('0000'), nqubits=4))

    CYGate = CGate(1, YGate(0))
    CYGate_matrix = Matrix(
        ((1, 0, 0, 0), (0, 1, 0, 0), (0, 0, 0, -I), (0, 0, I, 0)))
    # Test 2 qubit controlled-Y gate decompose method.
    assert represent(CYGate.decompose(), nqubits=2) == CYGate_matrix

    CZGate = CGate(0, ZGate(1))
    CZGate_matrix = Matrix(
        ((1, 0, 0, 0), (0, 1, 0, 0), (0, 0, 1, 0), (0, 0, 0, -1)))
    assert qapply(CZGate*Qubit('11')) == -Qubit('11')
    assert matrix_to_qubit(represent(CZGate*Qubit('11'), nqubits=2)) == \
        -Qubit('11')
    # Test 2 qubit controlled-Z gate decompose method.
    assert represent(CZGate.decompose(), nqubits=2) == CZGate_matrix

    CPhaseGate = CGate(0, PhaseGate(1))
    assert qapply(CPhaseGate*Qubit('11')) == \
        I*Qubit('11')
    assert matrix_to_qubit(represent(CPhaseGate*Qubit('11'), nqubits=2)) == \
        I*Qubit('11')

    # Test that the dagger, inverse, and power of CGate is evaluated properly
    assert Dagger(CZGate) == CZGate
    assert pow(CZGate, 1) == Dagger(CZGate)
    assert Dagger(CZGate) == CZGate.inverse()
    assert Dagger(CPhaseGate) != CPhaseGate
    assert Dagger(CPhaseGate) == CPhaseGate.inverse()
    assert Dagger(CPhaseGate) == pow(CPhaseGate, -1)
    assert pow(CPhaseGate, -1) == CPhaseGate.inverse()


def test_UGate_CGate_combo():
    a, b, c, d = symbols('a,b,c,d')
    uMat = Matrix([[a, b], [c, d]])
    cMat = Matrix([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, a, b], [0, 0, c, d]])

    # Test basic case where gate exists in 1-qubit space.
    u1 = UGate((0,), uMat)
    cu1 = CGate(1, u1)
    assert represent(cu1, nqubits=2) == cMat
    assert qapply(cu1*Qubit('10')) == a*Qubit('10') + c*Qubit('11')
    assert qapply(cu1*Qubit('11')) == b*Qubit('10') + d*Qubit('11')
    assert qapply(cu1*Qubit('01')) == Qubit('01')
    assert qapply(cu1*Qubit('00')) == Qubit('00')

    # Test case where gate exists in a larger space.
    u2 = UGate((1,), uMat)
    u2Rep = represent(u2, nqubits=2)
    for i in range(4):
        assert u2Rep*qubit_to_matrix(IntQubit(i, 2)) == \
            qubit_to_matrix(qapply(u2*IntQubit(i, 2)))

def test_UGate_OneQubitGate_combo():
    v, w, f, g = symbols('v w f g')
    uMat1 = ImmutableMatrix([[v, w], [f, g]])
    cMat1 = Matrix([[v, w + 1, 0, 0], [f + 1, g, 0, 0], [0, 0, v, w + 1], [0, 0, f + 1, g]])
    u1 = X(0) + UGate(0, uMat1)
    assert represent(u1, nqubits=2) == cMat1

    uMat2 = ImmutableMatrix([[1/sqrt(2), 1/sqrt(2)], [I/sqrt(2), -I/sqrt(2)]])
    cMat2_1 = Matrix([[Rational(1, 2) + I/2, Rational(1, 2) - I/2],
                      [Rational(1, 2) - I/2, Rational(1, 2) + I/2]])
    cMat2_2 = Matrix([[1, 0], [0, I]])
    u2 = UGate(0, uMat2)
    assert represent(H(0)*u2, nqubits=1) == cMat2_1
    assert represent(u2*H(0), nqubits=1) == cMat2_2

def test_represent_hadamard():
    """Test the representation of the hadamard gate."""
    circuit = HadamardGate(0)*Qubit('00')
    answer = represent(circuit, nqubits=2)
    # Check that the answers are same to within an epsilon.
    assert answer == Matrix([sqrt2_inv, sqrt2_inv, 0, 0])


def test_represent_xgate():
    """Test the representation of the X gate."""
    circuit = XGate(0)*Qubit('00')
    answer = represent(circuit, nqubits=2)
    assert Matrix([0, 1, 0, 0]) == answer


def test_represent_ygate():
    """Test the representation of the Y gate."""
    circuit = YGate(0)*Qubit('00')
    answer = represent(circuit, nqubits=2)
    assert answer[0] == 0 and answer[1] == I and \
        answer[2] == 0 and answer[3] == 0


def test_represent_zgate():
    """Test the representation of the Z gate."""
    circuit = ZGate(0)*Qubit('00')
    answer = represent(circuit, nqubits=2)
    assert Matrix([1, 0, 0, 0]) == answer


def test_represent_phasegate():
    """Test the representation of the S gate."""
    circuit = PhaseGate(0)*Qubit('01')
    answer = represent(circuit, nqubits=2)
    assert Matrix([0, I, 0, 0]) == answer


def test_represent_tgate():
    """Test the representation of the T gate."""
    circuit = TGate(0)*Qubit('01')
    assert Matrix([0, exp(I*pi/4), 0, 0]) == represent(circuit, nqubits=2)


def test_represent_sdggate():
    """Test the representation of the Sdg gate."""
    circuit = SdgGate(0)*Qubit('01')
    answer = represent(circuit, nqubits=2)
    assert Matrix([0, -I, 0, 0]) == answer


def test_represent_tdggate():
    """Test the representation of the Tdg gate."""
    circuit = TdgGate(0)*Qubit('01')
    assert Matrix([0, exp(-I*pi/4), 0, 0]) == represent(circuit, nqubits=2)


def test_represent_sxgate():
    """Test the representation of the SX gate."""
    circuit = SXGate(0)*Qubit('01')
    answer = represent(circuit, nqubits=2)
    # SX applied to qubit 0 which is in state |1>
    # SX|1> = ((1-i)/2)|0> + ((1+i)/2)|1>
    # Basis order is |00>, |01>, |10>, |11>
    # Qubit('01') is |01> where qubit 1 is |0> and qubit 0 is |1>
    # Applying SX to qubit 0: (|0> ⊗ SX)|01> = |0> ⊗ (SX|1>)
    # = |0> ⊗ [((1-i)/2)|0> + ((1+i)/2)|1>]
    # = ((1-i)/2)|00> + ((1+i)/2)|01>
    assert Matrix([(1-I)/2, (1+I)/2, 0, 0]) == answer


def test_represent_phaseshiftgate():
    """Test the representation of the parameterized phase shift gate."""
    # Test with pi/2 (should be same as S gate)
    circuit = PhaseShiftGate(0, pi/2)*Qubit('01')
    answer = represent(circuit, nqubits=2)
    assert Matrix([0, I, 0, 0]) == answer

    # Test with pi/4 (should be same as T gate)
    circuit = PhaseShiftGate(0, pi/4)*Qubit('01')
    answer = represent(circuit, nqubits=2)
    assert Matrix([0, exp(I*pi/4), 0, 0]) == answer


def test_represent_rxgate():
    """Test the representation of the RX gate."""
    # Test with pi (should flip |0> to |1> similar to X gate)
    circuit = RXGate(0, pi)*Qubit('00')
    answer = represent(circuit, nqubits=2)
    # RX(pi)|0> = -i|1>
    # Qubit('00') means both qubits are |0>
    # Applying RX(pi) to qubit 0: |0> ⊗ RX(pi)|0> = |0> ⊗ (-i|1>) = -i|01>
    assert Matrix([0, -I, 0, 0]) == answer


def test_represent_rygate():
    """Test the representation of the RY gate."""
    # Test with pi/2
    circuit = RYGate(0, pi/2)*Qubit('00')
    answer = represent(circuit, nqubits=2)
    # RY(pi/2)|0> = (|0> + |1>)/sqrt(2)
    assert Matrix([1/sqrt(2), 1/sqrt(2), 0, 0]) == answer


def test_represent_rzgate():
    """Test the representation of the RZ gate."""
    # Test with pi (should be similar to Z gate up to global phase)
    circuit = RZGate(0, pi)*Qubit('01')
    answer = represent(circuit, nqubits=2)
    # RZ(pi)|1> = exp(i*pi/2)|1> = i|1>
    assert Matrix([0, I, 0, 0]) == answer


def test_sdggate_qapply():
    """Test SdgGate gate with qapply on single-qubit states."""
    # Test with |0> state
    assert qapply(SdgGate(0)*Qubit('0')) == Qubit('0')

    # Test with |1> state
    assert qapply(SdgGate(0)*Qubit('1')) == -I*Qubit('1')


def test_tdggate_qapply():
    """Test TdgGate gate with qapply on single-qubit states."""
    # Test with |0> state
    assert qapply(TdgGate(0)*Qubit('0')) == Qubit('0')

    # Test with |1> state
    assert qapply(TdgGate(0)*Qubit('1')) == exp(-I*pi/4)*Qubit('1')


def test_sxgate_qapply():
    """Test SXGate gate with qapply on single-qubit states."""
    # Test with |0> state
    # SX|0> = (1+I)/2|0> + (1-I)/2|1> = 1/2|0> + I/2|0> + 1/2|1> - I/2|1>
    assert qapply(SXGate(0)*Qubit('0')) == Qubit('0')/2 + I*Qubit('0')/2 + Qubit('1')/2 - I*Qubit('1')/2

    # Test with |1> state
    # SX|1> = (1-I)/2|0> + (1+I)/2|1> = 1/2|0> - I/2|0> + 1/2|1> + I/2|1>
    assert qapply(SXGate(0)*Qubit('1')) == Qubit('0')/2 - I*Qubit('0')/2 + Qubit('1')/2 + I*Qubit('1')/2


def test_phaseshiftgate_qapply():
    """Test PhaseShiftGate gate with qapply on single-qubit states."""
    from sympy.core.symbol import symbols
    theta = symbols('theta', real=True)

    # Test with symbolic angle
    assert qapply(PhaseShiftGate(0, theta)*Qubit('0')) == Qubit('0')
    assert qapply(PhaseShiftGate(0, theta)*Qubit('1')) == exp(I*theta)*Qubit('1')

    # Test with concrete angle value pi/2
    assert qapply(PhaseShiftGate(0, pi/2)*Qubit('0')) == Qubit('0')
    assert qapply(PhaseShiftGate(0, pi/2)*Qubit('1')) == I*Qubit('1')


def test_rxgate_qapply():
    """Test RXGate gate with qapply on single-qubit states."""
    from sympy.core.symbol import symbols
    from sympy.functions.elementary.trigonometric import cos, sin
    theta = symbols('theta', real=True)

    # Test with symbolic angle
    assert qapply(RXGate(0, theta)*Qubit('0')) == cos(theta/2)*Qubit('0') - I*sin(theta/2)*Qubit('1')
    assert qapply(RXGate(0, theta)*Qubit('1')) == -I*sin(theta/2)*Qubit('0') + cos(theta/2)*Qubit('1')

    # Test with concrete angle value pi
    assert qapply(RXGate(0, pi)*Qubit('0')) == -I*Qubit('1')
    assert qapply(RXGate(0, pi)*Qubit('1')) == -I*Qubit('0')


def test_rygate_qapply():
    """Test RYGate gate with qapply on single-qubit states."""
    from sympy.core.symbol import symbols
    from sympy.functions.elementary.trigonometric import cos, sin
    theta = symbols('theta', real=True)

    # Test with symbolic angle
    assert qapply(RYGate(0, theta)*Qubit('0')) == cos(theta/2)*Qubit('0') + sin(theta/2)*Qubit('1')
    assert qapply(RYGate(0, theta)*Qubit('1')) == -sin(theta/2)*Qubit('0') + cos(theta/2)*Qubit('1')

    # Test with concrete angle value pi/2
    assert qapply(RYGate(0, pi/2)*Qubit('0')) == sqrt(2)/2*Qubit('0') + sqrt(2)/2*Qubit('1')
    assert qapply(RYGate(0, pi/2)*Qubit('1')) == -sqrt(2)/2*Qubit('0') + sqrt(2)/2*Qubit('1')


def test_rzgate_qapply():
    """Test RZGate gate with qapply on single-qubit states."""
    from sympy.core.symbol import symbols
    theta = symbols('theta', real=True)

    # Test with symbolic angle
    assert qapply(RZGate(0, theta)*Qubit('0')) == exp(-I*theta/2)*Qubit('0')
    assert qapply(RZGate(0, theta)*Qubit('1')) == exp(I*theta/2)*Qubit('1')

    # Test with concrete angle value pi
    assert qapply(RZGate(0, pi)*Qubit('0')) == -I*Qubit('0')
    assert qapply(RZGate(0, pi)*Qubit('1')) == I*Qubit('1')


def test_compound_gates():
    """Test a compound gate representation."""
    circuit = YGate(0)*ZGate(0)*XGate(0)*HadamardGate(0)*Qubit('00')
    answer = represent(circuit, nqubits=2)
    assert Matrix([I/sqrt(2), I/sqrt(2), 0, 0]) == answer


def test_cnot_gate():
    """Test the CNOT gate."""
    circuit = CNotGate(1, 0)
    assert represent(circuit, nqubits=2) == \
        Matrix([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 0, 1], [0, 0, 1, 0]])
    circuit = circuit*Qubit('111')
    assert matrix_to_qubit(represent(circuit, nqubits=3)) == \
        qapply(circuit)

    circuit = CNotGate(1, 0)
    assert Dagger(circuit) == circuit
    assert Dagger(Dagger(circuit)) == circuit
    assert circuit*circuit == 1


def test_gate_sort():
    """Test gate_sort."""
    for g in (X, Y, Z, H, S, T):
        assert gate_sort(g(2)*g(1)*g(0)) == g(0)*g(1)*g(2)
    e = gate_sort(X(1)*H(0)**2*CNOT(0, 1)*X(1)*X(0))
    assert e == H(0)**2*CNOT(0, 1)*X(0)*X(1)**2
    assert gate_sort(Z(0)*X(0)) == -X(0)*Z(0)
    assert gate_sort(Z(0)*X(0)**2) == X(0)**2*Z(0)
    assert gate_sort(Y(0)*H(0)) == -H(0)*Y(0)
    assert gate_sort(Y(0)*X(0)) == -X(0)*Y(0)
    assert gate_sort(Z(0)*Y(0)) == -Y(0)*Z(0)
    assert gate_sort(T(0)*S(0)) == S(0)*T(0)
    assert gate_sort(Z(0)*S(0)) == S(0)*Z(0)
    assert gate_sort(Z(0)*T(0)) == T(0)*Z(0)
    assert gate_sort(Z(0)*CNOT(0, 1)) == CNOT(0, 1)*Z(0)
    assert gate_sort(S(0)*CNOT(0, 1)) == CNOT(0, 1)*S(0)
    assert gate_sort(T(0)*CNOT(0, 1)) == CNOT(0, 1)*T(0)
    assert gate_sort(X(1)*CNOT(0, 1)) == CNOT(0, 1)*X(1)
    # This takes a long time and should only be uncommented once in a while.
    # nqubits = 5
    # ngates = 10
    # trials = 10
    # for i in range(trials):
    #     c = random_circuit(ngates, nqubits)
    #     assert represent(c, nqubits=nqubits) == \
    #            represent(gate_sort(c), nqubits=nqubits)


def test_gate_simp():
    """Test gate_simp."""
    e = H(0)*X(1)*H(0)**2*CNOT(0, 1)*X(1)**3*X(0)*Z(3)**2*S(4)**3
    assert gate_simp(e) == H(0)*CNOT(0, 1)*S(4)*X(0)*Z(4)
    assert gate_simp(X(0)*X(0)) == 1
    assert gate_simp(Y(0)*Y(0)) == 1
    assert gate_simp(Z(0)*Z(0)) == 1
    assert gate_simp(H(0)*H(0)) == 1
    assert gate_simp(T(0)*T(0)) == S(0)
    assert gate_simp(S(0)*S(0)) == Z(0)
    assert gate_simp(Integer(1)) == Integer(1)
    assert gate_simp(X(0)**2 + Y(0)**2) == Integer(2)


def test_swap_gate():
    """Test the SWAP gate."""
    swap_gate_matrix = Matrix(
        ((1, 0, 0, 0), (0, 0, 1, 0), (0, 1, 0, 0), (0, 0, 0, 1)))
    assert represent(SwapGate(1, 0).decompose(), nqubits=2) == swap_gate_matrix
    assert qapply(SwapGate(1, 3)*Qubit('0010')) == Qubit('1000')
    nqubits = 4
    for i in range(nqubits):
        for j in range(i):
            assert represent(SwapGate(i, j), nqubits=nqubits) == \
                represent(SwapGate(i, j).decompose(), nqubits=nqubits)


def test_one_qubit_commutators():
    """Test single qubit gate commutation relations."""
    # For gates with a parameter (such as an angle) add them here
    # using a partial.
    # Use two different angles to ensure commutators work with different parameters
    alpha, beta = symbols('alpha beta', real=True)

    # Create two gate sets with different angles
    gates_alpha = (
        IdentityGate, X, Y, Z, H, T, S, Sdg, Tdg, SX,
        lambda t: P(t, alpha),
        lambda t: RX(t, alpha),
        lambda t: RY(t, alpha),
        lambda t: RZ(t, alpha)
    )
    gates_beta = (
        IdentityGate, X, Y, Z, H, T, S, Sdg, Tdg, SX,
        lambda t: P(t, beta),
        lambda t: RX(t, beta),
        lambda t: RY(t, beta),
        lambda t: RZ(t, beta)
    )

    for g1 in gates_alpha:
        for g2 in gates_beta:
            e = Commutator(g1(0), g2(0))
            a = matrix_to_zero(represent(e, nqubits=1, format='sympy'))
            b = matrix_to_zero(represent(e.doit(), nqubits=1, format='sympy'))
            # For symbolic parameters, we need to simplify before comparing
            from sympy import simplify
            from sympy.functions.elementary.exponential import exp
            diff = a - b
            # Handle both scalar 0 and zero matrix by checking if all elements simplify to zero
            if diff == 0:
                pass  # Success - both are scalar 0
            elif hasattr(diff, 'shape'):
                # It's a matrix - check each element
                # Use rewrite(exp) then simplify to handle trig/exp identities
                for i in range(diff.shape[0]):
                    for j in range(diff.shape[1]):
                        elem = simplify(diff[i, j].rewrite(exp))
                        assert elem == 0, f"Element [{i},{j}] is {elem}, not 0"
            else:
                assert simplify(diff.rewrite(exp)) == 0

            e = Commutator(g1(0), g2(1))
            assert e.doit() == 0


def test_one_qubit_anticommutators():
    """Test single qubit gate anticommutation relations."""
    for g1 in (IdentityGate, X, Y, Z, H):
        for g2 in (IdentityGate, X, Y, Z, H):
            e = AntiCommutator(g1(0), g2(0))
            a = matrix_to_zero(represent(e, nqubits=1, format='sympy'))
            b = matrix_to_zero(represent(e.doit(), nqubits=1, format='sympy'))
            assert a == b
            e = AntiCommutator(g1(0), g2(1))
            a = matrix_to_zero(represent(e, nqubits=2, format='sympy'))
            b = matrix_to_zero(represent(e.doit(), nqubits=2, format='sympy'))
            assert a == b


def test_cnot_commutators():
    """Test commutators of involving CNOT gates."""
    assert Commutator(CNOT(0, 1), Z(0)).doit() == 0
    assert Commutator(CNOT(0, 1), T(0)).doit() == 0
    assert Commutator(CNOT(0, 1), S(0)).doit() == 0
    assert Commutator(CNOT(0, 1), X(1)).doit() == 0
    assert Commutator(CNOT(0, 1), CNOT(0, 1)).doit() == 0
    assert Commutator(CNOT(0, 1), CNOT(0, 2)).doit() == 0
    assert Commutator(CNOT(0, 2), CNOT(0, 1)).doit() == 0
    assert Commutator(CNOT(1, 2), CNOT(1, 0)).doit() == 0


def test_random_circuit():
    c = random_circuit(10, 3)
    assert isinstance(c, Mul)
    m = represent(c, nqubits=3)
    assert m.shape == (8, 8)
    assert isinstance(m, Matrix)


def test_hermitian_XGate():
    x = XGate(1, 2)
    x_dagger = Dagger(x)

    assert (x == x_dagger)


def test_hermitian_YGate():
    y = YGate(1, 2)
    y_dagger = Dagger(y)

    assert (y == y_dagger)


def test_hermitian_ZGate():
    z = ZGate(1, 2)
    z_dagger = Dagger(z)

    assert (z == z_dagger)


def test_unitary_XGate():
    x = XGate(1, 2)
    x_dagger = Dagger(x)

    assert (x*x_dagger == 1)


def test_unitary_YGate():
    y = YGate(1, 2)
    y_dagger = Dagger(y)

    assert (y*y_dagger == 1)


def test_unitary_ZGate():
    z = ZGate(1, 2)
    z_dagger = Dagger(z)

    assert (z*z_dagger == 1)
