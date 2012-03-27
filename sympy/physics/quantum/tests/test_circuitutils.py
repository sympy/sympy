from sympy import Wild, Integer
from sympy.physics.quantum.circuitutils import *
from sympy.physics.quantum.gate import (X, Y, Z, H, S, T, CNOT,
        CGate)

def create_gate_sequence(qubit=0):
    gates = (X(qubit), Y(qubit), Z(qubit), H(qubit))
    return gates

def test_kmp_table():
    word = ('a', 'b', 'c', 'd', 'a', 'b', 'd')
    expected_table = [-1, 0, 0, 0, 0, 1, 2]
    assert expected_table == kmp_table(word)

    word = ('P', 'A', 'R', 'T', 'I', 'C', 'I', 'P', 'A', 'T', 'E', ' ',
            'I', 'N', ' ', 'P', 'A', 'R', 'A', 'C', 'H', 'U', 'T', 'E')
    expected_table = [-1, 0, 0, 0, 0, 0, 0, 0, 1, 2, 0, 0,
                       0, 0, 0, 0, 1, 2, 3, 0, 0, 0, 0, 0]
    assert expected_table == kmp_table(word)

    x = X(0)
    y = Y(0)
    z = Z(0)
    h = H(0)
    word = (x, y, y, x, z)
    expected_table = [-1, 0, 0, 0, 1]
    assert expected_table == kmp_table(word)

    word = (x, x, y, h, z)
    expected_table = [-1, 0, 1, 0, 0]
    assert expected_table == kmp_table(word)

def test_find_subcircuit_with_seq():
    x = X(0)
    y = Y(0)
    z = Z(0)
    h = H(0)
    x1 = X(1)
    y1 = Y(1)

    i0 = Wild('i0')
    x_i0 = X(i0)
    y_i0 = Y(i0)
    z_i0 = Z(i0)

    circuit = (x, y, z)

    assert find_subcircuit_with_seq(circuit, (x,)) == 0
    assert find_subcircuit_with_seq(circuit, (x1,)) == -1
    assert find_subcircuit_with_seq(circuit, (y,)) == 1
    assert find_subcircuit_with_seq(circuit, (h,)) == -1
    assert find_subcircuit_with_seq(circuit, (x, h)) == -1
    assert find_subcircuit_with_seq(circuit, (x, y, z)) == 0
    assert find_subcircuit_with_seq(circuit, (y, z)) == 1
    assert find_subcircuit_with_seq(circuit, (x, y, z, h)) == -1
    assert find_subcircuit_with_seq(circuit, (z, y, x)) == -1
    assert find_subcircuit_with_seq(circuit, (x,), start=2, end=1) == -1

    circuit = (x, y, x, y, z)
    assert find_subcircuit_with_seq(circuit, (x, y, z)) == 2
    assert find_subcircuit_with_seq(circuit, (x,), start=1) == 2
    assert find_subcircuit_with_seq(circuit, (x, y), start=1, end=2) == -1
    assert find_subcircuit_with_seq(circuit, (x, y), start=1, end=3) == -1
    assert find_subcircuit_with_seq(circuit, (x, y), start=1, end=4) == 2
    assert find_subcircuit_with_seq(circuit, (x, y), start=2, end=4) == 2

    circuit = (x, y, z, x1, x, y, z, h, x, y, x1,
               x, y, z, h, y1, h)
    assert find_subcircuit_with_seq(circuit, (x, y, z, h, y1)) == 11

    circuit = (x, y, x_i0, y_i0, z_i0, z)
    assert find_subcircuit_with_seq(circuit, (x_i0, y_i0, z_i0)) == 2

def test_replace_subcircuit_with_seq():
    x = X(0)
    y = Y(0)
    z = Z(0)
    h = H(0)
    cnot = CNOT(1,0)
    cgate_z = CGate((0,), Z(1))

    # Standard cases
    circuit = (z, y, x, x)
    remove = (z, y, x)
    assert replace_subcircuit_with_seq(circuit, remove) == (x,)
    assert replace_subcircuit_with_seq(circuit, remove + (x,)) == ()
    assert replace_subcircuit_with_seq(circuit, remove, pos=1) == circuit
    assert replace_subcircuit_with_seq(circuit, remove, pos=0) == (x,)
    assert replace_subcircuit_with_seq(circuit, (x, x), pos=2) == (z, y)
    assert replace_subcircuit_with_seq(circuit, (h,)) == circuit

    circuit = (x, y, x, y, z)
    remove = (x, y, z)
    assert replace_subcircuit_with_seq(circuit, remove) == (x, y)
    remove = (x, y, x, y)
    assert replace_subcircuit_with_seq(circuit, remove) == (z,)

    circuit = (x, h, cgate_z, h, cnot)
    remove = (x, h, cgate_z)
    assert replace_subcircuit_with_seq(circuit, remove, pos=-1) == (h, cnot)
    assert replace_subcircuit_with_seq(circuit, remove, pos=1) == circuit
    remove = (h, h)
    assert replace_subcircuit_with_seq(circuit, remove) == circuit
    remove = (h, cgate_z, h, cnot)
    assert replace_subcircuit_with_seq(circuit, remove) == (x,)

    replace = (h, x)
    actual = replace_subcircuit_with_seq(circuit, remove,
                     replace=replace)
    assert actual == (x, h, x)

    circuit = (x, y, h, x, y, z)
    remove = (x, y)
    replace = (cnot, cgate_z)
    actual = replace_subcircuit_with_seq(circuit, remove,
                     replace=replace)
    assert actual == (cnot, cgate_z, h, x, y, z)

    actual = replace_subcircuit_with_seq(circuit, remove,
                     replace=replace, pos=1)
    assert actual == (x, y, h, cnot, cgate_z, z)

def test_conv2_symb_indices_with_seq():
    (x, y, z, h) = create_gate_sequence()

    i0 = Wild('i0')
    exp_map = {i0 : Integer(0)}
    actual, act_map, sndx = conv2_symbolic_qubits_with_seq(x)
    assert actual == (X(i0),)
    assert act_map == exp_map

    expected = (X(i0), Y(i0), Z(i0), H(i0))
    exp_map = {i0 : Integer(0)}
    actual, act_map, sndx = conv2_symbolic_qubits_with_seq(x, y, z, h)
    assert actual == expected
    assert exp_map == act_map

    (x1, y1, z1, h1) = create_gate_sequence(1)
    i1 = Wild('i1')

    expected = (X(i0), Y(i0), Z(i0), H(i0))
    exp_map = {i0 : Integer(1)}
    actual, act_map, sndx = conv2_symbolic_qubits_with_seq(x1, y1, z1, h1)
    assert actual == expected
    assert act_map == exp_map

    expected = (X(i0), Y(i0), Z(i0), H(i0), X(i1), Y(i1), Z(i1), H(i1))
    exp_map = {i0 : Integer(0), i1 : Integer(1)}
    actual, act_map, sndx = conv2_symbolic_qubits_with_seq(x, y, z, h,
                                    x1, y1, z1, h1)
    assert actual == expected
    assert act_map == exp_map

    exp_map = {i0 : Integer(1), i1 : Integer(0)}
    actual, act_map, sndx = conv2_symbolic_qubits_with_seq(x1, y1, z1, h1,
                                    x, y, z, h)
    assert actual == expected
    assert act_map == exp_map

    expected = (X(i0), X(i1), Y(i0), Y(i1), Z(i0), Z(i1), H(i0), H(i1))
    exp_map = {i0 : Integer(0), i1 : Integer(1)}
    actual, act_map, sndx = conv2_symbolic_qubits_with_seq(x, x1, y, y1,
                                    z, z1, h, h1)
    assert actual == expected
    assert act_map == exp_map

    exp_map = {i0 : Integer(1), i1 : Integer(0)}
    actual, act_map, sndx = conv2_symbolic_qubits_with_seq(x1, x, y1, y,
                                    z1, z, h1, h)
    assert actual == expected
    assert act_map == exp_map

    cnot_10 = CNOT(1, 0)
    cnot_01 = CNOT(0, 1)
    cgate_z_10 = CGate(1, Z(0))
    cgate_z_01 = CGate(0, Z(1))

    expected = (X(i0), X(i1), Y(i0), Y(i1), Z(i0), Z(i1),
                H(i0), H(i1), CNOT(i1, i0), CNOT(i0, i1),
                CGate(i1, Z(i0)), CGate(i0, Z(i1)))
    exp_map = {i0 : Integer(0), i1 : Integer(1)}
    args = (x, x1, y, y1, z, z1, h, h1, cnot_10, cnot_01,
            cgate_z_10, cgate_z_01)
    actual, act_map, sndx = conv2_symbolic_qubits_with_seq(*args)
    assert actual == expected
    assert act_map == exp_map

    args = (x1, x, y1, y, z1, z, h1, h, cnot_10, cnot_01,
            cgate_z_10, cgate_z_01)
    expected = (X(i0), X(i1), Y(i0), Y(i1), Z(i0), Z(i1),
                H(i0), H(i1), CNOT(i0, i1), CNOT(i1, i0),
                CGate(i0, Z(i1)), CGate(i1, Z(i0)))
    exp_map = {i0 : Integer(1), i1 : Integer(0)}
    actual, act_map, sndx = conv2_symbolic_qubits_with_seq(*args)
    assert actual == expected
    assert act_map == exp_map

    args = (cnot_10, h, cgate_z_01, h)
    expected = (CNOT(i0, i1), H(i1), CGate(i1, Z(i0)), H(i1))
    exp_map = {i0 : Integer(1), i1 : Integer(0)}
    actual, act_map, sndx = conv2_symbolic_qubits_with_seq(*args)
    assert actual == expected
    assert act_map == exp_map

    args = (cnot_01, h1, cgate_z_10, h1)
    exp_map = {i0 : Integer(0), i1 : Integer(1)}
    actual, act_map, sndx = conv2_symbolic_qubits_with_seq(*args)
    assert actual == expected
    assert act_map == exp_map

    args = (cnot_10, h1, cgate_z_01, h1)
    expected = (CNOT(i0, i1), H(i0), CGate(i1, Z(i0)), H(i0))
    exp_map = {i0 : Integer(1), i1 : Integer(0)}
    actual, act_map, sndx = conv2_symbolic_qubits_with_seq(*args)
    assert actual == expected
    assert act_map == exp_map

    i2 = Wild('i2')
    ccgate_z = CGate(0, CGate(1, Z(2)))
    ccgate_x = CGate(1, CGate(2, X(0)))
    args = (ccgate_z, ccgate_x)

    expected = (CGate(i0, CGate(i1, Z(i2))), CGate(i1, CGate(i2, X(i0))))
    exp_map = {i0 : Integer(0), i1 : Integer(1), i2 : Integer(2)}
    actual, act_map, sndx = conv2_symbolic_qubits_with_seq(*args)
    assert actual == expected
    assert act_map == exp_map

    ndx_map = {i0 : Integer(0)}
    actual, act_map, sndx = conv2_symbolic_qubits_with_seq(*args,
                                    qubit_map=ndx_map,
                                    start=i0)
    assert actual == expected
    assert act_map == exp_map

def test_conv2_real_qubits_with_seq():
    i0 = Wild('i0')
    i1 = Wild('i1')

    (x, y, z, h) = create_gate_sequence()

    x_i0 = X(i0)
    y_i0 = Y(i0)
    z_i0 = Z(i0)
    
    qubit_map = {i0 : Integer(0)}
    args = (z_i0, y_i0, x_i0)
    expected = (z, y, x)
    actual = conv2_real_qubits_with_seq(*args, qubit_map=qubit_map)
    assert actual == expected

    cnot_10 = CNOT(1,0)
    cnot_01 = CNOT(0,1)
    cgate_z_10 = CGate(1, Z(0))
    cgate_z_01 = CGate(0, Z(1))

    cnot_i1_i0 = CNOT(i1, i0)
    cnot_i0_i1 = CNOT(i0, i1)
    cgate_z_i1_i0 = CGate(i1, Z(i0))

    qubit_map = {i0 : Integer(0), i1 : Integer(1)}
    args = (cnot_i1_i0,)
    expected = (cnot_10,)
    actual = conv2_real_qubits_with_seq(*args, qubit_map=qubit_map)
    assert actual == expected

    args = (cgate_z_i1_i0,)
    expected = (cgate_z_10,)
    actual = conv2_real_qubits_with_seq(*args, qubit_map=qubit_map)
    assert actual == expected

    args = (cnot_i0_i1,)
    expected = (cnot_01,)
    actual = conv2_real_qubits_with_seq(*args, qubit_map=qubit_map)
    assert actual == expected

    qubit_map = {i0 : Integer(1), i1 : Integer(0)}
    args = (cgate_z_i1_i0,)
    expected = (cgate_z_01,)
    actual = conv2_real_qubits_with_seq(*args, qubit_map=qubit_map)
    assert actual == expected

    i2 = Wild('i2')
    ccgate_z = CGate(i0, CGate(i1, Z(i2)))
    ccgate_x = CGate(i1, CGate(i2, X(i0)))

    qubit_map = {i0 : Integer(0), i1 : Integer(1), i2 : Integer(2)}
    args = (ccgate_z, ccgate_x)
    expected = (CGate(0, CGate(1, Z(2))), CGate(1, CGate(2, X(0))))
    actual = conv2_real_qubits_with_seq(*args, qubit_map=qubit_map)
    assert actual == expected

    qubit_map = {i0 : Integer(1), i2 : Integer(0), i1 : Integer(2)}
    args = (ccgate_x, ccgate_z)
    expected = (CGate(2, CGate(0, X(1))), CGate(1, CGate(2, Z(0))))
    actual = conv2_real_qubits_with_seq(*args, qubit_map=qubit_map)
    assert actual == expected
