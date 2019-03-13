from sympy.external import import_module
from sympy import Mul, Integer
from sympy.physics.quantum.dagger import Dagger
from sympy.physics.quantum.gate import (X, Y, Z, H, CNOT,
        IdentityGate, CGate, PhaseGate, TGate)
from sympy.physics.quantum.identitysearch import (generate_gate_rules,
        generate_equivalent_ids, GateIdentity, bfs_identity_search,
        is_scalar_sparse_matrix,
        is_scalar_nonsparse_matrix, is_degenerate, is_reducible)
from sympy.utilities.pytest import skip, XFAIL


def create_gate_sequence(qubit=0):
    gates = (X(qubit), Y(qubit), Z(qubit), H(qubit))
    return gates


def test_generate_gate_rules_1():
    # Test with tuples
    (x, y, z, h) = create_gate_sequence()
    ph = PhaseGate(0)
    cgate_t = CGate(0, TGate(1))

    assert generate_gate_rules((x,)) == {((x,), ())}

    gate_rules = set([((x, x), ()),
                      ((x,), (x,))])
    assert generate_gate_rules((x, x)) == gate_rules

    gate_rules = set([((x, y, x), ()),
                      ((y, x, x), ()),
                      ((x, x, y), ()),
                      ((y, x), (x,)),
                      ((x, y), (x,)),
                      ((y,), (x, x))])
    assert generate_gate_rules((x, y, x)) == gate_rules

    gate_rules = set([((x, y, z), ()), ((y, z, x), ()), ((z, x, y), ()),
                      ((), (x, z, y)), ((), (y, x, z)), ((), (z, y, x)),
                      ((x,), (z, y)), ((y, z), (x,)), ((y,), (x, z)),
                      ((z, x), (y,)), ((z,), (y, x)), ((x, y), (z,))])
    actual = generate_gate_rules((x, y, z))
    assert actual == gate_rules

    gate_rules = set(
        [((), (h, z, y, x)), ((), (x, h, z, y)), ((), (y, x, h, z)),
         ((), (z, y, x, h)), ((h,), (z, y, x)), ((x,), (h, z, y)),
         ((y,), (x, h, z)), ((z,), (y, x, h)), ((h, x), (z, y)),
         ((x, y), (h, z)), ((y, z), (x, h)), ((z, h), (y, x)),
         ((h, x, y), (z,)), ((x, y, z), (h,)), ((y, z, h), (x,)),
         ((z, h, x), (y,)), ((h, x, y, z), ()), ((x, y, z, h), ()),
         ((y, z, h, x), ()), ((z, h, x, y), ())])
    actual = generate_gate_rules((x, y, z, h))
    assert actual == gate_rules

    gate_rules = set([((), (cgate_t**(-1), ph**(-1), x)),
                      ((), (ph**(-1), x, cgate_t**(-1))),
                      ((), (x, cgate_t**(-1), ph**(-1))),
                      ((cgate_t,), (ph**(-1), x)),
                      ((ph,), (x, cgate_t**(-1))),
                      ((x,), (cgate_t**(-1), ph**(-1))),
                      ((cgate_t, x), (ph**(-1),)),
                      ((ph, cgate_t), (x,)),
                      ((x, ph), (cgate_t**(-1),)),
                      ((cgate_t, x, ph), ()),
                      ((ph, cgate_t, x), ()),
                      ((x, ph, cgate_t), ())])
    actual = generate_gate_rules((x, ph, cgate_t))
    assert actual == gate_rules

    gate_rules = set([(Integer(1), cgate_t**(-1)*ph**(-1)*x),
                      (Integer(1), ph**(-1)*x*cgate_t**(-1)),
                      (Integer(1), x*cgate_t**(-1)*ph**(-1)),
                      (cgate_t, ph**(-1)*x),
                      (ph, x*cgate_t**(-1)),
                      (x, cgate_t**(-1)*ph**(-1)),
                      (cgate_t*x, ph**(-1)),
                      (ph*cgate_t, x),
                      (x*ph, cgate_t**(-1)),
                      (cgate_t*x*ph, Integer(1)),
                      (ph*cgate_t*x, Integer(1)),
                      (x*ph*cgate_t, Integer(1))])
    actual = generate_gate_rules((x, ph, cgate_t), return_as_muls=True)
    assert actual == gate_rules


def test_generate_gate_rules_2():
    # Test with Muls
    (x, y, z, h) = create_gate_sequence()
    ph = PhaseGate(0)
    cgate_t = CGate(0, TGate(1))

    # Note: 1 (type int) is not the same as 1 (type One)
    expected = {(x, Integer(1))}
    assert generate_gate_rules((x,), return_as_muls=True) == expected

    expected = {(Integer(1), Integer(1))}
    assert generate_gate_rules(x*x, return_as_muls=True) == expected

    expected = {((), ())}
    assert generate_gate_rules(x*x, return_as_muls=False) == expected

    gate_rules = set([(x*y*x, Integer(1)),
                      (y, Integer(1)),
                      (y*x, x),
                      (x*y, x)])
    assert generate_gate_rules(x*y*x, return_as_muls=True) == gate_rules

    gate_rules = set([(x*y*z, Integer(1)),
                      (y*z*x, Integer(1)),
                      (z*x*y, Integer(1)),
                      (Integer(1), x*z*y),
                      (Integer(1), y*x*z),
                      (Integer(1), z*y*x),
                      (x, z*y),
                      (y*z, x),
                      (y, x*z),
                      (z*x, y),
                      (z, y*x),
                      (x*y, z)])
    actual = generate_gate_rules(x*y*z, return_as_muls=True)
    assert actual == gate_rules

    gate_rules = set([(Integer(1), h*z*y*x),
                      (Integer(1), x*h*z*y),
                      (Integer(1), y*x*h*z),
                      (Integer(1), z*y*x*h),
                      (h, z*y*x), (x, h*z*y),
                      (y, x*h*z), (z, y*x*h),
                      (h*x, z*y), (z*h, y*x),
                      (x*y, h*z), (y*z, x*h),
                      (h*x*y, z), (x*y*z, h),
                      (y*z*h, x), (z*h*x, y),
                      (h*x*y*z, Integer(1)),
                      (x*y*z*h, Integer(1)),
                      (y*z*h*x, Integer(1)),
                      (z*h*x*y, Integer(1))])
    actual = generate_gate_rules(x*y*z*h, return_as_muls=True)
    assert actual == gate_rules

    gate_rules = set([(Integer(1), cgate_t**(-1)*ph**(-1)*x),
                      (Integer(1), ph**(-1)*x*cgate_t**(-1)),
                      (Integer(1), x*cgate_t**(-1)*ph**(-1)),
                      (cgate_t, ph**(-1)*x),
                      (ph, x*cgate_t**(-1)),
                      (x, cgate_t**(-1)*ph**(-1)),
                      (cgate_t*x, ph**(-1)),
                      (ph*cgate_t, x),
                      (x*ph, cgate_t**(-1)),
                      (cgate_t*x*ph, Integer(1)),
                      (ph*cgate_t*x, Integer(1)),
                      (x*ph*cgate_t, Integer(1))])
    actual = generate_gate_rules(x*ph*cgate_t, return_as_muls=True)
    assert actual == gate_rules

    gate_rules = set([((), (cgate_t**(-1), ph**(-1), x)),
                      ((), (ph**(-1), x, cgate_t**(-1))),
                      ((), (x, cgate_t**(-1), ph**(-1))),
                      ((cgate_t,), (ph**(-1), x)),
                      ((ph,), (x, cgate_t**(-1))),
                      ((x,), (cgate_t**(-1), ph**(-1))),
                      ((cgate_t, x), (ph**(-1),)),
                      ((ph, cgate_t), (x,)),
                      ((x, ph), (cgate_t**(-1),)),
                      ((cgate_t, x, ph), ()),
                      ((ph, cgate_t, x), ()),
                      ((x, ph, cgate_t), ())])
    actual = generate_gate_rules(x*ph*cgate_t)
    assert actual == gate_rules


def test_generate_equivalent_ids_1():
    # Test with tuples
    (x, y, z, h) = create_gate_sequence()

    assert generate_equivalent_ids((x,)) == {(x,)}
    assert generate_equivalent_ids((x, x)) == {(x, x)}
    assert generate_equivalent_ids((x, y)) == {(x, y), (y, x)}

    gate_seq = (x, y, z)
    gate_ids = set([(x, y, z), (y, z, x), (z, x, y), (z, y, x),
                    (y, x, z), (x, z, y)])
    assert generate_equivalent_ids(gate_seq) == gate_ids

    gate_ids = set([Mul(x, y, z), Mul(y, z, x), Mul(z, x, y),
                    Mul(z, y, x), Mul(y, x, z), Mul(x, z, y)])
    assert generate_equivalent_ids(gate_seq, return_as_muls=True) == gate_ids

    gate_seq = (x, y, z, h)
    gate_ids = set([(x, y, z, h), (y, z, h, x),
                    (h, x, y, z), (h, z, y, x),
                    (z, y, x, h), (y, x, h, z),
                    (z, h, x, y), (x, h, z, y)])
    assert generate_equivalent_ids(gate_seq) == gate_ids

    gate_seq = (x, y, x, y)
    gate_ids = {(x, y, x, y), (y, x, y, x)}
    assert generate_equivalent_ids(gate_seq) == gate_ids

    cgate_y = CGate((1,), y)
    gate_seq = (y, cgate_y, y, cgate_y)
    gate_ids = {(y, cgate_y, y, cgate_y), (cgate_y, y, cgate_y, y)}
    assert generate_equivalent_ids(gate_seq) == gate_ids

    cnot = CNOT(1, 0)
    cgate_z = CGate((0,), Z(1))
    gate_seq = (cnot, h, cgate_z, h)
    gate_ids = set([(cnot, h, cgate_z, h), (h, cgate_z, h, cnot),
                    (h, cnot, h, cgate_z), (cgate_z, h, cnot, h)])
    assert generate_equivalent_ids(gate_seq) == gate_ids


def test_generate_equivalent_ids_2():
    # Test with Muls
    (x, y, z, h) = create_gate_sequence()

    assert generate_equivalent_ids((x,), return_as_muls=True) == {x}

    gate_ids = {Integer(1)}
    assert generate_equivalent_ids(x*x, return_as_muls=True) == gate_ids

    gate_ids = {x*y, y*x}
    assert generate_equivalent_ids(x*y, return_as_muls=True) == gate_ids

    gate_ids = {(x, y), (y, x)}
    assert generate_equivalent_ids(x*y) == gate_ids

    circuit = Mul(*(x, y, z))
    gate_ids = set([x*y*z, y*z*x, z*x*y, z*y*x,
                    y*x*z, x*z*y])
    assert generate_equivalent_ids(circuit, return_as_muls=True) == gate_ids

    circuit = Mul(*(x, y, z, h))
    gate_ids = set([x*y*z*h, y*z*h*x,
                    h*x*y*z, h*z*y*x,
                    z*y*x*h, y*x*h*z,
                    z*h*x*y, x*h*z*y])
    assert generate_equivalent_ids(circuit, return_as_muls=True) == gate_ids

    circuit = Mul(*(x, y, x, y))
    gate_ids = {x*y*x*y, y*x*y*x}
    assert generate_equivalent_ids(circuit, return_as_muls=True) == gate_ids

    cgate_y = CGate((1,), y)
    circuit = Mul(*(y, cgate_y, y, cgate_y))
    gate_ids = {y*cgate_y*y*cgate_y, cgate_y*y*cgate_y*y}
    assert generate_equivalent_ids(circuit, return_as_muls=True) == gate_ids

    cnot = CNOT(1, 0)
    cgate_z = CGate((0,), Z(1))
    circuit = Mul(*(cnot, h, cgate_z, h))
    gate_ids = set([cnot*h*cgate_z*h, h*cgate_z*h*cnot,
                    h*cnot*h*cgate_z, cgate_z*h*cnot*h])
    assert generate_equivalent_ids(circuit, return_as_muls=True) == gate_ids


@XFAIL
def test_bfs_identity_search_xfail():
    s = PhaseGate(0)
    t = TGate(0)
    gate_list = [Dagger(s), t]
    id_set = {GateIdentity(Dagger(s), t, t)}
    assert bfs_identity_search(gate_list, 1, max_depth=3) == id_set
