from sympy import Symbol, ImmutableMatrix as Matrix
from sympy.physics.mechanics import (Body, PinJoint, SphericalJoint,
                                     JointsMethod)


def test_joint_method_simple():
    parent = Body('parent')
    child = Body('child')
    joints = []
    gravity = Symbol('gravity')
    child.apply_force(gravity * parent.frame.y)
    pin_joint = PinJoint('pinjoint', parent, child)
    joints.append(pin_joint)

    JM = JointsMethod(joints, parent)
    assert JM.rhs() == Matrix([
        [pin_joint.speeds[0]], [0]])
    assert JM.forcelist == child.loads
    assert JM.forcelist == [(child.masscenter, gravity*parent.frame.y)]
    assert JM.bodylist == [child, parent]
    assert JM.q == pin_joint.coordinates
    assert JM.u == pin_joint.speeds
    assert JM.kd == pin_joint.kds
    assert JM.forcing_full == Matrix([
        [pin_joint.speeds[0]], [0]])
    assert JM.mass_matrix_full == Matrix([
        [1,         0],
        [0, Symbol('child_ixx')]])
    assert hasattr(JM, '_qdot')
    assert hasattr(JM, '_udot')
    assert hasattr(JM, '_uaux')
