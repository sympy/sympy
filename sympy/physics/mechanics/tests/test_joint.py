from sympy import Symbol, sin, cos, ImmutableMatrix as Matrix
from sympy.physics.vector import Vector, cross
from sympy.physics.mechanics import dynamicsymbols, Body, PinJoint
from sympy.physics.mechanics.joint import Joint
from sympy.utilities.pytest import raises


def test_Joint():
    parent = Body('parent')
    child = Body('child')
    raises(TypeError, lambda: Joint('J', parent, child))


def test_pin_joint():
    parent = Body('parent')
    child = Body('child')
    l = Symbol('l')
    pin_joint = PinJoint('pin_joint', parent, child, child_point_pos=(0, l, 0))

    assert pin_joint.name == 'pin_joint'
    assert pin_joint.parent == parent
    assert pin_joint.child == child
    assert pin_joint._parent_joint_location == Vector(0)
    assert pin_joint._child_joint_location == l * child.frame.y
    assert pin_joint.parent_axis == parent.frame.x
    assert pin_joint.child_axis == child.frame.x

    theta = dynamicsymbols('pin_joint_theta')
    thetad = dynamicsymbols('pin_joint_theta', 1)
    omega = dynamicsymbols('pin_joint_omega')
    assert pin_joint.coordinates == [theta]
    assert pin_joint.speeds == [omega]
    assert pin_joint.kds == [thetad - omega]

    child_joint_point = pin_joint.child_joint_point
    parent_joint_point = pin_joint.parent_joint_point
    assert parent_joint_point.name == 'pin_joint_parent_joint'
    assert child_joint_point.name == 'pin_joint_child_joint'
    assert child_joint_point.pos_from(parent_joint_point) == 0
    assert child_joint_point.pos_from(child.masscenter) == l * child.frame.y
    assert parent_joint_point.pos_from(parent.masscenter) == 0
    assert child.masscenter.pos_from(parent.masscenter) == - l * child.frame.y
    assert child.frame.ang_vel_in(parent.frame) == omega * parent.frame.x

    assert child.masscenter.vel(parent.frame) == - l * omega * child.frame.z
    assert parent.masscenter.vel(parent.frame) == 0
    assert pin_joint.parent_joint_point.vel(parent.frame) == 0
    assert pin_joint.child_joint_point.vel(parent.frame) == 0
    assert parent.frame.dcm(child.frame) == Matrix([
        [1,          0,           0],
        [0, cos(theta), -sin(theta)],
        [0, sin(theta),  cos(theta)]])


def test_pin_joint_arguments():
    parent = Body('parent')
    child = Body('child')
    l1 = Symbol('l1')
    l2 = Symbol('l2')
    pin_joint = PinJoint('pin_joint', parent, child,
                         parent_point_pos=(0, l1, 0),
                         child_point_pos=(0, l2, 0),
                         parent_axis='y',
                         child_axis='y')
    assert pin_joint.name == 'pin_joint'
    assert pin_joint.parent == parent
    assert pin_joint.child == child
    assert pin_joint.parent_axis == parent.frame.y
    assert pin_joint.child_axis == child.frame.y
    assert pin_joint._parent_joint_location == l1 * parent.frame.y
    assert pin_joint._child_joint_location == l2 * child.frame.y
