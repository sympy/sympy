from sympy import Symbol
from sympy.physics.vector import Vector
from sympy.physics.mechanics import (dynamicsymbols, Body, PinJoint,
                                     SlidingJoint, CylindricalJoint)


def test_pin_joint():
    parent = Body('parent')
    child = Body('child')
    l = Symbol('l')
    pin_joint = PinJoint('pin_joint', parent, child, child_point_pos=(0, l, 0))

    assert pin_joint.name == 'pin_joint'
    assert pin_joint.parent == parent
    assert pin_joint.child == child
    assert pin_joint.parent_joint_vector == Vector(0)
    assert pin_joint.child_joint_vector == l * child.frame.y
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
    assert child_joint_point.pos_from(parent_joint_point) == Vector(0)
    assert child_joint_point.pos_from(child.masscenter) == l * child.frame.y
    assert parent_joint_point.pos_from(parent.masscenter) == Vector(0)

    assert child.frame.ang_vel_in(parent.frame) == omega * parent.frame.x


def test_sliding_joint():
    parent = Body('parent')
    child = Body('child')
    l = Symbol('l')
    sliding_joint = SlidingJoint('sliding_joint', parent, child,
                                 child_point_pos=(l, 0, 0))

    assert sliding_joint.name == 'sliding_joint'
    assert sliding_joint.parent == parent
    assert sliding_joint.child == child
    assert sliding_joint.parent_joint_vector == Vector(0)
    assert sliding_joint.child_joint_vector == l * child.frame.x
    assert sliding_joint.parent_axis == parent.frame.x
    assert sliding_joint.child_axis == child.frame.x

    dis = dynamicsymbols('sliding_joint_dis')
    disd = dynamicsymbols('sliding_joint_dis', 1)
    vel = dynamicsymbols('sliding_joint_vel')
    assert sliding_joint.coordinates == [dis]
    assert sliding_joint.speeds == [vel]
    assert sliding_joint.kds == [disd - vel]

    child_joint_point = sliding_joint.child_joint_vector
    parent_joint_point = sliding_joint.parent_joint_point
    assert parent_joint_point.name == 'sliding_joint_parent_joint'
    assert child_joint_point.name == 'sliding_joint_child_joint'
    assert child_joint_point.pos_from(parent_joint_point) == dis * parent.frame.x
    assert child_joint_point.pos_from(child.masscenter) == l * child.frame.x
    assert parent_joint_point.pos_from(parent.masscenter) == Vector(0)

    assert parent_joint_point.vel(parent.frame) == Vector(0)
    assert child_joint_point.vel(child.frame) == Vector(0)
    assert child_joint_point.vel(parent.frame) == vel * parent.frame.x
    assert child.masscenter.vel(parent.frame) == vel * parent.frame.x


def test_cylindrical_joint():
    parent = Body('parent')
    child = Body('child')
    l = Symbol('l')
    cylindrical_joint = CylindricalJoint('cylindrical_joint', parent, child,
                                         child_point_pos=(l, l, 0))

    assert cylindrical_joint.name == 'cylindrical_joint'
    assert cylindrical_joint.parent == parent
    assert cylindrical_joint.child == child
    assert cylindrical_joint.parent_joint_vector == Vector(0)
    assert cylindrical_joint.child_joint_vector == l * (child.frame.x + child.frame.y)
    assert cylindrical_joint.parent_axis == parent.frame.x
    assert cylindrical_joint.child_axis == child.frame.x
    
    dis = dynamicsymbols('cylindrical_joint_dis')
    disd = dynamicsymbols('cylindrical_joint_dis', 1)
    vel = dynamicsymbols('cylindrical_joint_vel')

    theta = dynamicsymbols('cylindrical_joint_theta')
    thetad = dynamicsymbols('cylindrical_joint_theta', 1)
    omega = dynamicsymbols('cylindrical_joint_omega')

    assert cylindrical_joint.coordinates == [dis, theta]
    assert cylindrical_joint.speeds == [vel, omega]
    assert cylindrical_joint.kds == [disd - vel, thetad - omega]

    child_joint_point = sliding_joint.child_joint_vector
    parent_joint_point = sliding_joint.parent_joint_point
    assert parent_joint_point.name == 'cylindrical_joint_parent_joint'
    assert child_joint_point.name == 'cylindrical_joint_child_joint'
    assert child_joint_point.pos_from(parent_joint_point) == dis * parent.frame.x
    child_joint_masscenter = l * (child.frame.x + l * child.frame.y)
    assert child_joint_point.pos_from(child.masscenter) == child_joint_masscenter
    assert parent_joint_point.pos_from(parent.masscenter) == Vector(0)
    assert child.frame.ang_vel_in(parent.frame) == omega * parent.frame.x
