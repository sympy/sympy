from sympy import Symbol, sin, cos, ImmutableMatrix as Matrix
from sympy.physics.vector import Vector, cross
from sympy.physics.mechanics import (dynamicsymbols, Body, Joint, PinJoint,
                                     SlidingJoint, CylindricalJoint, PlanarJoint,
                                     SphericalJoint)
from sympy.utilities.pytest import raises


def test_joint_abstract_class():
    parent = Body('parent')
    child = Body('child')
    raises(NotImplementedError, lambda: Joint('joint', parent, child))


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


def test_sliding_joint():
    parent = Body('parent')
    child = Body('child')
    l = Symbol('l')
    sliding_joint = SlidingJoint('sliding_joint', parent, child,
                                 child_point_pos=(l, 0, 0))

    assert sliding_joint.name == 'sliding_joint'
    assert sliding_joint.parent == parent
    assert sliding_joint.child == child
    assert sliding_joint._parent_joint_location == Vector(0)
    assert sliding_joint._child_joint_location == l * child.frame.x
    assert sliding_joint.parent_axis == parent.frame.x
    assert sliding_joint.child_axis == child.frame.x

    x = dynamicsymbols('sliding_joint_x')
    xd = dynamicsymbols('sliding_joint_x', 1)
    v = dynamicsymbols('sliding_joint_v')
    assert sliding_joint.coordinates == [x]
    assert sliding_joint.speeds == [v]
    assert sliding_joint.kds == [xd - v]

    child_joint_point = sliding_joint.child_joint_point
    parent_joint_point = sliding_joint.parent_joint_point
    assert parent_joint_point.name == 'sliding_joint_parent_joint'
    assert child_joint_point.name == 'sliding_joint_child_joint'
    assert child_joint_point.pos_from(parent_joint_point) == x * parent.frame.x
    assert child_joint_point.pos_from(child.masscenter) == l * child.frame.x
    assert parent_joint_point.pos_from(parent.masscenter) == Vector(0)
    assert child.masscenter.pos_from(parent.masscenter) == (-l * child.frame.x +
                                                            x * parent.frame.x)

    assert parent_joint_point.vel(parent.frame) == Vector(0)
    assert child_joint_point.vel(child.frame) == Vector(0)
    assert child_joint_point.vel(parent.frame) == v * parent.frame.x
    assert child.masscenter.vel(parent.frame) == v * parent.frame.x
    assert child.frame.ang_vel_in(parent.frame) == 0
    assert cross(parent.frame.x, child.frame.x) == 0
    assert parent.frame.dcm(child.frame) == Matrix([
        [1, 0, 0],
        [0, 1, 0],
        [0, 0, 1]])


def test_sliding_joint_arguments():
    parent = Body('parent')
    child = Body('child')
    l1 = Symbol('l1')
    l2 = Symbol('l2')
    sliding_joint = SlidingJoint('sliding_joint', parent, child,
                                 parent_point_pos=(0, l1, 0),
                                 child_point_pos=(0, l2, 0),
                                 parent_axis='z', child_axis='y')
    assert sliding_joint.name == 'sliding_joint'
    assert sliding_joint.parent == parent
    assert sliding_joint.child == child
    assert sliding_joint._parent_joint_location == l1 * parent.frame.y
    assert sliding_joint._child_joint_location == l2 * child.frame.y
    assert sliding_joint.parent_axis == parent.frame.z
    assert sliding_joint.child_axis == child.frame.y


def test_cylindrical_joint():
    parent = Body('parent')
    child = Body('child')
    l = Symbol('l')
    cylindrical_joint = CylindricalJoint('cylindrical_joint', parent, child,
                                         child_point_pos=(l, l, 0))

    assert cylindrical_joint.name == 'cylindrical_joint'
    assert cylindrical_joint.parent == parent
    assert cylindrical_joint.child == child
    assert cylindrical_joint._parent_joint_location == Vector(0)
    assert cylindrical_joint._child_joint_location == l * (child.frame.x +
                                                        child.frame.y)
    assert cylindrical_joint.parent_axis == parent.frame.x
    assert cylindrical_joint.child_axis == child.frame.x

    x = dynamicsymbols('cylindrical_joint_x')
    xd = dynamicsymbols('cylindrical_joint_x', 1)
    v = dynamicsymbols('cylindrical_joint_v')

    theta = dynamicsymbols('cylindrical_joint_theta')
    thetad = dynamicsymbols('cylindrical_joint_theta', 1)
    omega = dynamicsymbols('cylindrical_joint_omega')

    assert cylindrical_joint.coordinates == [x, theta]
    assert cylindrical_joint.speeds == [v, omega]
    assert cylindrical_joint.kds == [xd - v, thetad - omega]

    child_joint_point = cylindrical_joint.child_joint_point
    parent_joint_point = cylindrical_joint.parent_joint_point
    assert parent_joint_point.name == 'cylindrical_joint_parent_joint'
    assert child_joint_point.name == 'cylindrical_joint_child_joint'
    assert child_joint_point.pos_from(parent_joint_point) == x * parent.frame.x
    child_joint_masscenter = l * (child.frame.x + child.frame.y)
    assert child_joint_point.pos_from(child.masscenter) == child_joint_masscenter
    assert parent_joint_point.pos_from(parent.masscenter) == Vector(0)
    assert child.frame.ang_vel_in(parent.frame) == omega * parent.frame.x
    assert child.masscenter.vel(parent.frame) == - l * omega * child.frame.z
    assert child_joint_point.vel(parent.frame) == v * parent.frame.x
    assert parent.frame.dcm(child.frame) == Matrix([
        [1,          0,           0],
        [0, cos(theta), -sin(theta)],
        [0, sin(theta),  cos(theta)]])


def test_cylindrical_joint_arguments():
    parent = Body('parent')
    child = Body('child')
    l1 = Symbol('l1')
    l2 = Symbol('l2')
    cylindrical_joint = CylindricalJoint('cylindrical_joint', parent, child,
                                         parent_point_pos=(0, l1, 0),
                                         child_point_pos=(0, l2, 0),
                                         parent_axis='y', child_axis='z')
    assert cylindrical_joint.name == 'cylindrical_joint'
    assert cylindrical_joint.parent == parent
    assert cylindrical_joint.child == child
    assert cylindrical_joint.parent_axis == parent.frame.y
    assert cylindrical_joint.child_axis == child.frame.z
    assert cylindrical_joint._parent_joint_location == l1 * parent.frame.y
    assert cylindrical_joint._child_joint_location == l2 * child.frame.y


def test_planar_joint():
    parent = Body('parent')
    child = Body('child')
    l = Symbol('l')
    planar_joint = PlanarJoint('planar_joint', parent, child, child_axis='y',
                               parent_distance=l)

    assert planar_joint.name == 'planar_joint'
    assert planar_joint.parent == parent
    assert planar_joint.child == child
    assert planar_joint._parent_joint_location == l * parent.frame.z
    assert planar_joint._child_joint_location == 0
    assert planar_joint.parent_axis == parent.frame.z
    assert planar_joint.child_axis == child.frame.y

    theta = dynamicsymbols('planar_joint_theta')
    thetad = dynamicsymbols('planar_joint_theta', 1)
    omega = dynamicsymbols('planar_joint_omega')
    x_x = dynamicsymbols('planar_joint_x_x')
    x_xd = dynamicsymbols('planar_joint_x_x', 1)
    v_x = dynamicsymbols('planar_joint_v_x')
    x_y = dynamicsymbols('planar_joint_x_y')
    x_yd = dynamicsymbols('planar_joint_x_y', 1)
    v_y = dynamicsymbols('planar_joint_v_y')

    assert planar_joint.coordinates == [theta, x_x, x_y]
    assert planar_joint.speeds == [omega, v_x, v_y]
    assert planar_joint.kds == [thetad - omega, x_xd - v_x, x_yd - v_y]

    child_joint_point = planar_joint.child_joint_point
    parent_joint_point = planar_joint.parent_joint_point
    assert parent_joint_point.name == 'planar_joint_parent_joint'
    assert child_joint_point.name == 'planar_joint_child_joint'
    assert child_joint_point.pos_from(parent_joint_point) == (x_x*parent.frame.x +
                                                              x_y*parent.frame.y)
    assert child_joint_point.pos_from(child.masscenter) == 0
    assert parent_joint_point.pos_from(parent.masscenter) == l * parent.frame.z
    assert child_joint_point.pos_from(parent.masscenter) == (x_x*parent.frame.x +
                                                             x_y*parent.frame.y +
                                                             l*parent.frame.z)
    assert child.masscenter.pos_from(parent.masscenter) == (x_x*parent.frame.x +
                                                            x_y*parent.frame.y +
                                                            l*parent.frame.z)

    assert parent.masscenter.vel(parent.frame) == 0
    assert child.frame.ang_vel_in(parent.frame) == omega * planar_joint.parent_axis
    assert child_joint_point.vel(parent.frame) == (v_x * parent.frame.x +
                                                   v_y * parent.frame.y)
    assert parent_joint_point.vel(parent.frame) == 0
    assert child.masscenter.vel(parent.frame) == (- omega*x_y*parent.frame.x +
                                                  omega*x_x*parent.frame.y)
    assert parent.frame.dcm(child.frame) == Matrix([
        [cos(theta), -sin(theta), 0],
        [sin(theta),  cos(theta), 0],
        [0, 0, 1]])


def test_planar_joint_arguments():
    parent = Body('parent')
    child = Body('child')
    l1 = Symbol('l1')
    l2 = Symbol('l2')
    planar_joint = PlanarJoint('planar_joint', parent, child, child_axis='y',
                               parent_distance=l1, child_distance=l2)
    assert planar_joint.name == 'planar_joint'
    assert planar_joint.parent == parent
    assert planar_joint.child == child
    assert planar_joint.parent_axis == parent.frame.z
    assert planar_joint.child_axis == child.frame.y
    assert planar_joint._parent_joint_location == l1 * parent.frame.z
    assert planar_joint._child_joint_location == l2 * planar_joint.child_axis


def test_spherical_joint():
    parent = Body('parent')
    child = Body('child')
    l = Symbol('l')
    spherical_joint = SphericalJoint('spherical_joint', parent, child,
                                     child_point_pos=(0, l, 0))

    assert spherical_joint.name == 'spherical_joint'
    assert spherical_joint.parent == parent
    assert spherical_joint.child == child
    assert spherical_joint._parent_joint_location == Vector(0)
    assert spherical_joint._child_joint_location == l * child.frame.y

    alpha = dynamicsymbols('spherical_joint_alpha')
    beta = dynamicsymbols('spherical_joint_beta')
    gamma = dynamicsymbols('spherical_joint_gamma')
    alphad = dynamicsymbols('spherical_joint_alpha', 1)
    betad = dynamicsymbols('spherical_joint_beta', 1)
    gammad = dynamicsymbols('spherical_joint_gamma', 1)
    omega_alpha = dynamicsymbols('spherical_joint_omega_alpha')
    omega_beta = dynamicsymbols('spherical_joint_omega_beta')
    omega_gamma = dynamicsymbols('spherical_joint_omega_gamma')

    assert spherical_joint.coordinates == [alpha, beta, gamma]
    assert spherical_joint.speeds == [omega_alpha, omega_beta, omega_gamma]
    assert spherical_joint.kds == [alphad - omega_alpha, betad - omega_beta,
                                   gammad - omega_gamma]

    child_joint_point = spherical_joint.child_joint_point
    parent_joint_point = spherical_joint.parent_joint_point
    assert parent_joint_point.name == 'spherical_joint_parent_joint'
    assert child_joint_point.name == 'spherical_joint_child_joint'
    assert child_joint_point.pos_from(parent_joint_point) == Vector(0)
    assert child_joint_point.pos_from(child.masscenter) == l * child.frame.y
    assert parent_joint_point.pos_from(parent.masscenter) == Vector(0)
    assert child.frame.ang_vel_in(parent.frame) == (omega_alpha * parent.frame.x +
                                                    omega_beta * parent.frame.y +
                                                    omega_gamma * parent.frame.z)
    assert child_joint_point.vel(child.frame) == 0
    assert child.masscenter.vel(child.frame) == 0
    assert parent_joint_point.vel(parent.frame) == 0
    assert parent.masscenter.vel(parent.frame) == 0
    assert child_joint_point.vel(parent.frame) == 0
    assert child.masscenter.vel(parent.frame) == (l*(omega_alpha*sin(beta) -
        omega_beta*sin(alpha)*cos(beta) +
        omega_gamma*cos(alpha)*cos(beta))*child.frame.x -
        l*((sin(alpha)*sin(gamma) -
        sin(beta)*cos(alpha)*cos(gamma))*omega_gamma +
        (sin(alpha)*sin(beta)*cos(gamma) +
         sin(gamma)*cos(alpha))*omega_beta +
        omega_alpha*cos(beta)*cos(gamma))*child.frame.z)
    assert parent.frame.dcm(child.frame) == Matrix([
        [cos(beta)*cos(gamma), -sin(gamma)*cos(beta), sin(beta)],
        [sin(alpha)*sin(beta)*cos(gamma) + sin(gamma)*cos(alpha),
         -sin(alpha)*sin(beta)*sin(gamma) + cos(alpha)*cos(gamma),
         -sin(alpha)*cos(beta)],
        [sin(alpha)*sin(gamma) - sin(beta)*cos(alpha)*cos(gamma),
         sin(alpha)*cos(gamma) + sin(beta)*sin(gamma)*cos(alpha),
         cos(alpha)*cos(beta)]])


def test_spherical_joint_arguments():
    parent = Body('parent')
    child = Body('child')
    l1 = Symbol('l1')
    l2 = Symbol('l2')
    spherical_joint = SphericalJoint('spherical_joint', parent, child,
                                     parent_point_pos=(0, l1, 0),
                                     child_point_pos=(0, l2, 0),
                                     rot_order='XYX')
    assert spherical_joint.name == 'spherical_joint'
    assert spherical_joint.parent == parent
    assert spherical_joint.child == child
    assert spherical_joint._parent_joint_location == l1 * parent.frame.y
    assert spherical_joint._child_joint_location == l2 * child.frame.y
