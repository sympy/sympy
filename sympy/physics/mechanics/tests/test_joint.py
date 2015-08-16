from sympy import Symbol
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
    assert child_joint_point.pos_from(parent_joint_point) == Vector(0)
    assert child_joint_point.pos_from(child.masscenter) == l * child.frame.y
    assert parent_joint_point.pos_from(parent.masscenter) == Vector(0)
    assert child.masscenter.pos_from(parent.masscenter) == - l * child.frame.y
    assert child.frame.ang_vel_in(parent.frame) == omega * parent.frame.x


def test_pin_joint_arguments():
    parent = Body('parent')
    child = Body('child')
    l1 = Symbol('l1')
    l2 = Symbol('l2')
    pin_joint = PinJoint('pin_joint', parent, child,
                         parent_point_pos=(0, l1, 0),
                         child_point_pos=(0, l2, 0),
                         parent_axis=parent.frame.y,
                         child_axis=child.frame.y)
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

    dis = dynamicsymbols('sliding_joint_dis')
    disd = dynamicsymbols('sliding_joint_dis', 1)
    vel = dynamicsymbols('sliding_joint_vel')
    assert sliding_joint.coordinates == [dis]
    assert sliding_joint.speeds == [vel]
    assert sliding_joint.kds == [disd - vel]

    child_joint_point = sliding_joint.child_joint_point
    parent_joint_point = sliding_joint.parent_joint_point
    assert parent_joint_point.name == 'sliding_joint_parent_joint'
    assert child_joint_point.name == 'sliding_joint_child_joint'
    assert child_joint_point.pos_from(parent_joint_point) == dis * parent.frame.x
    assert child_joint_point.pos_from(child.masscenter) == l * child.frame.x
    assert parent_joint_point.pos_from(parent.masscenter) == Vector(0)
    assert child.masscenter.pos_from(parent.masscenter) == (-l * child.frame.x +
                                                            dis * parent.frame.x)

    assert parent_joint_point.vel(parent.frame) == Vector(0)
    assert child_joint_point.vel(child.frame) == Vector(0)
    assert child_joint_point.vel(parent.frame) == vel * parent.frame.x
    assert child.masscenter.vel(parent.frame) == vel * parent.frame.x


def test_sliding_joint_arguments():
    parent = Body('parent')
    child = Body('child')
    l1 = Symbol('l1')
    l2 = Symbol('l2')
    sliding_joint = SlidingJoint('sliding_joint', parent, child,
                                 parent_point_pos=(0, l1, 0),
                                 child_point_pos=(0, l2, 0),
                                 parent_axis=parent.frame.z,
                                 child_axis=child.frame.y)
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

    dis = dynamicsymbols('cylindrical_joint_dis')
    disd = dynamicsymbols('cylindrical_joint_dis', 1)
    vel = dynamicsymbols('cylindrical_joint_vel')

    theta = dynamicsymbols('cylindrical_joint_theta')
    thetad = dynamicsymbols('cylindrical_joint_theta', 1)
    omega = dynamicsymbols('cylindrical_joint_omega')

    assert cylindrical_joint.coordinates == [dis, theta]
    assert cylindrical_joint.speeds == [vel, omega]
    assert cylindrical_joint.kds == [disd - vel, thetad - omega]

    child_joint_point = cylindrical_joint.child_joint_point
    parent_joint_point = cylindrical_joint.parent_joint_point
    assert parent_joint_point.name == 'cylindrical_joint_parent_joint'
    assert child_joint_point.name == 'cylindrical_joint_child_joint'
    # assert child_joint_point.pos_from(parent_joint_point) == dis * parent.frame.x
    child_joint_masscenter = l * (child.frame.x + child.frame.y)
    assert child_joint_point.pos_from(child.masscenter) == child_joint_masscenter
    assert parent_joint_point.pos_from(parent.masscenter) == Vector(0)
    assert child.frame.ang_vel_in(parent.frame) == omega * parent.frame.x


def test_cylindrical_joint_arguments():
    parent = Body('parent')
    child = Body('child')
    l1 = Symbol('l1')
    l2 = Symbol('l2')
    cylindrical_joint = CylindricalJoint('cylindrical_joint', parent, child,
                                         parent_point_pos=(0, l1, 0),
                                         child_point_pos=(0, l2, 0),
                                         parent_axis=parent.frame.y,
                                         child_axis=parent.frame.z)
    assert cylindrical_joint.name == 'cylindrical_joint'
    assert cylindrical_joint.parent == parent
    assert cylindrical_joint.child == child
    assert cylindrical_joint.parent_axis == parent.frame.y
    assert cylindrical_joint.child_axis == parent.frame.z
    assert cylindrical_joint._parent_joint_location == l1 * parent.frame.y
    assert cylindrical_joint._child_joint_location == l2 * child.frame.y


def test_planar_joint():
    parent = Body('parent')
    child = Body('child')
    l = Symbol('l')
    planar_joint = PlanarJoint('planar_joint', parent, child,
                               child_point_pos=(0, l, 0),
                               parent_axis=parent.frame.y,
                               child_axis=child.frame.y)

    assert planar_joint.name == 'planar_joint'
    assert planar_joint.parent == parent
    assert planar_joint.child == child
    assert planar_joint._parent_joint_location == Vector(0)
    assert planar_joint._child_joint_location == l * child.frame.y
    assert planar_joint.parent_axis == parent.frame.y

    theta = dynamicsymbols('planar_joint_theta')  # rotation around z axis.
    thetad = dynamicsymbols('planar_joint_theta', 1)
    omega = dynamicsymbols('planar_joint_omega')
    disx = dynamicsymbols('planar_joint_disx')  # translation along x axis.
    disxd = dynamicsymbols('planar_joint_disx', 1)
    velx = dynamicsymbols('planar_joint_velx')
    disy = dynamicsymbols('planar_joint_disy')  # translation along y axis.
    disyd = dynamicsymbols('planar_joint_disy', 1)
    vely = dynamicsymbols('planar_joint_vely')

    assert planar_joint.coordinates == [theta, disx, disy]
    assert planar_joint.speeds == [omega, velx, vely]
    assert planar_joint.kds == [thetad - omega, disxd - velx, disyd - vely]

    child_joint_point = planar_joint.child_joint_point
    parent_joint_point = planar_joint.parent_joint_point
    assert parent_joint_point.name == 'planar_joint_parent_joint'
    assert child_joint_point.name == 'planar_joint_child_joint'
    assert child_joint_point.pos_from(parent_joint_point) == Vector(0)
    assert child_joint_point.pos_from(child.masscenter) == l * child.frame.y
    assert parent_joint_point.pos_from(parent.masscenter) == Vector(0)

    vecx = cross(planar_joint.parent_axis, planar_joint._child_joint_location)
    assert planar_joint.vecx == vecx
    vecy = cross(planar_joint.parent_axis, vecx)
    assert planar_joint.vecy == vecy
    assert planar_joint.vecz == planar_joint.parent_axis

    assert child.frame.ang_vel_in(parent.frame) == omega * planar_joint.parent_axis
    assert child_joint_point.vel(parent.frame) == disx * vecx + disy * vecy


def test_planar_joint_arguments():
    parent = Body('parent')
    child = Body('child')
    l1 = Symbol('l1')
    l2 = Symbol('l2')
    planar_joint = PlanarJoint('planar_joint', parent, child,
                               parent_point_pos=(0, l1, 0),
                               child_point_pos=(0, l2, 0),
                               parent_axis=parent.frame.y,
                               child_axis=child.frame.y)
    assert planar_joint.name == 'planar_joint'
    assert planar_joint.parent == parent
    assert planar_joint.child == child
    assert planar_joint.parent_axis == parent.frame.y
    assert planar_joint.child_axis == child.frame.y
    assert planar_joint._parent_joint_location == l1 * parent.frame.y
    assert planar_joint._child_joint_location == l2 * child.frame.y


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

    thetax = dynamicsymbols('spherical_joint_thetax')
    thetay = dynamicsymbols('spherical_joint_thetay')
    thetaz = dynamicsymbols('spherical_joint_thetaz')
    thetaxd = dynamicsymbols('spherical_joint_thetax', 1)
    thetayd = dynamicsymbols('spherical_joint_thetay', 1)
    thetazd = dynamicsymbols('spherical_joint_thetaz', 1)
    omegax = dynamicsymbols('spherical_joint_omegax')
    omegay = dynamicsymbols('spherical_joint_omegay')
    omegaz = dynamicsymbols('spherical_joint_omegaz')

    assert spherical_joint.coordinates == [thetax, thetay, thetaz]
    assert spherical_joint.speeds == [omegax, omegay, omegaz]
    assert spherical_joint.kds == [thetaxd - omegax, thetayd - omegay,
                                   thetazd - omegaz]

    child_joint_point = spherical_joint.child_joint_point
    parent_joint_point = spherical_joint.parent_joint_point
    assert parent_joint_point.name == 'spherical_joint_parent_joint'
    assert child_joint_point.name == 'spherical_joint_child_joint'
    assert child_joint_point.pos_from(parent_joint_point) == Vector(0)
    assert child_joint_point.pos_from(child.masscenter) == l * child.frame.y
    assert parent_joint_point.pos_from(parent.masscenter) == Vector(0)

    assert child.frame.ang_vel_in(parent.frame) == omegax * parent.frame.x + \
                                                   omegay * parent.frame.y + \
                                                   omegaz * parent.frame.z


def test_spherical_joint_arguments():
    parent = Body('parent')
    child = Body('child')
    l1 = Symbol('l1')
    l2 = Symbol('l2')
    spherical_joint = SphericalJoint('spherical_joint', parent, child,
                                     parent_point_pos=(0, l1, 0),
                                     child_point_pos=(0, l2, 0))
    assert spherical_joint.name == 'spherical_joint'
    assert spherical_joint.parent == parent
    assert spherical_joint.child == child
    assert spherical_joint._parent_joint_location == l1 * parent.frame.y
    assert spherical_joint._child_joint_location == l2 * child.frame.y
