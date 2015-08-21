from sympy import Symbol, symbols, ImmutableMatrix as Matrix
from sympy.physics.mechanics import *


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


def test_joint_method_simple_pendulum():
    # simple pendulum
    theta = dynamicsymbols('theta')
    thetad = dynamicsymbols('theta', 1)
    omega = dynamicsymbols('omega')
    length, gravity = symbols('length gravity')
    child_mass, parent_mass = symbols('child_mass parent_mass')

    parent_frame = ReferenceFrame('parent_frame')
    child_frame = parent_frame.orientnew('child_frame', 'Axis',
                                         [theta, parent_frame.z])

    child_frame.set_ang_vel(parent_frame, omega * parent_frame.z)

    parent_masscenter= Point('parent_masscenter')
    child_masscenter = parent_masscenter.locatenew('child_masscenter',
                                                   length * child_frame.x)

    parent_masscenter.set_vel(parent_frame, 0)
    child_masscenter.v2pt_theory(parent_masscenter, parent_frame, child_frame)

    parent_ixx, parent_iyy, parent_izz = symbols('parent_ixx parent_iyy parent_izz')
    parent_izx, parent_ixy, parent_iyz = symbols('parent_izx parent_ixy parent__iyz')
    parent_inertia = (inertia(parent_frame, parent_ixx, parent_iyy, parent_izz,
                      parent_ixy, parent_iyz, parent_izx), parent_masscenter)

    child_ixx, child_iyy, child_izz = symbols('child_ixx child_iyy child_izz')
    child_izx, child_ixy, child_iyz = symbols('child_izx child_ixy child__iyz')
    child_inertia = (inertia(child_frame, child_ixx, child_iyy, child_izz,
                     child_ixy, child_iyz, child_izx), child_masscenter)

    parent = RigidBody('parent', parent_masscenter, parent_frame, parent_mass,
                       parent_inertia)
    child = RigidBody('child', parent_masscenter, child_frame, parent_mass,
                      child_inertia)

    kd = [thetad - omega]
    FL = [(child_masscenter, child_mass * gravity * parent_frame.x)]
    BL = [parent, child]

    KM = KanesMethod(parent_frame, q_ind=[theta], u_ind=[omega], kd_eqs=kd)
    KM.kanes_equations(FL, BL)

    # Simple pendulum using joints method.
    parent = Body('parent')
    child = Body('child')
    length = Symbol('length')
    gravity = Symbol('gravity')
    child.apply_force(child.mass * gravity * child.frame.x)
    pin_joint = PinJoint('pin_joint', parent, child,
                         child_point_pos=(0, length, 0),
                         parent_axis='Z', child_axis='Z')
    JM = JointsMethod([pin_joint], parent)

    assert KM.rhs() == JM.rhs()
