from sympy import acos, cos, eye, Matrix, pprint, simplify, sin, sqrt, symbols, zeros
from sympy.physics.mechanics import (dynamicsymbols, Body, Joint, PinJoint,
                                     SlidingJoint)
from sympy.utilities.pytest import raises


def test_no_apply_joint_method():
    x = dynamicsymbols('x')
    y = x.diff()
    parent = Body('parent')
    child = Body('child')
    with raises(NotImplementedError):
        Joint('joint', parent, child, [x], [y])


def test_no_spatial_info():
    class JointChild(Joint):
        def apply_joint(self):
            # Some basic definitions so that the joint will be created
            # correctly. This code resembles PinJoint's apply_joint method
            tempname = self.name + "_child_joint_frame"
            tempaxis = self.parent_joint_frame.x
            temp = self.parent_joint_frame.orientnew(tempname, "Axis",
                                                     [0, tempaxis])
            self.child_joint_frame = temp
            self.child_joint_frame.set_ang_vel(self.parent_joint_frame,
                                               self.speeds[0] *
                                               self.parent_joint_frame.z)
            self.child_joint_point.set_pos(self.parent_joint_point, 0)
            self.child_joint_point.set_vel(self.parent_joint_frame, 0)

    x = dynamicsymbols('x')
    y = x.diff()
    parent = Body('parent')
    child = Body('child')
    joint = JointChild('joint', parent, child, [x], [y])
    with raises(NotImplementedError):
        joint.spatial_info()


def test_pinjoint_attributes():
    """Test that the attributes of PinJoint match the correct defaults"""
    x = dynamicsymbols('x')
    y = x.diff()
    parent = Body('parent')
    child = Body('child')
    joint = PinJoint('Joint', parent, child, x, y)

    # Set up the expected return results
    name = 'Joint'
    child_joints = [joint]
    parent_joint = joint
    coordinates = Matrix([x])
    speeds = Matrix([y])
    kin_diff = [x - y]

    # Test that the attributes return the above expected results
    assert joint.name == name
    assert joint.parent == parent
    assert joint.child == child
    assert joint.parent.child_joints == child_joints
    assert joint.child.parent_joint == parent_joint
    assert joint.coordinates == coordinates
    assert joint.speeds == speeds
    assert joint.kin_diff == kin_diff


def test_pinjoint_spatial():
    """Test the spatial components of a default PinJoint"""
    x = dynamicsymbols('x')
    y = x.diff()
    parent = Body('parent')
    child = Body('child')
    joint = PinJoint('Joint', parent, child, x, y)

    # Set up the expected return results
    joint_transform = Matrix([[cos(x),  sin(x), 0,    0,      0,    0],
                              [-sin(x), cos(x), 0,    0,      0,    0],
                              [0,         0,    1,    0,      0,    0],
                              [0,         0,    0,  cos(x), sin(x), 0],
                              [0,         0,    0, -sin(x), cos(x), 0],
                              [0,         0,    0,    0,      0,    1]])
    motion_subspace = Matrix([0, 0, 1, 0, 0, 0])
    spat_vel = Matrix([0, 0, y, 0, 0, 0])

    # Test the returned spatial results against the expected results
    assert joint.joint_transform == joint_transform
    assert joint.motion_subspace == motion_subspace
    out1, out2, out3 = joint.spatial_info()
    assert [out1, out2, out3] == [joint_transform, motion_subspace, spat_vel]


def test_pinjoint_defaults_points():
    """Test the relations between the 4 points in a default PinJoint"""
    x = dynamicsymbols('x')
    y = x.diff()
    parent = Body('parent')
    child = Body('child')
    joint = PinJoint('Joint', parent, child, x, y)

    # Collect the points for testing
    pmasscenter = joint.parent.masscenter
    pjointpoint = joint.parent_joint_point
    cjointpoint = joint.child_joint_point
    cmasscenter = joint.child.masscenter

    # Collect the frames for testing
    pframe = joint.parent.frame
    pjointframe = joint.parent_joint_frame
    cjointframe = joint.child_joint_frame
    cframe = joint.child.frame

    # Tests for the parent.masscenter
    assert pmasscenter.pos_from(pjointpoint) == 0
    assert pmasscenter.pos_from(cjointpoint) == 0
    assert pmasscenter.pos_from(cmasscenter) == 0

    assert pmasscenter.vel(pframe) == 0
    assert pmasscenter.vel(pjointframe) == 0
    assert pmasscenter.vel(cjointframe) == 0
    assert pmasscenter.vel(cframe) == 0

    # Tests for the parent_joint_point
    assert pjointpoint.name == 'Joint_parent_joint_point'

    assert pjointpoint.pos_from(cjointpoint) == 0
    assert pjointpoint.pos_from(cmasscenter) == 0

    assert pjointpoint.vel(pframe) == 0
    assert pjointpoint.vel(pjointframe) == 0
    assert pjointpoint.vel(cjointframe) == 0
    assert pjointpoint.vel(cframe) == 0

    # Tests for the child_joint_point
    assert cjointpoint.name == 'Joint_child_joint_point'

    assert cjointpoint.pos_from(cmasscenter) == 0

    assert cjointpoint.vel(pframe) == 0
    assert cjointpoint.vel(pjointframe) == 0
    assert cjointpoint.vel(cjointframe) == 0
    assert cjointpoint.vel(cframe) == 0

    # Tests for the child.mass
    assert cmasscenter.vel(pframe) == 0
    assert cmasscenter.vel(pjointframe) == 0
    assert cmasscenter.vel(cjointframe) == 0
    assert cmasscenter.vel(cframe) == 0


def test_pinjoint_defaults_frames():
    """Test the relations between the 4 frames in a default PinJoint"""
    x = dynamicsymbols('x')
    y = x.diff()
    parent = Body('parent')
    child = Body('child')
    joint = PinJoint('Joint', parent, child, x, y)

    # Collect the frames for testing
    pframe = joint.parent.frame
    pjointframe = joint.parent_joint_frame
    cjointframe = joint.child_joint_frame
    cframe = joint.child.frame

    # Tests for the parent.frame
    pframe2pjointframe = eye(3)

    pframe2cjointframe = Matrix([[cos(x), -sin(x), 0],
                                 [sin(x),  cos(x), 0],
                                 [0,            0, 1]])

    pframe2cframe = Matrix([[cos(x), -sin(x), 0],
                            [sin(x),  cos(x), 0],
                            [0,            0, 1]])

    assert pframe.dcm(pjointframe) == pframe2pjointframe
    assert pframe.dcm(cjointframe) == pframe2cjointframe
    assert pframe.dcm(cframe) == pframe2cframe

    assert pframe.ang_vel_in(pjointframe) == 0
    assert pframe.ang_vel_in(cjointframe) == -y * pjointframe.z
    assert pframe.ang_vel_in(cframe) == -y * pjointframe.z

    # Tests for the parent_joint_frame
    assert pjointframe.name == 'Joint_parent_joint_frame'

    pjointframe2cjointframe = Matrix([[cos(x), -sin(x), 0],
                                      [sin(x),  cos(x), 0],
                                      [0,            0, 1]])

    pjointframe2cframe = Matrix([[cos(x), -sin(x), 0],
                                 [sin(x),  cos(x), 0],
                                 [0,            0, 1]])

    assert pjointframe.dcm(cjointframe) == pjointframe2cjointframe
    assert pjointframe.dcm(cframe) == pjointframe2cframe

    assert pjointframe.ang_vel_in(pframe) == 0
    assert pjointframe.ang_vel_in(cjointframe) == -y * pjointframe.z
    assert pjointframe.ang_vel_in(cframe) == -y * pjointframe.z

    # Tests for the child_joint_frame
    assert cjointframe.name == 'Joint_child_joint_frame'

    cjointframe2cframe = eye(3)

    assert cjointframe.dcm(cframe) == cjointframe2cframe

    assert cjointframe.ang_vel_in(pframe) == y * pjointframe.z
    assert cjointframe.ang_vel_in(pjointframe) == y * pjointframe.z
    assert cjointframe.ang_vel_in(cframe) == 0

    # Tests for the child.frame
    assert cframe.ang_vel_in(pframe) == y * pjointframe.z
    assert cframe.ang_vel_in(pjointframe) == y * pjointframe.z
    assert cframe.ang_vel_in(cjointframe) == 0


def test_pinjoint_defaults_XTchild():
    """Test the output of XT_child for a default PinJoint"""
    x, y = dynamicsymbols('x y')
    u = x.diff()
    v = y.diff()
    parent = Body('parent')
    child = Body('child')
    child2 = Body('child2')
    joint1 = PinJoint('Joint1', parent, child, x, u)
    joint2 = PinJoint('Joint2', child, child2, y, v)

    # Test the XT_child of the first joint
    XTchild_exp1 = eye(6)
    assert joint1.XT_child() == XTchild_exp1

    # Test the XT_child of the first joint
    XTchild_exp2 = eye(6)
    assert joint2.XT_child() == XTchild_exp2


def test_pinjoint_fullargs_points():
    """Test the relations between the 4 points in a PinJoint for which all args
    were specified"""
    qx, qy, qz = symbols('qx qy qz')
    parent_point_pos = (qx, qy, qz)

    rx, ry, rz = symbols('rx ry rz')
    child_point_pos = (rx, ry, rz)

    sx, sy, sz = symbols('sx sy sz')
    parent_axis = (sx, sy, sz)

    tx, ty, tz = symbols('tx ty tz')
    child_axis = (tx, ty, tz)

    theta = dynamicsymbols('theta')
    omega = theta.diff()

    parent = Body('parent')
    child = Body('child')

    joint = PinJoint('Joint', parent, child, theta, omega, parent_point_pos,
                     child_point_pos, parent_axis, child_axis)

    # Collect the points for testing
    pmasscenter = joint.parent.masscenter
    pjointpoint = joint.parent_joint_point
    cjointpoint = joint.child_joint_point
    cmasscenter = joint.child.masscenter

    # Collect the frames for testing
    pframe = joint.parent.frame
    pjointframe = joint.parent_joint_frame
    cjointframe = joint.child_joint_frame
    cframe = joint.child.frame

    # Create the position vectors for testing
    q = qx*pframe.x + qy*pframe.y + qz*pframe.z
    r = rx*cframe.x + ry*cframe.y + rz*cframe.z

    # Tests for the parent.masscenter
    pmasscenter2pjointpoint = q
    pmasscenter2cjointpoint = q
    pmasscenter2cmasscenter = q - r

    assert pjointpoint.pos_from(pmasscenter) == pmasscenter2pjointpoint
    assert cjointpoint.pos_from(pmasscenter) == pmasscenter2cjointpoint
    assert cmasscenter.pos_from(pmasscenter) == pmasscenter2cmasscenter

    # assert pmasscenter.vel(pframe) == 0
    # assert pmasscenter.vel(pjointframe) == 0
    # assert pmasscenter.vel(cjointframe) == 0
    # assert pmasscenter.vel(cframe) == 0

    # Tests for the parent_joint_point
    assert pjointpoint.name == 'Joint_parent_joint_point'

    assert pjointpoint.pos_from(cjointpoint) == 0
    assert pjointpoint.pos_from(cmasscenter) == -r

    # assert pjointpoint.vel(pframe) == 0
    # assert pjointpoint.vel(pjointframe) == 0
    # assert pjointpoint.vel(cjointframe) == 0
    # assert pjointpoint.vel(cframe) == 0

    # Tests for the child_joint_point
    assert cjointpoint.name == 'Joint_child_joint_point'

    assert cjointpoint.pos_from(cmasscenter) == -r

    # assert cjointpoint.vel(pframe) == 0
    # assert cjointpoint.vel(pjointframe) == 0
    # assert cjointpoint.vel(cjointframe) == 0
    # assert cjointpoint.vel(cframe) == 0

    # Tests for the child.mass
    # assert cmasscenter.vel(pframe) == 0
    # assert cmasscenter.vel(pjointframe) == 0
    # assert cmasscenter.vel(cjointframe) == 0
    # assert cmasscenter.vel(cframe) == 0


def test_pinjoint_fullargs_frames():
    """Test the relations between the 4 frames in a PinJoint for which all args
    were specified"""
    qx, qy, qz = symbols('qx qy qz')
    parent_point_pos = (qx, qy, qz)

    rx, ry, rz = symbols('rx ry rz')
    child_point_pos = (rx, ry, rz)

    sx, sy, sz = symbols('sx sy sz')
    parent_axis = (sx, sy, sz)

    tx, ty, tz = symbols('tx ty tz')
    child_axis = (tx, ty, tz)

    theta = dynamicsymbols('theta')
    omega = theta.diff()

    parent = Body('parent')
    child = Body('child')

    joint = PinJoint('Joint', parent, child, theta, omega, parent_point_pos,
                     child_point_pos, parent_axis, child_axis)

    # Collect the points for testing
    pmasscenter = joint.parent.masscenter
    pjointpoint = joint.parent_joint_point
    cjointpoint = joint.child_joint_point
    cmasscenter = joint.child.masscenter

    # Collect the frames for testing
    pframe = joint.parent.frame
    pjointframe = joint.parent_joint_frame
    cjointframe = joint.child_joint_frame
    cframe = joint.child.frame

    # Create the position vectors for testing
    q = qx*pframe.x + qy*pframe.y + qz*pframe.z
    r = rx*cframe.x + ry*cframe.y + rz*cframe.z

    # Create the expected rotation matrices
    phi_P = acos(sz / sqrt(sx**2 + sy**2 + sz**2))
    phi_C = acos(tz / sqrt(tx**2 + ty**2 + tz**2))

    cP = cos(phi_P)
    sP = sin(phi_P)
    cC = cos(phi_C)
    sC = sin(phi_C)

    PF_R_PJ = Matrix([[(sy/sqrt(sy**2+sx**2))**2*(1 - cP) + cP, -sx/sqrt(sy**2+sx**2)*sy/sqrt(sy**2+sx**2)*(1 - cP),      sx/sqrt(sy**2+sx**2)*sP],
                      [-sx/sqrt(sy**2+sx**2)*sy/sqrt(sy**2+sx**2)*(1 - cP),      (sx/sqrt(sy**2+sx**2))**2*(1 - cP) + cP, sy/sqrt(sy**2+sx**2)*sP],
                      [-sx/sqrt(sy**2+sx**2)*sP,               -sy/sqrt(sy**2+sx**2)*sP,                 cP]])
    PJ_R_CJ = Matrix([[cos(theta), -sin(theta), 0],
                      [sin(theta),  cos(theta), 0],
                      [0,                0,     1]])
    CF_R_CJ = Matrix([[(ty/sqrt(ty**2+tx**2))**2*(1 - cC) + cC, -tx/sqrt(ty**2+tx**2)*ty/sqrt(ty**2+tx**2)*(1 - cC),      tx/sqrt(ty**2+tx**2)*sC],
                      [-tx/sqrt(ty**2+tx**2)*ty/sqrt(ty**2+tx**2)*(1 - cC),      (tx/sqrt(ty**2+tx**2))**2*(1 - cC) + cC, ty/sqrt(ty**2+tx**2)*sC],
                      [-tx/sqrt(ty**2+tx**2)*sC,               -ty/sqrt(ty**2+tx**2)*sC,                 cC]])

    # Tests for the parent.frame
    pframe2pjointframe = PF_R_PJ
    pframe2cjointframe = PF_R_PJ * PJ_R_CJ

    assert simplify(pframe.dcm(pjointframe) - pframe2pjointframe) == zeros(3)
    assert simplify(pframe.dcm(cjointframe) - pframe2cjointframe) == zeros(3)

    assert pframe.ang_vel_in(pjointframe) == 0
    assert pframe.ang_vel_in(cjointframe) == -omega * pjointframe.z
    assert pframe.ang_vel_in(cframe) == -omega * pjointframe.z

    # Tests for the parent_joint_frame
    assert pjointframe.name == 'Joint_parent_joint_frame'

    pjointframe2cjointframe = PJ_R_CJ
    pjointframe2cframe = PJ_R_CJ * CF_R_CJ.transpose()
    cjointframe2cframe = CF_R_CJ.transpose()

    pprint(pframe.dcm(pjointframe))
    # pprint(simplify(cjointframe.dcm(cframe) - cjointframe2cframe))

    # pprint(simplify(cframe.dcm(cjointframe)))
    # pprint(simplify(pjointframe.dcm(cframe) - pjointframe2cframe))

    assert simplify(pjointframe.dcm(cjointframe) - pjointframe2cjointframe) == zeros(3)
    #assert simplify(pjointframe.dcm(cframe) - pjointframe2cframe) == zeros(3)

    assert pjointframe.ang_vel_in(pframe) == 0
    assert pjointframe.ang_vel_in(cjointframe) == -y * pjointframe.z
    assert pjointframe.ang_vel_in(cframe) == -y * pjointframe.z

    # Tests for the child_joint_frame
    assert cjointframe.name == 'Joint_child_joint_frame'

    cjointframe2cframe = CF_R_CJ.transpose()


    assert cjointframe.dcm(cframe) == cjointframe2cframe

    assert cjointframe.ang_vel_in(pframe) == y * pjointframe.z
    assert cjointframe.ang_vel_in(pjointframe) == y * pjointframe.z
    assert cjointframe.ang_vel_in(cframe) == 0

    # Tests for the child.frame
    assert cframe.ang_vel_in(pframe) == y * pjointframe.z
    assert cframe.ang_vel_in(pjointframe) == y * pjointframe.z
    assert cframe.ang_vel_in(cjointframe) == 0


def test_pinjoint_fullargs_XTchild():
    """Test the output of XT_child for a PinJoint for which all args were
    specified"""


def test_pinjoint_input_errors():
    """Test ways that input arguments can be incorrectly specified for a
    PinJoint"""


def test_slidingjoint_attributes():
    """Test that the attributes of SlidingJoint match the correct defaults"""
    x = dynamicsymbols('x')
    y = x.diff()
    parent = Body('parent')
    child = Body('child')
    joint = SlidingJoint('Joint', parent, child, x, y)


def test_slidingjoint_spatial():
    """Test the spatial components of a default SlidingJoint"""


def test_slidingjoint_defaults_points():
    """Test the relations between the 4 points in a default SlidingJoint"""


def test_slidingjoint_defaults_frames():
    """Test the relations between the 4 frames in a default SlidingJoint"""


def test_slidingjoint_defaults_XTchild():
    """Test the output of XT_child for a default SlidingJoint"""


def test_slidingjoint_fullargs_points():
    """Test the relations between the 4 points in a SlidingJoint for which all
    args were specified"""


def test_slidingjoint_fullargs_frames():
    """Test the relations between the 4 frames in a SlidingJoint for which all
    args were specified"""


def test_slidingjoint_fullargs_XTchild():
    """Test the output of XT_child for a SlidingJoint for which all args were
    specified"""


def test_slidingjoint_input_errors():
    """Test ways that input arguments can be incorrectly specified for a
    SlidingJoint"""
