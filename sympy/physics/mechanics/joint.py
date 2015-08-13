from sympy import acos
from sympy.physics.vector import cross, Vector, dot, dynamicsymbols
from sympy.physics.mechanics.functions import convert_tuple_to_vector

__all__ = ['Joint', 'PinJoint', 'SlidingJoint', 'CylindricalJoint',
           'SphericalJoint', 'PlanarJoint']


class Joint(object):
    """Abstract Base class for all specific joints.

    Joints connects two bodies (a parent and child) by adding degrees of freedom
    to child w.r.t parent. Different joints adds different degrees of freedom.
    This is the base class for all specific joints and holds all common methods
    acting as an interface for all joints. Custom joint can also be created
    by creating a subclass of this class and overriding 'apply_joint()'
    method to add the dynamics of the joint.

    Parameters
    ----------
    name: String
        Name of the joint which makes it unique. Should be different from other
        joints.
    parent: Body
        Instance of Body class serving as the parent in the joint.
    child: Body
        Instance of Body class serving as the child in the joint.
    parent_point_pos: 3 Tuple, Optional
        Defines the joint's point where the parent will be connected to child.
        3 Tuple defines the values of x, y and z directions w.r.t parent's
        frame. If it is not supplied, center of mass is added as default.
    child_point_pos: 3 Tuple, Optional
        Defines the joint's point where the child will be connected to parent.
        3 Tuple defines the values of x, y and z directions w.r.t child's frame.
        If it is not supplied, center of mass is added as default.
    parent_axis: Vector, Optional
        Defines the orientation as a vector which must be aligned with child's
        axis before adding joint. If it is not passed, default is x axis in
        parent's frame.
    child_axis: Vector, Optional
        Defines the orientation as a vector which must be aligned with parent's
        axis before adding joint. If it is not passed, default is x axis in
        child's frame.
    """

    def __init__(self, name, parent, child, parent_point_pos=None,
                 child_point_pos=None):
        self.name = name
        self.parent = parent
        self.child = child
        self.coordinates = []
        self.speeds = []
        self.kds = []

        if parent_point_pos is None:
            parent_point_pos = (0, 0, 0)  # Parent's Center of mass

        if child_point_pos is None:
            child_point_pos = (0, 0, 0)  # Child's Center of mass

        self.parent_joint_vector = convert_tuple_to_vector(
            parent.frame, parent_point_pos)
        self.child_joint_vector = convert_tuple_to_vector(
            child.frame, child_point_pos)

        self._locate_joint_point()
        self.apply_joint()

    def _locate_joint_point(self):
        self.parent_joint_point = self.parent.masscenter.locatenew(
            self.name + '_parent_joint',
            self.parent_joint_vector)
        self.child_joint_point = self.child.masscenter.locatenew(
            self.name + '_child_joint',
            self.child_joint_vector)

    def _align_axes(self, parent_axis, child_axis):
        """Rotates child_frame so that child_axis is aligned to parent_axis."""
        mag1 = parent_axis.magnitude()
        mag2 = child_axis.magnitude()
        angle = acos(dot(parent_axis, child_axis)/(mag1 * mag2))
        axis = cross(child_axis, parent_axis)
        if axis != Vector(0):
            self.child.frame.orient(
                self.parent.frame, 'Axis', [angle, axis])

    def apply_joint(self):
        """To create a custom joint, this method should add degrees of freedom
        to the bodies in subclass of Joint"""
        raise NotImplementedError("To define a custom pydy.Joint, you need to" +
                                  " override apply_joint method in Joint's" +
                                  " subclass.")


class PinJoint(Joint):
    """
    Pin (Revolute) Joint.

    Provides one rotational degree of freedom to the child w.r.t parent. Single
    generalized coordinate, 'theta', is the rotational angle and generalized
    speed, 'omega', is the angular velocity. Generalized speed is the time derivative of
    the generalized coordinate.

    Parameters
    ----------
    name: String
        Name of the joint which makes it unique. Should be different from other
        joints.
    parent: Body
        Instance of Body class serving as the parent in the joint.
    child: Body
        Instance of Body class serving as the child in the joint.
    parent_point_pos: 3 Tuple, Optional
        Defines the joint's point where the parent will be connected to child.
        3 Tuple defines the values of x, y and z directions w.r.t parent's
        frame. If it is not supplied, center of mass is added as default.
    child_point_pos: 3 Tuple, Optional
        Defines the joint's point where the child will be connected to parent.
        3 Tuple defines the values of x, y and z directions w.r.t child's frame.
        If it is not supplied, center of mass is added as default.
    parent_axis: Vector, Optional
        Defines the orientation as a vector which must be aligned with child's
        axis before adding joint. If it is not passed, default is x axis in
        parent's frame.
    child_axis: Vector, Optional
        Defines the orientation as a vector which must be aligned with parent's
        axis before adding joint. If it is not passed, default is x axis in
        child's frame.

    Example
    -------
    Adding a Pin Joint which connects center of mass of parent to a point
    pointed by child.frame.x + child.frame.y . Gravity is along y axis of
    parent.
    >>> from sympy import Symbol
    >>> from sympy.physics.mechanics import Body, PinJoint
    >>> parent = Body('parent')
    >>> child = Body('child')
    >>> gravity = Symbol('gravity')
    >>> child.add_force(child.mass * gravity * parent.frame.y, child.masscenter)
    >>> pin_joint = PinJoint('pin_joint', parent, \
                             child, child_point_pos=(1, 1, 0))
    >>> pin_joint.coordinates
    [theta(t)]
    >>> pin_joint.speeds
    [omega(t)]

    """
    def __init__(self, name, parent, child, parent_point_pos=None,
                 child_point_pos=None, parent_axis=None, child_axis=None):

        if parent_axis is None:
            self.parent_axis = parent.frame.x
        else:
            self.parent_axis = convert_tuple_to_vector(parent.frame,
                                                        parent_axis)
        if child_axis is None:
            self.child_axis = child.frame.x
        else:
            self.child_axis = convert_tuple_to_vector(child.frame,
                                                       child_axis)

        super(PinJoint, self).__init__(name, parent, child, parent_point_pos,
                                       child_point_pos)

    def apply_joint(self):
        theta = dynamicsymbols(self.name + '_theta')
        thetad = dynamicsymbols(self.name + '_theta', 1)
        omega = dynamicsymbols(self.name + '_omega')

        self.coordinates.append(theta)
        self.speeds.append(omega)
        self.kds.append(thetad - omega)

        self.child.frame.orient(self.parent.frame, 'Axis',
                                [theta, self.parent_axis])
        self.child.frame.set_ang_vel(self.parent.frame,
                                     omega * self.parent_axis)
        self._align_axes(self.parent_axis, self.child_axis)
        self.child_joint_point.set_pos(self.parent_joint_point, 0)
        self.child.masscenter.v2pt_theory(self.parent.masscenter,
                                          self.parent.frame, self.child.frame)


class SlidingJoint(Joint):
    """
    Sliding (Prismatic) Joint.

    Provides one translational degree of freedom to the child w.r.t parent. Single
    generalized coordinate, 'dis', is the displacement and generalized
    speed, 'vel', is the velocity. Generalized speed is the time derivative of
    the generalized coordinate.

    Parameters
    ----------
    name: String
        Name of the joint which makes it unique. Should be different from other
        joints.
    parent: Body
        Instance of Body class serving as the parent in the joint.
    child: Body
        Instance of Body class serving as the child in the joint.
    parent_point_pos: 3 Tuple, Optional
        Defines the joint's point where the parent will be connected to child.
        3 Tuple defines the values of x, y and z directions w.r.t parent's
        frame. If it is not supplied, center of mass is added as default.
    child_point_pos: 3 Tuple, Optional
        Defines the joint's point where the child will be connected to parent.
        3 Tuple defines the values of x, y and z directions w.r.t child's frame.
        If it is not supplied, center of mass is added as default.
    parent_axis: Vector, Optional
        Defines the orientation as a vector which must be aligned with child's
        axis before adding joint. If it is not passed, default is x axis in
        parent's frame.
    child_axis: Vector, Optional
        Defines the orientation as a vector which must be aligned with parent's
        axis before adding joint. If it is not passed, default is x axis in
        child's frame.

    Example
    -------
    Adds sliding Joint between parent's masscenter and a point located at unit
    distance in x axis of child.
    >>> from sympy import
    >>> from sympy.physics.mechanics import Body, SlidingJoint
    >>> parent = Body('parent')
    >>> child = Body('child')
    >>> sliding_joint = SlidingJoint('slidingjoint', parent, child, \
                                     child_point_pos=(1, 0, 0))
    >>> sliding_joint.coordinates
    [dis(t)]
    >>> sliding_joint.speeds
    [vel(t)]

    """
    def __init__(self, name, parent, child, parent_point_pos, child_point_pos,
                 parent_axis=None, child_axis=None):

        if parent_axis is None:
            self.parent_axis = parent.frame.x
        else:
            self.parent_axis = convert_tuple_to_vector(parent.frame,
                                                        parent_axis)
        if child_axis is None:
            self.child_axis = child.frame.x
        else:
            self.child_axis = convert_tuple_to_vector(child.frame,
                                                       child_axis)

        super(SlidingJoint, self).__init__(name, parent, child,
                                           parent_point_pos, child_point_pos)

    def apply_joint(self):
        dis = dynamicsymbols(self.name + '_dis')
        disd = dynamicsymbols(self.name + '_dis', 1)
        vel = dynamicsymbols(self.name + '_vel')

        self.coordinates.append(dis)
        self.speeds.append(vel)
        self.kds.append(disd - vel)

        self.child.frame.orient(self.parent.frame, 'Axis',
                                [0, self.parent.frame.z])
        self._align_axes(self.parent_axis, self.child_axis)
        self._locate_joint_point()

        self.parent_joint_point.set_vel(self.parent.frame, 0)
        self.child_joint_point.set_vel(self.child.frame, 0)

        self.child_joint_point.set_pos(self.parent_joint_point,
                                       dis * self.parent_axis)
        self.child_joint_point.set_vel(self.parent.frame,
                                       vel * self.parent_axis)
        self.child.masscenter.set_vel(self.parent.frame,
                                      vel * self.parent_axis)


class CylindricalJoint(Joint):
    """
    Cylindrical Joint.

    Provides one translational and one rotational degree of freedom to the
    child w.r.t parent. It is a combination of Pin Joint and Sliding Joint.
    Adds two generalized coordinates, theta - rotational angle and dis -
    displacement. Two generalized speeds are omega - angular velocity and
    vel - linear velocity. Generalized speeds are the time derivatives of
    respective generalized coordinates.

    Parameters
    ----------
    name: String
        Name of the joint which makes it unique. Should be different from other
        joints.
    parent: Body
        Instance of Body class serving as the parent in the joint.
    child: Body
        Instance of Body class serving as the child in the joint.
    parent_point_pos: 3 Tuple, Optional
        Defines the joint's point where the parent will be connected to child.
        3 Tuple defines the values of x, y and z directions w.r.t parent's
        frame. If it is not supplied, center of mass is added as default.
    child_point_pos: 3 Tuple, Optional
        Defines the joint's point where the child will be connected to parent.
        3 Tuple defines the values of x, y and z directions w.r.t child's frame.
        If it is not supplied, center of mass is added as default.
    parent_axis: Vector, Optional
        Defines the orientation as a vector which must be aligned with child's
        axis before adding joint. If it is not passed, default is x axis in
        parent's frame.
    child_axis: Vector, Optional
        Defines the orientation as a vector which must be aligned with parent's
        axis before adding joint. If it is not passed, default is x axis in
        child's frame.

    Example
    -------
    Adds cylindrical Joint between parent's masscenter and a point located by
    child.frame.x + child.frame.y in child.
    >>> from sympy.physics.mechanics import Body, CylindricalJoint
    >>> parent = Body('parent')
    >>> child = Body('child')
    >>> cylindrical_joint = CylindricalJoint('cylindricaljoint', parent, child, \
                                             child_point_pos=(1, 1, 0))
    >>> cylindrical_joint.coordinates
    [theta(t), dis(t)]
    >>>cylindrical_joint.speeds
    [omega(t), vel(t)]

    """
    def __init__(self, name, parent, child, parent_point_pos=None,
                 child_point_pos=None, parent_axis=None, child_axis=None):

        if parent_axis is None:
            self.parent_axis = parent.frame.x
        else:
            self.parent_axis = convert_tuple_to_vector(parent.frame,
                                                        parent_axis)
        if child_axis is None:
            self.child_axis = child.frame.x
        else:
            self.child_axis = convert_tuple_to_vector(child.frame,
                                                       child_axis)

        super(CylindricalJoint, self).__init__(name, parent, child,
                                               parent_point_pos, child_point_pos)

    def apply_joint(self):
        dis = dynamicsymbols(self.name + '_dis')
        disd = dynamicsymbols(self.name + '_dis', 1)
        vel = dynamicsymbols(self.name + '_vel')

        theta = dynamicsymbols(self.name + '_theta')
        thetad = dynamicsymbols(self.name + '_theta', 1)
        omega = dynamicsymbols(self.name + '_omega')

        self.coordinates.append(dis)
        self.speeds.append(vel)
        self.kds.append(disd - vel)

        self.coordinates.append(theta)
        self.speeds.append(omega)
        self.kds.append(thetad - omega)

        self.child.frame.orient(self.parent.frame, 'Axis',
                                [0, self.parent.frame.z])
        self._align_axes(self.parent_axis, self.child_axis)
        self._locate_joint_point()
        self.child.frame.orient(self.parent.frame, 'Axis',
                                [theta, self.parent_axis])
        self.child.frame.set_ang_vel(self.parent.frame,
                                     omega * self.parent_axis)

        self.parent_joint_point.set_vel(self.parent.frame, 0)
        self.child_joint_point.set_vel(self.child.frame, 0)

        self.child_joint_point.set_pos(self.parent_joint_point,
                                       dis * self.parent_axis)
        self.child_joint_point.set_vel(self.parent.frame,
                                       vel * self.parent_axis)
        self.child.masscenter.set_vel(self.parent.frame,
                                      vel * self.parent_axis)
        self.child_joint_point.set_pos(self.parent_joint_point, 0)
        self.child.masscenter.v2pt_theory(self.parent.masscenter,
                                          self.parent.frame,
                                          self.child.frame)


class PlanarJoint(Joint):
    """
    Planar Joint.

    Provides two translational and one rotational degrees of freedom to the
    child w.r.t parent. Translational is in a plane or along two perpendicular
    vectors. The generalized coordinates are theta - the rotational angle,
    disx - displacement along one of the perpendicular vector and disy -
    displacement along the other perpendicular vector . The generalized speeds
    are omega - angular velocity perpendicular to the plane, velx - linear
    velocity along disx and vely - linear velocity along disy.

    Parameters
    ----------
    name: String
        Name of the joint which makes it unique. Should be different from other
        joints.
    parent: Body
        Instance of Body class serving as the parent in the joint.
    child: Body
        Instance of Body class serving as the child in the joint.
    parent_point_pos: 3 Tuple, Optional
        Defines the joint's point where the parent will be connected to child.
        3 Tuple defines the values of x, y and z directions w.r.t parent's
        frame. If it is not supplied, center of mass is added as default.
    child_point_pos: 3 Tuple, Optional
        Defines the joint's point where the child will be connected to parent.
        3 Tuple defines the values of x, y and z directions w.r.t child's frame.
        If it is not supplied, center of mass is added as default.
    parent_axis: Vector, Optional
        Defines the orientation as a vector which must be aligned with child's
        axis before adding joint. If it is not passed, default is x axis in
        parent's frame.
    child_axis: Vector, Optional
        Defines the orientation as a vector which must be aligned with parent's
        axis before adding joint. If it is not passed, default is x axis in
        child's frame.

    Example
    -------
    Adds planar Joint between parent's masscenter and a point at unit distance
    along child's y axis.
    >>> from sympy.physics.mechanics import Body, PlanarJoint
    >>> parent = Body('parent')
    >>> child = Body('child')
    >>> planar_joint = PlanarJoint('planarjoint', parent, child, \
                                   child_point_pos=(0, 1, 0))
    >>> planar_joint.coordinates
    [theta(t), disx(t), disy(t)]
    >>> planar_joint.speeds
    [omega(t), velx(t), vely(t)]

    """
    def __init__(self, name, parent, child, parent_point_pos=None,
                 child_point_pos=None, parent_axis=None, child_axis=None):

        if parent_axis is None:
            self.parent_axis = parent.frame.x
        else:
            self.parent_axis = convert_tuple_to_vector(parent.frame,
                                                        parent_axis)
        if child_axis is None:
            self.child_axis = child.frame.x
        else:
            self.child_axis = convert_tuple_to_vector(child.frame,
                                                       child_axis)

        super(PlanarJoint, self).__init__(name, parent, child,
                                          parent_point_pos, child_point_pos)

    def apply_joint(self):
        # generalized coordinates in specific order.
        theta = dynamicsymbols(self.name + '_theta')  # rotation around z axis.
        thetad = dynamicsymbols(self.name + '_theta', 1)
        omega = dynamicsymbols(self.name + '_omega')
        disx = dynamicsymbols(self.name + '_disx')  # translation along x axis.
        disxd = dynamicsymbols(self.name + '_disx', 1)
        velx = dynamicsymbols(self.name + '_velx')
        disy = dynamicsymbols(self.name + '_disy')  # translation along y axis.
        disyd = dynamicsymbols(self.name + '_disy', 1)
        vely = dynamicsymbols(self.name + '_vely')

        self.coordinates.append(theta)
        self.speeds.append(omega)
        self.kds.append(thetad - omega)

        self.coordinates.append(disx)
        self.speeds.append(velx)
        self.kds.append(disxd - velx)

        self.coordinates.append(disy)
        self.speeds.append(vely)
        self.kds.append(disyd - vely)
        
        self.child.frame.orient(self.parent.frame, 'Axis',
                                [0, self.parent.frame.x])
        self.child_joint_point.set_pos(self.parent_joint_point, 0)
        self._align_axes(self.parent_axis, self.child_axis)
        self._locate_joint_point()

        # Adding rotation
        self.child.frame.orient(self.parent.frame, 'Axis',
                                [theta, self.parent.frame.z])
        self.child.frame.set_ang_vel(self.parent.frame,
                                     omega * self.parent.frame.z)
        
        # Adding translation along x axis.
        self.child_joint_point.set_pos(self.parent_joint_point,
                                       disx * self.parent.frame.x)
        self.child_joint_point.set_vel(self.parent.frame,
                                       velx * self.parent.frame.x)
        
        # Adding translation along y axis
        self.child_joint_point.set_pos(self.parent_joint_point,
                                       disy * self.parent.frame.y)
        self.child_joint_point.set_vel(self.parent.frame,
                                       vely * self.parent.frame.y)
        
        self.child.masscenter.v2pt_theory(self.parent.masscenter,
                                          self.parent.frame, self.child.frame)


class SphericalJoint(Joint):
    """
    Spherical (Ball) Joint.

    Provides three rotational degrees of freedom to the child w.r.t parent.
    The three generalized coordinates (thetax, thetay and thetax) are always
    the three angles along three perpendicular vectors and the three
    generalized speeds (omegax, omegay and omegaz) are the angular velocities
    along those vectors respectively.

    Parameters
    ----------
    name: String
        Name of the joint which makes it unique. Should be different from other
        joints.
    parent: Body
        Instance of Body class serving as the parent in the joint.
    child: Body
        Instance of Body class serving as the child in the joint.
    parent_point_pos: 3 Tuple, Optional
        Defines the joint's point where the parent will be connected to child.
        3 Tuple defines the values of x, y and z directions w.r.t parent's
        frame. If it is not supplied, center of mass is added as default.
    child_point_pos: 3 Tuple, Optional
        Defines the joint's point where the child will be connected to parent.
        3 Tuple defines the values of x, y and z directions w.r.t child's frame.
        If it is not supplied, center of mass is added as default.
    parent_axis: Vector, Optional
        Defines the orientation as a vector which must be aligned with child's
        axis before adding joint. If it is not passed, default is x axis in
        parent's frame.
    child_axis: Vector, Optional
        Defines the orientation as a vector which must be aligned with parent's
        axis before adding joint. If it is not passed, default is x axis in
        child's frame.

    Example
    -------
    Adds spherical Joint between parent's masscenter and a point at unit distance
    along child's y axis.
    >>> from sympy.physics.mechanics import Body, SphericalJoint
    >>> parent = Body('parent')
    >>> child = Body('child')
    >>> spherical_joint = SphericalJoint('sphericaljoint', parent, child, \
                                         child_point_pos=(0, 1, 0))
    >>> spherical_joint.coordinates
    [thetax(t), thetay(t), thetaz(t)]
    >>> spherical_joint.speeds
    [omegax(t), omegay(t), omegaz(t)]

    """
    def __init__(self, name, parent, child, parent_point_pos=None,
                 child_point_pos=None):

        super(SphericalJoint, self).__init__(name, parent, child,
                                             parent_point_pos, child_point_pos)

    def apply_joint(self):
        thetax = dynamicsymbols(self.name + '_thetax')
        thetay = dynamicsymbols(self.name + '_thetay')
        thetaz = dynamicsymbols(self.name + '_thetaz')
        thetaxd = dynamicsymbols(self.name + '_thetax', 1)
        thetayd = dynamicsymbols(self.name + '_thetay', 1)
        thetazd = dynamicsymbols(self.name + '_thetaz', 1)
        omegax = dynamicsymbols(self.name + '_omegax')
        omegay = dynamicsymbols(self.name + '_omegay')
        omegaz = dynamicsymbols(self.name + '_omegaz')

        self.coordinates.append(thetax)
        self.speeds.append(omegax)
        self.kds.append(thetaxd - omegax)

        self.coordinates.append(thetay)
        self.speeds.append(omegay)
        self.kds.append(thetayd - omegay)

        self.coordinates.append(thetaz)
        self.speeds.append(omegaz)
        self.kds.append(thetazd - omegaz)

        self._locate_joint_point()
        self.child_joint_point.set_pos(self.parent_joint_point, 0)
        
        self.child.frame.orient(self.parent.frame, 'Axis',
                                [thetax, self.parent.frame.x])
        self.child.frame.set_ang_vel(self.parent.frame,
                                     omegax * self.parent.frame.x)
        self.child.frame.orient(self.parent.frame, 'Axis',
                                [thetay, self.parent.frame.y])
        self.child.frame.set_ang_vel(self.parent.frame,
                                     omegay * self.parent.frame.y)
        self.child.frame.orient(self.parent.frame, 'Axis',
                                [thetaz, self.parent.frame.z])
        self.child.frame.set_ang_vel(self.parent.frame,
                                     omegaz * self.parent.frame.z)
        self.child.masscenter.v2pt_theory(self.parent.masscenter,
                                          self.parent.frame,
                                          self.child.frame)
