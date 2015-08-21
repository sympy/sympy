from sympy.physics.mechanics import KanesMethod

__all__ = ['JointsMethod']


class JointsMethod(object):
    """
    Joints method object.

    This object is used to form the equations of motion of a system created
    using Joints.

    Parameters
    ----------
    joints: List of Joints
        List of all the joints created must be passed as argument.
    root_body: Body
        Root body w.r.t which equations of motion much be generated.
        Is necessary to remove the  ambiguity since user may not
        connect all bodies using joints and thus, there is not
        way to find the root_body.

    Examples
    --------
    >>> from sympy import Symbol
    >>> from sympy.physics.mechanics import (Body, PinJoint, SphericalJoint, \
                                             JointsMethod)
    >>> root_body = Body('root_body')
    >>> parent = Body('parent')
    >>> child = Body('child')
    >>> l1 = Symbol('l1')
    >>> joints = []
    >>> pin_joint = PinJoint('pin_joint' root_body, parent, \
                             child_point_pos=(0, l1, 0))
    >>> joints.append(pin_joint)
    >>> l2 = Symbol('l2')
    >>> spherical_joint = SphericalJoint('spherical_joint', parent, child, \
                                         child_point_pos=(0, l2, 0))
    >>> joints.append(spherical_joint)
    >>> JM = JointsMethod(joints, root_body)
    >>> JM.rhs()

    """
    def __init__(self, joints, root_body):
        self.joints = joints
        self.root_body = root_body
        self._bodylist = self._generate_bodylist()
        self._forcelist = self._generate_forcelist()
        self._q = self._generate_q()
        self._u = self._generate_u()
        self._kds = self._generate_kds()
        self._set_kanes()

    @property
    def bodylist(self):
        return self._bodylist

    def _generate_bodylist(self):
        bodies = []
        for joint in self.joints:
            if joint.child not in bodies:
                bodies.append(joint.child)
            if joint.parent not in bodies:
                bodies.append(joint.parent)
        return bodies

    @property
    def forcelist(self):
        return self._forcelist

    def _generate_forcelist(self):
        force_list = []
        for body in self.bodylist:
            for force in body.force_list:
                force_list.append(force)
        return force_list

    @property
    def q(self):
        return self._q

    def _generate_q(self):
        q_ind = []
        for joint in self.joints:
            coordinates = joint.coordinates
            for coordinate in coordinates:
                q_ind.append(coordinate)
        return q_ind

    @property
    def u(self):
        return self._u

    def _generate_u(self):
        u_ind = []
        for joint in self.joints:
            speeds = joint.speeds
            for speed in speeds:
                u_ind.append(speed)
        return u_ind

    @property
    def kd(self):
        return self._kds

    def _generate_kds(self):
        kd_ind = []
        for joint in self.joints:
            kds = joint.kds
            for kd in kds:
                kd_ind.append(kd)
        return kd_ind

    @property
    def forcing_full(self):
        return self._KM.forcing_full

    @property
    def mass_matrix_full(self):
        return self._KM.mass_matrix_full

    def _set_kanes(self):
        self._KM = KanesMethod(self.root_body.get_frame(), q_ind=self.q, u_ind=self.u,
                               kd_eqs=self.kd)
        self._KM.kanes_equations(self.forcelist, self.bodylist)
        # TODO Removing call to private attributes in pydy.System and fix this.
        self._qdot = self._KM._qdot
        self._udot = self._KM._udot
        self._uaux = self._KM._uaux

    def rhs(self, inv_method=None):
        return self._KM.rhs(inv_method)
