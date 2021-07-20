from sympy.physics.mechanics import KanesMethod, Body

__all__ = ['JointsMethod']


class JointsMethod(object):
    """ Joints Method

    Explanation
    ===========

    A joints method object is used to form the equations of
    motion of a system created with joints.

    Parameters
    ==========

    body : Body or ReferenceFrame
        The inertial(Newtonion) frame.
    joints : Joint
        The joints in the system

    Attributes
    ==========

    q, u : Matrix
        Matrices of the generalized coordinates and speeds
    bodylist : iterable
        Iterable of Body objects in the system.
    forcelist : iterable
        Iterable of (Point, vector) or (ReferenceFrame, vector) tuples
        describing the forces on the system.
    auxiliary_eqs : Matrix
        If applicable, the set of auxiliary Kane's
        equations used to solve for non-contributing
        forces.
    mass_matrix : Matrix
        The system's mass matrix
    forcing : Matrix
        The system's forcing vector
    mass_matrix_full : Matrix
        The "mass matrix" for the u's and q's
    forcing_full : Matrix
        The "forcing vector" for the u's and q's
    kane : KanesMethod
        Kane's method object used.

    """

    def __init__(self, body, *joints):
        if isinstance(body, Body):
            self.frame = body.frame
        else:
            self.frame = body

        self._joints = joints
        self._bodylist = self._generate_bodylist()
        self._loadlist = self._generate_loadlist()
        self._q = self._generate_q()
        self._u = self._generate_u()
        self._kdes = self._generate_kdes()
        self._kane = KanesMethod(self.frame, q_ind=self.q,
                                u_ind=self.u, kd_eqs=self.kdes)

    @property
    def bodylist(self):
        return self._bodylist

    @property
    def loadlist(self):
        return self._loadlist

    @property
    def q(self):
        return self._q

    @property
    def u(self):
        return self._u

    @property
    def kdes(self):
        return self._kdes

    @property
    def forcing_full(self):
        return self._kane.forcing_full

    @property
    def mass_matrix_full(self):
        return self._kane.mass_matrix_full

    @property
    def fr(self):
        return self._fr

    @property
    def frstar(self):
        return self._frstar

    @property
    def mass_matrix(self):
        return self._kane.mass_matrix

    @property
    def forcing(self):
        return self._kane.forcing

    @property
    def auxiliary_eqs(self):
        return self._kane.auxiliary_eqs

    @property
    def kane(self):
        return self._kane

    def _generate_bodylist(self):
        bodies = []
        for joint in self._joints:
            if joint.child not in bodies:
                bodies.append(joint.child)
            if joint.parent not in bodies:
                bodies.append(joint.parent)
        return bodies

    def _generate_loadlist(self):
        load_list = []
        for body in self.bodylist:
            load_list.extend(body.loads)
        return load_list

    def _generate_q(self):
        q_ind = []
        for joint in self._joints:
            q_ind.extend(joint.coordinates)
        return q_ind

    def _generate_u(self):
        u_ind = []
        for joint in self._joints:
            u_ind.extend(joint.speeds)
        return u_ind

    def _generate_kdes(self):
        kd_ind = []
        for joint in self._joints:
            kd_ind.extend(joint.kdes)
        return kd_ind

    def rhs(self, inv_method=None):
        """Returns the system's equations of motion in first order form. The
        output is the right hand side of::

           x' = |q'| =: f(q, u, r, p, t)
                |u'|

        The right hand side is what is needed by most numerical ODE
        integrators.

        Parameters
        ==========

        inv_method : str
            The specific sympy inverse matrix calculation method to use. For a
            list of valid methods, see
            :meth:`~sympy.matrices.matrices.MatrixBase.inv`

        """

        return self._kane.rhs(inv_method)

    def kanes_equations(self):
        """ Method to form Kane's equations, Fr + Fr* = 0.

        Explanation
        ===========

        Returns (Fr, Fr*). In the case where auxiliary generalized speeds are
        present (say, s auxiliary speeds, o generalized speeds, and m motion
        constraints) the length of the returned vectors will be o - m + s in
        length. The first o - m equations will be the constrained Kane's
        equations, then the s auxiliary Kane's equations. These auxiliary
        equations can be accessed with the auxiliary_eqs().

        """

        if self.loadlist == []:
           self._fr, self._frstar = self._kane.kanes_equations(self.bodylist)
           return self.fr, self.frstar
        self._fr, self._frstar = self._kane.kanes_equations(self.bodylist, self.loadlist)
        return self.fr, self.frstar
