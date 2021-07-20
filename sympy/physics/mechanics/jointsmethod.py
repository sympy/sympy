from sympy.physics.mechanics import KanesMethod, Body

__all__ = ['JointsMethod']


class JointsMethod(object):
    """Method for formulating the equations of motion using a set of interconnected bodies with joints.

    Parameters
    ==========

    body : Body or ReferenceFrame
        The inertial(Newtonion) frame.
    *joints : Joint
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
    mass_matrix : Matrix
        The system's mass matrix
    forcing : Matrix
        The system's forcing vector
    mass_matrix_full : Matrix
        The "mass matrix" for the u's and q's
    forcing_full : Matrix
        The "forcing vector" for the u's and q's
    method : KanesMethod or Lagrange's method
        Method's object.

    Examples
    ========

    This is a simple example for a one degree of freedom translational
    spring-mass-damper.

    >>> from sympy import symbols
    >>> from sympy.physics.mechanics import Body, JointsMethod, PrismaticJoint
    >>> from sympy.physics.vector import dynamicsymbols
    >>> c, k = symbols('c k')
    >>> q, u = dynamicsymbols('q u')
    >>> W = Body('W')
    >>> B = Body('B')
    >>> J = PrismaticJoint('J', W, B, coordinates=q, speeds=u)
    >>> W.apply_force(c*u*W.x, reaction_body=B)
    >>> W.apply_force(k*q*W.x, reaction_body=B)
    >>> method = JointsMethod(W, J)
    >>> (fr, frstar) = method.form_eoms()
    >>> MM = method.mass_matrix
    >>> forcing = method.forcing
    >>> rhs = MM.inv() * forcing
    >>> rhs
    Matrix([[(-c*u(t) - k*q(t))/B_mass]])

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
        return self.method.forcing_full

    @property
    def mass_matrix_full(self):
        return self.method.mass_matrix_full

    @property
    def mass_matrix(self):
        return self.method.mass_matrix

    @property
    def forcing(self):
        return self.method.forcing

    @property
    def method(self):
        return self._method

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

    def _form_kanes_equations(self):
        self._method = KanesMethod(self.frame, q_ind=self.q,
                                u_ind=self.u, kd_eqs=self.kdes)
        if self.loadlist == []:
           fr, frstar = self.method.kanes_equations(self.bodylist)
           return fr, frstar
        fr, frstar = self.method.kanes_equations(self.bodylist, self.loadlist)
        return fr, frstar

    def _form_lagranges_equations(self):
        pass

    def form_eoms(self, method='kane'):
        """ Method to form system's equation of motions.

        Parameters
        ==========

        method : String
            - `kane` : Form equations of motions using kane's method. Default method.
            - `lagrange` : Form equations of motions using lagrange's method.

        """

        if method == 'kane':
            (fr, frstar) = self._form_kanes_equations()
            return fr, frstar
        if method == 'lagrange':
            lagr_eqns = self._form_lagranges_equations()
            return lagr_eqns
