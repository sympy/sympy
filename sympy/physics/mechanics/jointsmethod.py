from sympy.physics.mechanics import Body, Lagrangian, KanesMethod, LagrangesMethod
from sympy.physics.mechanics.method import _Methods

__all__ = ['JointsMethod']


class JointsMethod(_Methods):
    """Method for formulating the equations of motion using a set of interconnected bodies with joints.

    Parameters
    ==========

    newtonion : Body or ReferenceFrame
        The newtonion(inertial) frame.
    *joints : Joint
        The joints in the system

    Attributes
    ==========

    q, u : Matrix
        Matrices of the generalized coordinates and speeds
    bodies : iterable
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
    >>> method.form_eoms()
    (Matrix([[-c*u(t) - k*q(t)]]), Matrix([[-B_mass*Derivative(u(t), t)]]))
    >>> MM = method.mass_matrix_full
    >>> forcing = method.forcing_full
    >>> rhs = MM.LUsolve(forcing)
    >>> rhs
    Matrix([
    [                     u(t)],
    [(-c*u(t) - k*q(t))/B_mass]])

    Notes
    =====

    ``JointsMethod`` currently only works with systems that do not have any
    configuration or motion constraints.

    """

    def __init__(self, newtonion, *joints):
        if isinstance(newtonion, Body):
            self.frame = newtonion.frame
        else:
            self.frame = newtonion

        self._joints = joints
        self._bodylist = self._generate_bodylist()
        self._loadlist = self._generate_loadlist()
        self._q = self._generate_q()
        self._u = self._generate_u()
        self._kdes = self._generate_kdes()

        self._method = None

    @property
    def bodies(self):
        return self._bodylist

    @property
    def forcelist(self):
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
        for body in self.bodies:
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

    def form_eoms(self, method=KanesMethod):
        """Method to form system's equation of motions.

        Parameters
        ==========

        method : Class
            Class name of method.

        Examples
        ========

        This is a simple example for a one degree of freedom translational
        spring-mass-damper.

        >>> from sympy.physics.mechanics import LagrangesMethod, dynamicsymbols, Body
        >>> from sympy.physics.mechanics import PrismaticJoint, JointsMethod
        >>> from sympy import symbols
        >>> q = dynamicsymbols('q')
        >>> qd = dynamicsymbols('q', 1)
        >>> m, k, b = symbols('m k b')
        >>> B = Body('B')
        >>> P = Body('P', mass=m)
        >>> P.potential_energy = k * q**2 / 2.0
        >>> J = PrismaticJoint('J', B, P, coordinates=q, speeds=qd)
        >>> B.apply_force(b * qd * B.x, reaction_body=P)
        >>> method = JointsMethod(B, J)
        >>> method.form_eoms(LagrangesMethod)
        Matrix([[b*Derivative(q(t), t) + 1.0*k*q(t) + m*Derivative(q(t), (t, 2))]])

        We can also solve for the states using the 'rhs' method.

        >>> method.rhs()
        Matrix([
        [                    Derivative(q(t), t)],
        [(-b*Derivative(q(t), t) - 1.0*k*q(t))/m]])

        """

        if issubclass(method, KanesMethod): #KanesMethod or similar
            self._method = method(self.frame, q_ind=self.q, u_ind=self.u, kd_eqs=self.kdes,
                                    forcelist=self.forcelist, bodies=self.bodies)
        if issubclass(method, LagrangesMethod): #LagrangesMethod or similar
            L = Lagrangian(self.frame, *self.bodies)
            self._method = method(L, self.q, self.forcelist, self.bodies, self.frame)
        soln = self.method._form_eoms()
        return soln

    def rhs(self, inv_method=None):
        """Returns equations that can be solved numerically.

        Parameters
        ==========

        inv_method : str
            The specific sympy inverse matrix calculation method to use. For a
            list of valid methods, see
            :meth:`~sympy.matrices.matrices.MatrixBase.inv`
        """

        return self.method.rhs(inv_method=inv_method)
