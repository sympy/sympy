from __future__ import print_function, division

from sympy import zeros, Matrix, diff
from sympy.core.compatibility import range
from sympy.utilities import default_sort_key
from sympy.physics.vector import (ReferenceFrame, dynamicsymbols,
                                  partial_velocity)
from sympy.physics.mechanics.particle import Particle
from sympy.physics.mechanics.rigidbody import RigidBody
from sympy.physics.mechanics.functions import (msubs, find_dynamicsymbols,
                                               _f_list_parser)
from sympy.physics.mechanics.linearize import Linearizer
from sympy.utilities.exceptions import SymPyDeprecationWarning
from sympy.utilities.iterables import iterable

__all__ = ['KanesMethod']


def _none_handler(x):
    """Returns and empty matrix if x is False."""
    if x:
        return Matrix(x)
    else:
        return Matrix()


class KanesMethod(object):
    """This class is used automatically form the equations of motion of a
    multibody system using Kane's method [Kane1985]_.

    The dynamical equations of motion are provided in Kane's form::

        fr + frstar = 0

    in the factored implicit form::

        mass_matrix * speeds.diff() = forcing

    and in the explicit form::

        speeds.diff() = mass_matrix.inv() * forcing

    where fr, frstar, the mass matrix, M, and the right hand side vector, f,
    are available as class attributes after Kane's equations are formed.

    Attributes
    ==========

    inertial_frame : ReferenceFrame
        The inertial reference frame the system is defined with respect to.
    n, num_coordinates, num_speeds : integer
        The total number of generalized coordinates and speeds.
    q, coordinates : Matrix, shape(n, 1)
        The generalized coordinates. These are ordered as [q_s, q_r]^T.
    o, num_dependent_coordinates : integer
        The number of dependent coordinates.
    q_s, independent_coodrinates : Matrix, shape(n-o, 1)
        The independent generalized coordinates.
    q_r, dependent_coodrinates : Matrix, shape(o, 1)
        The dependent generalized coordinates.
    u, speeds : Matrix, shape(n, 1)
        The generalized speeds. For a nonholomic system these are ordered as
        [u_s, u_r]^T.
    p, num_indepenent_speeds : integer
        The number of independen speeds.
    u_s, independent_speeds : Matrix, shape(p, 1)
        The independent generalized speeds.
    m, num_dependent_speeds : integer
        The number of dependent speeds.
    u_r, dependent_speeds : Matrix, shape(m, 1)
        The dependent generalized speeds.
    holonomic : boolean
        True if the system is a holonomic system.
    nonholonomic : boolean
        True if the system is a nonholonomic system.
    bodies, bodylist : list
        The Particle and RigidBody objects in the system.
    loads, forcelist : list
        The (Point, Vector) or (ReferenceFrame, Vector) tuples describing the
        forces and torques on the system.
    auxiliary : Matrix
        If applicable, the set of auxiliary Kane's equations used to solve for
        non-contributing forces.
    fr : Matrix, shape(n, 1)
        The generalized active forces acting on the system.
    frstar : Matrix, shape(n, 1)
        The genealized inertia forcies acting on the system.
    mass_matrix : Matrix, shape(n, n)
        The "mass matrix" of the implicit dynamic differential equations.
    forcing : Matrix, shape(n, 1)
        The right hand side or "forcing vector" of the implicit dynamic
        differential equations.
    mass_matrix_full : Matrix, shape(2n, 2n)
        The "mass matrix" of the full first order form of the implicit dynamic
        and kinematic differential equations.
    forcing_full : Matrix, shape(2n, 1)
        The right hand side or "forcing vector" of the full first order form of
        the implicit dynamic and kinematic differential equations.

    Examples
    ========

    This is a simple example for a one degree of freedom translational
    spring-mass-damper.

    In this example, we first need to do the kinematics.  This involves
    creating generalized speeds and coordinates and their derivatives.  Then we
    create a point and set its velocity in a frame.

        >>> from sympy import symbols
        >>> from sympy.physics.mechanics import dynamicsymbols, ReferenceFrame
        >>> from sympy.physics.mechanics import Point, Particle, KanesMethod
        >>> q, u = dynamicsymbols('q u')
        >>> qd, ud = dynamicsymbols('q u', 1)
        >>> m, c, k = symbols('m c k')
        >>> N = ReferenceFrame('N')
        >>> P = Point('P')
        >>> P.set_vel(N, u * N.x)

    Next we need to arrange/store information in the way that KanesMethod
    requires.  The kinematic differential equations need to be stored in a
    dict.  A list of forces/torques must be constructed, where each entry in
    the list is a (Point, Vector) or (ReferenceFrame, Vector) tuple, where the
    Vectors represent the Force or Torque.
    Next a particle needs to be created, and it needs to have a point and mass
    assigned to it.
    Finally, a list of all bodies and particles needs to be created.

    >>> kd = [qd - u]
    >>> FL = [(P, (-k * q - c * u) * N.x)]
    >>> pa = Particle('pa', P, m)
    >>> BL = [pa]

    Finally we can generate the equations of motion.
    First we create the KanesMethod object and supply an inertial frame,
    coordinates, generalized speeds, and the kinematic differential equations.
    Additional quantities such as configuration and motion constraints,
    dependent coordinates and speeds, and auxiliary speeds are also supplied
    here (see the online documentation).
    Next we form FR* and FR to complete: Fr + Fr* = 0.
    We have the equations of motion at this point.
    It makes sense to rearrnge them though, so we calculate the mass matrix and
    the forcing terms, for E.o.M. in the form: [MM] udot = forcing, where MM is
    the mass matrix, udot is a vector of the time derivatives of the
    generalized speeds, and forcing is a vector representing "forcing" terms.

        >>> KM = KanesMethod(N, q_ind=[q], u_ind=[u], kd_eqs=kd)
        >>> (fr, frstar) = KM.kanes_equations(BL, FL)
        >>> MM = KM.mass_matrix
        >>> forcing = KM.forcing
        >>> rhs = MM.inv() * forcing
        >>> rhs
        Matrix([[(-c*u(t) - k*q(t))/m]])
        >>> KM.linearize(A_and_B=True, new_method=True)[0]
        Matrix([
        [   0,    1],
        [-k/m, -c/m]])

    Please look at the documentation pages for more information on how to
    perform linearization and how to deal with dependent coordinates & speeds,
    and how do deal with bringing non-contributing forces into evidence.

    References
    ==========

    .. [Kane1985] Kane, T., Levinson, D. Dynamics Theory and Applications. 1985
       McGraw-Hill

    """

    def __init__(self, frame, q_ind, u_ind, kd_eqs=None, q_dependent=None,
                 configuration_constraints=None, u_dependent=None,
                 velocity_constraints=None, acceleration_constraints=None,
                 u_auxiliary=None):
        """Initializes a KanesMethod object.

        Parameters
        ==========

        frame : ReferenceFrame
            The inertial reference frame with which to form the equations of
            motion with respect to. All bodies and particles must have their
            velocity defined with respect to this frame, either directly or
            indirectly.
        q_ind : ordered iterable of functions of time, len(n)
            The n independent generalized coordinates.
        u_ind : ordered iterable of functions of time, len(p)
            The p independent generalized speeds.
        kd_eqs : iterable of SymPy expressions, optional
            The left hand side of the zero equal kinematical differential
            equations. (TODO: check if following statement is true) If None,
            the simple definition of u = q' is used.
        q_dependent : ordered iterable functions of time, len(o), optional
            The o dependent generalized coordinates.
        configuration_constraints : iterable of expressions, optional
            The left hand side of the zero equal holonomic constraints.
        u_dependent : ordered iterable, len(m), optional
            The m dependent generalized speeds.
        velocity_constraints : iterable of expressions, optional
            This iterable should contain the expressions that represent the
            left hand side of zero equal nonholomic constraints.
        acceleration_constraints : iterable of expressions, optional
            This iterable should contain the derivatives of the velocity
            constraints. If not provided they will be computed automatically if
            velocity constraints are provided.
        u_auxilary : ordered iterable of functions of time, optional
            The auxiliary speeds that are present in the bodies' velocity
            expressions.

        Notes
        =====

        n : number of independent generalized coordinates
        m : number of dependent generalized speeds and nonholomic constraints
        o : number of dependent generalized coordinates and holonomic constraints
        p : number of independent generalized speeds

        """

        if not isinstance(frame, ReferenceFrame):
            raise TypeError('A ReferenceFrame must be supplied.')
        self._inertial = frame

        self._fr = None
        self._frstar = None

        self._forcelist = None
        self._bodylist = None

        self._initialize_vectors(q_ind, q_dependent, u_ind, u_dependent,
                                 u_auxiliary)

        self._initialize_kindiffeq_matrices(kd_eqs)

        self._initialize_constraint_matrices(configuration_constraints,
                                             velocity_constraints,
                                             acceleration_constraints)

    def _initialize_vectors(self, q_ind, q_dep, u_ind, u_dep, u_aux):
        """Initialize the coordinate and speed vectors and sets the private
        class attributes that contain them. The results are publicly accessible
        through the attributes described in the class docstring."""

        # Initialize generalized coordinates
        q_dep = _none_handler(q_dep)
        if not iterable(q_dep):
            msg = 'Dependent generalized coordinates must be an iterable.'
            raise TypeError(msg)
        if not iterable(q_ind):
            msg = 'Independent generalized coordinates must be an iterable.'
            raise TypeError(msg)
        q_ind = Matrix(q_ind)

        self._qind = q_ind
        self._qdep = q_dep
        self._q = Matrix([q_ind, q_dep])
        self._qdot = self.q.diff(dynamicsymbols._t)

        # Initialize generalized speeds
        u_dep = _none_handler(u_dep)
        if not iterable(u_dep):
            raise TypeError('Generalized dependent speeds must be an iterable.')
        if not iterable(u_ind):
            raise TypeError('Generalized indepenent speeds must be an iterable.')
        u_ind = Matrix(u_ind)

        self._uind = u_ind
        self._udep = u_dep
        self._u = Matrix([u_ind, u_dep])
        self._udot = self.u.diff(dynamicsymbols._t)
        self._uaux = _none_handler(u_aux)

        self._n = len(self._u)
        self._m = len(u_dep)
        self._o = len(q_dep)
        self._p = len(u_ind)

    def _initialize_constraint_matrices(self, config, vel, acc):
        """Initializes constraint matrices.

        Given::

            config = f_h(q, t) = 0
            vel = K_nh(q, t) * u + f_nh(q, t) = 0
            acc = K_dnh(q, t) * u' + f_dnh(q, u, t) = 0

        this function finds and stores::

            f_h, K_nh, f_nh, K_dnh, f_dnh, and A_rs

        where::

            u = [u_s, u_r]^T
            u_r = A_rs * u_s + B_r

        Parameters
        ==========

        config : iterable of expressions
            The configuration constraint expressions.
        vel : iterable of expressions
            The velocity constraint expressions.
        acc : iterable of expressions
            The acceleration constraint expressions.

        """

        # Initialize configuration constraints
        config = _none_handler(config)
        if len(config) != self.num_dep_coords:
            raise ValueError('There must be an equal number of dependent '
                             'coordinates and configuration constraints.')
        self._f_h = _none_handler(config)

        # Initialize velocity and acceleration constraints
        vel = _none_handler(vel)
        if len(vel) != self.num_dependent_speeds:
            raise ValueError('There must be an equal number of dependent '
                             'speeds and velocity constraints.')
        acc = _none_handler(acc)
        if acc and (len(acc) != self.num_dependent_speeds):
            raise ValueError('There must be an equal number of dependent '
                             'speeds and acceleration constraints.')
        if vel:
            u_zero = dict((i, 0) for i in self.u)
            udot_zero = dict((i, 0) for i in self._udot)

            # When calling kanes_equations, another class instance will be
            # created if auxiliary u's are present. In this case, the
            # computation of kinetic differential equation matrices will be
            # skipped as this was computed during the original KanesMethod
            # object, and the qd_u_map will not be available.
            if self._qdot_u_map is not None:
                vel = vel.xreplace(self._qdot_u_map)

            # mxn * nx1 + mx1 = 0
            # K_nh(q, t) * u + f_nh(q, t) = 0
            self._f_nh = vel.xreplace(u_zero)
            self._k_nh = vel.jacobian(self.u)

            # If no acceleration constraints given, calculate them.
            # d[K_nh(q, t) * u + f_nh(q, t)] / dt = 0
            # K_nh(q, t) * u' + K_nh'(q, t) * u + f_nh'(q, t) = 0
            # K_dnh(q, t) * u' + f_dnh(q, q', u, t) = 0
            if not acc:
                self._f_dnh = (self._k_nh.diff(dynamicsymbols._t) * self.u +
                               self._f_nh.diff(dynamicsymbols._t))
                self._k_dnh = self._k_nh
            else:
                if self._qdot_u_map is not None:
                    acc = acc.xreplace(self._qdot_u_map)
                self._f_dnh = acc.xreplace(udot_zero)
                self._k_dnh = acc.jacobian(self._udot)

            # The dependent speeds can be defined in terms of the independent
            # speeds using the nonholomic constraints. First, the
            # non-holonmic matrix needs to be partitioned wrt to the speeds:
            #
            # K_nh(q, t) = |B_s|
            #              |B_r|
            #
            # so the nonholomic constraints can be written as:
            #
            # B_r * u_r + B_s * u_s + f_nh(q, t) = 0
            #
            # The dependent speeds are solved for as such:
            #
            # B_r * u_r = -B_s * u_s - f_nh(q, t)
            # u_r = -B_r^-1 * B_s * u_s - B_r^-1 * f_nh(q, t)
            #
            # Kane defines the linear coefficient matrix that relates the
            # dependent speeds to the independent speeds, A_rs (mxp), as:
            #
            # u_r = A_rs * u_s + B_r (Kane Page 43, equation 1)
            #
            # So,
            #
            # A_rs = -B_r^-1 * B_s
            # B_r = - B_r^-1 * f_nh(q, t)
            #
            # with these sizes:
            #
            # mxm * mxp  =  mxp
            # B_r * A_rs = -B_s

            B_s = self._k_nh[:, :self.num_independent_speeds]
            B_r = self._k_nh[:, self.num_independent_speeds:]
            self._A_rs = -B_r.LUsolve(B_s)

        else:  # this is a holonomic system
            self._f_nh = Matrix()
            self._k_nh = Matrix()
            self._f_dnh = Matrix()
            self._k_dnh = Matrix()
            self._A_rs = Matrix()

    def _initialize_kindiffeq_matrices(self, kdeqs):
        """Initialize the kinematic differential equation matrices.

        The kinematic differential equations are in this form::

            K_ku(q, t) * u + K_kqdot(q, t) * q' + f_k(q, t) = 0

        This method finds and stores K_ku, K_kqdot, and f_k. It also solves the
        equations for q' and stores the solution.

        """

        if kdeqs:
            if len(kdeqs) != self.num_coordinates:
                msg = ('There must be an equal number of kinematic '
                       'differential equations and generalized coordinates.')
                raise ValueError(msg)

            kdeqs = Matrix(kdeqs)

            # Substitution dictionaries for setting variables to zero.
            u_zero = dict((i, 0) for i in self.u)
            uaux_zero = dict((i, 0) for i in self._uaux)
            qdot_zero = dict((i, 0) for i in self._qdot)

            # K_ku(q, t) * u + K_kq(q, t) * q' + f_k(q, t) = 0
            f_k = kdeqs.xreplace(u_zero).xreplace(qdot_zero)
            k_ku = kdeqs.jacobian(self.u)
            k_kqdot = kdeqs.jacobian(self._qdot)

            # K_kq(q, t) * q' = -K_ku(q, t) * u - f_k(q, t)
            sol = k_kqdot.LUsolve(-k_ku * self.u - f_k)
            self._qdot_u_map = dict(zip(self._qdot, sol))

            # NOTE : Not sure why this is necessary.
            self._f_k = f_k.xreplace(uaux_zero)
            self._k_ku = k_ku.xreplace(uaux_zero)
            self._k_kqdot = k_kqdot
        else:
            self._qdot_u_map = None
            self._f_k = Matrix()
            self._k_ku = Matrix()
            self._k_kqdot = Matrix()

    def _form_fr(self, forces):
        """Returns the holonomic or nonholonomic generalized active forces."""

        if forces is not None and (len(forces) == 0 or not iterable(forces)):
            raise ValueError('Force pairs must be supplied in an non-empty '
                             'iterable or it should be None.')

        # Pull out relevant velocities for constructing partial velocities.
        vel_list, f_list = _f_list_parser(forces, self._inertial)
        # Make sure all velocities are strictly a function of the speeds.
        vel_list = [msubs(v, self._qdot_u_map) for v in vel_list]

        # Fill Fr with dot product of partial velocities and forces.
        nu = len(f_list)
        Fr = zeros(self.n, 1)
        partials = partial_velocity(vel_list, self.u, self._inertial)
        for i in range(self.n):
            Fr[i] = sum(partials[j][i].dot(f_list[j]) for j in range(nu))

        # In case there are dependent speeds
        # Kane Page 99, Equation 3: F_r_tilde = F_r + F_s * A_sr
        if self.nonholonomic:
            Frtilde = Fr[:self.p, 0]
            Frold = Fr[self.p:self.n, 0]
            Frtilde += self._A_rs.T * Frold
            Fr = Frtilde

        self._forcelist = forces
        self._fr = Fr
        return Fr

    def _form_frstar(self, bl):
        """Returns the holonomic or nonholonomic generalized inertia forces."""

        if not iterable(bl):
            raise TypeError('Bodies must be supplied in an iterable.')

        t = dynamicsymbols._t
        N = self._inertial
        # Dicts setting things to zero
        udot_zero = dict((i, 0) for i in self._udot)
        uaux_zero = dict((i, 0) for i in self._uaux)
        uauxdot = [diff(i, t) for i in self._uaux]
        uauxdot_zero = dict((i, 0) for i in uauxdot)
        # Dictionary of q' and q'' to u and u'
        q_ddot_u_map = dict((k.diff(t), v.diff(t)) for (k, v) in
                self._qdot_u_map.items())
        q_ddot_u_map.update(self._qdot_u_map)

        # Fill up the list of partials: format is a list with num elements
        # equal to number of entries in body list. Each of these elements is a
        # list - either of length 1 for the translational components of
        # particles or of length 2 for the translational and rotational
        # components of rigid bodies. The inner most list is the list of
        # partial velocities.
        def get_partial_velocity(body):
            if isinstance(body, RigidBody):
                vlist = [body.masscenter.vel(N), body.frame.ang_vel_in(N)]
            elif isinstance(body, Particle):
                vlist = [body.point.vel(N),]
            else:
                raise TypeError('The body list may only contain either '
                                'RigidBody or Particle as list elements.')
            v = [msubs(vel, self._qdot_u_map) for vel in vlist]
            return partial_velocity(v, self.u, N)
        partials = [get_partial_velocity(body) for body in bl]

        # Compute fr_star in two components:
        # fr_star = -(MM*u' + nonMM)
        o = len(self.u)
        MM = zeros(o, o)
        nonMM = zeros(o, 1)
        zero_uaux = lambda expr: msubs(expr, uaux_zero)
        zero_udot_uaux = lambda expr: msubs(msubs(expr, udot_zero), uaux_zero)
        for i, body in enumerate(bl):
            if isinstance(body, RigidBody):
                M = zero_uaux(body.mass)
                I = zero_uaux(body.central_inertia)
                vel = zero_uaux(body.masscenter.vel(N))
                omega = zero_uaux(body.frame.ang_vel_in(N))
                acc = zero_udot_uaux(body.masscenter.acc(N))
                inertial_force = (M.diff(t) * vel + M * acc)
                inertial_torque = zero_uaux((I.dt(body.frame) & omega) +
                    msubs(I & body.frame.ang_acc_in(N), udot_zero) +
                    (omega ^ (I & omega)))
                for j in range(o):
                    tmp_vel = zero_uaux(partials[i][0][j])
                    tmp_ang = zero_uaux(I & partials[i][1][j])
                    for k in range(o):
                        # translational
                        MM[j, k] += M * (tmp_vel & partials[i][0][k])
                        # rotational
                        MM[j, k] += (tmp_ang & partials[i][1][k])
                    nonMM[j] += inertial_force & partials[i][0][j]
                    nonMM[j] += inertial_torque & partials[i][1][j]
            else:
                M = zero_uaux(body.mass)
                vel = zero_uaux(body.point.vel(N))
                acc = zero_udot_uaux(body.point.acc(N))
                inertial_force = (M.diff(t) * vel + M * acc)
                for j in range(o):
                    temp = zero_uaux(partials[i][0][j])
                    for k in range(o):
                        MM[j, k] += M * (temp & partials[i][0][k])
                    nonMM[j] += inertial_force & partials[i][0][j]
        # Compose fr_star out of MM and nonMM
        MM = zero_uaux(msubs(MM, q_ddot_u_map))
        nonMM = msubs(msubs(nonMM, q_ddot_u_map), udot_zero, uauxdot_zero,
                      uaux_zero)
        fr_star = -(MM * msubs(Matrix(self._udot), uauxdot_zero) + nonMM)

        # If there are dependent speeds, we need to find fr_star_tilde
        if self._udep:
            p = o - len(self._udep)
            fr_star_ind = fr_star[:p, 0]
            fr_star_dep = fr_star[p:o, 0]
            fr_star = fr_star_ind + (self._A_rs.T * fr_star_dep)
            # Apply the same to MM
            MMi = MM[:p, :]
            MMd = MM[p:o, :]
            MM = MMi + (self._A_rs.T * MMd)

        self._bodylist = bl
        self._frstar = fr_star
        self._k_d = MM
        self._f_d = -msubs(self._fr + self._frstar, udot_zero)
        return fr_star

    def to_linearizer(self):
        """Returns an instance of the Linearizer class, initiated from the
        data in the KanesMethod class. This may be more desirable than using
        the linearize class method, as the Linearizer object will allow more
        efficient recalculation (i.e. about varying operating points)."""

        if (self._fr is None) or (self._frstar is None):
            raise ValueError('Need to compute Fr, Fr* first.')

        # Get required equation components. The Kane's method class breaks
        # these into pieces. Need to reassemble
        f_c = self._f_h
        if self._f_nh and self._k_nh:
            f_v = self._f_nh + self._k_nh*Matrix(self.u)
        else:
            f_v = Matrix()
        if self._f_dnh and self._k_dnh:
            f_a = self._f_dnh + self._k_dnh*Matrix(self._udot)
        else:
            f_a = Matrix()
        # Dicts to sub to zero, for splitting up expressions
        u_zero = dict((i, 0) for i in self.u)
        ud_zero = dict((i, 0) for i in self._udot)
        qd_zero = dict((i, 0) for i in self._qdot)
        qd_u_zero = dict((i, 0) for i in Matrix([self._qdot, self.u]))
        # Break the kinematic differential eqs apart into f_0 and f_1
        f_0 = msubs(self._f_k, u_zero) + self._k_kqdot*Matrix(self._qdot)
        f_1 = msubs(self._f_k, qd_zero) + self._k_ku*Matrix(self.u)
        # Break the dynamic differential eqs into f_2 and f_3
        f_2 = msubs(self._frstar, qd_u_zero)
        f_3 = msubs(self._frstar, ud_zero) + self._fr
        f_4 = zeros(len(f_2), 1)

        # Get the required vector components
        q = self.q
        u = self.u
        if self._qdep:
            q_i = q[:-len(self._qdep)]
        else:
            q_i = q
        q_d = self._qdep
        if self._udep:
            u_i = u[:-len(self._udep)]
        else:
            u_i = u
        u_d = self._udep

        # Form dictionary to set auxiliary speeds & their derivatives to 0.
        uaux = self._uaux
        uauxdot = uaux.diff(dynamicsymbols._t)
        uaux_zero = dict((i, 0) for i in Matrix([uaux, uauxdot]))

        # Checking for dynamic symbols outside the dynamic differential
        # equations; throws error if there is.
        sym_list = set(Matrix([q, self._qdot, u, self._udot, uaux, uauxdot]))
        if any(find_dynamicsymbols(i, sym_list) for i in [self._k_kqdot,
                self._k_ku, self._f_k, self._k_dnh, self._f_dnh, self._k_d]):
            raise ValueError('Cannot have dynamicsymbols outside dynamic \
                             forcing vector.')

        # Find all other dynamic symbols, forming the forcing vector r.
        # Sort r to make it canonical.
        r = list(find_dynamicsymbols(msubs(self._f_d, uaux_zero), sym_list))
        r.sort(key=default_sort_key)

        # Check for any derivatives of variables in r that are also found in r.
        for i in r:
            if diff(i, dynamicsymbols._t) in r:
                raise ValueError('Cannot have derivatives of specified \
                                 quantities when linearizing forcing terms.')
        return Linearizer(f_0, f_1, f_2, f_3, f_4, f_c, f_v, f_a, q, u, q_i,
                          q_d, u_i, u_d, r)

    def linearize(self, **kwargs):
        """ Linearize the equations of motion about a symbolic operating point.

        If kwarg A_and_B is False (default), returns M, A, B, r for the
        linearized form, M*[q', u']^T = A*[q_ind, u_ind]^T + B*r.

        If kwarg A_and_B is True, returns A, B, r for the linearized form
        dx = A*x + B*r, where x = [q_ind, u_ind]^T. Note that this is
        computationally intensive if there are many symbolic parameters. For
        this reason, it may be more desirable to use the default A_and_B=False,
        returning M, A, and B. Values may then be substituted in to these
        matrices, and the state space form found as
        A = P.T*M.inv()*A, B = P.T*M.inv()*B, where P = Linearizer.perm_mat.

        In both cases, r is found as all dynamicsymbols in the equations of
        motion that are not part of q, u, q', or u'. They are sorted in
        canonical form.

        The operating points may be also entered using the ``op_point`` kwarg.
        This takes a dictionary of {symbol: value}, or a an iterable of such
        dictionaries. The values may be numberic or symbolic. The more values
        you can specify beforehand, the faster this computation will run.

        As part of the deprecation cycle, the new method will not be used unless
        the kwarg ``new_method`` is set to True. If the kwarg is missing, or set
        to false, the old linearization method will be used. After next release
        the need for this kwarg will be removed.

        For more documentation, please see the ``Linearizer`` class."""

        if 'new_method' not in kwargs or not kwargs['new_method']:
            # User is still using old code.
            SymPyDeprecationWarning('The linearize class method has changed '
                    'to a new interface, the old method is deprecated. To '
                    'use the new method, set the kwarg `new_method=True`. '
                    'For more information, read the docstring '
                    'of `linearize`.').warn()
            return self._old_linearize()
        # Remove the new method flag, before passing kwargs to linearize
        kwargs.pop('new_method')
        linearizer = self.to_linearizer()
        result = linearizer.linearize(**kwargs)
        return result + (linearizer.r,)

    def _old_linearize(self):
        """Old method to linearize the equations of motion. Returns a tuple of
        (f_lin_A, f_lin_B, y) for forming [M]qudot = [f_lin_A]qu + [f_lin_B]y.

        Deprecated in favor of new method using Linearizer class. Please change
        your code to use the new `linearize` method."""

        if (self._fr is None) or (self._frstar is None):
            raise ValueError('Need to compute Fr, Fr* first.')

        # Note that this is now unneccessary, and it should never be
        # encountered; I still think it should be in here in case the user
        # manually sets these matrices incorrectly.
        for i in self.q:
            if self._k_kqdot.diff(i) != 0 * self._k_kqdot:
                raise ValueError('Matrix K_kqdot must not depend on any q.')

        t = dynamicsymbols._t
        uaux = self._uaux
        uauxdot = [diff(i, t) for i in uaux]
        # dictionary of auxiliary speeds & derivatives which are equal to zero
        subdict = dict(zip(uaux[:] + uauxdot[:],
                           [0] * (len(uaux) + len(uauxdot))))

        # Checking for dynamic symbols outside the dynamic differential
        # equations; throws error if there is.
        insyms = set(self.q[:] + self._qdot[:] + self.u[:] + self._udot[:] +
                     uaux[:] + uauxdot)
        if any(find_dynamicsymbols(i, insyms) for i in [self._k_kqdot,
                self._k_ku, self._f_k, self._k_dnh, self._f_dnh, self._k_d]):
            raise ValueError('Cannot have dynamicsymbols outside dynamic \
                             forcing vector.')
        other_dyns = list(find_dynamicsymbols(msubs(self._f_d, subdict), insyms))

        # make it canonically ordered so the jacobian is canonical
        other_dyns.sort(key=default_sort_key)

        for i in other_dyns:
            if diff(i, dynamicsymbols._t) in other_dyns:
                raise ValueError('Cannot have derivatives of specified '
                                 'quantities when linearizing forcing terms.')

        o = len(self.u)  # number of speeds
        n = len(self.q)  # number of coordinates
        l = len(self._qdep)  # number of configuration constraints
        m = len(self._udep)  # number of motion constraints
        qi = Matrix(self.q[: n - l])  # independent coords
        qd = Matrix(self.q[n - l: n])  # dependent coords; could be empty
        ui = Matrix(self.u[: o - m])  # independent speeds
        ud = Matrix(self.u[o - m: o])  # dependent speeds; could be empty
        qdot = Matrix(self._qdot)  # time derivatives of coordinates

        # with equations in the form MM udot = forcing, expand that to:
        # MM_full [q,u].T = forcing_full. This combines coordinates and
        # speeds together for the linearization, which is necessary for the
        # linearization process, due to dependent coordinates. f1 is the rows
        # from the kinematic differential equations, f2 is the rows from the
        # dynamic differential equations (and differentiated nonholomic
        # constraints).
        f1 = self._k_ku * Matrix(self.u) + self._f_k
        f2 = self._f_d
        # Only want to do this if these matrices have been filled in, which
        # occurs when there are dependent speeds
        if m != 0:
            f2 = self._f_d.col_join(self._f_dnh)
            fnh = self._f_nh + self._k_nh * Matrix(self.u)
        f1 = msubs(f1, subdict)
        f2 = msubs(f2, subdict)
        fh = msubs(self._f_h, subdict)
        fku = msubs(self._k_ku * Matrix(self.u), subdict)
        fkf = msubs(self._f_k, subdict)

        # In the code below, we are applying the chain rule by hand on these
        # things. All the matrices have been changed into vectors (by
        # multiplying the dynamic symbols which it is paired with), so we can
        # take the jacobian of them. The basic operation is take the jacobian
        # of the f1, f2 vectors wrt all of the q's and u's. f1 is a function of
        # q, u, and t; f2 is a function of q, qdot, u, and t. In the code
        # below, we are not considering perturbations in t. So if f1 is a
        # function of the q's, u's but some of the q's or u's could be
        # dependent on other q's or u's (qd's might be dependent on qi's, ud's
        # might be dependent on ui's or qi's), so what we do is take the
        # jacobian of the f1 term wrt qi's and qd's, the jacobian wrt the qd's
        # gets multiplied by the jacobian of qd wrt qi, this is extended for
        # the ud's as well. dqd_dqi is computed by taking a taylor expansion of
        # the holonomic constraint equations about q*, treating q* - q as dq,
        # separating into dqd (depedent q's) and dqi (independent q's) and the
        # rearranging for dqd/dqi. This is again extended for the speeds.

        # First case: configuration and motion constraints
        if (l != 0) and (m != 0):
            fh_jac_qi = fh.jacobian(qi)
            fh_jac_qd = fh.jacobian(qd)
            fnh_jac_qi = fnh.jacobian(qi)
            fnh_jac_qd = fnh.jacobian(qd)
            fnh_jac_ui = fnh.jacobian(ui)
            fnh_jac_ud = fnh.jacobian(ud)
            fku_jac_qi = fku.jacobian(qi)
            fku_jac_qd = fku.jacobian(qd)
            fku_jac_ui = fku.jacobian(ui)
            fku_jac_ud = fku.jacobian(ud)
            fkf_jac_qi = fkf.jacobian(qi)
            fkf_jac_qd = fkf.jacobian(qd)
            f1_jac_qi = f1.jacobian(qi)
            f1_jac_qd = f1.jacobian(qd)
            f1_jac_ui = f1.jacobian(ui)
            f1_jac_ud = f1.jacobian(ud)
            f2_jac_qi = f2.jacobian(qi)
            f2_jac_qd = f2.jacobian(qd)
            f2_jac_ui = f2.jacobian(ui)
            f2_jac_ud = f2.jacobian(ud)
            f2_jac_qdot = f2.jacobian(qdot)

            dqd_dqi = - fh_jac_qd.LUsolve(fh_jac_qi)
            dud_dqi = fnh_jac_ud.LUsolve(fnh_jac_qd * dqd_dqi - fnh_jac_qi)
            dud_dui = - fnh_jac_ud.LUsolve(fnh_jac_ui)
            dqdot_dui = - self._k_kqdot.inv() * (fku_jac_ui +
                                                fku_jac_ud * dud_dui)
            dqdot_dqi = - self._k_kqdot.inv() * (fku_jac_qi + fkf_jac_qi +
                    (fku_jac_qd + fkf_jac_qd) * dqd_dqi + fku_jac_ud * dud_dqi)
            f1_q = f1_jac_qi + f1_jac_qd * dqd_dqi + f1_jac_ud * dud_dqi
            f1_u = f1_jac_ui + f1_jac_ud * dud_dui
            f2_q = (f2_jac_qi + f2_jac_qd * dqd_dqi + f2_jac_qdot * dqdot_dqi +
                    f2_jac_ud * dud_dqi)
            f2_u = f2_jac_ui + f2_jac_ud * dud_dui + f2_jac_qdot * dqdot_dui
        # Second case: configuration constraints only
        elif l != 0:
            dqd_dqi = - fh.jacobian(qd).LUsolve(fh.jacobian(qi))
            dqdot_dui = - self._k_kqdot.inv() * fku.jacobian(ui)
            dqdot_dqi = - self._k_kqdot.inv() * (fku.jacobian(qi) +
                fkf.jacobian(qi) + (fku.jacobian(qd) + fkf.jacobian(qd)) *
                dqd_dqi)
            f1_q = (f1.jacobian(qi) + f1.jacobian(qd) * dqd_dqi)
            f1_u = f1.jacobian(ui)
            f2_jac_qdot = f2.jacobian(qdot)
            f2_q = (f2.jacobian(qi) + f2.jacobian(qd) * dqd_dqi +
                    f2.jac_qdot * dqdot_dqi)
            f2_u = f2.jacobian(ui) + f2_jac_qdot * dqdot_dui
        # Third case: motion constraints only
        elif m != 0:
            dud_dqi = fnh.jacobian(ud).LUsolve(- fnh.jacobian(qi))
            dud_dui = - fnh.jacobian(ud).LUsolve(fnh.jacobian(ui))
            dqdot_dui = - self._k_kqdot.inv() * (fku.jacobian(ui) +
                                                fku.jacobian(ud) * dud_dui)
            dqdot_dqi = - self._k_kqdot.inv() * (fku.jacobian(qi) +
                    fkf.jacobian(qi) + fku.jacobian(ud) * dud_dqi)
            f1_jac_ud = f1.jacobian(ud)
            f2_jac_qdot = f2.jacobian(qdot)
            f2_jac_ud = f2.jacobian(ud)
            f1_q = f1.jacobian(qi) + f1_jac_ud * dud_dqi
            f1_u = f1.jacobian(ui) + f1_jac_ud * dud_dui
            f2_q = (f2.jacobian(qi) + f2_jac_qdot * dqdot_dqi + f2_jac_ud
                    * dud_dqi)
            f2_u = (f2.jacobian(ui) + f2_jac_ud * dud_dui + f2_jac_qdot *
                    dqdot_dui)
        # Fourth case: No constraints
        else:
            dqdot_dui = - self._k_kqdot.inv() * fku.jacobian(ui)
            dqdot_dqi = - self._k_kqdot.inv() * (fku.jacobian(qi) +
                    fkf.jacobian(qi))
            f1_q = f1.jacobian(qi)
            f1_u = f1.jacobian(ui)
            f2_jac_qdot = f2.jacobian(qdot)
            f2_q = f2.jacobian(qi) + f2_jac_qdot * dqdot_dqi
            f2_u = f2.jacobian(ui) + f2_jac_qdot * dqdot_dui
        f_lin_A = -(f1_q.row_join(f1_u)).col_join(f2_q.row_join(f2_u))
        if other_dyns:
            f1_oths = f1.jacobian(other_dyns)
            f2_oths = f2.jacobian(other_dyns)
            f_lin_B = -f1_oths.col_join(f2_oths)
        else:
            f_lin_B = Matrix()
        return (f_lin_A, f_lin_B, Matrix(other_dyns))

    def kanes_equations(self, bodies, loads=None):
        """Returns the components of Kane's equations: Fr + Fr* = 0.

        In the case where auxiliary generalized speeds are present (say, s
        auxiliary speeds, n generalized speeds, and m motion constraints) the
        length of the returned vectors will be n - m + s in length. The first n
        - m equations will be the constrained Kane's equations, then the s
        auxiliary Kane's equations. These auxiliary equations can be accessed
        with the auxiliary_eqs().

        Parameters
        ==========

        bodies : iterable of Particle and RigidBody
            An iterable of all RigidBody's and Particle's in the system. A
            system must have at least one particle or body.
        loads : iterable
            Takes in an iterable of (Particle, Vector) or (ReferenceFrame,
            Vector) tuples which represent the force at a point or torque on a
            frame, respectively. Must be either a non-empty iterable of tuples
            or None which corresponds to a system with no loads.

        """
        if (bodies is None and loads is not None) or isinstance(bodies[0], tuple):
            # This switches the order if they use the old way.
            bodies, loads = loads, bodies
            SymPyDeprecationWarning(value='The API for kanes_equations() has changed such '
                    'that the loads (forces and torques) are now the second argument '
                    'and is optional with None being the default.',
                    feature='The kanes_equation() argument order',
                    useinstead='switched argument order to update your code, For example: '
                    'kanes_equations(loads, bodies) > kanes_equations(bodies, loads).',
                    issue=10945, deprecated_since_version="1.1").warn()

        # TODO : Why is it optional to pass in the kinematic differential
        # equations? I'm not sure what the use case is without them.
        if not self._k_kqdot:
            raise AttributeError('Create an instance of KanesMethod with '
                                 'kinematic differential equations to use this '
                                 'method.')

        # TODO : Seems self._fr and self._frstar are set inside the following
        # methods. Make it consistent with setting them in this method below.
        fr = self._form_fr(loads)
        frstar = self._form_frstar(bodies)

        if self._uaux:
            if self.holonomic:
                km = KanesMethod(self._inertial, self.q, self._uaux,
                                 u_auxiliary=self._uaux)
            else:
                km = KanesMethod(self._inertial, self.q, self._uaux,
                                 u_auxiliary=self._uaux,
                                 u_dependent=self._udep,
                                 velocity_constraints=(self._k_nh * self.u +
                                                       self._f_nh))
            km._qdot_u_map = self._qdot_u_map
            self._km = km
            fraux = km._form_fr(loads)
            frstaraux = km._form_frstar(bodies)
            self._aux_eq = fraux + frstaraux
            self._fr = fr.col_join(fraux)
            self._frstar = frstar.col_join(frstaraux)

        return (self._fr, self._frstar)

    def rhs(self, inv_method=None):
        """Returns the system's equations of motion in first order form.

        The output of this will be the right hand side of:

        [qdot, udot].T = f(q, u, t)

        Or, the equations of motion in first order form.  The right hand side
        is what is needed by most numerical ODE integrators.

        Parameters
        ==========

        inv_method : str
            The specific sympy inverse matrix calculation method to use. For a
            list of valid methods, see
            :meth:`~sympy.matrices.matrices.MatrixBase.inv`

        """
        if inv_method is None:
            self._rhs = self.mass_matrix_full.LUsolve(self.forcing_full)
        else:
            self._rhs = (self.mass_matrix_full.inv(inv_method,
                         try_block_diag=True) * self.forcing_full)
        return self._rhs

    def kindiffdict(self):
        """Returns a dictionary mapping q' to functions of u."""
        # TODO : This should be a property not a method since it is already
        # computed.
        if not self._qdot_u_map:
            msg = ('Create an instance of KanesMethod with kinematic '
                   'differential equations to use this method.')
            raise AttributeError(msg)
        return self._qdot_u_map

    @property
    def auxiliary_eqs(self):
        """A matrix containing the auxiliary equations."""
        if not self._fr or not self._frstar:
            raise ValueError('Need to compute Fr, Fr* first.')
        if not self._uaux:
            raise ValueError('No auxiliary speeds have been declared.')
        return self._aux_eq

    @property
    def mass_matrix(self):
        """The mass matrix of the following dynamic differential equations::

            M * u' = f

        where u contains the independent and dependent generalized speeds::

            u = [u_ind]
                [u_dep]

        The equations can be formed with these class attributes::

            >>> lhs = mass_matrix * generalized_accelerations
            >>> rhs = forcing
            >>> Eq(lhs, rhs)

        """
        if not self._fr or not self._frstar:
            msg = ('kanes_equations() must be run before the mass matrix is '
                   'available.')
            raise ValueError(msg)
        return Matrix([self._k_d, self._k_dnh])

    @property
    def mass_matrix_full(self):
        """The mass matrix of the system, augmented by the kinematic
        differential equations."""
        if not self._fr or not self._frstar:
            raise ValueError('Need to compute Fr, Fr* first.')
        o = len(self.u)
        n = len(self.q)
        return ((self._k_kqdot).row_join(zeros(n, o))).col_join((zeros(o,
                n)).row_join(self.mass_matrix))

    @property
    def forcing(self):
        """The right hand side of the following dynamic equations of motion::

            M * u' = f

        where u contains the independent and dependent generalized speeds:

            u = [u_ind]
                [u_dep]

        The equations can be formed with these class attributes::

            >>> lhs = mass_matrix * generalized_accelerations
            >>> rhs = forcing
            >>> Eq(lhs, rhs)

        """
        if not self._fr or not self._frstar:
            raise ValueError('Need to compute Fr, Fr* first.')
        return -Matrix([self._f_d, self._f_dnh])

    @property
    def forcing_full(self):
        """The forcing vector of the system, augmented by the kinematic
        differential equations."""
        if not self._fr or not self._frstar:
            raise ValueError('Need to compute Fr, Fr* first.')
        f1 = self._k_ku * Matrix(self.u) + self._f_k
        return -Matrix([f1, self._f_d, self._f_dnh])

    # Read only attributes are below.

    @property
    def num_coordinates(self):
        """Returns the number of generalized coordinates."""
        return len(self._q)

    @property
    def q(self):
        """Returns a column matrix of the generalized coordinates."""
        return self._q

    @property
    def coordinates(self):
        """Returns a column matrix of the generalized coordinates."""
        return self._q

    @property
    def num_dep_coords(self):
        """Returns the number of dependent generalized coordinates."""
        return self._o

    @property
    def num_ind_coords(self):
        """Returns the number of independent generalized coordinates."""
        return self._n - self._o

    @property
    def n(self):
        """Returns the number of generalized speeds."""
        return self._n

    @property
    def num_speeds(self):
        """Returns the total number of generalized speeds."""
        return self._n

    @property
    def u(self):
        """Returns the column matrix of generalized speeds."""
        return self._u

    @property
    def speeds(self):
        """Returns the column matrix of generalized speeds."""
        return self._u

    @property
    def p(self):
        return self._p

    @property
    def num_independent_speeds(self):
        """Returns the number of independent generalized speeds."""
        return self._p

    @property
    def independent_speeds(self):
        """Returns a matrix of independent generalized speeds."""
        return self._uind

    @property
    def num_dependent_speeds(self):
        """Returns the number of dependent generalized speeds."""
        return self._m

    @property
    def dependent_speeds(self):
        """Returns a matrix of dependent generalized speeds."""
        return self._udep

    @property
    def num_auxiliary_speeds(self):
        """Returns the number of auxiliary speeds."""
        return len(self._uaux)

    @property
    def auxiliary_speeds(self):
        """Returns a matrix of auxiliary generalized speeds."""
        return self._uaux

    @property
    def bodylist(self):
        """Returns a list of the Partical and Body objects."""
        return self._bodylist

    @property
    def bodies(self):
        """Returns a list of the Particle and Body objects."""
        return self._bodylist

    @property
    def forcelist(self):
        """Returns a list of the load tuples."""
        return self._forcelist

    @property
    def loads(self):
        """Returns a list of the load tuples."""
        return self._forcelist

    @property
    def holonomic(self):
        if self.num_dependent_speeds > 0:
            return False
        else:
            return True

    @property
    def nonholonomic(self):
        return not self.holonomic
