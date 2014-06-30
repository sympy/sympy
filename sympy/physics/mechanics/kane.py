from __future__ import print_function, division

__all__ = ['KanesMethod']

from sympy import Symbol, zeros, Matrix, diff, solve_linear_system_LU, eye
from sympy.core.compatibility import reduce
from sympy.utilities import default_sort_key
from sympy.physics.vector import ReferenceFrame, dynamicsymbols, \
     Point, partial_velocity
from sympy.physics.mechanics.particle import Particle
from sympy.physics.mechanics.rigidbody import RigidBody
from sympy.physics.mechanics.functions import inertia_of_point_mass, _mat_inv_mul

class KanesMethod(object):
    """Kane's method object.

    This object is used to do the "book-keeping" as you go through and form
    equations of motion in the way Kane presents in:
    Kane, T., Levinson, D. Dynamics Theory and Applications. 1985 McGraw-Hill

    The attributes are for equations in the form [M] udot = forcing.

    Attributes
    ==========

    auxiliary : Matrix
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

    Examples
    ========

    This is a simple example for a one degree of freedom translational
    spring-mass-damper.

    In this example, we first need to do the kinematics.
    This involves creating generalized speeds and coordinates and their
    derivatives.
    Then we create a point and set its velocity in a frame::

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
    Finally, a list of all bodies and particles needs to be created::

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
    generalized speeds, and forcing is a vector representing "forcing" terms::

        >>> KM = KanesMethod(N, q_ind=[q], u_ind=[u], kd_eqs=kd)
        >>> (fr, frstar) = KM.kanes_equations(FL, BL)
        >>> MM = KM.mass_matrix
        >>> forcing = KM.forcing
        >>> rhs = MM.inv() * forcing
        >>> rhs
        Matrix([[(-c*u(t) - k*q(t))/m]])
        >>> KM.linearize()[0]
        Matrix([
        [ 0,  1],
        [-k, -c]])

    Please look at the documentation pages for more information on how to
    perform linearization and how to deal with dependent coordinates & speeds,
    and how do deal with bringing non-contributing forces into evidence.

    """

    simp = True
    ___KDEqError = AttributeError('Create an instance of KanesMethod with' +
                                  'kinematic differential equations to use' +
                                  'this method.')

    def __init__(self, frame, q_ind, u_ind, kd_eqs=None, q_dependent=[],
            configuration_constraints=[], u_dependent=[],
            velocity_constraints=[], acceleration_constraints=None,
            u_auxiliary=[]):

        """Please read the online documentation. """
        # Big storage things
        if not isinstance(frame, ReferenceFrame):
            raise TypeError('An intertial ReferenceFrame must be supplied')
        self._inertial = frame
        self._forcelist = None
        self._bodylist = None
        self._fr = None
        self._frstar = None
        self._rhs = None
        self._aux_eq = None

        # States
        self._q = None
        self._qdep = []
        self._qdot = None
        self._u = None
        self._udep = []
        self._udot = None
        self._uaux = None

        # Differential Equations Matrices and Map
        self._k_d = None
        self._f_d = None
        self._k_kqdot = None
        self._k_ku = None
        self._f_k = None
        self._qdot_u_map = None

        # Constraint Matrices
        self._f_h = Matrix([])
        self._k_nh = Matrix([])
        self._f_nh = Matrix([])
        self._k_dnh = Matrix([])
        self._f_dnh = Matrix([])

        self._coords(q_ind, q_dependent, configuration_constraints)
        self._speeds(u_ind, u_dependent, velocity_constraints,
                acceleration_constraints, u_auxiliary)
        if kd_eqs is not None:
            self._kindiffeq(kd_eqs)

    def _find_dynamicsymbols(self, inlist, insyms=[]):
        """Finds all non-supplied dynamicsymbols in the expressions."""
        from sympy.core.function import AppliedUndef, Derivative
        t = dynamicsymbols._t
        return reduce(set.union, [set([i]) for j in inlist
            for i in j.atoms(AppliedUndef, Derivative)
            if i.free_symbols == set([t])], set()) - insyms

        temp_f = set().union(*[i.atoms(AppliedUndef) for i in inlist])
        temp_d = set().union(*[i.atoms(Derivative) for i in inlist])
        set_f = set([a for a in temp_f if a.args == (t,)])
        set_d = set([a for a in temp_d if ((a.args[0] in set_f) and all([i == t
                     for i in a.variables]))])
        return list(set.union(set_f, set_d) - set(insyms))

    def _find_othersymbols(self, inlist, insyms=[]):
        """Finds all non-dynamic symbols in the expressions."""
        return list(reduce(set.union, [i.free_symbols for i in inlist]) -
                    set(insyms))


    def _coords(self, qind, qdep=[], coneqs=[]):
        """Supply all the generalized coordinates in a list.

        If some coordinates are dependent, supply them as part of qdep. Their
        dependent nature will only show up in the linearization process though.

        Parameters
        ==========

        qind : list
            A list of independent generalized coords
        qdep : list
            List of dependent coordinates
        coneq : list
            List of expressions which are equal to zero; these are the
            configuration constraint equations
        """

        if not isinstance(qind, (list, tuple)):
            raise TypeError('Generalized coords. must be supplied in a list.')
        self._q = qind + qdep
        self._qdot = [diff(i, dynamicsymbols._t) for i in self._q]

        if not isinstance(qdep, (list, tuple)):
            raise TypeError('Dependent coordinates and constraints must each be '
                            'provided in their own list.')
        if len(qdep) != len(coneqs):
            raise ValueError('There must be an equal number of dependent '
                             'coordinates and constraints.')
        coneqs = Matrix(coneqs)
        self._qdep = qdep
        self._f_h = coneqs

    def _speeds(self, uind, udep=[], coneqs=[], diffconeqs=None, u_auxiliary=[]):
        """Supply all the generalized speeds in a list.

        If there are motion constraints or auxiliary speeds, they are provided
        here as well (as well as motion constraints).

        Parameters
        ==========

        uind : list
            A list of independent generalized speeds
        udep : list
            Optional list of dependent speeds
        coneqs : list
            Optional List of constraint expressions; these are expressions
            which are equal to zero which define a speed (motion) constraint.
        diffconeqs : list
            Optional, calculated automatically otherwise; list of constraint
            equations; again equal to zero, but define an acceleration
            constraint.
        u_auxiliary : list
            An optional list of auxiliary speeds used for brining
            non-contributing forces into evidence

        """

        if not hasattr(uind, '__iter__'):
            raise TypeError('Supply generalized speeds in an iterable.')
        self._u = uind + udep
        self._udot = [diff(i, dynamicsymbols._t) for i in self._u]
        self._uaux = u_auxiliary

        if not hasattr(udep, '__iter__'):
            raise TypeError('Supply dependent speeds in an iterable.')
        if len(udep) != len(coneqs):
            raise ValueError('There must be an equal number of dependent '
                             'speeds and constraints.')
        if diffconeqs is not None:
            if len(udep) != len(diffconeqs):
                raise ValueError('There must be an equal number of dependent '
                                 'speeds and constraints.')
        if len(udep) != 0:
            u = self._u
            uzero = dict(list(zip(u, [0] * len(u))))
            coneqs = Matrix(coneqs)
            udot = self._udot
            udotzero = dict(list(zip(udot, [0] * len(udot))))

            self._udep = udep
            self._f_nh = coneqs.subs(uzero)
            self._k_nh = (coneqs - self._f_nh).jacobian(u)
            # if no differentiated non holonomic constraints were given, calculate
            if diffconeqs is None:
                self._k_dnh = self._k_nh
                self._f_dnh = (self._k_nh.diff(dynamicsymbols._t) * Matrix(u) +
                               self._f_nh.diff(dynamicsymbols._t))
            else:
                self._f_dnh = diffconeqs.subs(udotzero)
                self._k_dnh = (diffconeqs - self._f_dnh).jacobian(udot)

            o = len(u)  # number of generalized speeds
            m = len(udep)  # number of motion constraints
            p = o - m  # number of independent speeds
            # For a reminder, form of non-holonomic constraints is:
            # B u + C = 0
            B = self._k_nh[:, :]
            C = self._f_nh[:, 0]

            # We partition B into indenpendent and dependent columns
            # Ars is then -Bdep.inv() * Bind, and it relates depedent speeds to
            # independent speeds as: udep = Ars uind, neglecting the C term here.
            self._depB = B
            self._depC = C
            mr1 = B[:, :p]
            ml1 = B[:, p:o]
            self._Ars = - _mat_inv_mul(ml1, mr1)

    def _partial_velocity(self, vlist, ulist, frame):
        """Returns the list of partial velocities, replacing qdot's in the
        velocity list if necessary.
        """
        if self._qdot_u_map is None:
            raise ___KDEqError
        v = [vel.subs(self._qdot_u_map) for vel in vlist]
        return partial_velocity(v, ulist, frame)

    def kindiffdict(self):
        """Returns the qdot's in a dictionary. """
        if self._qdot_u_map is None:
            raise ___KDEqError
        return self._qdot_u_map

    def _kindiffeq(self, kdeqs):
        """Supply all the kinematic differential equations in a list.

        They should be in the form [Expr1, Expr2, ...] where Expri is equal to
        zero

        Parameters
        ==========

        kdeqs : list (of Expr)
            The listof kinematic differential equations

        """
        if len(self._q) != len(kdeqs):
            raise ValueError('There must be an equal number of kinematic '
                             'differential equations and coordinates.')

        uaux = self._uaux
        # dictionary of auxiliary speeds which are equal to zero
        uaz = dict(list(zip(uaux, [0] * len(uaux))))

        #kdeqs = Matrix(kdeqs).subs(uaz)
        kdeqs = Matrix(kdeqs)

        qdot = self._qdot
        qdotzero = dict(list(zip(qdot, [0] * len(qdot))))
        u = self._u
        uzero = dict(list(zip(u, [0] * len(u))))

        f_k = kdeqs.subs(uzero).subs(qdotzero)
        k_kqdot = (kdeqs.subs(uzero) - f_k).jacobian(Matrix(qdot))
        k_ku = (kdeqs.subs(qdotzero) - f_k).jacobian(Matrix(u))

        self._k_ku = _mat_inv_mul(k_kqdot, k_ku)
        self._f_k = _mat_inv_mul(k_kqdot, f_k)
        self._k_kqdot = eye(len(qdot))
        self._qdot_u_map = solve_linear_system_LU(Matrix([self._k_kqdot.T,
            -(self._k_ku * Matrix(self._u) + self._f_k).T]).T, self._qdot)

        self._k_ku = _mat_inv_mul(k_kqdot, k_ku).subs(uaz)
        self._f_k = _mat_inv_mul(k_kqdot, f_k).subs(uaz)

    def _form_fr(self, fl):
        """Form the generalized active force.

        Computes the vector of the generalized active force vector.
        Used to compute E.o.M. in the form Fr + Fr* = 0.

        Parameters
        ==========

        fl : list
            Takes in a list of (Point, Vector) or (ReferenceFrame, Vector)
            tuples which represent the force at a point or torque on a frame.

        """

        if not hasattr(fl, '__iter__'):
            raise TypeError('Force pairs must be supplied in an iterable.')

        N = self._inertial
        self._forcelist = fl[:]
        u = self._u
        o = len(u)  # number of gen. speeds
        b = len(fl)  # number of forces

        FR = zeros(o, 1)

        # pull out relevant velocities for constructing partial velocities
        vel_list = []
        f_list = []
        for i in fl:
            if isinstance(i[0], ReferenceFrame):
                vel_list += [i[0].ang_vel_in(N)]
            elif isinstance(i[0], Point):
                vel_list += [i[0].vel(N)]
            else:
                raise TypeError('First entry in pair must be point or frame.')
            f_list += [i[1]]
        partials = self._partial_velocity(vel_list, u, N)

        # Fill Fr with dot product of partial velocities and forces
        for i in range(o):
            for j in range(b):
                FR[i] += partials[j][i] & f_list[j]

        # In case there are dependent speeds
        m = len(self._udep)  # number of dependent speeds
        if m != 0:
            p = o - m
            FRtilde = FR[:p, 0]
            FRold = FR[p:o, 0]
            FRtilde += self._Ars.T * FRold
            FR = FRtilde

        self._fr = FR
        return FR

    def _form_frstar(self, bl):
        """Form the generalized inertia force.

        Computes the vector of the generalized inertia force vector.
        Used to compute E.o.M. in the form Fr + Fr* = 0.

        Parameters
        ==========

        bl : list
            A list of all RigidBody's and Particle's in the system.

        """

        if not hasattr(bl, '__iter__'):
            raise TypeError('Bodies must be supplied in an iterable.')
        t = dynamicsymbols._t
        N = self._inertial
        self._bodylist = bl
        u = self._u  # all speeds
        udep = self._udep  # dependent speeds
        o = len(u)
        m = len(udep)
        p = o - m
        udot = self._udot
        udotzero = dict(list(zip(udot, [0] * o)))
        # auxiliary speeds
        uaux = self._uaux
        uauxdot = [diff(i, t) for i in uaux]
        # dictionary of auxiliary speeds which are equal to zero
        uaz = dict(list(zip(uaux, [0] * len(uaux))))
        uadz = dict(list(zip(uauxdot, [0] * len(uauxdot))))
        # dictionary of qdot's to u's
        qdots = dict(list(zip(list(self._qdot_u_map.keys()),
                         list(self._qdot_u_map.values()))))
        for k, v in list(qdots.items()):
            qdots[k.diff(t)] = v.diff(t)

        MM = zeros(o, o)
        nonMM = zeros(o, 1)
        partials = []
        # Fill up the list of partials: format is a list with no. elements
        # equal to number of entries in body list. Each of these elements is a
        # list - either of length 1 for the translational components of
        # particles or of length 2 for the translational and rotational
        # components of rigid bodies. The inner most list is the list of
        # partial velocities.
        for v in bl:
            if isinstance(v, RigidBody):
                partials += [self._partial_velocity([v.masscenter.vel(N),
                                                     v.frame.ang_vel_in(N)],
                                                    u, N)]
            elif isinstance(v, Particle):
                partials += [self._partial_velocity([v.point.vel(N)], u, N)]
            else:
                raise TypeError('The body list needs RigidBody or '
                                'Particle as list elements.')

        # This section does 2 things - computes the parts of Fr* that are
        # associated with the udots, and the parts that are not associated with
        # the udots. This happens for RigidBody and Particle a little
        # differently, but similar process overall.
        for i, v in enumerate(bl):
            if isinstance(v, RigidBody):
                M = v.mass.subs(uaz).doit()
                vel = v.masscenter.vel(N).subs(uaz).doit()
                acc = v.masscenter.acc(N).subs(udotzero).subs(uaz).doit()
                inertial_force = (M.diff(t) * vel + M * acc)
                omega = v.frame.ang_vel_in(N).subs(uaz).doit()
                I = v.central_inertia.subs(uaz).doit()
                inertial_torque = ((I.dt(v.frame) & omega).subs(uaz).doit() +
                    (I & v.frame.ang_acc_in(N)).subs(udotzero).subs(uaz).doit() +
                    (omega ^ (I & omega)).subs(uaz).doit())
                for j in range(o):
                    tmp_vel = partials[i][0][j].subs(uaz).doit()
                    tmp_ang = (I & partials[i][1][j].subs(uaz).doit())
                    for k in range(o):
                        # translational
                        MM[j, k] += M * (tmp_vel & partials[i][0][k])
                        # rotational
                        MM[j, k] += (tmp_ang & partials[i][1][k])
                    nonMM[j] += inertial_force & partials[i][0][j]
                    nonMM[j] += inertial_torque & partials[i][1][j]

            if isinstance(v, Particle):
                M = v.mass.subs(uaz).doit()
                vel = v.point.vel(N).subs(uaz).doit()
                acc = v.point.acc(N).subs(udotzero).subs(uaz).doit()
                inertial_force = (M.diff(t) * vel + M * acc)
                for j in range(o):
                    temp = partials[i][0][j].subs(uaz).doit()
                    for k in range(o):
                        MM[j, k] += M * (temp & partials[i][0][k])
                    nonMM[j] += inertial_force & partials[i][0][j]
        # Negate FRSTAR since Kane defines the inertia forces/torques
        # to be negative and we didn't do so above.
        MM = MM.subs(qdots).subs(uaz).doit()
        nonMM = nonMM.subs(qdots).subs(udotzero).subs(uadz).subs(uaz).doit()
        FRSTAR = -(MM * Matrix(udot).subs(uadz) + nonMM)

        # For motion constraints, m is the number of constraints
        # Really, one should just look at Kane's book for descriptions of this
        # process
        if m != 0:
            FRSTARtilde = FRSTAR[:p, 0]
            FRSTARold = FRSTAR[p:o, 0]
            FRSTARtilde += self._Ars.T * FRSTARold
            FRSTAR = FRSTARtilde

            MMi = MM[:p, :]
            MMd = MM[p:o, :]
            MM = MMi + self._Ars.T * MMd
        self._frstar = FRSTAR

        zeroeq = -(self._fr + self._frstar)
        zeroeq = zeroeq.subs(udotzero)

        self._k_d = MM
        self._f_d = zeroeq
        return FRSTAR

    def kanes_equations(self, FL, BL):
        """ Method to form Kane's equations, Fr + Fr* = 0.

        Returns (Fr, Fr*). In the case where auxiliary generalized speeds are
        present (say, s auxiliary speeds, o generalized speeds, and m motion
        constraints) the length of the returned vectors will be o - m + s in
        length. The first o - m equations will be the constrained Kane's
        equations, then the s auxiliary Kane's equations. These auxiliary
        equations can be accessed with the auxiliary_eqs().

        Parameters
        ==========

        FL : list
            Takes in a list of (Point, Vector) or (ReferenceFrame, Vector)
            tuples which represent the force at a point or torque on a frame.
        BL : list
            A list of all RigidBody's and Particle's in the system.

        """

        if (self._q is None) or (self._u is None):
            raise ValueError('Speeds and coordinates must be supplied first.')
        if (self._k_kqdot is None):
            raise __KDEqError


        fr = self._form_fr(FL)
        frstar = self._form_frstar(BL)
        if self._uaux != []:
            if self._udep == []:
                km = KanesMethod(self._inertial, self._q, self._uaux,
                             u_auxiliary=self._uaux)
            else:
                km = KanesMethod(self._inertial, self._q, self._uaux,
                u_auxiliary=self._uaux, u_dependent=self._udep,
                velocity_constraints=(self._k_nh * Matrix(self._u) + self._f_nh))
            km._qdot_u_map = self._qdot_u_map
            self._km = km
            fraux = km._form_fr(FL)
            frstaraux = km._form_frstar(BL)
            self._aux_eq = fraux + frstaraux
            self._fr = fr.col_join(fraux)
            self._frstar = frstar.col_join(frstaraux)
            return (self._fr, self._frstar)
        else:
            return (fr, frstar)

    def linearize(self):
        """ Method used to generate linearized equations.

        Note that for linearization, it is assumed that time is not perturbed,
        but only coordinates and positions. The "forcing" vector's jacobian is
        computed with respect to the state vector in the form [Qi, Qd, Ui, Ud].
        This is the "f_lin_A" matrix.

        It also finds any non-state dynamicsymbols and computes the jacobian of
        the "forcing" vector with respect to them. This is the "f_lin_B"
        matrix; if this is empty, an empty matrix is created.

        Consider the following:
        If our equations are: [M]qudot = f, where [M] is the full mass matrix,
        qudot is a vector of the deriatives of the coordinates and speeds, and
        f in the full forcing vector, the linearization process is as follows:
        [M]qudot = [f_lin_A]qu + [f_lin_B]y, where qu is the state vector,
        f_lin_A is the jacobian of the full forcing vector with respect to the
        state vector, f_lin_B is the jacobian of the full forcing vector with
        respect to any non-speed/coordinate dynamicsymbols which show up in the
        full forcing vector, and y is a vector of those dynamic symbols (each
        column in f_lin_B corresponds to a row of the y vector, each of which
        is a non-speed/coordinate dynamicsymbol).

        To get the traditional state-space A and B matrix, you need to multiply
        the f_lin_A and f_lin_B matrices by the inverse of the mass matrix.
        Caution needs to be taken when inverting large symbolic matrices;
        substituting in numerical values before inverting will work better.

        A tuple of (f_lin_A, f_lin_B, other_dynamicsymbols) is returned.

        """

        if (self._fr is None) or (self._frstar is None):
            raise ValueError('Need to compute Fr, Fr* first.')

        # Note that this is now unneccessary, and it should never be
        # encountered; I still think it should be in here in case the user
        # manually sets these matrices incorrectly.
        for i in self._q:
            if self._k_kqdot.diff(i) != 0 * self._k_kqdot:
                raise ValueError('Matrix K_kqdot must not depend on any q.')

        t = dynamicsymbols._t
        uaux = self._uaux
        uauxdot = [diff(i, t) for i in uaux]
        # dictionary of auxiliary speeds & derivatives which are equal to zero
        subdict = dict(list(zip(uaux + uauxdot, [0] * (len(uaux) + len(uauxdot)))))

        # Checking for dynamic symbols outside the dynamic differential
        # equations; throws error if there is.
        insyms = set(
            self._q + self._qdot + self._u + self._udot + uaux + uauxdot)
        if any(self._find_dynamicsymbols(i, insyms) for i in [self._k_kqdot,
                                                              self._k_ku,
                                                              self._f_k,
                                                              self._k_dnh,
                                                              self._f_dnh,
                                                              self._k_d]):
            raise ValueError('Cannot have dynamicsymbols outside dynamic '
                             'forcing vector.')
        other_dyns = list(self._find_dynamicsymbols(self._f_d.subs(subdict),
                                             insyms))

        # make it canonically ordered so the jacobian is canonical
        other_dyns.sort(key=default_sort_key)

        for i in other_dyns:
            if diff(i, dynamicsymbols._t) in other_dyns:
                raise ValueError('Cannot have derivatives of specified '
                                 'quantities when linearizing forcing terms.')

        o = len(self._u)  # number of speeds
        n = len(self._q)  # number of coordinates
        l = len(self._qdep)  # number of configuration constraints
        m = len(self._udep)  # number of motion constraints
        qi = Matrix(self._q[: n - l])  # independent coords
        qd = Matrix(self._q[n - l: n])  # dependent coords; could be empty
        ui = Matrix(self._u[: o - m])  # independent speeds
        ud = Matrix(self._u[o - m: o])  # dependent speeds; could be empty
        qdot = Matrix(self._qdot)  # time derivatives of coordinates

        # with equations in the form MM udot = forcing, expand that to:
        # MM_full [q,u].T = forcing_full. This combines coordinates and
        # speeds together for the linearization, which is necessary for the
        # linearization process, due to dependent coordinates. f1 is the rows
        # from the kinematic differential equations, f2 is the rows from the
        # dynamic differential equations (and differentiated non-holonomic
        # constraints).
        f1 = self._k_ku * Matrix(self._u) + self._f_k
        f2 = self._f_d
        # Only want to do this if these matrices have been filled in, which
        # occurs when there are dependent speeds
        if m != 0:
            f2 = self._f_d.col_join(self._f_dnh)
            fnh = self._f_nh + self._k_nh * Matrix(self._u)
        f1 = f1.subs(subdict)
        f2 = f2.subs(subdict)
        fh = self._f_h.subs(subdict)
        fku = (self._k_ku * Matrix(self._u)).subs(subdict)
        fkf = self._f_k.subs(subdict)

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

            dqd_dqi = - _mat_inv_mul(fh_jac_qd, fh_jac_qi)
            dud_dqi = _mat_inv_mul(fnh_jac_ud, (fnh_jac_qd *
                                        dqd_dqi - fnh_jac_qi))
            dud_dui = - _mat_inv_mul(fnh_jac_ud, fnh_jac_ui)
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
            dqd_dqi = - _mat_inv_mul(fh.jacobian(qd), fh.jacobian(qi))
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
            dud_dqi = _mat_inv_mul(fnh.jacobian(ud), - fnh.jacobian(qi))
            dud_dui = - _mat_inv_mul(fnh.jacobian(ud), fnh.jacobian(ui))
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
            f_lin_B = Matrix([])
        return (f_lin_A, f_lin_B, Matrix(other_dyns))

    def rhs(self, inv_method=None):
        """ Returns the system's equations of motion in first order form.

        The output of this will be the right hand side of:

        [qdot, udot].T = f(q, u, t)

        Or, the equations of motion in first order form.  The right hand side
        is what is needed by most numerical ODE integrators.

        Parameters
        ==========

        inv_method : str
            The specific sympy inverse matrix calculation method to use. For a
            list of valid methods, see :py:method:
            `~sympy.matrices.matrices.MatrixBase.inv`

        """
        if inv_method is None:
            self._rhs = _mat_inv_mul(self.mass_matrix_full,
                                          self.forcing_full)
        else:
            self._rhs = (self.mass_matrix_full.inv(inv_method,
                         try_block_diag=True) * self.forcing_full)
        return self._rhs

    @property
    def auxiliary_eqs(self):
        if (self._fr is None) or (self._frstar is None):
            raise ValueError('Need to compute Fr, Fr* first.')
        if self._uaux == []:
            raise ValueError('No auxiliary speeds have been declared.')
        return self._aux_eq

    @property
    def mass_matrix(self):
        # Returns the mass matrix, which is augmented by the differentiated non
        # holonomic equations if necessary
        if (self._frstar is None) & (self._fr is None):
            raise ValueError('Need to compute Fr, Fr* first.')
        return Matrix([self._k_d, self._k_dnh])

    @property
    def mass_matrix_full(self):
        # Returns the mass matrix from above, augmented by kin diff's k_kqdot
        if (self._frstar is None) & (self._fr is None):
            raise ValueError('Need to compute Fr, Fr* first.')
        o = len(self._u)
        n = len(self._q)
        return ((self._k_kqdot).row_join(zeros(n, o))).col_join((zeros(o,
                n)).row_join(self.mass_matrix))

    @property
    def forcing(self):
        # Returns the forcing vector, which is augmented by the differentiated
        # non holonomic equations if necessary
        if (self._frstar is None) & (self._fr is None):
            raise ValueError('Need to compute Fr, Fr* first.')
        return -Matrix([self._f_d, self._f_dnh])

    @property
    def forcing_full(self):
        # Returns the forcing vector, which is augmented by the differentiated
        # non holonomic equations if necessary
        if (self._frstar is None) & (self._fr is None):
            raise ValueError('Need to compute Fr, Fr* first.')
        f1 = self._k_ku * Matrix(self._u) + self._f_k
        return -Matrix([f1, self._f_d, self._f_dnh])
