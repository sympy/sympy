__all__ = ['Kane']

from sympy import Symbol, zeros, Matrix, diff, solve_linear_system_LU, eye
from sympy.physics.mechanics.essential import ReferenceFrame, dynamicsymbols
from sympy.physics.mechanics.particle import Particle
from sympy.physics.mechanics.point import Point
from sympy.physics.mechanics.rigidbody import RigidBody

class Kane(object):
    """Kane's method object.

    This object is used to do the "book-keeping" as you go through and form
    equations of motion in the way Kane presents in:
    Kane, T., Levinson, D. Dynamics Theory and Applications. 1985 McGraw-Hill

    The attributes are for equations in the form [M] udot = forcing


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
    simp : Boolean
        Flag determining whether simplification of symbolic matrix
        inversion can occur or not
    mass_matrix_full : Matrix
        The "mass matrix" for the u's and q's
    forcing_full : Matrix
        The "forcing vector" for the u's and q's

    Examples
    ========

    This is a simple example for a one defree of freedom translational
    spring-mass-damper.

    In this example, we first need to do the kinematics.
    This involves creating generalized speeds and coordinates and their
    derivatives.
    Then we create a point and set its velocity in a frame::

        >>> from sympy import symbols
        >>> from sympy.physics.mechanics import dynamicsymbols, ReferenceFrame
        >>> from sympy.physics.mechanics import Point, Particle, Kane
        >>> q, u = dynamicsymbols('q u')
        >>> qd, ud = dynamicsymbols('q u', 1)
        >>> m, c, k = symbols('m c k')
        >>> N = ReferenceFrame('N')
        >>> P = Point('P')
        >>> P.set_vel(N, u * N.x)

    Next we need to arrange/store information in the way the Kane requires.
    The kinematic differential equations need to be stored in a dict.
    A list of forces/torques must be constructed, where each entry in the list
    is a (Point, Vector) or (ReferenceFrame, Vector) tuple, where the Vectors
    represent the Force or Torque.
    Next a particle needs to be created, and it needs to have a point and mass
    assigned to it.
    Finally, a list of all bodies and particles needs to be created::

        >>> kd = [qd - u]
        >>> FL = [(P, (-k * q - c * u) * N.x)]
        >>> pa = Particle()
        >>> pa.mass = m
        >>> pa.point = P
        >>> BL = [pa]

    Finally we can generate the equations of motion.
    First we create the Kane object and supply an inertial frame.
    Next we pass it the generalized speeds.
    Then we pass it the kinematic differential equation dict.
    Next we form FR* and FR to complete: Fr + Fr* = 0.
    We have the equations of motion at this point.
    It makes sense to rearrnge them though, so we calculate the mass matrix and
    the forcing terms, for E.o.M. in the form: [MM] udot = forcing, where MM is
    the mass matrix, udot is a vector of the time derivatives of the
    generalized speeds, and forcing is a vector representing "forcing" terms::

        >>> KM = Kane(N)
        >>> KM.coords([q])
        >>> KM.speeds([u])
        >>> KM.kindiffeq(kd)
        >>> (fr, frstar) = KM.kanes_equations(FL, BL)
        >>> MM = KM.mass_matrix
        >>> forcing = KM.forcing
        >>> rhs = MM.inv() * forcing
        >>> rhs
        [-(c*u(t) + k*q(t))/m]
        >>> KM.linearize()[0]
        [0, 1]
        [k, c]

    Please look at the documentation pages for more information on how to
    perform linearization and how to deal with dependent coordinates & speeds,
    and how do deal with bringing non-contributing forces into evidence.

    """

    simp = False

    def __init__(self, frame):
        """Supply the inertial frame for Kane initialization. """
        # Big storage things
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

        # Differential Equations Matrices
        self._k_d = None
        self._f_d = None
        self._k_kqdot = None
        self._k_ku = None
        self._f_k = None

        # Constraint Matrices
        self._f_h = Matrix([])
        self._k_nh = Matrix([])
        self._f_nh = Matrix([])
        self._k_dnh = Matrix([])
        self._f_dnh = Matrix([])

    def _find_dynamicsymbols(self, inlist, insyms=[]):
        """Finds all non-supplied dynamicsymbols in the expressions."""
        from sympy.core.function import UndefinedFunction, Derivative
        t = dynamicsymbols._t

        def _deeper(iexpr):
            oli = []
            if isinstance(type(iexpr), UndefinedFunction):
                if iexpr.args == (t,):
                    oli += [iexpr]
            elif isinstance(iexpr, Derivative):
                if (bool([i == t for i in iexpr.variables]) &
                    isinstance(type(iexpr.args[0]), UndefinedFunction)):
                    ol = str(iexpr.args[0].func)
                    for i, v in enumerate(iexpr.variables):
                        ol += '\''
                    oli += [iexpr]
            else:
                for i, v in enumerate(iexpr.args):
                    oli += _deeper(v)
            return oli

        ol = []
        for i in list(inlist):
            ol += _deeper(i)
        seta = {}
        map(seta.__setitem__, ol, [])
        ol = seta.keys()
        for i, v in enumerate(insyms):
            if ol.__contains__(v):
                ol.remove(v)
        return ol

    def _find_othersymbols(self, inlist, insyms=[]):
        """Finds all non-dynamic symbols in the expressions."""
        def _deeper(iexpr):
            oli = []
            if isinstance(iexpr, Symbol):
                oli += [iexpr]
            else:
                for i, v in enumerate(iexpr.args):
                    oli += _deeper(v)
            return oli

        ol = []
        for i in list(inlist):
            ol += _deeper(i)
        seta = {}
        map(seta.__setitem__, ol, [])
        ol = seta.keys()
        for i, v in enumerate(insyms):
            if ol.__contains__(v):
                ol.remove(v)
        return ol

    def _mat_inv_mul(self, A, B):
        """Internal Function

        Computes A^-1 * B symbolically w/ substitution, where B is not
        necessarily a vector, but can be a matrix.

        """

        # Note: investigate difficulty in only creating symbols for non-zero
        # entries; this could speed things up, perhaps?

        r1, c1 = A.shape
        r2, c2 = B.shape
        temp1 = Matrix(r1, c1, lambda i, j: Symbol('x' + str(j + r1 * i)))
        temp2 = Matrix(r2, c2, lambda i, j: Symbol('y' + str(j + r2 * i)))
        for i in range(len(temp1)):
            if A[i] == 0:
                temp1[i] = 0
        for i in range(len(temp2)):
            if B[i] == 0:
                temp2[i] = 0
        temp3 = []
        for i in range(c2):
            temp3.append(temp1.LUsolve(temp2.extract(range(r2), [i])))
        temp3 = Matrix([i.T for i in temp3]).T
        if Kane.simp == True:
            temp3.simplify()
        return temp3.subs(dict(zip(temp1, A))).subs(dict(zip(temp2, B)))

    def coords(self, qind, qdep=[], coneqs=[]):
        """Supply all the generalized coordiantes in a list.

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
            raise TypeError('Generalized coords. must be supplied in a list')
        self._q = qind + qdep
        self._qdot = [diff(i, dynamicsymbols._t) for i in self._q]

        if not isinstance(qdep, (list, tuple)):
            raise TypeError('Dependent speeds and constraints must each be '
                            'provided in their own list')
        if len(qdep) != len(coneqs):
            raise ValueError('There must be an equal number of dependent '
                             'speeds and constraints')
        coneqs = Matrix(coneqs)
        self._qdep = qdep
        self._f_h = coneqs

    def speeds(self, uind, udep=[], coneqs=[], diffconeqs=None, u_auxiliary=[]):
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

        if not isinstance(uind, (list, tuple)):
            raise TypeError('Generalized speeds must be supplied in a list')
        self._u = uind + udep
        self._udot = [diff(i, dynamicsymbols._t) for i in self._u]
        self._uaux = u_auxiliary

        if not isinstance(udep, (list, tuple)):
            raise TypeError('Dependent speeds and constraints must each be '
                            'provided in their own list')
        if len(udep) != len(coneqs):
            raise ValueError('There must be an equal number of dependent '
                             'speeds and constraints')
        if diffconeqs != None:
            if len(udep) != len(diffconeqs):
                raise ValueError('There must be an equal number of dependent '
                                 'speeds and constraints')
        if len(udep) != 0:
            u = self._u
            uzero = dict(zip(u, [0] * len(u)))
            coneqs = Matrix(coneqs)
            udot = self._udot
            udotzero = dict(zip(udot, [0] * len(udot)))

            self._udep = udep
            self._f_nh = coneqs.subs(uzero)
            self._k_nh = (coneqs - self._f_nh).jacobian(u)
            # if no differentiated non holonomic constraints were given, calculate
            if diffconeqs == None:
                self._k_dnh = self._k_nh
                self._f_dnh = (self._k_nh.diff(dynamicsymbols._t) * Matrix(u) +
                               self._f_nh.diff(dynamicsymbols._t))
            else:
                self._f_dnh = diffconeqs.subs(udotzero)
                self._k_dnh = (diffconeqs - self._f_dnh).jacobian(udot)

            o = len(u) # number of generalized speeds
            m = len(udep) # number of motion constraints
            p = o - m # number of independent speeds
            # For a reminder, form of non-holonomic constraints is:
            # B u + C = 0
            B = self._k_nh.extract(range(m), range(o))
            C = self._f_nh.extract(range(m), [0])

            # We partition B into indenpendent and dependent columns
            # Ars is then -Bdep.inv() * Bind, and it relates depedent speeds to
            # independent speeds as: udep = Ars uind, neglecting the C term here.
            self._depB = B
            self._depC = C
            mr1 = B.extract(range(m), range(p))
            ml1 = B.extract(range(m), range(p, o))
            self._Ars = - self._mat_inv_mul(ml1, mr1)

    def kindiffdict(self):
        """Returns the qdot's in a dictionary. """
        if self._k_kqdot == None:
            raise ValueError('Kin. diff. eqs need to be supplied first')
        sub_dict = solve_linear_system_LU(Matrix([self._k_kqdot.T,
            -(self._k_ku * Matrix(self._u) + self._f_k).T]).T, self._qdot)
        return sub_dict

    def kindiffeq(self, kdeqs):
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
        uaz = dict(zip(uaux, [0] * len(uaux)))

        kdeqs = Matrix(kdeqs).subs(uaz)

        qdot = self._qdot
        qdotzero = dict(zip(qdot, [0] * len(qdot)))
        u = self._u
        uzero = dict(zip(u, [0] * len(u)))

        f_k = kdeqs.subs(uzero).subs(qdotzero)
        k_kqdot = (kdeqs.subs(uzero) - f_k).jacobian(Matrix(qdot))
        k_ku = (kdeqs.subs(qdotzero) - f_k).jacobian(Matrix(u))

        self._k_ku = self._mat_inv_mul(k_kqdot, k_ku)
        self._f_k = self._mat_inv_mul(k_kqdot, f_k)
        self._k_kqdot = eye(len(qdot))

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

        if not isinstance(fl, (list, tuple)):
            raise TypeError('Forces must be supplied in a list of: lists or '
                            'tuples.')
        N = self._inertial
        self._forcelist = fl[:]
        u = self._u
        o = len(u)

        FR = zeros(o, 1)
        # goes through each Fr (where this loop's i is r)
        for i, v in enumerate(u):
            # does this for each force pair in list (pair is w)
            for j, w in enumerate(fl):
                if isinstance(w[0], ReferenceFrame):
                    speed = w[0].ang_vel_in(N)
                    FR[i] += speed.diff(v, N) & w[1]
                elif isinstance(w[0], Point):
                    speed = w[0].vel(N)
                    FR[i] += speed.diff(v, N) & w[1]
                else:
                    raise TypeError('First entry in force pair is a point or'
                                    ' frame')
        # for dependent speeds
        if len(self._udep) != 0:
            m = len(self._udep)
            p = o - m
            FRtilde = FR.extract(range(p), [0])
            FRold = FR.extract(range(p, o), [0])
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

        if not isinstance(bl, (list, tuple)):
            raise TypeError('Bodies must be supplied in a list.')
        if self._fr == None:
            raise ValueError('Calculate Fr first, please')
        t = dynamicsymbols._t
        N = self._inertial
        self._bodylist = bl
        u = self._u # all speeds
        udep = self._udep # dependent speeds
        o = len(u)
        p = o - len(udep)
        udot = self._udot
        udotzero = dict(zip(udot, [0] * len(udot)))
        uaux = self._uaux
        uauxdot = [diff(i, t) for i in uaux]
        # dictionary of auxiliary speeds which are equal to zero
        uaz = dict(zip(uaux, [0] * len(uaux)))
        # dictionary of derivatives of auxiliary speeds which are equal to zero
        uadz = dict(zip(uauxdot, [0] * len(uauxdot)))

        # Form R*, T* for each body or particle in the list
        # This is stored as a list of tuples [(r*, t*),...]
        # Each tuple is for a body or particle
        # Within each rs is a tuple and ts is a tuple
        # These have the same structure: ([list], value)
        # The list is the coefficients of rs/ts wrt udots, value is everything
        # else in the expression
        # Partial velocities are stored as a list of tuple; a tuple for each
        # body
        # Each tuple has two elements, lists which represent the partial
        # velocity for each ur; The first list is translational partial
        # velocities, the second list is rotational translational velocities
        MM = zeros(o, o)
        nonMM = zeros(o, 1)
        rsts = []
        partials = []
        for i, v in enumerate(bl): # go through list of bodies, particles
            if isinstance(v, RigidBody):
                om = v.frame.ang_vel_in(N).subs(uadz).subs(uaz) # ang velocity
                omp = v.frame.ang_vel_in(N) # ang velocity, for partials
                alp = v.frame.ang_acc_in(N).subs(uadz).subs(uaz) # ang acc
                ve = v.mc.vel(N).subs(uadz).subs(uaz) # velocity
                vep = v.mc.vel(N) # velocity, for partials
                acc = v.mc.acc(N).subs(uadz).subs(uaz) # acceleration
                m = (v.mass).subs(uadz).subs(uaz)
                I, P = v.inertia
                I = I.subs(uadz).subs(uaz)
                if P != v.mc:
                    # redefine I about mass center
                    # have I S/O, want I S/S*
                    # I S/O = I S/S* + I S*/O; I S/S* = I S/O - I S*/O

                    # This block of code needs to have a test written for it
                    print('This functionality has not yet been tested yet, '
                          'use at your own risk')
                    f = v.frame
                    d = v.mc.pos_from(P)
                    I -= m * (((f.x | f.x) + (f.y | f.y) + (f.z | f.z)) *
                              (d & d) - (d | d))
                templist = []
                # One could think of r star as a collection of coefficients of
                # the udots plus another term. What we do here is get all of
                # these coefficients and store them in a list, then we get the
                # "other" term and put the list and other term in a tuple, for
                # each body/particle. The same is done for t star. The reason
                # for this is to not let the expressions get too large; so we
                # keep them seperate for as long a possible
                for j, w in enumerate(udot):
                    templist.append(-m * acc.diff(w, N))
                other = -m.diff(t) * ve - m * acc.subs(udotzero)
                rs = (templist, other)
                templist = []
                # see above note
                for j, w in enumerate(udot):
                    templist.append(-I & alp.diff(w, N))
                other = -((I.dt(v.frame) & om) + (I & alp.subs(udotzero))
                          + (om ^ (I & om)))
                ts = (templist, other)
                tl1 = []
                tl2 = []
                # calculates the partials only once and stores them for later
                for j, w in enumerate(u):
                    tl1.append(vep.diff(w, N))
                    tl2.append(omp.diff(w, N))
                partials.append((tl1, tl2))

            elif isinstance(v, Particle):
                ve = v.point.vel(N).subs(uadz).subs(uaz)
                vep = v.point.vel(N)
                acc = v.point.acc(N).subs(uadz).subs(uaz)
                m = v.mass.subs(uadz).subs(uaz)
                templist = []
                # see above note
                for j, w in enumerate(udot):
                    templist.append(-m * acc.diff(w, N))
                other = -m.diff(t) * ve - m * acc.subs(udotzero)
                rs = (templist, other)
                # We make an empty t star here so that way the later code
                # doesn't care whether its operating on a body or particle
                ts = ([0] * len(u), 0)
                tl1 = []
                tl2 = []
                # calculates the partials only once, makes 0's for angular
                # partials so the later code is body/particle indepedent
                for j, w in enumerate(u):
                    tl1.append(vep.diff(w, N))
                    tl2.append(0)
                partials.append((tl1, tl2))
            else:
                raise TypeError('The body list needs RigidBody or '
                                'Particle as list elements')
            rsts.append((rs, ts))

        # Use R*, T* and partial velocities to form FR*
        FRSTAR = zeros(o, 1)
        # does this for each body in the list
        for i, v in enumerate(rsts):
            rs, ts = v # unpact r*, t*
            vps, ops = partials[i] # unpack vel. partials, ang. vel. partials
            # Computes the mass matrix entries from r*, there are from the list
            # in the rstar tuple
            ii = 0
            for x in vps:
                for w in rs[0]:
                    MM[ii] += w & x
                    ii += 1
            # Computes the mass matrix entries from t*, there are from the list
            # in the tstar tuple
            ii = 0
            for x in ops:
                for w in ts[0]:
                    MM[ii] += w & x
                    ii += 1
            # Non mass matrix entries from rstar, from the other in the rstar
            # tuple
            for j, w in enumerate(vps):
                nonMM[j] += w & rs[1]
            # Non mass matrix entries from tstar, from the other in the tstar
            # tuple
            for j, w in enumerate(ops):
                nonMM[j] += w & ts[1]
        FRSTAR = MM * Matrix(udot) + nonMM

        # For motion constraints, m is the number of constraints
        # Really, one should just look at Kane's book for descriptions of this
        # process
        if len(self._udep) != 0:
            FRSTARtilde = FRSTAR.extract(range(p), [0])
            FRSTARold = FRSTAR.extract(range(p, o), [0])
            FRSTARtilde += self._Ars.T * FRSTARold
            FRSTAR = FRSTARtilde

            MMi = MM.extract(range(p), range(o))
            MMd = MM.extract(range(p, o), range(o))
            MM = MMi + self._Ars.T * MMd
        self._frstar = FRSTAR

        zeroeq = self._fr + self._frstar
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

        if (self._q == None) or (self._u == None):
            raise ValueError('Speeds and coordinates must be supplied first')
        if (self._k_kqdot == None):
            raise ValueError('Supply kinematic differential equations please')


        fr = self._form_fr(FL)
        frstar = self._form_frstar(BL)
        if self._uaux != []:
            km = Kane(self._inertial)
            km.coords(self._q)
            km.speeds(self._uaux, u_auxiliary=self._uaux)
            fraux = km._form_fr(FL)
            frstaraux = km._form_frstar(BL)
            self._aux_eq = fraux + frstaraux
            self._fr = fr.col_join(fraux)
            self._frstar = frstar.col_join(frstaraux)
            return (self._fr, self._frstar)
        else:
            return (fr, frstar)

    @property
    def auxiliary_eqs(self):
        if (self._fr == None) or (self._frstar == None):
            raise ValueError('Need to compute Fr, Fr* first')
        if self._uaux == []:
            raise ValueError('No auxiliary speeds have been declared')
        return self._aux_eq

    def linearize(self):
        """ Method used to generate linearized equations.

        Note that for linearization, it is assumed that time is not perturbed,
        but only coordinates and positions. The "forcing" vector's jacobian is
        computed with respect to the state vector in the form [Qi, Qd, Ui, Ud].
        This is the "A" matrix.

        It also finds any non-state dynamicsymbols and computes the jacobian of
        the "forcing" vector with respect to them. This is the "B" matrix; if
        this is empty, an empty matrix is created

        A tuple of ("A", "B") is returned.

        """

        if (self._fr == None) or (self._frstar == None):
            raise ValueError('Need to compute Fr, Fr* first')

        # Note that this is now unneccessary, and it should never be
        # encountered; I still think it should be in here in case the user
        # manually sets these matrices incorrectly.
        for i in self._q:
            if self._k_kqdot.diff(i) != 0 * self._k_kqdot:
                raise ValueError('Matrix K_kqdot must not depend on any q')

        t = dynamicsymbols._t
        uaux = self._uaux
        uauxdot = [diff(i, t) for i in uaux]
        # dictionary of auxiliary speeds which are equal to zero
        uaz = dict(zip(uaux, [0] * len(uaux)))
        # dictionary of derivatives of auxiliary speeds which are equal to zero
        uadz = dict(zip(uauxdot, [0] * len(uauxdot)))

        # Checking for dynamic symbols outside the dynamic differential
        # equations; throws error if there is.
        insyms = self._q + self._qdot + self._u + self._udot + uaux + uauxdot
        bad_oths= self._find_dynamicsymbols(self._k_kqdot, insyms)
        bad_oths += self._find_dynamicsymbols(self._k_ku, insyms)
        bad_oths += self._find_dynamicsymbols(self._f_k, insyms)
        bad_oths += self._find_dynamicsymbols(self._k_dnh, insyms)
        bad_oths += self._find_dynamicsymbols(self._f_dnh, insyms)
        bad_oths += self._find_dynamicsymbols(self._k_d, insyms)
        if bad_oths != []:
            raise ValueError('Cannot have dynamic symbols outside dynamic ' +
                             'forcing vector')
        oth_dyns = self._find_dynamicsymbols(self._f_d.subs(uadz).subs(uaz),
                                             insyms)
        for i in oth_dyns:
            if diff(i, dynamicsymbols._t) in oth_dyns:
                raise ValueError('Cannot have derivatives of forcing terms ' +
                                 'when linearizing forcing terms')

        o = len(self._u) # number of speeds
        n = len(self._q) # number of coordinates
        l = len(self._qdep) # number of configuration constraints
        m = len(self._udep) # number of motion constraints
        qi = Matrix(self._q[0: n - l]) # independent coords
        qd = Matrix(self._q[n - l: n]) # dependent coords; could be empty
        ui = Matrix(self._u[0: o - m]) # independent speeds
        ud = Matrix(self._u[o - m: o]) # dependent speeds; could be empty

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
        f1 = f1.subs(uadz).subs(uaz)
        f2 = f2.subs(uadz).subs(uaz)
        fh = self._f_h.subs(uadz).subs(uaz)
        fku = (self._k_ku * Matrix(self._u)).subs(uadz).subs(uaz)
        fkf = self._f_k.subs(uadz).subs(uaz)

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
        # seperating into dqd (depedent q's) and dqi (independent q's) and the
        # rearranging for dqd/dqi. This is again extended for the speeds.

        # First case: configuration and motion constraints
        if (l != 0) and (m != 0):
            dqd_dqi = - self._mat_inv_mul(fh.jacobian(qd), fh.jacobian(qi))
            dud_dqi = self._mat_inv_mul(fnh.jacobian(ud), (fnh.jacobian(qd) *
                                        dqd_dqi - fnh.jacobian(qi)))
            dud_dui = - self._mat_inv_mul(fnh.jacobian(ud), fnh.jacobian(ui))
            dqdot_dui = - self._k_kqdot.inv() * (fku.jacobian(ui) +
                                                fku.jacobian(ud) * dud_dui)
            dqdot_dqi = - self._k_kqdot.inv() * (fku.jacobian(qi) +
                    fkf.jacobian(qi) + (fku.jacobian(qd) + fkf.jacobian(qd)) *
                    dqd_dqi + fku.jacobian(ud) * dud_dqi)
            f1_q = (f1.jacobian(qi) + f1.jacobian(qd) * dqd_dqi +
                    f1.jacobian(ud) * dud_dqi)
            f1_u = (f1.jacobian(ui) + f1.jacobian(ud) * dud_dui)
            f2_q = (f2.jacobian(qi) + f2.jacobian(qd) * dqd_dqi +
                    f2.jacobian(self._qdot) * dqdot_dqi + f2.jacobian(ud) *
                    dud_dqi)
            f2_u = (f2.jacobian(ui) + f2.jacobian(ud) * dud_dui +
                    f2.jacobian(self._qdot) * dqdot_dui)
        # Second case: configuration constraints only
        elif l != 0:
            dqd_dqi = - self._mat_inv_mul(fh.jacobian(qd), fh.jacobian(qi))
            dqdot_dui = - self._k_kqdot.inv() * fku.jacobian(ui)
            dqdot_dqi = - self._k_kqdot.inv() * (fku.jacobian(qi) +
                fkf.jacobian(qi) + (fku.jacobian(qd) + fkf.jacobian(qd)) *
                    dqd_dqi)
            f1_q = (f1.jacobian(qi) + f1.jacobian(qd) * dqd_dqi)
            f1_u = f1.jacobian(ui)
            f2_q = (f2.jacobian(qi) + f2.jacobian(qd) * dqd_dqi +
                    f2.jacobian(self._qdot) * dqdot_dqi)
            f2_u = (f2.jacobian(ui) + f2.jacobian(self._qdot) * dqdot_dui)
        # Third case: motion constraints only
        elif m != 0:
            dud_dqi = self._mat_inv_mul(fnh.jacobian(ud), - fnh.jacobian(qi))
            dud_dui = - self._mat_inv_mul(fnh.jacobian(ud), fnh.jacobian(ui))
            dqdot_dui = - self._k_kqdot.inv() * (fku.jacobian(ui) +
                                                fku.jacobian(ud) * dud_dui)
            dqdot_dqi = - self._k_kqdot.inv() * (fku.jacobian(qi) +
                    fkf.jacobian(qi) + fku.jacobian(ud) * dud_dqi)
            f1_q = (f1.jacobian(qi) + f1.jacobian(ud) * dud_dqi)
            f1_u = (f1.jacobian(ui) + f1.jacobian(ud) * dud_dui)
            f2_q = (f2.jacobian(qi) + f2.jacobian(self._qdot) * dqdot_dqi +
                    f2.jacobian(ud) * dud_dqi)
            f2_u = (f2.jacobian(ui) + f2.jacobian(ud) * dud_dui +
                    f2.jacobian(self._qdot) * dqdot_dui)
        # Fourth case: No constraints
        else:
            dqdot_dui = - self._k_kqdot.inv() * fku.jacobian(ui)
            dqdot_dqi = - self._k_kqdot.inv() * (fku.jacobian(qi) +
                    fkf.jacobian(qi))
            f1_q = f1.jacobian(qi)
            f1_u = f1.jacobian(ui)
            f2_q = f2.jacobian(qi) + f2.jacobian(self._qdot) * dqdot_dqi
            f2_u = f2.jacobian(ui) + f2.jacobian(self._qdot) * dqdot_dui
        Amat = -(f1_q.row_join(f1_u)).col_join(f2_q.row_join(f2_u))
        if oth_dyns != []:
            f1_oths = f1.jacobian(oth_dyns)
            f2_oths = f2.jacobian(oth_dyns)
            Bmat = -f1_oths.col_join(f2_oths)
        else:
            Bmat = Matrix([])
        return (Amat, Bmat)

    @property
    def mass_matrix(self):
        # Returns the mass matrix, which is augmented by the differentiated non
        # holonomic equations if necessary
        if (self._frstar == None) & (self._fr == None):
            raise ValueError('Need to compute Fr, Fr* first')
        return Matrix([self._k_d, self._k_dnh])

    @property
    def mass_matrix_full(self):
        # Returns the mass matrix from above, augmented by kin diff's k_kqdot
        if (self._frstar == None) & (self._fr == None):
            raise ValueError('Need to compute Fr, Fr* first')
        o = len(self._u)
        n = len(self._q)
        return ((self._k_kqdot).row_join(zeros(n, o))).col_join((zeros(o,
                n)).row_join(self.mass_matrix))

    @property
    def forcing(self):
        # Returns the forcing vector, which is augmented by the differentiated
        # non holonomic equations if necessary
        if (self._frstar == None) & (self._fr == None):
            raise ValueError('Need to compute Fr, Fr* first')
        return -Matrix([self._f_d, self._f_dnh])

    @property
    def forcing_full(self):
        # Returns the forcing vector, which is augmented by the differentiated
        # non holonomic equations if necessary
        if (self._frstar == None) & (self._fr == None):
            raise ValueError('Need to compute Fr, Fr* first')
        f1 = self._k_ku * Matrix(self._u) + self._f_k
        return -Matrix([f1, self._f_d, self._f_dnh])

