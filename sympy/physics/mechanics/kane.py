__all__ = ['Kane']

from sympy import Symbol, zeros, simplify, Matrix, diff
from sympy.physics.mechanics.essential import ReferenceFrame
from sympy.physics.mechanics.point import Point
from sympy.physics.mechanics.dynamicsymbol import DynamicSymbol
from sympy.physics.mechanics.rigidbody import RigidBody
from sympy.physics.mechanics.particle import Particle

class Kane(object):
    """Kane's method object.

    This object is used to do the "book-keeping" as you go through and form
    equations of motion in the way Kane presents in:
    Kane, T., Levinson, D. Dynamics Theory and Applications. 1985 McGraw-Hill

    The attributes are for equations in the form [M] udot = forcing

    Attributes
    ==========
    mass_matrix : Matrix
        The system's mass matrix
    forcing : Matrix
        The system's forcing vector


    A simple example for a one degree of freedom translational
    spring-mass-damper follows

    Example
    =======

    In this example, we first need to do the kinematics.
    This involves creating generalized speeds and coordinates and their
    derivatives.
    Then we create a point and set its velocity in a frame.

    >>> q, qd, u, ud = dynamicsymbols('q qd u ud')
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
    Finally, a list of all bodies and particles needs to be created.

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
    generalized speeds, and forcing is a vector representing "forcing" terms.

    >>> KM = Kane(N)
    >>> KM.coords([q])
    >>> KM.speeds([u])
    >>> KM.kindiffeq(kd)
    >>> frstar = KM.form_frstar(BL)
    >>> fr = KM.form_fr(FL)
    >>> MM = KM.mass_matrix()
    >>> forcing = KM.forcing()
    >>> rhs = MM.inv() * forcing
    >>> rhs
    [(-q*k - u*c)/m]

    """

    def __init__(self, frame):
        """Supply the inertial frame for Kane initialization. """
        # Big storage things
        self._inertial = frame
        self._forcelist = None
        self._bodylist = None
        self._fr = None
        self._frstar = None
        self._rhs = None

        # States
        self._q = None
        self._qdep = None
        self._qdot = None
        self._u = None
        self._udep = []
        self._udot = None

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

    def _find_others(self, inlist, insyms):
        """Finds all non-supplied DynamicSymbols in the list of expressions."""
        def _deeper(inexpr):
            oli = []
            try:
                for i, v in enumerate(inexpr.args):
                    oli += _deeper(v)
            except:
                if isinstance(inexpr, DynamicSymbol):
                    oli += inexpr
            return oli

        ol = []
        for i, v in enumerate(inlist):
            try:
                for i2, v2 in enumerate(v.args):
                    ol += _deeper(v2)
            except:
                if isinstance(v, DynamicSymbol):
                    ol += inexpr
        seta = {}
        map(seta.__setitem__, ol, [])
        ol = seta.keys()
        for i, v in enumerate(insyms):
            if ol.__contains__(v):
                ol.remove(v)
        return ol

    def coords(self, inlist):
        """Supply all the generalized coordinates in a list. """
        if not isinstance(inlist, (list, tuple)):
            raise TypeError('Generalized coords. must be supplied in a list')
        self._q = inlist
        self._qdot = [diff(i, Symbol('t')) for i in inlist]

    def speeds(self, inlist):
        """Supply all the generalized speeds in a list.

        If there are motion constraints, the order needs to be:
        [Ui, Ud], where Ui is independent speeds & Ud is dependent speeds.

        Parameters
        ==========
        inlist : list
            A list of generalized speeds

        """

        if not isinstance(inlist, (list, tuple)):
            raise TypeError('Generalized speeds must be supplied in a list')
        self._u = inlist
        self._udot = [diff(i, Symbol('t')) for i in inlist]

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

        kdeqs = Matrix(kdeqs)

        qdot = self._qdot
        qdotzero = dict(zip(qdot, [0] * len(qdot)))
        u = self._u
        uzero = dict(zip(u, [0] * len(u)))

        self._f_k = kdeqs.subs(uzero).subs(qdotzero)
        self._k_kqdot = (kdeqs.subs(uzero) - self._f_k).jacobian(Matrix(qdot))
        self._k_ku = (kdeqs.subs(qdotzero) - self._f_k).jacobian(Matrix(u))

    def dependent_speeds(self, udep, coneq, diffconeq=None):
        """This is used to when dealing with systems with motion constraints.

        There will be much more documentation for this in the future.

        """

        if not isinstance(udep, (list, tuple)):
            raise TypeError('Dependent speeds and constraints must each be '
                            'provided in their own list')
        if len(udep) != len(coneq):
            raise ValueError('There must be an equal number of dependent '
                             'speeds and constraints')
        if diffconeq != None:
            if len(udep) != len(diffconeq):
                raise ValueError('There must be an equal number of dependent '
                                 'speeds and constraints')
        for i in range(len(udep)):
            if self._u[-(i + 1)] != udep[-(i + 1)]:
                raise ValueError('The order of dependent speeds does not match'
                                 ' previously specified speeds')

        u = self._u
        uzero = dict(zip(u, [0] * len(u)))
        coneq = Matrix(coneq)
        udot = self._udot
        udotzero = dict(zip(udot, [0] * len(udot)))

        self._udep = udep
        self._f_nh = coneq.subs(uzero)
        self._k_nh = (coneq - self._f_nh).jacobian(u)
        if diffconeq == None:
            self._k_dnh = self._k_nh
            self._f_dnh = (self._k_nh.diff(Symbol('t')) * Matrix(u) +
                           self._f_nh.diff(Symbol('t')))
        else:
            self._f_dnh = diffconeq.subs(udotzero)
            self._k_dnh = (diffconeq - self._f_dnh).jacobian(udot)

        o = len(u) # number of generalized speeds
        m = len(udep) # number of motion constraints
        p = o - m # number of independent speeds
        # puts independent speeds first
        B = self._k_nh.extract(range(m), range(o))
        C = self._f_nh.extract(range(m), [0])

        mr1 = B.extract(range(m), range(p))
        ml1 = B.extract(range(m), range(p, o))
        self._depB = B
        self._depC = C
        ml1i = ml1.inverse_ADJ()
        self._Ars = ml1i * mr1

    def form_fr(self, fl):
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

        FR = zeros((o, 1))
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

    def form_frstar(self, bl):
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
        N = self._inertial
        self._bodylist = bl
        u = self._u # all speeds
        udep = self._udep # dependent speeds
        o = len(u)
        p = o - len(udep)
        udot = self._udot
        udots = []
        udotzero = dict(zip(udot, [0] * len(udot)))

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
        MM = zeros((o, o))
        nonMM = zeros((o, 1))
        rsts = []
        partials = []
        for i, v in enumerate(bl):
            if isinstance(v, RigidBody):
                om = v.frame.ang_vel_in(N)
                alp = v.frame.ang_acc_in(N)
                ve = v.mc.vel(N)
                acc = v.mc.acc(N)
                m = v.mass
                I, P = v.inertia
                if P != v.mc:
                    # redefine I about mass center
                    # have I S/O, want I S/S*
                    # I S/O = I S/S* + I S*/O; I S/S* = I S/O - I S*/O
                    f = v.frame
                    d = v.mc.pos_from(p)
                    I -= m * (((f.x | f.x) + (f.y | f.y) + (f.z | f.z)) *
                              (d & d) - (d | d))
                templist = []
                for j, w in enumerate(udot):
                    templist.append(-m * acc.diff(w, N))
                other = -m.diff(Symbol('t')) * ve - m * acc
                rs = (templist, other)
                templist = []
                for j, w in enumerate(udot):
                    templist.append(-I & alp.diff(w, N))
                other = -((I.dt(v.frame) & om) + (I & alp) + (om ^ (I & om)))
                ts = (templist, other)
                tl1 = []
                tl2 = []
                for j, w in enumerate(u):
                    tl1.append(ve.diff(w, N))
                    tl2.append(om.diff(w, N))
                partials.append((tl1, tl2))

            elif isinstance(v, Particle):
                ve = v.point.vel(N)
                acc = v.point.acc(N)
                m = v.mass
                templist = []
                for j, w in enumerate(udot):
                    templist.append(-m * acc.diff(w, N))
                other = -m.diff(Symbol('t')) * ve - m * acc
                rs = (templist, other)
                ts = ([0] * len(u), 0)
                tl1 = []
                tl2 = []
                for j, w in enumerate(u):
                    tl1.append(ve.diff(w, N))
                    tl2.append(0)
                partials.append((tl1, tl2))
            else:
                raise TypeError('The body list needs RigidBody or '
                                'Particle as list elements')
            rsts.append((rs, ts))

        # Use R*, T* and partial velocities to form FR*
        FRSTAR = zeros((o, 1))
        # does this for each body in the list
        for i, v in enumerate(rsts):
            rs, ts = v
            vps, ops = partials[i] # velocity partials, angular vel. partials
            ii = 0
            for j, w in enumerate(rs[0]):
                for k, x in enumerate(vps):
                    MM[ii] += w & x
                    ii += 1
            ii = 0
            for j, w in enumerate(ts[0]):
                for k, x in enumerate(ops):
                    MM[ii] += w & x
                    ii += 1
            for j, w in enumerate(vps):
                nonMM[j] += w & rs[1]
            for j, w in enumerate(ops):
                nonMM[j] += w & ts[1]
        FRSTAR = MM * Matrix(udot) + nonMM

        # For motion constraints, m is the number of constraints
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
        self._f_d = zeroeq.subs(udotzero)
        return FRSTAR

    @property
    def mass_matrix(self):
        if (self._frstar == None) & (self._fr == None):
            raise ValueError('Need to compute Fr, Fr* first')
        return Matrix([self._k_d, self._k_dnh])

    @property
    def forcing(self):
        if (self._frstar == None) & (self._fr == None):
            raise ValueError('Need to compute Fr, Fr* first')
        return -Matrix([self._f_d, self._f_dnh])

if __name__ == "__main__":
    import doctest
    from sympy import symbols
    from sympy.physics.mechanics import (dynamicsymbols, ReferenceFrame, Point,
                                         Particle, Kane)
    global_dict = {'symbols': symbols,
                   'dynamicsymbols': dynamicsymbols,
                   'ReferenceFrame': ReferenceFrame,
                   'Point': Point,
                   'Particle': Particle,
                   'Kane': Kane}
    doctest.testmod(globs=global_dict)
