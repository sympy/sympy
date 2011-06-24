__all__ = ['Kane']

from sympy import Symbol, zeros, simplify, Matrix, solve_linear_system_LU
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

    >>> kd = {qd: u}
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
    >>> KM.gen_speeds([u])
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
        self._inertial = frame
        self._us = None
        self._qs = None
        self._qdots = None
        self._kd = dict()
        self._forcelist = None
        self._bodylist = None
        self._fr = None
        self._frstar = None
        self._uds = []
        self._dep_cons = []
        self._forcing = None
        self._massmatrix = None
        self._rhs = None

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

    def gen_speeds(self, inlist):
        """Supply all the generalized speeds in a list. """
        if not isinstance(inlist, (list, tuple)):
            raise TypeError('Generalized Speeds must be supplied in a list')
        self._us = inlist

    def kindiffeq(self, indict):
        """Supply all the kinematic differential equations in a dict form. """
        if not isinstance(indict, dict):
            raise TypeError('Kinematic Differential Equations must be '
                            'supplied in a dictionary')
        self._kd = indict
        newlist = indict.keys()
        for i, v in enumerate(newlist):
            newlist[i] = DynamicSymbol(v.name[:-1])
        self._qs = newlist
        self._qdots = indict.keys()

    def dependent_speeds(self, uds, conl):
        """This is used to when dealing with systems with motion constraints.

        There will be much more documentation for this in the future.

        """

        if not isinstance(uds, (list, tuple)):
            raise TypeError('Dependent speeds and constraints must each be '
                            'provided in their own list')
        if not uds.__len__() == conl.__len__():
            raise ValueError('There must be an equal number of dependent '
                             'speeds and constraints')
        self._uds = uds
        for i, v in enumerate(conl):
            conl[i] = v.subs(self._kd)
        self._dep_cons = conl
        uis = self._us

        # non coord/speed DynamicSymbols
        oth = self._find_others(conl, uis + uds + self._qs + self._qdots)

        n = len(uis)
        m = len(uds)
        p = n - m
        # puts independent speeds first
        for i, v in enumerate(uds):
            uis.remove(v)
        uis += uds
        self._us = uis
        B = zeros((m, n))
        C = zeros((m, 1))
        uzeros = {}
        for i, v in enumerate(uis):
            uzeros.update({v: 0})
        ii = 0
        for i1, v1 in enumerate(uds):
            for i2, v2 in enumerate(uis):
                B[ii] = conl[i1].diff(uis[i2])
                ii += 1
        for i1, v1 in enumerate(uds):
                C[i1] = conl[i1].subs(uzeros)
        mr1 = B.extract(range(m), range(p))
        ml1 = B.extract(range(m), range(p, n))
        self._depB = B
        self._depC = C
        ml1i = ml1.inverse_ADJ()
        self._Ars = ml1i * mr1
        self._Brs = ml1i * C

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
        uis = self._us
        uds = self._uds
        n = len(uis)

        FR = zeros((n, 1))
        # goes through each Fr (where this loop's i is r)
        for i, v in enumerate(uis):
            # does this for each force pair in list (pair is w)
            for j, w in enumerate(fl):
                if isinstance(w[0], ReferenceFrame):
                    speed = w[0].ang_vel_in(N).subs(self._kd)
                    FR[i] += speed.diff(v, N) & w[1].subs(self._kd)
                elif isinstance(w[0], Point):
                    speed = w[0].vel(N).subs(self._kd)
                    FR[i] += speed.diff(v, N) & w[1].subs(self._kd)
                else:
                    raise TypeError('First entry in force pair is a point or'
                                    ' frame')
        # for dependent speeds
        if len(uds) != 0:
            m = len(uds)
            p = n - m
            FRtilde = FR.extract(range(p), [0])
            FRold = FR.extract(range(p, n), [0])
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
            A list of all RigidBody's and Point's in the system.

        """

        if not isinstance(bl, (list, tuple)):
            raise TypeError('Bodies must be supplied in a list.')
        N = self._inertial
        self._bodylist = bl[:]
        uis = self._us
        uds = self._uds
        n = len(uis)
        udots = []
        udotszero = {}
        for i, v in enumerate(uis):
            udots.append(v.diff(Symbol('t')))
            udotszero.update({udots[i]: 0})
        # Form R*, T* for each body or particle in the list
        MM = zeros((n, n))
        nonMM = zeros((n, 1))
        rsts = []
        partials = []
        for i, v in enumerate(bl):
            if isinstance(v, RigidBody):
                om = v.frame.ang_vel_in(N).subs(self._kd)
                alp = v.frame.ang_acc_in(N).subs(self._kd)
                ve = v.mc.vel(N).subs(self._kd)
                acc = v.mc.acc(N).subs(self._kd)
                m = v.mass
                I, p = v.inertia
                if p != v.mc:
                    # redefine I about mass center
                    # have I S/O, want I S/S*
                    # I S/O = I S/S* + I S*/O; I S/S* = I S/O - I S*/O
                    f = v.frame
                    d = v.mc.pos_from(p)
                    I -= m * (((f.x | f.x) + (f.y | f.y) + (f.z | f.z)) *
                              (d & d) - (d | d))
                templist = []
                for j, w in enumerate(udots):
                    templist.append(-m * acc.diff(w, N))
                other = -m.diff(Symbol('t')) * ve - m * acc.subs(udotszero)
                rs = (templist, other)
                templist = []
                for j, w in enumerate(udots):
                    templist.append(-I & alp.diff(w, N))
                other = -((I.dt(v.frame) & om) + (I & alp.subs(udotszero)) + (om
                          ^ (I & om)))
                ts = (templist, other)
                tl1 = []
                tl2 = []
                for j, w in enumerate(uis):
                    tl1.append(ve.diff(w, N))
                    tl2.append(om.diff(w, N))
                partials.append((tl1, tl2))

            elif isinstance(v, Particle):
                ve = v.point.vel(N).subs(self._kd)
                acc = v.point.acc(N).subs(self._kd)
                m = v.mass
                templist = []
                for j, w in enumerate(udots):
                    templist.append(-m * acc.diff(w, N))
                other = -m.diff(Symbol('t')) * ve - m * acc.subs(udotszero)
                rs = (templist, other)
                ts = (udotszero.values(),0)
                tl1 = []
                tl2 = []
                for j, w in enumerate(uis):
                    tl1.append(ve.diff(w, N))
                    tl2.append(0)
                partials.append((tl1, tl2))
            else:
                raise TypeError('The body list needs RigidBody or '
                                'Particle as list elements')
            rsts.append((rs, ts))

        # Use R*, T* and partial velocities to form FR*
        FRSTAR = zeros((n, 1))
        # does this for each body in the list
        for i, v in enumerate(bl):
            rs, ts = rsts[i]
            vps, ops = partials[i]
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
        FRSTAR = MM * Matrix(udots) + nonMM

        if len(uds) != 0:
            m = len(uds)
            p = n - m
            FRSTARtilde = FRSTAR.extract(range(p), [0])
            FRSTARold = FRSTAR.extract(range(p, n), [0])
            FRSTARtilde += self._Ars.T * FRSTARold
            FRSTAR = FRSTARtilde
        self._frstar = FRSTAR
        self._massmatrix = -MM
        return FRSTAR

    def mass_matrix(self):
        """Computes the system's mass matrix.

        The mass matrix for the system represents a generalized m, in a
        multibody interpretation of f = m a.

        When working this way, the E.o.M. in the form: [MM] udot = forcing,
        where MM is the mass matrix, udot is a vector of the time derivatives
        of the generalized speeds, and forcing is a vector representing
        various "forcing" terms.

        """

        if (self._frstar == None) & (self._fr == None):
            raise ValueError('Need to compute FR* first')
        if len(self._uds) == 0:
            return self._massmatrix
        uis = self._us
        uds = self._uds
        n = len(uis)
        m = len(uds)
        p = n - m
        MM = self._massmatrix
        MMi = MM.extract(range(p), range(n))
        MMd = MM.extract(range(p, n), range(n))
        MM = MMi + self._Ars.T * MMd
        MM = Matrix([MM, self._depB])
        self._massmatrix = MM
        return MM

    def forcing(self):
        """Computes the system's forcing terms

        The mass matrix for the system represents a generalized f, in a
        multibody interpretation of f = m a.

        When working this way, the E.o.M. in the form: [MM] udot = forcing,
        where MM is the mass matrix, udot is a vector of the time derivatives
        of the generalized speeds, and forcing is a vector representing
        various "forcing" terms.

        """

        if (self._frstar == None) & (self._fr == None):
            raise ValueError('Need to compute FR* first')
        t = Symbol('t')
        uis = self._us
        uds = self._uds
        n = len(uis)
        m = len(uds)
        p = n - m
        udotszero = {}
        for i, v in enumerate(uis):
            udotszero.update({v.diff(t): 0})
        zeroeq = self._fr + self._frstar
        zeroeq = zeroeq.subs(udotszero)
        if len(self._uds) == 0:
            return zeroeq
        extra = - self._depB.diff(t) * Matrix(uis) - self._depC.diff(t)
        self._forcing = None
        return Matrix([zeroeq, extra])

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
