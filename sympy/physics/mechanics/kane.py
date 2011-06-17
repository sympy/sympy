__all__ = ['Kane']

from sympy import Symbol, zeros, simplify, expand, Matrix
from sympy.physics.mechanics.essential import ReferenceFrame
from sympy.physics.mechanics.point import Point
from sympy.physics.mechanics.dynamicsymbol import DynamicSymbol
from sympy.physics.mechanics.rigidbody import RigidBody

class Kane(object):
    """Kane's method object. """

    def __init__(self, frame):
        """Supply the inertial frame. """
        self._inertial = frame
        self._us = None
        self._udots = None
        self._qs = None
        self._qdots = None
        self._kd = dict()
        self._forcelist = None
        self._bodylist = None
        self._fr = None
        self._frstar = None
        self._uds = []
        self._dep_cons = []

    def _find_others(self, inlist, insyms):
        """Finds all non-supplied DynamicSymbols in the list of expressions"""

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
        if not isinstance(inlist, (list, tuple)):
            raise TypeError('Generalized Speeds must be supplied in a list')
        self._us = inlist
        ol = []
        for i, v in enumerate(inlist):
            ol.append(v.diff(Symbol('t')))
        self._udots = ol

    def kindiffeq(self, indict):
        if not isinstance(indict, dict):
            raise TypeError('Kinematic Differential Equations must be '
                            'supplied in a dictionary')
        self._kd = indict
        newlist = indict.keys()
        for i, v in enumerate(newlist):
            newlist[i] = DynamicSymbol(v.name[:-1])
        self._qs = newlist
        self._qdots = indict.keys()

    def dependent_speeds(self, speedl, conl):
        if not isinstance(speedl, (list, tuple)):
            raise TypeError('Dependent speeds and constraints must each be '
                            'provided in their own list')
        if not speedl.__len__() == conl.__len__():
            raise ValueError('There must be an equal number of dependent '
                             'speeds and constraints')
        self._uds = speedl
        for i, v in enumerate(conl):
            conl[i] = v.subs(self._kd)
        self._dep_cons = conl

        uis = self._us[:]
        uds = self._uds[:]
        oth = self._find_others(conl, uis + uds + self._qs + self._qdots)
        n = len(uis)
        m = len(uds)
        p = n - m
        # puts independent speeds first
        for i, v in enumerate(uds):
            uis.remove(v)
        uis += uds
        B = zeros((m, n))
        C = zeros((m, len(oth)))

        ii = 0
        for i1, v1 in enumerate(uds):
            for i2, v2 in enumerate(uis):
                B[ii] = conl[i1].diff(uis[i2])
                ii += 1
        #REDO THIS PART BY SUBTRACTING THE ABOVE FROM conl[i]
        ii = 0
        for i1, v1 in enumerate(uds):
            for i2, v2 in enumerate(oth):
                C[ii] = conl[i1].diff(oth[i2])
                ii += 1
        for i, v in enumerate(uds):
            uis.remove(v)
        mr1 = B.extract(range(m), range(p))
        ml1 = B.extract(range(m), range(p, n))
        mr2 = C * Matrix(oth)
        ml1i = ml1.inverse_ADJ()
        self._Ars = ml1i * mr1
        self._Brs = ml1i * mr2

    def form_fr(self, fl):
        #dict is body, force
        if not isinstance(fl, (list, tuple)):
            raise TypeError('Forces must be supplied in a list of: lists or '
                            'tuples.')
        N = self._inertial
        self._forcelist = fl[:]
        uis = self._us
        uds = self._uds
        for i, v in enumerate(uds):
            uis.remove(v)
        uis += uds
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
        if not isinstance(bl, (list, tuple)):
            raise TypeError('Bodies must be supplied in a list.')
        N = self._inertial
        self._bodylist = bl[:]
        uis = self._us
        uds = self._uds
        for i, v in enumerate(uds):
            uis.remove(v)
        uis += uds
        n = len(uis)

        # Form R*, T* for each body or particle in the list
        rsts = []
        for i, v in enumerate(bl):
            if isinstance(v, RigidBody):
                om = v.frame.ang_vel_in(N).subs(self._kd)
                ve = v.mc.vel(N).subs(self._kd)
                if v.mass.diff(Symbol('t')) != 0:
                    r = (v.mass * ve).dt(N)
                else:
                    r = v.mass * v.mc.acc(N).subs(self._kd)
                I, p = v.inertia
                if p != v.mc:
                    pass
                    #redefine I
                if I.dt(v.frame) != 0:
                    t = (I & om).diff(Symbol('t'), N)
                else:
                    t = ((v.frame.ang_acc_in(N).subs(self._kd) & I) +
                         ((om ^ I) & om))
            elif isinstance(v, Particle):
                t = 0
                ve = v.point.vel(N).subs(self._kd)
                if v.mass.diff(Symbol('t')) != 0:
                    r = (v.mass * ve).dt(N)
                else:
                    r = w.mass * v.point.acc(N).subs(self._kd)
            else:
                raise TypeError('The body list needs RigidBody or '
                                'Particle as list elements')
            rsts.append((r, t))

        # Use R*, T* and partial velocities to form FR*
        FRSTAR = zeros((n, 1))
        # goes through each Fr (where this loop's i is r)
        for i, v in enumerate(uis):
            # does this for each body in the list (w)
            for j, w in enumerate(bl):
                if isinstance(w, RigidBody):
                    om = w.frame.ang_vel_in(N).subs(self._kd)
                    ve = w.mc.vel(N).subs(self._kd)
                    r = rsts[j][0] & ve.diff(v, N)
                    t = rsts[j][1] & om.diff(v, N)
                    FRSTAR[i] -= (r + t)
                elif isinstance(w, Particle):
                    ve = w.mc.vel(N).subs(self._kd)
                    FRSTAR[i] -= rsts[j][0] & ve.diff(v, N)
        if len(uds) != 0:
            m = len(uds)
            p = n - m
            FRSTARtilde = FRSTAR.extract(range(p), [0])
            FRSTARold = FRSTAR.extract(range(p, n), [0])
            FRSTARtilde += self._Ars.T * FRSTARold
            FRSTAR = FRSTARtilde
        self._frstar = FRSTAR
        return FRSTAR

    def mass_matrix(self):
        uis = self._us
        udots = self._udots
        uds = self._uds
        for i, v in enumerate(uds):
            uis.remove(v)
        p = len(uis)
        zeroeq = self._fr + self._frstar
        MM = zeros((p, p))
        # strange counter is because sympy matrix only has 1 index even if 2d
        ii = 0
        for i in range(p):
            for j in range(p):
                MM[ii] = zeroeq[i].diff(udots[j])
                ii += 1
        self._massmatrix = MM
        return MM

    def rhs(self):
        uis = self._us
        udots = self._udots
        uds = self._uds
        for i, v in enumerate(uds):
            uis.remove(v)
        p = len(uis)
        zeroeq = self._fr + self._frstar
        MM = self._massmatrix
        # strange counter is because sympy matrix only has 1 index even if 2d
        ii = 0
        for i in range(p):
            for j in range(p):
                zeroeq[i] -= MM[ii] * udots[j]
                ii += 1
            zeroeq[i] = -simplify(expand(zeroeq[i]))
        self._rhs = zeroeq
        return zeroeq

if __name__ == "__main__":
    import doctest
    doctest.testmod()
