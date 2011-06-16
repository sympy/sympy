__all__ = ['Kane']

from sympy import Symbol, zeros, simplify
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
        self._uds = None
        self._qs = None
        self._qds = None
        self._kd = dict()
        self._time_varying = []
        self._forcelist = None
        self._bodylist = None
        self._fr = None
        self._frstar = None
        self._udeps = []
        self._dep_cons = []

    def _find_others(self, inlist, insyms):
        """Finds all non-supplied DynamicSymbols in the list of expressions"""

        def _deeper(inexpr):
            oli = []
            try:
                for i, v in enumerate(inexpr.args):
                    oli.append(_deeper(v))
            except:
                if isinstance(inexpr, DynamicSymbol):
                    oli.append(inexpr)
            return oli

        ol = []
        for i, v in enumerate(inlist):
            try:
                for i2, v2 in enumerate(v.args):
                    ol.append(_deeper(v2))
            except:
                if isinstance(v, DynamicSymbol):
                    ol.append(inexpr)
        seta = {} 
        map(seta.__setitem__, ol, []) 
        ol = seta.keys()
        for i, v in enumerate(insyms):
            ol.remove(v)
        return ol

    def gen_speeds(self, inlist):
        if not isinstance(inlist, (list, tuple)):
            raise TypeError('Generalized Speeds must be supplied in a list')
        self._us = inlist
        ol = []
        for i, v in enumerate(inlist):
            ol.append(v.diff(Symbol('t')))
        self._uds = ol

    def kindiffeq(self, indict):
        if not isinstance(indict, dict):
            raise TypeError('Kinematic Differential Equations must be '
                            'supplied in a dictionary')
        self._kd = indict
        newlist = indict.keys()
        for i, v in enumerate(newlist):
            newlist[i] = DynamicSymbol(v.name[:-1])
        self._qs = newlist
        self._qds = indict.keys()

    def dependent_speeds(self, speedl, conl):
        if not isinstance(speedl, (list, tuple)):
            raise TypeError('Dependent speeds and constraints must each be '
                            'provided in their own list')
        if not speedl.__len__() == conl.__len__():
            raise ValueError('There must be an equal number of dependent '
                             'speeds and constraints')
        self._udeps = speedl
        for i, v in enumerate(conl):
            conl[i] = v.subs(self._kd)
        self._dep_cons = conl

        uis = self._us[:]
        uds = self._udeps[:]
        oth = self._find_others(conl, uis + uds + self._qs + self._qds)
        n = len(uis)
        m = len(uds)
        p = n - m
        for i, v in enumerate(uds):
            uis.remove(v)
        uis += udsv
        B = zeros((m, n))
        C = zeros((m, len(oth)))

        ii = 0
        for i1, v1 in enumerate(uds):
            for i2, v2 in enumerate(uis):
                A[ii] = conl[i1].coeff(uis[i2])
                ii += 1
        #REDO THIS PART BY SUBTRACTING THE ABOVE FROM conl[i]
        ii = 0
        for i1, v1 in enumerate(uds):
            for i2, v2 in enumerate(oth):
                B[ii] = conl[i1].coeff(oth[i2])
                ii += 1
        for i, v in enumerate(uds):
            uis.remove(v)
        mr1 = B.extract([range(m)], [range(p)])
        ml1 = B.extract([range(m)], [range(p, n)])
        mr2 = C * Matrix(oth)
        self._Ars = ml1.inv() * mr1
        self._Brs = ml1.inv() * mr2

    def form_fr(self, fl):
        #dict is body, force
        if not isinstance(fl, (list, tuple)):
            raise TypeError('Forces must be supplied in a list of: lists or '
                            'tuples.')
        N = self._inertial
        self._forcelist = fl[:]
        uis = self._us
        uds = self._udeps
        for i, v in enumerate(uds):
            uis.remove(v)
        uis += uds
        n = len(uis)

        FR = zeros((n, 1))
        for i, v in enumerate(uis):
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
        if len(uds) != 0:
            m = len(uds)
            p = n - m
            FRtilde = FR.extract([range(p), 0])
            FRold = FR.extract([range(p, n), 0])
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
        uds = self._udeps
        for i, v in enumerate(uds):
            uis.remove(v)
        uis += uds
        n = len(uis)

        FRSTAR = zeros((n, 1))
        for i, v in enumerate(uis):
            for j, w in enumerate(bl):
                if isinstance(w, RigidBody):
                    om = w.frame.ang_vel_in(N).subs(self._kd)
                    ve = w.mc.vel(N).subs(self._kd)
                    if w.mass.diff(Symbol('t')) != 0:
                        r = (w.mass * ve).dt(N)
                    else:
                        r = w.mass * ve
                    r = r & ve.diff(v, N)
                    I, p = w.inertia
                    if p != w.mc:
                        pass
                        #redefine I
                    if I.dt(w.frame) != 0:
                        t = (I & om).diff(Symbol('t'), N)
                    else:
                        t = ((w.frame.ang_acc_in(N).subs(self._kd) & I) +
                             ((om ^ I) & om))
                    t = t & om.diff(v, N)
                    FRSTAR[i] -= (r + t)
                elif isinstance(w, Particle):
                    traspeed = w.mc.vel(N).subs(self._kd).diff(v)
                    r = w.mass * w.mc.vel(N).subs(self._kd)
                    r = r.diff(Symbol('t'))
                    r = r & traspeed
                    FRSTAR[i] -= traspeed & r
                else:
                    raise TypeError('The body list needs RigidBody or '
                                    'Particle as list elements')
        if len(uds) != 0:
            m = len(uds)
            p = n - m
            FRSTARtilde = FRSTAR.extract([range(p), 0])
            FRSTARold = FRSTAR.extract([range(p, n), 0])
            FRSTARtilde += self._Ars.T * FRSTARold
            FRSTAR = FRSTARtilde
        self._frstar = FRSTAR
        return FRSTAR

    def mass_matrix(self):
        uis = self._us
        udots = self._uds
        uds = self._udeps
        for i, v in enumerate(uds):
            uis.remove(v)
        p = len(uis)
        zeroeq = self._fr + self._frstar
        MM = zeros((p, p))
        ii = 0
        for i in range(p):
            for j in range(p):
                MM[ii] = zeroeq[i].diff(udots[j])
                ii += 1
        self._massmatrix = MM
        return MM

    def rhs(self):
        uis = self._us
        udots = self._uds
        uds = self._udeps
        for i, v in enumerate(uds):
            uis.remove(v)
        p = len(uis)
        zeroeq = self._fr + self._frstar
        MM = self._massmatrix
        ii = 0
        for i in range(p):
            for j in range(p):
                zeroeq[i] -= MM[ii] * udots[j]
                ii += 1
            zeroeq[i] = -simplify(simplify(zeroeq[i]))
        self._rhs = zeroeq
        return zeroeq

if __name__ == "__main__":
    import doctest
    doctest.testmod()
