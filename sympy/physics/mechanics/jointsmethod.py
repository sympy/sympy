from sympy.physics.mechanics import KanesMethod

__all__ = ['JointsMethod']


class JointsMethod(object):

    def __init__(self, root_body, *joints):
        self._joints = joints
        self.root_body = root_body
        self._bodylist = self._generate_bodylist()
        self._loadlist = self._generate_loadlist()
        self._q = self._generate_q()
        self._u = self._generate_u()
        self._kdes = self._generate_kdes()
        self._kane = KanesMethod(self.root_body.frame, q_ind=self.q,
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
        return self._kane.rhs(inv_method)

    def kanes_equations(self):
        if self.loadlist == []:
           self._fr, self._frstar = self._kane.kanes_equations(self.bodylist)
           return self.fr, self.frstar
        self._fr, self._frstar = self._kane.kanes_equations(self.bodylist, self.loadlist)
        return self.fr, self.frstar
