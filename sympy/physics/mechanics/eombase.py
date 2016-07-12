# Mock up of class to hold the equations of motion (not a parent class)
from sympy import Matrix, zeros


class EOMBase(object):
    def __init__(self, coordinates, speeds, forcing_vector, mass_matrix,
                 k_kqdot, k_ku, f_k, loads=(), bodies=()):

        # Dynamical differential equation
        self._mass_matrix = mass_matrix
        self._forcing = forcing_vector

        # Kinematical differential equation
        # Using the equation form defined in Kane's method documentation
        # k_kqdot * qdot + k_ku * u + f_k = 0
        self._k_kqdot = k_kqdot
        self._k_ku = k_ku
        self._f_k = f_k

        # Coordinates and speeds
        self._q = coordinates
        self._u = speeds

        # Default values
        self._bodylist = bodies
        self._forcelist = loads

    def eom_check(self):
        try:
            self._mass_matrix
        except AttributeError:
            raise ValueError('Need to compute the equations of motion first')

    @property
    def bodylist(self):
        return self._bodylist

    @property
    def forcelist(self):
        return self._forcelist

    @property
    def forcing(self):
        self.eom_check()
        return self._forcing

    @property
    def forcing_full(self):
        self.eom_check()
        f1 = self._k_ku * Matrix(self.u) + self._f_k
        return Matrix([f1, self.forcing])

    @property
    def mass_matrix(self):
        self.eom_check()
        return self._mass_matrix

    @property
    def mass_matrix_full(self):
        self.eom_check()
        o = len(self.u)
        n = len(self.q)
        return ((self._k_kqdot).row_join(zeros(n, o))).col_join((zeros(o,
                n)).row_join(self.mass_matrix))

    @property
    def q(self):
        return self._q

    @property
    def u(self):
        return self._u
