# Mock up of class to hold the equations of motion (not a parent class)
from sympy import Matrix, zeros


class EOM(object):
    def __init__(self, coordinates, speeds, forcing_vector, mass_matrix,
                 k_kqdot, k_ku, f_k, loads=(), bodies=()):

        # Dynamical differential equation
        self._k_d = mass_matrix
        self._f_d = forcing_vector

        # Kinematical differential equation
        # Using the equation form defined in Kane's method documentation
        # k_kqdot * qdot + k_ku * u + f_k = 0
        self._k_kqdot = k_kqdot
        self._k_ku = k_ku
        self._f_k = f_k

        # Coordinates and speeds
        self._coordinates = coordinates
        self._speeds = speeds

        # Default values
        self._bodies = bodies
        self._loads = loads

    @property
    def bodies(self):
        return self._bodies

    @property
    def loads(self):
        return self._loads

    @property
    def forcing(self):
        return self._f_d

    @property
    def forcing_full(self):
        f1 = self._k_ku * Matrix(self.speeds) + self._f_k
        return Matrix([f1, self.forcing])

    @property
    def mass_matrix(self):
        return self._k_d

    @property
    def mass_matrix_full(self):
        o = len(self.speeds)
        n = len(self.coordinates)
        return ((self._k_kqdot).row_join(zeros(n, o))).col_join((zeros(o,
                n)).row_join(self.mass_matrix))

    @property
    def coordinates(self):
        return self._coordinates

    @property
    def speeds(self):
        return self._speeds
