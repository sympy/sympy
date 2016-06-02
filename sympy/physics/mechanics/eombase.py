# Mock up of class to hold the equations of motion (not a parent class)


class EOM(object):
    def __init__(self, coordinates, speeds, forcing_vector, mass_matrix,
                 loads=(), bodies=()):
        self._k_d = mass_matrix
        self._f_d = forcing_vector
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
        pass

    @property
    def mass_matrix(self):
        return self._k_d

    @property
    def mass_matrix_full(self):
        pass

    @property
    def coordinates(self):
        return self._coordinates

    @property
    def speeds(self):
        return self._speeds
