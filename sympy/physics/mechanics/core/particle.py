from sympy.physics.mechanics import MovingRefFrame


class Particle(object):
    """
    A particle with non-zero mass and no spatial extension.

    Values need to be supplied at initialization only.
    """

    def __init__(self, name, pos_vector=None, trans_vel=None, trans_acc=None, \
                 mass=0, pot_energy=0, frame=None, **kwargs):
        """
        Initializer for Particle class
        """
        self._mass = sympify(mass)
        if frame is None or not isinstance(frame, MovingRefFrame):
            raise ValueError("Valid frame needs to be entered")
        self._frame = MovingRefFrame(name + "_frame", pos_vector, trans_vel, \
                                     trans_acc, parentframe = frame, **kwargs)
        self._pe = pot_energy

    def __str__(self):
        return name

    @property
    def frame(self):
        """
        Returns frame associated with this Particle
        """
        return self._frame

    @property
    def mass(self):
        """
        Returns value of particle's mass
        """
        return self._mass

    def locatenew(self, name, pos_vect=0, mass=0, pe=0):
        """
        Initialize a new Particle at a certain position wrt this Particle.

        Parameters
        ==========

        name : String
            Name of the new particle

        pos_vect : vector
            The position vector of the new particle wrt this one

        mass : sympifiable
            The mass of the new particle

        pe : sympifiable
            The potential energy of the new particle

        Examples
        ========

        """

        return Particle(name, pos_vector = pos_vect, mass = mass,
                        pot_energy = pe, frame = self.frame)

    def trans_vel_wrt(self, other):
        """
        The translational velocity of this particle wrt a specified frame
        or Particle.

        Parameters
        ==========

        other : MovingRefFrame/Particle
            The frame/particle wrt which this particle's translational vel is
            to be calculated

        Examples
        ========

        """

        if isinstance(other, MovingRefFrame):
            return self.frame.trans_vel_in(other)
        elif isinstance(other, Particle):
            return self.frame.trans_vel_in(other.frame)
        else:
            raise TypeError("Wrong type of argument - " + \
                            str(type(other)))

    def trans_acc_wrt(self, frame):
        """
        The translational acceleration of this particle wrt a specified frame
        or Particle.

        Parameters
        ==========

        other : MovingRefFrame/Particle
            The frame/particle wrt which this particle's translational accln is
            to be calculated

        Examples
        ========

        """

        if isinstance(other, MovingRefFrame):
            return self.frame.trans_acc_in(other)
        elif isinstance(other, Particle):
            return self.frame.trans_acc_in(other.frame)
        else:
            raise TypeError("Wrong type of argument - " + \
                            str(type(other)))

    def linear_momentum(self, frame):
        """
        Linear momentum of the particle in the specified frame.

        The linear momentum L, of a particle P, with respect to frame N is
        given by

        L = m * v

        where m is the mass of the particle, and v is the velocity of the
        particle in the frame N.

        Parameters
        ==========

        frame : MovingRefFrame
            The frame to express the linear momentum of this particle in

        Examples
        ========
        
        """
        return self.mass * self.trans_vel_in(frame)

    def angular_momentum(self, pos_vector, frame):
        """
        Angular momentum of the particle wrt a certain point in the specified
        frame.

        The angular momentum H, about some point O of a particle, P, is given
        by:

        H = r x m * v

        where r is the position vector from point O to the particle P, m is
        the mass of the particle, and v is the velocity of the particle in
        the inertial frame, N.

        Parameters
        ==========

        pos_vector : vector
            Position vector of the point wrt which the angular momentum is to be
            calculated

        frame : MovingRefFrame
            The frame in which the point's position vector is defined

        Examples
        ========

        """
        
        pos_vector = self.frame.convert_pos_vector(pos_vector, frame)
        return (-1 * pos_vector) ^ (self.mass * self.trans_vel_in(frame))

    def kinetic_energy(self, frame):
        """Kinetic energy of the particle in the specified frame

        The kinetic energy, T, of a particle, P, is given by

        'T = 1/2 m v^2'

        where m is the mass of particle P, and v is the velocity of the
        particle in the supplied ReferenceFrame.

        Parameters
        ==========

        frame : MovingRefFrame
            The Particle's velocity is typically defined with respect to
            an inertial frame but any relevant frame in which the velocity is
            known can be supplied.

        Examples
        ========

        """

        return (self.mass / sympify(2) * self.trans_vel_in(frame) &
                self.trans_vel_in(frame))

    @property
    def potential_energy(self):
        """
        The potential energy of this particle
        """
        return self._pe
