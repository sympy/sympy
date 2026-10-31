from sympy import pi, sympify, S
from sympy.physics.mechanics import Vector, Particle, _check_frame, \
     get_motion_params
from sympy.physics.em import gradient, scalar_potential, divergence, laplacian

#The electromagnetic 'k' constant = 1/(4*pi*eo)
k = Symbol('k')


class ParticleCharge(Particle):
    """
    Class to represent a charged particle in space

    Parameters
    ==========

    name : string
        The name for the ParticleCharge instance

    point : Point
        The Point instance corresponding to this ParticleCharge

    mass : sympyfiable
        The mass of the charged particle

    charge : sympyfiable
        The charge of the charged particle

    """

    def __init__(self, name, point, mass, charge=S(0)):
        super(ParticleCharge, self).__init__(name, point, mass)
        self._charge = charge

    @property
    def charge(self):
        """ Charge of this particle """
        return self._charge

    def set_motion(self, frame, **kwargs):
        """
        Set motion of this ParticleCharge according to given parameters,
        in the user-specified frame.

        See the docs of functions.get_motion_params for more info-

        """
        params = get_motion_params(frame, **kwargs)
        self._point.set_pos(frame, params[2])
        self._point.set_vel(frame, params[1])
        self._point.set_acc(frame, params[0])
        pos_vector = frame.express(params[2], variables=True)
        pos_vector = frame[0]*frame.x + frame[1]*frame.y + frame[2]*frame.z -\
                     pos_vector
        self._pot_func = k*self.charge / \
                        (frame.x.dot(pos_vector)**2 + \
                         frame.y.dot(pos_vector)**2 + \
                         frame.z.dot(pos_vector)**2)**0.5
        self._field = electrostatic_field(self._pot_func, frame)

    set_motion.__doc__ = get_motion_params.__doc__

    def electrostatic_potential(self, frame, point=None):
        """
        Scalar potential function of the particle in given field.
        If a Point is provided, the potential at that position is returned.

        Parameters
        ==========

        frame : ReferenceFrame
            The field to express the potential function in

        point : Point
            Point to calculate the potential at

        Examples
        ========
        
        """

        _check_frame(frame)
        if point is None:
            return frame.express(self._pot_func)
        else:
            point = frame.express(point, variables=True)
            subs_dict = {}
            for i, x in enumerate(frame):
                subs_dict[frame[i]] = x.dot(point)
            return frame.express(self._pot_func).subs(subs_dict)

    def electrostatic_field(self, frame, point=None):
        """
        Scalar potential function of the particle in given field
        If a point(position vector) is provided, the field at that
        point is returned.

        Parameters
        ==========

        frame : ReferenceFrame
            The field to express the potential function in

        point(optional) : Point
            Point to calculate the field at

        Examples
        ========
        
        """

        if self._frame is None:
            raise ValueError("Motion has not been set for Particle- "+\
                             str(self))
        _check_frame(frame)
        if point is None:
            return frame.express(self._field, variables=True)
        else:
            point = frame.express(point, variables=True)
            subs_dict = {}
            for i, x in enumerate(frame):
                subs_dict[frame[i]] = x.dot(point)
            return frame.express(self._field, variables=True).subs(subs_dict)

    def electrostatic_force(self, *efields):
        """
        The force experience by this particle under the influence of given
        electrostatic fields

        Parameters
        ==========

        efields : list(of Vectors)
            List of electrostatic fields to consider

        Examples
        ========

        """

        total_force = 0
        for x in efields:
            if not isinstance(x, Vector):
                raise TypeError(str(x)+ " is not a vector field")
            total_force = self.charge * x
        return total_force


def electrostatic_field(potential, frame):
    """
    The electric field corresponding to a scalar electrostatic potential \
    in given frame

    Parameters
    ==========

    potential : sympyfiable scalar
        The electrostatic potential function

    frame : ReferenceFrame
        The frame to do the calculations in

    Examples
    ========

    """
    return -1 * gradient(potential, frame)


def electrostatic_potential(electric_field, frame):
    """
    The electric scalar potential corresponding to an electrostatic field \
    in given frame

    Parameters
    ==========

    electric_field : Vector
        The electric field whose potential is to be calculated

    frame : ReferenceFrame
        The frame to do the calculations in

    Examples
    ========

    """
    return -1 * scalar_potential(electric_field, frame)


def charge_assembly_energy(*chargedparticles):
    """
    The work done to assemble a system of charged particles from infinity
    in free space.

    Parameters
    ==========

    chargedparticles : list(of ParticleCharge instances)
        The system whose energy is to be calculated

    Examples
    ========

    """

    total_work = 0
    for i in range(0, len(chargedparticles)-1):
        for j in range(i+1, len(chargedparticles)):
            total_work += chargedparticles[i].charge * \
                          chargedparticles[j].electrostatic_potential(
                              chargedparticles[j]._frame, \
                              chargedparticles[i].pos_vector_wrt(
                                  chargedparticles[j]))
    return total_work


def charge_density(field, frame):
    """
    The charge density wrt the given field or scalar potential in a frame.

    Parameters
    ==========

    field : Vector/sympyfiable scalar
        The vector electrostatic field or scalar potential field
        under consideration

    frame : ReferenceFrame
        The frame to do the calculations in

    Examples
    ========

    """

    _check_frame(frame)
    if isinstance(field, Vector):
        return divergence(field, frame) / (4*pi*k)
    else:
        field = sympify(field)
        return -1 * divergence(gradient(field, frame), frame) / (4*pi*k)
