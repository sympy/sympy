__all__ = ['EField', 'ParticleCharge']

from sympy.physics.mechanics import Vector, Particle
from sympy.physics.mechanics.essential import _check_vector
from sympy import integrate
from math import pi

e0 = 8.854187817620 * 10 ** -12


class EField:
    """
    Represents a purely electric vector field in space. It must be defined in
    such a way that it can be expressed entirely in one frame.
    """
    
    def __init__(self, vector):
        """
        Constructor for EField class.

        Parameters
        ==========

        vector : Vector

        A Vector defining the EField. The specified Vector must be a function
        of x ,y and z only. 'vector' must be conservative in nature.

        Examples
        ========

        >>> from sympy.physics.electrical import EField
        >>> from sympy.physics.mechanics import Vector, ReferenceFrame
        >>> from sympy import symbols
        >>> N = ReferenceFrame('N')
        >>> x, y, v = symbols('x y v')
        >>> v = - y /(x ** 2 + y ** 2) * N.x + x /(x ** 2 + y ** 2) * N.y
        >>> field = EField(v)
        >>> print field
        - y/(x**2 + y**2)*N.x + x/(x**2 + y**2)*N.y

        """
        self.set_value(vector)

    def get_vector(self):
        """Returns the Vector corresponding to this EField."""
        return self.vector

    def set_value(self, vector):
        """Set the vector value of an EField."""
        _check_vector(vector)
        #Check whether the field can be expressed entirely in one of the frames
        #it is defined in.
        vector.express(vector.args[0][1])
        #Check whether field is conservative.
        if (_check_conservative(vector)):
            self.vector = vector
            #Store a reference frame to do all the relevant calculations in.
            self.frame = vector.args[0][1]
        else:
            raise ValueError("Given vector must define conservative field.")

    def __add__(self, other_field):
        """Addition operator for EField."""
        return EField(self.vector + other_field.get_vector())

    def __sub__(self, other_field):
        """Subtraction operator for EField."""
        return EField(self.vector - other_field.get_vector())

    def __mul__(self, expr):
        """
        Returns the EField whose Vector is the product of this field's
        Vector and expr.
        'expr' must be a Sympifiable expression.
        """
        return EField(expr * self.vector)

    def __div__(self, expr):
        """Division operator for EField."""
        return EField(self.vector / expr)

    def __eq__(self, other):
        """To test the equality of two EFields."""
        return self.vector == other.get_vector()

    def __str__(self):
        """Printing methid."""
        return (self.vector).__str__()
    
    def potential_difference(self, origin, point1, point2):
        """
        Calculates the electrostatic potential difference between two points for the field
        and the specified origin. The position vectors of the points with respect to the
        origin must be such that they can be expressed entirely in the reference frames
        contained in the definition of the electric field.

        The potential difference is irrespective of the reference frame used for calculations.

        Examples
        ========
        >>> from sympy.physics.electrical import EField
        >>> from sympy.physics.mechanics import Vector, ReferenceFrame, Point
        >>> from sympy import Symbol
        >>> N = ReferenceFrame('N')
        >>> x = Symbol('x')
        >>> field = EField(3*N.x)
        >>> o = Point('o')
        >>> p = Point('p')
        >>> p.set_pos(o,3*N.x+4*N.y)
        >>> q = Point('q')
        >>> q.set_pos(o,2*N.x+7*N.y)
        >>> field.potential_difference(o,p,q)
        3
        
        """
        from sympy.abc import x, y, z
        variables = [x, y, z]
        #Express every vector in one frame.
        pos_vector1 = (point1.pos_from(origin)).express(self.frame)
        pos_vector2 = (point2.pos_from(origin)).express(self.frame)
        temp = (self.vector).express(self.frame)
        initial_values = []
        final_values = []
        #Store limits of integration in initial_values and final_values
        for i in range(3):
            initial_values.append(pos_vector1.args[0][0][i])
            final_values.append(pos_vector2.args[0][0][i])
        pot_difference = 0
        #Integrate electric field to get potential difference.
        for i in range(3):
            function = integrate(temp.args[0][0][i], variables[i])
            pot_difference = pot_difference - function.subs(zip(variables, final_values)) + function.subs(zip(variables, initial_values))
        return pot_difference


class ParticleCharge(Particle):
    """
    Represents a charged particle of zero volume. May have non-zero mass and
    electric charge.
    """

    def __init__(self, name, point, charge = 0, mass = 0):
        Particle.__init__(self, name, point, mass)
        self.set_charge(charge)

    def set_charge(self, newcharge):
        self.charge = newcharge

    def get_charge(self):
        return self.charge

    def field_at(self, point):
        pos_vector = point.pos_from(self.get_point())
        magnitude = (self.charge)/(4 * pi * e0 * (pos_vector.magnitude()) ** 2)
        return EField(magnitude * (pos_vector.normalize()))

    def potential_at(self,vpoint):
        pos_vector = point.pos_from(self.get_point())
        return (self.charge)/(4 * pi * e0 * abs(pos_vector.magnitude()))


def _check_conservative(vector):
    """
    Checks if 'vector' defines a field that is conservative.
    If yes, return True.
    Else, return False.
    """

    from sympy import symbols
    x, y, z = symbols('x y z')
    dfdx = ((vector.args[0][0][0]).diff(y)).diff(z)
    dfdy = ((vector.args[0][0][1]).diff(x)).diff(z)
    dfdz = ((vector.args[0][0][2]).diff(x)).diff(y)
    if (dfdx == dfdy == dfdz):
        return True
    else:
        return False
