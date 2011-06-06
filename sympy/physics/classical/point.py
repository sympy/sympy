__all__ = ['Point']

from sympy.physics.classical.essential import Vector, ReferenceFrame

class Point(object):
    """This object represents a point in a dynamic system.
    It stores the: position, velocity, and acceleration of a point.
    The position is a vector defined as the vector distance from a parent
    point to this point.

    """

    def __init__(self, name):
        """Initialization of a Point object.  Takes in a name, sets
        attributes to zero.

        """

        self.name = name
        self._pos = None
        self._pos_par = None
        self._vel = None
        self._vel_frame = None
        self._acc = None
        self._acc_frame = None

    def set_pos(self, value, point = None):
        """Used to set the position of this point w.r.t. another point.

        """

        if value != 0:
            if not isinstance(value, Vector):
                raise TypeError('Position is a Vector')
        if point != None:
            if not isinstance(point, Point):
                raise TypeError('Need to supply a parent point')
        self._pos = value
        self._pos_par = point

    def pos(self, otherpoint = None):
        """Returns a Vector distance between this Point and the other Point.

        If no other Point is given, the value of this Point's position is
        returned.

        """

        if type(otherpoint) == type(None):
            return self._pos
        common_pos_par = self._common_pos_par(otherpoint)
        leg1 = 0
        ptr = self
        while ptr != common_pos_par:
            leg1 += ptr._pos
            ptr = ptr._pos_par
        leg2 = 0
        ptr = 0
        while ptr != common_pos_par:
            leg2 -= ptr._pos
            ptr = ptr._pos_par
        return leg1 + leg2

    def _common_pos_par(self,other):
        """This returns the first common parent between two ReferenceFrames."""
        leg1 = [self]
        ptr = self
        while ptr._pos_par != None:
            ptr = ptr._pos_par
            leg1.append(ptr)
        leg2 = [other]
        ptr = other
        while ptr._pos_par != None:
            ptr = ptr._pos_par
            leg2.append(ptr)
        for i,v1 in enumerate(leg1):
            for j, v2 in enumerate(leg2):
                if v1 == v2:
                    return v1
        raise ValueError('No Common Position Parent')

    def set_vel(self, value, frame):
        """
        Used to set the velocity Vector of this Point in a ReferenceFrame.
        """
        if not isinstance(frame, ReferenceFrame):
            raise TypeError('Velocity is defined in a ReferenceFrame')
        if value != 0:
            if not isinstance(value, Vector):
                raise TypeError('Velocity is a Vector')
        self._vel_frame = frame
        self._vel = value

    def vel(self, frame):
        """Returns the velocity of this Point in the ReferenceFrame, as a Vector.

        """

        if not isinstance(frame, ReferenceFrame):
            raise TypeError('Velocity is described in a frame')
        if frame != self._vel_frame:
            raise ValueError('Velocity has not been defined in '
                             'that ReferenceFrame; redefine it first')
        return self._vel

    def set_acc(self, value, frame):
        """Used to set the acceleration of this Point in a ReferenceFrame.

        """

        if not isinstance(frame, ReferenceFrame):
            raise TypeError('Acceleration is defined in a ReferenceFrame')
        if value != 0:
            if not isinstance(value, Vector):
                raise TypeError('Acceleration is a Vector')
        self._acc_frame = frame
        self._acc = value

    def acc(self, frame):
        """
        Returns the acceleration of this Point in a ReferenceFrame, as a
        Vector.
        """
        if not isinstance(frame, ReferenceFrame):
            raise TypeError('Velocity is described in a frame')
        if frame != self._acc_frame:
            raise ValueError('Acceleration has not been defined in that '
                             'ReferenceFrame; redefine it first')
        return self._acc

