__all__ = ['Point']

from sympy.physics.classical.essential import Vector, ReferenceFrame

class Point(object):
    """
    This object represents a point in a dynamic system.
    It stores the: position, velocity, and acceleration of a point.
    The position is a vector defined as the vector distance from a parent
    point to this point.
    """

    def __init__(self, name):
        """
        Initialization of a Point object.  Takes in a name, sets
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
        """
        Used to set the position of this point with respect to another point.
        """
        if value != 0:
            assert isinstance(value, Vector), ('Need to supply a Vector for ',
                                               'the position of this Point')
        if point != None:
            assert isinstance(point, Point), 'Need to supply a parent point'
        self._pos = value
        self._pos_par = point

    def pos(self, otherpoint = None):
        """
        Returns a Vector distance between this Point and the other Point.
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
        """
        This returns the first common parent between two ReferenceFrames.
        Takes in another ReferenceFrame, and returns a ReferenceFrame.
        """
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
        assert isinstance(frame, ReferenceFrame), 'Velocity is defined \
                in a specific ReferenceFrame'
        if value != 0:
            assert isinstance(value, Vector), 'Velocity is a Vector'
        self._vel_frame = frame
        self._vel = value

    def vel(self, frame):
        """
        Returns the velocity of this Point in the ReferenceFrame, as a Vector.
        """
        assert isinstance(frame, ReferenceFrame), 'Velocity is described \
                in a frame'
        assert frame == self._vel_frame, 'Velocity has not been defined in \
                that ReferenceFrame; redefine it first'
        return self._vel

    def set_acc(self, value, frame):
        """
        Used to set the acceleration of this Point in a ReferenceFrame.
        """
        assert isinstance(frame, ReferenceFrame), 'Acceleration is defined \
                in a specific ReferenceFrame'
        if value != 0:
            assert isinstance(value, Vector), 'Acceleration is a Vector'
        self._acc_frame = frame
        self._acc = value

    def acc(self, frame):
        """
        Returns the acceleration of this Point in a ReferenceFrame, as a
        Vector.
        """
        assert isinstance(frame, ReferenceFrame), 'Velocity is described \
                in a frame'
        assert frame == self._acc_frame, 'Acceleration has not been defined \
                in that ReferenceFrame; redefine it first'
        return self._acc

