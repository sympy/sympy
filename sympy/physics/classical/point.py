__all__ = ['Point']

from sympy.physics.classical.essential import Vector, ReferenceFrame

class Point(object):
    """This object represents a point in a dynamic system.
    It stores the: position, velocity, and acceleration of a point.
    The position is a vector defined as the vector distance from a parent
    point to this point.

    """

    def __init__(self, name):
        """Initialization of a Point object. """
        self.name = name
        self._pos_dict = {}
        self._vel_dict = {}
        self._acc_dict = {}
        self._dlist = [self._pos_dict, self._vel_dict, self._acc_dict]

    def _check_frame(self, other):
        if not isinstance(other, ReferenceFrame):
            raise TypeError('A ReferenceFrame must be supplied')

    def _check_point(self, other):
        if not isinstance(other, Point):
            raise TypeError('A Point must be supplied')

    def _check_vector(self, other):
        if isintance(other, int):
            if other == 0:
                return
        if not isinstance(other, Vector):
            raise TypeError('A Vector must be supplied')

    def _dict_list(self, other, num):
        """Creates a list from self to other using _dcm_dict. """
        outlist = [[self]]
        oldlist = [[]]
        while outlist != oldlist:
            oldlist = outlist[:]
            for i, v in enumerate(outlist):
                templist = v[-1]._dlist[num].keys()
                for i2, v2 in enumerate(templist):
                    if not v.__contains__(v2):
                        littletemplist = v + [v2]
                        if not outlist.__contains__(littletemplist):
                            outlist.append(littletemplist)
        for i, v in enumerate(oldlist):
            if v[-1] != other:
                outlist.remove(v)
        outlist.sort(key = len)
        if len(outlist) != 0:
            return outlist[0]
        raise ValueError('No Connecting Path Found')

    def set_pos(self, value, otherpoint):
        """Used to set the position of this point w.r.t. another point.

        Parameters
        ==========
        value : Vector
            The vector which defines the location of this point
        point : Point
            The other point which this point's location is defined relative to

        Examples
        ========

        >>> from sympy.physics.classical import Point, ReferenceFrame
        >>> N = ReferenceFrame('N')
        >>> p1 = Point('p1')
        >>> p2 = Point('p2')
        >>> p1.set_pos(10 * N.x, p2)
        >>> p1.pos(p2)
        (10)*nx>

        """

        self._check_vector(value)
        self._check_point(otherpoint)
        self._pos_dict.update({otherpoint: value})
        otherpoint._pos_dict.update({self: -value})

    def pos(self, otherpoint):
        """Returns a Vector distance between this Point and the other Point.

        Parameters
        ==========
        otherpoint : Point
            The otherpoint we are locating this one relative to

        Examples
        ========

        >>> from sympy.physics.classical import Point, ReferenceFrame
        >>> N = ReferenceFrame('N')
        >>> p1 = Point('p1')
        >>> p2 = Point('p2')
        >>> p1.set_pos(10 * N.x, p2)
        >>> p1.pos(p2)
        (10)*nx>

        """

        outvec = 0
        plist = self._dict_list(otherpoint, 0)
        for i in range(len(plist) - 1):
            outvec += plist[i]._pos_dict[plist[i + 1]]
        return outvec

# add dict search func

    def set_vel(self, value, frame):
        """Sets the velocity Vector of this Point in a ReferenceFrame.

        Parameters
        ==========
        value : Vector
            The vector value of this point's velocity in the frame
        frame : ReferenceFrame
            The frame in which this point's velocity is defined

        Examples
        ========

        >>> from sympy.physics.classical import Point, ReferenceFrame
        >>> N = ReferenceFrame('N')
        >>> p1 = Point('p1')
        >>> p1.set_vel(10 * N.x, N)
        >>> p1.vel(N)
        (10)*nx>

        """

        self._check_vector(value)
        self._check_frame(frame)
        self._vel_dict.update({frame: value})

    def vel(self, frame):
        """The velocity Vector of this Point in the ReferenceFrame.

        Parameters
        ==========
        frame : ReferenceFrame
            The frame in which the returned velocity vector will be defined in

        Examples
        ========

        >>> from sympy.physics.classical import Point, ReferenceFrame
        >>> N = ReferenceFrame('N')
        >>> p1 = Point('p1')
        >>> p1.set_vel(10 * N.x, N)
        >>> p1.vel(N)
        (10)*nx>

        """

        self._check_frame(frame)
        if not self._vel_dict.has_key(frame):
            raise ValueError('Velocity has not been defined in '
                             'that ReferenceFrame; define it first')
        return self._vel_dict[frame]

    def set_acc(self, value, frame):
        """Used to set the acceleration of this Point in a ReferenceFrame.

        Parameters
        ==========
        value : Vector
            The vector value of this point's acceleration in the frame
        frame : ReferenceFrame
            The frame in which this point's acceleration is defined

        Examples
        ========

        >>> from sympy.physics.classical import Point, ReferenceFrame
        >>> N = ReferenceFrame('N')
        >>> p1 = Point('p1')
        >>> p1.set_acc(10 * N.x, N)
        >>> p1.acc(N)
        (10)*nx>

        """

        self._check_vector(value)
        self._check_frame(frame)
        self._acc_dict.update({frame: value})

    def acc(self, frame):
        """The acceleration Vector of this Point in a ReferenceFrame.

        Parameters
        ==========
        frame : ReferenceFrame
            The frame in which the returned acceleration vector will be defined in

        Examples
        ========

        >>> from sympy.physics.classical import Point, ReferenceFrame
        >>> N = ReferenceFrame('N')
        >>> p1 = Point('p1')
        >>> p1.set_acc(10 * N.x, N)
        >>> p1.acc(N)
        (10)*nx>

        """

        self._check_frame(frame)
        if not self._acc_dict.has_key(frame):
            raise ValueError('Acceleration has not been defined in '
                             'that ReferenceFrame; define it first')
        return self._acc_dict[frame]
