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
        self._pdlist = [self._pos_dict, self._vel_dict, self._acc_dict]

    def __str__(self):
        return self.name

    __repr__ = __str__

    def _check_frame(self, other):
        if not isinstance(other, ReferenceFrame):
            raise TypeError('A ReferenceFrame must be supplied')

    def _check_point(self, other):
        if not isinstance(other, Point):
            raise TypeError('A Point must be supplied')

    def _check_vector(self, other):
        if isinstance(other, int):
            if other == 0:
                return
        if not isinstance(other, Vector):
            raise TypeError('A Vector must be supplied')

    def _pdict_list(self, other, num):
        """Creates a list from self to other using _dcm_dict. """
        outlist = [[self]]
        oldlist = [[]]
        while outlist != oldlist:
            oldlist = outlist[:]
            for i, v in enumerate(outlist):
                templist = v[-1]._pdlist[num].keys()
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
        raise ValueError('No Connecting Path found between ' + other.name +
                         ' and ' + self.name)

    def a1pt(self, otherpoint, outframe, interframe):
        """Sets the acceleration of this point with the 1-point theory.

        The 1-point theory for point acceleration looks like this:

        :math:`^{N} \vec{a} ^{P} = ^{B} \vec{a} ^{P} + ^{N} \vec{a} ^{O} +
        ^{N} \vec{\alpha} ^{B} \times \vec{r} ^{OP} + ^{N} \vec{\omega} ^{B}
        \times (^{N} \vec{\omega} ^{B} \times \vec{r} ^{OP}) + 2 ^{N}
        \vec{\omega}^{B} \times ^{B} \vec{v} ^{P} `

        where Q is a point fixed in B, P is a point moving in B, and B is
        rotating in frame N.

        Parameters
        ==========
        otherpoint : Point
            The first point of the 1-point theory
        outframe : ReferenceFrame
            The frame we want this point's acceleration defined in
        fixedframe : ReferenceFrame
            The intermediate frame in this calculation

        Examples
        ========
        >>> from sympy.physics.classical import Point, ReferenceFrame,\
                DynamicSymbol
        >>> q = DynamicSymbol('q')
        >>> q2 = DynamicSymbol('q2')
        >>> qd = DynamicSymbol('qd')
        >>> q2d = DynamicSymbol('q2d')
        >>> N = ReferenceFrame('N')
        >>> B = ReferenceFrame('B')
        >>> B.set_ang_vel(N, 5 * B.y)
        >>> O = Point('O')
        >>> P = O.newpoint('P', q * B.x)
        >>> P.set_vel(qd * B.x + q2d * B.y, B)
        >>> O.set_vel(0, N)
        >>> P.a1pt(O, N, B)
        (-25*q + qdd)*bx> + (q2dd)*by> + (-10*qd)*bz>

        """

        self._check_frame(outframe)
        self._check_frame(interframe)
        self._check_point(otherpoint)
        dist = self.pos_from(otherpoint)
        v = self.vel(interframe)
        a1 = otherpoint.acc(outframe)
        a2 = self.acc(interframe)
        omega = interframe.ang_vel_in(outframe)
        alpha = interframe.ang_acc_in(outframe)
        self.set_acc(a2 + 2 * (omega ^ v) + a1 + (alpha ^ dist) + (omega ^
                     (omega ^ dist)), outframe)
        return self.acc(outframe)

    def a2pt(self, otherpoint, outframe, fixedframe):
        """Sets the acceleration of this point with the 2-point theory.

        The 2-point theory for point acceleration looks like this:

        :math:`^{N} \vec{a} ^{P} = ^{N} \vec{a} ^{O} + ^{N} \vec{\alpha} ^{B}
         \times \vec{r} ^{OP} + ^{N} \vec{\omega} ^{B} \times (
         ^{N} \vec{\omega} ^{B} \times \vec{r} ^{OP})`

        where Q and P are both points fixed in frame B, which is rotating in
        frame N.

        Parameters
        ==========
        otherpoint : Point
            The first point of the 2-point theory
        outframe : ReferenceFrame
            The frame we want this point's acceleration defined in
        fixedframe : ReferenceFrame
            The frame in which both this point and the other point are fixed

        Examples
        ========

        >>> from sympy.physics.classical import Point, ReferenceFrame,\
                DynamicSymbol
        >>> q = DynamicSymbol('q')
        >>> qd = DynamicSymbol('qd')
        >>> N = ReferenceFrame('N')
        >>> B = N.orientnew('B', 'Simple', q, 3)
        >>> O = Point('O')
        >>> P = O.newpoint('P', 10 * B.x)
        >>> O.set_vel(5 * N.x, N)
        >>> P.a2pt(O, N, B)
        (-10*qd**2)*bx> + (10*qdd)*by>

        """

        self._check_frame(outframe)
        self._check_frame(fixedframe)
        self._check_point(otherpoint)
        dist = self.pos_from(otherpoint)
        v = otherpoint.vel(outframe)
        a = otherpoint.acc(outframe)
        omega = fixedframe.ang_vel_in(outframe)
        alpha = fixedframe.ang_acc_in(outframe)
        self.set_acc(a + (alpha ^ dist) + (omega ^ (omega ^ dist)), outframe)
        return self.acc(outframe)

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
            if self._vel_dict[frame] != 0:
                return (self._vel_dict[frame]).dt(frame)
            else:
                return 0
        return self._acc_dict[frame]

    def newpoint(self, name, value):
        """Creates a new point with a position defined from this point.

        Parameters
        ==========
        name : str
            The name for the new point
        value : Vector
            The position of the new point relative to this point

        Examples
        ========

        >>> from sympy.physics.classical import ReferenceFrame, Point
        >>> N = ReferenceFrame('N')
        >>> P1 = Point('P1')
        >>> P2 = P1.newpoint('P2', 10 * N.x)

        """

        if not isinstance(name, str):
            raise TypeError('Must supply a valid name')
        self._check_vector(value)
        p = Point(name)
        p.set_pos(-value, self)
        self.set_pos(-value, p)
        return p

    def pos_from(self, otherpoint):
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
        >>> p1.pos_from(p2)
        (10)*nx>

        """

        outvec = 0
        plist = self._pdict_list(otherpoint, 0)
        for i in range(len(plist) - 1):
            outvec += plist[i]._pos_dict[plist[i + 1]]
        return outvec

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
        >>> p1.pos_from(p2)
        (10)*nx>

        """

        self._check_vector(value)
        self._check_point(otherpoint)
        self._pos_dict.update({otherpoint: value})
        otherpoint._pos_dict.update({self: -value})

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

    def v1pt(self, otherpoint, outframe, interframe):
        """Sets the velocity of this point with the 1-point theory.

        The 1-point theory for point velocity looks like this:

        :math:`^{N} \vec{v} ^{P} = ^{B} \vec{v} ^{P} + ^{N} \vec{v} ^{O} +
        ^{N} \vec{\omega} ^{B} \times \vec{r} ^{OP}`

        where Q is a point fixed in B, P is a point moving in B, and B is
        rotating in frame N.

        Parameters
        ==========
        otherpoint : Point
            The first point of the 2-point theory
        outframe : ReferenceFrame
            The frame we want this point's velocity defined in
        interframe : ReferenceFrame
            The intermediate frame in this calculation

        Examples
        ========

        >>> from sympy.physics.classical import Point, ReferenceFrame,\
                DynamicSymbol
        >>> q = DynamicSymbol('q')
        >>> q2 = DynamicSymbol('q2')
        >>> qd = DynamicSymbol('qd')
        >>> q2d = DynamicSymbol('q2d')
        >>> N = ReferenceFrame('N')
        >>> B = ReferenceFrame('B')
        >>> B.set_ang_vel(N, 5 * B.y)
        >>> O = Point('O')
        >>> P = O.newpoint('P', q * B.x)
        >>> P.set_vel(qd * B.x + q2d * B.y, B)
        >>> O.set_vel(0, N)
        >>> P.v1pt(O, N, B)
        (qd)*bx> + (q2d)*by> + (-5*q)*bz>

        """

        self._check_frame(outframe)
        self._check_frame(interframe)
        self._check_point(otherpoint)
        dist = self.pos_from(otherpoint)
        v1 = self.vel(interframe)
        v2 = otherpoint.vel(outframe)
        omega = interframe.ang_vel_in(outframe)
        self.set_vel(v1 + v2 + (omega ^ dist), outframe)
        return self.vel(outframe)

    def v2pt(self, otherpoint, outframe, fixedframe):
        """Sets the velocity of this point with the 2-point theory.

        The 2-point theory for point velocity looks like this:

        :math:`^{N} \vec{v} ^{P} = ^{N} \vec{v} ^{O} + ^{N} \vec{\omega} ^{B}
         \times \vec{r} ^{OP}`

        where Q and P are both points fixed in frame B, which is rotating in
        frame N.

        Parameters
        ==========
        otherpoint : Point
            The first point of the 2-point theory
        outframe : ReferenceFrame
            The frame we want this point's velocity defined in
        fixedframe : ReferenceFrame
            The frame in which both this point and the other point are fixed

        Examples
        ========

        >>> from sympy.physics.classical import Point, ReferenceFrame,\
                DynamicSymbol
        >>> q = DynamicSymbol('q')
        >>> qd = DynamicSymbol('qd')
        >>> N = ReferenceFrame('N')
        >>> B = N.orientnew('B', 'Simple', q, 3)
        >>> O = Point('O')
        >>> P = O.newpoint('P', 10 * B.x)
        >>> O.set_vel(5 * N.x, N)
        >>> P.v2pt(O, N, B)
        (5)*nx> + (10*qd)*by>

        """

        self._check_frame(outframe)
        self._check_frame(fixedframe)
        self._check_point(otherpoint)
        dist = self.pos_from(otherpoint)
        v = otherpoint.vel(outframe)
        omega = fixedframe.ang_vel_in(outframe)
        self.set_vel(v + (omega ^ dist), outframe)
        return self.vel(outframe)

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
            raise ValueError('Velocity of point ' + self.name + ' has not been'
                             ' defined in ReferenceFrame ' + frame.name)
        return self._vel_dict[frame]

if __name__ == "__main__":
    import doctest
    doctest.testmod()
