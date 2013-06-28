from sympy import Symbol, diff, sympify, trigsimp
from sympy.core.cache import cacheit
from sympy.vector import CoordSys, Vector, VectAdd, VectMul, BaseScalar
from sympy.physics.mechanics.core import get_motion_acc, get_motion_vel, get_motion_pos


class MovingRefFrame(CoordSysRect):
    """
    A moving frame of reference in classical mechanics.

    It subclasses the static CoordSys class of vector module
    and adds motion-related functionality to it.

    This frame can have translational and angular motion with
    respect to another frame, which affects the expression of
    vectors/points defined in this frame in other frames, and
    vice versa.
    """

    def __init__(name, pos_vector=None, trans_vel=None, trans_acc=None, orient_type=None,
                 orient_amount=None, orient_order=None, ang_vel=None, ang_acc=None, parentframe=None, **kwargs):
        """
        Initializer for the MovingRefFrame class.

        For a class whose parentframe == None, the 'timevar' keyword arg must be given to initialize a Symbol for
        time. If not specified, it falls back to default Symbol('t').

        Translational motion can be specified either by specifying the position vector(pos_vector), velocity
        (trans_vel) or acceleration(trans_acc) as functions of time. Required boundary conditions can be specified
        as keyword arguments. The relevant keyword-args are -
        
        For velocity- 'pos_vector', 't' - the position at the value of time t
        For acceleration - 'trans_vel', 'pos_vector', 't1', 't2' - the velocity and position boundary conditions at
        times t1 and t2 respectively

        If more than one arguments are entered, preference is given in order- position, velocity, acceleration.
        Hence, entering both velocity and acceleration will result in motion being set as per the velocity conditions.

        Rotational motion can be specified either by specifying the orientation(refer CoordSys docs) or the angular
        velocity(ang_vel) or the angular acceleration(ang_acc). Boundary conditions are required, similar to the
        translational motion params.

        For angular velocity- 'rotation', 'rt' - the rotation at the value of time rt
        For angular acceleration - 'ang_vel', 'rotation', 'rt1', 'rt2' - the velocity and position
        boundary conditions at times rt1 and rt2 respectively

        If more than one arguments are entered, preference is given in order- orientation, velocity, acceleration.
        Hence, entering both angular velocity and acceleration will result in motion being set as per the velocity
        conditions.

        Any of the time-dependent arguments(apart from pos_vector and orientation) can be specified as a vector,
        or as a tuple of the form - ([v1, v2, v3], V) where V is a vector defined wrt parent frame, and v1, v2,
        v3 are components of a vector V' defined wrt this frame itself. The total vectorial value of the param
        would be V + V'.

        Boundary conditions must be expressed entirely in the parentframe. If one or more of the boundary condition
        parameters are not entered, they are taken to be zero automatically.
        """

        #Special case of no parent frame
        if parentframe is None:
            self.parent = None
            self._root = self
            if 'timevar' not in kwargs:
                #Default time variable
                self.time = Symbol('t')
            else:
                #Set global time
                kwargs['timevar'] = sympify(kwargs['timevar'])
                if not kwargs['timevar'].is_Symbol:
                    raise ValueError("timevar must be a Symbol")
                self.time = kwargs['timevar']
            #All motion params zero
            self._pos_vector = 0
            self._trans_vel = 0
            self._trans_acc = 0
            self._ang_vel = 0
            self._ang_acc = 0
            super(MovingRefFrame, self).__init__(name, dim=3, wrt = None)
        else:
            self.parent = parentframe
            self._root = parentframe._root
            #Take time variable from parent frame
            self.time = parentframe.time
            #Initial orientation and positioning in case of tuple args
            flag = False
            for x in (trans_vel, trans_acc, ang_vel, ang_acc):
                if type(x) == tuple:
                    flag = True
                    break
            if flag:
                for x in ('t', 'rt', 't1', 't2', 'rt1', 'rt2'):
                    if x in kwargs:
                        if kwargs[x] != 0:
                            raise ValueError("initial values(t=0) needed for initialization with tuple args")
                if pos_vector is not None:
                    kwargs['pos_vector'] = pos_vector
                elif 'pos_vector' not in kwargs:
                    kwargs['pos_vector'] = 0
                if 'rotation' in kwargs:
                    orient_type_temp = 'Axis'
                    orient_amount_temp = [kwargs['rotation'].magnitude(), kwargs['rotation'].normalize()]
                super(MovingRefFrame, self).__init__(name, dim=3, kwargs['pos_vector'], position_coord='rect',\
                                                     orient_type_temp, orient_amount_temp, orient_order)
            #Set translational params as functions of time
            if pos_vector is not None:
                self._trans_acc, self._trans_vel, self._pos_vector = get_motion_pos(pos_vector, parentframe)
            elif trans_vel is not None:
                trans_vel = self._to_vector(trans_vel)
                for x in ('posvalue', 't'):
                    if x not in kwargs:
                        kwargs[x] = 0
                self._trans_acc, self._trans_vel, self._pos_vector = \
                                 get_motion_vel(trans_vel, kwargs['posvalue'], kwargs['t'], parentframe)
            elif trans_acc is not None:
                trans_acc = self._to_vector(trans_acc)
                for x in ('velvalue', 'posvalue', 't1', 't2'):
                    if x not in kwargs:
                        kwargs[x] = 0
                self._trans_acc, self._trans_vel, self._pos_vector = \
                                 get_motion_vel(trans_acc, kwargs['velvalue'], kwargs['posvalue'], \
                                                kwargs['t1'], kwargs['t2'], parentframe)
            else:
                self._trans_acc, self._trans_vel, self._pos_vector = get_motion_pos(0, parentframe)
            #Set rotational params as functions of time
            if orient_type is not None or (orient_type is None and ang_vel is None and ang_acc is None):
                super(MovingRefFrame, self).__init__(name, dim=3, self._pos_vector, position_coord='rect',\
                                                     orient_type, orient_amount, orient_order, wrt)
                #Set angular velocity and angular accln params by time-differentiation of DCM
                dcm2diff = self.dcm(parentframe)
                diffed = dcm2diff.diff(self.time)
                angvelmat = diffed * dcm2diff.T
                w1 = trigsimp(expand(angvelmat[7]), recursive=True)
                w2 = trigsimp(expand(angvelmat[2]), recursive=True)
                w3 = trigsimp(expand(angvelmat[3]), recursive=True)
                self._ang_vel = w1 * parentframe.basis(0) + w2 * parentframe.basis(1) + \
                                w3 * parentframe.basis(2)
                self._ang_acc = parentframe.time_derivative(self._ang_vel)
            elif ang_vel is not None:
                ang_vel = self._to_vector(ang_vel)
                for x in ('rotation', 'rt'):
                    if x not in kwargs:
                        kwargs[x] = 0
                self._ang_acc, self._ang_vel, rotation = \
                                 get_motion_vel(ang_vel, kwargs['rotation'], kwargs['rt'], parentframe)
                angle = rotation.magnitude()
                axis = rotation.normalize()
                super(MovingRefFrame, self).__init__(name, dim=3, self._pos_vector, \
                                                     position_coord='rect', orient_type = 'Axis', \
                                                     orient_amount = [angle, axis], wrt = parentframe)
            elif ang_acc is not None:
                ang_acc = self._to_vector(ang_acc)
                for x in ('ang_vel', 'rotation', 'rt1', 'rt2'):
                    if x not in kwargs:
                        kwargs[x] = 0
                self._ang_acc, self._ang_vel, rotation = \
                                 get_motion_acc(ang_vel, kwargs['ang_vel'], kwargs['rotation'], \
                                                kwargs['rt1'], kwargs['rt2'], parentframe)
                angle = rotation.magnitude()
                axis = rotation.normalize()
                super(MovingRefFrame, self).__init__(name, dim=3, self._pos_vector, \
                                                     position_coord='rect', orient_type = 'Axis', \
                                                     orient_amount = [angle, axis], wrt = parentframe)

    def _to_vector(self, inputparam):
        """
        Converts input given by the user into a vector.

        Input may be a vector, or a tuple defining a vector using two components - one in own's basis vectors,
        and the other in some other frame
        """
        if type(inputparam) != tuple:
            return inputparam
        else:
            #Process tuple to get a vector
            outvect = sum([inputparam[0][i] * self.basis(i) for i in range(self.dim)])
            outvect += inputparam[1]
            return outvect
        
    def _frame_path(self, otherframe):
        """
        Calculates 'path' of frames starting from this frame to the other,
        along with the index of the common root

        Returns index, list pair
        """
        if self._root != otherframe._root:
            raise ValueError("No connecting path between the two frames-"+ str(self) + \
                             " and " + str(otherframe))
        other_path = []
        frame = otherframe
        while frame.parent is not None:
            other_path.append(frame)
            frame = frame.parent
        other_path.append(frame)
        frameset = set(other_path)
        self_path = []
        frame = self
        while frame not in frameset:
            self_path.append(frame)
            frame = frame.parent
        index = len(self_path)
        i = other_path.index(frame)
        while i >= 0:
            self_path.append(other_path[i])
            i -= 1
        return index, self_path
            
            
    def convert_pos_vector(self, pos_vector, frame=None):
        """
        Convert a position vector defined in another frame to this frame

        The position vector must be defined in a single frame.

        If pos_vector = 0, frame will need to be specified. In this case, the other
        frame's origin's position vector wrt this frame will be returned.

        Parameters
        ==========

        pos_vector : vector
            The position vector to be converted

        frame : MovingRefFrame
            Frame whose origin's position vector has to be calculated

        Examples
        ========

        >>> from sympy.physics.mechanics import MovingRefFrame
        >>> R1 = MovingRefFrame('R1', parentframe=None)
        >>> R2 = MovingRefFrame('R2', parentframe=R1, pos_vector = 2 * R1.basis(0), ang_vel = R1.basis(2))
        >>> pos_vector = 5 * R1.basis(0) + 3 * R1.basis(1) + 4 * R1.basis(2)
        >>> R2.convert_pos_vector(R1)
        ...
        """

        if pos_vector == 0 and type(frame) != MovingRefFrame:
            raise ValueError("Valid frame has to be specified for zero vector")
        pos_vector = sympify(pos_vector)
        #Check if pos_vector is entirely defined in a single frame
        if pos_vector.is_vector:
            if type(pos_vector) == Vector:
                frame = pos_vector.coord_sys
            elif pos_vector != 0:
                if type(pos_vector) == VectAdd:
                    frame = pos_vector.args[0].coord_sys
                else:
                    frame = pos_vector.coord_sys
                for x in condition.atoms():
                    if type(x) == Vector or type(x) == BaseScalar:
                        if x.coord_sys != frame:
                            raise ValueError("Position vector must be defined in a single frame")
        else:
            raise ValueError("pos_vector must be a valid vector")
        #Convert
        pos_vector = self.express(pos_vector)
        return self.express(frame.pos_vector_in(self)) + pos_vector

    @cacheit
    def pos_vector_in(self, otherframe):
        """
        Returns the relative position vector of this frame's origin in
        another frame.

        Parameters
        ==========

        otherframe : MovingRefFrame
            The frame to calculate the position vector in

        Examples
        ========

        >>> from sympy.physics.mechanics import MovingRefFrame
        >>> R1 = MovingRefFrame('R1', parentframe=None)
        >>> R2 = MovingRefFrame('R2', parentframe=R1, ang_vel = R1.basis(2))
        >>> R3 = MovingRefFrame('R3', parentframe = R2, pos_vector = R2.basis(0))
        >>> R3.pos_vector_in(R1)
        #Something equivalent to R2.x
        """

        if otherframe == self:
            return 0
        elif otherframe == self.parent:
            return self._pos_vector
        rootindex, path = self._frame_path(otherframe)
        result = 0
        i = -1
        for i in range(rootindex):
            result += path[i]._pos_vector
        i += 2
        while i < len(path):
            result -= path[i]._pos_vector
            i += 1
        return result
        
    @cacheit
    def trans_vel_in(self, otherframe):
        """
        Returns the relative translational velocity vector of this frame's origin in
        another frame.

        Parameters
        ==========

        otherframe : MovingRefFrame
            The frame to calculate the relative velocity in

        Examples
        ========

        >>> from sympy.physics.mechanics import MovingRefFrame
        >>> R1 = MovingRefFrame('R1', parentframe=None)
        >>> R2 = MovingRefFrame('R2', parentframe=R1, ang_vel = R1.basis(2))
        >>> R3 = MovingRefFrame('R3', parentframe = R2, pos_vector = R2.basis(0))
        >>> R3.trans_vel_in(R1)
        #Something equivalent to R2.y
        """
        if otherframe == self:
            return 0
        elif otherframe == self.parent:
            return self._trans_vel
        return otherframe.time_derivative(self.pos_vector_in(otherframe))
        
    @cacheit
    def trans_acc_in(self, otherframe):
        """
        Returns the relative translational acceleration vector of this frame's origin in
        another frame.

        Parameters
        ==========

        otherframe : MovingRefFrame
            The frame to calculate the relative acceleration in

        Examples
        ========

        >>> from sympy.physics.mechanics import MovingRefFrame
        >>> R1 = MovingRefFrame('R1', parentframe=None)
        >>> R2 = MovingRefFrame('R2', parentframe=R1, ang_vel = R1.basis(2))
        >>> R3 = MovingRefFrame('R3', parentframe = R2, pos_vector = R2.basis(0))
        >>> R3.trans_vel_in(R1)
        #Something equivalent to -R1.x
        """
        if otherframe == self:
            return 0
        elif otherframe == self.parent:
            return self._trans_acc
        return otherframe.time_derivative(self.pos_vector_in(otherframe), 2)
    
    @cacheit
    def ang_vel_in(self, otherframe):
        """
        Returns the relative angular velocity vector of this frame in
        another frame.

        Parameters
        ==========

        otherframe : MovingRefFrame
            The frame to calculate the relative angular velocity in

        Examples
        ========

        ToBeDone
        """
        if otherframe == self:
            return 0
        elif otherframe == self.parent:
            return self._ang_vel
        rootindex, path = self._frame_path(otherframe)
        result = 0
        for i in range(rootindex):
            result += path[i]._ang_vel
        i += 2
        while i < len(path):
            result -= path[i]._ang_vel
            i += 1
        return result
    
    @cacheit
    def ang_acc_in(self, otherframe):
        """
        Returns the relative angular acceleration vector of this frame in
        another frame.

        Parameters
        ==========

        otherframe : MovingRefFrame
            The frame to calculate the relative angular acceleration in

        Examples
        ========

        ToBeDone
        """
        if otherframe == self:
            return 0
        elif otherframe == self.parent:
            return self._ang_acc
        rootindex, path = self._frame_path(otherframe)
        result = 0
        for i in range(rootindex):
            result += path[i]._ang_acc
        i += 2
        while i < len(path):
            result -= path[i]._ang_acc
            i += 1
        return result
    
    def time_derivative(self, expr, order=1):
        """
        Calculate the time derivative of a field function in this frame.

        References
        ==========

        http://en.wikipedia.org/wiki/Rotating_reference_frame#Time_derivatives_in_the_two_frames

        Parameters
        ==========

        expr : vector/scalar Expr
            The field whose time derivative is to be calculated

        order : integer
            The order of the derivative to be calculated

        Examples
        ========

        ToBeDone
        """
        
        expr = sympify(expr)
        if order == 0:
            return expr
        if order%1 != 0 or order < 0:
            raise ValueError("Unsupported value of order entered")
        if expr.is_vector:
            #Decompose vector into its constituents in each frame
            frame_dict = {}
            if type(expr) == Vector:
                frame_dict[expr.coord_sys] = expr
            elif type(expr) == VectMul:
                frame_dict[expr.coord_sys] = expr.coord_sys.express(expr)
            else:
                for x in expr.args:
                    if x.coord_sys not in frame_dict:
                        frame_dict[x.coord_sys] = 0
                    frame_dict[x.coord_sys] += x.coord_sys.express(expr)
            #Process each constituent separately, and add to get result
            dt = 0
            for frame in frame_dict:
                if frame == self:
                    dt += diff(frame_dict[frame], self.time)
                else:
                    dt += diff(frame_dict[frame], self.time) + \
                          frame.ang_vel_in(self).cross(frame_dict[frame])
            return self.time_derivative(dt, order-1)
        else:
            return diff(self.express(expr), self.time, order)
