from sympy import Symbol, diff, sympify, trigsimp
from sympy.core.cache import cacheit
from sympy.vector import CoordSys, Vector, VectAdd, VectMul, BaseScalar
from sympy.physics.mechanics import get_motion_acc, get_motion_vel, \
     get_motion_pos, dynamicsymbols


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

    def __init__(name, pos_vector=None, trans_vel=None, trans_acc=None, \
                 orient_type=None, orient_amount=None, orient_order=None, \
                 ang_vel=None, ang_acc=None, parentframe=None, **kwargs):
        """
        Initializer for the MovingRefFrame class.

        Note
        ====

        For a class whose parentframe == None, the 'timevar' keyword arg must be
        given to initialize a Symbol forvtime. If not specified, it falls back to
        default Symbol('t').

        Translational motion can be specified either by specifying the position
        vector(pos_vector), velocity (trans_vel) or acceleration(trans_acc) as
        functions of time. Required boundary conditions can be specified as keyword
        arguments. The relevant keyword-args are -
        
        For velocity- 'pos_vector_b', 't' - the position at the value of time t
        For acceleration - 'trans_vel_b', 'pos_vector_b', 't1', 't2' - the velocity
        and position boundary conditions at times t1 and t2 respectively

        If more than one arguments are entered, preference is given in order-
        position, velocity, acceleration.
        Hence, entering both velocity and acceleration will result in motion being
        set as per the velocity conditions.

        Rotational motion can be specified either by specifying the orientation(refer
        CoordSys docs) or the angular velocity(ang_vel) or the angular acceleration
        (ang_acc). Boundary conditions are required, similar to the translational
        motion params.

        For angular velocity- 'rotation_b', 'rt' - the rotation at the value of
        time rt.
        For angular acceleration - 'ang_vel_b', 'rotation_b', 'rt1', 'rt2' -
        the velocity and position boundary conditions at times rt1 and rt2
        respectively

        If more than one arguments are entered, preference is given in order-
        orientation, velocity, acceleration.
        Hence, entering both angular velocity and acceleration will result in motion
        being set as per the velocity conditions.

        Any of the time-dependent arguments(apart from pos_vector and orientation)
        can be specified as a vector, or as a tuple of the form - ([v1, v2, v3], V)
        where V is a vector defined wrt parent frame, and v1, v2, v3 are components of
        a vector V' defined wrt this frame itself. The total vectorial value of the
        param would be V + V'.

        Boundary conditions must be expressed entirely in the parentframe. If one
        or more of the boundary condition parameters are not entered, they are taken
        to be zero automatically.

        If any boundary condition is dependent on time, the corresponding value of time
        will be substituted to convert it to the time-independent form.

        Parameters
        ==========

        parentframe : MovingRefFrame
            The parent of the frame to be initialized. If it's the first frame,
            parentframe = None

        name : String
            The name of the frame

        pos_vector : vector
            The position vector of the new frame wrt its parent

        trans_vel : vector/tuple
            The translational velocity of the new frame wrt its parent

        trans_acc : vector/tuple
            The translational acceleration of the new frame wrt its parent

        orient_type , orient_amount, orient_order : orientation params
            Refer CoordSys docs

        ang_vel : vector/tuple
            The angular velocity of the ne frame in its parent

        ang_acc : vector/tuple
            The angular acceleration of the new frame wrt its parent

        kwargs
        ------

        'timevar' : Symbol
            The desired Symbol for global time

        't', 'rt', 't1', 't2', 'rt1', rt2' : sympifiable
            The times at which boundary conditions are specified. Must be
            independent of time variable
        
        'pos_vector_b' : vector
            Boundary condition for position

        'trans_vel_b' : vector
            Boundary condition for velocity

        'rotation_b' : vector
            Boundary condition for rotation

        'ang_vel_b' : vector
            Boundary condition for angular velocity

        Examples
        ========

        >>> #Initialise first frame
        >>> F1 = MovingRefFrame('F1')
        >>> #Initialise F2 having trans vel of 3 * F1.x and initial pos (2, 0, 0) in F1
        >>> F2 = MovingRefFrame(parentframe = F1, name = 'F2',
            trans_vel = 3 * F1.basis(0), pos_vector_b = 2 * F1.basis(0))
        
        """

        #Special case of no parent frame
        if parentframe is None:
            self._parent = None
            self._root = self
            if 'timevar' not in kwargs:
                #Default time variable is the Symbol used by
                #dynamicsymbols
                self._time = dynamicsymbols._t
            else:
                #Set global time
                kwargs['timevar'] = sympify(kwargs['timevar'])
                if not kwargs['timevar'].is_Symbol:
                    raise TypeError("timevar must be a Symbol")
                self._time = kwargs['timevar']
            #All motion params zero
            self._pos_vector = 0
            self._trans_vel = 0
            self._trans_acc = 0
            self._ang_vel = 0
            self._ang_acc = 0
            super(MovingRefFrame, self).__init__(name, dim=3, wrt = None)
        else:
            if not isinstance(parentframe, MovingRefFrame):
                raise TypeError("parentframe should be of type MovingRefFrame")
            self._parent = parentframe
            self._root = parentframe._root
            #Take time variable from parent frame
            self._time = parentframe.time
            #Initial orientation and positioning in case of tuple args
            flag = False
            for x in (trans_vel, trans_acc, ang_vel, ang_acc):
                if type(x) == tuple:
                    flag = True
                    break
            if flag:
                #Time values for boundary conditions must all be zero if tuple args
                #are to be processed. If not, raise exception.
                for x in ('t', 'rt', 't1', 't2', 'rt1', 'rt2'):
                    if x in kwargs:
                        if kwargs[x] != 0:
                            raise ValueError("initial values(t=0) needed for\
                                             initialization with tuple args")
                #Fix initial position from given args/ kwargs
                if pos_vector is not None:
                    kwargs['pos_vector_b'] = pos_vector
                elif 'pos_vector_b' not in kwargs:
                    kwargs['pos_vector_b'] = 0
                self._pos_vector = kwargs['pos_vector_b']
                #Fix initial orientation from given args/ kwargs
                orient_type_temp = None
                orient_amount_temp = None
                if orient_type is not None:
                    orient_type_temp = orient_type
                    orient_amount_temp = orient_amount
                elif 'rotation_b' in kwargs:
                    orient_type_temp = 'Axis'
                    orient_amount_temp = [kwargs['rotation'].magnitude(),
                                          kwargs['rotation'].normalize()]
                #Calling the superclass' __init__ function for the first time sets
                #the initial position and orientation of this frame wrt parent frame.
                #This allows usage of this frame's basis vectors in setting velocites/
                #acceleration
                super(MovingRefFrame, self).__init__(name, dim=3,
                                                     position=kwargs['pos_vector'],
                                                     position_coord='rect',
                                                     orient_type=orient_type_temp,
                                                     orient_amount=orient_amount_temp,
                                                     orient_order=orient_order)
            #Set translational params as functions of time
            if pos_vector is not None:
                #User has provided pos_vector, hence set motion according to that
                self._trans_acc, self._trans_vel, self._pos_vector = \
                                 get_motion_pos(pos_vector, parentframe)
            elif trans_vel is not None:
                #User has provided trans_vel. Process rest of translational
                #motion params from this
                trans_vel = self._to_vector(trans_vel)
                for x in ('pos_vector_b', 't'):
                    if x not in kwargs:
                        kwargs[x] = 0
                self._trans_acc, self._trans_vel, self._pos_vector = \
                                 get_motion_vel(trans_vel, kwargs['pos_vector_b'],
                                                kwargs['t'], parentframe)
            elif trans_acc is not None:
                #User has provided trans_acc. Process rest of translational motion params
                #using this
                trans_acc = self._to_vector(trans_acc)
                for x in ('trans_vel_b', 'pos_vector_b', 't1', 't2'):
                    if x not in kwargs:
                        kwargs[x] = 0
                self._trans_acc, self._trans_vel, self._pos_vector = \
                                 get_motion_acc(trans_acc, kwargs['trans_vel_b'],
                                                kwargs['pos_vector_b'],
                                                kwargs['t1'], kwargs['t2'], parentframe)
            else:
                #None of the params are provided. This means this frame has no
                #translational motion wrt parent
                self._trans_acc, self._trans_vel, self._pos_vector = (0, 0, 0)
            #Set rotational params as functions of time
            if orient_type is not None or (orient_type is None and ang_vel is None and \
                                           ang_acc is None):
                super(MovingRefFrame, self).__init__(name, dim=3,
                                                     position=self._pos_vector,
                                                     position_coord='rect',
                                                     orient_type, orient_amount,
                                                     orient_order, wrt)
                #Set angular velocity and angular accln params by
                #time-differentiation of DCM
                dcm2diff = self.dcm(parentframe)
                diffed = dcm2diff.diff(self._time)
                angvelmat = diffed * dcm2diff.T
                w1 = trigsimp(expand(angvelmat[7]), recursive=True)
                w2 = trigsimp(expand(angvelmat[2]), recursive=True)
                w3 = trigsimp(expand(angvelmat[3]), recursive=True)
                self._ang_vel = -w1 * parentframe.basis(0) - w2 * parentframe.basis(1) -\
                                w3 * parentframe.basis(2)
                self._ang_acc = parentframe.dt(self._ang_vel)
            elif ang_vel is not None:
                #ang_vel is provided. Process other rotation params using that
                #(and boundary conditions)
                ang_vel = self._to_vector(ang_vel)
                for x in ('rotation_b', 'rt'):
                    if x not in kwargs:
                        kwargs[x] = 0
                self._ang_acc, self._ang_vel, rotation = \
                                 get_motion_vel(ang_vel, kwargs['rotation_b'],
                                                kwargs['rt'], parentframe)
                angle = rotation.magnitude()
                axis = rotation.normalize()
                super(MovingRefFrame, self).__init__(name, dim=3,
                                                     position=self._pos_vector,
                                                     position_coord='rect',
                                                     orient_type = 'Axis',
                                                     orient_amount = [angle, axis],
                                                     wrt = parentframe)
            elif ang_acc is not None:
                #ang_acc is provided. Process other rotation params using that
                #and boundary conditions
                ang_acc = self._to_vector(ang_acc)
                for x in ('ang_vel_b', 'rotation_b', 'rt1', 'rt2'):
                    if x not in kwargs:
                        kwargs[x] = 0
                self._ang_acc, self._ang_vel, rotation = \
                                 get_motion_acc(ang_vel, kwargs['ang_vel_b'],
                                                kwargs['rotation_b'],
                                                kwargs['rt1'], kwargs['rt2'],
                                                parentframe)
                angle = rotation.magnitude()
                axis = rotation.normalize()
                super(MovingRefFrame, self).__init__(name, dim=3,
                                                     position=self._pos_vector,
                                                     position_coord='rect',
                                                     orient_type = 'Axis', 
                                                     orient_amount = [angle, axis],
                                                     wrt = parentframe)

    @property
    def parent(self):
        return self._parent

    @property
    def time(self):
        return self._time

    def _to_vector(self, inputparam):
        """
        Converts input given by the user into a vector.

        Input may be a vector, or a tuple defining a vector using two components -
        one in own's basis vectors, and the other in some other frame
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
            raise ValueError("No connecting path between the two frames-"+ \
                             str(self) + " and " + str(otherframe))
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
            
    def convert_pos_vector(self, pos_vector, frame):
        """
        Convert a position vector defined in another frame to this frame

        Parameters
        ==========

        pos_vector : vector
            The position vector to be converted

        frame : MovingRefFrame
            The frame the given position vector is defined in

        Examples
        ========

        >>> from sympy.physics.mechanics import MovingRefFrame
        >>> R1 = MovingRefFrame('R1', parentframe=None)
        >>> R2 = MovingRefFrame('R2', parentframe=R1, pos_vector = \
                                2 * R1.basis(0), ang_vel = R1.basis(2))
        >>> pos_vector = 5 * R1.basis(0) + 3 * R1.basis(1) + 4 * R1.basis(2)
        >>> R2.convert_pos_vector(R1)
        ...
        """
        
        pos_vector = sympify(pos_vector)
        if pos_vector.is_Vector or pos_vector == 0:
            #Add the given position vector and the specified frame's
            #position vector in this frame
            return pos_vector + frame.pos_vector_in(self)
        else:
            raise TypeError("pos_vector must be a valid vector")

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
        Returns the relative translational velocity vector of this frame's
        origin in another frame.

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
        elif otherframe.parent == self:
            return -1 * self._trans_vel
        return otherframe.dt(self.pos_vector_in(otherframe))
        
    @cacheit
    def trans_acc_in(self, otherframe):
        """
        Returns the relative translational acceleration vector of this frame's
        origin in another frame.

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
        elif otherframe.parent == self:
            return -1 * self._trans_acc
        return otherframe.dt(self.trans_vel_in(otherframe))
    
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
        elif otherframe.parent == self:
            return -1 * self._ang_vel
        rootindex, path = self._frame_path(otherframe)
        result = 0
        i = -1
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
        elif otherframe.parent == self:
            return -1 * self._ang_acc
        return otherframe.dt(self.ang_vel_in(otherframe))
    
    def dt(self, expr, order=1):
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
            frame_dict = vector.separate()
            #Process each constituent separately, and add to get result
            dt = 0
            for frame in frame_dict:
                if frame == self:
                    dt += diff(frame_dict[frame], self.time)
                else:
                    dt += diff(frame_dict[frame], self.time) + \
                          frame.ang_vel_in(self).cross(frame_dict[frame])
            return self.dt(dt, order-1)
        else:
            return diff(self.express(expr), self.time, order)
