from sympy import diff, sympify
from sympy.core.cache import cacheit
from sympy.vector import CoordSys, Vector, VectAdd, VectMul, BaseScalar

class MovingRefFrame(CoordSysRect): #For now, I have subclassed CoordSysRect
    """
    A moving frame of reference in classical mechanics.

    It subclasses the static CoordSys class of vector module
    and adds motion-related functionality to it.

    This frame can have translational and angular motion with
    respect to another frame, which affects the expression of
    vectors/points defined in this frame in other frames, and
    vice versa.
    """
    global_time = None

    @classmethod
    def set_time(cls, time_var):
        """
        Set the global time variable

        All the frames use this common Symbol as the time variable
        for motion related calculations.

        Parameters
        ==========

        time_var : Symbol
            The time variable to be shared by all frames
        """
        
        if cls.global_time is None:
            cls.global_time = time_var
        else:
            #Global time is already defined
            #Raise error to maintain immutability
            raise ValueError("Global time already defined!")

    def __init__(name, dim, position_coord, translation=[0, None, None], orient_type,
                 orient_amount, orient_order, rotation=[None, None], parentframe=None):
        """
        Initializer for the MovingRefFrame class.
        """

        if self.global_time is None:
            raise ValueError("Time variable has not been set!")
        #motion parameters will always be stored wrt a global 'static' frame
        if parentframe is None:
            #Frame is effectively the 'global fixed frame'
            self._parent = None
            self._pos = 0
            self._trans_vel = 0
            self._trans_acc = 0
            self._ang_vel = 0
            self._ang_acc = 0
            super(MovingCoordSys, self).__init__(name, dim, wrt=None)
        else:
            #Check if args are right
            if type(parentframe) != MovingRefFrame:
                raise ValueError("Parent frame must be a MovingRefFrame")
            if len(translation) != 3:
                raise ValueError("'translation' must be a list of 3 vector values!")
            if translation[0] is None:
                translation[0] = 0
            for x in translation:
                x = sympify(x)
                try:
                    if x is not None:
                        if x != 0 and !(x.is_vector):
                            raise ValueError
                except:
                    raise ValueError("Unsupported parameters in 'translation'")
            if len(rotation) != 2:
                raise ValueError("'rotation' must be a list of 2 vector values!")
            for x in rotation:
                x = sympify(x)
                try:
                    if x is not None:
                        if x != 0 and !(x.is_vector):
                            raise ValueError
                except:
                    raise ValueError("Unsupported parameters in 'rotation'")
            #Call superclass initializer
            super(MovingCoordSys, self).__init__(name, dim, translation[0], position_coord,
                                                 orient_type, orient_amount, orient_order wrt=parentframe)
            #Motion initialization
            self._parent = parentframe
            superparent = parentframe
            while superparent._parent is not None:
                superparent = superparent._parent
            self._pos = superparent.express(translation[0]) + parentframe._pos
            self._trans_vel = parentframe._trans_vel
            if translation[1] is None:
                self._trans_vel += diff(self._pos, global_time)
            else:
                self._trans_vel += superparent.express(translation[1])
            self._trans_acc = parentframe._trans_acc
            if translation[2] is None:
                self._trans_acc += diff(self._trans_vel, global_time)
            else:
                self._trans_acc += superparent.express(translation[2])
            self._ang_vel = parentframe._ang_vel
            if rotation[0] is None:
                #TODO: Find angular velocity using orientation params
            else:
                self._ang_vel += superparent.express(rotation[0])
            self._ang_acc = parentframe._trans_acc
            if rotation[1] is None:
                self._ang_acc += diff(self._ang_vel, global_time)
            else:
                self._ang_accln += superparent.express(rotation[1])

    @cacheit
    def trans_vel_in(self, otherframe):
        """
        Returns the relative translational velocity vector of this frame in
        otherframe, expressed in otherframe.

        Parameters
        ==========

        otherframe : MovingRefFrame
            The frame to calculate the relative velocity in

        Examples
        ========

        ToBeDone
        """
        return otherframe.express(self._trans_vel - otherframe._trans_vel)
        
    @cacheit
    def trans_acc_in(self, otherframe):
        """
        Returns the relative translational acceleration vector of this frame in
        otherframe, expressed in otherframe.

        Parameters
        ==========

        otherframe : MovingRefFrame
            The frame to calculate the relative velocity in

        Examples
        ========

        ToBeDone
        """
        return otherframe.express(self._trans_acc - otherframe._trans_acc)
    
    @cacheit
    def ang_vel_in(self, otherframe):
        """
        Returns the relative angular velocity vector of this frame in
        otherframe, expressed in otherframe.

        Parameters
        ==========

        otherframe : MovingRefFrame
            The frame to calculate the relative velocity in

        Examples
        ========

        ToBeDone
        """
        return otherframe.express(self._ang_vel - otherframe._ang_vel)
    
    @cacheit
    def ang_acc_in(self, otherframe):
        """
        Returns the relative angular acceleration vector of this frame in
        otherframe, expressed in otherframe.

        Parameters
        ==========

        otherframe : MovingRefFrame
            The frame to calculate the relative velocity in

        Examples
        ========

        ToBeDone
        """
        return otherframe.express(self._ang_acc - otherframe._ang_acc)
    
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
        
        if order == 0:
            return expr
        if order%1 != 1 or order < 0:
            raise ValueError("Unsupported value of order entered")
        if expr.is_vector:
            #Decompose vector into its constituents in each frame
            frame_dict = {}
            if type(expr) == Vector:
                frame_dict[expr.system] = expr
            elif type(expr) == VectMul:
                frame_dict[expr.system] = expr.system.express(expr)
            else:
                for x in expr.args:
                    if x.system not in frame_dict:
                        frame_dict[x.system] = 0
                    frame_dict[x.system] += x.system.express(expr)
            #Process each constituent separately, and add to get result
            dt = 0
            for frame in frame_dict:
                dt += diff(frame_dict[frame], global_time) + \
                      frame.ang_vel_in(self).cross(frame_dict[frame])
            return self.time_derivative(dt, order-1)
        else:
            return diff(self.express(expr), global_time, order)
