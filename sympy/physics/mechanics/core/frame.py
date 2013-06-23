from sympy import Symbol, diff, sympify
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

    def __init__(name, dim, pos_vector=None, trans_vel=None, trans_acc=None, orient_type,
                 orient_amount, orient_order, rot_vel=None, rot_acc=None, parentframe=None, **kwargs):
        """
        Initializer for the MovingRefFrame class.
        """

        #WIP

    def convert_pos_vector(self, pos_vector, frame=None):
        """
        Convert a position vector defined in another frame to this frame

        The position vector must be defined in a single frame.

        If pos_vector = 0, frame will need to be specified. In this case, the other
        frame's position vector wrt this frame will be returned.

        Parameters
        ==========

        pos_vector : vector
            The position vector to be converted

        frame : MovingRefFrame
            Frame whose origin's position vector has to be calculated

        Examples
        ========
        
        """

        if pos_vector == 0 and type(frame) != MovingRefFrame:
            raise ValueError("Valid frame has to be specified for zero vector")
        pos_vector = sympify(pos_vector)
        #Check if pos_vector is entirely defined in a single frame
        if pos_vector.is_vector:
            if type(pos_vector) == Vector:
                frame = pos_vector.system
            elif pos_vector != 0:
                if type(pos_vector) == VectAdd:
                    frame = pos_vector.args[0].system
                else:
                    frame = pos_vector.system
                for x in condition.atoms():
                    if type(x) == Vector or type(x) == BaseScalar:
                        if x.system != frame:
                            raise ValueError("Position vector must be defined in a single frame")
        else:
            raise ValueError("pos_vector must be a valid vector")
        #Convert
        pos_vector = self.express(pos_vector)
        return frame.pos_vector_in(self) + pos_vector

    @cacheit
    def pos_vector_in(self, otherframe):
        """
        Returns the relative position vector of this frameis origin in
        otherframe, expressed in otherframe.

        Parameters
        ==========

        otherframe : MovingRefFrame
            The frame to calculate the position vector in

        Examples
        ========

        ToBeDone
        """
        return otherframe.express(self._pos_vector - otherframe._pos_vector)

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
            The frame to calculate the relative acceleration in

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
            The frame to calculate the relative angular velocity in

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
            The frame to calculate the relative angular acceleration in

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
        
        expr = sympify(expr)
        if order == 0:
            return expr
        if order%1 != 0 or order < 0:
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
