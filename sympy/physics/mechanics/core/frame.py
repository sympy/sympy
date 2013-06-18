from sympy.vector import CoordSys

class MovingRefFrame(CoordSys):
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
        
        if global_time == None:
            cls.global_time = time_var
        else:
            #Global time is already defined
            #Raise error to maintain immutability
            raise ValueError("Global time already defined!")

    def __new__(cls, ...):
        #Check if global time variable has been defined
        #If not, raise error to maintain immutability of instances
        if cls.global_time == None:
            raise ValueError("Global time variable has not been defined!")
        else:
            super.__new__(cls, ...)

    def __init__(name, dim, bases, variables, trans_vel=(0, None), trans_acc=(0, None),
                 ang_vel=(0, None), ang_acc=(0, None)):
        """
        Initializer for the MovingRefFrame class.

        All the motion parameters must be specified as a (vectorial value,frame) couple
        """

        #translational parameters will always be stored wrt a global 'static' frame
        if trans_vel[1] == None:
            self._trans_vel = 0
        else:
            self._trans_vel = trans_vel[0] + trans_vel[1]._trans_vel
        #If translational acceleration is not provided, it is taken as the time derivative
        #of the translational velocity of this frame in the frame it is defined wrt
        if trans_acc[1] == None:
            if trans_vel[1] == None:
                self._trans_acc = 0
            else:
                self._trans_acc = trans_vel[1].time_derivative(self._trans_vel)
        else:
            self._trans_acc = trans_vel[0] + trans_vel[1]._trans_vel
        #Will depend on CoordSys API
    
    def trans_vel_in(self, otherframe):
        """
        Returns the relative translational velocity vector of this frame in
        otherframe.

        Parameters
        ==========

        otherframe : MovingRefFrame
            The frame to calculate the relative velocity in

        Examples
        ========
        """
        
    def trans_acc_in(self, otherframe):
        """
        Returns the relative translational acceleration vector of this frame in
        otherframe.

        Parameters
        ==========

        otherframe : MovingRefFrame
            The frame to calculate the relative velocity in

        Examples
        ========
        """
    
    def ang_vel_in(self, otherframe):
        """
        Returns the relative angular velocity vector of this frame in
        otherframe.

        Parameters
        ==========

        otherframe : MovingRefFrame
            The frame to calculate the relative velocity in

        Examples
        ========
        """
    
    def ang_acc_in(self, otherframe):
        """
        Returns the relative angular acceleration vector of this frame in
        otherframe.

        Parameters
        ==========

        otherframe : MovingRefFrame
            The frame to calculate the relative velocity in

        Examples
        ========
        """
    
    def time_derivative(self, vector):
        """
        Calculate the time derivative of a vector in this frame.

        References
        ==========

        http://en.wikipedia.org/wiki/Rotating_reference_frame#Time_derivatives_in_the_two_frames

        Parameters
        ==========

        vector : Vector/VectAdd/VectMul
        """
        #Decompose VectAdd/VectMul into its components in each constituent frame
        #Then express each constituent in this frame, then add
