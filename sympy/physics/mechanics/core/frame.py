from sympy.vector import CoordSys, Vector, VectAdd, VectMul

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
        
        if cls.global_time is None:
            cls.global_time = time_var
        else:
            #Global time is already defined
            #Raise error to maintain immutability
            raise ValueError("Global time already defined!")

    def __init__(name, dim, position, position_coord, orient_type, orient_amount,
                 orient_order, wrt, trans_vel=(0, None), trans_acc=(0, None),
                 ang_vel=(0, None), ang_acc=(0, None)):
        """
        Initializer for the MovingRefFrame class.

        All the motion parameters must be specified as a (vectorial value,frame) couple.

        User must take care not to create inconsistent conditions through the parameters
        passed.
        """

        if self.global_time is None:
            raise ValueError("Time variable has not been set!")
        #The line below will change, have to incorporate the possibility of different
        #kinds of coordinate systems
        super(MovingCoordSys, self).__init__(name, dim, position, position_coord,
                                             orient_type, orient_amount, orient_order, wrt)
        #translational parameters will always be stored wrt a global 'static' frame
        #Maybe check whther only one of the three-position, velocity, and acceleration are defined?
        #The above idea may avoid any inconsistencies
        if trans_vel[1] is None:
            self._trans_vel = self.time_derivative(self.position) + wrt._trans_vel
        else:
            self._trans_vel = trans_vel[0] + trans_vel[1]._trans_vel
        #If translational acceleration is not provided, it is taken as the time derivative
        #of the translational velocity of this frame in the frame it is defined wrt
        if trans_acc[1] is None:
            if trans_vel[1] is None:
                self._trans_acc = self.time_derivative(self.position, 2) + wrt._trans_acc
            else:
                self._trans_acc = trans_vel[1].time_derivative(self._trans_vel)
        else:
            self._trans_acc = trans_vel[0] + trans_vel[1]._trans_vel

    def map_variables(self, otherframe):
        """
        Returns a dictionary mapping the base scalars of this frame to the base scalars
        of another frame

        Every base scalar of this frame is expressed as a function of the base scalars of
        the other frame, and the global time variable

        Parameters
        ==========

        otherframe : MovingRefFrame
            The other frame to express the base scalars of this frame in

        Examples
        ========

        """

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
        return self._trans_vel - otherframe._trans_vel
        
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
        return self._trans_acc - otherframe._trans_acc
    
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
    
    def time_derivative(self, vector, order=1):
        """
        Calculate the time derivative of a vector in this frame.

        References
        ==========

        http://en.wikipedia.org/wiki/Rotating_reference_frame#Time_derivatives_in_the_two_frames

        Parameters
        ==========

        vector : Vector/VectAdd/VectMul
            The vector whose time derivative is to be calculated

        order : integer
            The order of the derivative to be calculated
        """
        #Decompose VectAdd/VectMul into its components in each constituent frame
        #Then express each constituent in this frame, then add
        if type(vector) == Vector:
            if vector.system() == self:
                return 0
            else:
                return (vector.system.ang_vel_in(self)).cross(vector)
        #VectMul and VectAdd cases yet to be done
            
