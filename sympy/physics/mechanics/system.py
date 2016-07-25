class SymbolicSystem(object):
    """SymbolicSystem is a class that contains all the information about a
    system in a symbolic format such as the equations of motions and the bodies
    and loads in the system.

    There are three ways that the equations of motion can be described for
    Symbolic System:


        [1] Explicit form where the kinematics and dynamics are combined
            x' = F(x, t, r, p)

        [2] Implicit form where the kinematics and dynamics are combined
            M(x, p) x' = F(x, t, r, p)

        [3] Implicit form where the kinematics and dynamics are separate
            M(q, p) u' = F(q, u, t, r, p)
            q' = G(q, u, t, r, p)

    where

        x : states, i.e. [q, u]
        t : time
        r : specified (exogenous) inputs
        p : constants
        q : generalized coordinates
        u : generalized speeds
        M : mass matrix or generalized inertia matrix (implicit dynamical
            equations or combined dynamics and kinematics)
        F : right hand side (implicit dynamical equations or combined dynamics
            and kinematics)
        G : right hand side of the kinematical differential equations



    Attributes
    ==========

    TODO - come up with the size of the matrices and add a notes section
    detailing any variables I used

    coordinates : Matrix
        This is a matrix containing the generalized coordinates of the system

    speeds : Matrix
        This is a matrix containing the generalized speeds of the system

    states : Matrix
        This is a matrix containing the state variables of the system

    alg_con : List
        This list contains the indices of the algebraic constraints in the
        combined equations of motion. The presence of these constraints requires
        that a DAE solver be used instead of an ODE solver

    dyn_implicit_mat : Matrix
        This is the M matrix in form [3] of the equations of motion (the mass
        matrix or generalized inertia matrix of the dynamical equations of
        motion in implicit form).

    dyn_implicit_rhs : Matrix
        This is the F vector in form [3] of the equations of motion (the right
        hand side of the dynamical equations of motion in implicit form).

    comb_implicit_mat : Matrix
        This is the M matrix in form [2] of the equations of motion. This matrix
        contains a block diagonal structure where the top left block (the first
        rows) represent the matrix in the implicit form of the kinematical
        equations and the bottom right block (the last rows) represent the
        matrix in the implicit form of the dynamical equations.

    comb_implicit_rhs : Matrix
        This is the F vector in form [2] of the equations of motion. The top
        part of the vector represents the right hand side of the implicit form
        of the kinemaical equations and the bottom of the vector represents the
        right hand side of the implicit form of the dynamical equations of
        motion.

    comb_explicit_rhs : Matrix
        This vector represents the right hand side of the combined equations of
        motion in explicit form (form [1] from above).

    kin_explicit_rhs : Matrix
        This is the right hand side of the explicit form of the kinematical
        equations of motion as can be seen in form [3] (the G matrix).

    dynamic_symbols : Matrix
        This matrix contains all of the symbols that have a time dependence in
        the equations of motion of the system

    constant_symbols : Matrix
        This matrix contains all of the symbols that do not have a time
        dependence in the equations of motion of the system

    output_eqns : Dictionary
        If output equations were given they are stored in a dictionary where the
        key corresponds to the name given for the specific equation and the
        value is the equation itself in symbolic form

    bodies : List
        If the bodies in the system were given they are stored in a list for
        future access

    loads : List
        If the loads in the system were given they are stored in a list for
        future access. This includes forces and torques where forces are given
        by (point of application, force vector) and torques are given by
        (reference frame acted upon, torque vector).

    Parameters
    ==========

    coord_states : ordered iterable of functions of time
        This input will either be a collection of the coordinates or states of
        the system depending on whether or not the speeds are also given. If
        speeds are specified this input will be assumed to be the coordinates
        otherwise this input will be assumed to be the states.

    right_hand_side : Matrix
        This variable is the right hand side of the equations of motion in any
        of the forms. The specific form will be assumed depending on whether a
        mass matrix or coordinate derivatives are given.

    speeds : ordered iterable of functions of time, optional
        This is a collection of the generalized speeds of the system. If given
        it will be assumed that the first argument (coord_states) will represent
        the generalized coordinates of the system.

    mass_matrix : Matrix, optional
        The matrix of the implicit forms of the equations of motion (forms [2]
        and [3]). The distinction between the forms is determined by whether or
        not the coordinate derivatives are passed in. If they are given form [3]
        will be assumed otherwise form [2] is assumed.

    coordinate_derivatives : Matrix, optional
        The right hand side of the kinematical equations in explicit form. If
        given it will be assumed that the equations of motion are being entered
        in form [3].

    alg_con : Iterable, optional
        The indexes of the rows in the equations of motion that contain
        algebraic constraints instead of differential equations. If the
        equations are input in form [3], it will be assumed the indexes are
        referencing the mass_matrix/right_hand_side combination and not the
        coordinate_derivatives.

    output_eqns : Dictionary, optional
        Any output equations that are desired to be tracked are stored in a
        dictionary where the key corresponds to the name given for the specific
        equation and the value is the equation itself in symbolic form

    coord_idxs : Iterable, optional
        If coord_states corresponds to the states rather than the coordinates
        this variable will tell SymbolicSystem which indexes of the states
        correspond to generalized coordinates.

    speed_idxs : Iterable, optional
        If coord_states corresponds to the states rather than the coordinates
        this variable will tell SymbolicSystem which indexes of the states
        correspond to generalized speeds.

    bodies : List, optional
        List containing the bodies of the system

    loads : List, optional
        List containing the loads of the system where forces are given by (point
        of application, force vector) and torques are given by (reference frame
        acting upon, torque vector). Ex [(point, force), (ref_frame, torque)]

    """

    def __init__(self):
        """init method"""

    @property
    def coordinates(self):
        """Returns the column matrix of the generalized coordinates"""

    @property
    def speeds(self):
        """Returns the column matrix of generalized speeds"""

    @property
    def states(self):
        """Returns the column matrix of the state variables"""

    @property
    def alg_con(self):
        """Returns a list with the indices of the rows containing algebraic
        constraints in the combined form of the equations of motion"""

    @property
    def dyn_implicit_mat(self):
        """Returns the matrix, M, corresponding to the dynamic equations in
        implicit form, M x' = F, where the kinematical equations are not
        included"""

    @property
    def dyn_implicit_rhs(self):
        """Returns the column matrix, F, corresponding to the dynamic equations
        in implicit form, M x' = F, where the kinematical equations are not
        included"""

    @property
    def comb_implicit_mat(self):
        """Returns the matrix, M, corresponding to the equations of motion in
        implicit form (form [2]), M x' = F, where the kinematical equations are
        included"""

    @property
    def comb_implicit_rhs(self):
        """Returns the column matrix, F, corresponding to the equations of
        motion in implicit form (form [2]), M x' = F, where the kinematical
        equations are included"""

    @property
    def comb_explicit_rhs(self):
        """Returns the right hand side of the equations of motion in explicit
        form, x' = F, where the kinematical equations are included"""

    @property
    def kin_explicit_rhs(self):
        """Returns the right hand side of the kinematical equations in explicit
        form, q' = G"""

    @property
    def dynamic_symbols(self):
        """Returns a column matrix containing all of the symbols in the system
        that depend on time"""

    @property
    def constant_symbols(self):
        """Returns a column matrix containing all of the symbols in the system
        that do not depend on time"""

    @property
    def output_eqns(self):
        """Returns a dictionary of the given output equations where the key's
        are user defined names for each equation and the values are the symbolic
        equations themselves"""

    @property
    def bodies(self):
        """Returns the bodies in the system"""

    @property
    def loads(self):
        """Returns the loads in the system"""
