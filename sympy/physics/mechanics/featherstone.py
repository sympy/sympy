# NOTE: THIS CODE WILL BE MAKING USE OF THE JOINTS CODE THAT STILL NEEDS TO BE
# FIXED AND ADDED

from sympy import Matrix, symbols, zeros
# from sympy.physics.mechanics import SymbolicSystem
from sympy.physics.vector.spatial import cross_m, cross_f


class FeatherstonesMethod(object):
    """Class docstring
    Discuss the class itself and method

    Will use articulated body algorithm described in chapter 7 and detailed in
    Table 7.1 on page 132
    Attributes
    ==========

    discuss necessary inputs
    Parameters
    ==========

    TODO: Should the class simply make Denavit-Hartenberg parameters for all
    of the joints where num_joints=num_bodies for an open chain system?
    coordinates : ordered iterables of functions of time

    lambda : list, len(n)
        parent array

    Notes
    =====

    TODO: determine the variable that featherstone uses for number of
    bodies, number of degrees of freedoms, number of joints

    Examples
    ========
    Double pendulum example

        >>> from sympy import symbols
        >>> from sympy.physics.mechanics import Body, FeatherstonesMethod, \
        ...     PinJoint

        >>> bodyA = Body('A')
        >>> bodyB = Body('B')
        >>> ground_body = Body('G')

        >>> joint1 = PinJoint('joint1', ground, bodyA,
        ...                   parent_point_pos=(0, 0, 0),
        ...                   child_point_pos=(0, 1, 0),
        ...                   parent_axis=ground.frame.z,
        ...                   child_axis=bodyA.frame.z)
        >>> joint2 = PinJoint('joint2', bodyA, bodyB,
        ...                   parent_point_pos=(0, 0, 0),
        ...                   child_point_pos=(0, 1, 0),
        ...                   parent_axis=bodyA.frame.z,
        ...                   child_axis=bodyB.frame.z)

        >>> featherstone = FeatherstonesMethod(ground)

    References
    ==========

    .. [Featherstone2007] Roy Featherstone, Rigid Body Dynamics Algorithms. 2007
       Springer
    """

    def __init__(self, base, tau, grav_accn=None):
        """init docstring"""

        # Build the kinematic tree from the base body.  Bodies hold joint
        # objects information, joints hold parent and child body information.
        # Should this be incorporated into pass 1?

        # Definition of variables
        #   vi - Spatial velocity of body i
        #   XJ - Joint transform
        #   Si - Motion subspace allowed by joint i
        #   vJ - Joint spatial velocity
        #   i_X_lam - Transform from frame attached to the parent of body to
        #             body i
        #   XT - Transform from frame for body i to the location of the joint on
        #        body i
        #   ci - Intermediate term for body motion (Sdot * qdot)
        #   i_X_o - Transform from the base to body i
        #   IA - Articulated body inertia
        #   I - Inertia of body i
        #   pA - Articulated body bias force
        #   i_X_o_star - Transform from the base to body i in the Force subspace
        #   U - Intermediate term
        #   D - Intermediate term
        #   u - Intermediate term
        #   fx - Extermal forces on body i (pg 130)
        #   tau - Force term from cannonical form of the equations of motion
        #   Ia - Apparent inertia of a massless handle
        #   pa - Apparent bias force of a massless handle

        bodies = [base]  # base is the input body
        parent_array = [0]  # The base does not have a parent object
        indx = 0
        while True:

            # If the last body in the bodies list does not have child joints
            # exit the while loop
            if bodies[indx].child_joints == [] and bodies[indx] == bodies[-1]:
                break

            # Add child bodies of current index to the bodies array and parent
            # arrays
            if bodies[indx].child_joints != []:
                for joint in bodies[indx].child_joints:
                    bodies.append(joint.child)
                    parent_array.append(indx)

            indx += 1


        # Initialize variables for the algorithm
        num_bodies = len(bodies)
        iters = list(range(num_bodies))
        v = [0 for i in iters]
        c = [0 for i in iters]
        S = [0 for i in iters]
        IA = [0 for i in iters]
        pA = [0 for i in iters]
        i_X_lam = [0 for i in iters]
        # i_X_o = [0 for i in iters]

        # Parent Defaults
        v[0] = Matrix([0, 0, 0, 0, 0, 0])

        for i in iters[1:]:
            # jcalc information should be obtained from the joint objects,
            # whether by attributes or a method
            XJ, S[i], vJ = bodies[i].parent_joint.spatial_info()
            i_X_lam[i] = XJ*bodies[i].parent_joint.XT_child()
            # if parent_array[i] != 0:
            #     i_X_o[i] = i_X_lam[i] * i_X_o[parent_array[i]]
            # else:
            #     i_X_o[i] = i_X_lam[i]
            v[i] = i_X_lam[i] * v[parent_array[i]] + vJ
            c[i] = cross_m(v[i])*vJ
            temp = bodies[i].parent_joint.child.masscenter
            vec = bodies[i].parent_joint.child_joint_point.pos_from(temp)
            IA[i] = bodies[i].spatial_inertia(vec)
            pA[i] = cross_f(v[i]) * bodies[i].spatial_inertia(vec) * v[i]

        # Pass 2 mock up (note only being shown here and not in combined_pass_1
        # because it should come out the same)

        # Initialize pass 2 variables
        U = [0 for i in iters]
        D = [0 for i in iters]
        u = [0 for i in iters]

        for i in reversed(iters[1:]):
            U[i] = IA[i] * S[i]
            D[i] = S[i].transpose() * U[i]
            u[i] = tau[i] - S[i].transpose()*pA[i]
            if parent_array[i] != 0:
                Ia = IA[i] - U[i]*D[i].inv()*U[i].transpose()
                pa = pA[i] + Ia*c[i] + U[i]*D[i].inv()*u[i]
                lam_X_i_star = i_X_lam[i].transpose()  # I believe
                IA[parent_array[i]] = IA[parent_array[i]] + lam_X_i_star[i] * \
                    Ia*i_X_lam[i]
                pA[parent_array[i]] = pA[parent_array[i]] + lam_X_i_star * pa

        # Pass 3 mock up (note only being shown here and not in combined_pass_1
        # because it should come out the same)

        # Initialize variables for pass 3
        a = [0 for i in iters]
        qddot = [0 for i in iters]
        # Default gravity in negative y direction? Give user option to turn
        # gravitational effect off?
        if grav_accn is None:
            g = symbols('g')
            a[0] = Matrix([0, 0, 0, 0, - g, 0])
        else:
            a[0] = Matrix([0, 0, 0, grav_accn[0], grav_accn[1], grav_accn[2]])

        for i in iters[1:]:
            a[i] = i_X_lam[i]*a[parent_array[i]] + c[i]
            qddot[i] = D[i].inv()*(u[i] - U[i].transpose()*a[i])
            a[i] += S[i]*qddot[i]

    def separate_pass_1(self, base, tau, grav_accn=None):
        """This method is being used to prototype what the algorithm looks
        like if the kinematic tree is built separately from the first pass
        in the algorithm

        base - root of the kinematic chain, link 0

        tau - vector of force variables

        grav_accn - 3x1 vector representing the gravitational acceleration in
        the base link frame"""

        bodies = [base]  # base is the input body
        parent_array = [0]  # The base does not have a parent object
        while True:
            indx = 0

            # If the last body in the bodies list does not have child joints
            # exit the while loop
            if bodies[indx].child_joints == [] and bodies[indx] == bodies[-1]:
                break

            # Add child bodies of current index to the bodies array and parent
            # arrays
            if bodies[indx].child_joints != []:
                for joint in bodies[indx].child_joints:
                    bodies.append(joint.child_body)
                    parent_array.append(indx)

            indx += 1

        # Initialize variables for the algorithm
        num_bodies = len(bodies)
        iters = [range(num_bodies)]
        v = Matrix(iters)
        c = Matrix(iters)
        S = Matrix(iters)
        IA = [0 for i in iters]
        pA = [0 for i in iters]
        i_X_lam = [0 for i in iters]
        # i_X_o = [0 for i in iters]

        v[0] = 0
        for i in iters:
            # jcalc information should be obtained from the joint objects,
            # whether by attributes or a method
            XJ, S[i], vJ = bodies[i].parent_joint.spatial_info
            i_X_lam[i] = XJ*bodies[i].parent_joint.XT_child
            # if parent_array[i] != 0:
            #     i_X_o[i] = i_X_lam[i] * i_X_o[parent_array[i]]
            # else:
            #     i_X_o[i] = i_X_lam[i]
            v[i] = i_X_lam[i] * v[parent_array[i]] + vJ
            c[i] = cross_m(v[i])*vJ
            IA[i] = bodies[i].I  # Need to convert body inertia to spatial inertia? (3x3 -> 6x6)
            pA[i] = cross_f(v[i]) * bodies[i].I * v[i]

        # Pass 2 mock up (note only being shown here and not in combined_pass_1
        # because it should come out the same)

        # Initialize pass 2 variables
        U = [0 for i in iters]
        D = [0 for i in iters]
        u = [0 for i in iters]

        for i in reversed(iters):
            U[i] = IA[i] * S[i]
            D[i] = S[i].transpose() * U[i]
            u[i] = tau[i] - S[i].transpose()*pA[i]
            if parent_array[i] != 0:
                Ia = IA[i] - U[i]*D[i].inv()*U[i].transpose()
                pa = pA[i] + Ia*c[i] + U[i]*D[i].inv()*u[i]
                lam_X_i_star = i_X_lam[i].transpose()  # I believe
                IA[parent_array[i]] = IA[parent_array[i]] + lam_X_i_star[i] * \
                    Ia*i_X_lam[i]
                pA[parent_array[i]] = pA[parent_array[i]] + lam_X_i_star * pa

        # Pass 3 mock up (note only being shown here and not in combined_pass_1
        # because it should come out the same)

        # Initialize variables for pass 3
        a = Matrix(iters)
        # Default gravity in negative y direction? Give user option to turn
        # gravitational effect off?
        if grav_accn is None:
            g = symbols('g')
            a[0] = Matrix([0, 0, 0, 0, - g, 0])
        else:
            a[0] = Matrix([0, 0, 0, grav_accn[0], grav_accn[1], grav_accn[2]])

        for i in iters:
            a[i] = i_X_lam[i]*a[parent_array[i]] + c[i]
            qddot = D[i].inv()*(u[i] - U[i].transpose()*a[i])
            a[i] += S[i]*qddot[i]

    def combined_pass_1(self, base):
        """This method is being used to prototype what the algorithm looks
        like if the kinematic tree is built alongside the first pass
        in the algorithm"""

        # Set up variables for the first iteration
        bodies = [base]  # base is the input body
        parent_array = [0]  # The base does not have a parent object
        S = [[0]]  # Will need a more proper zero subspace matrix?
        i_X_lam = [zeros(6)]
        i_X_o = [zeros(6)]
        v = [0]
        c = [0]
        IA = [base.inertia]
        pA = [0]

        while True:
            indx = 0

            # If the last body in the bodies list does not have child joints
            # exit the while loop
            if bodies[indx].child_joints == [] and bodies[indx] == bodies[-1]:
                break

            # Add child bodies of current index to the bodies array and parent
            # arrays
            for joint in bodies[indx].child_joints:
                bodies.append(joint.child_body)
                parent_array.append(indx)

                # jcalc information should be obtained from the joint objects,
                # whether by attributes or a method
                XJ, S_temp, vJ = joint.spatial_info
                S.append(S_temp)
                i_X_lam.append(XJ * joint.XT_child)
                if joint.parent_body != base:
                    i_X_o.append(i_X_lam[-1] * i_X_o[parent_array[-1]])
                else:
                    i_X_o.append(i_X_lam)
                v.append(i_X_lam[-1] * v[parent_array[-1]] + vJ)
                c.append( cross_m(v[-1])*vJ)
                IA.append(joint.child_body.I)
                # Need to find i_X_o_star from i_X_o
                i_X_o_star = i_X_o[-1]  # Need actual equation for transformation
                pA.append(cross_f(v[-1])*joint.child_body.I*v[-1] -
                          i_X_o_star*joint.child_body.fx)

            indx += 1

    def to_system(self):
        """This will return a SymbolicSystem containing the equations of motion
        information along with other system information such as the body
        objects"""

        # Algorithm forms combined explicit right hand side or just dynamic
        # explicit right hand side. If the latter what are the kinematic's rhs?

        coordinates = Matrix([i.coordinates[0] for i in bodies])
        return SymbolicSystem(coordinates, right_hand_side, bodies=self.bodies)
