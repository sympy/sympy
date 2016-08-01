# NOTE: THIS CODE WILL BE MAKING USE OF THE JOINTS CODE THAT STILL NEEDS TO BE
# FIXED AND ADDED

from sympy import Matrix, zeros
from sympy.physics.vector import cross_m, cross_f


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

    def __init__(self, base):
        """init docstring"""

        # Build the kinematic tree from the base body.  Bodies hold joint
        # objects information, joints hold parent and child body information.
        # Should this be incorporated into pass 1?

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
            for joint in bodies[indx].child_joints:
                bodies.append(joint.child_body)
                parent_array.append(indx)

            indx += 1

        num_bodies = len(bodies)

        # Definition of variables
        #   vi - Spatial velocity of body i
        #   XJ - Joint transform
        #   Si - Motion subspace allowed by joint i
        #   vJ - Joint spatial velocity
        #   cJ - Intermediate term for joint motion (Sdot * qdot)
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
        #   fx - Extermal forces on body i
        #   U - Intermediate term
        #   D - Intermediate term
        #   u - Intermediate term
        #   tau - Force term from cannonical form of the equations of motion
        #   Ia - Apparent inertia of a massless handle
        #   pa - Apparent bias force of a massless handle

        # iters = [range(num_bodies)]

        # Algorithm as presented in Featherstone's book on page 132
        # # Pass 1

        # # Will probably need to build an idea of the kinematic tree in the
        # # first loop
        # v[0] = 0
        # for i in iters:
        #     # jcalc information should be obtained from the joint objects,
        #     # whether by attributes or a method
        #     XJ, S[i], vJ, cJ = jcalc(jtype, q[i], qdot[i])
        #     i_X_lam[i] = XJ*XT[i]
        #     if lam[i] != 0:
        #         i_X_o = i_X_lam[i] * i_X_o[lam[i]]
        #     v[i] = i_X_lam[i]*vlam + vJ
        #     c[i] = cJ + cross_m(v[i], vJ)
        #     IA[i] = I[i]
        #     # Need to find i_X_o_star from i_X_o
        #     pA[i] = cross_f(v[i], I[i]*v[i]) - i_X_o_star*fx[i]

        # # Pass 2

        # for i in reversed(iters):
        #     U[i] = IA[i] * S[i]
        #     D[i] = S[i].transpose() * U[i]
        #     u[i] = tau[i] - S[i].transpose()*pA[i]
        #     if lam[i] != 0:
        #         Ia = IA[i] - U[i]*D[i].inv()*U[i].transpose()
        #         pa = pA[i] + Ia*c[i] + U[i]*D[i].inv()*u[i]
        #         IA[lam[i]] = IA[lam[i]] + lam_X_i_star[i]*Ia*i_X_lam[i]
        #         pA[lam[i]] = pA[lam[i]] + lam_X_i_star * pa

        # # Pass 3

        # a[0] = -ag
        # for i in iters:
        #     a_prime = i_X_lam[i]*alam + c[i]
        #     qddot = D[i].inv()*(u[i] - U[i].transpose()*a_prime)
        #     a[i] = a_prime + S[i]*qddot[i]

    def separate_pass_1(self, base):
        """This method is being used to prototype what the algorithm looks
        like if the kinematic tree is built separately from the first pass
        in the algorithm"""

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
        i_X_o = [0 for i in iters]

        v[0] = 0
        for i in iters:
            # jcalc information should be obtained from the joint objects,
            # whether by attributes or a method
            XJ, S[i], vJ, cJ = bodies[i].parent_joint.spatial_info
            i_X_lam[i] = XJ*bodies[i].parent_joint.XT_child
            if parent_array[i] != 0:
                i_X_o[i] = i_X_lam[i] * i_X_o[parent_array[i]]
            else:
                i_X_o[i] = i_X_lam[i]
            v[i] = i_X_lam[i] * v[parent_array[i]] + vJ
            c[i] = cJ + cross_m(v[i], vJ)
            IA[i] = bodies[i].I
            # Need to find i_X_o_star from i_X_o
            i_X_o_star = i_X_o[i]  # Need actual equation for transformation
            pA[i] = cross_f(v[i], bodies[i].I*v[i]) - i_X_o_star*bodies[i].fx

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
                XJ, S_temp, vJ, cJ = joint.spatial_info
                S.append(S_temp)
                i_X_lam.append(XJ * joint.XT_child)
                if joint.parent_body != base:
                    i_X_o.append(i_X_lam[-1] * i_X_o[parent_array[-1]])
                else:
                    i_X_o.append(i_X_lam)
                v.append(i_X_lam[-1] * v[parent_array[-1]] + vJ)
                c.append(cJ + cross_m(v[-1], vJ))
                IA.append(joint.child_body.I)
                # Need to find i_X_o_star from i_X_o
                i_X_o_star = i_X_o[-1]  # Need actual equation for transformation
                pA.append(cross_f(v[-1], joint.child_body.I*v[-1]) -
                          i_X_o_star*joint.child_body.fx)

            indx += 1
