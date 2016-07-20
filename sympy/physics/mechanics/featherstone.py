# This file will serve as a view of how featherstone method will look if the
# featherstone class init initialized a system model and the equations of motion
# are found using a forward_dynamics method just like is done in
# http://cctbx.sourceforge.net/scitbx_rigid_body_essence/featherstone.py

# NOTE: THIS CODE WILL BE MAKING USE OF THE JOINTS CODE THAT STILL NEEDS TO BE
# FIXED AND ADDED


class Featherstone(object):
    """Class docstring
    Discuss the class itself and method
    Attributes
    ==========

    Examples
    ========
    Double pendulum example

        >>> from sympy import symbols
        >>> from sympy.physics.mechanics import Body, FeatherstonesMethod,
        ...     inertia, PinJoint, Point, ReferenceFrame

        >>> mA, mB = symbols('mA mB')
        >>> mass_center_A = Point('moA')
        >>> mass_center_B = Point('moB')
        >>> ground_point = Point('O')

        >>> frameA = ReferenceFrame('A')
        >>> frameB = ReferenceFrame('B')
        >>> ground = ReferenceFrame('G')

        >>> body_A_inertia, body_B_inertia = symbols('I_Az I_Bz')
        >>> body_A_inertia_dyadic = inertia(frameA, 0, 0, body_A_inertia)
        >>> body_B_inertia_dyadic = inertia(frameB, 0, 0, body_B_inertia)
        >>> ground_inertia_dyadic = inertia(ground, 0, 0, 0)

        >>> bodyA = Body('A', mass_center_A, mA, frameA, body_A_inertia_dyadic)
        >>> bodyB = Body('B', mass_center_B, mB, frameB, body_B_inertia_dyadic)
        >>> ground_body = Body('G', ground_point, 0, ground,
        ...                    ground_inertia_dyadic)

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
    """

    def __init__(self):
        """init docstring
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
        """
