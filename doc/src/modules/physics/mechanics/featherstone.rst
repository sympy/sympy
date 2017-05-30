Simple 2 link pendulum Example
==============================

::

    >>> from sympy import symbols
    >>> from sympy.physics.mechanics import Body, FeatherstonesMethod, PinJoint
    >>> bodyA = Body('A')
    >>> bodyB = Body('B')
    >>> ground_body = Body('G')
    >>> joint1 = PinJoint('joint1', ground, bodyA, parent_point_pos=(0, 0, 0),
    ...                   child_point_pos=(0, 1, 0), parent_axis=ground.frame.z,
    ...                   child_axis=bodyA.frame.z)
    >>> joint2 = PinJoint('joint2', bodyA, bodyB, parent_point_pos=(0, 0, 0),
    ...                   child_point_pos=(0, 1, 0), parent_axis=bodyA.frame.z,
    ...                   child_axis=bodyB.frame.z)
    >>> featherstone = FeatherstonesMethod(ground)

Some expected return values. ::

    >>> featherstone.coordinates
    Matrix([joint1_theta, joint2_theta])
    >>> featherstone.speeds
    Matrix([joint1_omega, joint2_omega])
    >>> symsystem = featherstone.to_system()
    >>> symsystem
    <SymbolicSystem object at ...>
